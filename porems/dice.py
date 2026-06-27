################################################################################
# Dice Class                                                                   #
#                                                                              #
"""Separation of a molecule object into smaller cubes for pair-search."""
################################################################################


import math
import multiprocessing as mp

import numpy as np
import porems.geometry as geometry


def _neighbor_ids(cube_id, count, size, is_pbc):
    """Return valid neighboring cube ids for *cube_id* (module-level for pickling)."""
    def step(idx, dim, delta):
        idx = list(idx)
        idx[dim] += delta
        if idx[dim] >= count[dim]:
            idx[dim] = 0 if is_pbc else None
        elif idx[dim] < 0:
            idx[dim] = count[dim] - 1 if is_pbc else None
        return None if None in idx else tuple(idx)

    neighbors = []
    for dz in (-1, 0, 1):
        nz = step(cube_id, 2, dz)
        if nz is None:
            continue
        for dy in (-1, 0, 1):
            ny = step(nz, 1, dy)
            if ny is None:
                continue
            for dx in (-1, 0, 1):
                nx = step(ny, 0, dx)
                if nx is not None:
                    neighbors.append(nx)
    return neighbors


def _find_bond_task(cube_list, atom_data, pointer, count, size, is_pbc, mol_box_arr, atom_type, distance):
    """Module-level worker for parallel bond search.

    Accepts only the minimal data needed to avoid pickling the full Molecule object.
    """
    d_lo, d_hi = distance
    t0, t1 = atom_type
    bond_list = []

    for cube_id in cube_list:
        neighbors = _neighbor_ids(cube_id, count, size, is_pbc)
        atoms = [aid for nb in neighbors for aid in pointer[nb]]

        for atom_id_a in pointer[cube_id]:
            if atom_data[atom_id_a][0] != t0:
                continue
            entry = [atom_id_a, []]
            pos_a = atom_data[atom_id_a][1]
            for atom_id_b in atoms:
                if atom_data[atom_id_b][0] != t1 or atom_id_a == atom_id_b:
                    continue
                bv = pos_a - atom_data[atom_id_b][1]
                mask = np.abs(bv) > 3 * size
                bv[mask] -= mol_box_arr[mask] * np.round(bv[mask] / mol_box_arr[mask])
                length = math.sqrt(bv[0]*bv[0] + bv[1]*bv[1] + bv[2]*bv[2])
                if d_lo <= length <= d_hi:
                    entry[1].append(atom_id_b)
            bond_list.append(entry)

    return bond_list


class Dice:
    """This class splits the molecule into smaller sub boxes and provides
    parallelized functions for pair-search.

    The aim is reducing the workload on search algorithms for atom pairs.
    Normally the computational effort would be

    .. math::

        \\mathcal O(n^2)

    with the number of atoms :math:`n`, because each atom has to be compared
    with all other atoms.

    Since a bond distance is fixed, the idea is reducing the search space by
    dividing the molecule box into cubes and only performing the  pair-search
    within these smaller boxes and their 26 immediate neighbors. Assuming that
    the grid structure is ideal in a geometrical sense, that all bond length and
    angles are constant, the number of atoms in each cube are a constant
    :math:`b`. The computational effort for each atom is thus a constant

    .. math::

        \\mathcal O(27\\cdot b^2)

    Therefore, the computational effort for an entire search scales linear with
    the number of cubes. For example, doubling the cristobalite block size only
    increases the effort eightfold.

    Furthermore, the search is easily parallelizable, since no communication is
    needed between the subprocesses that each cover a set of cubes. The effort
    therefore has an ideal speedup.

    Note that the cube size must be strictly greater than the intended
    bond length searches.

    Parameters
    ----------
    mol : Molecule
        Molecule to be divided
    size : float
        Cube edge size
    is_pbc : bool
        True if periodic boundary conditions are needed
    """
    def __init__(self, mol, size, is_pbc):
        # Initialize
        self._dim = 3
        self._np = mp.cpu_count()

        self._mol = mol
        self._size = size
        self._is_pbc = is_pbc

        self._atom_data = {atom_id: [atom.get_atom_type(), np.asarray(self._mol.pos(atom_id))] for atom_id, atom in enumerate(self._mol.get_atom_list())}
        self._mol_box = self._mol.get_box()
        self._mol_box_arr = np.asarray(self._mol_box)

        # Split molecule box into cubes and fill them with atom ids
        self._split()
        self._fill()


    ##############
    # Management #
    ##############
    def _split(self):
        """Here the number of cubes is calculated for each dimension for the
        defined cube size and molecule dimension. A dictionary of cubes is
        generated containing the coordinates of the origin point of each cube.
        Furthermore, an empty list for each cube is added, that will contain
        atom ids of atom objects.

        Cube ids are tuples containing three elements with the x-axis as the
        first entry, y-axis as the second entry and the z-axis as the third
        entry

        .. math::

            \\text{id}=\\begin{pmatrix}x&y&z\\end{pmatrix}.
        """
        # Calculate number of cubes in each dimension
        self._count = [math.floor(box/self._size) for box in self._mol_box]

        # Fill cube origins
        self._origin = {}
        self._pointer = {}
        for i in range(self._count[0]):
            for j in range(self._count[1]):
                for k in range(self._count[2]):
                    self._origin[(i, j, k)] = [self._size*x for x in [i, j, k]]
                    self._pointer[(i, j, k)] = []

    def _pos_to_index(self, position):
        """Calculate the cube index for a given position.

        Parameters
        ----------
        position : list
            Three-dimensional coordinates

        Returns
        -------
        index : tuple
            Cube index
        """
        return tuple([math.floor(pos/self._size) if math.floor(pos/self._size)<self._count[dim] else self._count[dim]-1 for dim, pos in enumerate(position)])

    def _fill(self):
        """Based on their coordinates, the atom ids, as defined in the molecule
        object, are filled into the cubes.
        """
        for atom_id in range(self._mol.get_num()):
            self._pointer[self._pos_to_index(self._mol.pos(atom_id))].append(atom_id)


    ############
    # Iterator #
    ############
    def _step(self, dim, step, index):
        """Helper function for iterating through the cubes. Optionally,
        periodic boundary conditions are applied.

        Parameters
        ----------
        dim : integer
            Stepping dimension
        step : integer
            Step to move
        index : list
            Cube index

        Returns
        -------
        index : list
            New cube index
        """
        # Step in intended dimension
        index = list(index)
        index[dim] += step

        # Periodicity
        if index[dim] >= self._count[dim]:
            index[dim] = 0 if self._is_pbc else None
        elif index[dim] < 0:
            index[dim] = self._count[dim]-1 if self._is_pbc else None

        return tuple(index)

    def _right(self, index):
        """Step one cube to the right considering the x-axis.

        Parameters
        ----------
        index : list
            Cube index

        Returns
        -------
        index : list
            New cube index
        """
        return self._step(0, 1, index)

    def _left(self, index):
        """Step one cube to the left considering the x-axis.

        Parameters
        ----------
        index : list
            Cube index

        Returns
        -------
        index : list
            New cube index
        """
        return self._step(0, -1, index)

    def _top(self, index):
        """Step one cube to the top considering the y-axis.

        Parameters
        ----------
        index : list
            Cube index

        Returns
        -------
        index : list
            New cube index
        """
        return self._step(1, 1, index)

    def _bot(self, index):
        """Step one cube to the bottom considering the y-axis.

        Parameters
        ----------
        index : list
            Cube index

        Returns
        -------
        index : list
            New cube index
        """
        return self._step(1, -1, index)

    def _front(self, index):
        """Step one cube to the front considering the z-axis.

        Parameters
        ----------
        index : list
            Cube index

        Returns
        -------
        index : list
            New cube index
        """
        return self._step(2, 1, index)

    def _back(self, index):
        """Step one cube to the back considering the z-axis.

        Parameters
        ----------
        index : list
            Cube index

        Returns
        -------
        index : list
            New cube index
        """
        return self._step(2, -1, index)

    def neighbor(self, cube_id, is_self=True):
        """Get the ids of the cubes surrounding the given one.

        Parameters
        ----------
        index : list
            Main cube index
        is_self : bool, optional
            True to add the main cube to the output

        Returns
        -------
        neighbor : list
            List of surrounding cube ids, optionally including given cube id
        """
        # Initialize
        neighbor = []

        # Find neighbors
        z = [self._back(cube_id), cube_id, self._front(cube_id)]
        y = [[self._top(i), i, self._bot(i)] for i in z]

        for i in range(len(z)):
            for j in range(len(y[i])):
                neighbor.append(self._left(y[i][j]))
                neighbor.append(y[i][j])
                neighbor.append(self._right(y[i][j]))

        return [n for n in neighbor if n is not None and (is_self or n != cube_id)]


    ##########
    # Search #
    ##########
    def find_bond(self, cube_list, atom_type, distance):
        """Single-threaded bond search used as a reference implementation.

        Searches for atom-pairs that fulfill the distance requirements within the
        given cube and all 26 surrounding ones. For parallel execution see
        :func:`find_parallel`, which calls the module-level :func:`_find_bond_task`
        to avoid pickling the full :class:`Molecule` object.

        Parameters
        ----------
        cube_list : list
            List of cube indices to search in, use an empty list for all cubes
        atom_type : list
            List of two atom types, e.g. ``["Si", "O"]``
        distance : list
            Bounds of allowed distance ``[lower, upper]`` in nm

        Returns
        -------
        bond_list : list
            Bond array containing lists of two atom ids
        """
        # Process input
        cube_list = cube_list if cube_list else list(self._pointer.keys())

        # Loop through all given cubes
        bond_list = []
        for cube_id in cube_list:
            # Get atom ids of surrounding cubes
            atoms = sum([self._pointer[x] for x in self.neighbor(cube_id) if None not in x], [])

            # Run through atoms in the main cube
            for atom_id_a in self._pointer[cube_id]:
                # Check type
                if self._atom_data[atom_id_a][0] == atom_type[0]:
                    entry = [atom_id_a, []]
                    # Search in all surrounding cubes for partners
                    for atom_id_b in atoms:
                        if self._atom_data[atom_id_b][0] == atom_type[1] and not atom_id_a == atom_id_b:
                            # Calculate bond vector with nearest-image convention
                            bv = self._atom_data[atom_id_a][1] - self._atom_data[atom_id_b][1]
                            mol_box = self._mol_box_arr
                            mask = np.abs(bv) > 3*self._size
                            bv[mask] -= mol_box[mask] * np.round(bv[mask] / mol_box[mask])
                            length = geometry.length(bv)
                            # Check if bond distance is within range
                            if length >= distance[0] and length <= distance[1]:
                                entry[1].append(atom_id_b)

                    # Add pairs to bond list
                    bond_list.append(entry)

        return bond_list

    def find_parallel(self, cube_list, atom_type, distance):
        """Parallelized bond search using :func:`_find_bond_task`.

        Divides the cube list across ``cpu_count()`` workers. Workers receive only
        the atom-data and pointer dicts rather than the full
        :class:`~porems.Molecule`, reducing IPC overhead. On single-threaded
        Unix processes a ``fork`` start-context is used so workers inherit the
        parent's memory via copy-on-write (zero-copy). ``spawn`` is used on
        Windows and inside multi-threaded processes (e.g. test runners) where
        ``fork`` risks deadlocks.

        Parameters
        ----------
        cube_list : list
            List of cube indices to search in, use an empty list for all cubes
        atom_type : list
            List of two atom types, e.g. ``["Si", "O"]``
        distance : list
            Bounds of allowed distance ``[lower, upper]`` in nm

        Returns
        -------
        bond_list : list
            Bond array containing lists of two atom ids
        """
        # Process input
        cube_list = cube_list if cube_list else list(self._pointer.keys())

        # Divide cubes on processors
        cube_num = math.floor(len(cube_list)/self._np)
        cube_np = [cube_list[cube_num*i:] if i == self._np-1 else cube_list[cube_num*i:cube_num*(i+1)] for i in range(self._np)]

        # Pass only the essential data to workers (avoids pickling the full Molecule)
        task_data = (self._atom_data, self._pointer, self._count, self._size,
                     self._is_pbc, self._mol_box_arr, atom_type, distance)
        tasks = [(cube_chunk,) + task_data for cube_chunk in cube_np]

        # fork context inherits parent memory (copy-on-write) — avoids pickling
        # molecule data entirely. Use spawn on Windows (no fork support), during
        # testing (pytest adds OS threads not visible to threading.active_count()),
        # and in any explicitly multi-threaded process.
        import sys, threading
        _in_test = "_pytest" in sys.modules or "pytest" in sys.modules
        if sys.platform == "win32" or threading.active_count() > 1 or _in_test:
            _start = "spawn"
        else:
            _start = "fork"
        ctx = mp.get_context(_start)
        pool = ctx.Pool(processes=self._np)
        results = pool.starmap(_find_bond_task, tasks)
        pool.terminate()
        bond_list = sum(results, [])

        return bond_list


    ##################
    # Setter methods #
    ##################
    def set_pbc(self, pbc):
        """Turn the periodic boundary conditions on or off.

        Parameters
        ----------
        bond : bool
            True to turn on periodic boundary conditions
        """
        self._is_pbc = pbc


    ##################
    # Getter methods #
    ##################
    def get_origin(self):
        """Return the origin positions of the cubes.

        Returns
        -------
        origin : dictionary
            Dictionary of origin positions for each cube
        """
        return self._origin

    def get_pointer(self):
        """Return the list of atoms in each cube.

        Returns
        -------
        pointer : dictionary
            Pointer dictionary
        """
        return self._pointer

    def get_count(self):
        """Return the number of cubes in each dimension.

        Returns
        -------
        count : list
            Number of cubes
        """
        return self._count

    def get_size(self):
        """Return the cubes size.

        Returns
        -------
        size : integer
            Cube size
        """
        return self._size

    def get_mol(self):
        """Return the molecule.

        Returns
        -------
        mol : Molecule
            Molecule
        """
        return self._mol
