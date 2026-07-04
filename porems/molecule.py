################################################################################
# Molecule Class                                                               #
#                                                                              #
"""All necessary function for creating and editing molecules."""
################################################################################

import numpy as np
import pandas as pd

import porems.database as db
import porems.geometry as geometry
from porems.atom import Atom


class Molecule:
    """This class defines a molecule object, which is basically a list of Atom
    objects combined with a numpy-backed positions list. Each atom object carries
    the atom type, name and residue; positions are stored separately for
    vectorised bulk operations.

    Functions have been provided for editing, moving and transforming the
    atom objects, either as a collective or specific part of the molecule.

    Parameters
    ----------
    name : string, optional
        Molecule name
    short : string, optional
        Molecule short name
    inp : None, string, list, optional
        None for empty Molecule, string to read a molecule from a
        specified file link, or a list of Molecule objects to concatenate.

    Examples
    --------
    Following example generates a benzene molecule without hydrogen atoms

    .. code-block:: python

        import porems as pms

        mol = pms.Molecule("benzene", "BEN")
        mol.add("C", [0,0,0])
        mol.add("C", 0, r=0.1375, theta= 60)
        mol.add("C", 1, r=0.1375, theta=120)
        mol.add("C", 2, r=0.1375, theta=180)
        mol.add("C", 3, r=0.1375, theta=240)
        mol.add("C", 4, r=0.1375, theta=300)

    """

    def __init__(self, name="molecule", short="MOL", inp=None):
        # Initialize
        self._dim = 3

        self._name = name
        self._short = short

        self._box = []
        self._charge = 0
        self._masses = []
        self._mass = 0

        self._atom_list = []
        self._pos_list = []  # list of [x, y, z] — fast to append, vectorised on bulk ops

        # Check data input
        if inp is None:
            pass
        elif isinstance(inp, str):
            self._read(inp, inp.split(".")[-1].upper())
        elif isinstance(inp, list) and inp:
            if isinstance(inp[0], Molecule):
                self._concat(inp)
            # Legacy: list of Atom objects (positions are lost — use Molecule(inp=[mol]) instead)
            elif isinstance(inp[0], Atom):
                self._atom_list = list(inp)
                self._pos_list = [[0.0, 0.0, 0.0] for _ in inp]

    ##################
    # Representation #
    ##################
    def __repr__(self):
        """Create a pandas table of the molecule data.

        Returns
        -------
        repr : string
            Pandas data frame string of the molecule object
        """
        columns = ["Residue", "Name", "Type", "x", "y", "z"]
        data = []
        for i, atom in enumerate(self._atom_list):
            p = self._pos_list[i]
            data.append(
                [
                    atom.get_residue(),
                    atom.get_name(),
                    atom.get_atom_type(),
                    p[0],
                    p[1],
                    p[2],
                ]
            )
        return pd.DataFrame(data, columns=columns).to_string()

    ##############
    # Management #
    ##############
    def _read(self, file_path, file_type):
        """Read a molecule from a file. Currently only **GRO**, **PDB** and
        **MOL2** files are supported.

        Parameters
        ----------
        file_path : string
            Link to requested file
        file_type : string
            File extension name
        """
        if file_type not in ["GRO", "PDB", "MOL2"]:
            print("Unsupported filetype.")
            return

        atom_list = []
        pos_list = []

        with open(file_path, "r") as file_in:
            for line_idx, line in enumerate(file_in):
                line_val = line.split()
                is_add = False

                # Gro file
                if file_type == "GRO":
                    if line_idx > 0 and len(line_val) > 3:
                        residue = int(line[0:5]) - 1
                        pos = [float(line_val[i]) for i in range(3, 5 + 1)]
                        name = line_val[1]
                        atom_type = "".join([i for i in line_val[1] if not i.isdigit()])
                        is_add = True

                # Pdb file
                elif file_type == "PDB":
                    if line_val[0] in ["ATOM", "HETATM"]:
                        residue = int(line[22:26]) - 1
                        pos = [float(line_val[i]) / 10 for i in range(6, 8 + 1)]
                        name = line_val[11]
                        atom_type = line_val[11]
                        is_add = True

                # Mol2 file
                elif file_type == "MOL2":
                    if len(line_val) > 8:
                        residue = 0
                        pos = [float(line_val[i]) / 10 for i in range(2, 4 + 1)]
                        name = line_val[1]
                        atom_type = "".join([i for i in line_val[1] if not i.isdigit()])
                        is_add = True

                if is_add:
                    atom_list.append(Atom(atom_type, name, residue))
                    pos_list.append(pos)

        self._atom_list = atom_list
        self._pos_list = pos_list

    def _concat(self, mol_list):
        """Concatenate a molecule list into the current object.

        Parameters
        ----------
        mol_list : list
            List of Molecule objects to be concatenated
        """
        self._atom_list = sum([mol._atom_list for mol in mol_list], [])
        self._pos_list = sum([mol._pos_list for mol in mol_list], [])

    def _temp(self, atoms):
        """Create a temporary molecule copy of specified atom ids.

        Returns the molecule copy and the index list so that callers can write
        back modified positions via :meth:`_writeback_temp`.

        Parameters
        ----------
        atoms : list
            List of atom ids to include

        Returns
        -------
        mol : Molecule
            Molecule copy
        indices : list
            Atom ids in *self* corresponding to the copy
        """
        mol = Molecule.__new__(Molecule)
        mol._dim = self._dim
        mol._name = self._name
        mol._short = self._short
        mol._box = []
        mol._charge = 0
        mol._masses = []
        mol._mass = 0
        mol._atom_list = [self._atom_list[x] for x in atoms]
        mol._pos_list = [list(self._pos_list[x]) for x in atoms]
        return mol, list(atoms)

    def _writeback_temp(self, temp, indices):
        """Write positions from a temporary molecule back into *self*.

        Parameters
        ----------
        temp : Molecule
            Temporary molecule whose positions should be written back
        indices : list
            Atom ids in *self* that correspond to *temp*
        """
        for i, atom_id in enumerate(indices):
            self._pos_list[atom_id] = temp._pos_list[i]

    def append(self, mol):
        """Append a given molecule to the current object.

        Parameters
        ----------
        mol : Molecule
            Molecule object
        """
        self._atom_list += mol._atom_list
        self._pos_list += mol._pos_list

    def column_pos(self):
        """Create column list of atom positions.

        Returns
        -------
        column : list
            Columns of all atom positions — ``[[x0, x1, ...], [y0, y1, ...], [z0, z1, ...]]``
        """
        if not self._pos_list:
            return [[], [], []]
        return np.array(self._pos_list).T.tolist()

    ############
    # Geometry #
    ############
    def _vector(self, pos_a, pos_b):
        """Calculate the vector between two positions as defined in
        :func:`porems.geometry.vector` with the addition to define the inputs
        as atom indices.

        Parameters
        ----------
        pos_a : integer, list
            First position
        pos_b : integer, list
            Second position

        Returns
        -------
        vector : numpy.ndarray
            Bond vector
        """
        if isinstance(pos_a, int) and isinstance(pos_b, int):
            pos_a = self.pos(pos_a)
            pos_b = self.pos(pos_b)
        elif not (
            isinstance(pos_a, (list, np.ndarray))
            and isinstance(pos_b, (list, np.ndarray))
        ):
            print("Vector: Wrong input...")
            return None

        if not len(pos_a) == self._dim:
            print("Vector: Wrong dimensions...")
            return None

        return geometry.vector(pos_a, pos_b)

    def _box_size(self):
        """Calculate the box size of the current molecule.

        Returns
        -------
        box : list
            Box length in all dimensions
        """
        if not self._pos_list:
            return [0.001] * self._dim
        maxes = np.array(self._pos_list).max(axis=0)
        return [float(m) if float(m) > 0 else 0.001 for m in maxes]

    ##############
    # Properties #
    ##############
    def pos(self, atom):
        """Get the position of an atom.

        Parameters
        ----------
        atom : integer
            Atom id

        Returns
        -------
        pos : list
            Position vector of the specified atom
        """
        return self._pos_list[atom]

    def bond(self, inp_a, inp_b):
        """Return the bond vector of a specified bond. The two inputs can either
        be atom indices or positional vectors.

        Parameters
        ----------
        inp_a : integer, list
            Either an atom id or a position vector
        inp_b : integer, list
            Either an atom id or a position vector

        Returns
        -------
        bond : numpy.ndarray
            Bond vector

        Examples
        --------
        .. code-block:: python

            mol.bond(0, 1)
            mol.bond(*[0, 1])
            mol.bond([1, 0, 0], [0, 0, 0])
        """
        return self._vector(inp_a, inp_b)

    def centroid(self):
        """Calculate the geometrical centre of mass.

        Returns
        -------
        centroid : list
            Geometrical centre of mass
        """
        return np.mean(np.array(self._pos_list, dtype=float), axis=0).tolist()

    def com(self):
        """Calculate the centre of mass.

        Returns
        -------
        com : list
            Centre of mass
        """
        arr = np.array(self._pos_list, dtype=float)
        masses = np.array(self.get_masses(), dtype=float)
        return ((arr * masses[:, np.newaxis]).sum(axis=0) / masses.sum()).tolist()

    #################
    # Basic Editing #
    #################
    def translate(self, vec):
        """Translate all atom positions along a vector.

        Parameters
        ----------
        vec : list
            Translation vector
        """
        arr = np.array(self._pos_list, dtype=float)
        arr += np.asarray(vec, dtype=float)
        self._pos_list = arr.tolist()

    def rotate(self, axis, angle, is_deg=True):
        """Rotate all atom positions around an axis using
        :func:`porems.geometry.rotate`.

        Parameters
        ----------
        axis : integer, string, list
            Rotation axis
        angle : float
            Angle
        is_deg : bool, optional
            True if the input is in degrees
        """
        if not self._pos_list:
            return
        arr = np.array(self._pos_list, dtype=float)  # (N, 3)
        result = geometry.rotate(arr.T, axis, angle, is_deg)  # (3, N)
        if result is not None:
            self._pos_list = np.asarray(result, dtype=float).T.tolist()

    def move(self, atom, pos):
        """Move the molecule so that a given atom sits at a specified position.

        Parameters
        ----------
        atom : integer
            Atom id to use as the dragging point
        pos : list
            Target position
        """
        self.translate(self._vector(self.pos(atom), pos))

    def zero(self, pos=[0, 0, 0]):
        """Shift the molecule so that its minimum coordinate equals ``pos``.

        Parameters
        ----------
        pos : list, optional
            Zero-point coordinates

        Returns
        -------
        vec : list
            Translation vector used
        """
        if not self._pos_list:
            return [0.0] * self._dim
        arr = np.array(self._pos_list, dtype=float)
        mins = arr.min(axis=0)
        vec = [float(pos[i]) - float(mins[i]) for i in range(self._dim)]
        self._box = []
        self.translate(vec)
        return vec

    def put(self, atom, pos):
        """Move a single atom to a given position.

        Parameters
        ----------
        atom : integer
            Atom id
        pos : list
            New position vector
        """
        self._pos_list[atom] = list(pos)

    ####################
    # Advanced Editing #
    ####################
    def part_move(self, bond, atoms, length, vec=None):
        """Change the length of a specified bond by translating a subset of atoms.

        Parameters
        ----------
        bond : list
            Two atom ids defining the bond to adjust
        atoms : integer, list
            Atom ids to translate
        length : float
            New bond length
        vec : list, optional
            Override the translation direction
        """
        if isinstance(atoms, int):
            atoms = [atoms]
        temp, indices = self._temp(atoms)

        length = abs(length - geometry.length(self.bond(*bond)))

        if not vec:
            vec = self._vector(bond[0], bond[1])
        vec = [v * length for v in geometry.unit(vec)]

        temp.translate(vec)
        self._writeback_temp(temp, indices)

    def part_rotate(self, bond, atoms, angle, zero):
        """Rotate a subset of atoms around a bond axis.

        Parameters
        ----------
        bond : list
            Two atom ids defining the rotation axis
        atoms : integer, list
            Atom ids to rotate
        angle : float
            Rotation angle
        zero : integer
            Atom id used to set the coordinate origin before rotating
        """
        self.move(zero, [0, 0, 0])
        if isinstance(atoms, int):
            atoms = [atoms]
        temp, indices = self._temp(atoms)

        temp.rotate([self.pos(bond[0]), self.pos(bond[1])], angle)
        self._writeback_temp(temp, indices)

    def part_angle(self, bond_a, bond_b, atoms, angle, zero):
        """Change the angle between two bonds by rotating a subset of atoms.

        Parameters
        ----------
        bond_a : list
            First bond (two atom ids or a 3-vector)
        bond_b : list
            Second bond (two atom ids or a 3-vector)
        atoms : integer, list
            Atom ids to rotate
        angle : float
            Rotation angle
        zero : integer
            Atom id used to set the coordinate origin before rotating
        """
        self.move(zero, [0, 0, 0])
        if isinstance(atoms, int):
            atoms = [atoms]
        temp, indices = self._temp(atoms)

        if len(bond_a) == len(bond_b):
            if len(bond_a) == 2:
                vec = geometry.cross_product(
                    self._vector(*bond_a), self._vector(*bond_b)
                )
            elif len(bond_a) == self._dim:
                vec = geometry.cross_product(bond_a, bond_b)
            else:
                print("Part_Angle : Wrong bond input...")
                return
        else:
            print("Part_Angle : Wrong bond dimensions...")
            return

        temp.rotate(vec, angle)
        self._writeback_temp(temp, indices)

    #########
    # Atoms #
    #########
    def add(
        self,
        atom_type,
        pos,
        bond=None,
        r=0,
        theta=0,
        phi=0,
        is_deg=True,
        name="",
        residue=0,
    ):
        """Add a new atom in polar coordinates. The ``pos`` input is either
        an atom id that determines the bond-start, or a vector for a specific
        position.

        Parameters
        ----------
        atom_type : string
            Atom type
        pos : integer, list
            Position of the atom
        bond : list, optional
            Bond axis
        r : float, optional
            Bond length
        theta : float, optional
            Azimuthal angle
        phi : float, optional
            Polar angle
        is_deg : bool, optional
            True if angles are given in degrees
        name : string, optional
            Unique atom name
        residue : integer, optional
            Residue number

        Examples
        --------
        .. code-block:: python

            mol.add("C", [0, 0, 0])
            mol.add("C", 0, r=0.153, theta=-135)
            mol.add("C", 1, [0, 1], r=0.153, theta= 135)
        """
        pos = self.pos(pos) if isinstance(pos, int) else pos
        vec = self._vector(*bond) if bond else geometry.main_axis("z")

        phi += geometry.angle_polar(vec, is_deg)
        theta += geometry.angle_azi(vec, is_deg)

        phi *= np.pi / 180 if is_deg else 1
        theta *= np.pi / 180 if is_deg else 1

        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        self._atom_list.append(Atom(atom_type, name, residue))
        self._pos_list.append([float(pos[0] + x), float(pos[1] + y), float(pos[2] + z)])

    def delete(self, atoms):
        """Delete specified atoms from the molecule.

        Parameters
        ----------
        atoms : integer, list
            Atom id or list of ids to delete
        """
        atoms = [atoms] if isinstance(atoms, int) else atoms
        for atom in sorted(atoms, reverse=True):
            self._atom_list.pop(atom)
            self._pos_list.pop(atom)

    def overlap(self, error=0.005):
        """Search for overlapping atoms using a KD-tree (Chebyshev metric).

        Parameters
        ----------
        error : float, optional
            Maximum coordinate-wise distance to consider atoms overlapping

        Returns
        -------
        duplicates : dictionary
            Dictionary mapping each representative atom id to the list of
            atoms that overlap with it
        """
        from scipy.spatial import KDTree

        if not self._pos_list:
            return {}
        arr = np.array(self._pos_list)
        pairs = sorted(KDTree(arr).query_pairs(error, p=np.inf))

        marked_as_b = set()
        duplicates = {}
        for a, b in pairs:
            if a not in marked_as_b:
                if a not in duplicates:
                    duplicates[a] = []
                duplicates[a].append(b)
                marked_as_b.add(b)

        return duplicates

    def switch_atom_order(self, atom_a, atom_b):
        """Swap two atoms in the atom list and positions.

        Parameters
        ----------
        atom_a : integer
            First atom id
        atom_b : integer
            Second atom id
        """
        self._atom_list[atom_a], self._atom_list[atom_b] = (
            self._atom_list[atom_b],
            self._atom_list[atom_a],
        )
        self._pos_list[atom_a], self._pos_list[atom_b] = (
            self._pos_list[atom_b],
            self._pos_list[atom_a],
        )

    def set_atom_type(self, atom, atom_type):
        """Change the atom type of a specified atom.

        Parameters
        ----------
        atom : integer
            Atom id
        atom_type : string
            New atom type
        """
        self._atom_list[atom].set_atom_type(atom_type)

    def set_atom_name(self, atom, name):
        """Change the atom name of a specified atom.

        Parameters
        ----------
        atom : integer
            Atom id
        name : string
            New atom name
        """
        self._atom_list[atom].set_name(name)

    def set_atom_residue(self, atom, residue):
        """Change the residue index of a specified atom.

        Parameters
        ----------
        atom : integer
            Atom id
        residue : integer
            New residue index
        """
        self._atom_list[atom].set_residue(residue)

    def get_atom_type(self, atom):
        """Return the atom type of the given atom id.

        Parameters
        ----------
        atom : integer
            Atom id

        Returns
        -------
        atom_type : string
            Atom type
        """
        return self._atom_list[atom].get_atom_type()

    def get_atom_list(self):
        """Return the atom list.

        Returns
        -------
        atom_list : list
            List of :class:`porems.atom.Atom` objects
        """
        return self._atom_list

    ##################
    # Setter Methods #
    ##################
    def set_name(self, name):
        """Set the molecule name.

        Parameters
        ----------
        name : string
            Molecule name
        """
        self._name = name

    def set_short(self, short):
        """Set the molecule short name.

        Parameters
        ----------
        short : string
            Molecule short name
        """
        self._short = short

    def set_box(self, box):
        """Set the box size.

        Parameters
        ----------
        box : list
            Box size in all dimensions
        """
        self._box = box

    def set_charge(self, charge):
        """Set the total charge of the molecule.

        Parameters
        ----------
        charge : float
            Total molecule charge
        """
        self._charge = charge

    def set_masses(self, masses=None):
        """Set the molar masses of the atoms.

        Parameters
        ----------
        masses : list, optional
            List of molar masses in g/mol
        """
        self._masses = (
            masses
            if masses
            else [db.get_mass(atom.get_atom_type()) for atom in self._atom_list]
        )

    def set_mass(self, mass=0):
        """Set the molar mass of the molecule.

        Parameters
        ----------
        mass : float, optional
            Molar mass in g/mol
        """
        self._mass = mass if mass else sum(self.get_masses())

    ##################
    # Getter Methods #
    ##################
    def get_name(self):
        """Return the molecule name.

        Returns
        -------
        name : string
            Molecule name
        """
        return self._name

    def get_short(self):
        """Return the molecule short name.

        Returns
        -------
        short : string
            Molecule short name
        """
        return self._short

    def get_box(self):
        """Return the box size of the molecule.

        Returns
        -------
        box : list
            Box size in all dimensions
        """
        return self._box if self._box else self._box_size()

    def get_num(self):
        """Return the number of atoms.

        Returns
        -------
        num : integer
            Number of atoms
        """
        return len(self._atom_list)

    def get_charge(self):
        """Return the total charge of the molecule.

        Returns
        -------
        charge : float
            Total charge
        """
        return self._charge

    def get_masses(self):
        """Return a list of molar masses of the atoms.

        Returns
        -------
        masses : list
            Masses in g/mol
        """
        if not self._masses:
            self.set_masses()
        return self._masses

    def get_mass(self):
        """Return the molar mass of the molecule.

        Returns
        -------
        mass : float
            Molar mass in g/mol
        """
        if not self._mass:
            self.set_mass()
        return self._mass
