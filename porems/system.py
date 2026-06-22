################################################################################
# Basic Pore System Classes                                                    #
#                                                                              #
"""Here basic pore system constructions are defined."""
################################################################################


import os
import math
import yaml
import pandas as pd
import porems as pms
import numpy as np


class PoreKit():
    """Pore construction kit.
    """
    # Key used in sites_sl_shape / sites_shape for binding sites not assigned to any specific shape
    _UNASSIGNED_KEY = 20

    def __init__(self):
        # Initialize
        self._sort_list = ["OM", "SI"]
        self._res = 0
        self._hydro = {"in": [], "ex": 0}
        self._shapes = []
        self._yml = {}

    def structure(self, structure):
        """Add structure to the Pore builder.

        Parameters
        ----------
        structure : Molecule
            Crystal structure given as a PoreMS Molecule object
        """
        # Globalize crystal structure
        self._block = structure
        self._box = self._block.get_box()
        self._centroid = self._block.centroid()

    def build(self, bonds=[0.155-1e-2, 0.155+1e-2]):
        """Process provided structure to build connectivity matrix.

        Parameters
        ----------
        bonds : list, optional
            Minimal ans maximal bond length range for Si-O bonds
        """
        # Dice up block
        dice = pms.Dice(self._block, 0.4, True)
        self._matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], bonds))

        # Create pore object
        self._pore = pms.Pore(self._block, self._matrix)
        self._pore.set_name("pore")

    ##########
    # Shapes #
    ##########
    def exterior(self, res, hydro=0):
        """Optionally create and adjust exterior surface.

        Parameters
        ----------
        res : float
            Reservoir length in nm
        hydro : float, optional
            Hydroxilation degree for exterior surface in
            :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`
            leave zero for no adjustment
        """
        self._pore.exterior()
        self._hydro["ex"] = hydro
        self._res = res

    def add_shape(self, shape, section={"x": [], "y": [], "z": []}, hydro=0):
        """Add shape to pore system for drilling.

        Parameters
        ----------
        shape : list
            Shape type (pos 0) and shape (pos 1)
        section : dictionary, optional
            Range of shape from start x,y,z-length to end x,y,z-length,
            leave empty for whole range - mainly used to assign
            silanol groups to the shapes
        hydro : float, optional
            Hydroxylation degree for interior surface in
            :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`
        """
        # Check shape type
        if shape[0] not in ["CYLINDER", "SLIT", "SPHERE", "CONE"]:
            print("Wrong shape type...")
            return

        # Process user input
        coord_id = {"x": 0, "y": 1, "z": 2}
        ranges = {}
        for coord, sec in section.items():
            ranges[coord_id[coord]] = sec if sec else [0, self._box[coord_id[coord]]]
        self._hydro["in"].append(hydro)

        # Append to shape list
        shape.append(ranges)
        self._shapes.append(shape)

    def shape_cylinder(self, diam, length=0, centroid=[], central=[0, 0, 1]):
        """Add cylindrical shape

        Parameters
        ----------
        diam : float
            Cylinder diameter
        length : float, optional
            length of cylindrical shape leave zero for full length
        centroid : list, optional
            Cylinder centroid - leave zero for system centroid
        central : list, optional
            Central axis for cylinder

        Returns
        -------
        shape : list
            Shape type (pos 0) and shape (pos 1)
        """
        # Process user input
        centroid = centroid if centroid else self.centroid()
        length = length if length else self._box[2]

        # Define shape
        cylinder = pms.Cylinder({"centroid": centroid, "central": central, "length": length, "diameter": diam-0.5})  # Preparation precaution

        return ["CYLINDER", cylinder]

    def shape_slit(self, height, length=0, centroid=[], central=[0, 0, 1]):
        """Add slit shape

        Parameters
        ----------
        height : float
            Cylinder diameter
        length : float, optional
            length of cylindrical shape
        centroid : list, optional
            Cylinder centroid - leave zero for system centroid
        central : list, optional
            Central axis for cylinder

        Returns
        -------
        shape : list
            Shape type (pos 0) and shape (pos 1)
        """
        # Process user input
        centroid = centroid if centroid else self.centroid()
        length = length if length else self._box[2]

        # Define shape
        cuboid = pms.Cuboid({"centroid": centroid, "central": central, "length": length, "width": self._box[0], "height": height-0.5})  # Preparation precaution

        return ["SLIT", cuboid]

    def shape_sphere(self, diameter, centroid=[], central=[0, 0, 1]):
        """Add sphere shape

        Parameters
        ----------
        diameter : float
            Sphere diameter
        centroid : list, optional
            Sphere centroid - leave zero for system centroid
        central : list, optional
            Central axis for Sphere

        Returns
        -------
        shape : list
            Shape type (pos 0) and shape (pos 1)
        """
        # Process user input
        centroid = centroid if centroid else self.centroid()

        # Define shape
        sphere = pms.Sphere({"centroid": centroid, "central": central, "diameter": diameter})

        return ["SPHERE", sphere]

    def shape_cone(self, diam_1, diam_2, length=0, centroid=[], central=[0, 0, 1]):
        """Add cone shape

        Parameters
        ----------
        diam_1 : float
            Cone first diameter
        diam_2 : float
            Cone second diameter
        length : float, optional
            length of cylindrical shape leave zero for full length
        centroid : list, optional
            Cone centroid - leave zero for system centroid
        central : list, optional
            Central axis for cone

        Returns
        -------
        shape : list
            Shape type (pos 0) and shape (pos 1)
        """
        # Process user input
        centroid = centroid if centroid else self.centroid()
        length = length if length else self._box[2]

        # Define shape
        cone = pms.Cone({"centroid": centroid, "central": central, "length": length, "diameter_1": diam_1-0.5, "diameter_2": diam_2-0.5})  # Preparation precaution

        return ["CONE", cone]

    def prepare(self):
        """Prepare pore surface, add siloxane bridges, assign sites to sections,
        and assign a unique normal vector to each site by updating the
        original pore object site list.
        """
        # Carve out shape
        del_list = []
        for shape in self._shapes:
            del_list += [atom_id for atom_id in range(self._block.get_num()) if shape[1].is_in(self._block.pos(atom_id))]
        self._matrix.strip(del_list)

        # Prepare pore surface
        self._pore.prepare()

        # Determine sites
        self._pore.sites()
        site_list = self._pore.get_sites()

        # Check if sections are given
        self._sections = [shape[2] for shape in self._shapes]
        self._is_section = False
        for section in self._sections:
            if not section=={0: [0, self._box[0]], 1: [0, self._box[1]], 2: [0, self._box[2]]}:
                self._is_section = True

        # Define sites and save binding site si positions and allocate to interior section
        self._site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="ex"]
        self._si_pos_ex = [self._block.pos(site_key) for site_key in self._site_ex]

        self._site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]
        self._si_pos_in = [[] for x in self._sections]

        # Add normal vector to interior pore site list
        for site in self._site_in:
            # Get geometric position
            pos = self._block.pos(site)

            # Multiple shapes
            if len(self._shapes)>1:
                # Check distance to shape centroid
                lengths = []
                for i in range(len(self._shapes)):
                    if self._shapes[i][1].get_inp()["central"]==[0,0,1]:
                        inp_i = self._shapes[i][1].get_inp()
                        if "diameter_1" in inp_i:
                            dia = (inp_i["diameter_1"] + inp_i["diameter_2"]) / 2
                            lengths.append(pms.geom.length(pms.geom.vector(inp_i["centroid"][:2], pos[:2])) / (0.5*dia))
                        else:
                            lengths.append(pms.geom.length(pms.geom.vector(inp_i["centroid"][:2], pos[:2])) / (0.5*inp_i["diameter"]))
                    if self._shapes[i][1].get_inp()["central"]==[0,1,0]:
                        centroid = self._shapes[i][1].get_inp()["centroid"]
                        lengths.append(pms.geom.length(pms.geom.vector([centroid[0], centroid[2]], [pos[0],pos[2]]))/(0.5*self._shapes[i][1].get_inp()["diameter"]))
                    if self._shapes[i][1].get_inp()["central"]==[1,0,0]:
                        centroid = self._shapes[i][1].get_inp()["centroid"]
                        lengths.append(pms.geom.length(pms.geom.vector([centroid[1], centroid[2]], [pos[0],pos[2]]))/(0.5*self._shapes[i][1].get_inp()["diameter"]))
                    if self._shapes[i][1].get_inp()["central"]==[1,1,0]:
                        centroid = self._shapes[i][1].get_inp()["centroid"]
                        lengths.append(pms.geom.length(pms.geom.vector([centroid[-1]], [pos[-1]]))/(0.5*self._shapes[i][1].get_inp()["diameter"]))
                
                # Fin minimal length id                            
                min_len_id = lengths.index(min(lengths))
  
                # If sections are given
                if self._is_section:
                    for i in range(len(self._shapes)):
                        # Check if position within section
                        for i, section in enumerate(self._sections):
                            if self._shapes[i][1].get_inp()["central"]==[0,0,1] :
                                is_in = []
                                for dim in range(3):
                                    if dim == 2:
                                        #is_in.append(pos[dim]>=self._shapes[i][1].get_inp()["centroid"][dim]-self._shapes[i][1].get_inp()["length"]/2 and pos[dim]<=self._shapes[i][1].get_inp()["centroid"][dim]+self._shapes[i][1].get_inp()["length"]/2)
                                        is_in.append(self._shapes[i][1].is_in(pos))
                                    else:
                                        inp_i = self._shapes[i][1].get_inp()
                                        if "diameter_1" in inp_i:
                                            dia = (inp_i["diameter_1"] + inp_i["diameter_2"]) / 2
                                        else:
                                            dia = inp_i["diameter"]
                                        is_in.append(inp_i["centroid"][dim]-dia/2 <= pos[dim] <= inp_i["centroid"][dim]+dia/2)
                                if sum(is_in)==3:
                                    min_len_id = i
            else:
                min_len_id = 0

            # Add position to global list
            self._si_pos_in[min_len_id].append(pos)

            # Add normal vector to pore site list
            site_list[site]["normal"] = self._shapes[min_len_id][1].normal

        # Add normal vector to exterior pore site list
        for site in self._site_ex:
            site_list[site]["normal"] = self._normal_ex

        # Sanity check
        num_site_err = sum([1 for site in site_list if "normal" not in site_list[site].keys()])
        if num_site_err > 0:
            print("%i"%num_site_err+" sites were not assigned to shapes. Consider adjusting section intervals.")

        # Siloxane bridges
        if self._hydro["in"]: 
            self._siloxane("in")
        if self._hydro["ex"]:
            self._siloxane("ex")

        # Update site list after siloxan bridges
        self._site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="ex"]
        self._site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"]=="in"]

        # Calculate free binding sites of the several pores
        # Initialize dictonaries
        self.sites_sl_shape = {}

        # Loop over the shapes
        for i,shape in enumerate(self._shapes):
            self.sites_sl_shape[i] = []
            self.sites = []

        # Bucket for binding sites not assigned to any specific shape
        self.sites_sl_shape[self._UNASSIGNED_KEY] = []

        # Loop over the free binding sites
        for site in self._site_in:
            # Get position of the free site
            p = self._pore.get_block().pos(site)
            for i, shape in enumerate(self._shapes):
                # Set properties of the shape
                centroid = self._shapes[i][1].get_inp()["centroid"]
                if not shape[0]=="SPHERE":
                    length = self._shapes[i][1].get_inp()["length"]
                radi = pms.geometry.length(pms.geometry.vector([centroid[0], centroid[1], p[2]], p))
                
                # if shape is cyclinder
                if shape[0]=="CYLINDER":
                    if (centroid[2]-length/2)<p[2]<(centroid[2]+length/2):
                        if radi < ((self._shapes[i][1].get_inp()["diameter"]*1.5)/2):
                            self.sites_sl_shape[i].append(site)
                            self.sites.append(site)
                
                # if shape is cone
                elif shape[0]=="CONE":
                    if (centroid[2]-length/2)<p[2]<(centroid[2]+(length)/2):
                        if radi < ((self._shapes[i][1].get_inp()["diameter_1"]*1.5)/2):
                            self.sites_sl_shape[i].append(site)
                            self.sites.append(site)
                
                # if shape is slit
                elif shape[0]=="SLIT":
                    self.sites_sl_shape[i].append(site)
                    self.sites.append(site)

                elif shape[0]=="SPHERE":
                    if (centroid[2]-(length)/2)<p[2]<(centroid[2]+(length)/2):
                        if radi < ((self._shapes[i][1].get_inp()["diameter"]*1.05)/2):
                            self.sites_sl_shape[i].append(site)
                            self.sites.append(site)


            # If no match to one shape
            if site not in self.sites:
                self.sites_sl_shape[self._UNASSIGNED_KEY].append(site)
        # Drop the unassigned bucket if it is empty
        if self.sites_sl_shape[self._UNASSIGNED_KEY] == []:
            del self.sites_sl_shape[self._UNASSIGNED_KEY]

               
        # Count the numbers of attached siloxane
        self._pore.sites_sl_shape = self.sites_sl_shape
        for i in self.sites_sl_shape:
            self._pore.sites_attach_mol[i]["SL"] = 0
            self._pore.sites_attach_mol[i]["SLG"] = 0
            for si in self.sites_sl_shape[i]:
                if (len(self._pore._sites[si]["o"])==1 and self._pore._sites[si]["type"]=="in"):
                    self._pore.sites_attach_mol[i]["SL"] +=1
                elif (len(self._pore._sites[si]["o"])==2 and self._pore._sites[si]["type"]=="in"):
                    self._pore.sites_attach_mol[i]["SLG"] +=1
        
        # Objectify grid
        non_grid = self._matrix.bound(1)+list(site_list.keys())
        bonded = self._matrix.bound(0, "gt")
        grid_atoms = [atom for atom in bonded if not atom in non_grid]
        self._pore.objectify(grid_atoms)

    ##############
    # Attachment #
    ##############
    def _normal_ex(self, pos):
        """Normal function for the exterior surface

        Parameters
        ----------
        pos : list
            Position on the surface

        Returns
        -------
        normal : list
            Vector perpendicular to surface of the given position
        """
        return [0, 0, -1] if pos[2] < self.centroid()[2] else [0, 0, 1]

    def _siloxane(self, site_type, slx_dist=[0.507-1e-2, 0.507+1e-2]):
        """Attach siloxane bridges using function
        :func:`porems.pore.Pore.siloxane`.

        Parameters
        ----------
        site_type : string
            Site type - interior **in**, exterior **ex**
        slx_dist : list, optional
            Silicon atom distance to search for partners in proximity
        """
        # Initialize
        site_list = self._pore.get_sites()
        sites = self._site_in if site_type=="in" else self._site_ex
        hydro = self._hydro["in"] if site_type=="in" else self._hydro["ex"]

        # Find free binding sites of Si atoms in the structure
        if site_type=="in":
            # Initialize dictonaries
            self.sites_shape = {}                   #Index of the binding sites
            self._pore.sites_attach_mol = {}        #number of attach molecules
            # Loop over the shapes
            for i,shapes in enumerate(self._shapes):
                self.sites_shape[i] = []
                self._pore.sites_attach_mol[i] = {}
            # Bucket for binding sites not assigned to any specific shape
            self._pore.sites_attach_mol[self._UNASSIGNED_KEY] = {}
            self.sites_shape[self._UNASSIGNED_KEY] = []

            # Loop over the free binding sites
            self.sites = []
            for site in sites:
                # Get position of the free site
                p = self._pore.get_block().pos(site)
                for i, shape in enumerate(self._shapes):
                    # Set properties of the shape
                    centroid = self._shapes[i][1].get_inp()["centroid"]
                    if not shape[0]=="SPHERE":
                        length = self._shapes[i][1].get_inp()["length"]
                    
                    radi = pms.geom.length(pms.geom.vector([centroid[0], centroid[1], p[2]], p))
                    
                    # if shape is cyclinder
                    if site not in self.sites:
                        if shape[0]=="CYLINDER":
                            if (centroid[2]-(length)/2)<p[2]<(centroid[2]+(length)/2):
                                if radi < ((self._shapes[i][1].get_inp()["diameter"]*1.5)/2):
                                    self.sites_shape[i].append(site)
                                    self.sites.append(site)
                        # if shape is cone
                        elif shape[0]=="CONE":
                            if (centroid[2]-(length)/2)<p[2]<(centroid[2]+(length)/2):
                                if radi < ((self._shapes[i][1].get_inp()["diameter_1"]*1.5)/2):
                                    self.sites_shape[i].append(site)
                                    self.sites.append(site)
                        # if shape is slit
                        elif shape[0]=="SLIT":
                            self.sites_shape[i].append(site)
                            self.sites.append(site)
                        elif shape[0]=="SPHERE":
                            if (centroid[2]-(length)/2)<p[2]<(centroid[2]+(length)/2):
                                if radi < ((self._shapes[i][1].get_inp()["diameter"]*1.05)/2):
                                    self.sites_shape[i].append(site)
                                    self.sites.append(site)

                # If no match to one shape
                if site not in self.sites:
                    self.sites_shape[self._UNASSIGNED_KEY].append(site)

            # Drop the unassigned bucket if it is empty
            if self.sites_shape[self._UNASSIGNED_KEY] == []:
               del self.sites_shape[self._UNASSIGNED_KEY]

            # If there are unassigned binding sites print a warning
            else:
                n = len(self.sites_shape[self._UNASSIGNED_KEY])
                print("Warning: %i interior sites could not be assigned to a shape. They will be filled with silanol/siloxane bridges." % n)

            # Amount - Connect two oxygen to one siloxane
            amount = []
            for hydro,surface,sites in zip(hydro,self.surface(is_sum=False)[site_type],self.sites_shape.values()):
                oh = len(sum([site_list[site]["o"] for site in sites], []))
                oh_goal = pms.utils.mumol_m2_to_mols(hydro,surface)
                amount.append(round((oh-oh_goal)/2))

            # Fill siloxane
            for amount, sites in zip(amount, self.sites_shape.items()):
                if amount > 0:
                    # Run attachment
                    mols = self._pore.siloxane(sites[1], amount, slx_dist=slx_dist, site_type=site_type)
                    for mol in mols:
                        # Count the numbers of attached siloxane
                        if mol.get_short() in self._pore.sites_attach_mol[sites[0]]:
                            self._pore.sites_attach_mol[sites[0]][mol.get_short()] +=1
                        else:
                            self._pore.sites_attach_mol[sites[0]][mol.get_short()] = 0
                            self._pore.sites_attach_mol[sites[0]][mol.get_short()] += 1

                    # Add to sorting list
                    for mol in mols:
                        if not mol.get_short() in self._sort_list:
                            self._sort_list.append(mol.get_short())
        else:
             # Amount - Connect two oxygen to one siloxane
            oh = len(sum([site_list[site]["o"] for site in sites], []))
            oh_goal = pms.utils.mumol_m2_to_mols(hydro, self.surface()[site_type])
            amount = round((oh-oh_goal)/2)

            # Fill siloxane
            if amount > 0:
                # Run attachment
                mols = self._pore.siloxane(sites, amount, slx_dist=slx_dist, site_type=site_type)
                
                # Add to sorting list
                for mol in mols:
                    if not mol.get_short() in self._sort_list:
                        self._sort_list.append(mol.get_short())
                        

    def attach(self, mol, mount, axis, amount, site_type="in", inp="num", shape="all", pos_list=[], scale=1, trials=1000, is_proxi=True, is_rotate=False, is_g=True):
        """Attach molecule on the surface.

        Parameters
        ----------
        mol : Molecule
            Molecule object to attach
        mount : integer
            Atom id of the molecule that is placed on the surface silicon atom
        axis : list
            List of two atom ids of the molecule that define the molecule axis
        amount : int
            Number of molecules to attach
        site_type : string, optional
            Use **in** for the interior surface and **ex** for the exterior
        inp : string, optional
            Input type: **num** - Number of molecules,
            **molar** - :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`,
            **percent** - :math:`\\%` of OH groups
        shape : string, optional
            Optional is "all" this means every shape will be functionalize.
            Otherwise specific the shape if you want to functionalize.
        pos_list : list, optional
            List of positions (Cartesian) to find nearest available binding site for
        scale : float, optional
            Circumference scaling around the molecule position
        trials : integer, optional
            Number of trials picking a random site
        is_proxi : bool, optional
            True to fill binding sites in proximity of filled binding site
        is_rotate : bool, optional
            True to randomly rotate molecule around own axis
        is_g : bool, optional
            Force to add molecules only on single binding sites
        """
        # Process input
        if site_type not in ["in", "ex"]:
            print("Pore: Wrong site_type input...")
            return

        if inp not in ["num", "molar", "percent"]:
            print("Pore: Wrong inp type...")
            return

        # Temporarily remove the unassigned-sites bucket so attach() only sees shape-specific sites
        if shape == "all":
            _saved_unassigned = self._pore.sites_sl_shape.pop(self._UNASSIGNED_KEY, None)

        # Amount of SL molecules
        # Input molar
        if inp=="molar":
            if site_type=="in" and shape!="all":
                amount = int(pms.utils.mumol_m2_to_mols(amount, self.surface(is_sum=False)[site_type][int(shape[-1])]))
            elif site_type=="in" and shape=="all":
                amount_list = {}
                for sites_shape_idx in self._pore.sites_sl_shape:
                    amount_list[sites_shape_idx] = int(pms.utils.mumol_m2_to_mols(amount, self.surface(is_sum=False)[site_type][sites_shape_idx]))
            elif site_type=="ex":
                 amount = int(pms.utils.mumol_m2_to_mols(amount, self.surface()[site_type]))
        
        # Input percent         
        elif inp=="percent":
            if site_type=="in" and shape!="all":
                sites = self._pore.sites_sl_shape[int(shape[-1])]
                num_oh = len(sites)
                num_oh += sum([1 for x in self._pore.get_sites().values() if len(x["o"])==2 and x["type"]==site_type])
                amount = int(amount/100*num_oh)
            elif site_type=="in" and shape=="all":
                amount_list = {}
                for sites_shape_idx in self._pore.sites_sl_shape:
                    sites = self._pore.sites_sl_shape[sites_shape_idx]
                    num_oh = len(sites)
                    num_oh += sum([1 for x in self._pore.get_sites().values() if len(x["o"])==2 and x["type"]==site_type])
                    amount_list[sites_shape_idx] = int(amount/100*num_oh)
            elif site_type=="ex":
                num_oh = len(sites)
                num_oh += sum([1 for x in self._pore.get_sites().values() if len(x["o"])==2 and x["type"]==site_type])
                amount = int(amount/100*num_oh)
        # Input number of molecules
        else:
            if site_type=="in" and shape!="all":
                amount = amount
            elif site_type=="in" and shape=="all":
                amount_list = {}
                for sites_shape_idx in self._pore.sites_sl_shape:
                    amount_list[sites_shape_idx] = int(amount/len(self._pore.sites_sl_shape))
            elif site_type=="ex":
                amount = amount
            

        # Check number of given positions
        if pos_list and not len(pos_list)==amount:
            print("Pore: Number of given positions does not match number of groups to attach...")
            return

        # Run attachment
        # Attachment on exterior
        if site_type=="ex":
            sites = self._site_ex
            mols = self._pore.attach(mol, mount, axis, sites, amount, scale, trials, pos_list=pos_list, site_type=site_type, is_proxi=is_proxi, is_random=True, is_rotate=is_rotate, is_g=is_g)
        # Attachment in interior (for a specific shape)
        elif site_type=="in" and shape!="all":
            sites = self._pore.sites_sl_shape[int(shape[-1])]
            mols = self._pore.attach(mol, mount, axis, sites, amount, scale, trials, pos_list=pos_list, site_type=site_type, is_proxi=is_proxi, is_random=True, is_rotate=is_rotate, is_g=is_g)
            # Remove Si sites which are no longer free
            for i in self._pore.sites_sl_shape:
                for si in self._pore.sites_sl_shape[i]:
                    state = self._pore._sites[si]["state"]
                    if state == False:
                        self._pore.sites_sl_shape[i].remove(si)
            # Count type and number of attach molecules
            for mol in mols:
                if mol.get_short() not in ["SL","SLG"]:
                    if mol.get_short() in self._pore.sites_attach_mol[int(shape[-1])]:
                        self._pore.sites_attach_mol[int(shape[-1])][mol.get_short()] +=1
                    else:                
                        self._pore.sites_attach_mol[int(shape[-1])][mol.get_short()] = 0
                        self._pore.sites_attach_mol[int(shape[-1])][mol.get_short()] += 1
        # Attachment in interior (in all shapes)    
        elif shape=="all" and site_type=="in":
            for sites_shape_idx in self._pore.sites_sl_shape:
                sites = self._pore.sites_sl_shape[sites_shape_idx]
                mols = self._pore.attach(mol, mount, axis, sites, amount_list[sites_shape_idx], scale, trials, pos_list=pos_list, site_type=site_type, is_proxi=is_proxi, is_random=True, is_rotate=is_rotate, is_g=is_g)
                # Remove Si sites which are no longer free
                for i in self._pore.sites_sl_shape:
                    for si in self._pore.sites_sl_shape[i]:
                        state = self._pore._sites[si]["state"]
                        if state == False:
                            self._pore.sites_sl_shape[i].remove(si)
                # Count type and number of attach molecules
                for mol in mols:
                    if mol.get_short() not in ["SL","SLG"]:
                        if mol.get_short() in self._pore.sites_attach_mol[i]:
                            self._pore.sites_attach_mol[i][mol.get_short()] +=1
                        else:
                            self._pore.sites_attach_mol[i][mol.get_short()] = 0
                            self._pore.sites_attach_mol[i][mol.get_short()] += 1
        
        # Add to sorting list
        for mol in mols:
            if not mol.get_short() in self._sort_list:
                self._sort_list.append(mol.get_short())

        # Restore unassigned-sites bucket
        if shape == "all" and _saved_unassigned is not None:
            self._pore.sites_sl_shape[self._UNASSIGNED_KEY] = _saved_unassigned

    def attach_special(self, mol, mount, axis, amount, scale=1, symmetry="point", is_proxi=True, is_rotate=False):
        """Attach molecules at geometrically symmetric positions along z-axis.

        Parameters
        ----------
        mol : Molecule
            Molecule to attach
        mount : integer
            Atom id used as mounting point
        axis : list
            Two atom ids defining the molecule axis
        amount : int
            Number of molecules to attach
        scale : float, optional
            Circumference scaling factor
        symmetry : string, optional
            ``"point"`` — alternating sides; ``"mirror"`` — same side
        is_proxi : bool, optional
            Fill proximity sites with silanol if True
        is_rotate : bool, optional
            Randomly rotate molecule around own axis if True
        """
        if symmetry not in ["point", "mirror"]:
            print("attach_special: unsupported symmetry type (use 'point' or 'mirror')")
            return

        dist = self._box[2] / amount if amount > 0 else 0
        start = dist / 2
        diam = self.diameter()
        radius = diam[0] / 2 if isinstance(diam, list) else diam / 2

        pos_list = []
        for i in range(amount):
            coeff = -1 if (symmetry == "point" and i % 2 == 0) else 1
            pos_list.append([
                self._centroid[0] + coeff * radius,
                self._centroid[1],
                start + dist * i,
            ])

        mols = self._pore.attach(
            mol, mount, axis, self._site_in, len(pos_list),
            scale=scale, pos_list=pos_list,
            site_type="in", is_proxi=is_proxi, is_random=False, is_rotate=is_rotate,
        )
        for m in mols:
            if m.get_short() not in self._sort_list:
                self._sort_list.append(m.get_short())

    ################
    # Finalization #
    ################
    def finalize(self):
        """Finalize pore system."""

        # Fill silanol on the exterior surface
        mols_ex = self._pore.fill_sites(self._site_ex, "ex") if self._site_ex else []
        for mol in mols_ex:
            if not mol.get_short() in self._sort_list:
                self._sort_list.append(mol.get_short())

        # Fill silanol on the interior surface
        for i,sites in self._pore.sites_sl_shape.items():
            mols_in = self._pore.fill_sites(sites, "in") if self._site_in else []
            for mol in mols_in:
                if not mol.get_short() in self._sort_list:
                    self._sort_list.append(mol.get_short())

        # Create reservoir
        self._pore.reservoir(self._res)

    def store(self, link="./", sort_list=[]):
        """Store pore system and all necessary files for simulation at given
        link.

        Parameters
        ----------
        link : string, optional
            Folder link for output
        """
        # Process input
        link = link if link[-1] == "/" else link+"/"

        # Set sort list
        sort_list = sort_list if sort_list else self._sort_list

        # Create store object
        store = pms.Store(self._pore, link, sort_list)

        # Save files
        store.gro(use_atom_names=True)
        store.obj()
        store.top()
        store.grid()
        pms.utils.save(self, link+self._pore.get_name()+"_system.obj")
        self.yml(link)

    def yml(self, link="./"):
        """Save yaml file with properties necessary for analysis.

        Parameters
        ----------
        link : string, optional
            Folder link for output
        """
        # Process input
        link = link if link[-1] == "/" else link+"/"

        # Fill system properties
        self._yml["system"] = {}
        self._yml["system"]["dimensions"] = self.box()
        self._yml["system"]["centroid"] = self.centroid()
        self._yml["system"]["reservoir"] = self.reservoir()
        self._yml["system"]["volume"] = self.volume()
        self._yml["system"]["surface"] = self.surface()

        # Calculate properties
        diameter = self.diameter()
        roughness = self.roughness()
        volume = self.volume(is_sum=False)
        surface = self.surface(is_sum=False)

        # Fill properties for each shape
        for i, shape in enumerate(self._shapes):
            shape_id = "shape_"+"%02i"%i
            self._yml[shape_id] = {}
            self._yml[shape_id]["shape"] = shape[0]
            self._yml[shape_id]["parameter"] = shape[1].get_inp().copy()
            if shape[0]=="CYLINDER":
                self._yml[shape_id]["parameter"]["diameter"] += 0.5
            elif shape[0]=="SLIT":
                self._yml[shape_id]["parameter"]["height"] += 0.5
            self._yml[shape_id]["diameter"] = diameter[i]
            self._yml[shape_id]["roughness"] = roughness["in"][i]
            self._yml[shape_id]["volume"] = volume[i]
            self._yml[shape_id]["surface"] = surface["in"][i]

        # Export
        yaml.Dumper.ignore_aliases = lambda *args : True
        with open(link+self._pore.get_name()+".yml", "w") as file_out:
            file_out.write(yaml.dump(self._yml, default_flow_style=False))

    ############
    # Analysis #
    ############
    def _radii_per_shape(self):
        """Compute the radial distances of binding-site Si atoms from their
        shape's central axis, grouped by shape.

        Returns
        -------
        radii : list
            One list of floats per shape in ``self._shapes``
        """
        radii = []
        pos_new = [0, 0, 0]

        for i, shape in enumerate(self._shapes):
            index_si = self.sites_shape.get(i, self._si_pos_in[i])
            centroid = self._shapes[i][1].get_inp()["centroid"]
            central  = self._shapes[i][1].get_inp()["central"]

            if shape[0] != "SPHERE":
                length = self._shapes[i][1].get_inp()["length"]
                z_min = centroid[2] - length/2 + 0.1
                z_max = centroid[2] + length/2 - 0.1
                centroid_new = [0, 0, 0]
                centroid_new[0] = centroid[0]*np.cos(-np.pi/4) - centroid[1]*np.sin(-np.pi/4)
                centroid_new[1] = centroid[0]*np.sin(-np.pi/4) + centroid[1]*np.cos(-np.pi/4)
                x_min = centroid_new[0] - 0.2
                x_max = centroid_new[0] + 0.2
            else:
                length = self._shapes[i][1].get_inp()["diameter"]
                z_min = centroid[2] - length/2
                z_max = centroid[2] + length/2

            radii_temp = []
            for index in index_si:
                pos = self._pore.get_block().pos(index) if isinstance(index, int) else index

                if shape[0] in ("CYLINDER", "CONE"):
                    if z_min < pos[2] < z_max and central == [0, 0, 1]:
                        radii_temp.append(pms.geom.length(pms.geom.vector([centroid[0], centroid[1], pos[2]], pos)))
                    elif central == [1, 1, 0] and shape[0] == "CYLINDER":
                        pos_new[0] = pos[0]*np.cos(-np.pi/4) - pos[1]*np.sin(-np.pi/4)
                        pos_new[1] = pos[0]*np.sin(-np.pi/4) + pos[1]*np.cos(-np.pi/4)
                        if x_min < pos_new[0] < x_max:
                            r = pms.geom.length(pms.geom.vector([pos_new[0], centroid_new[1], centroid_new[2]], pos_new))
                            diameter_inp = self._shapes[i][1].get_inp()["diameter"] + 0.5
                            if (diameter_inp/2)*0.9 < r < (diameter_inp/2)*1.1:
                                radii_temp.append(r)
                elif shape[0] == "SLIT":
                    radii_temp.append(pms.geom.length(pms.geom.vector([pos[0], centroid[1], pos[2]], pos)))
                elif shape[0] == "SPHERE":
                    if z_min < pos[2] < z_max and central == [0, 0, 1]:
                        radii_temp.append(pms.geom.length(pms.geom.vector(centroid, pos)))
            radii.append(radii_temp)
        return radii

    def diameter(self):
        """Calculate true diameter after drilling and preparation.

        Returns
        -------
        diameter : list
            List of shape diameters after preparation
        """
        radii = self._radii_per_shape()
        r_bar = [sum(r)/len(r) if r else 0 for r in radii]
        return [2*r for r in r_bar]

    def roughness(self):
        """Calculate surface roughness (RMS) for interior and exterior surfaces.

        Returns
        -------
        roughness : dict
            ``{"in": [rq_shape_0, ...], "ex": rq_exterior}``
        """
        # Interior
        radii_in = self._radii_per_shape()
        r_bar_in = [sum(r)/len(r) if r else 0 for r in radii_in]
        r_q_in = [
            math.sqrt(sum((r_i - r_bar)**2 for r_i in r) / len(r)) if r else 0
            for r, r_bar in zip(radii_in, r_bar_in)
        ]

        # Exterior
        size = 0
        if self._res and self._si_pos_ex:
            temp_mol = pms.Molecule()
            for pos in self._si_pos_ex:
                temp_mol.add("Si", pos)
            temp_mol.zero()
            size = temp_mol.get_box()[2]
        r_ex = [pos[2] if size == 0 or pos[2] < size/2 else abs(pos[2]-size) for pos in self._si_pos_ex]
        r_bar_ex = sum(r_ex)/len(r_ex) if r_ex else 0
        r_q_ex = math.sqrt(sum((r_i - r_bar_ex)**2 for r_i in r_ex)/len(r_ex)) if r_ex else 0

        return {"in": r_q_in, "ex": r_q_ex}

    def volume(self, is_sum=True):
        """Calculate pore volume. This is done by defining a new shape object
        with system sizes after pore preparation and using the volume functions.

        Notes
        -----
        Sphere volume is fully calculated.

        Parameters
        ----------
        is_sum : bool, optional
            True to return sum of all sections

        Returns
        -------
        volume : float, list
            Total pore volume if is_sum is True. Otherwise a list of volumes
        """
        # Get diameters
        diam = self.diameter()
        diam = diam if isinstance(diam, list) else [diam]

        # Calculate volumes
        volume = []
        for i, shape in enumerate(self._shapes):
            centroid = self._shapes[i][1].get_inp()["centroid"]
            if not shape[0]=="SPHERE":
                length = self._shapes[i][1].get_inp()["length"]
            if shape[0]=="CYLINDER":
                volume.append(pms.Cylinder({"centroid": centroid, "central": [0, 0, 1], "length": length, "diameter": diam[i]}).volume())
            elif shape[0]=="SLIT":
                volume.append(pms.Cuboid({"centroid": centroid, "central": [0, 0, 1], "length": length, "width": self._box[0], "height": diam[i]}).volume())
            elif shape[0]=="SPHERE":
                volume.append(pms.Sphere({"centroid": centroid, "central": [0, 0, 1], "diameter": diam[i]}).volume())
            if shape[0]=="CONE":
                diam_cone =  [self._shapes[i][1].get_inp()["diameter_1"],self._shapes[i][1].get_inp()["diameter_2"]]
                volume.append(pms.Cone({"centroid": centroid, "central": [0, 0, 1], "length": length, "diameter_1": diam_cone[0], "diameter_2": diam_cone[1]}).volume())

        return sum(volume) if is_sum else volume

    def surface(self, is_sum=True):
        """Calculate pore surface and exterior surface. This is done by defining
        a new shape object with system sizes after pore preparation and using
        the surface functions.

        Notes
        -----
        Sphere surface is fully calculated and outer surface assumes
        that spheres intersect at the center.

        Parameters
        ----------
        is_sum : bool, optional
            True to return sum of all sections

        Returns
        -------
        surface : dictionary
            Pore surface of interior and exterior
        """
        # Get diameters
        diam = self.diameter()
        diam = diam if isinstance(diam, list) else [diam]

        # Interior surface
        surf_in = []
        for i, shape in enumerate(self._shapes):
            centroid = self._shapes[i][1].get_inp()["centroid"]
            if not shape[0]=="SPHERE":
                length = self._shapes[i][1].get_inp()["length"]
            if shape[0]=="CYLINDER":
                surf_in.append(pms.Cylinder({"centroid": centroid, "central": [0, 0, 1], "length": length, "diameter": diam[i]}).surface())
            elif shape[0]=="SLIT":
                surf_in.append(pms.Cuboid({"centroid": centroid, "central": [0, 0, 1], "length": length, "width": self._box[0], "height": diam[i]}).surface()/2)
            elif shape[0]=="SPHERE":
                surf_in.append(pms.Sphere({"centroid": centroid, "central": [0, 0, 1], "diameter": diam[i]}).surface())
            if shape[0]=="CONE":
                diam_cone =  [self._shapes[i][1].get_inp()["diameter_1"],self._shapes[i][1].get_inp()["diameter_2"]]
                surf_in.append(pms.Cone({"centroid": centroid, "central": [0, 0, 1],"length": length, "diameter_1": diam_cone[0], "diameter_2": diam_cone[1]}).surface())

        # Exterior surface
        sections_ex = [0, 0]
        for i, section in enumerate(self._sections):
            if section[2][0]-0<=1e-2:
                sections_ex[0] = i
            if section[2][1]-self._box[2]<=1e-2:
                sections_ex[1] = i

        surf_ex = []
        for section in sections_ex:
            if self._shapes[section][0]=="CYLINDER":
                surf_ex.append(self._box[0]*self._box[1]-math.pi*(diam[section]/2)**2)
            elif self._shapes[section][0]=="SLIT":
                surf_ex.append(self._box[0]*(self._box[1]-diam[section]))
            elif self._shapes[section][0]=="SPHERE":
                surf_ex.append(self._box[0]*self._box[1]-math.pi*(diam[section]/2)**2)
            if self._shapes[section][0]=="CONE":
                surf_ex.append(0)

        return {"in": sum(surf_in) if is_sum else surf_in,"ex": sum(surf_ex) if is_sum else surf_ex}

    def allocation(self):
        """Calculate molecule allocation on the surface. Using interior and
        exterior surfaces, the allocation rates can be determined by counting
        the number of used molecules on the surfaces.

        Using the conversion function :func:`porems.utils.mols_to_mumol_m2`, the
        number of molecules is converted to concentration in
        :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`.

        Returns
        -------
        alloc : dictionary
            Dictionary containing the surface allocation of all molecules in
            number of molecules, :math:`\\frac{\\text{mols}}{\\text{nm}^2}` and
            :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`
        """
        # Get surfaces
        surf = self.surface()
        site_dict = self._pore.get_site_dict()

        # Calculate allocation for all molecules
        alloc = {}
        for mol in sorted(self._sort_list):
            for site_type in ["in", "ex"]:
                if mol in site_dict[site_type]:
                    if not mol in alloc:
                        alloc[mol] = {"in": [0, 0, 0], "ex": [0, 0, 0]}
                    # Number of molecules
                    alloc[mol][site_type][0] = len(site_dict[site_type][mol])

                    # Molecules per nano meter
                    alloc[mol][site_type][1] = len(site_dict[site_type][mol])/surf[site_type] if surf[site_type]>0 else 0

                    # Micromolar per meter
                    alloc[mol][site_type][2] = pms.utils.mols_to_mumol_m2(len(site_dict[site_type][mol]), surf[site_type]) if surf[site_type]>0 else 0

        # OH allocation
        alloc["OH"] = {"in": [0, 0, 0], "ex": [0, 0, 0]}
        for site_type in ["in", "ex"]:
            num_oh = len(sum([x["o"] for x in self._pore.get_sites().values() if x["type"]==site_type], []))
            for mol in site_dict[site_type].keys():
                num_oh -= len(site_dict[site_type][mol]) if mol not in ["SL", "SLG", "SLX"] else 0

            # num_oh = num_oh-num_in_ex if site_type=="ex" else num_oh+num_in_ex

            alloc["OH"][site_type][0] = num_oh
            alloc["OH"][site_type][1] = num_oh/surf[site_type] if surf[site_type]>0 else 0
            alloc["OH"][site_type][2] = pms.utils.mols_to_mumol_m2(num_oh, surf[site_type]) if surf[site_type]>0 else 0

        # Hydroxylation - Total number of binding sites
        alloc["Hydro"] = {"in": [0, 0, 0], "ex": [0, 0, 0]}
        for site_type in ["in", "ex"]:
            num_tot = len(sum([x["o"] for x in self._pore.get_sites().values() if x["type"]==site_type], []))

            # num_tot = num_tot-num_in_ex if site_type=="ex" else num_tot+num_in_ex

            alloc["Hydro"][site_type][0] = num_tot
            alloc["Hydro"][site_type][1] = num_tot/surf[site_type] if surf[site_type]>0 else 0
            alloc["Hydro"][site_type][2] = pms.utils.mols_to_mumol_m2(num_tot, surf[site_type]) if surf[site_type]>0 else 0

        return alloc

    def reservoir(self):
        """Return the reservoir length.

        Returns
        -------
        res : float
            Reservoir length
        """
        return self._res

    def box(self):
        """Return the box size of the pore block.

        Returns
        -------
        box : list
            Box size in all dimensions
        """
        return self._pore.get_box()

    def centroid(self):
        """Return pore centroid.

        Returns
        -------
        centroid : list
            Centroid of the pore
        """
        return self._centroid

    def shape(self):
        """Return the pore shape for analysis using PoreAna.

        Returns
        -------
        pore_shape : dictionary
            Pore shape
        """
        return self._shapes

    #########
    # Table #
    #########
    def table(self, decimals=3):
        """Create properties as pandas table for easy viewing.

        Parameters
        ----------
        decimals : integer, optional
            Number of decimals to be rounded to

        Returns
        -------
        tables : DataFrame
            Pandas table of all properties
        """
        # Initialize
        form = "%."+str(decimals)+"f"

        # Get allocation data
        allocation = self.allocation()

        # Get properties
        surf = self.surface()
        roughness = self.roughness()

        # Save data
        data = {"Interior": {}, "Exterior": {}}

        data["Interior"]["Silica block xyz-dimensions (nm)"] = " "
        data["Exterior"]["Silica block xyz-dimensions (nm)"] = "["+form%self.box()[0]+", "+form%self.box()[1]+", "+form%(self.box()[2]-2*self.reservoir())+"]"
        data["Interior"]["Simulation box xyz-dimensions (nm)"] = " "
        data["Exterior"]["Simulation box xyz-dimensions (nm)"] = "["+form%self.box()[0]+", "+form%self.box()[1]+", "+form%self.box()[2]+"]"
        data["Interior"]["Surface roughness (nm)"] = [form%val for val in roughness["in"]] if "in" in roughness else form%0
        data["Exterior"]["Surface roughness (nm)"] = form%roughness["ex"] if "ex" in roughness else form%0
        for i, val in enumerate(self.diameter()):
            data["Interior"]["Pore "+ str(i+1) +" diameter (nm)"] = form%val
            data["Exterior"]["Pore "+ str(i+1) +" diameter (nm)"] = " "
        data["Interior"]["Solvent reservoir z-dimension (nm)"] = " "
        data["Exterior"]["Solvent reservoir z-dimension (nm)"] = form%self.reservoir()
        for i,val in enumerate(self.volume(is_sum=False)):
            data["Interior"]["Pore "+ str(i+1) +" volume (nm^3)"] = form%val
            data["Exterior"]["Pore "+ str(i+1) +" volume (nm^3)"] = " "
        data["Interior"]["Pore volume (nm^3)"] = form%self.volume()
        data["Exterior"]["Pore volume (nm^3)"] = " "
        data["Interior"]["Solvent reservoir volume (nm^3)"] = " "
        data["Exterior"]["Solvent reservoir volume (nm^3)"] = "2 * "+form%(self.box()[0]*self.box()[1]*self.reservoir())
        for i,val in enumerate(self.surface(is_sum=False)["in"]):
            data["Interior"]["Surface "+ str(i+1) +" area (nm^2)"] = form%val
            data["Exterior"]["Surface "+ str(i+1) +" area (nm^2)"] = " "
        data["Interior"]["Surface area (nm^2)"] = form%surf["in"]
        data["Exterior"]["Surface area (nm^2)"] = "2 * "+form%(surf["ex"]/2)

        data["Interior"]["Surface chemistry - Before Functionalization"] = " "
        data["Exterior"]["Surface chemistry - Before Functionalization"] = " "
        data["Interior"]["    Number of single silanol groups"] = "%i"%sum([1 for x in self._pore.get_sites().values() if len(x["o"])==1 and x["type"]=="in"])
        data["Exterior"]["    Number of single silanol groups"] = "%i"%sum([1 for x in self._pore.get_sites().values() if len(x["o"])==1 and x["type"]=="ex"])
        data["Interior"]["    Number of geminal silanol groups"] = "%i"%sum([1 for x in self._pore.get_sites().values() if len(x["o"])==2 and x["type"]=="in"])
        data["Exterior"]["    Number of geminal silanol groups"] = "%i"%sum([1 for x in self._pore.get_sites().values() if len(x["o"])==2 and x["type"]=="ex"])
        data["Interior"]["    Number of siloxane bridges"] = "%i"%allocation["SLX"]["in"][0] if "SLX" in allocation else "0"
        data["Exterior"]["    Number of siloxane bridges"] = "%i"%allocation["SLX"]["ex"][0] if "SLX" in allocation else "0"
        data["Interior"]["    Total number of OH groups"] = "%i"%allocation["Hydro"]["in"][0]
        data["Exterior"]["    Total number of OH groups"] = "%i"%allocation["Hydro"]["ex"][0]
        data["Interior"]["    Overall hydroxylation (mumol/m^2)"] = form%allocation["Hydro"]["in"][2]
        data["Exterior"]["    Overall hydroxylation (mumol/m^2)"] = form%allocation["Hydro"]["ex"][2]

        try:
            self._pore.sites_attach_mol
            for i in  self._pore.sites_attach_mol:
                if i != 20:
                    data["Interior"]["Surface chemistry - Before Functionalization (Pore " + str(i+1) +")"] = " "
                    data["Exterior"]["Surface chemistry - Before Functionalization (Pore " + str(i+1) +")"] = " "
                    data["Interior"]["    Pore " + str(i+1) + " Number of single silanol groups"] = "%i"%self._pore.sites_attach_mol[i]["SL"]
                    data["Exterior"]["    Pore " + str(i+1) + " Number of single silanol groups"] = " "
                    data["Interior"]["    Pore " + str(i+1) + " Number of geminal silanol groups"] = "%i"%self._pore.sites_attach_mol[i]["SLG"]
                    data["Exterior"]["    Pore " + str(i+1) + " Number of geminal silanol groups"] = " "
                    data["Interior"]["    Pore " + str(i+1) + " Number of siloxane bridges"] = "%i"%self._pore.sites_attach_mol[i]["SLX"] if "SLX" in self._pore.sites_attach_mol[i] else "0"
                    data["Exterior"]["    Pore " + str(i+1) + " Number of siloxane bridges"] = " "
                    data["Interior"]["    Pore " + str(i+1) + " Total number of OH groups"] = "%i"%(self._pore.sites_attach_mol[i]["SL"]+2*self._pore.sites_attach_mol[i]["SLG"])
                    data["Exterior"]["    Pore " + str(i+1) + " Total number of OH groups"] = " "
                    data["Interior"]["    Pore " + str(i+1) + " Overall hydroxylation (mumol/m^2)"] = form%(pms.utils.mols_to_mumol_m2(self._pore.sites_attach_mol[i]["SL"]+2*self._pore.sites_attach_mol[i]["SLG"],self.surface(is_sum=False)["in"][i]))
                    data["Exterior"]["    Pore " + str(i+1) + " Overall hydroxylation (mumol/m^2)"] = " "
                elif i == 20:
                    data["Interior"]["Surface chemistry - Before Functionalization (Unassigned binding sites)"] = " "
                    data["Exterior"]["Surface chemistry - Before Functionalization (Unassigned binding sites)"] = " "
                    data["Interior"]["    Unassigned binding sites" + " Number of single silanol groups"] = "%i"%self._pore.sites_attach_mol[i]["SL"]
                    data["Exterior"]["    Unassigned binding sites" + " Number of single silanol groups"] = " "
                    data["Interior"]["    Unassigned binding sites" + " Number of geminal silanol groups"] = "%i"%self._pore.sites_attach_mol[i]["SLG"]
                    data["Exterior"]["    Unassigned binding sites" + " Number of geminal silanol groups"] = " "
                    data["Interior"]["    Unassigned binding sites" + " Number of siloxane bridges"] = "%i"%self._pore.sites_attach_mol[i]["SLX"] if "SLX" in self._pore.sites_attach_mol[i] else "0"
                    data["Exterior"]["    Unassigned binding sites" + " Number of siloxane bridges"] = " "
                    data["Interior"]["    Unassigned binding sites" + " Total number of OH groups"] = "%i"%(self._pore.sites_attach_mol[i]["SL"]+2*self._pore.sites_attach_mol[i]["SLG"])
                    data["Exterior"]["    Unassigned binding sites" + " Total number of OH groups"] = " "
                    data["Interior"]["    Unassigned binding sites" + " Overall hydroxylation (mumol/m^2)"] = form%(pms.utils.mols_to_mumol_m2(self._pore.sites_attach_mol[i]["SL"]+2*self._pore.sites_attach_mol[i]["SLG"],self.surface(is_sum=False)["in"][i]))
                    data["Exterior"]["    Unassigned binding sites" + " Overall hydroxylation (mumol/m^2)"] = " "
        except AttributeError:
            pass


        data["Interior"]["Surface chemistry - After Functionalization"] = " "
        data["Exterior"]["Surface chemistry - After Functionalization"] = " "
        for mol in allocation.keys():
            if mol not in ["SL", "SLG", "SLX", "Hydro", "OH"]:
                data["Interior"]["    Number of "+mol+" groups"] = "%i"%allocation[mol]["in"][0]
                data["Exterior"]["    Number of "+mol+" groups"] = "%i"%allocation[mol]["ex"][0]
                data["Interior"]["    "+mol+" density (mumol/m^2)"] = form%allocation[mol]["in"][2]
                data["Exterior"]["    "+mol+" density (mumol/m^2)"] = form%allocation[mol]["ex"][2]
        data["Interior"]["    Bonded-phase density (mumol/m^2)"] = form%(allocation["Hydro"]["in"][2]-allocation["OH"]["in"][2])
        data["Exterior"]["    Bonded-phase density (mumol/m^2)"] = form%(allocation["Hydro"]["ex"][2]-allocation["OH"]["ex"][2])
        data["Interior"]["    Number of residual OH groups"] = "%i"%allocation["OH"]["in"][0]
        data["Exterior"]["    Number of residual OH groups"] = "%i"%allocation["OH"]["ex"][0]
        data["Interior"]["    Residual hydroxylation (mumol/m^2)"] = form%allocation["OH"]["in"][2]
        data["Exterior"]["    Residual hydroxylation (mumol/m^2)"] = form%allocation["OH"]["ex"][2]

        try:
            self._pore.sites_attach_mol
            for i in  self._pore.sites_attach_mol:
                for mol in self._pore.sites_attach_mol[i]:
                    if (mol not in ["SL", "SLG", "SLX", "Hydro", "OH"]):
                        if i == 20:
                            data["Interior"]["Surface chemistry - After Functionalization (Unassigned binding sites)"] = " "
                            data["Exterior"]["Surface chemistry - After Functionalization (Unassigned binding sites)"] = " "
                            data["Interior"]["    Unassigned binding sites" + " Number of "+mol+" groups"] = "%i"%self._pore.sites_attach_mol[i][mol]
                            data["Exterior"]["    Unassigned binding sites" + " Number of "+mol+" groups"] = " "
                            data["Interior"]["    Unassigned binding sites" + " "+mol+" density (mumol/m^2)"] = form%(pms.utils.mols_to_mumol_m2(self._pore.sites_attach_mol[i][mol],self.surface(is_sum=False)["in"][i]))
                            data["Exterior"]["    Unassigned binding sites" + " "+mol+" density (mumol/m^2)"] = " "
                        else:
                            data["Interior"]["Surface chemistry - After Functionalization (Pore " + str(i+1) +")"] = " "
                            data["Exterior"]["Surface chemistry - After Functionalization (Pore " + str(i+1) +")"] = " "
                            data["Interior"]["    Pore " + str(i+1) + " Number of "+mol+" groups"] = "%i"%self._pore.sites_attach_mol[i][mol]
                            data["Exterior"]["    Pore " + str(i+1) + " Number of "+mol+" groups"] = " "
                            data["Interior"]["    Pore " + str(i+1) + " "+mol+" density (mumol/m^2)"] = form%(pms.utils.mols_to_mumol_m2(self._pore.sites_attach_mol[i][mol],self.surface(is_sum=False)["in"][i]))
                            data["Exterior"]["    Pore " + str(i+1) + " "+mol+" density (mumol/m^2)"] = " "
        except AttributeError:
            pass
        return pd.DataFrame.from_dict(data)


class PoreCylinder(PoreKit):
    """This class carves a cylindrical pore system out of a
    :math:`\\beta`-cristobalite block.

    Parameters
    ----------
    size : list
        Size of the silicon-oxygen-grid
    diam : float
        Cylinder diameter
    res : float, optional
        Reservoir size on each side
    hydro : list, optional
        Hydroxilation degree for interior and exterior of the pore in
        :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`

    Examples
    --------
    Following example generates a cylindrical pore with a diameter of 4nm,
    reservoirs of 5nm on each side and a surface functionalized with TMS

    .. code-block:: python

        import porems as pms

        pore = pms.PoreCylinder([6, 6, 6], 4, 5)

        pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in")
        pore.attach(pms.gen.tms(), 0, [0, 1], 20, "ex")

        pore.finalize()

        pore.store("output/")
    """
    def __init__(self, size, diam, res=5, hydro=[0, 0]):
        # Call super class
        super(PoreCylinder, self).__init__()

        # Create structure
        self.structure(pms.BetaCristobalit().generate(size, "z"))
        self.build()

        # Create reservoir
        self.exterior(res, hydro=hydro[1])

        # Add pore shape
        self.add_shape(self.shape_cylinder(diam), hydro=hydro[0])
        self.prepare()


class PoreSlit(PoreKit):
    """This class carves a slit-pore out of a :math:`\\beta`-cristobalite block.

    Parameters
    ----------
    size : list
        Size of the silicon-oxygen-grid
    height : float
        Pore height
    res : float, optional
        Reservoir size on each side
    hydro: list, optional
        Hydroxilation degree for interior and exterior of the pore in
        :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`

    Examples
    --------
    Following example generates a slit-pore with a height of 4nm functionalized
    with TMS

    .. code-block:: python

        import porems as pms

        pore = pms.PoreSlit([6, 6, 6], 4)

        pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in")

        pore.finalize()

        pore.store("output/")
    """
    def __init__(self, size, height, res=0, hydro=[0, 0]):
        # Call super class
        super(PoreSlit, self).__init__()

        # Create structure
        self.structure(pms.BetaCristobalit().generate(size, "z"))
        self.build()

        # Create reservoir
        self.exterior(res, hydro=hydro[1])

        # Add pore shape
        self.add_shape(self.shape_slit(height), hydro=hydro[0])
        self.prepare()


class PoreCapsule(PoreKit):
    """This class carves a capsule pore system out of a
    :math:`\\beta`-cristobalite block.

    Parameters
    ----------
    size : list
        Size of the silicon-oxygen-grid
    diam : float
        Cylinder diameter
    sep : float
        Separation length between capsule
    res : float, optional
        Reservoir size on each side
    hydro: list, optional
        Hydroxilation degree for interior and exterior of the pore in
        :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`

    Examples
    --------
    Following example generates a capsule pore with a diameter of 4nm, a
    separation distance between the capsules of 2nm, reservoirs of 5nm on each
    side and a surface functionalized with TMS

    .. code-block:: python

        import porems as pms

        pore = pms.PoreCapsule([6, 6, 12], 4, 2)

        pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in")
        pore.attach(pms.gen.tms(), 0, [0, 1], 20, "ex")

        pore.finalize()

        pore.store("output/")
    """
    def __init__(self, size, diam, sep, res=5, hydro=[0, 0]):
        # Call super class
        super(PoreCapsule, self).__init__()

        # Create structure
        self.structure(pms.BetaCristobalit().generate(size, "z"))
        self.build()

        # Create reservoir
        self.exterior(res, hydro=hydro[1])

        # Add pore shape
        center = [size[0]/2, size[1]/2]
        len_cyl = (size[2]-diam-sep)/2

        centroids = []
        centroids.append(center+[0+len_cyl/2])
        centroids.append(center+[len_cyl])
        centroids.append(center+[len_cyl+diam+sep])
        centroids.append(center+[size[2]-len_cyl/2])

        sections = []
        sections.append({"x": [], "y": [], "z": [0, len_cyl]})
        sections.append({"x": [], "y": [], "z": [len_cyl, len_cyl+diam/2+sep/2]})
        sections.append({"x": [], "y": [], "z": [len_cyl+diam/2+sep/2, len_cyl+diam+sep]})
        sections.append({"x": [], "y": [], "z": [size[2]-len_cyl, size[2]]})

        self.add_shape(self.shape_cylinder(diam, len_cyl, centroids[0]), section=sections[0], hydro=hydro[0])
        self.add_shape(self.shape_sphere(  diam,          centroids[1]), section=sections[1], hydro=hydro[0])
        self.add_shape(self.shape_sphere(  diam,          centroids[2]), section=sections[2], hydro=hydro[0])
        self.add_shape(self.shape_cylinder(diam, len_cyl, centroids[3]), section=sections[3], hydro=hydro[0])

        self.prepare()


class PoreAmorphCylinder(PoreKit):
    """This class carves a cylindrical pore system out of an amorph
    :math:`\\beta`-cristobalite block from
    `Vink et al. <http://doi.org/10.1103/PhysRevB.67.245201>`_ with dimensions
    [9.605, 9.605, 9.605] (x, y, z).

    Parameters
    ----------
    diam : float
        Cylinder diameter
    res : float, optional
        Reservoir size on each side
    hydro: list, optional
        Hydroxilation degree for interior and exterior of the pore in
        :math:`\\frac{\\mu\\text{mol}}{\\text{m}^2}`

    Examples
    --------
    Following example generates a cylindrical pore with a diameter of 4nm,
    reservoirs of 5nm on each side and a surface functionalized with TMS

    .. code-block:: python

        import porems as pms

        pore = pms.PoreAmorphCylinder(4, 5)

        pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in")
        pore.attach(pms.gen.tms(), 0, [0, 1], 20, "ex")

        pore.finalize()

        pore.store("output/")
    """
    def __init__(self, diam, res=5, hydro=[0, 0]):
        # Call super class
        super(PoreAmorphCylinder, self).__init__()

        # Create structure
        self.structure(pms.Molecule(inp=os.path.split(__file__)[0]+"/templates/amorph.gro"))
        self.build(bonds=[0.160-0.02, 0.160+0.02])
        self._matrix.split(57790, 2524)

        # Create reservoir
        self.exterior(res, hydro=hydro[1])

        # Add pore shape
        self.add_shape(self.shape_cylinder(diam), section={"x": [], "y": [], "z": [-1, 10]}, hydro=hydro[0])
        self.prepare()


class PoreMultiChannel(PoreKit):
    """Convenience class that carves *N* parallel cylindrical channels into a
    :math:`\\beta`-cristobalite block.

    The channels are arranged in a regular 1-D row along the x-axis and share
    the same y-coordinate as the block centroid.  Their centre-to-centre
    spacing is ``spacing`` (defaulting to ``diam + 0.5`` nm).

    Parameters
    ----------
    size : list
        ``[nx, ny, nz]`` repeat counts for the BetaCristobalit block (nm).
    n_channels : int
        Number of parallel channels to carve.
    diam : float
        Channel diameter in nm.
    spacing : float, optional
        Centre-to-centre channel spacing in nm.  Defaults to ``diam + 0.5``.
    res : float, optional
        Reservoir length on each side in nm (0 = no reservoir).
    hydro : list, optional
        Hydroxylation density ``[interior, exterior]`` in
        :math:`\\mu\\text{mol}\\,\\text{m}^{-2}`.

    Examples
    --------
    .. code-block:: python

        import porems as pms

        pore = pms.PoreMultiChannel([8, 4, 8], 3, 2.0, spacing=3.0, res=5)
        pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in")
        pore.finalize()
        pore.store("output/multi/")
    """
    def __init__(self, size, n_channels, diam, spacing=None, res=5, hydro=[0, 0]):
        super(PoreMultiChannel, self).__init__()

        if spacing is None:
            spacing = diam + 0.5

        # Build crystal block
        self.structure(pms.BetaCristobalit().generate(size, "z"))
        self.build()

        # Reservoir / exterior surface
        self.exterior(res, hydro=hydro[1])

        # Place channels in a row along x, centred on the block
        cy = self._box[1] / 2
        cz = self._box[2] / 2
        total_width = (n_channels - 1) * spacing
        x_start = self._box[0] / 2 - total_width / 2

        for i in range(n_channels):
            cx = x_start + i * spacing
            centroid = [cx, cy, cz]
            self.add_shape(
                self.shape_cylinder(diam, centroid=centroid),
                hydro=hydro[0],
            )

        self.prepare()


class PoreCone(PoreKit):
    """Convenience class that carves a conical pore into a β-cristobalite block.

    The cone transitions linearly from ``diam_in`` at one end to ``diam_out``
    at the other.  Both end diameters are on the interior surface.

    Parameters
    ----------
    size : list
        ``[nx, ny, nz]`` block repeat counts (nm).
    diam_in : float
        Diameter at the narrow end (nm).
    diam_out : float
        Diameter at the wide end (nm).
    res : float, optional
        Reservoir length on each side (nm). Default 5.
    hydro : list, optional
        Hydroxylation density ``[interior, exterior]`` in μmol m⁻².

    Examples
    --------
    .. code-block:: python

        import porems as pms

        pore = pms.PoreCone([8, 8, 10], 2.0, 4.0, res=5)
        pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in")
        pore.finalize()
        pore.store("output/cone/")
    """
    def __init__(self, size, diam_in, diam_out, res=5, hydro=[0, 0]):
        super(PoreCone, self).__init__()
        self.structure(pms.BetaCristobalit().generate(size, "z"))
        self.build()
        self.exterior(res, hydro=hydro[1])
        self.add_shape(self.shape_cone(diam_in, diam_out), hydro=hydro[0])
        self.prepare()
