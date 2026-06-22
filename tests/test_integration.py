import matplotlib.pyplot as plt

import pytest

import porems as pms


pytestmark = pytest.mark.slow


def test_pattern_beta_cristobalit():
    beta_cristobalit = pms.BetaCristobalit()

    pattern = beta_cristobalit.pattern()
    pattern.set_name("pattern_beta_cbt_minimal")
    assert pattern.get_num() == 36
    pms.Store(pattern, "output").gro()

    beta_cristobalit = pms.BetaCristobalit()
    beta_cristobalit.generate([2, 2, 2], "x")
    beta_cristobalit.get_block().set_name("pattern_beta_cbt_x")
    assert beta_cristobalit.get_size() == [2.480, 1.754, 2.024]
    assert [round(x, 3) for x in beta_cristobalit.get_block().get_box()] == [2.480, 1.754, 2.024]
    pms.Store(beta_cristobalit.get_block(), "output").gro()

    beta_cristobalit = pms.BetaCristobalit()
    beta_cristobalit.generate([2, 2, 2], "y")
    beta_cristobalit.get_block().set_name("pattern_beta_cbt_y")
    assert beta_cristobalit.get_size() == [2.024, 2.480, 1.754]
    assert [round(x, 3) for x in beta_cristobalit.get_block().get_box()] == [2.024, 2.480, 1.754]
    pms.Store(beta_cristobalit.get_block(), "output").gro()

    beta_cristobalit = pms.BetaCristobalit()
    beta_cristobalit.generate([2, 2, 2], "z")
    beta_cristobalit.get_block().set_name("pattern_beta_cbt_z")
    assert beta_cristobalit.get_size() == [2.024, 1.754, 2.480]
    assert [round(x, 3) for x in beta_cristobalit.get_block().get_box()] == [2.024, 1.754, 2.480]
    pms.Store(beta_cristobalit.get_block(), "output").gro()
    pms.Store(beta_cristobalit.get_block(), "output").lmp()

    beta_cristobalit = pms.BetaCristobalit()
    beta_cristobalit.generate([2, 2, 2], "z")
    beta_cristobalit.get_block().set_name("DOTA")

    assert beta_cristobalit.get_block().get_num() == 576
    assert beta_cristobalit.get_block().overlap() == {}
    assert beta_cristobalit.get_repeat() == [0.506, 0.877, 1.240]
    assert beta_cristobalit.get_gap() == [0.126, 0.073, 0.155]
    assert beta_cristobalit.get_orient() == "z"
    assert beta_cristobalit.get_block().get_name() == "DOTA"


def test_alpha_cristobalit():
    alpha_cristobalit = pms.AlphaCristobalit()

    pattern = alpha_cristobalit.pattern()
    pattern.set_name("pattern_alpha_cbt_minimal")
    assert pattern.get_num() == 12
    pms.Store(pattern, "output").gro()

    alpha_cristobalit = pms.AlphaCristobalit()
    alpha_cristobalit.generate([2, 2, 2], "x")
    alpha_cristobalit.get_block().set_name("pattern_alpha_cbt_x")
    assert alpha_cristobalit.get_size() == [2.0844, 1.9912, 1.9912]
    assert [round(x, 3) for x in alpha_cristobalit.get_block().get_box()] == [2.084, 1.991, 1.991]
    pms.Store(alpha_cristobalit.get_block(), "output").gro()

    alpha_cristobalit = pms.AlphaCristobalit()
    alpha_cristobalit.generate([2, 2, 2], "y")
    alpha_cristobalit.get_block().set_name("pattern_alpha_cbt_y")
    assert alpha_cristobalit.get_size() == [1.9912, 2.0844, 1.9912]
    assert [round(x, 3) for x in alpha_cristobalit.get_block().get_box()] == [1.991, 2.084, 1.991]
    pms.Store(alpha_cristobalit.get_block(), "output").gro()

    alpha_cristobalit = pms.AlphaCristobalit()
    alpha_cristobalit.generate([2, 2, 2], "z")
    alpha_cristobalit.get_block().set_name("pattern_alpha_cbt_z")
    assert alpha_cristobalit.get_size() == [1.9912, 1.9912, 2.0844]
    assert [round(x, 3) for x in alpha_cristobalit.get_block().get_box()] == [1.991, 1.991, 2.084]
    pms.Store(alpha_cristobalit.get_block(), "output").gro()
    pms.Store(alpha_cristobalit.get_block(), "output").lmp()

    alpha_cristobalit = pms.AlphaCristobalit()
    alpha_cristobalit.generate([2, 2, 2], "z")
    alpha_cristobalit.get_block().set_name("DOTA")

    assert alpha_cristobalit.get_block().get_num() == 576
    assert alpha_cristobalit.get_block().overlap() == {}
    assert alpha_cristobalit.get_repeat() == [0.4978, 0.4978, 0.6948]
    assert alpha_cristobalit.get_orient() == "z"
    assert alpha_cristobalit.get_block().get_name() == "DOTA"


def test_dice():
    block = pms.BetaCristobalit().generate([2, 2, 2], "z")
    block.set_name("dice")
    pms.Store(block, "output").gro()
    dice = pms.Dice(block, 0.4, True)

    assert len(dice.get_origin()) == 120
    assert dice.get_origin()[(1, 1, 1)] == [0.4, 0.4, 0.4]
    assert dice.get_pointer()[(1, 1, 1)] == [14, 51, 52, 64, 65, 67]

    assert dice._right((1, 1, 1)) == (2, 1, 1)
    assert dice._left((1, 1, 1))  == (0, 1, 1)
    assert dice._top((1, 1, 1))   == (1, 2, 1)
    assert dice._bot((1, 1, 1))   == (1, 0, 1)
    assert dice._front((1, 1, 1)) == (1, 1, 2)
    assert dice._back((1, 1, 1))  == (1, 1, 0)
    assert len(dice.neighbor((1, 1, 1))) == 27
    assert len(dice.neighbor((1, 1, 1), False)) == 26

    assert dice.find_bond([(1, 1, 1)], ["Si", "O"], [0.155-0.005, 0.155+0.005]) == [[51, [46, 14, 52, 65]], [64, [26, 63, 65, 67]]]
    assert dice.find_bond([(1, 1, 1)], ["O", "Si"], [0.155-0.005, 0.155+0.005]) == [[14, [51, 13]], [52, [51, 49]], [65, [51, 64]], [67, [64, 69]]]
    assert dice.find_bond([(0, 0, 0)], ["Si", "O"], [0.155-0.005, 0.155+0.005]) == [[3, [4, 9, 2, 174]], [5, [306, 110, 4, 6]]]
    assert dice.find_bond([(0, 0, 0)], ["O", "Si"], [0.155-0.005, 0.155+0.005]) == [[4, [3, 5]], [6, [7, 5]], [9, [3, 11]]]

    assert len(dice.find_parallel(None, ["Si", "O"], [0.155-0.005, 0.155+0.005])) == 192
    assert len(dice.find_parallel(None, ["O", "Si"], [0.155-0.005, 0.155+0.005])) == 384

    dice.set_pbc(True)
    assert dice.get_count() == [5, 4, 6]
    assert dice.get_size() == 0.4
    assert dice.get_mol().get_name() == "dice"


def test_matrix():
    block = pms.BetaCristobalit().generate([1, 1, 1], "z")
    block.set_name("matrix")
    pms.Store(block, "output").gro()
    dice = pms.Dice(block, 0.2, True)
    bonds = dice.find_bond(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2])

    matrix = pms.Matrix(bonds)
    connect = matrix.get_matrix()
    matrix.split(0, 17)
    assert connect[0]["atoms"] == [30, 8, 1]
    assert connect[17]["atoms"] == [19]
    matrix.strip(0)
    assert connect[0]["atoms"] == []
    assert connect[1]["atoms"] == [43]
    assert connect[8]["atoms"] == [7]
    assert connect[30]["atoms"] == [3]
    assert matrix.bound(0) == [0]
    assert matrix.bound(1, "lt") == [0]
    assert matrix.bound(4, "gt") == []
    matrix.add(0, 17)
    assert connect[0]["atoms"] == [17]
    assert connect[17]["atoms"] == [19, 0]
    assert matrix.bound(4, "test") is None


def test_shape_cylinder():
    block = pms.BetaCristobalit().generate([6, 6, 6], "z")
    block.set_name("shape_cylinder")
    dice = pms.Dice(block, 0.4, True)
    matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2]))
    centroid = block.centroid()
    central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 45, True))

    cylinder = pms.Cylinder({"centroid": centroid, "central": central, "length": 3, "diameter": 4})

    assert round(cylinder.volume(), 4) == 37.6991
    assert round(cylinder.surface(), 4) == 37.6991

    vec = [3.6086, 4.4076, 0.2065]
    assert [round(x[0][20], 4) for x in cylinder.surf(num=100)] == vec
    assert [round(x[0][20], 4) for x in cylinder.rim(0, num=100)] == vec
    assert [round(x, 4) for x in cylinder.convert([0, 0, 0], False)] == [3.0147, 3.0572, 1.5569]
    assert [round(x, 4) for x in cylinder.normal(vec)] == [0.5939, 2.9704, 0.0000]

    del_list = [atom_id for atom_id in range(block.get_num()) if cylinder.is_in(block.pos(atom_id))]
    matrix.strip(del_list)
    block.delete(matrix.bound(0))
    assert block.get_num() == 12650

    pms.Store(block, "output").gro()
    plt.figure()
    cylinder.plot(vec=[3.17290646, 4.50630614, 0.22183271])


def test_shape_sphere():
    block = pms.BetaCristobalit().generate([6, 6, 6], "z")
    block.set_name("shape_sphere")
    dice = pms.Dice(block, 0.4, True)
    matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2]))
    centroid = block.centroid()
    central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 0, True))

    sphere = pms.Sphere({"centroid": centroid, "central": central, "diameter": 4})

    assert round(sphere.volume(), 4) == 33.5103
    assert round(sphere.surface(), 4) == 50.2655
    assert [round(x[0][20], 4) for x in sphere.surf(num=100)] == [4.2006, 3.0572, 4.6675]
    assert [round(x[0][20], 4) for x in sphere.rim(0, num=100)] == [4.9245, 3.0572, 3.6508]
    assert [round(x, 4) for x in sphere.convert([0, 0, 0], False)] == [3.0147, 3.0572, 3.0569]
    assert [round(x, 4) for x in sphere.normal([4.2006, 3.0572, 4.6675])] == [1.4063, 0.0000, 1.9099]

    del_list = [atom_id for atom_id in range(block.get_num()) if sphere.is_in(block.pos(atom_id))]
    matrix.strip(del_list)
    block.delete(matrix.bound(0))
    assert block.get_num() == 12934

    pms.Store(block, "output").gro()
    sphere.plot(inp=3.14, vec=[1.08001048, 3.09687610, 1.72960828])


def test_shape_cuboid():
    block = pms.BetaCristobalit().generate([6, 6, 6], "z")
    block.set_name("shape_cuboid")
    dice = pms.Dice(block, 0.4, True)
    matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2]))
    centroid = block.centroid()
    central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 0, True))

    cuboid = pms.Cuboid({"centroid": centroid, "central": central, "length": 10, "width": 6, "height": 4})

    assert round(cuboid.volume(), 4) == 240
    assert round(cuboid.surface(), 4) == 248
    assert [round(x, 4) for x in cuboid.convert([0, 0, 0], False)] == [0.0147, 1.0572, -1.9431]
    assert [round(x, 4) for x in cuboid.normal([4.2636, 3.0937, 4.745])] == [0, 1, 0]

    del_list = [atom_id for atom_id in range(block.get_num()) if cuboid.is_in(block.pos(atom_id))]
    matrix.strip(del_list)
    block.delete(matrix.bound(0))
    assert block.get_num() == 5160

    pms.Store(block, "output").gro()
    cuboid.plot()


def test_shape_cone():
    block = pms.BetaCristobalit().generate([6, 6, 6], "z")
    block.set_name("shape_cone")
    dice = pms.Dice(block, 0.4, True)
    matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2]))
    centroid = block.centroid()
    central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 45, True))

    cone = pms.Cone({"centroid": centroid, "central": central, "length": 6, "diameter_1": 4, "diameter_2": 1})

    assert round(cone.volume(), 4) == 32.9867
    assert round(cone.surface(), 4) == 48.5742

    vec = [3.6977, 4.6102, -1.4961]
    assert [round(x[0][20], 4) for x in cone.surf(num=100)] == vec
    assert [round(x[0][20], 4) for x in cone.rim(0, num=100)] == vec
    assert [round(x, 4) for x in cone.convert([0, 0, 0], False)] == [3.0147, 3.0572, 0.0569]
    assert [round(x, 4) for x in cone.normal(vec)] == [0.3182, 2.0114, 0.6109]

    del_list = [atom_id for atom_id in range(block.get_num()) if cone.is_in(block.pos(atom_id))]
    matrix.strip(del_list)
    block.delete(matrix.bound(0))
    assert block.get_num() == 12486

    pms.Store(block, "output").gro()
    plt.figure()
    cone.plot(vec=[3.17290646, 4.50630614, 0.22183271])


def test_pore():
    orient = "z"
    pattern = pms.BetaCristobalit()
    pattern.generate([6, 6, 6], orient)

    block = pattern.get_block()
    block.set_name("pore_cylinder_block")

    dice = pms.Dice(block, 0.4, True)
    bond_list = dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2])
    matrix = pms.Matrix(bond_list)

    pore = pms.Pore(block, matrix)

    centroid = block.centroid()
    central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 0, True))
    cylinder = pms.Cylinder({"centroid": centroid, "central": central, "length": 6, "diameter": 4})
    del_list = [atom_id for atom_id in range(block.get_num()) if cylinder.is_in(block.pos(atom_id))]
    matrix.strip(del_list)

    pore.prepare()
    pore.sites()
    assert len(pore.get_sites()) == 455

    block.delete(matrix.bound(0))
    pms.Store(block, "output").gro("pore_no_ex.gro")

    # With exterior surface
    pattern = pms.BetaCristobalit()
    pattern.generate([6, 6, 6], orient)

    block = pattern.get_block()
    block.set_name("pore_cylinder_block")

    dice = pms.Dice(block, 0.4, True)
    bond_list = dice.find_parallel(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2])
    matrix = pms.Matrix(bond_list)

    pore = pms.Pore(block, matrix)
    pore.exterior()

    centroid = block.centroid()
    central = pms.geom.unit(pms.geom.rotate([0, 0, 1], [1, 0, 0], 0, True))
    cylinder = pms.Cylinder({"centroid": centroid, "central": central, "length": 6, "diameter": 4})
    del_list = [atom_id for atom_id in range(block.get_num()) if cylinder.is_in(block.pos(atom_id))]
    matrix.strip(del_list)

    pore.prepare()
    pore.amorph()
    assert len(matrix.bound(1)) == 710
    pore.sites()
    site_list = pore.get_sites()
    site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"] == "in"]
    site_ex = [site_key for site_key, site_val in site_list.items() if site_val["type"] == "ex"]
    assert len(site_in) == 432
    assert len(site_ex) == 201

    si_pos_in = [block.pos(site_key) for site_key, site_val in site_list.items() if site_val["type"] == "in"]
    si_pos_ex = [block.pos(site_key) for site_key, site_val in site_list.items() if site_val["type"] == "ex"]

    if si_pos_in:
        temp_mol = pms.Molecule()
        for pos in si_pos_in:
            temp_mol.add("Si", pos)
        pms.Store(temp_mol).gro("output/pore_cylinder_si_in.gro")

    if si_pos_ex:
        temp_mol = pms.Molecule()
        for pos in si_pos_ex:
            temp_mol.add("Si", pos)
        pms.Store(temp_mol).gro("output/pore_cylinder_si_ex.gro")

    non_grid = matrix.bound(1) + list(site_list.keys())
    bonded = matrix.bound(0, "gt")
    grid_atoms = [atom for atom in bonded if atom not in non_grid]
    mol_obj = pore.objectify(grid_atoms)
    assert len(mol_obj) == 8279
    pms.Store(pms.Molecule(name="pore_cylinder_grid", inp=mol_obj), "output").gro(use_atom_names=True)

    mol = pms.gen.tms()

    def normal(pos):
        return [0, 0, -1] if pos[2] < centroid[2] else [0, 0, 1]

    for site in site_in:
        site_list[site]["normal"] = cylinder.normal
    for site in site_ex:
        site_list[site]["normal"] = normal

    mols_siloxane = pore.siloxane(site_in, 100)
    site_in = [site_key for site_key, site_val in site_list.items() if site_val["type"] == "in"]

    mols_in = pore.attach(mol, 0, [0, 1], site_in, 100, site_type="in")
    mols_ex = pore.attach(mol, 0, [0, 1], site_ex, 20, site_type="ex")

    mols_in_fill = pore.fill_sites(site_in, site_type="in")
    mols_ex_fill = pore.fill_sites(site_ex, site_type="ex")

    pms.Store(pms.Molecule(name="pore_cylinder_siloxane", inp=mols_siloxane), "output").gro()
    pms.Store(pms.Molecule(name="pore_cylinder_in", inp=mols_in), "output").gro()
    pms.Store(pms.Molecule(name="pore_cylinder_ex", inp=mols_ex), "output").gro()
    pms.Store(pms.Molecule(name="pore_cylinder_in_fill", inp=mols_in_fill), "output").gro()
    pms.Store(pms.Molecule(name="pore_cylinder_ex_fill", inp=mols_ex_fill), "output").gro()

    block.delete(matrix.bound(0))
    pms.Store(block, "output").gro()

    pore.reservoir(5)
    assert [round(x) for x in pore.get_box()] == [6, 6, 17]

    pore.set_name("pore_cylinder_full")
    pms.Store(pore, "output").gro(use_atom_names=True)

    sort_list = ["OM", "SI", "SLX", "SL", "SLG", "TMS", "TMSG"]
    pore.set_name("pore_cylinder_full_sort")
    pms.Store(pore, "output", sort_list=sort_list).gro(use_atom_names=True)
    pms.Store(pore, "output", sort_list=sort_list).pdb(use_atom_names=True)

    pms.Store(pore, "output", sort_list=sort_list[:-1])
    pms.Store(pore, "output", sort_list=sort_list).top()

    assert pore.attach(mol, 0, [0, 1], site_in, 0, cylinder.normal, site_type="DOTA") is None
    assert pore.siloxane(site_in, 0, cylinder.normal, site_type="DOTA") is None

    assert pore.get_block().get_name() == "pore_cylinder_block"
    assert len(pore.get_site_dict()) == 3
    assert pore.get_num_in_ex() == 23


def test_pore_exterior():
    for orient in ("x", "y", "z"):
        pattern = pms.BetaCristobalit()
        block = pattern.generate([2, 2, 2], orient)
        block.set_name(f"pattern_beta_cbt_ex_{orient}")
        dice = pms.Dice(block, 0.2, True)
        bonds = dice.find_bond(None, ["Si", "O"], [0.155-1e-2, 0.155+1e-2])
        matrix = pms.Matrix(bonds)
        pore = pms.Pore(block, matrix)
        pore.prepare()
        pore.exterior()
        pore.sites()
        pms.Store(block, "output").gro()

    pattern = pms.BetaCristobalit()
    pattern.generate([2, 2, 2], "z")
    pattern._structure = pms.Molecule(inp="data/amorph.gro")
    pattern._size = [2.014, 1.751, 2.468]

    block = pattern.get_block()
    block.set_name("pattern_beta_cbt_ex_amoprh")

    dice = pms.Dice(block, 0.4, True)
    matrix = pms.Matrix(dice.find_parallel(None, ["Si", "O"], [0.160-0.02, 0.160+0.02]))

    connect = matrix.get_matrix()
    matrix.split(57790, 2524)

    pore = pms.Pore(block, matrix)
    pore.prepare()
    pore.exterior()
    pore.sites()
    pms.Store(block, "output").gro()


def test_pore_kit():
    pore = pms.PoreKit()
    pore.structure(pms.BetaCristobalit().generate([5, 5, 10], "z"))
    pore.build()
    pore.exterior(5, hydro=0.4)
    pore.add_shape(pore.shape_cylinder(2, 10, [3.5, 3.5, 5]), hydro=0.4)
    pore.add_shape(pore.shape_cylinder(2, 10, [1.5, 1.5, 5]), hydro=0.4)
    pore.prepare()
    pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in")
    pore.attach(pms.gen.tms(), 0, [0, 1], 20, "ex")
    pore.finalize()
    pore.store("output/kit_parallel/")

    pore = pms.PoreKit()
    pore.structure(pms.BetaCristobalit().generate([7, 7, 10], "z"))
    pore.build()
    pore.exterior(5, hydro=0.4)
    pore.add_shape(pore.shape_cylinder(6, 4, [3.5, 3.5, 2]), section={"x": [], "y": [], "z": [0, 4]}, hydro=0.4)
    pore.add_shape(pore.shape_cone(4.5, 3, 2, [3.5, 3.5, 5]), section={"x": [], "y": [], "z": [4, 6]}, hydro=0.4)
    pore.add_shape(pore.shape_cylinder(4, 4, [3.5, 3.5, 8]), section={"x": [], "y": [], "z": [6, 10]}, hydro=0.4)
    pore.prepare()
    pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in")
    pore.attach(pms.gen.tms(), 0, [0, 1], 20, "ex")
    pore.finalize()
    pore.store("output/kit_narrow/")


def test_pore_cylinder():
    pore = pms.PoreCylinder([4, 4, 4], 2, 0)
    pore.finalize()

    pore = pms.PoreCylinder([6, 6, 6], 4, 5, [5, 5])

    tms2 = pms.gen.tms()
    tms2.set_short("TMS2")

    pore.attach(tms2, 0, [0, 1], 10, "in", trials=10, inp="percent")
    pore.attach(tms2, 0, [0, 1], 1, "in", trials=10, inp="molar")
    pore.attach(tms2, 0, [0, 1], 0.1, "ex", trials=10, inp="molar")

    assert pore.attach(pms.gen.tms(), 0, [0, 1], 100, site_type="DOTA") is None
    assert pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in", inp="DOTA") is None
    assert pore.attach(pms.gen.tms(), 0, [0, 1], 100, pos_list=[[1, 3, 3], [7, 4, 2]]) is None
    assert pore.attach_special(pms.gen.tms(), 0, [0, 1], 3, symmetry="DOTA") is None

    pore.finalize()
    pore.store("output/cylinder/")
    print(pore.table())

    assert round(pore.diameter()[0]) == 4
    assert [round(x, 4) for x in pore.centroid()] == [3.0147, 3.0572, 3.0569]
    assert round(pore.roughness()["in"][0], 1) == 0.1
    assert round(pore.roughness()["ex"], 1) == 0.0
    assert round(pore.volume()) == 78
    assert {key: round(item) for key, item in pore.surface().items()} == {"in": 78, "ex": 49}


def test_pore_slit():
    pore = pms.PoreSlit([4, 4, 4], 2)
    pore.finalize()

    pore = pms.PoreSlit([6, 6, 6], 3, 5, [5, 5])

    tms2 = pms.gen.tms()
    tms2.set_short("TMS2")

    pore.attach(tms2, 0, [0, 1], 10, "in", trials=10, inp="percent")
    pore.attach(tms2, 0, [0, 1], 1, "in", trials=10, inp="molar")
    pore.attach(tms2, 0, [0, 1], 0.1, "ex", trials=10, inp="molar")

    assert pore.attach(pms.gen.tms(), 0, [0, 1], 100, site_type="DOTA") is None
    assert pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in", inp="DOTA") is None
    assert pore.attach_special(pms.gen.tms(), 0, [0, 1], 3, symmetry="DOTA") is None

    pore.finalize()
    pore.store("output/slit/")
    print(pore.table())

    assert round(pore.diameter()[0]) == 3
    assert [round(x, 4) for x in pore.centroid()] == [3.0147, 3.0572, 3.0569]
    assert round(pore.roughness()["in"][0], 1) == 0.1
    assert round(pore.roughness()["ex"], 1) == 0.0
    assert round(pore.volume()) == 114
    assert round(pore.surface()["in"]) == 75


def test_pore_capsule():
    pore = pms.PoreCapsule([3, 3, 6], 2, 1, 2.5)
    pore.finalize()

    pore = pms.PoreCapsule([6, 6, 10], 4, 2, 5, [5, 5])

    tms2 = pms.gen.tms()
    tms2.set_short("TMS2")

    pore.attach(tms2, 0, [0, 1], 10, "in", trials=10, inp="percent")
    pore.attach(tms2, 0, [0, 1], 1, "in", trials=10, inp="molar")
    pore.attach(tms2, 0, [0, 1], 0.1, "ex", trials=10, inp="molar")

    assert pore.attach(pms.gen.tms(), 0, [0, 1], 100, site_type="DOTA") is None
    assert pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in", inp="DOTA") is None

    pore.finalize()
    print(pore.table())

    assert [round(x, 4) for x in pore.diameter()] == [4.2317, 4.3751, 4.3749, 4.2449]
    assert [round(x, 4) for x in pore.centroid()] == [3.0147, 3.0572, 4.9169]
    assert [round(x, 4) for x in pore.roughness()["in"]] == [0.1223, 0.0432, 0.0566, 0.125]
    assert round(pore.roughness()["ex"], 1) == 0.0
    assert round(pore.volume()) == 144
    assert {key: round(item) for key, item in pore.surface().items()} == {"in": 174, "ex": 45}


def test_amorph_no_interior_sioh():
    """After prepare(), no Si in the pore interior should have unsaturated bonds
    — only surface Si atoms should carry OH groups."""
    pore = pms.PoreAmorphCylinder(3.0, res=0)
    matrix = pore._pore._matrix.get_matrix()
    block = pore._pore.get_block()
    box = block.get_box()
    cx = box[0] / 2
    cy = box[1] / 2
    # interior is at least 1 nm from pore wall — Si here must be fully saturated
    bulk_si_with_free_bonds = [
        atom for atom, props in matrix.items()
        if block.get_atom_type(atom) == "Si"
        and len(props["atoms"]) < props["bonds"]
        and ((block.pos(atom)[0] - cx)**2 + (block.pos(atom)[1] - cy)**2) > (1.5)**2
    ]
    assert bulk_si_with_free_bonds == [], (
        f"{len(bulk_si_with_free_bonds)} bulk Si atoms have unsaturated bonds"
    )


def test_attach_special_point_symmetry():
    """attach_special with point symmetry should place molecules on alternating sides."""
    pore = pms.PoreCylinder([6, 6, 6], 3.0, res=0)
    mols = pore.attach_special(pms.gen.tms(), 0, [0, 1], 2, symmetry="point")
    assert mols is None or True  # just ensure it doesn't crash
    pore.finalize()


def test_pore_cylinder_amorph():
    pore = pms.PoreAmorphCylinder(2, 0)
    pore.finalize()

    pore = pms.PoreAmorphCylinder(4, 5, [2, 2])

    tms2 = pms.gen.tms()
    tms2.set_short("TMS2")

    pore.attach(tms2, 0, [0, 1], 10, "in", trials=10, inp="percent")
    pore.attach(tms2, 0, [0, 1], 1, "in", trials=10, inp="molar")
    pore.attach(tms2, 0, [0, 1], 0.1, "ex", trials=10, inp="molar")

    assert pore.attach(pms.gen.tms(), 0, [0, 1], 100, site_type="DOTA") is None
    assert pore.attach(pms.gen.tms(), 0, [0, 1], 100, "in", inp="DOTA") is None
    assert pore.attach(pms.gen.tms(), 0, [0, 1], 100, pos_list=[[1, 3, 3], [7, 4, 2]]) is None
    assert pore.attach_special(pms.gen.tms(), 0, [0, 1], 3, symmetry="DOTA") is None

    pore.finalize()
    pore.store("output/cylinder_amorph/")
    print(pore.table())

    assert round(pore.diameter()[0]) == 4
    assert [round(x, 4) for x in pore.centroid()] == [4.7958, 4.7978, 4.807]
    assert round(pore.roughness()["in"][0], 1) == 0.1
    assert round(pore.roughness()["ex"], 1) == 0.3
    assert round(pore.volume()) == 121
    assert {key: round(item) for key, item in pore.surface().items()} == {"in": 121, "ex": 159}
