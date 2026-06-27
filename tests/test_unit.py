import matplotlib.pyplot as plt

import porems as pms


def test_utils():
    file_link = "output/test/test.txt"

    pms.utils.mkdirp("output/test")

    with open(file_link, "w") as file_out:
        file_out.write("TEST")
    pms.utils.copy(file_link, file_link + "t")
    pms.utils.replace(file_link + "t", "TEST", "DOTA")
    with open(file_link + "t", "r") as file_in:
        for line in file_in:
            assert line == "DOTA\n"

    assert pms.utils.column([[1, 1, 1], [2, 2, 2]]) == [[1, 2], [1, 2], [1, 2]]

    pms.utils.save([1, 1, 1], file_link)
    assert pms.utils.load(file_link) == [1, 1, 1]

    assert round(pms.utils.mumol_m2_to_mols(3, 100), 4) == 180.66
    assert round(pms.utils.mols_to_mumol_m2(180, 100), 4) == 2.989
    assert round(pms.utils.mmol_g_to_mumol_m2(0.072, 512), 2) == 0.14
    assert round(pms.utils.mmol_l_to_mols(30, 1000), 4) == 18.066
    assert round(pms.utils.mols_to_mmol_l(18, 1000), 4) == 29.8904

    pms.utils.toc(pms.utils.tic(), message="Test", is_print=True)
    assert round(pms.utils.toc(pms.utils.tic(), is_print=True)) == 0


def test_geometry():
    vec_a = [1, 1, 2]
    vec_b = [0, 3, 2]

    # Correctness
    assert round(pms.geom.dot_product(vec_a, vec_b), 4) == 7
    assert round(pms.geom.length(vec_a), 4) == 2.4495
    assert [round(x, 4) for x in pms.geom.vector(vec_a, vec_b)] == [-1, 2, 0]
    assert pms.geom.vector([0, 1], [0, 0, 0]) is None          # dimension mismatch → None
    assert [round(x, 4) for x in pms.geom.unit(vec_a)] == [0.4082, 0.4082, 0.8165]
    assert [round(x, 4) for x in pms.geom.cross_product(vec_a, vec_b)] == [-4, -2, 3]
    assert round(pms.geom.angle(vec_a, vec_b), 4) == 37.5714
    assert round(pms.geom.angle_polar(vec_a), 4) == 0.7854
    assert round(pms.geom.angle_azi(vec_b), 4) == 0.9828
    assert round(pms.geom.angle_azi([0, 0, 0]), 4) == 1.5708   # zero-length → π/2
    assert [round(x, 4) for x in pms.geom.main_axis(1)] == [1, 0, 0]
    assert [round(x, 4) for x in pms.geom.main_axis(2)] == [0, 1, 0]
    assert [round(x, 4) for x in pms.geom.main_axis(3)] == [0, 0, 1]
    assert [round(x, 4) for x in pms.geom.main_axis("x")] == [1, 0, 0]
    assert [round(x, 4) for x in pms.geom.main_axis("y")] == [0, 1, 0]
    assert [round(x, 4) for x in pms.geom.main_axis("z")] == [0, 0, 1]
    assert pms.geom.main_axis("h") == "Wrong axis definition..."
    assert pms.geom.main_axis(100) == "Wrong axis definition..."
    assert pms.geom.main_axis(0.1) == "Wrong axis definition..."
    assert [round(x, 4) for x in pms.geom.rotate(vec_a, "x", 90, True)] == [1.0, -2.0, 1.0]
    assert pms.geom.rotate(vec_a, [0, 1, 2, 3], 90, True) is None
    assert pms.geom.rotate(vec_a, "h", 90, True) is None

    # Return types: all scalar functions return plain Python types (not numpy)
    import math as _math
    assert isinstance(pms.geom.dot_product(vec_a, vec_b), (int, float))
    assert isinstance(pms.geom.length(vec_a), float)
    assert isinstance(pms.geom.vector(vec_a, vec_b), list)
    assert isinstance(pms.geom.unit(vec_a), list)
    assert isinstance(pms.geom.cross_product(vec_a, vec_b), list)
    assert isinstance(pms.geom.angle(vec_a, vec_b), float)
    # rotate: list for single vector, ndarray for array data
    import numpy as _np
    assert isinstance(pms.geom.rotate(vec_a, "z", 0, True), list)
    assert isinstance(pms.geom.rotate(_np.tile(vec_a, (4, 1)).T, "z", 0, True), _np.ndarray)

    # Dimension-generic: length and vector work for 2D inputs
    assert round(pms.geom.length([3, 4]), 4) == 5.0
    assert pms.geom.vector([1, 2], [4, 6]) == [3, 4]


def test_database():
    assert pms.db.get_mass("H") == 1.0079
    assert pms.db.get_mass("DOTA") is None


def test_atom():
    atom = pms.Atom("O", "O", 5)

    atom.set_atom_type("H")
    atom.set_name("HO1")
    atom.set_residue(0)

    assert atom.get_atom_type() == "H"
    assert atom.get_name() == "HO1"
    assert atom.get_residue() == 0
    assert atom.__str__() == "   Residue Name Type\n0        0  HO1    H"


def test_molecule_loading():
    mol_gro = pms.Molecule(inp="data/benzene.gro")
    mol_pdb = pms.Molecule(inp="data/benzene.pdb")
    mol_mol2 = pms.Molecule(inp="data/benzene.mol2")

    mol_copy = pms.Molecule(inp=[mol_mol2])
    mol_concat = pms.Molecule(inp=[mol_gro, mol_pdb])

    mol_append = pms.Molecule(inp="data/benzene.gro")
    mol_append.append(mol_gro)

    pos_gro = [[round(x, 4) for x in col] for col in mol_gro.column_pos()]
    pos_pdb = [[round(x, 4) for x in col] for col in mol_pdb.column_pos()]
    pos_mol2 = [[round(x, 4) for x in col] for col in mol_mol2.column_pos()]
    pos_copy = [[round(x, 4) for x in col] for col in mol_copy.column_pos()]
    pos_concat = [[round(x, 4) for x in col] for col in mol_concat.column_pos()]
    pos_append = [[round(x, 4) for x in col] for col in mol_append.column_pos()]

    assert pos_gro == pos_pdb
    assert pos_gro == pos_mol2
    assert pos_gro == pos_copy
    assert [col + col for col in pos_gro] == pos_concat
    assert [col + col for col in pos_gro] == pos_append
    assert pms.Molecule(inp="data/benzene.DOTA").get_num() == 0


def test_molecule_properties():
    mol = pms.Molecule(inp="data/benzene.gro")

    assert mol.pos(0) == [0.0935, 0.0000, 0.3143]
    assert [round(x, 4) for x in mol.bond(0, 1)] == [0.1191, 0.0, 0.0687]
    assert list(mol.bond([1, 0, 0], [0, 0, 0])) == [-1, 0, 0]
    assert mol.get_box() == [0.4252, 0.001, 0.491]
    assert [round(x, 4) for x in mol.centroid()] == [0.2126, 0.0, 0.2455]
    assert [round(x, 4) for x in mol.com()] == [0.2126, 0.0, 0.2455]


def test_molecule_editing():
    mol = pms.Molecule(inp="data/benzene.gro")

    mol.translate([0, 0.1, 0.2])
    assert [round(x, 4) for x in mol.pos(3)] == [0.3317, 0.1000, 0.3768]
    mol.rotate("x", 45)
    assert [round(x, 4) for x in mol.pos(3)] == [0.3317, -0.1957, 0.3371]
    mol.move(0, [1, 1, 1])
    assert [round(x, 4) for x in mol.pos(3)] == [1.2382, 1.0972, 0.9028]
    mol.zero()
    assert [round(x, 4) for x in mol.pos(3)] == [0.3317, 0.2222, 0.1250]
    mol.put(3, [0, 0, 0])
    assert [round(x, 4) for x in mol.pos(3)] == [0.0000, 0.0000, 0.0000]
    mol.part_move([0, 1], [2, 3, 4], 0.5)
    mol.part_move([0, 1], 1, 0.5)
    assert [round(x, 4) for x in mol.pos(3)] == [0.3140, -0.1281, 0.1281]
    mol.part_rotate([0, 1], [2, 3, 4], 45, 1)
    mol.part_rotate([0, 1], 1, 45, 1)
    assert [round(x, 4) for x in mol.pos(3)] == [-0.1277, 0.0849, -0.3176]
    mol.part_angle([0, 1], [1, 2], [1, 2, 3, 4], 45, 1)
    assert [round(x, 4) for x in mol.pos(3)] == [-0.1360, -0.1084, -0.3068]
    mol.part_angle([0, 0, 1], [0, 1, 0], 1, 45, 1)
    assert [round(x, 4) for x in mol.pos(3)] == [-0.1360, -0.1084, -0.3068]

    assert mol._vector(0.1, 0.1) is None
    assert mol._vector([0, 0], [0, 0]) is None
    assert mol.part_angle([0, 0, 1, 0], [0, 1, 0, 0], 1, 45, 1) is None
    assert mol.part_angle([0, 0], [0, 1, 2], 1, 45, 1) is None


def test_molecule_creation():
    mol = pms.Molecule()

    mol.add("C", [0, 0.1, 0.2])
    mol.add("C", 0, r=0.1, theta=90)
    mol.add("C", 1, [0, 1], r=0.1, theta=90)
    mol.add("C", 2, [0, 2], r=0.1, theta=90, phi=45)
    assert [round(x, 4) for x in mol.pos(3)] == [0.0500, 0.0500, 0.0293]
    mol.delete(2)
    assert [round(x, 4) for x in mol.pos(2)] == [0.0500, 0.0500, 0.0293]
    mol.add("C", [0, 0.1, 0.2])
    assert mol.overlap() == {0: [3]}
    mol.switch_atom_order(0, 2)
    assert [round(x, 4) for x in mol.pos(0)] == [0.0500, 0.0500, 0.0293]
    mol.set_atom_type(0, "R")
    assert mol.get_atom_list()[0].get_atom_type() == "R"
    assert mol.get_atom_type(0) == "R"
    mol.set_atom_name(0, "RuX")
    assert mol.get_atom_list()[0].get_name() == "RuX"
    mol.set_atom_residue(0, 1)
    assert mol.get_atom_list()[0].get_residue() == 1


def test_molecule_set_get():
    mol = pms.Molecule()

    mol.set_name("test_mol")
    mol.set_short("TMOL")
    mol.set_box([1, 1, 1])
    mol.set_charge(1.5)
    mol.set_masses([1, 2, 3])

    assert mol.get_name() == "test_mol"
    assert mol.get_short() == "TMOL"
    assert mol.get_box() == [1, 1, 1]
    assert mol.get_num() == 0
    assert mol.get_charge() == 1.5
    assert mol.get_masses() == [1, 2, 3]
    assert mol.get_mass() == 6


def test_molecule_representation():
    mol = pms.Molecule()
    mol.add("H", [0.0, 0.1, 0.2], name="HO1")

    assert mol.__str__() == "   Residue Name Type    x    y    z\n0        0  HO1    H  0.0  0.1  0.2"


def test_generic():
    assert [round(x, 4) for x in pms.gen.alkane(10, "decane", "DEC").pos(5)] == [0.0472, 0.1028, 0.7170]
    assert [round(x, 4) for x in pms.gen.alkane(1, "methane", "MET").pos(0)] == [0.0514, 0.0890, 0.0363]
    assert [round(x, 4) for x in pms.gen.alcohol(10, "decanol", "DCOL").pos(5)] == [0.0363, 0.1028, 0.7170]
    assert [round(x, 4) for x in pms.gen.alcohol(1, "methanol", "MEOL").pos(0)] == [0.0715, 0.0890, 0.0363]
    assert [round(x, 4) for x in pms.gen.ketone(10, 5, "decanone", "DCON").pos(5)] == [0.0472, 0.1028, 0.7170]
    assert [round(x, 4) for x in pms.gen.tms(separation=30).pos(5)] == [0.0273, 0.0472, 0.4525]
    assert [round(x, 4) for x in pms.gen.tms(is_si=False).pos(5)] == [0.0273, 0.0472, 0.4976]
    assert [round(x, 4) for x in pms.gen.silanol().pos(0)] == [0.000, 0.000, 0.000]
    assert pms.gen.ketone(2, 0) is None


def test_store():
    mol = pms.Molecule(inp="data/benzene.gro")

    mol.set_atom_residue(1, 1)

    pms.Store(mol, "output").job("store_job", "store_master.job")
    pms.Store(mol, "output").obj("store_obj.obj")
    pms.Store(mol, "output").gro("store_gro.gro", True)
    pms.Store(mol, "output").pdb("store_pdb.pdb", True)
    pms.Store(mol, "output").xyz("store_xyz.xyz")
    pms.Store(mol, "output").lmp("store_lmp.lmp")
    pms.Store(mol, "output").grid("store_grid.itp")

    pms.Store({})
    assert pms.Store(mol).top() is None
