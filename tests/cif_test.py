from unittest import TestCase

from fr3d.cif.reader import CIF

READER = CIF("files/1GID.cif")


class ReaderStructureTest(TestCase):

    def setUp(self):
        self.structures = READER.structures()

    def test_loads_all_structures(self):
        val = len(self.structures)
        ans = 1
        self.assertEqual(ans, val)

    def test_assigns_pdb_id(self):
        val = self.structures[0].pdb
        ans = '1GID'
        self.assertEqual(ans, val)


class ReaderResidueTest(TestCase):

    def setUp(self):
        self.residues = list(READER.structures()[0].residues())

    def test_assigns_atoms_to_residues(self):
        self.residues.sort(key=lambda r: r.number)
        val = len(self.residues[49].atoms())
        ans = 20
        self.assertEqual(ans, val)

    def test_assigns_numbers_correctly(self):
        self.residues.sort(key=lambda r: r.number)
        val = self.residues[0].number
        ans = '103'
        self.assertEqual(ans, val)

    def test_assigns_pdb(self):
        val = self.residues[0].pdb
        ans = '1GID'
        self.assertEqual(ans, val)

    def test_assigns_model(self):
        val = self.residues[0].model
        ans = '1'
        self.assertEqual(ans, val)

    def test_assigns_chain(self):
        self.residues.sort(key=lambda r: '%s%s' % (r.chain, r.number))
        val = self.residues[0].chain
        ans = 'A'
        self.assertEqual(ans, val)

    def test_assigns_symmetry(self):
        val = self.residues[0].symmetry
        ans = '1_555'
        self.assertEqual(ans, val)

    def test_assigns_ins_code(self):
        val = self.residues[0].ins_code
        ans = '?'
        self.assertEqual(ans, val)

    def test_assigns_sequence(self):
        self.residues.sort(key=lambda r: r.number)
        val = self.residues[0].sequence
        ans = 'G'
        self.assertEqual(ans, val)

    def test_can_generate_unit_id(self):
        self.residues.sort(key=lambda r: '%s%s' % (r.chain, r.number))
        val = self.residues[0].unit_id()
        ans = '1GID|1|A|G|103'
        self.assertEqual(ans, val)


class ReaderAtomTest(TestCase):

    def setUp(self):
        self.atoms = []
        for residue in READER.structures()[0].residues():
            self.atoms.extend(residue.atoms())

    def test_loads_all_atoms(self):
        val = len(self.atoms)
        ans = 6824
        self.assertEqual(ans, val)
