from unittest import TestCase

from fr3d.cif.reader import CIF


class ReaderStructureTest(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.reader = CIF('files/1GID.cif')

    def setUp(self):
        self.structures = self.reader.structures()

    def test_loads_all_structures(self):
        val = len(self.structures)
        ans = 1
        self.assertEqual(ans, val)

    def test_assigns_pdb_id(self):
        val = self.structures[0].pdb
        ans = '1GID'
        self.assertEqual(ans, val)


class ReaderResidueTest(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.reader = CIF('files/1GID.cif')

    def setUp(self):
        self.residues = self.reader.structures()[0].residues()

    def test_assigns_atoms_to_residues(self):
        self.residues.sort(key=lambda r: r.number)
        val = len(self.residues[49].atoms())
        ans = 20
        self.assertEqual(ans, val)


class ReaderAtomTest(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.reader = CIF('files/1GID.cif')

    def setUp(self):
        self.atoms = []
        for residue in self.reader.structures()[0].residues():
            self.atoms.extend(residue.atoms())

    def test_loads_all_atoms(self):
        val = len(self.atoms)
        ans = 6824
        self.assertEqual(ans, val)
