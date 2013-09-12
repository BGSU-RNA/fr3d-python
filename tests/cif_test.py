from unittest import TestCase

from fr3d.cif.reader import CIF


with open('files/1GID.cif', 'rb') as raw:
    reader = CIF(raw)
    STRUCTURES = reader.structures()


class ReaderStructureTest(TestCase):
    def setUp(self):
        self.structures = STRUCTURES

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
        self.residues = STRUCTURES[0].residues()

    def test_assigns_atoms_to_residues(self):
        self.residues.sort(key=lambda r: r.number)
        val = len(self.residues[49].atoms())
        ans = 20
        self.assertEqual(ans, val)


class ReaderAtomTest(TestCase):
    def setUp(self):
        self.atoms = []
        for residue in STRUCTURES[0].residues():
            self.atoms.extend(residue.atoms())

    def test_loads_all_atoms(self):
        val = len(self.atoms)
        ans = 6824
        self.assertEqual(ans, val)
