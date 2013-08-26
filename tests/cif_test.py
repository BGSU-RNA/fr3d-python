from unittest import TestCase

from pdbx.reader.PdbxParser import PdbxReader as Reader

from fr3d.cif.reader import CIF


class CifTest(TestCase):
    def setUp(self):
        self.data = []
        with open('files/1GID.cif', 'rb') as raw:
            reader = Reader(raw)
            reader.read(self.data)
        self.data = self.data[0]

    def test_loads_a_file(self):
        self.assertIsNot(self.data, [])


class FR3DReaderTest(TestCase):
    def setUp(self):
        with open('files/1GID.cif', 'rb') as raw:
            self.reader = CIF(raw)

    def test_generates_all_unit_ids(self):
        self.fail()
