from unittest import TestCase

from pdbx.reader.PdbxParser import PdbxReader as Reader


class CifTest(TestCase):
    def setUp(self):
        self.data = []
        with open('files/1GID.cif', 'rb') as raw:
            reader = Reader(raw)
            reader.read(self.data)
        self.data = self.data[0]

    def test_loads_a_file(self):
        self.assertIsNot(self.data, [])
