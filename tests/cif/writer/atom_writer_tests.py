import cStringIO as sio
import unittest

from fr3d.cif.reader import Cif as Reader
from fr3d.cif.writer import CifAtom as Writer


class BasicTest(unittest.TestCase):
    def setUp(self):
        self.handle = sio.StringIO()
        self.writer = Writer(self.handle)
        with open('files/1FAT.cif', 'rb') as raw:
            self.data = Reader(raw).structure()

    def write(self):
        self.writer(self.data)
        return self.handle.getvalue()

    def test_can_write_all_atoms(self):
        val = len(self.write().split("\n"))
        self.assertTrue(val > 231)

    def test_will_include_unit_ids(self):
        val = self.write().split('\n')
        line = val[-3].split()
        assert '1FAT|1|D|HOH|312' == line[-1]

    def test_will_not_include_unit_ids_if_requested(self):
        self.writer.unit_ids = False
        val = self.write().split('\n')
        line = val[-3].split()
        assert '1' == line[-1]
