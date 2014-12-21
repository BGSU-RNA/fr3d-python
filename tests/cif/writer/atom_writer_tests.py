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
