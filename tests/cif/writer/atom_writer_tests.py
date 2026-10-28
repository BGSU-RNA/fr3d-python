import io
import os
import unittest

from fr3d.cif.reader import Cif as Reader
from fr3d.cif.writer import CifAtom as Writer


class BasicTest(unittest.TestCase):
    def setUp(self):
        self.handle = io.StringIO()
        self.writer = Writer(self.handle)

        base_dir = os.path.dirname(os.path.dirname(__file__))
        candidates = [
            os.path.join(os.path.dirname(base_dir), 'cif_files', '1FAT.cif'),
            os.path.join(base_dir, 'files', '1FAT.cif'),
            os.path.join(os.path.dirname(base_dir), 'files', '1FAT.cif'),
            os.path.join('tests', 'cif_files', '1FAT.cif'),
            os.path.join('files', '1FAT.cif'),
        ]
        filename = next((path for path in candidates if os.path.exists(path)), None)
        if filename is None:
            self.skipTest('Missing CIF fixture for 1FAT')

        with open(filename, 'rt') as raw:
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
