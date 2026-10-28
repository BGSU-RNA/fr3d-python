import os
from unittest import SkipTest
from unittest import TestCase

from fr3d.cif.reader import Cif


class ReaderTest(TestCase):
    name = None

    @classmethod
    def setUpClass(cls):
        base_dir = os.path.dirname(__file__)
        candidates = [
            os.path.join(os.path.dirname(base_dir), 'cif_files', cls.name + '.cif'),
            os.path.join(base_dir, 'files', cls.name + '.cif'),
            os.path.join(os.path.dirname(base_dir), 'files', cls.name + '.cif'),
            os.path.join('tests', 'cif_files', cls.name + '.cif'),
            os.path.join('files', cls.name + '.cif'),
        ]

        filename = next((path for path in candidates if os.path.exists(path)), None)
        if filename is None:
            raise SkipTest('Missing CIF fixture for %s' % cls.name)

        with open(filename, 'rt') as raw:
            cls.cif = Cif(raw)
            cls.structure = cls.cif.structure()

    def setUp(self):
        self.cif = self.__class__.cif
        self.structure = self.__class__.structure
