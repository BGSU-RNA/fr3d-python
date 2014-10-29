import os
from unittest import TestCase

from fr3d.cif.reader import CIF


class ReaderTest(TestCase):
    name = None

    @classmethod
    def setUpClass(cls):
        with open(os.path.join('files', cls.name + '.cif'), 'rb') as raw:
            cls.cif = CIF(raw)
            cls.structure = cls.cif.structure()

    def setUp(self):
        self.cif = self.__class__.cif
        self.structure = self.__class__.structure
