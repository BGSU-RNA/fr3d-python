from unittest import TestCase
import numpy
from numpy import array
from numpy.testing import assert_almost_equal
from fr3d.cif.reader import CIF
from fr3d.definitions import *

from fr3d.geometry.inferhydrogen import infer_all_hydrogen_coordinates
from fr3d.geometry.inferhydrogen import infer_hydrogen_coordinates


class InferhydrogenTest(TestCase):

    def setUp(self):  
        with open('files/1GID.cif', 'rb') as raw:
            reader = CIF(raw)
            self.structures = reader.structures()[0]
            
    def test_hydrogen_infers_on_residue(self):
        heavyatoms = RNAbaseheavyatoms
        basehydrogens = RNAbasehydrogens
        basecoordinates = RNAbasecoordinates 
        
        residue = self.structures.residues(sequence=['A'])[0]
        infer_hydrogen_coordinates(residue, heavyatoms, basehydrogens, 
                                   basecoordinates)
        
        atoms = list(residue.atoms(name=['H2', 'H8', 'H9', '1H6', '2H6']))
        
        self.assertEqual(len(atoms), 5)


