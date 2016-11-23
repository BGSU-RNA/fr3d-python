# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 16:58:14 2016

@author: Poorna
"""
import os
import sys

here = os.path.abspath(os.path.dirname(__file__))
fr3d = os.path.abspath(os.path.join(here, ".."))
sys.path.insert(0, fr3d)

from fr3d.cif.reader import Cif
from fr3d.classifiers.base_aafg import Classifier

def get_structure(filename):
    with open(filename, 'rb') as raw:
        structure = Cif(raw).structure()
        """All RNA bases are placed in the standard orientation. All Hydrogen
 atoms are inferred. Rotation matrix is calculated for each base."""
        structure.infer_hydrogens()
        return structure


PDB_List = ['1FJG']

if __name__=="__main__":
    
    for PDB in PDB_List:
        structure = get_structure('E:\\Leontis\\Python scripts\\CIF\\%s.cif' % PDB)
        classifier = Classifier()
        print classifier.classify(structure)
        