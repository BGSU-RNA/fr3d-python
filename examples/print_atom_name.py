# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 15:08:09 2014

@author: Poorna
"""

from fr3d.cif.reader import Cif

def main(filename):
    with open(filename, 'rb') as raw:
        structure = Cif(raw).structure()
        print('Iterating over atoms')
    for residue in structure.residues(chain='A', sequence = 'C'):
        for atom in residue.atoms():
            print residue.unit_id()            
            print(atom.name)
main('E:\\Leontis\\Python scripts\\CIF\\2AW7.cif')
        