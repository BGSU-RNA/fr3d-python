# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 16:58:14 2016

@author: Poorna
"""
from fr3d.cif.reader import Cif
from fr3d.classifiers.base_aafg import classification
from fr3d.classifiers.base_aafg import detect_edge
def get_structure(filename):
    with open(filename, 'rb') as raw:
        structure = Cif(raw).structure()
        """All RNA bases are placed in the standard orientation. All Hydrogen
 atoms are inferred. Rotation matrix is calculated for each base."""
        structure.infer_hydrogens()
        return structure


PDB_List = ['4YBB']
base_seq_list = ['A', 'U', 'C', 'G']
aa_list = ['ALA','VAL','ILE','LEU','ARG','LYS','HIS','ASP','GLU','ASN','GLN','THR','SER','TYR','TRP','PHE','PRO','CYS','MET']
if __name__=="__main__":
    
    for PDB in PDB_List:
        structure = get_structure('E:\\Leontis\\Python scripts\\CIF\\%s.cif' % PDB)
                
        aa_part = 'aa_fg'
        base_part = 'base'
                                               
                 
        bases = structure.residues(chain = "AA", sequence= base_seq_list)
        amino_acids = structure.residues(sequence=aa_list)
        
        for base in bases:
            for aa in amino_acids:
                interaction = classification(base, aa)
                if interaction == "pseudopair":
                    edge = detect_edge(base, aa)
                else:
                    edge == "N/A"
                    
                print base, aa, interaction, edge