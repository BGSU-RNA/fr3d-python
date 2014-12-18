"""Iterate over parts of a cif file."""



from fr3d.cif.reader import Cif
import numpy as np
import matplotlib.pyplot as plt
#read cif file and store structure
new_aa_x = []
new_aa_y = []

def main(filename):
    with open(filename, 'rb') as raw:
        structure = Cif(raw).structure()
    
    for residue in structure.residues(chain='A', sequence = ['G']):
        base_center = residue.centers['base']
        for aa_residue in structure.residues(sequence = 'LYS'):
                         
            for atom in aa_residue.atoms():
                aa_center = aa_residue.centers['aa_sidechain']
                dist = np.linalg.norm(base_center - aa_center)
                if dist <= 6:
                    new_aa_center = aa_center - base_center
                    new_aa_x.append(new_aa_center[0])
                    new_aa_y.append(new_aa_center[1])
                    
main('E:\\Leontis\\Python scripts\\2AW7.cif')
print(new_aa_x, new_aa_y)
plt.scatter(new_aa_x,new_aa_y)