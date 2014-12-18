"""Iterate over parts of a cif file."""



from fr3d.cif.reader import Cif
import numpy as np
import matplotlib.pyplot as plt
#read cif file and store structure
def main(filename):
    with open(filename, 'rb') as raw:
        structure = Cif(raw).structure()
    
    for residue in structure.residues(chain='A', sequence = ['G']):
        base_cent = residue.centers['base']
        for aa_residue in structure.residues(chain = 'H'):
            for atom in aa_residue.atoms():
                atom_coord = atom.coordinates()
                dist = np.linalg.norm(base_cent-atom_coord)
                if dist <= 5:
                    amino_cent = atom_coord - base_cent
                    plt.scatter(amino_cent[0],amino_cent[1])
                    plt.show()
main('E:\\Leontis\\Python scripts\\2AW7.cif')
