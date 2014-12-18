"""Iterate over parts of a cif file."""



from fr3d.cif.reader import Cif
import numpy as np
#read cif file and store structure
def main(filename):
    with open(filename, 'rb') as raw:
        structure = Cif(raw).structure()
    
    for residue in structure.residues(chain='A', sequence = ['G']):
        base_cent = residue.centers['base']
        for aa_residue in structure.residues():
            if aa_residue.chain == 'A':
                continue
            for atom in aa_residue.atoms():
                atom_coord = np.array(atom.coordinates())
                dist = np.linalg.norm(base_cent-atom_coord)
                if dist <= 5:
                    print('  residue: %s' % aa_residue.unit_id())
main('E:\\Leontis\\Python scripts\\2AW7.cif')
