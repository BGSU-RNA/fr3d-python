"""Detect and plot RNA base- amino acid interactions."""



from fr3d.cif.reader import Cif
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from fr3d.definitions import RNAbasecoordinates
from fr3d.definitions import RNAbaseheavyatoms
#read cif file and store structure
new_aa_x = []
new_aa_y = []
new_aa_z = []
#creates lists of rotated base coordinates
new_base_x = []
new_base_y = []
new_base_z = []
def main(filename):
    with open(filename, 'rb') as raw:
        structure = Cif(raw).structure()
        """All RNA bases are placed in the standard orientation. All Hydrogen
 atoms are inferred. Rotation matrix is calculated for each base 
 and is appended as an attribute."""
        structure.infer_hydrogens()
         # Input base of interest
    sequence = 'G'
    
    for residue in structure.residues(chain='A', sequence = sequence):
        base_center = residue.centers['base']
        for aa_residue in structure.residues(sequence = 'GLN'):
            aa_sidechain_center = aa_residue.centers['aa_sidechain']
            c2c_vector = np.subtract(base_center, aa_sidechain_center)
            c2c_matrix = np.matrix(c2c_vector)
            c2c_matrix = c2c_matrix.transpose()
            c2c_dist = np.linalg.norm(c2c_vector)
            if c2c_dist <= 10:
                new_aa_center = residue.rotation_matrix * c2c_matrix
                new_aa_x.append(new_aa_center[0])
                new_aa_y.append(new_aa_center[1])
                new_aa_z.append(new_aa_center[2])
    #creates a list of base atoms and stores 
    baseatoms = RNAbaseheavyatoms[sequence]
    for atomname in baseatoms:
        coord = RNAbasecoordinates[sequence][atomname]
        new_base_x.append(coord[0])
        new_base_y.append(coord[1])
        new_base_z.append(coord[2])

main('E:\\Leontis\\Python scripts\\2AW7.cif')
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(new_aa_x,new_aa_y,new_aa_z, c= 'r', marker = 'o')
ax.scatter(new_base_x, new_base_y, new_base_z, label= 'Base')
ax.scatter(0, 0, 0, c='b', marker='o')
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')
plt.title('C with Tyr sidechain')
plt.show()