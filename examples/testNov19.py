"""Detect and plot RNA base- amino acid interactions."""

from fr3d.cif.reader import Cif
from fr3d.definitions import RNAbasecoordinates
#from fr3d.definitions import RNAbaseheavyatoms
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
#creates lists of rotated amino acid coordinates
rotated_aa_x = []
rotated_aa_y = []
rotated_aa_z = []
#creates lists of rotated base coordinates
new_base_x = []
new_base_y = []
new_base_z = []
#reads cif file and store structure
def main(filename):
    with open(filename, 'rb') as raw:
        structure = Cif(raw).structure()
        """All RNA bases are placed in the standard orientation. All Hydrogen
 atoms are inferred. Rotation matrix is calculated for each base."""
        structure.infer_hydrogens()
        
    # Input base of interest
    sequence = 'A'
    
    for residue in structure.residues(chain='A', sequence=sequence):
        base_center = residue.centers['base']
        for aa_residue in structure.residues(sequence = 'SER'):
            for atom in aa_residue.atoms():
                aa_center = aa_residue.centers['aa_sidechain']
                dist_vector = np.subtract(base_center, aa_center)
                dist_matrix = np.matrix(dist_vector)
                transposed_dist_vector = dist_matrix.transpose()
                dist_scalar = np.linalg.norm(dist_vector)
                if dist_scalar <= 10:
                    rotated_aa_center = residue.rotation_matrix * transposed_dist_vector
                    rotated_aa_x.append(rotated_aa_center[0])
                    rotated_aa_y.append(rotated_aa_center[1])
                    rotated_aa_z.append(rotated_aa_center[2])

    #creates a list of base atoms and stores 
    #baseatoms = RNAbaseheavyatoms[sequence]
    connections = ['N1','C6','C6','N6','C6','C5','C5','C4','C5','N7','N7','C8','C8','N9','N9','C4','C4','C5','C4','N3','N3','C2','C2','N1']  
    for atomname in connections:
        coord = RNAbasecoordinates[sequence][atomname]
        new_base_x.append(coord[0])
        new_base_y.append(coord[1])
        new_base_z.append(coord[2])
   
main('E:\\Leontis\\Python scripts\\2AW7.cif')
# 3D plots of base-aa interactions
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(rotated_aa_x,rotated_aa_y,rotated_aa_z, c= 'r', marker = 'o')
ax.scatter(0, 0, 0, c='b', marker='o')
#ax.scatter(new_base_x,new_base_y,new_base_z, c= 'b', marker = 'o')
ax.plot(new_base_x, new_base_y, new_base_z, label= 'Base')
ax.set_xlabel('X Axis')
ax.set_ylabel('Y Axis')
ax.set_zlabel('Z Axis')
plt.title('Adenosine with Ser sidechain')
plt.show()