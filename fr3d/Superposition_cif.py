import numpy
import os

# We need to figure out how to set the working directory within Spyder so
# we don't need to do this, which is different on each computer!
#os.chdir('C:/Users/Kevin/Documents/GitHub/fr3d-python')
os.chdir('C:/Users/zirbel/Documents/GitHub/fr3d-python')

from fr3d.cif.reader import CIF
from fr3d.geometry.superpositions import besttransformation

with open('files/1GID.cif', 'rb') as raw:
    reader = CIF(raw)
    STRUCTURES = reader.structures()

from fr3d.definitions import *
struct = STRUCTURES[0] 
# (above: just renaming for typing convenience, recall this is one cif file)

for residue in struct.residues(sequence=['A', 'C', 'G', 'U']):
# above begins loop for all residues which is classified by base
# as with either A, C, G, or U. (DNA would be DA, DC, DG, DT per PDB standard)
    R = [] # empty list to be filled w/ [x,y,z] of this residue's
           # heavy atoms in crystal structure of 1GIF.cif
    S = [] # empty list to be filled w/ [x,y,z] of this residue's
           # heavy atoms in standard base
    baseheavy = RNAbaseheavyatoms[residue.sequence]  # imports the 
           # base of the residue of interest see:
           # struct.residues(sequence=['A','C','U','G'])[0].sequence
           # Above is the base of the first entry in struct.residues(sequence=['A','C','U','G'])
           # This exports ...
           # RNAbaseheavyatoms is a dictionary in definitions.py
           # Each entry in this dictionary is a list,
           # a list of the heavy atoms in the RNA base.
           # RNAbaseheavyatoms[struct.residues(sequence=['A','C','U','G'])[0].sequence]

    #Now I must get the coordinates of the heavy atoms within each residue
    #in the crystal structure of 1GID.cif
    for atom in residue.atoms(name=baseheavy):  
        #Now pertaining to each heavy atom, we will record [x,y,z] and
        #add it to the list R for crystal structure, and S for standard
        coordinates = atom.coordinates()
        #In coordinates, the [x,y,z] coordinate of the heavy atom in current
        #loop status is stored.         
        R.append(coordinates)
        # R now contains this [x,y,z] coordinate as an element in its list.        
        S.append(RNAbasecoordinates[residue.sequence][atom.name])
        # Now, from definitions find the standard [x,y,z] coordinate 
        # by stating that we are considering the specific residue's 
        # base, residue.sequence,  and the atom.name which is the name for the 
        # position of the base within a molecular structure.
        # Thus, each nth element R is in correspondence to the nth element in S

    
    # Now, R and S are filled, and are what we wanted.  Thus, we can 
    # superimpose the coordinates of R onto S.  Recall R and S are a 
    # collection of the [x,y,z] coordinates of heavy atoms 
    # for a residue of the crystal, and standard, respectively.  
    # R will be, hypothetically, a transformation of the standard's coordinates
    # in a three dimensional space.  Thus, we want the crystal to be mapped to 
    # the standard.
    R=numpy.array(R)
    R = R.astype(np.float)
    S=numpy.array(S)
#    R=R.tolist()
#    S=S.tolist()
    rotation_matrix, fitted, base_center, rmsd = besttransformation(R, S)

    print residue.atoms(name='P')[0].pdb, residue.atoms(name='P')[0].chain, residue.sequence, int(residue.number)+102
    print base_center
    print rotation_matrix
    print rmsd
    