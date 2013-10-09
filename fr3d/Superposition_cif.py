import numpy
import os

# We need to figure out how to set the working directory within Spyder so
# we don't need to do this, which is different on each computer!
os.chdir('C:/Users/Kevin/Documents/GitHub/fr3d-python')
#os.chdir('C:/Users/zirbel/Documents/GitHub/fr3d-python')

from fr3d.cif.reader import CIF
from fr3d.geometry.superpositions import besttransformation
from fr3d.data import Atom

with open('files/1GID.cif', 'rb') as raw:
    reader = CIF(raw)
    STRUCTURES = reader.structures()

from fr3d.definitions import *
struct = STRUCTURES[0] 

# (above: just renaming for typing convenience, recall this is one cif file)

def infer_all_hydrogen_coordinates(struct, heavyatoms, basehydrogens, 
                                   basecoordinates):
                                       
    for residue in struct.residues(sequence=['A', 'C', 'G', 'U']):
        R = []
        S = []
        baseheavy = heavyatoms[residue.sequence]
    
        for atom in residue.atoms(name=baseheavy):
            coordinates = atom.coordinates()
            R.append(coordinates)
            S.append(basecoordinates[residue.sequence][atom.name])
    
        R = numpy.array(R)
        R = R.astype(np.float)
        S = numpy.array(S)
        rotation_matrix, fitted, base_center, rmsd = besttransformation(R, S)
        hydrogens = basehydrogens[residue.sequence]
    
        for hydrogenatom in hydrogens:
            hydrogencoordinates = basecoordinates[residue.sequence][hydrogenatom]
            newcoordinates = base_center + \
                numpy.dot(hydrogencoordinates, numpy.transpose(rotation_matrix))
            residue._atoms.append(Atom({'name': hydrogenatom,
                                    'x': newcoordinates[0,0], 
                                    'y': newcoordinates[0,1], 
                                    'z': newcoordinates[0,2]}))
    return struct

def infer_hydrogen_coordinates(residue, heavyatoms, basehydrogens, 
                               basecoordinates):
    R = []
    S = []
    baseheavy = heavyatoms[residue.sequence]  # imports the 
    
    for atom in residue.atoms(name=baseheavy):
        coordinates = atom.coordinates()
        R.append(coordinates)
        S.append(basecoordinates[residue.sequence][atom.name])
    
    R = numpy.array(R)
    R = R.astype(np.float)
    S = numpy.array(S)
    rotation_matrix, fitted, base_center, rmsd = besttransformation(R, S)
    hydrogens = basehydrogens[residue.sequence]
    
    for hydrogenatom in hydrogens:
        hydrogencoordinates = basecoordinates[residue.sequence][hydrogenatom]
        newcoordinates = base_center + \
            numpy.dot(hydrogencoordinates, numpy.transpose(rotation_matrix))
        residue._atoms.append(Atom({'name': hydrogenatom,
                                    'x': newcoordinates[0,0], 
                                    'y': newcoordinates[0,1], 
                                    'z': newcoordinates[0,2]}))
    
    return residue    
