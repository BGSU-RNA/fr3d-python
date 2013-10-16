import numpy
from fr3d.geometry.superpositions import besttransformation
from fr3d.data import Atom
from fr3d.definitions import RNAbaseheavyatoms
from fr3d.definitions import RNAbasehydrogens
from fr3d.definitions import RNAbasecoordinates


def infer_hydrogen_coordinates(residue, heavyatoms, basehydrogens,
                               basecoordinates):
    R = []
    S = []
    baseheavy = heavyatoms[residue.sequence]

    for atom in residue.atoms(name=baseheavy):
        coordinates = atom.coordinates()
        R.append(coordinates)
        S.append(basecoordinates[residue.sequence][atom.name])

    R = numpy.array(R)
    R = R.astype(numpy.float)
    S = numpy.array(S)
    rotation_matrix, fitted, base_center, rmsd = besttransformation(R, S)
    hydrogens = basehydrogens[residue.sequence]

    for hydrogenatom in hydrogens:
        hydrogencoordinates = basecoordinates[residue.sequence][hydrogenatom]
        newcoordinates = base_center + \
            numpy.dot(hydrogencoordinates, numpy.transpose(rotation_matrix))
        residue._atoms.append(Atom({'name': hydrogenatom,
                                    'x': newcoordinates[0, 0],
                                    'y': newcoordinates[0, 1],
                                    'z': newcoordinates[0, 2]}))


def infer_all_hydrogen_coordinates(struct):

    for residue in struct.residues(sequence=['A', 'C', 'G', 'U']):
        infer_hydrogen_coordinates(residue, RNAbaseheavyatoms,
                                   RNAbasehydrogens, RNAbasecoordinates)
