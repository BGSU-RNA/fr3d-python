# -*- coding: utf-8 -*-
"""
    This program reads one or more CIF files and produces annotations
    the glycosidic bond orientation.

"""

import numpy as np
from fr3d.modified_parent_mapping import modified_nucleotides

from fr3d.cif.reader import Cif
from fr3d.definitions import RNAconnections
from fr3d.definitions import NAbaseheavyatoms
from fr3d.definitions import NAbasehydrogens
from fr3d.definitions import nt_sugar
from fr3d.definitions import nt_phosphate
from fr3d.definitions import Ribophos_connect
from fr3d.definitions import aa_connections
from fr3d.definitions import aa_backconnect
from fr3d.definitions import aa_hydrogen_connections
from fr3d.definitions import aa_fg
from fr3d.definitions import aa_linker
from fr3d.definitions import aa_backbone
from fr3d.definitions import tilt_cutoff
from fr3d.definitions import planar_atoms
from fr3d.definitions import HB_donors
from fr3d.definitions import HB_weak_donors
from fr3d.definitions import HB_acceptors

from discrepancy import matrix_discrepancy
import csv
import urllib
import pickle
import math
import sys

import matplotlib.pyplot as plt
from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D
# note that fr3d.localpath does not synchronize with Git, so you can change it locally to point to your own directory structure
from fr3d.localpath import outputText
from fr3d.localpath import outputNAPairwiseInteractions
from fr3d.localpath import contact_list_file
from fr3d.localpath import inputPath
from fr3d.localpath import outputHTML
from fr3d.data.base import EntitySelector

from fr3d.ordering.greedyInsertion import orderWithPathLengthFromDistanceMatrix

#from fr3d.classifiers.base_aafg import distance_metrics
from datetime import datetime
from math import floor
import os
from os import path

from time import time

from class_limits import nt_nt_cutoffs

HB_donor_hydrogens = {}
HB_donor_hydrogens['A'] = {"N6":["1H6","2H6"], "C2":["H2"], "C8":["H8"], "O2'":[]}
HB_donor_hydrogens['G'] = {"N1":["H1"], "N2":["2H2","1H2"], "C8":["H8"], "O2'":[]}
HB_donor_hydrogens['C'] = {"N4":["1H4","2H4"], "C5":["H5"], "C6":["H6"], "O2'":[]}
HB_donor_hydrogens['U'] = {"N3":["H3"], "C5":["H5"], "C6":["H6"], "O2'":[]}

def myTimer(state,data={}):

    # add elapsed time to the current state of the timer
    if "currentState" in data:
        currentState = data["currentState"]
        data[currentState] += time() - data["lastTime"]

    if state == "summary":
        print("Summary of time taken:")
        for state in data["allStates"]:
            if not state == "lastTime" and not state == "currentState":
                print("%-31s: %10.3f seconds %10.3f minutes" % (state,data[state],data[state]/60))
    elif not state in data:
        data[state] = 0
        # keep track of states and the order in which they were seen
        if "allStates" in data:
            data["allStates"].append(state)
        else:
            data["allStates"] = [state]

    # change to the state just starting now
    data["currentState"] = state
    data["lastTime"] = time()

    return data

def get_structure(filename,PDB):

    if ".pdb" in filename:
        filename = filename.replace(".cif","")

    # if not available locally, download from PDB
    if not os.path.exists(filename):
        print("  Downloading %s from https://files.rcsb.org/download/%s.cif" % (PDB,PDB))
        if sys.version_info[0] < 3:
            urllib.urlretrieve("http://files.rcsb.org/download/%s.cif" % PDB, filename)  # python 2
        else:
            urllib.request.urlretrieve("http://files.rcsb.org/download/%s.cif" % PDB, filename)  # python 3

    with open(filename, 'rb') as raw:
        print("  Loading " + filename)
        structure = Cif(raw).structure()
        """All RNA bases are placed in the standard orientation.
        Rotation matrix is calculated for each base."""

        structure.infer_hydrogens()  # add hydrogens to NA bases and amino acids; slow

        return structure

def build_atom_to_unit_part_list():

    atom_to_part_list = defaultdict(lambda: "unknown")

    for base in NAbaseheavyatoms.keys():
        for atom in NAbaseheavyatoms[base]:
            atom_to_part_list[(base,atom)] = "base"
        for atom in NAbasehydrogens[base]:
            atom_to_part_list[(base,atom)] = "base"
        for atom in nt_phosphate[base]:
            atom_to_part_list[(base,atom)] = "nt_phosphate"
        for atom in nt_sugar[base]:
            atom_to_part_list[(base,atom)] = "nt_sugar"

    for aa in aa_backbone.keys():
        for atom in aa_backbone[aa]:
            atom_to_part_list[(aa,atom)] = "aa_backbone"
        for atom in aa_linker[aa]:
            atom_to_part_list[(aa,atom)] = "aa_linker"
        for atom in aa_fg[aa]:
            atom_to_part_list[(aa,atom)] = "aa_fg"

#    print(atom_to_part_list)

    return atom_to_part_list


def find_atom_atom_contacts(bases,amino_acids,atom_atom_min_distance):
    """Find all atoms within atom_atom_min_distance of each other,
    one from a nt, one from an aa, and report these"""
    # build a set of cubes and record which bases are in which cube
    # also record which other cubes are neighbors of each cube

    atom_atom_min_distance_squared = atom_atom_min_distance**2
    contact_list = []

    atom_to_part_list = build_atom_to_unit_part_list()

    # build a set of cubes and record which nt atoms are in which cube
    ntAtomCubeList = {}
    ntAtomCubeNeighbors = {}
    for base in bases:
        for atom in base.atoms():
            center = atom.coordinates()
            if len(center) == 3:
                x = floor(center[0]/atom_atom_min_distance)
                y = floor(center[1]/atom_atom_min_distance)
                z = floor(center[2]/atom_atom_min_distance)
                key = "%d,%d,%d" % (x,y,z)
                entry = (center[0],center[1],center[2],atom.name,base.unit_id(),atom_to_part_list[(base.sequence,atom.name)])
                if key in ntAtomCubeList:
                    ntAtomCubeList[key].append(entry)
                else:
                    ntAtomCubeList[key] = [entry]
                    ntAtomCubeNeighbors[key] = []
                    for a in [-1,0,1]:
                        for b in [-1,0,1]:
                            for c in [-1,0,1]:
                                k = "%d,%d,%d" % (x+a,y+b,z+c)
                                ntAtomCubeNeighbors[key].append(k)

    # build a set of cubes and record which amino acids are in which cube
    aaAtomCubeList = {}
    for aa in amino_acids:
        for atom in aa.atoms():
            center = atom.coordinates()
            if len(center) == 3:
                x = floor(center[0]/atom_atom_min_distance)
                y = floor(center[1]/atom_atom_min_distance)
                z = floor(center[2]/atom_atom_min_distance)
                key = "%d,%d,%d" % (x,y,z)
                entry = (center[0],center[1],center[2],atom.name,aa.unit_id(),atom_to_part_list[(aa.sequence,atom.name)])
                if key in aaAtomCubeList:
                    aaAtomCubeList[key].append(entry)
                else:
                    aaAtomCubeList[key] = [entry]

    # loop over nt atoms and aa atoms and find those within atom_atom_min_distance
    for key in ntAtomCubeList:                           # one nt atom cube
        for aakey in ntAtomCubeNeighbors[key]:           # neighboring cube
            if aakey in aaAtomCubeList:                  # if an aa atom lies in the neighbor cube
                for ntAtom in ntAtomCubeList[key]:       # loop over nt atoms in first cube
                    for aaAtom in aaAtomCubeList[aakey]: # and over aa atoms in neighbor cube
                        x = abs(ntAtom[0] - aaAtom[0])   # coordinates are stored in 0, 1, 2 of entry
                        y = abs(ntAtom[1] - aaAtom[1])   # coordinates are stored in 0, 1, 2 of entry
                        z = abs(ntAtom[2] - aaAtom[2])   # coordinates are stored in 0, 1, 2 of entry

                        if x > atom_atom_min_distance or \
                           y > atom_atom_min_distance or \
                           x*x + y*y + z*z > atom_atom_min_distance_squared:
                            continue

                        distance = math.sqrt(x*x + y*y + z*z)

                        contact_list.append("%s\t%s\t%s\t%s\t%s\t%s\t%8.4f\n" % (ntAtom[4],ntAtom[5],ntAtom[3],aaAtom[4],aaAtom[5],aaAtom[3],distance))

    return contact_list

def make_nt_cubes(bases, screen_distance_cutoff, nt_reference="base"):
    """Builds cubes with side length screen_distance_cutoff
    using nt_reference as the point for each nucleotide.
    """

    # build a set of cubes and record which bases are in which cube
    # also record which other cubes are neighbors of each cube
    baseCubeList = {}
    baseCubeNeighbors = {}

    # build a set of cubes and record which bases are in which cube
    for base in bases:
        center = base.centers[nt_reference]  # chosen reference point
        if len(center) == 3:
            x = floor(center[0]/screen_distance_cutoff)
            y = floor(center[1]/screen_distance_cutoff)
            z = floor(center[2]/screen_distance_cutoff)
            key = "%d,%d,%d" % (x,y,z)
            if key in baseCubeList:
                baseCubeList[key].append(base)
            else:
                baseCubeList[key] = [base]
                baseCubeNeighbors[key] = []
                for a in [-1,0,1]:
                    for b in [-1,0,1]:
                        for c in [-1,0,1]:
                            k = "%d,%d,%d" % (x+a,y+b,z+c)
                            baseCubeNeighbors[key].append(k)

    return baseCubeList, baseCubeNeighbors

# produce a list of interacting atoms marked by their group
def get_interacting_atoms(residue1,residue2):

    interacting_atoms = defaultdict(list)

    atom_atom_min_distance_squared = atom_atom_min_distance**2
    nt_parts = ['base','nt_sugar','nt_phosphate']
    aa_parts = ['aa_fg','aa_backbone','aa_linker']
    nt_parts = ['base']
    aa_parts = ['aa_fg']
    nt_atoms = []
    aa_atoms = []

    for nt_part in nt_parts:
        if nt_part == "base" and residue1.sequence in NAbaseheavyatoms:
            nt_atoms += NAbaseheavyatoms[residue1.sequence]
            if residue1.sequence in NAbasehydrogens:
                nt_atoms += NAbasehydrogens[residue1.sequence]
        elif nt_part == "nt_sugar" and residue1.sequence in nt_sugar:
            nt_atoms += nt_sugar[residue1.sequence]
        elif residue1.sequence in nt_phosphate:
            nt_atoms += nt_phosphate[residue1.sequence]
        else:
            continue

        for aa_part in aa_parts:
            if aa_part == "aa_fg" and residue2.sequence in aa_fg:
                aa_atoms += aa_fg[residue2.sequence]
            elif aa_part == "aa_backbone" and residue2.sequence in aa_backbone:
                aa_atoms += aa_backbone[residue2.sequence]
            elif residue2.sequence in aa_linker:
                aa_atoms += aa_linker[residue2.sequence]
            else:
                continue

            for nt_atom in residue1.atoms(name=nt_atoms):
                nt = nt_atom.coordinates()
                for aa_atom in residue2.atoms(name=aa_atoms):
                    aa = aa_atom.coordinates()
                    x = abs(nt[0]-aa[0])
                    y = abs(nt[1]-aa[1])
                    z = abs(nt[2]-aa[2])

                    if x <= atom_atom_min_distance and \
                       y <= atom_atom_min_distance and \
                       x**2 + y**2 + z**2 <= atom_atom_min_distance_squared:
                        distance = math.sqrt(x**2 + y**2 + z**2)
                        interacting_atoms[(nt_part,aa_part)].append((nt_atom,aa_atom,distance))

#    if len(interacting_atoms) > 0:
#        print(interacting_atoms)

    return interacting_atoms

def reverse_edges(interaction):

    if len(interaction) == 3:
        return interaction[0] + interaction[2] + interaction[1]
    elif len(interaction) == 4:
        return interaction[0] + interaction[1] + interaction[3] + interaction[2]
    else:
        return None

def annotate_nt_nt_interactions(bases, screen_distance_cutoff, baseCubeList, baseCubeNeighbors):

    # loop through nt cubes, loop through neighboring nt cubes,
    # then loop through bases in the two cubes,
    # screening distances between them, then annotating interactions

    count_pair = 0

    pair_to_bp_type = {}

    list_base_coord = []
    list_nt_nt = []
    contact_list = []
    output = []

    max_screen_distance = 0     # record the largest screening distance for which an interaction is found

    for nt1key in baseCubeList:                                # key to first cube
        for nt2key in baseCubeNeighbors[nt1key]:               # key to each potential neighboring cube, including the first
            if nt2key in baseCubeList:                         # if this cube was actually made
                for nt1 in baseCubeList[nt1key]:               # first nt of a potential pair
                    parent1 = get_parent(nt1.sequence)
                    gly1 = get_glycosidic_atom_coordinates(nt1,parent1)

                    if len(nt1.centers["base"]) < 3:
                        print("Missing base center for %s" % nt1.unit_id())
                        print(nt1.centers["base"])
                        continue
                    for nt2 in baseCubeList[nt2key]:           # second nt of a potential pair
                        if len(nt2.centers["base"]) < 3:
                            print("Missing base center for %s" % nt2.unit_id())
                            print(nt2.centers["base"])
                            continue
                        displacement = abs(nt2.centers["base"]-nt1.centers["base"]) # center-center

                        if displacement[0] > screen_distance_cutoff or \
                           displacement[1] > screen_distance_cutoff:
                            continue

                        screen_distance = np.linalg.norm(displacement)

                        if screen_distance > screen_distance_cutoff:
                            continue

                        if screen_distance < 2:       # some structures have overlapping nucleotides
                            continue

#                        print("  Checking for an interaction between %-18s and %-18s center-center distance %7.4f" % (nt1.unit_id(),nt2.unit_id(),screen_distance))

                        parent2 = get_parent(nt2.sequence)
                        parent_pair = parent1 + "," + parent2

#                       Note that AA, CC, GG, UU are being checked twice already!
#                       need to add T and have a plan for DNA nucleotides as well
                        if parent_pair in ['A,A','A,C','A,G','A,U','C,C','G,C','C,U','G,G','G,U','U,U']:
#                        if parent_pair in ['A,U']:
                            cutoffs = nt_nt_cutoffs[parent1+","+parent2]

                            gly2 = get_glycosidic_atom_coordinates(nt2,parent2)
                            glycosidic_displacement = np.subtract(gly2,gly1)

                            pair_data = check_coplanar(nt1,nt2,glycosidic_displacement)

                            base_ribose_stack = check_base_ribose_stack(nt1,nt2,pair_data)

                            interaction = check_basepair_cutoffs(nt1,nt2,pair_data,cutoffs)
                            if False and len(interaction) > 1:
                                print("  Identified parents as %s and %s" % (parent1,parent2))
                                print("  Found %s interaction between %-18s and %-18s" % (interaction,nt1.unit_id(),nt2.unit_id()))
                                print("  Gap value %0.8f" % pair_data["gap12"])
                                print("  Coplanar value %0.8f" % pair_data["coplanar_value"])
                                #print(pair_data)
                                print("  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt1.unit_id(),nt2.unit_id()))
                                print("")

                            if len(interaction) > 0:
                                count_pair += 1

                                inter = interaction[0][0]

                                pair_to_bp_type[(nt1.unit_id(),nt2.unit_id())] = inter
                                pair_to_bp_type[(nt2.unit_id(),nt1.unit_id())] = reverse_edges(inter)


    #print(pair_to_bp_type)

    print("  Found %d nucleotide-amino acid pairs" % count_pair)
    print("  Recorded %d nucleotide-amino acid pairs" % len(list_nt_nt))
    print("  Maximum screen distance for actual contacts is %8.4f" % max_screen_distance)

    return pair_to_bp_type, list_nt_nt, list_base_coord

def get_parent(sequence):
    """ Look up parent sequence for RNA, DNA, and modified nucleotides
    """

    if sequence in ['A','C','G','U']:
        return sequence
    elif sequence in ['DA','DC','DG','DT']:
        return sequence[1]
    elif sequence in modified_nucleotides:
        return modified_nucleotides[sequence]["standard"]
    else:
        return None

def check_base_ribose_stack(nt1,nt2,pair_data):

    base_ribose_stack = ""

    # rotate nt2 to origin
    standard_nt2 = nt1.translate_rotate_component(nt2)

    ribose_atom = standard_nt2.centers["O5'"]

    if len(ribose_atom) == 3:

        x = ribose_atom[0]
        y = ribose_atom[1]
        z = ribose_atom[2]
        r = math.sqrt(x*x + y*y)
        q = math.sqrt(x*x + y*y + 2*(abs(z)-2.9)**2)           # rough measure of quality of location

        if r < 2.0 and abs(z) < 3.5:
            #print("%s with %s has O4' coordinates %0.4f,%0.4f,%0.4f" % (nt1.unit_id(),nt2.unit_id(),ribose_atom[0],ribose_atom[1],ribose_atom[2]))
            print('%s\t%s\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t=hyperlink("http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s")' % (nt1.unit_id(),nt2.unit_id(),x,y,z,r,q,nt1.unit_id(),nt2.unit_id()))

            #print("%s with %s has O4' coordinates " % (nt1.unit_id(),nt2.unit_id() ))
            #print(ribose_atom)

def check_coplanar(nt1,nt2,glycosidic_displacement):
    """ Given nt1 and nt2, check specific criteria to say that
    they are enough in the same plane to be called coplanar.
    If so, True, and return a number from 0 to 1 to measure the
    degree of coplanarity with 1 being best.
    Criteria for being coplanar or near coplanar:
      Pair.Gap must be < 97th percentile among basepairs (1.5179 Angstroms)
      min_distance must be < 97th percentile among basepairs (2.4589 A)
      Angle between center-center vector and normals must be > 70.2388 degrees
      Angle between normal vectors must be < 39.1315 degrees
    Return data that will be helpful for basepair classification.
    """

    pair_data = {}
    pair_data["coplanar"] = False
    pair_data["coplanar_value"] = -1         # 0 to 1 is coplanar, 1 is the best

    displ12 = np.matmul(glycosidic_displacement,nt1.rotation_matrix)  # vector from origin to nt2 when standardized
    pair_data["displ12"] = displ12

    gap12, points2 = calculate_basepair_gap(nt1,nt2)
    pair_data["gap12"] = gap12

    if gap12 >= 1.5179:
        return pair_data

    min_distance = 100
    points1 = []
    for atom in nt1.atoms():                  # nt1 atoms
        if not "'" in atom.name and not atom.name in ["P","OP1","OP2"]:  # exclude backbone atoms
            q = [atom.x, atom.y, atom.z]      # nt1 base atoms
            points1.append(q)                 # save for later gap21 calculation
            for p in points2:                 # nt2 base atoms
                d = np.linalg.norm(np.subtract(p,q))
                if d < min_distance:
                    min_distance = d

    pair_data["min_distance"] = min_distance

    # modified nucleotides don't have hydrogens, so be more flexible with them
    if min_distance >= 3.4589:
        return pair_data

    # if working with regular bases, insist on close contact
    if min_distance >= 2.4589 and nt1.sequence in ['A','C','G','U'] and nt2.sequence in ['A','C','G','U']:
        return pair_data

    center_displ = np.subtract(nt1.centers["base"],nt2.centers["base"])
    center_displ = center_displ / np.linalg.norm(center_displ) # normalize

    # calculate angle between center_displ and normal vectors to bases
    dot1 = abs(np.matmul(center_displ,nt1.rotation_matrix[:,2]))[0,0]
    if dot1 >= 0.3381:
        return pair_data

    dot2 = abs(np.matmul(center_displ,nt2.rotation_matrix[:,2]))[0,0]
    if dot2 >= 0.3381:
        return pair_data

    # calculate angle between normal vectors to the bases
    dot3 = abs(np.matmul(nt1.rotation_matrix[:,2].T,nt2.rotation_matrix[:,2]))
    if dot3 <= 0.7757:
        return pair_data

    gap21, points1 = calculate_basepair_gap(nt1,nt2,points1)

    if gap12 <  0.5062:             # 70th percentile
      Gap1Val = 1
    elif gap12 <  0.9775:           # 90th percentile
      Gap1Val = 1+(gap12- 0.5062)*(-1.0609)
    elif gap12 <  1.5179:           # 97th percentile
      Gap1Val = 0.5+(gap12- 0.9775)*(-0.9252)
    else:
      Gap1Val = 0

    if gap21 <  0.5062:             # 70th percentile
      Gap2Val = 1
    elif gap21 <  0.9775:           # 90th percentile
      Gap2Val = 1+(gap21- 0.5062)*(-1.0609)
    elif gap21 <  1.5179:           # 97th percentile
      Gap2Val = 0.5+(gap21- 0.9775)*(-0.9252)
    else:
      Gap2Val = 0

    if dot1 <  0.1139:              # 70th percentile
      dot1Val = 1
    elif dot1 <  0.2193:            # 90th percentile
      dot1Val = 1+(dot1- 0.1139)*(-4.7408)
    elif dot1 <  0.3381:            # 97th percentile
      dot1Val = 0.5+(dot1- 0.2193)*(-4.2103)
    else:
      dot1Val = 0

    if dot2 <  0.1139:              # 70th percentile
      dot2Val = 1
    elif dot2 <  0.2193:            # 90th percentile
      dot2Val = 1+(dot2- 0.1139)*(-4.7408)
    elif dot2 <  0.3381:            # 97th percentile
      dot2Val = 0.5+(dot2- 0.2193)*(-4.2103)
    else:
      dot2Val = 0

    if -dot3 < -0.9509:             # 70th percentile
      dot3Val = 1
    elif -dot3 < -0.8835:           # 90th percentile
      dot3Val = 1+(-dot3-(-0.9509))*(-7.4217)
    elif -dot3 < -0.7757:           # 97th percentile
      dot3Val = 0.5+(-dot3-(-0.8835))*(-4.6390)
    else:
      dot3Val = 0

    if min_distance <  1.8982:      # 70th percentile
      MinDistVal = 1
    elif min_distance <  2.1357:    # 90th percentile
      MinDistVal = 1+(min_distance- 1.8982)*(-2.1050)
    elif min_distance <  2.4859:    # 97th percentile
      MinDistVal = 0.5+(min_distance- 2.1357)*(-1.4280)
    else:
      MinDistVal = 0

    # Pair.Coplanar is 1 if all are within the 70th percentile
    # Pair.Coplanar is 0.5 if all are within the 90th percentile
    # Pair.Coplanar is > 0 if all are within the 97th percentile
    # Between these, it decreases linearly

    pair_data["coplanar"] = True
    pair_data["coplanar_value"] = min([Gap1Val, Gap2Val, dot1Val, dot2Val, dot3Val, MinDistVal])

    return pair_data

def calculate_basepair_gap(nt1,nt2,points2=None):

    displacements = []
    distances = []

    if points2:
        for p in points2:
            v = np.subtract(p,nt1.centers["base"])
            d = np.linalg.norm(v)
            displacements.append(v)
            distances.append(d)

    else:
        points2 = []
        for atom in nt2.atoms():
            if not "'" in atom.name and not atom.name in ["P","OP1","OP2"]:  # exclude backbone atoms
                p = [atom.x, atom.y, atom.z]
                v = np.subtract(p,nt1.centers["base"])
                d = np.linalg.norm(v)
                points2.append(p)
                displacements.append(v)
                distances.append(d)

    indices = np.argsort(distances)

    gap12 = 100
    for k in range(0,3):              # 3 nearest points
        p = displacements[indices[k]]
        z = abs(np.matmul(p,nt1.rotation_matrix[:,2])[0,0])  # distance out of plane of nt1
        if z < gap12:
            gap12 = z                 # gap is smallest z value

    return gap12, points2


def check_basepair_cutoffs(nt1,nt2,pair_data,cutoffs):
    """ Given nt1 and nt2 and the dictionary of cutoffs
    for that pair of nucleotides, check cutoffs for each
    basepair interaction type
    """

    displ = pair_data["displ12"]  # vector from origin to nt2 when standardized

    if abs(displ[0,2]) > 3.6:                            # too far out of plane
        return []

    ok_displacement_screen = []

    for interaction in cutoffs:                          # cWW, tWW, etc.
        for subcategory in cutoffs[interaction].keys():  # 0, 1, etc.
            cut = cutoffs[interaction][subcategory]
            if displ[0,2] < cut['zmin']:
                continue
            if displ[0,2] > cut['zmax']:
                continue
            if displ[0,0] < cut['xmin']:
                continue
            if displ[0,0] > cut['xmax']:
                continue
            if displ[0,1] < cut['ymin']:
                continue
            if displ[0,1] > cut['ymax']:
                continue
            ok_displacement_screen.append((interaction,subcategory)) # ("cWW",0), etc.

    if len(ok_displacement_screen) == 0:
        return []

#    return ok_displacement_screen

    rotation_1_to_2 = np.matmul(np.transpose(nt1.rotation_matrix), nt2.rotation_matrix)

    normal_Z = rotation_1_to_2[2,2]   # z component of normal vector to second base

#    print("%s with %s normal_Z is %0.4f" % (nt1.unit_id(),nt2.unit_id(),normal_Z))

    ok_normal = []

    for interaction,subcategory in ok_displacement_screen:
        cut = cutoffs[interaction][subcategory]
        if normal_Z < cut['normalmin']:
            continue
        if normal_Z > cut['normalmax']:
            continue
        ok_normal.append((interaction,subcategory))

    if len(ok_normal) == 0:
        return []

#    return ok_normal

    angle_in_plane = math.atan2(rotation_1_to_2[1,1],rotation_1_to_2[1,0])*57.29577951308232 - 90

    if angle_in_plane <= -90:
        angle_in_plane += 360

    ok_angle_in_plane = []

    for interaction,subcategory in ok_normal:
        cut = cutoffs[interaction][subcategory]
        if cut['anglemin'] < cut['anglemax']:     # for ranges in -90 to 270 like 50 to 120
            if angle_in_plane <= cut['anglemin']:
                continue
            if angle_in_plane >= cut['anglemax']:
                continue
        else:                                     # for ranges straddling 270 like 260 to -75
            if angle_in_plane >= cut['anglemax'] and angle_in_plane <= cut['anglemin']:
                continue

        ok_angle_in_plane.append((interaction,subcategory))

    if len(ok_angle_in_plane) == 0:
        return []

    ok_gap = []

    for interaction,subcategory in ok_angle_in_plane:
        cut = cutoffs[interaction][subcategory]
        if cut['gapmax'] > 0.1:
            if pair_data["gap12"] > cut['gapmax']:
                continue

        ok_gap.append((interaction,subcategory))

    return ok_gap





    interactions = {}
    return interactions

    """
      ro = N1.Rot'*N2.Rot;                       % rotation matrix from 1 to 2
  Pair.Normal = ro(:,3)';                    % normal to second plane

  if ro(3,3) > 0,                            % depending on orientation of 2,
    [ax,ang] = zAxisAngle(ro);               % rotation angle without a flip
  else
    [ax,ang] = zAxisAngle(ro*diag([-1 1 -1])); % flip base 2 first
  end

  Pair.Rot      = ro;
  Pair.RotAx    = ax';
  Pair.Ang      = ang;

  Pair.PlaneAng = acos(abs(N1.Rot(:,3)'*N2.Rot(:,3)))*57.29577951308232;
                                             % angle between planes

  % Calculate angle to rotate glycosidic bond of N2 counterclockwise in the plane to align with glycosidic of N1
  angle_in_plane = atan2(ro(2,2),ro(2,1))*57.29577951308232 - 90;

  if angle_in_plane <= -90,
    angle_in_plane = angle_in_plane + 360;
  end

%fprintf('zAnalyzepair:  Angles %8.4f  %8.4f\n', Pair.Ang, angle_in_plane);

  Pair.angle_in_plane = angle_in_plane;

  Lim(2,:) = [15 13 16 12];     % total number of atoms, including hydrogen

  d = zDistance(N2.Fit(1:Lim(2,Code2),:), N1.Center);
                                           % distances to base 1 center
  [sorteddistances,sortedindices] = sort(d);

  Pair.Gap = min(abs(N1.Rot(:,3)'*(N2.Fit(sortedindices(1:3),:)-ones(3,1)*N1.Center)'));

  min_distance = min(min(zDistance(N1.Fit,N2.Fit)));

  if isfield(N1,'Useangle_in_plane') && N1.Useangle_in_plane == 0,
    a = zCheckCutoffs(Pair.Displ,Pair.Normal,Pair.Ang,Pair.Gap,CL(:,:,Pair.Paircode));  % for comparison with historical
  else
    a = zCheckCutoffs(Pair.Displ,Pair.Normal,Pair.angle_in_plane,Pair.Gap,CL(:,:,Pair.Paircode));
  end
"""

def get_glycosidic_atom_coordinates(nt,parent):

    gly = None

    if nt.sequence in ['A','G','DA','DG']:
        gly = nt.centers["N9"]
    elif nt.sequence in ['C','U','DC','DT']:
        gly = nt.centers["N1"]
    elif nt.sequence in modified_nucleotides:
        if parent in ['A','G','DA','DG']:
            gly = nt.centers[modified_nucleotides[nt.sequence]["atoms"]["N9"]]
        elif parent in ['C','U','DC','DT']:
            gly = nt.centers[modified_nucleotides[nt.sequence]["atoms"]["N1"]]

    return gly

def get_axis_angle_from_rotation_matrix(rotation):
    """
    Turn a 3x3 rotation matrix into an axis of rotation and angle of rotation
    """

    values, vectors = np.linalg.eig(rotation) # get eigenvectors and eigenvalues of rotation

    imag0 = abs(values[0].imag)
    imag1 = abs(values[1].imag)
    imag2 = abs(values[2].imag)

    min_imag = np.argsort(np.absolute(np.imag(values)))[0]

    """
    if imag0 < imag1:
        if imag0 < imag2:
            min_imag = 0
        elif imag2 <= imag0:
            min_imag = 2
    else:
        if imag1 < imag2:
            min_imag = 1
        elif imag2 <= imag1:
            min_imag = 2
    """

    axis = np.real(np.array(vectors[:,min_imag]))  # column vector for axis

    angle = None
    i = np.argsort(np.absolute(axis),axis=0)      # find two largest entries of axis

    b = np.zeros((3,1))                # column vector of zeros
    b[i[1],0] = axis[i[2],0]
    b[i[2],0] = -axis[i[1],0]

    """
    print(values)
    print(vectors)
    print("    Eigenvalue with smallest imaginary part is %d" % min_imag)

    print("Axis:")
    print(axis)
    print(np.absolute(axis))
    print(i)
    print(b)
    print(b.T)
    print((b.T * rotation * b)[0,0])
    print(b.T.dot(b))

    print(np.matmul(np.matmul(b.T,rotation),b)[0,0])
    print(np.matmul(b.T,b)[0,0])
    """

    angle = math.acos(np.matmul(np.matmul(b.T,rotation),b)[0,0] / np.matmul(b.T,b)[0,0]).real
    angle = angle * np.sign(np.linalg.det(np.concatenate((b,np.matmul(rotation,b),axis),axis=1)))
    angle = angle * 57.29577951308232

    if angle <= -90:
        angle += 360

    return axis,angle


def type_of_interaction(base_residue, aa_residue, aa_coordinates, standard_aa_center, base_atoms):
    """ This function works with base and aa in standard position """
    squared_xy_dist_list = []

    """Defines different sets of amino acids"""
    planar_aa = set (["ARG", "ASN", "ASP", "GLU", "GLN", "HIS", "PHE", "TRP", "TYR"])
    stacked_aliphatic = set(["ALA", "CYS", "ILE", "LEU", "MET", "PRO", "SER", "THR", "VAL"])
    # Note:  LYS and GLY are not in the previous lists

    edge_to_edge_aa = set (["ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "LYS", "PHE", "SER", "THR", "TYR", "TRP"])
    shb_aa = set (["ARG", "ASN", "ASP", "GLU", "GLN", "HIS", "LYS", "SER", "THR", "TYR"])

    # calculate distances from aa atoms to base center
    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key  = aa_atom.name
        aa_x = aa_coordinates[key][0]
        aa_y = aa_coordinates[key][1]

        squared_xy_dist = (aa_x**2) + (aa_y**2)
        squared_xy_dist_list.append(squared_xy_dist)

    #print base_residue.unit_id(), aa_residue.unit_id(), min(squared_xy_dist_list), mean_z

    # for a stacking interaction, the x,y coordinate of at least one atom of the amino acid group needs to be
    # within sqrt(5) = 2.236 of the base center 0,0
    min_dist = np.sqrt(min(squared_xy_dist_list))

    if min_dist <= 2.236:
        #print base_residue.unit_id(), aa_residue.unit_id(), min(squared_xy_dist_list), mean_z
        if aa_residue.sequence in planar_aa:
            return stacking_planar_annotation(base_residue, aa_residue, min_dist)

        elif aa_residue.sequence in stacked_aliphatic:
            return stacking_non_planar_annotation(aa_residue, aa_coordinates, min_dist)

        else:
            return ("other-stack",{"dist-xy-from-center":min_dist})

    # check for interactions in the plane of the base
    mean_z = standard_aa_center[2]
    (num_hydrogen_bonds,hydrogen_bond_list) = count_hydrogen_bonds(base_residue, aa_residue, base_atoms)

    if -1.8 <= mean_z < 1.8:
        if aa_residue.sequence in edge_to_edge_aa:
            angle = calculate_angle_between_planar_residues(base_residue, aa_residue)
            if angle:
                if 0 <= angle <= 45 and num_hydrogen_bonds >= 2:
                    return ("pseudopair",{"hydrogen-bonds":hydrogen_bond_list,"angle-between-planes":angle})
                elif 45 <= angle:
                    return ("perpendicular-edge",{"hydrogen-bonds":hydrogen_bond_list,"angle-between-planes":angle})

        if aa_residue.sequence in shb_aa:
#            base_seq = base_residue.sequence
#            base_atoms = NAbaseheavyatoms[base_seq]
            if num_hydrogen_bonds >= 1:
                return ("SHB",{"hydrogen-bonds":hydrogen_bond_list})

        return ("other-edge",{"hydrogen-bonds":hydrogen_bond_list})

    return ("other",{"dist-xy-from-center":min_dist,"hydrogen-bonds":hydrogen_bond_list})

def count_hydrogen_bonds(base_residue, aa_residue, base_atoms):
    """Calculates number of Hydrogen bonds between amino acid part and base_part
    and returns the number and a list of pairs of atoms from (base,aa)
    """

    n = 0                            # number of independent hydrogen bonds being counted toward pseudopairs
    hydrogen_bond_list = []

    aa_key = aa_residue.sequence

    if not (aa_key in HB_donors or aa_key in HB_acceptors or aa_key in HB_weak_donors):
        return (n,hydrogen_bond_list)

    n0 = 0
    hb0 = []

    hb_angle_cutoff = 100
    screen_distance = 4.3                 # around 85% are higher than 4; 4 is a good first screen

    used_base_atoms = []
    used_aa_atoms = []
    carboxylate_donor_used = False   # track the two oxygens on ASP and GLU; only one can be a donor
    HIS_acceptor_used = False        # track the two nitrogens on HIS; only one can be an acceptor

    base_key = base_residue.sequence
    base_donors = HB_donor_hydrogens[base_key].keys()
    base_acceptors = HB_acceptors[base_key]
    base_HB_atoms = list(set(base_donors + base_acceptors))  # don't list atoms twice

    # these amino acids can be mis-modeled with the functional group flipped 180 degrees
    if aa_key in ["ASN","GLN"]:
        num_flip_states = 2
    else:
        num_flip_states = 1

    # check some amino acids in original and flipped orientations
    for flip in range(0,num_flip_states):
        n = 0
        hydrogen_bond_list = []
        if flip == 0:
            aa_donors = HB_donors[aa_key]
            if aa_key in HB_weak_donors:
                base_donors += HB_weak_donors[aa_key]
            aa_acceptors = HB_acceptors[aa_key]
            flip_name = ""
        else:
            # pretend that these three amino acids are flipped; they are often mis-modeled
            if aa_key == "ASN":
                aa_donors = ["OD1"]
                aa_acceptors = ["ND2"]
            elif aa_key == "GLN":
                aa_donors = ["OE1"]
                aa_acceptors = ["NE2"]
            # for now, don't try to recognize flipped HIS, just check for donors/acceptors
#                elif aa_key == "HIS":
#                    aa_donors = ['CD2', 'CE1']
#                    aa_acceptors = ['ND1', 'NE2']
            flip_name = "f"

        for base_atom in base_residue.atoms(name=base_HB_atoms):
            for aa_atom in aa_residue.atoms(name=aa_fg[aa_key]):
                difference = np.subtract(base_atom.coordinates(), aa_atom.coordinates())
                distance = np.linalg.norm(difference)

                if distance > screen_distance:
                    continue

                # check the distance between heavy atoms, depending on the two interacting atoms and the residues
                if ("O" in base_atom.name or "N" in base_atom.name) and ("O" in aa_atom.name or "N" in aa_atom.name):
                    h_bond_ideal_distance  = distance_limit_coulombic[base_residue.sequence][base_atom.name]
                    h_bond_ideal_distance += distance_limit_coulombic[aa_residue.sequence][aa_atom.name]
                    h_bond_over_distance = distance - h_bond_ideal_distance

                    hbondtext = base_residue.sequence +"\t"+ base_atom.name +"\t"+ str(distance_limit_coulombic[base_residue.sequence][base_atom.name]) +"\t"
                    hbondtext += aa_residue.sequence +"\t"+ aa_atom.name +"\t"+ str(distance_limit_coulombic[aa_residue.sequence][aa_atom.name])
                    hbondtext += "\t"+ str(h_bond_ideal_distance) +"\t"+ str(distance) +"\t"+ str(h_bond_over_distance)

                else:
                    h_bond_ideal_distance  = distance_limit_vdw[base_residue.sequence][base_atom.name]
                    h_bond_ideal_distance += distance_limit_vdw[aa_residue.sequence][aa_atom.name]
                    h_bond_over_distance = distance - h_bond_ideal_distance

                    hbondtext = base_residue.sequence +"\t"+ base_atom.name +"\t"+ str(distance_limit_vdw[base_residue.sequence][base_atom.name]) +"\t"
                    hbondtext += aa_residue.sequence +"\t"+ aa_atom.name +"\t"+ str(distance_limit_vdw[aa_residue.sequence][aa_atom.name])
                    hbondtext += "\t"+ str(h_bond_ideal_distance) +"\t"+ str(distance) +"\t"+ str(h_bond_over_distance)

#                print("hbond\t"+hbondtext)

                if distance > h_bond_ideal_distance + 0.4:
                    continue
                    h_bond_over_distance = distance - h_bond_ideal_distance

                if base_atom.name in base_donors and aa_atom.name in aa_acceptors:
                    # loop through hydrogens whose locations are known
                    for hydrogen_atom in base_residue.atoms(name=HB_donor_hydrogens[base_key][base_atom.name]):
                        hb_angle = calculate_hb_angle(base_atom.coordinates(),hydrogen_atom.coordinates(),aa_atom.coordinates())
#                            print("hb_angle %10.8f" % hb_angle)
                        if hb_angle > hb_angle_cutoff:
                            if not base_atom.name in used_base_atoms and not aa_atom.name in used_aa_atoms:
                                if aa_key == "HIS":
                                    if not HIS_acceptor_used:
                                        n = n + 1
                                    HIS_acceptor_used = True
                                else:
                                    n = n + 1
                            hydrogen_bond_list.append((base_atom.name,hydrogen_atom.name,aa_atom.name+flip_name,h_bond_ideal_distance,distance,hb_angle))
                            used_base_atoms.append(base_atom.name)
                            used_aa_atoms.append(aa_atom.name)

                    # for O2', coordinates of H are not known, so check angles separately
                    if base_atom.name == "O2'":
                        hb_angle = calculate_hb_angle(base_residue.centers["C2'"],base_atom.coordinates(),aa_atom.coordinates())
#                        print("C2'-O2'-A angle %10.8f" % hb_angle)
                        if hb_angle > 80:
                            if not base_atom.name in used_base_atoms and not aa_atom.name in used_aa_atoms:
                                if aa_key == "HIS":
                                    if not HIS_acceptor_used:
                                        n = n + 1
                                    HIS_acceptor_used = True
                                else:
                                    n = n + 1
                            hydrogen_bond_list.append((base_atom.name,"O2'H",aa_atom.name+flip_name,h_bond_ideal_distance,distance,hb_angle))
                            used_base_atoms.append(base_atom.name)
                            used_aa_atoms.append(aa_atom.name)


                elif base_atom.name in base_acceptors and aa_atom.name in aa_donors:
                    # check OH group on certain amino acids; coordinates of H are not known
                    if (aa_key == "SER" and aa_atom.name == "OG") or \
                       (aa_key == "THR" and aa_atom.name == "OG1") or \
                       (aa_key == "TYR" and aa_atom.name == "OH"):

                        if aa_key == "TYR":
                            carbon = "CZ"
                        else:
                            carbon = "CB"

                        hb_angle = calculate_hb_angle(aa_residue.centers[carbon],aa_atom.coordinates(),base_atom.coordinates())
#                        if hb_angle:
#                            print("%s-OH-%s angle %10.8f ------------ http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (carbon,base_atom.name,hb_angle,base_residue.unit_id(),aa_residue.unit_id()))

                        if hb_angle > 80:
                            if not base_atom.name in used_base_atoms and not aa_atom.name in used_aa_atoms:
                                n = n + 1
                            hydrogen_bond_list.append((base_atom.name,"OHH",aa_atom.name+flip_name,h_bond_ideal_distance,distance,hb_angle))
                            used_base_atoms.append(base_atom.name)
                            used_aa_atoms.append(aa_atom.name)

                    if not base_atom.name in used_base_atoms and not aa_atom.name in used_aa_atoms:
                        # aa with carboxylate, only count one oxygen as a donor
                        if aa_key == "ASP" or aa_key == "GLU":
                            if not carboxylate_donor_used:
                                n = n + 1
                            carboxylate_donor_used = True
                        else:
                            n = n+1
                    hydrogen_bond_list.append((base_atom.name,"H?",aa_atom.name+flip_name,h_bond_ideal_distance,distance,0))
                    used_base_atoms.append(base_atom.name)
                    used_aa_atoms.append(aa_atom.name)
                if flip == 0 and num_flip_states == 2:
                    n0 = n
                    hb0 = hydrogen_bond_list

    if num_flip_states == 2:
        if n0 >= n:       # second flip is not strictly better
            n = n0        # use the first set of hydrogen bonds
            hydrogen_bond_list = hb0
        else:
            print("  Found a flipped amino acid "+aa_residue.unit_id()+" "+base_key+" "+aa_key+" "+str(n)+" $$$$$$$$$$$$$$$$$")
#            print(hydrogen_bond_list)

    return (n,hydrogen_bond_list)

def stacking_planar_annotation (base_residue, aa_residue, min_dist):
    """ For planar amino acids, determine the stacking classification
    according to the angle between the plane of the amino acid and
    the plane of the base """

    angle = calculate_angle_between_planar_residues(base_residue, aa_residue)

    # cation is about the type of amino acid.  List them ... HIS is positive sometimes.
    #
    perpendicular_stack_aa = set(["HIS", "PHE", "TRP", "TYR"])
    perpendicular_aa = set (["HIS", "ARG", "LYS", "ASN", "GLN"])

    if angle:
        if angle <= 45:
            return ("pi-pi-stacking",{"angle-between-planes":angle})
        elif angle > 45:
            if aa_residue.sequence in perpendicular_stack_aa:
                return ("perpendicular-stacking",{"angle-between-planes":angle,"dist-xy-from-center":min_dist})
            elif aa_residue.sequence in perpendicular_aa:
                return ("cation-pi",{"angle-between-planes":angle,"dist-xy-from-center":min_dist})

    return ("other-stack",{"angle-between-planes":None,"dist-xy-from-center":min_dist})

def stacking_non_planar_annotation(aa_residue, aa_coordinates, min_dist):
    """ For non-planar amino acids, determine the stacking type
    by looking at how spread out the atoms are above/below the
    plane of the base """
    baa_dist_list = []

    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key = aa_atom.name
        aa_z = aa_coordinates[key][2]
        baa_dist_list.append(aa_z)
    max_baa = max(baa_dist_list)
    min_baa = min(baa_dist_list)
    diff = max_baa - min_baa
    #print aa_residue.unit_id(), diff
    if diff <= tilt_cutoff[aa_residue.sequence]:
        return ("stacked",{"stacking-diff":diff,"dist-xy-from-center":min_dist})

    return ("other-stack",{"stacking-diff":diff,"dist-xy-from-center":min_dist})

def calculate_angle_between_planar_residues (base_residue, aa_residue):
    vec1 = normal_vector_calculation(base_residue)
    vec2 = normal_vector_calculation(aa_residue)

    if len(vec1) == 3 and len(vec2) == 3:
        angle = smaller_angle_between_vectors(vec1, vec2)
        return angle
    else:
        print("Missing a normal vector %s %s" % (base_residue.unit_id(),aa_residue.unit_id()))
        return None

def normal_vector_calculation(residue):
    key = residue.sequence
    P1 = residue.centers[planar_atoms[key][0]]
    P2 = residue.centers[planar_atoms[key][1]]
    P3 = residue.centers[planar_atoms[key][2]]
#    print key, residue.unit_id(), P1, P2, P3

    if len(P1) == 3 and len(P2) == 3 and len(P3) == 3:
        normal_vector = np.cross((P2 - P1),(P3 - P1))
        return normal_vector
    else:
        return []

# this function calculates the angle made from A to B to C from 0 to 180 degrees
def calculate_hb_angle(A,B,C):
    if len(A) == 3 and len(B) == 3 and len(C) == 3:
        return angle_between_vectors(np.subtract(A,B),np.subtract(C,B))

# This function calculates an angle from 0 to 90 degrees between two vectors
def smaller_angle_between_vectors(vec1, vec2):
    if len(vec1) == 3 and len(vec2) == 3:
        # the following line sometimes causes "RuntimeWarning: invalid value encountered in double_scalars" on 5JTE
        cosang = abs(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)))
        angle = np.arccos(cosang)
        return 180*abs(angle)/np.pi
    else:
        return None

# This function calculates an angle from 0 to 180 degrees between two vectors
def angle_between_vectors(vec1, vec2):
    if len(vec1) == 3 and len(vec2) == 3:
        cosang = np.dot(vec1, vec2)
        sinang = np.linalg.norm(np.cross(vec1, vec2))
        angle = np.arctan2(sinang, cosang)
        return 180*angle/np.pi
    else:
        return None

# This function calculates an angle from 0 to 180 degrees between two vectors
def angle_between_three_points(P1,P2,P3):
    if len(P1) == 3 and len(P2) == 3 and len(P3) == 3:
        return angle_between_vectors(P1-P2,P3-P2)
    else:
        return None

# This function calculates an angle from 0 to 180 degrees between two vectors
def distance_between_vectors(vec1, vec2):
    if len(vec1) == 3 and len(vec2) == 3:
        return np.linalg.norm(np.subtract(vec1,vec2))
    else:
        return None

def detect_base_edge(base_residue, base_coordinates, aa_residue, aa_coordinates):
    aa_x = []
    aa_y = []
    base_x = []
    base_y = []
    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key = aa_atom.name
        aa_x.append(aa_coordinates[key][0])
        aa_y.append(aa_coordinates[key][1])

    aa_center_x = np.mean(aa_x)
    aa_center_y = np.mean(aa_y)

    NAsixsidedringatoms = ['N1','C2','N3','C4','C5','C6']

#    for base_atom in base_residue.atoms(name=NAbaseheavyatoms[base_residue.sequence]):
    for base_atom in base_residue.atoms(name=NAsixsidedringatoms):
        key = base_atom.name
        base_x.append(base_coordinates[key][0])
        base_y.append(base_coordinates[key][1])

    base_center_x = np.mean(base_x)
    base_center_y = np.mean(base_y)

    """
    reference_atom = {'A':'C2', 'C':'O2', 'G':'N2', 'U':'O2'}    # for WC versus Sugar edge
    reference_atom = {'A':'N6', 'C':'N4', 'G':'O6', 'U':'O4'}    # for WC versus Hoogsteen edge
    x = base_coordinates[reference_atom[base_residue.sequence]][0] - base_center_x
    y = base_coordinates[reference_atom[base_residue.sequence]][1] - base_center_y
    angle_aa = np.arctan2(y,x)         # values -pi to pi
    angle_deg = (180*angle_aa)/np.pi # values -180 to 180
    print("Base %s has S/WC reference angle %10.8f" % (base_residue.sequence,angle_deg))

    reference_atom = {'A':'N9', 'C':'N1', 'G':'N9', 'U':'N1'}    # for WC versus Hoogsteen edge
    print("Base %s has N1/N9 x value %10.8f" % (base_residue.sequence,base_coordinates[reference_atom[base_residue.sequence]][0]))
    """

    x = aa_center_x - base_center_x
    y = aa_center_y - base_center_y
    angle_aa = np.arctan2(y,x)         # values -pi to pi
    angle_deg = (180*angle_aa)/np.pi # values -180 to 180

    purine = set(["A", "G", "DA", "DG"])
    pyrimidine = set(["C", "U", "DC", "DT"])

    if base_residue.sequence in purine:
        if -14 <= angle_deg <= 104:
            return ("fgWC",angle_deg)
        elif 104 < angle_deg or aa_center_x < -1.3:
            return ("fgH",angle_deg)
        else:
            return ("fgS",angle_deg)

    elif base_residue.sequence in pyrimidine:
        if -32 <= angle_deg <= 86:
            return ("fgWC",angle_deg)
        elif 86 < angle_deg or aa_center_x < -0.37:
            return ("fgH",angle_deg)
        else:
            return ("fgS",angle_deg)

def detect_face(aa_residue, aa_coordinates):
    aa_z =[]

    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key = aa_atom.name
        aa_z.append(aa_coordinates[key][2])

    mean_z = np.mean(aa_z)
    if mean_z <= 0:
        return ("fgs5",mean_z)
    else:
        return ("fgs3",mean_z)

def unit_vector(v):
    return v / np.linalg.norm(v)


def load_basepair_annotations(filename,all_pair_types):

    with open(filename,'rb') as opener:
        basepairs = pickle.load(opener)

    pair_to_bp_type={}

    for bp_type in all_pair_types:
      for (u1,u2,c) in basepairs[bp_type]:
         pair_to_bp_type[(u1,u2)]=bp_type

      near_bp_type="n"+bp_type
#      for (u1,u2,c) in basepairs[near_bp_type]:
#         pair_to_bp_type[(u1,u2)]=near_bp_type

    return pair_to_bp_type


def compare_annotations(pairs_1,pairs_2,all_pair_types):

    allkeys = set(pairs_1.keys()) | set(pairs_1.keys())

    with open("comparison.txt", mode='w') as file:

        file.write("Unit id 1"+"\t"+"Unit id 2"+"\t"+"Base Combination"+"\t"+"Annotation 1"+"\t"+"Annotation 2"+"\t"+"Match Count"+"\t"+"url"+"\n")

        total_annotations  = 0
        count_agreements = 0
        for key in allkeys:
            #print(key)
            x = key[0].split("|")
            y = key[1].split("|")
            base_combination = x[3]+y[3]

            pair_type=""
            if key in pairs_1:
                pair_type = pairs_1[key].replace("n","")

            elif key in pairs_2:
                pair_type = pairs_2[key]

            #print([pair_type])
            if pair_type in ['cHW', 'tHW', 'cSH', 'tSH', 'cSW', 'tSW']:
                continue

            if pair_type in ['cWW', 'tWW','cSS', 'tSS', 'cHH', 'tHH', '0  ']:
                if base_combination in ['UA', 'CG', 'UG', 'GA', 'CA', 'UC']:
                    continue
                if base_combination in ['AA', 'CC', 'GG', 'UU']:
                    if int(x[4])>=int(y[4]):
                        continue

            if key in pairs_1:
                ann1 = pairs_1[key]
            else:
                ann1 = ""
            if key in pairs_2:
                ann2 = pairs_2[key]
            else:
                ann2 = ""

            ## Keep track of how often ann1 and ann2 agree over all pairs
            ## New Python annotations use upper and lowercase to indicate
            ## details of which part of the edge is used with cWw and tHh and such,
            ## so change to lowercase to compare

            if len(ann1) > 0 or len(ann2) > 0:
                total_annotations += 1

            if ann1.lower() == ann2.lower() and len(ann2) > 0:
                count_agreements +=  1

            ## Counting the number of times the annotations match with each other
            ## for this particular pair of neucleotides
            count = 0
            if ann1.lower() == ann2.lower()  and len(ann2) > 0:
                count = 1
            elif ann2.lower() in ann1.lower() and len(ann2) > 0:
                count += 0.6

            url= "http://rna.bgsu.edu/rna3dhub/display3D/unitid/" + key[0] + "," +  key[1]

            if count < 1:
                file.write(key[0]+"\t"+key[1]+"\t"+x[3]+y[3]+"\t"+ann1+"\t"+ann2+"\t"+str(count)+"\t"+url+"\n")

    print("  Number of times ann1 and ann2 agree %d" % count_agreements)
    print("  Total number of ann1 or ann2 annotations %d" % total_annotations)
    print("  Percentage of times ann1 and ann2 agree %0.4f" % (count_agreements*100/total_annotations))


def draw_base(base_seq, ax):
    """Connects atoms to draw neighboring bases and amino acids for 3D plots"""
     #creates lists of rotated base coordinates
    for basecoord_list in list_base_coord:
        new_base_x = []
        new_base_y = []
        new_base_z = []

        back_base_x = []
        back_base_y = []
        back_base_z = []


        try:
            for atomname in RNAconnections[base_seq]:
                coord_base = []
                coord_base= basecoord_list[atomname]
                new_base_x.append(coord_base[0])
                new_base_y.append(coord_base[1])
                new_base_z.append(coord_base[2])
            base_lines= ax.plot(new_base_x, new_base_y, new_base_z, label= 'Base')
            #ax.scatter(basecenter[0], basecenter[1], basecenter[2], zdir='y', color='b', marker='o')
            #ax.scatter(x = 0, y= 0, z= 0, color='b', marker='o')
            plt.setp(base_lines, 'color', 'b', 'linewidth', 1.0)

            for atomname in Ribophos_connect[base_seq]:
                back_base=[]
                back_base= basecoord_list[atomname]
                back_base_x.append(back_base[0])
                back_base_y.append(back_base[1])
                back_base_z.append(back_base[2])
            base_lines= ax.plot(back_base_x, back_base_y, back_base_z, label= 'Base')
            plt.setp(base_lines, 'color', 'g', 'linewidth', 1.0)
            #ax.text(9, 1, 1, base_residue)
        except:
            print "Missing residues"
            continue


def text_output(result_list):
    with open(outputText % PDB, 'wb') as target:
        for result in result_list:
            target.write(str(result))
            target.write("\r\n")
            target.close

def csv_output(result_list):
    with open(outputNAPairwiseInteractions % PDB, 'wb') as csvfile:
        fieldnames = ['RNA ID', 'AA ID', 'RNA Chain ID', 'RNA residue','RNA residue number','Protein Chain ID', 'AA residue','AA residue number', 'Interaction', 'Edge', 'Param']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for base_residue, aa_residue, interaction, edge, standard_aa_center, param in result_list:
            base = base_residue.unit_id()
            aa = aa_residue.unit_id()
            #print base, aa, interaction
            base_component = str(base).split("|")
            aa_component = str(aa).split("|")
            writer.writerow({'RNA ID': base, 'AA ID': aa, 'RNA Chain ID': base_component[2], \
                'RNA residue':base_component[3],'RNA residue number': base_component[4],\
                'Protein Chain ID':aa_component[2],'AA residue': aa_component[3],\
                'AA residue number': aa_component[4], 'Interaction': interaction, 'Edge': edge, 'Param': param})

        """for base_residue, aa_residue,interaction in result_list:
                    base_component = str(base_residue).split("|")
                    aa_component = str(aa_residue).split("|")
                    writer.writerow({'RNA Chain ID': base_component[2], 'RNA residue':base_component[3],\
                    'RNA residue number': base_component[4],'Protein Chain ID':ChainNames[PDB][aa_component[2]],\
                    'AA residue': aa_component[3],'AA residue number': aa_component[4], 'Interaction': interaction})"""

def WriteProteinUnits(PDB,amino_acids):
    all_aa_list = []
    for aa in amino_acids:
        fields = aa.unit_id().split("|")
        my_tuple = (aa.unit_id(),fields[2],aa.index,aa.centers['aa_fg'])
        all_aa_list.append(my_tuple)

    new_all_aa_list = sorted(all_aa_list, key=lambda x: (x[1],x[2]))

    ids = []
    chainPositions = []
    centers = []
    rotations = []
    for my_tuple in new_all_aa_list:
        if len(my_tuple[3]) == 3:
#            print(my_tuple)
            ids.append(my_tuple[0])
            chainPositions.append(my_tuple[2])
            centers.append(my_tuple[3])
            rotations.append(np.zeros((0,3,3)))
        else:
            print("missing center",my_tuple)

    dataPathUnits = "C:/Users/zirbel/Dropbox/2018 FR3D intersecting pairs/data/units/"

    with open(dataPathUnits + PDB + '_protein.pickle', 'wb') as handle:
        pickle.dump((ids,chainPositions,centers,rotations), handle, protocol = 2)  # protocol 2 is safe for Python 2.7


def writeInteractionsHTML(allInteractionDictionary,outputHTML,version):

    SERVER = True
    SERVER = False
    if not SERVER:
        JS1 = '<script src="./js/JSmol.min.nojq.js"></script>'
        JS2 = '<script src="./js/jquery.jmolToolsRNAProtein.js"></script>'
        JS3 = '<script src="./js/imagehandlinglocal.js"></script>'
        JS4 = '<script src="./js/jmolplugin.js" type="text/javascript"></script>'
        JS5 = '<script type="text/javascript" src="./js/heatmap.js"></script>'
        JS6 = '<script src="./js/scroll_table_addon.js"></script>'   # scrollable table
        JS7 = '<script src="./js/table_js.js"></script>'   # scrollable table
        CSS1 = '<link rel="stylesheet" type="text/css" href="./css/table_style.css">'
        CSS2 = '<link rel="stylesheet" type="text/css" href="./css/scroll_table_addon.css"/>'

    else:
        JS1 = '<script src="http://rna.bgsu.edu/rna3dhub/js/jsmol/JSmol.min.nojq.js"></script>'
        JS2 = '<script src="http://rna.bgsu.edu/rna3dhub/js/jquery.jmolTools.js"></script>'
        JS2 = '<script src="./js/jquery.jmolToolsRNAProtein.js"></script>'   # special code to superimpose bases
        JS3 = '<script src="http://rna.bgsu.edu/webfr3d/js/imagehandling.js"></script>'
        JS4 = '<script src="http://rna.bgsu.edu/webfr3d/js/jmolplugin.js" type="text/javascript"></script>'
        JS5 = '<script type="text/javascript" src="http://rna.bgsu.edu/webfr3d/js/heatmap.js"></script>'
        JS6 = '<script src="./js/scroll_table_addon.js"></script>'   # scrollable table
        JS7 = '<script src="./js/table_js.js"></script>'   # scrollable table
        CSS1 = '<link rel="stylesheet" type="text/css" href="./css/table_style.css">'
        CSS2 = '<link rel="stylesheet" type="text/css" href="./css/scroll_table_addon.css"/>'

    # not working yet, so omit:
    JS6 = ""
    JS7 = ""

    count_pair = 0

    for key in allInteractionDictionary:
        pagetitle = key.replace(" ","-")
        htmlfilename = key.replace(" ","-") + version

#        print("Writing HTML file for "+key+" in "+htmlfilename+".html, found "+ str(len(allInteractionDictionary[key])) + " instances")

        fields = key.split("_")
        print(fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t"+str(len(allInteractionDictionary[key])))
        count_pair += len(allInteractionDictionary[key])

        # limit the number of instances shown, to be able to compute and display discrepancy
        numForDiscrepancy = min(300,len(allInteractionDictionary[key]))

        # calculate discrepancies between all instances, up to 300
        discrepancy = np.zeros((numForDiscrepancy,numForDiscrepancy))
        for i in range(0,numForDiscrepancy):
            instance_1 = allInteractionDictionary[key][i]
            aa_1 = instance_1[4]
            for j in range(i+1,numForDiscrepancy):
                instance_2 = allInteractionDictionary[key][j]
                aa_2 = instance_2[4]
                s = 0
                for atom_name in aa_fg[aa_1.sequence]:
                    d = distance_between_vectors(aa_1.centers[atom_name],aa_2.centers[atom_name])
                    if d:
                        s += d**2

                discrepancy[i][j] = np.sqrt(s)/len(aa_fg[aa_1.sequence])  # average distance between corresponding atoms
                discrepancy[j][i] = discrepancy[i][j]
                # discrepancy[j][i] = np.linalg.norm(standard_aa_center_1 - standard_aa_center_2)


        # base_aa = (base_residue, aa_residue, interaction, edge, standard_aa_center, param)

        # use greedy insertion 100 times to find a decent ordering of the instances
        newOrder, bestPathLength, distances = orderWithPathLengthFromDistanceMatrix(discrepancy,10)

        # rewrite the list of instances, limiting it to numForDiscrepancy
        newList = []
        for i in range(0,len(newOrder)):
            newList.append(allInteractionDictionary[key][newOrder[i]])
        allInteractionDictionary[key] = newList

        # write out text for radio boxes to display each individual interaction
        i = 1
        queryNote = "<h2>"+pagetitle+"</h2>\n"

        candidatelist = '<table style="white-space:nowrap;" id="table">\n'
        candidatelist += '<thead><tr><th><span onclick="sortTable(1)">Number</span></th>'
        candidatelist += '<th>View</th>'
        candidatelist += '<th><span onclick="sortTable(3)">Nucleotide</span></th>'
        candidatelist += '<th><span onclick="sortTable(4)">Amino acid</span></th>'
        candidatelist += '<th><span onclick="sortTable(5)">Interaction</span></th>'
        candidatelist += '<th><span onclick="sortTable(6)">Edge</span></th>'
        candidatelist += '<th><span onclick="sortTable(7)">a.a. x</span></th>'
        candidatelist += '<th><span onclick="sortTable(8)">a.a. y</span></th>'
        candidatelist += '<th><span onclick="sortTable(9)">a.a. z</span></th>'
        param = allInteractionDictionary[key][0][5]
        param_list = sorted(list(param.keys()))
        col = 9
        for header in param_list:
            col += 1
            candidatelist += '<th><span onclick="sortTable(%d)">%s</span></th>' % (col,header)
        candidatelist += "</tr></thead>\n"
        candidatelist += '<tbody id="table_rows">\n'
        for base_id, aa_id, interaction, edge, standard_aa, param in allInteractionDictionary[key]:
            candidatelist += '<tr><td>'+str(i)+'.</td><td><label><input type="checkbox" id="'+str(i-1)+'" class="jmolInline" data-coord="'
            candidatelist += base_id +","+ aa_id
            candidatelist += '">&nbsp</td>'
            candidatelist += '<td>%s</td>' % base_id
            candidatelist += '<td>%s</td>' % aa_id
            candidatelist += '<td>%s</td>' % interaction
            candidatelist += '<td>%s</td>' % edge
            for header in param_list:
                if isinstance(param[header],float):
                    candidatelist += '<td>%0.4f</td>' % param[header]
                elif isinstance(param[header],list):
                    for hbond in param[header]:
                        candidatelist += '<td>%s-%s-%s,%0.2fA,%0.2fA,%0.2fd</td>' % hbond
                else:
                    candidatelist += '<td>Error</td>'

            candidatelist += '</tr>\n'
            i += 1
        candidatelist += '</tbody></table>\n'
        candidatelist = candidatelist

        # write out text to tell what values to put in the heat map
        discrepancyText = ''
        for c in range(0,numForDiscrepancy):
            instance1 = allInteractionDictionary[key][c][0]  # id of base
            for d in range(0,numForDiscrepancy):
                instance2 = allInteractionDictionary[key][d][0]  # id of base

                discrepancyText += '{"discrepancy": ' + str(discrepancy[newOrder[c]][newOrder[d]])
                discrepancyText += ', "ife1": "' + str(c+1) + "-" + instance1 + '", "ife1_index": ' + str(c)
                discrepancyText += ', "ife2": "' + str(d+1) + "-" + instance2 + '", "ife2_index": ' + str(d) + '}'
                if c < numForDiscrepancy-1 or d < numForDiscrepancy-1:
                    discrepancyText += ',\n'

        # read template.html into one string
#        with open(outputHTML+'/localtemplate.html', 'r') as myfile:
        with open('template.html', 'r') as myfile:
            template = myfile.read()

        # replace ###PAGETITLE### with pagetitle
        template = template.replace("###PAGETITLE###",pagetitle)

        # replace ###CANDIDATELIST### with candidatelist
        template = template.replace("###CANDIDATELIST###",candidatelist)

        # replace ###DISCREPANCYDATA### with discrepancyText
        discrepancyText = "var data =  [\n" + discrepancyText + "]"
        discrepancyText = '<script type="text/javascript">\n' + discrepancyText + '\n</script>'
        template = template.replace("###DISCREPANCYDATA###",discrepancyText)

        template = template.replace("###QUERYNAME###",queryNote)
        template = template.replace("###JS1###",JS1)
        template = template.replace("###JS2###",JS2)
        template = template.replace("###JS3###",JS3)
        template = template.replace("###JS4###",JS4)
        template = template.replace("###REFRESH###","")
        template = template.replace("###JS5###",JS5)    # include heatmap.js code
        template = template.replace("###JS6###",JS6)
        template = template.replace("###JS7###",JS7)
        template = template.replace("###CSS1###",CSS1)
        template = template.replace("###CSS2###",CSS2)

        # write htmlfilename
        with open(outputHTML+'/'+htmlfilename+'.html', 'w') as myfile:
            myfile.write(template)

    print("Wrote out %d pairwise interactions" % count_pair)

        # upload the files to /var/www/html/RNAprotein

def write_unit_data_file(PDB,unit_data_path,structure):
    """
    Write out data file(s) of nucleotide centers and rotation matrices,
    primarily for use by the FR3D motif search tool.
    If unit_data_path is empty, no files are written.
    One file for each chain.
    """

    print("Writing data file")
    print(structure.pdb)
    print(structure.model)
    print(structure.symmetry)
    print(structure.chain)

    if len(unit_data_path) > 0:
        # get a list of RNA and DNA chains in the structure
        # loop over the chains
        # get the nucleotides in the chain in sequence order
        # write out data for each nucleotide including its sequence position

        filename = unit_data_path + "/" + PDB + "_NA_base_rotation.pickle"
        if not os.path.exists(filename):

            units = []
            order = []
            cntrs = []
            rttns = []

            for nt in nucleotides:
                units.append(nt.unit_id())
                order.append(nt.index)
                cntrs.append(nt.centers["base"])
                rttns.append(nt.rotation_matrix)

            rsset = [units, order, cntrs, rttns]

            with open(filename, 'wb') as fh:
                # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
                pickle.dump(rsset, fh, 2)

#=======================================================================

PDB_List = ['4TNA']
PDB_List = ['4ARC']
PDB_List = ['4ARC','4V9F']

ReadPickleFile = True                  # when true, just read the .pickle file from a previous run
ReadPickleFile = False                 # when true, just read the .pickle file from a previous run

base_seq_list = ['DA','DT','DC','DG']  # for DNA
base_seq_list = ['A','U','C','G']      # for RNA
base_seq_list = ['A','U','C','G']      # for RNA
base_seq_list = []                     # for all nucleic acids, modified or not

nt_reference_point = "base"
nt_nt_screen_distance = 10.5    # minimum distance between base centers to screen for interaction
atom_atom_min_distance = 5    # minimum distance between atoms in nts to consider them interacting

# plot one instance of each of the pairwise interactions
PlotPair = True
PlotPair = False
AlreadyPlotted = {}

ShowStructureReadingErrors = True

unit_data_path = "C:/Users/zirbel/Documents/GitHub/fr3d-python/data/units"

if __name__=="__main__":

    timerData = myTimer("start")

    allInteractionDictionary = defaultdict(list)

    result_nt_nt = []               # for accumulating a complete list over all PDB files

    timerData = myTimer("Making PDB list",timerData)

    allOutputDataFile = outputNAPairwiseInteractions + "AllInteractions.pickle"

    if ReadPickleFile:
        print("Reading " + outputDataFile)
        timerData = myTimer("Reading pickle file",timerData)
        allInteractionDictionary,allAATwoBaseDictionary,PDB_List = pickle.load(open(allOutputDataFile,'rb'))
        writeInteractionsHTML(allInteractionDictionary,outputHTML,version)
    else:
        PDB_IFE_Dict = defaultdict(str)      # accumulate PDB-IFE pairs
        for PDB in PDB_List:
            if "nrlist" in PDB and "NR_" in PDB:
                                          # referring to an equivalence class online
                                          # download the entire representative set,
                                          # then find the right line for the equivalence class
                                          # then extract the list
                f = urllib.urlopen(PDB)
                myfile = f.read()
                alltbody = myfile.split("tbody")
                alllines = alltbody[1].split("<a class='pdb'>")
                del alllines[0]
                for line in alllines:
                    fields = line.split("</a>")
                    if len(fields[0]) > 1:
                        newIFE = fields[0].replace(" ","")   # remove spaces
                        newPDB = newIFE[0:4]
                        PDB_IFE_Dict[newPDB] += "+" + newIFE

            elif "nrlist" in PDB:           # referring to a representative set online
                f = urllib.urlopen(PDB)
                myfile = f.read()
                alllines = myfile.split("\n")
                for line in alllines:
                    fields = line.split(",")

                    if len(fields) > 1 and len(fields[1]) > 4:
                        newPDB = fields[1][1:5]       # use only PDB identifier, ignore IFE for now
                        PDB_IFE_Dict[newPDB] += "+" + fields[1].replace('"','')

            elif "+" in PDB:                      # in case multiple chains in an IFE
                newPDB = PDB.split("|")[0]        # in case model, chain is indicated
                PDB_IFE_Dict[newPDB] = PDB

            elif "|" in PDB:                      # in case model, chain is indicated
                newPDB = PDB.split("|")[0]
                PDB_IFE_Dict[newPDB] = PDB

            else:
                PDB_IFE_Dict[PDB] = ""            # indicates to process the whole PDB file

        print("PDB_IFE_Dict is %s" % PDB_IFE_Dict)

        counter = 0
        count_pair = 0

        # loop through 3D structures and annotate interactions
        PDBs = PDB_IFE_Dict.keys()
        #PDBs = PDBs[::-1]  # reverse the order of the list, for debugging

        for PDB in PDBs:
            counter += 1

            outputDataFileCSV =    outputNAPairwiseInteractions + PDB + ".csv"
            outputDataFilePickle = outputNAPairwiseInteractions + PDB + ".pickle"

            print("Reading file " + PDB + ", which is number "+str(counter)+" out of "+str(len(PDB_IFE_Dict)))
            timerData = myTimer("Reading CIF files",timerData)

            if ShowStructureReadingErrors:
                # do this to make sure to see any error messages
                structure = get_structure(inputPath % PDB,PDB)
            else:
                # do it this way to suppress error messages
                try:
                    structure = get_structure(inputPath % PDB,PDB)
                except:
                    print("Could not load structure %s" % PDB)
                    continue





