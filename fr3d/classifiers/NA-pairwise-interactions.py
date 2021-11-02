# -*- coding: utf-8 -*-
"""
    This program reads one or more CIF files and produces annotations
    of nucleotide-nucleotide interactions.
    Basepairs are annotated with Leontis-Westhof annotations like cWW, tHS.
    Basepairs might also be annotated as "near" with ncWW, ntHS.
    A few basepair categories have "alternative" geometries like acWW, ctWW.
    Alternative geometries are not checked for hydrogen bonds.
"""

"""
    Developer notes:
        CC acHS needs a tighter gap requirement, avoid http://rna.bgsu.edu/rna3dhub/display3D/unitid/4V9F|1|0|C|2309,4V9F|1|0|C|2281
        AC cWW needs tigher requirements, avoid http://rna.bgsu.edu/rna3dhub/display3D/unitid/4V9F|1|0|A|2465,4V9F|1|0|C|2396
        AG cWS has near that should be true, see big cluster at http://rna.bgsu.edu/webfr3d/Results/60ca1ea55ae61/60ca1ea55ae61.html

    When fr3d is changed, python setup.py install
"""

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
from fr3d.modified_parent_mapping import modified_nucleotides

from discrepancy import matrix_discrepancy
import numpy as np
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
from fr3d.localpath import outputNAPickleInteractions
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

nt_nt_screen_distance = 12

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


def reverse_edges(inter):

    if inter[0] == "n":
        rev = inter[0] + inter[1] + inter[3] + inter[2] + inter[4:]
    else:
        rev = inter[0] + inter[2] + inter[1] + inter[3:]

    #print("reversed %s and got %s" % (inter,rev))

    return rev


def annotate_nt_nt_interactions(bases, screen_distance_cutoff, baseCubeList, baseCubeNeighbors):

    # loop through nt cubes, loop through neighboring nt cubes,
    # then loop through bases in the two cubes,
    # screening distances between them, then annotating interactions

    screen_distance_cutoff = 12    # minimum distance between base centers to screen for interaction

    count_pair = 0

    pair_to_interaction = {}
    interaction_to_pair_list = defaultdict(list)

    list_nt_nt = []
    contact_list = []
    output = []

    max_screen_distance = 0     # record the largest screening distance for which an interaction is found

    for nt1key in baseCubeList:                                # key to first cube
        for nt2key in baseCubeNeighbors[nt1key]:               # key to each potential neighboring cube, including the first
            if nt2key in baseCubeList:                         # if this cube was actually made
                for nt1 in baseCubeList[nt1key]:               # first nt of a potential pair
                    if len(nt1.centers["base"]) < 3:
                        print("  Missing base center for %s" % nt1.unit_id())
                        print(nt1.centers["base"])
                        continue

                    parent1 = get_parent(nt1.sequence)
                    gly1 = get_glycosidic_atom_coordinates(nt1,parent1)

                    if len(gly1) < 3:
                        print("  Missing glycosidic atom for %s" % nt1.unit_id())
                        continue

                    for nt2 in baseCubeList[nt2key]:           # second nt of a potential pair
                        if len(nt2.centers["base"]) < 3:
                            print("  Missing base center for %s" % nt2.unit_id())
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

                        interaction = check_base_backbone_O_stack(nt1,nt2)

                        if len(interaction) > 0:

                            # these interactions are currently experimental
                            # labeling them as such makes it possible to compare to previous ones
                            interaction += "_exp"

                            pair_to_interaction[(nt1.unit_id(),nt2.unit_id())] = interaction
                            interaction_to_pair_list[interaction].append((nt1.unit_id(),nt2.unit_id()))
                            max_screen_distance = max(max_screen_distance,screen_distance)

                        parent2 = get_parent(nt2.sequence)
                        parent_pair = parent1 + "," + parent2

#                       Note that AA, CC, GG, UU are being checked twice already!
#                       need to add T and have a plan for DNA nucleotides as well
                        if parent_pair in ['A,A','A,C','A,G','A,U','C,C','G,C','C,U','G,G','G,U','U,U']:
#                        if parent_pair in ['A,U']:
                            cutoffs = nt_nt_cutoffs[parent1+","+parent2]

                            gly2 = get_glycosidic_atom_coordinates(nt2,parent2)

                            if len(gly2) < 3:
                                print("  Missing glycosidic atom for %s" % nt2.unit_id())
                                continue


                            glycosidic_displacement = np.subtract(gly2,gly1)

                            pair_data = check_coplanar(nt1,nt2,glycosidic_displacement)

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
                                max_screen_distance = max(max_screen_distance,screen_distance)

                                inter = interaction[0][0]

                                # these interactions are currently experimental
                                # labeling them as such makes it possible to compare to previous ones
                                inter += "_exp"

                                pair_to_interaction[(nt1.unit_id(),nt2.unit_id())] = inter
                                interaction_to_pair_list[inter].append((nt1.unit_id(),nt2.unit_id()))

                                if interaction[0] in ["c","t"] or interaction in ["s33","s35","s53","s55"]:
                                    pair_to_interaction[(nt2.unit_id(),nt1.unit_id())] = reverse_edges(inter)

    print("  Found %d nucleotide-nucleotide pairs" % count_pair)
    print("  Maximum screen distance for actual contacts is %8.4f" % max_screen_distance)

    interaction_to_triple_list = calculate_crossing_numbers(bases,interaction_to_pair_list)

    return interaction_to_triple_list, pair_to_interaction, list_nt_nt

def calculate_crossing_numbers(bases,interaction_to_pair_list):
    # Identify which cWW pairs are nested
    # Then for each interaction, calculate the number of nested
    # cWW pairs it crosses

    # map unit_id to chain and sequence index
    unit_id_to_index = {}
    chain_to_max_index = defaultdict(lambda: 0)
    for nt in bases:
        unit_id = nt.unit_id()
        fields = unit_id.split("|")
        chain = fields[2]
        unit_id_to_index[unit_id] = (chain,nt.index)
        #print("%s\t%s" % (nt.index,nt.unit_id()))
        chain_to_max_index[chain] = max(chain_to_max_index[chain],nt.index)

    chain_to_cWW_pairs = defaultdict(list)   # separate list for each chain

    # find cWW basepairs within each chain
    for cWW in ['cWW_exp','cWw_exp','cwW_exp','acWW_exp','acWw_exp','acwW_exp']:
        for u1,u2 in interaction_to_pair_list[cWW]:
            chain1, index1 = unit_id_to_index[u1]
            chain2, index2 = unit_id_to_index[u2]

            # record cWW pairs within each chain
            if chain1 == chain2:
                if index1 < index2:
                    chain_to_cWW_pairs[chain1].append((index1,index2))
                else:
                    chain_to_cWW_pairs[chain1].append((index2,index1))

    chain_nested_cWW_endpoints = {}

    # within each chain, sort nested cWW by distance between them
    # started with the shortest-range pairs, record nested cWW pairs
    # by mapping one index to the other in chain_nested_cWW_endpoints
    for chain in chain_to_cWW_pairs.keys():
        cWW_pairs = sorted(chain_to_cWW_pairs[chain], key=lambda p: (p[1]-p[0],p[0]))

        nested_cWW_endpoints = []

        # at first, each index maps to itself
        for i in range(0,chain_to_max_index[chain]+1):
            nested_cWW_endpoints.append(i)

        # loop over cWW pairs and if no conflict, record as being nested
        for index1,index2 in cWW_pairs:
            # loop over indices within this pair, see if they map outside this pair
            i = index1+1
            while i < index2 and nested_cWW_endpoints[i] > index1 and nested_cWW_endpoints[i] < index2:
                i += 1

            if i == index2:
                # this pair is nested, so record the endpoints

                #print("index1",index1)
                #print("index2",index2)
                #print("list length",len(nested_cWW_endpoints))

                nested_cWW_endpoints[index1] = index2
                nested_cWW_endpoints[index2] = index1
            else:
                #print("cWW pair %s,%s is not nested" % (index1,index2))
                pass

        # record the nested cWW endpoints for this chain
        # might this only result in a pointer and so mixed up lists?
        chain_nested_cWW_endpoints[chain] = nested_cWW_endpoints


    #print('chain_to_max_index.keys()',chain_to_max_index.keys())
    #print('chain_to_cWW_pairs.keys()',chain_to_cWW_pairs.keys())
    #print('chain_nested_cWW_endpoints.keys()',chain_nested_cWW_endpoints.keys())
    interaction_to_triple_list = defaultdict(list)

    # loop over pairs, calculate crossing number
    # record interacting pairs and their crossing number as triples
    for interaction in interaction_to_pair_list.keys():

        for u1,u2 in interaction_to_pair_list[interaction]:
            chain1,index1 = unit_id_to_index[u1]
            chain2,index2 = unit_id_to_index[u2]

            crossing = 0

            # interactions within the same chain can have non-zero crossing number
            # some chains may not have any cWW pairs, then all interactions are nested
            if chain1 == chain2 and chain1 in chain_nested_cWW_endpoints:
                # put indices in increasing order
                index1,index2 = sorted([index1,index2])

                # count nested cWW that reach outside of [index1,index2]
                for i in range(index1+1,index2):

                    j = chain_nested_cWW_endpoints[chain1][i]

                    if j < index1 or j > index2:
                        crossing += 1

                if False and crossing > 0:
                    print("%-20s and %-20s make %s and have crossing number %d" % (u1,u2,interaction,crossing))

            interaction_to_triple_list[interaction].append((u1,u2,crossing))

            # duplicate certain pairs in reversed order
            if interaction[0] in ["c","t"] or interaction in ["s33","s35","s53","s55"]:
                interaction_to_triple_list[reverse_edges(interaction)].append((u2,u1,crossing))
            elif interaction[0:2] in ["nc","nt"] or interaction in ["ns33","ns35","ns53","ns55"]:
                interaction_to_triple_list[reverse_edges(interaction)].append((u2,u1,crossing))

    return interaction_to_triple_list

def annotate_nt_nt_in_structure(structure):
    """
    This function can be called from the pipeline to annotate a structure
    """

    bases = structure.residues(type = ["RNA linking","DNA linking"])  # load all RNA/DNA nucleotides
    print("  Building nucleotide cubes in " + PDB)
    baseCubeList, baseCubeNeighbors = make_nt_cubes(bases, nt_nt_screen_distance, nt_reference_point)

    # annotate nt-nt interactions
    print("  Annotating interactions")
    interaction_to_triple_list, pair_to_interaction, list_nt_nt = annotate_nt_nt_interactions(bases, nt_nt_screen_distance, baseCubeList, baseCubeNeighbors)

    return interaction_to_triple_list, pair_to_interaction


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

def check_base_backbone_O_stack(nt1,nt2):

    interaction = ""

    # rotate nt2 to origin
    standard_nt2 = nt1.translate_rotate_component(nt2)

    oxygens = ["O2'","O3'","O4'","O5'","OP1","OP2"]
    qmin = 6
    qmax = 0      # track to see when to bail out for large q

    for oxygen in oxygens:

        oxygen_atom = standard_nt2.centers[oxygen]

        if len(oxygen_atom) == 3:

            x = oxygen_atom[0]
            y = oxygen_atom[1]
            z = oxygen_atom[2]
            r = math.sqrt(x*x + y*y)
            q = x*x + y*y + 8*(abs(z)-2.9)**2    # measure of quality of location

            if q < qmin:
                qmin = q
                xmin = x
                ymin = y
                zmin = z
                rmin = r
                min_oxygen = oxygen

    # require r < 2 and abs(z)-2.9 < sqrt(1/2)=0.707, and elliptical between those two
    if qmin < 4:
        if z > 0:
            interaction = "s3" + min_oxygen
        else:
            interaction = "s5" + min_oxygen

    # require r < sqrt(6)=.4495 and abs(z)-2.9 < sqrt(6/9)=0.8165 and elliptical between
    elif qmin < 6:
        if z > 0:
            interaction = "ns3" + min_oxygen
        else:
            interaction = "ns5" + min_oxygen

    if False and len(interaction) > 0:
        print('%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t\t=hyperlink("http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s")' % (nt1.unit_id(),nt2.unit_id(),interaction,xmin,ymin,zmin,rmin,qmin,nt1.unit_id(),nt2.unit_id()))

    return interaction

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


def unit_vector(v):
    return v / np.linalg.norm(v)


def load_basepair_annotations(filename,all_pair_types):

    with open(filename,'rb') as opener:
        basepairs = pickle.load(opener)

    pair_to_interaction={}

    for bp_type in all_pair_types:
      for (u1,u2,c) in basepairs[bp_type]:
         pair_to_interaction[(u1,u2)]=bp_type

      near_bp_type="n"+bp_type
#      for (u1,u2,c) in basepairs[near_bp_type]:
#         pair_to_interaction[(u1,u2)]=near_bp_type

    return pair_to_interaction


def compare_annotations(pairs_1,pairs_2,all_pair_types,filename):

    allkeys = set(pairs_1.keys()) | set(pairs_1.keys())

    with open(filename, mode='w') as file:

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
            print("Missing residues")
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

PDB_List = ['5AJ3']
PDB_List = ['6hiv']
PDB_List = ['3QRQ','5J7L']
PDB_List = ['4V9F','4YBB','4Y4O','6AZ3','4P95']
PDB_List = ['3BT7']
PDB_List = ['5I4A']
PDB_List = ['6A2H']
PDB_List = ['3JB9']
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.72/3.0A/csv']
version = "_3.72_3.0"
PDB_List = ['1OCT']
PDB_List = ['4v9fFH.pdb']
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/2.5A/csv']
version = "_3.48_2.5"
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.74/4.0A/csv']
version = "_3.74_4.0"
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']
version = "_3.48_3.0"
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/view/NR_4.0_56726.45']
PDB_List = ['5KCR', '4WOI', '6C4I', '5JC9', '5L3P', '5KPW', '3J9Y', '3J9Z', '6BU8', '5WF0', '4V55', '4V54', '4V57', '4V56', '4V50', '4V53', '4V52', '4WF1', '5H5U', '4V5B', '5WFS', '5O2R', '5WFK', '5LZD', '5LZA', '6O9J', '6O9K', '6ORL', '6ORE', '3R8O', '3R8N', '4V85', '5MDV', '5MDW', '4V80', '4U27', '4U26', '4U25', '4U24', '4U20', '5KPS', '6GXM', '5KPX', '4U1U', '3JBU', '4V9P', '3JBV', '6Q9A', '6DNC', '4U1V', '6GXO', '5IQR', '5NWY', '4V9C', '6OSK', '4V9D', '4V9O', '5MGP', '6Q97', '3JCJ', '5J91', '3JCD', '3JCE', '6I7V', '6GXN', '4V64', '5J7L', '5AFI', '6BY1', '6ENU', '4V7V', '4V7U', '4V7T', '4V7S', '3JA1', '6ENF', '6OUO', '6ENJ', '5JU8', '5J8A', '6GWT', '4YBB', '5NP6', '5J88', '5U9G', '5U9F', '4V6D', '4V6E', '4V6C', '5JTE', '6OT3', '5J5B', '4WWW', '6OSQ', '5U4J', '5MDZ', '5U4I', '6NQB', '5UYQ', '5UYP', '5MDY', '5WDT', '6H4N', '5UYK', '4V89', '5UYM', '5UYL', '5UYN', '5WE6', '5WE4', '5KCS', '4V4Q', '4V4H', '5IT8']
PDB_List = ['4V51','4V9K']
PDB_List = ['6WJR']
PDB_List = ['6TPQ']
PDB_List = ['4KTG']
version = "_3.160_2.5"
PDB_List = ['5KCR']
PDB_List = ['7ECF']  # DNA quadruplex
PDB_List = ['4TNA']
PDB_List = ['5J7L']
PDB_List = ['4V9F','5J7L','4ARC']
PDB_List = ['4ARC']
PDB_List = ['4V9F']
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.201/2.5A/csv']
PDB_List = ['4ARC']
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.201/all/csv']

ReadPickleFile = True                  # when true, just read the .pickle file from a previous run
ReadPickleFile = False                 # when true, just read the .pickle file from a previous run

base_seq_list = ['DA','DT','DC','DG']  # for DNA
base_seq_list = ['A','U','C','G']      # for RNA
base_seq_list = ['A','U','C','G']      # for RNA
base_seq_list = []                     # for all nucleic acids, modified or not

nt_reference_point = "base"
atom_atom_min_distance = 5    # minimum distance between atoms in nts to consider them interacting

# plot one instance of each of the pairwise interactions
PlotPair = False
PlotPair = True
AlreadyPlotted = {}

ShowStructureReadingErrors = True
ShowStructureReadingErrors = False

unit_data_path = "C:/Users/zirbel/Documents/GitHub/fr3d-python/data/units"

# Just do the annotation, no diagnostics or anything
test_for_pipeline = True

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
            outputDataFilePickle = outputNAPickleInteractions + PDB + "_RNA_pairs_exp.pickle"

            if test_for_pipeline:

                if not os.path.exists(outputDataFilePickle):

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
                            print("  Could not load structure %s" % PDB)
                            continue

                    # write out data file of nucleotide centers and rotations that can be used by FR3D for searches
                    # need to be able to identify each chain that is available
                    # write_unit_data_file(PDB,unit_data_path,structure)

                    timerData = myTimer("Test for pipeline",timerData)

                    interaction_to_triple_list, pair_to_interaction = annotate_nt_nt_in_structure(structure)
                    print("  Annotated these interactions: %s" % interaction_to_triple_list.keys())



                    pickle.dump(interaction_to_triple_list,open(outputDataFilePickle,"wb"),2)

            else:

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

                # extract nucleotides to analyze
                IFE = PDB_IFE_Dict[PDB]          #
                if len(IFE) == 0:                # use the whole PDB file
                    if base_seq_list:
                        bases = structure.residues(sequence = base_seq_list)  # load just the types of bases in base_seq_list
                    else:
                        bases = structure.residues(type = ["RNA linking","DNA linking"])  # load all RNA/DNA nucleotides
                else:                            # use specific chains only
                    chain_ids = []
                    print("  Keeping only bases in chains %s" % IFE)
                    chains = IFE.split("+")
                    for chain in chains[1:]:            #skip element zero, leading +
                        fields = chain.split("|")
                        chain_ids.append(fields[2])
                    if base_seq_list:
                        bases = structure.residues(chain = chain_ids, sequence = base_seq_list)  # load just the types of bases in base_seq_list
                    else:
                        bases = structure.residues(chain = chain_ids)  # load all bases

                # ??? record which RNA/DNA chains are actually present
                # count nucleotides
                numBases = 0
                for base in bases:
                    numBases += 1
                print("  Found " + str(numBases) + " bases in " + PDB)

                # build cubes to be able to find potential pairs quickly
                timerData = myTimer("Building cubes",timerData)
                print("  Building nucleotide cubes in " + PDB)
                baseCubeList, baseCubeNeighbors = make_nt_cubes(bases, nt_nt_screen_distance, nt_reference_point)

                # annotate nt-nt interactions
                timerData = myTimer("Annotating interactions",timerData)
                Python_pairs, list_nt_nt = annotate_nt_nt_interactions(bases, nt_nt_screen_distance, baseCubeList, baseCubeNeighbors)

                timerData = myTimer("Recording interactions",timerData)





                # compare to previous annotations, may be machine specific

                pathAndFileName = "C:/Users/zirbel/Documents/FR3D/Python FR3D/data/pairs/%s_RNA_pairs.pickle" % PDB

                if not os.path.exists(pathAndFileName):
                    pairsFileName = PDB + '_RNA_pairs.pickle'
                    print("Downloading "+pairsFileName)
                    if sys.version_info[0] < 3:
                        urllib.urlretrieve("http://rna.bgsu.edu/pairs/"+pairsFileName, pathAndFileName)  # python 2
                    else:
                        urllib.request.urlretrieve("http://rna.bgsu.edu/pairs/"+pairsFileName, pathAndFileName)  # python 3


                all_pair_types= ['cWW', 'tWW', 'cWH', 'tWH', 'cWS', 'tWS', 'cHH', 'tHH', 'cHS', 'tHS', 'cSS', 'tSS']
                Matlab_pairs = load_basepair_annotations(pathAndFileName,all_pair_types)

                comparison_filename = "C:/Users/zirbel/Documents/FR3D/Python FR3D/comparison/comparison_%s.txt" % PDB
                compare_annotations(Matlab_pairs,Python_pairs,all_pair_types,comparison_filename)


                # write out pairs in the format that WebFR3D reads


                # accumulate list of interacting units by base, interaction type, and edges
                for nt1, nt2, interaction, edge, standard_aa, param in list_nt_nt:
                    base = base_residue.unit_id()
                    # skip symmetry operated instances; generally these are just duplicates anyway
                    if not "||||" in str(base):
                        aa = aa_residue.unit_id()
                        base_component = str(base).split("|")
                        aa_component = str(aa).split("|")
                        key = base_component[3]+"_"+aa_component[3]+"_"+interaction+"_"+edge
                        count_pair += 1
                        allInteractionDictionary[key].append((base,aa,interaction,edge,standard_aa,param))  # store tuples

                """ 3D plots of base-aa interactions
                for base, aa, interaction in list_nt_nt:
                    base_seq = base.sequence
                    aa= aa.sequence

                    draw_base(base_seq, ax)
                    draw_aa(aa, ax)

                    ax.set_xlabel('X Axis')
                    ax.set_ylabel('Y Axis')
                    ax.set_zlabel('Z Axis')
                    ax.set_xlim3d(10, -15)
                    ax.set_ylim3d(10, -15)
                    ax.set_zlim3d(10, -15)
                    plt.title('%s with ' % base_seq +'%s' % aa + ' %s' % "null")
                    plt.show()
                              """
                #accumulate a full list of resultant RNA-aa pairs
        #        result_nt_nt.extend(list_nt_nt)

                #writing out output files
                #text_output(result_nt_nt)

     #           csv_output(list_nt_nt)
    #            print("  Wrote output to " + outputNAPairwiseInteractions + PDB)

                myTimer("summary",timerData)

#        print("Recorded %d pairwise interactions" % count_pair)

        # when appropriate, write out HTML files

        """
        if len(PDB_IFE_Dict) > 100:
            print("Writing " + outputDataFile)
            timerData = myTimer("Writing HTML files",timerData)
            pickle.dump((allInteractionDictionary,allAATwoBaseDictionary,PDB_List),open(outputDataFile,"wb"))
            writeInteractionsHTML(allInteractionDictionary,outputHTML,version)
        """

myTimer("summary",timerData)
