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

from urllib import request
from fr3d.cif.reader import Cif
from fr3d.definitions import RNAconnections
from fr3d.definitions import NAbaseheavyatoms
from fr3d.definitions import NAbasehydrogens
from fr3d.definitions import nt_sugar
from fr3d.definitions import nt_phosphate
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
if sys.version_info[0] < 3:
    from urllib import urlopen
else:
    from urllib.request import urlopen 

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

fr3d_classification_version = 'v1'   # temporary, for changes here



def myTimer(state,data={}):

    # add elapsed time to the current state of the timer
    if "currentState" in data:
        currentState = data["currentState"]
        data[currentState] += time() - data["lastTime"]

    if state == "summary":
        total = 0.0
        for state in data["allStates"]:
            if not state == "lastTime" and not state == "currentState":
                total += data[state]

        print("Summary of time taken:")
        for state in data["allStates"]:
            if not state == "lastTime" and not state == "currentState":
                print("%-31s: %10.3f seconds %10.3f minutes %10.3f%% of total" % (state,data[state],data[state]/60,100*data[state]/total))

        print("%-31s: %10.3f seconds %10.3f minutes %10.3f%% of total" % ("Total",total,total/60,100))


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

    # if not available locally, download from PDB and save locally
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


def annotate_nt_nt_interactions(bases, center_center_distance_cutoff, baseCubeList, baseCubeNeighbors, timerData):
    """
    loop through nt cubes, loop through neighboring nt cubes,
    then loop through bases in the two cubes,
    screening distances between them, then annotating interactions
    """


    get_datapoint = True           # collect data about each pair to pass back

    count_pair = 0

    pair_to_interaction = {}       # note:  this will need to return a list, not a text string,
                                   # because one pair of nucleotides could make more than one interaction
    interaction_to_pair_list = defaultdict(list)

    pair_to_data = defaultdict(dict)

    list_nt_nt = []
    contact_list = []
    output = []

    max_center_center_distance = 0     # record the largest screening distance for which an interaction is found

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

                        # vector displacement between base centers
                        displacement = abs(nt2.centers["base"]-nt1.centers["base"]) # center-center

                        # quick screens for base centers being too far apart
                        if displacement[0] > center_center_distance_cutoff or \
                           displacement[1] > center_center_distance_cutoff or \
                           displacement[2] > center_center_distance_cutoff:
                            continue

                        # calculate actual center-center distance, screen
                        center_center_distance = np.linalg.norm(displacement)

                        # base centers are too far apart to interact, screen them out
                        if center_center_distance > center_center_distance_cutoff:
                            continue

                        # some structures have overlapping nucleotides, screen, those out
                        if center_center_distance < 2:
                            continue

                        # screen for the same nucleotide with alternate coordinates?
                        # There are very few of those, and they can be screened out in other places

#                        print("  Checking for an interaction between %-18s and %-18s center-center distance %7.4f" % (nt1.unit_id(),nt2.unit_id(),center_center_distance))

                        unit_id_pair = (nt1.unit_id(),nt2.unit_id())  # tuple for these nucleotides in this order

                        parent2 = get_parent(nt2.sequence)
                        parent_pair = parent1 + "," + parent2

                        if get_datapoint:
                            datapoint = {}
                            datapoint['center_center_distance'] = center_center_distance
                            datapoint['nt1_seq'] = nt1.sequence
                            datapoint['nt2_seq'] = nt2.sequence
                            datapoint['nt1_parent'] = parent1
                            datapoint['nt2_parent'] = parent2
                            datapoint['url'] = "http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt1.unit_id(),nt2.unit_id())
                        else:
                            datapoint = None

                        # check base-oxygen stack
                        timerData = myTimer("Check base oxygen stack",timerData)
                        interaction, datapoint = check_base_oxygen_stack_rings(nt1,nt2,parent1,datapoint)

                        if len(interaction) > 0:
                            pair_to_interaction[unit_id_pair] = interaction  # fragile, only stores one interaction per pair
                            interaction_to_pair_list[interaction].append(unit_id_pair)
                            max_center_center_distance = max(max_center_center_distance,center_center_distance)  # for setting optimally

                        # check coplanar and basepairing for bases in specific orders
                        # AA, CC, GG, UU will be checked in both nucleotide orders, that's OK
                        # need to add T and have a plan for DNA nucleotides as well
                        if parent_pair in ['A,A','A,C','A,G','A,U','C,C','G,C','C,U','G,G','G,U','U,U']:
                            gly2 = get_glycosidic_atom_coordinates(nt2,parent2)

                            if len(gly2) < 3:
                                print("  Missing glycosidic atom for %s" % nt2.unit_id())
                                continue

                            glycosidic_displacement = np.subtract(gly2,gly1)

                            timerData = myTimer("Get basepair parameters",timerData)

                            pair_data, datapoint = get_basepair_parameters(nt1,nt2,glycosidic_displacement,datapoint)

                            timerData = myTimer("Check basepairing",timerData)

                            cutoffs = nt_nt_cutoffs[parent1+","+parent2]
                            interaction, datapoint = check_basepair_cutoffs(nt1,nt2,pair_data,cutoffs,datapoint)

                            check_base_base_stacking(nt1, nt2, pair_data)

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
                                max_center_center_distance = max(max_center_center_distance,center_center_distance)

                                inter = interaction[0][0]

                                # these interactions are currently experimental
                                # labeling them as such makes it possible to compare to previous ones
                                inter += "_exp"

                                pair_to_interaction[unit_id_pair] = inter
                                interaction_to_pair_list[inter].append(unit_id_pair)

                                # TODO this should also be done for near interactions
                                if interaction[0] in ["c","t"] or interaction in ["s33","s35","s53","s55"]:
                                    pair_to_interaction[(nt2.unit_id(),nt1.unit_id())] = reverse_edges(inter)

                        pair_to_data[unit_id_pair] = datapoint

    print("  Found %d nucleotide-nucleotide pairs" % count_pair)
    print("  Maximum screen distance for actual contacts is %8.4f" % max_center_center_distance)

    # calculate and save crossing numbers for each annoated interaction
    timerData = myTimer("Calculate crossing",timerData)
    interaction_to_triple_list = calculate_crossing_numbers(bases,interaction_to_pair_list)

    return interaction_to_triple_list, pair_to_interaction, pair_to_data, timerData

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

def annotate_nt_nt_in_structure(structure,timerData):
    """
    This function can be called from the pipeline to annotate a structure
    structure is an output from
    """

    bases = structure.residues(type = ["RNA linking","DNA linking"])  # load all RNA/DNA nucleotides
    print("  Building nucleotide cubes in " + PDB)
    timerData = myTimer("Building cubes",timerData)
    baseCubeList, baseCubeNeighbors = make_nt_cubes(bases, nt_nt_screen_distance, nt_reference_point)

    # annotate nt-nt interactions
    print("  Annotating interactions")
    timerData = myTimer("Annotating interactions",timerData)
    interaction_to_triple_list, pair_to_interaction, pair_to_data, timerData = annotate_nt_nt_interactions(bases, nt_nt_screen_distance, baseCubeList, baseCubeNeighbors, timerData)

    return interaction_to_triple_list, pair_to_interaction, pair_to_data, timerData


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

def translate_rotate_point(nt,point):
    """
    Use the rotation matrix and center of nt to move point into standard position
    """

    translated_coord = np.subtract(point, nt.centers["base"])
    translated_coord_matrix = np.matrix(translated_coord)
    rotated_coord = translated_coord_matrix * nt.rotation_matrix
    coord_array = np.array(rotated_coord)
    a = coord_array.flatten()
    new_point = a.tolist()

    return new_point

"""
def check_base_oxygen_stack(nt1,nt2,get_datapoint=False):
    '''
    Does one of the backbone oxygens of nt2 stack on the base of nt1?
    '''

    interaction = ""

    oxygens = ["O2'","O3'","O4'","O5'","OP1","OP2"]

    qmin = 6      # q measures the quality of the interaction
    qmax = 0      # track to see when to bail out for large q

    for oxygen in oxygens:

        oxygen_point = nt2.centers[oxygen]

        if len(oxygen_point) == 3:

            x,y,z = translate_rotate_point(nt1,oxygen_point)

            r = math.sqrt(x*x + y*y)
            q = x*x + y*y + 8*(abs(z)-2.9)**2    # measure of quality of location

            if q < qmin:
                qmin = q
                xmin = x
                ymin = y
                zmin = z
                rmin = r
                oxygenmin = oxygen

    # require r < 2 and abs(z)-2.9 < sqrt(1/2)=0.707, and elliptical between those two
    if qmin < 4:
        if zmin > 0:
            interaction = "s3" + oxygenmin
        else:
            interaction = "s5" + oxygenmin

    # require r < sqrt(6)=.4495 and abs(z)-2.9 < sqrt(6/9)=0.8165 and elliptical between
    elif qmin < 6:
        if zmin > 0:
            interaction = "ns3" + oxygenmin
        else:
            interaction = "ns5" + oxygenmin

    if False and len(interaction) > 0:
        print('%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t\t=hyperlink("http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s")' % (nt1.unit_id(),nt2.unit_id(),interaction,xmin,ymin,zmin,rmin,qmin,nt1.unit_id(),nt2.unit_id()))


    if get_datapoint:
        datapoint = {}
        if len(interaction) > 0:
            datapoint['sOx'] = xmin
            datapoint['sOy'] = ymin
            datapoint['sOz'] = zmin
            datapoint['sOq'] = qmin
            datapoint['sOr'] = rmin
            datapoint['sOinteraction'] = interaction

    return interaction, datapoint
"""

def check_base_oxygen_stack_rings(nt1,nt2,parent1,datapoint):
    '''
    Does one of the backbone oxygens of nt2 stack inside a ring on the base of nt1?
    '''

    true_z_cutoff = 3.5
    near_z_cutoff = 3.6

    interaction = ""

    oxygens = ["O2'","O3'","O4'","O5'","OP1","OP2"]
    oxygen_points = []  # list of translated rotated points

    true_found = False
    near_found = False

    zmin = 999    # keep track of oxygen closest to the plane and over a ring

    for oxygen in oxygens:

        oxygen_point = nt2.centers[oxygen]

        if len(oxygen_point) == 3:   # avoid atoms with missing coordinates

            x,y,z = translate_rotate_point(nt1,oxygen_point)  # put into standard orientation

            oxygen_points.append([x,y,z,oxygen])  # store for checking near interactions later

            # exclude impossibly close stacking, for example, from alternate locations of nt atoms
            if abs(z) < 2:
                continue

            ring5 = False
            ring6 = False

            # check z component, then check if projected point is inside a ring
            if abs(z) < near_z_cutoff:
                if parent1 == 'A' or parent1 == 'DA':
                    if -1.302671*x + -0.512161*y + -0.512114 > 0:  # Left of C4-C5
                        if -0.014382*x + -1.379291*y +  0.382370 > 0:  # Left of C5-N7
                            if  1.286593*x + -0.316949*y +  2.517358 > 0:  # Left of N7-C8
                                if  0.833587*x +  1.089911*y +  2.912966 > 0:  # Left of C8-N9
                                    if -0.803127*x +  1.118490*y +  1.147479 > 0:  # Left of N9-C4
                                        ring5 = True
                    else:
                        if  0.363524*x +  1.290539*y +  1.313698 > 0:  # Left of C4-N3
                            if -1.076359*x +  0.793555*y +  2.495722 > 0:  # Left of N3-C2
                                if -1.308429*x + -0.337740*y +  2.633517 > 0:  # Left of C2-N1
                                    if -0.319116*x + -1.301200*y +  1.862429 > 0:  # Left of N1-C6
                                        if  1.037709*x + -0.957315*y +  0.793620 > 0:  # Left of C6-C5
                                            ring6 = True
                elif parent1 == 'C' or parent1 == 'DC':
                    if -0.599253*x +  1.289335*y +  1.686062 > 0:  # Left of N1-C2
                        if -1.378522*x +  0.022802*y +  1.272927 > 0:  # Left of C2-N3
                            if -0.676851*x + -1.128767*y +  1.187225 > 0:  # Left of N3-C4
                                if  0.596389*x + -1.312333*y +  1.653099 > 0:  # Left of C4-C5
                                    if  1.359882*x + -0.033090*y +  2.071781 > 0:  # Left of C5-C6
                                        if  0.698355*x +  1.162053*y +  1.990943 > 0:  # Left of C6-N1
                                            ring6 = True
                elif parent1 == 'G' or parent1 == 'DG':
                    if -1.306197*x + -0.492373*y + -0.896488 > 0:  # Left of C4-C5
                        if -0.023230*x + -1.376606*y +  0.510698 > 0:  # Left of C5-N7
                            if  1.278249*x + -0.337248*y +  2.960145 > 0:  # Left of N7-C8
                                if  0.841883*x +  1.088640*y +  3.089984 > 0:  # Left of C8-N9
                                    if -0.790705*x +  1.117587*y +  0.761380 > 0:  # Left of N9-C4
                                        ring5 = True
                    else:
                        if  0.449709*x +  1.286231*y +  1.337347 > 0:  # Left of C4-N3
                            if -0.992445*x +  0.855594*y +  2.112909 > 0:  # Left of N3-C2
                                if -1.324604*x + -0.362005*y +  2.250906 > 0:  # Left of C2-N1
                                    if -0.533023*x + -1.330285*y +  2.026599 > 0:  # Left of N1-C6
                                        if  1.094166*x + -0.941908*y +  1.272410 > 0:  # Left of C6-C5
                                            ring6 = True
                elif parent1 == 'U':
                    if -0.589251*x +  1.260286*y +  1.716262 > 0:  # Left of N1-C2
                        if -1.384641*x + -0.064970*y +  1.232961 > 0:  # Left of C2-N3
                            if -0.834465*x + -1.135313*y +  1.246706 > 0:  # Left of N3-C4
                                if  0.745842*x + -1.256133*y +  1.824059 > 0:  # Left of C4-C5
                                    if  1.352820*x +  0.018369*y +  2.049668 > 0:  # Left of C5-C6
                                        if  0.709695*x +  1.177761*y +  2.015286 > 0:  # Left of C6-N1
                                            ring6 = True
                elif parent1 == 'DT':
                    if -0.675137*x +  1.198579*y +  2.053967 > 0:  # Left of N1-C2
                        if -1.365448*x + -0.109817*y +  1.633725 > 0:  # Left of C2-N3
                            if -0.742906*x + -1.165341*y +  1.298813 > 0:  # Left of N3-C4
                                if  0.767749*x + -1.221287*y +  1.359137 > 0:  # Left of C4-C5
                                    if  1.338191*x +  0.092630*y +  1.600513 > 0:  # Left of C5-C6
                                        if  0.677551*x +  1.205236*y +  1.959719 > 0:  # Left of C6-N1
                                            ring6 = True

            if ring5 or ring6:
                if abs(z) < true_z_cutoff:
                    true_found = True
                else:
                    near_found = True

                if abs(z) < abs(zmin):       # better than any previous stacking
                    xmin = x
                    ymin = y
                    zmin = z
                    oxygenmin = oxygen
                    if ring5:
                        ringmin = "ring5"
                    else:
                        ringmin = "ring6"

    if true_found:  # over base ring and z value is OK
        if zmin > 0:
            interaction = "s3" + oxygenmin
        else:
            interaction = "s5" + oxygenmin

    elif near_found:  # over a base ring, but z value too large for true
        if zmin > 0:
            interaction = "ns3" + oxygenmin
        else:
            interaction = "ns5" + oxygenmin

    else:            # not over a base ring, but maybe close enough

        r2min = 999   # keep track of distance to base center, use the minimum

        for x,y,z,oxygen in oxygen_points:

            nearring5 = False
            nearring6 = False

            if abs(z) < true_z_cutoff:
                if parent1 == 'A' or parent1 == 'DA':
                    if -1.302671*x + -0.512161*y + -0.512114 > 0:  # Left of C4-C5
                        if -0.014382*x + -1.379291*y +  1.072053 > 0:  # Within  0.500000 Angstroms of being left of C5-N7
                            if  1.286593*x + -0.316949*y +  3.179887 > 0:  # Within  0.500000 Angstroms of being left of N7-C8
                                if  0.833587*x +  1.089911*y +  3.599037 > 0:  # Within  0.500000 Angstroms of being left of C8-N9
                                    if -0.803127*x +  1.118490*y +  1.835962 > 0:  # Within  0.500000 Angstroms of being left of N9-C4
                                        nearring5 = True
                    else:
                        if  0.363524*x +  1.290539*y +  1.984079 > 0:  # Within  0.500000 Angstroms of being left of C4-N3
                            if -1.076359*x +  0.793555*y +  3.164355 > 0:  # Within  0.500000 Angstroms of being left of N3-C2
                                if -1.308429*x + -0.337740*y +  3.309175 > 0:  # Within  0.500000 Angstroms of being left of C2-N1
                                    if -0.319116*x + -1.301200*y +  2.532309 > 0:  # Within  0.500000 Angstroms of being left of N1-C6
                                        if  1.037709*x + -0.957315*y +  1.499540 > 0:  # Within  0.500000 Angstroms of being left of C6-C5
                                            nearring6 = True
                elif parent1 == 'C' or parent1 == 'DC':
                    if -0.599253*x +  1.289335*y +  2.396957 > 0:  # Within  0.500000 Angstroms of being left of N1-C2
                        if -1.378522*x +  0.022802*y +  1.962283 > 0:  # Within  0.500000 Angstroms of being left of C2-N3
                            if -0.676851*x + -1.128767*y +  1.845298 > 0:  # Within  0.500000 Angstroms of being left of N3-C4
                                if  0.596389*x + -1.312333*y +  2.373845 > 0:  # Within  0.500000 Angstroms of being left of C4-C5
                                    if  1.359882*x + -0.033090*y +  2.751923 > 0:  # Within  0.500000 Angstroms of being left of C5-C6
                                        if  0.698355*x +  1.162053*y +  2.668820 > 0:  # Within  0.500000 Angstroms of being left of C6-N1
                                            nearring6 = True
                elif parent1 == 'G' or parent1 == 'DG':
                    if -1.306197*x + -0.492373*y + -0.896488 > 0:  # Left of C4-C5
                        if -0.023230*x + -1.376606*y +  1.199099 > 0:  # Within  0.500000 Angstroms of being left of C5-N7
                            if  1.278249*x + -0.337248*y +  3.621140 > 0:  # Within  0.500000 Angstroms of being left of N7-C8
                                if  0.841883*x +  1.088640*y +  3.778080 > 0:  # Within  0.500000 Angstroms of being left of C8-N9
                                    if -0.790705*x +  1.117587*y +  1.445889 > 0:  # Within  0.500000 Angstroms of being left of N9-C4
                                        nearring5 = True
                    else:
                        if  0.449709*x +  1.286231*y +  2.018638 > 0:  # Within  0.500000 Angstroms of being left of C4-N3
                            if -0.992445*x +  0.855594*y +  2.768079 > 0:  # Within  0.500000 Angstroms of being left of N3-C2
                                if -1.324604*x + -0.362005*y +  2.937496 > 0:  # Within  0.500000 Angstroms of being left of C2-N1
                                    if -0.533023*x + -1.330285*y +  2.743148 > 0:  # Within  0.500000 Angstroms of being left of N1-C6
                                        if  1.094166*x + -0.941908*y +  1.994280 > 0:  # Within  0.500000 Angstroms of being left of C6-C5
                                            nearring6 = True
                elif parent1 == 'DT':
                    if -0.675137*x +  1.198579*y +  2.741790 > 0:  # Within  0.500000 Angstroms of being left of N1-C2
                        if -1.365448*x + -0.109817*y +  2.318653 > 0:  # Within  0.500000 Angstroms of being left of C2-N3
                            if -0.742906*x + -1.165341*y +  1.989814 > 0:  # Within  0.500000 Angstroms of being left of N3-C4
                                if  0.767749*x + -1.221287*y +  2.080417 > 0:  # Within  0.500000 Angstroms of being left of C4-C5
                                    if  1.338191*x +  0.092630*y +  2.271210 > 0:  # Within  0.500000 Angstroms of being left of C5-C6
                                        if  0.677551*x +  1.205236*y +  2.651034 > 0:  # Within  0.500000 Angstroms of being left of C6-N1
                                            nearring6 = True
                elif parent1 == 'U':
                    if -0.589251*x +  1.260286*y +  2.411880 > 0:  # Within  0.500000 Angstroms of being left of N1-C2
                        if -1.384641*x + -0.064970*y +  1.926043 > 0:  # Within  0.500000 Angstroms of being left of C2-N3
                            if -0.834465*x + -1.135313*y +  1.951203 > 0:  # Within  0.500000 Angstroms of being left of N3-C4
                                if  0.745842*x + -1.256133*y +  2.554496 > 0:  # Within  0.500000 Angstroms of being left of C4-C5
                                    if  1.352820*x +  0.018369*y +  2.726141 > 0:  # Within  0.500000 Angstroms of being left of C5-C6
                                        if  0.709695*x +  1.177761*y +  2.702815 > 0:  # Within  0.500000 Angstroms of being left of C6-N1
                                            nearring6 = True

            """
            if abs(z) < true_z_cutoff:
                if parent1 == 'A' or parent1 == 'DA':
                    if -1.302671*x + -0.512161*y + -0.512114 > 0:  # Left of C4-C5
                        if -0.014382*x + -1.379291*y +  0.934116 > 0:  # Within  0.400000 Angstroms of being left of C5-N7
                            if  1.286593*x + -0.316949*y +  3.047381 > 0:  # Within  0.400000 Angstroms of being left of N7-C8
                                if  0.833587*x +  1.089911*y +  3.461823 > 0:  # Within  0.400000 Angstroms of being left of C8-N9
                                    if -0.803127*x +  1.118490*y +  1.698266 > 0:  # Within  0.400000 Angstroms of being left of N9-C4
                                        nearring5 = True
                    else:
                        if  0.363524*x +  1.290539*y +  1.850003 > 0:  # Within  0.400000 Angstroms of being left of C4-N3
                            if -1.076359*x +  0.793555*y +  3.030628 > 0:  # Within  0.400000 Angstroms of being left of N3-C2
                                if -1.308429*x + -0.337740*y +  3.174043 > 0:  # Within  0.400000 Angstroms of being left of C2-N1
                                    if -0.319116*x + -1.301200*y +  2.398333 > 0:  # Within  0.400000 Angstroms of being left of N1-C6
                                        if  1.037709*x + -0.957315*y +  1.358356 > 0:  # Within  0.400000 Angstroms of being left of C6-C5
                                            nearring6 = True
                elif parent1 == 'C' or parent1 == 'DC':
                    if -0.599253*x +  1.289335*y +  2.254778 > 0:  # Within  0.400000 Angstroms of being left of N1-C2
                        if -1.378522*x +  0.022802*y +  1.824412 > 0:  # Within  0.400000 Angstroms of being left of C2-N3
                            if -0.676851*x + -1.128767*y +  1.713684 > 0:  # Within  0.400000 Angstroms of being left of N3-C4
                                if  0.596389*x + -1.312333*y +  2.229696 > 0:  # Within  0.400000 Angstroms of being left of C4-C5
                                    if  1.359882*x + -0.033090*y +  2.615895 > 0:  # Within  0.400000 Angstroms of being left of C5-C6
                                        if  0.698355*x +  1.162053*y +  2.533245 > 0:  # Within  0.400000 Angstroms of being left of C6-N1
                                            nearring6 = True
                elif parent1 == 'G' or parent1 == 'DG':
                    if -1.306197*x + -0.492373*y + -0.896488 > 0:  # Left of C4-C5
                        if -0.023230*x + -1.376606*y +  1.061418 > 0:  # Within  0.400000 Angstroms of being left of C5-N7
                            if  1.278249*x + -0.337248*y +  3.488941 > 0:  # Within  0.400000 Angstroms of being left of N7-C8
                                if  0.841883*x +  1.088640*y +  3.640461 > 0:  # Within  0.400000 Angstroms of being left of C8-N9
                                    if -0.790705*x +  1.117587*y +  1.308988 > 0:  # Within  0.400000 Angstroms of being left of N9-C4
                                        nearring5 = True
                    else:
                        if  0.449709*x +  1.286231*y +  1.882380 > 0:  # Within  0.400000 Angstroms of being left of C4-N3
                            if -0.992445*x +  0.855594*y +  2.637045 > 0:  # Within  0.400000 Angstroms of being left of N3-C2
                                if -1.324604*x + -0.362005*y +  2.800178 > 0:  # Within  0.400000 Angstroms of being left of C2-N1
                                    if -0.533023*x + -1.330285*y +  2.599839 > 0:  # Within  0.400000 Angstroms of being left of N1-C6
                                        if  1.094166*x + -0.941908*y +  1.849906 > 0:  # Within  0.400000 Angstroms of being left of C6-C5
                                            nearring6 = True
                elif parent1 == 'DT':
                    if -0.675137*x +  1.198579*y +  2.604225 > 0:  # Within  0.400000 Angstroms of being left of N1-C2
                        if -1.365448*x + -0.109817*y +  2.181667 > 0:  # Within  0.400000 Angstroms of being left of C2-N3
                            if -0.742906*x + -1.165341*y +  1.851614 > 0:  # Within  0.400000 Angstroms of being left of N3-C4
                                if  0.767749*x + -1.221287*y +  1.936161 > 0:  # Within  0.400000 Angstroms of being left of C4-C5
                                    if  1.338191*x +  0.092630*y +  2.137070 > 0:  # Within  0.400000 Angstroms of being left of C5-C6
                                        if  0.677551*x +  1.205236*y +  2.512771 > 0:  # Within  0.400000 Angstroms of being left of C6-N1
                                            nearring6 = True
                elif parent1 == 'U':
                    if -0.589251*x +  1.260286*y +  2.272756 > 0:  # Within  0.400000 Angstroms of being left of N1-C2
                        if -1.384641*x + -0.064970*y +  1.787427 > 0:  # Within  0.400000 Angstroms of being left of C2-N3
                            if -0.834465*x + -1.135313*y +  1.810304 > 0:  # Within  0.400000 Angstroms of being left of N3-C4
                                if  0.745842*x + -1.256133*y +  2.408409 > 0:  # Within  0.400000 Angstroms of being left of C4-C5
                                    if  1.352820*x +  0.018369*y +  2.590846 > 0:  # Within  0.400000 Angstroms of being left of C5-C6
                                        if  0.709695*x +  1.177761*y +  2.565310 > 0:  # Within  0.400000 Angstroms of being left of C6-N1
                                            nearring6 = True
            """
            """
            if abs(z) < near_z_cutoff:
                if parent1 == 'A' or parent1 == 'DA':
                    if -1.302671*x + -0.512161*y + -0.512114 > 0:  # Left of C4-C5
                        if -0.014382*x + -1.379291*y +  0.658243 > 0:  # Within  0.200000 Angstroms of being left of C5-N7
                            if  1.286593*x + -0.316949*y +  2.782370 > 0:  # Within  0.200000 Angstroms of being left of N7-C8
                                if  0.833587*x +  1.089911*y +  3.187395 > 0:  # Within  0.200000 Angstroms of being left of C8-N9
                                    if -0.803127*x +  1.118490*y +  1.422873 > 0:  # Within  0.200000 Angstroms of being left of N9-C4
                                        nearring5 = True
                    else:
                        if  0.363524*x +  1.290539*y +  1.581851 > 0:  # Within  0.200000 Angstroms of being left of C4-N3
                            if -1.076359*x +  0.793555*y +  2.763175 > 0:  # Within  0.200000 Angstroms of being left of N3-C2
                                if -1.308429*x + -0.337740*y +  2.903780 > 0:  # Within  0.200000 Angstroms of being left of C2-N1
                                    if -0.319116*x + -1.301200*y +  2.130381 > 0:  # Within  0.200000 Angstroms of being left of N1-C6
                                        if  1.037709*x + -0.957315*y +  1.075988 > 0:  # Within  0.200000 Angstroms of being left of C6-C5
                                            nearring6 = True
                elif parent1 == 'C' or parent1 == 'DC':
                    if -0.599253*x +  1.289335*y +  1.970420 > 0:  # Within  0.200000 Angstroms of being left of N1-C2
                        if -1.378522*x +  0.022802*y +  1.548670 > 0:  # Within  0.200000 Angstroms of being left of C2-N3
                            if -0.676851*x + -1.128767*y +  1.450454 > 0:  # Within  0.200000 Angstroms of being left of N3-C4
                                if  0.596389*x + -1.312333*y +  1.941398 > 0:  # Within  0.200000 Angstroms of being left of C4-C5
                                    if  1.359882*x + -0.033090*y +  2.343838 > 0:  # Within  0.200000 Angstroms of being left of C5-C6
                                        if  0.698355*x +  1.162053*y +  2.262094 > 0:  # Within  0.200000 Angstroms of being left of C6-N1
                                            nearring6 = True
                elif parent1 == 'G' or parent1 == 'DG':
                    if -1.306197*x + -0.492373*y + -0.896488 > 0:  # Left of C4-C5
                        if -0.023230*x + -1.376606*y +  0.786058 > 0:  # Within  0.200000 Angstroms of being left of C5-N7
                            if  1.278249*x + -0.337248*y +  3.224543 > 0:  # Within  0.200000 Angstroms of being left of N7-C8
                                if  0.841883*x +  1.088640*y +  3.365222 > 0:  # Within  0.200000 Angstroms of being left of C8-N9
                                    if -0.790705*x +  1.117587*y +  1.035184 > 0:  # Within  0.200000 Angstroms of being left of N9-C4
                                        nearring5 = True
                    else:
                        if  0.449709*x +  1.286231*y +  1.609864 > 0:  # Within  0.200000 Angstroms of being left of C4-N3
                            if -0.992445*x +  0.855594*y +  2.374977 > 0:  # Within  0.200000 Angstroms of being left of N3-C2
                                if -1.324604*x + -0.362005*y +  2.525542 > 0:  # Within  0.200000 Angstroms of being left of C2-N1
                                    if -0.533023*x + -1.330285*y +  2.313219 > 0:  # Within  0.200000 Angstroms of being left of N1-C6
                                        if  1.094166*x + -0.941908*y +  1.561158 > 0:  # Within  0.200000 Angstroms of being left of C6-C5
                                            nearring6 = True
                elif parent1 == 'DT':
                    if -0.675137*x +  1.198579*y +  2.329096 > 0:  # Within  0.200000 Angstroms of being left of N1-C2
                        if -1.365448*x + -0.109817*y +  1.907696 > 0:  # Within  0.200000 Angstroms of being left of C2-N3
                            if -0.742906*x + -1.165341*y +  1.575214 > 0:  # Within  0.200000 Angstroms of being left of N3-C4
                                if  0.767749*x + -1.221287*y +  1.647649 > 0:  # Within  0.200000 Angstroms of being left of C4-C5
                                    if  1.338191*x +  0.092630*y +  1.868792 > 0:  # Within  0.200000 Angstroms of being left of C5-C6
                                        if  0.677551*x +  1.205236*y +  2.236245 > 0:  # Within  0.200000 Angstroms of being left of C6-N1
                                            nearring6 = True
                elif parent1 == 'U':
                    if -0.589251*x +  1.260286*y +  1.994509 > 0:  # Within  0.200000 Angstroms of being left of N1-C2
                        if -1.384641*x + -0.064970*y +  1.510194 > 0:  # Within  0.200000 Angstroms of being left of C2-N3
                            if -0.834465*x + -1.135313*y +  1.528505 > 0:  # Within  0.200000 Angstroms of being left of N3-C4
                                if  0.745842*x + -1.256133*y +  2.116234 > 0:  # Within  0.200000 Angstroms of being left of C4-C5
                                    if  1.352820*x +  0.018369*y +  2.320257 > 0:  # Within  0.200000 Angstroms of being left of C5-C6
                                        if  0.709695*x +  1.177761*y +  2.290298 > 0:  # Within  0.200000 Angstroms of being left of C6-N1
                                            nearring6 = True
            """

            if nearring5 or nearring6:
                near_found = True

                r2 = x**2 + y**2

                if r2 < r2min:       # closer to the base center than other near interactions
                    r2min = r2
                    xmin = x
                    ymin = y
                    zmin = z
                    oxygenmin = oxygen
                    if nearring5:
                        ringmin = "ring5"
                    else:
                        ringmin = "ring6"

        if near_found:
            if zmin > 0:
                interaction = "ns3" + oxygenmin
            else:
                interaction = "ns5" + oxygenmin

    if False and len(interaction) > 0:
        print('%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\t\t=hyperlink("http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s")' % (nt1.unit_id(),nt2.unit_id(),interaction,xmin,ymin,zmin,nt1.unit_id(),nt2.unit_id()))

    if datapoint:
        if len(interaction) > 0:
            datapoint['sOx'] = xmin
            datapoint['sOy'] = ymin
            datapoint['sOz'] = zmin
            datapoint['sOring'] = ringmin
            datapoint['sOoxygen'] = oxygenmin
            datapoint['sOinteraction'] = interaction

    return interaction, datapoint

def check_base_base_stacking(nt1, nt2, pair_data):
    """Pass in two nucleotides. Base stacking."""
    true_z_cutoff = 3.5 #maybe this should be 3.4? 
    near_z_cutoff = 3.6
    interaction = "" 
    true_found = False
    near_found = False
    #Outermost Atoms Of NT Bases
    convexGAtoms = ['N2','N1','O6','N7','C8','N9','N3']
    convexCAtoms = ['N4','C5','C6','N1','O2','N3']
    convexAAtoms = ['N6','N1','C2','N3','N9','C8','N7']
    convexUAtoms = ['O4','N3','O2','N1','C6','C5']
    #Add Atoms for DT

    nt1BaseZcenter = nt1.centers["base"][2] #Dont think I need this.
    nt2BaseZcenter = nt2.centers["base"][2] #z coordinates for center base
    #print(list(nt1.sequence)) #Sequence holds the nt type A C G or U 
    x, y, z  = translate_rotate_point(nt1, nt2.centers["base"]) #Puts into standard orientation?
    #nt1normal = normal_vector_calculation(nt1) #This function takes in an entire residue...Does this normal vector actually represent the normal vector from the base? There might be useful stuff in the check coplanar function
    #nt2normal = normal_vector_calculation(nt2) #Dont think I actually need these either...
    #print("Z" + str(z))
    #draw_base(nt1,2)
    #draw_base(nt2,2)
    #print("NT2 Center Contents: ")
    #print(nt2.centers)
    #print("NT2 ATOMS:")
    #print(nt2._atoms)
    #for atom in nt2.atoms["base"]:
    #    print(str(atom))
    #print("PAIR DATA: ")
    #print(pair_data)
    #displ12 vector from origin to nt2 when standardized
    #for atom in nt2.atoms()
    #Check the pair data to see if nt1 and nt2 are coplanar. If they are not coplanar, they are not stacking. 
    for atom in convexAAtoms:
        x, y, z  = translate_rotate_point(nt1, nt2.centers[convexAAtoms]) #Puts into standard orientation?
        print("I made it into for loop")
        print(nt2.sequence)
        print(atom)
        print(z)
        if abs(z) < near_z_cutoff:
            if  nt1.sequence == 'A':
                if  0.439603*x + -2.409029*y + -3.904243 > 0:  # Left of N3-N9
                    if -0.833587*x + -1.089911*y + -2.912966 > 0:  # Left of N9-C8
                            if -1.286593*x +  0.316949*y + -2.517358 > 0:  # Left of C8-N7
                                    if -2.326854*x +  1.976066*y + -4.969193 > 0:  # Left of N7-N6
                                            if  1.622643*x +  1.661740*y + -4.510171 > 0:  # Left of N6-N1
                                                    if  1.308429*x +  0.337740*y + -2.633517 > 0:  # Left of N1-C2
                                                            if  1.076359*x + -0.793555*y + -2.495722 > 0:  # Left of C2-N3
                                                                print("INSIDE FOR LOOP IN SIDE END OF A ") 
    if pair_data["coplanar"] == True:
        #print("Nt1 and Nt2Are Coplanar")
        garbage =3
    else:
        return None
    #Checks ti see if the atoms from nt2 are within the convex hull of nt1     
    if abs(z) < near_z_cutoff:
        if  nt1.sequence == 'A':
            if  0.439603*x + -2.409029*y + -3.904243 > 0:  # Left of N3-N9
                if -0.833587*x + -1.089911*y + -2.912966 > 0:  # Left of N9-C8
                        if -1.286593*x +  0.316949*y + -2.517358 > 0:  # Left of C8-N7
                                if -2.326854*x +  1.976066*y + -4.969193 > 0:  # Left of N7-N6
                                        if  1.622643*x +  1.661740*y + -4.510171 > 0:  # Left of N6-N1
                                                if  1.308429*x +  0.337740*y + -2.633517 > 0:  # Left of N1-C2
                                                        if  1.076359*x + -0.793555*y + -2.495722 > 0:  # Left of C2-N3
                                                            print("A IM HERE") 
        elif nt1.sequence == 'C':
            if  2.031427*x +  1.030781*y + -2.400765 > 0:  # Left of N4-N3
                    if  2.098558*x +  0.957313*y + -2.427068 > 0:  # Left of N3-O2
                            if -0.120783*x + -2.269450*y + -3.415154 > 0:  # Left of O2-N1
                                    if -0.698355*x + -1.162053*y + -1.990943 > 0:  # Left of N1-C6
                                            if -1.359882*x +  0.033090*y + -2.071781 > 0:  # Left of C6-C5
                                                    if -1.950965*x +  1.410319*y + -3.754099 > 0:  # Left of C5-N4
                                                        something = 5
        elif nt1.sequence == 'DT':
            if  2.031505*x +  1.412190*y +  3.691439 > 0:  # Left of C7-C6
                    if  0.677551*x +  1.205236*y +  1.959719 > 0:  # Left of C6-N1
                            if -0.116119*x +  2.281277*y +  3.818305 > 0:  # Left of N1-O2
                                    if -1.924466*x + -1.192515*y +  2.687361 > 0:  # Left of O2-N3
                                            if -1.969872*x + -1.161301*y +  2.728760 > 0:  # Left of N3-O4
                                                    if  1.301401*x + -2.544887*y +  5.949759 > 0:  # Left of O4-C7
                                                        something = 5
        elif nt1.sequence == 'G':
            if  1.736168*x +  1.500421*y + -3.920978 > 0:  # Left of O6-N1
                    if  1.594663*x +  1.698560*y + -3.904578 > 0:  # Left of N1-N2
                            if  0.722386*x + -2.192149*y + -3.689399 > 0:  # Left of N2-N3
                                    if  0.340996*x + -2.403818*y + -3.618345 > 0:  # Left of N3-N9
                                            if -0.841883*x + -1.088640*y + -3.089984 > 0:  # Left of N9-C8
                                                    if -1.278249*x +  0.337248*y + -2.960145 > 0:  # Left of C8-N7
                                                            if -2.274081*x +  2.148378*y + -5.898397 > 0:  # Left of N7-O6
                                                                something = 5
        elif nt1.sequence == 'U':
            if  1.957545*x + -1.373286*y +  3.733149 > 0:  # Left of O4-C5
                    if  1.352820*x +  0.018369*y +  2.049668 > 0:  # Left of C5-C6
                            if  0.709695*x +  1.177761*y +  2.015286 > 0:  # Left of C6-N1
                                    if  0.048394*x +  2.292490*y +  3.487594 > 0:  # Left of N1-O2
                                            if -2.022286*x + -1.097174*y +  2.261275 > 0:  # Left of O2-N3
                                                    if -2.046168*x + -1.018160*y +  2.245721 > 0:  # Left of N3-O4"""
                                                        something = 5
    """Code to check that an (x,y) point is close to being inside a base ring
    A
        if  0.439603*x + -2.409029*y + -2.679838 > 0:  # Within  0.500000 Angstroms of being left of N3-N9
                if -0.833587*x + -1.089911*y + -2.226895 > 0:  # Within  0.500000 Angstroms of being left of N9-C8
                        if -1.286593*x +  0.316949*y + -1.854829 > 0:  # Within  0.500000 Angstroms of being left of C8-N7
                                if -2.326854*x +  1.976066*y + -3.442834 > 0:  # Within  0.500000 Angstroms of being left of N7-N6
                                        if  1.622643*x +  1.661740*y + -3.348884 > 0:  # Within  0.500000 Angstroms of being left of N6-N1
                                                if  1.308429*x +  0.337740*y + -1.957859 > 0:  # Within  0.500000 Angstroms of being left of N1-C2
                                                        if  1.076359*x + -0.793555*y + -1.827089 > 0:  # Within  0.500000 Angstroms of being left of C2-N3
C
        if  2.031427*x +  1.030781*y + -1.261774 > 0:  # Within  0.500000 Angstroms of being left of N4-N3
                if  2.098558*x +  0.957313*y + -1.273769 > 0:  # Within  0.500000 Angstroms of being left of N3-O2
                        if -0.120783*x + -2.269450*y + -2.278823 > 0:  # Within  0.500000 Angstroms of being left of O2-N1
                                if -0.698355*x + -1.162053*y + -1.313067 > 0:  # Within  0.500000 Angstroms of being left of N1-C6
                                        if -1.359882*x +  0.033090*y + -1.391639 > 0:  # Within  0.500000 Angstroms of being left of C6-C5
                                                if -1.950965*x +  1.410319*y + -2.550431 > 0:  # Within  0.500000 Angstroms of being left of C5-N4
DT
        if  2.031505*x +  1.412190*y +  4.928502 > 0:  # Within  0.500000 Angstroms of being left of C7-C6
                if  0.677551*x +  1.205236*y +  2.651034 > 0:  # Within  0.500000 Angstroms of being left of C6-N1
                        if -0.116119*x +  2.281277*y +  4.960420 > 0:  # Within  0.500000 Angstroms of being left of N1-O2
                                if -1.924466*x + -1.192515*y +  3.819357 > 0:  # Within  0.500000 Angstroms of being left of O2-N3
                                        if -1.969872*x + -1.161301*y +  3.872112 > 0:  # Within  0.500000 Angstroms of being left of N3-O4
                                                if  1.301401*x + -2.544887*y +  7.378928 > 0:  # Within  0.500000 Angstroms of being left of O4-C7
G
        if  1.736168*x +  1.500421*y + -2.773639 > 0:  # Within  0.500000 Angstroms of being left of O6-N1
                if  1.594663*x +  1.698560*y + -2.739669 > 0:  # Within  0.500000 Angstroms of being left of N1-N2
                        if  0.722386*x + -2.192149*y + -2.535345 > 0:  # Within  0.500000 Angstroms of being left of N2-N3
                                if  0.340996*x + -2.403818*y + -2.404403 > 0:  # Within  0.500000 Angstroms of being left of N3-N9
                                        if -0.841883*x + -1.088640*y + -2.401888 > 0:  # Within  0.500000 Angstroms of being left of N9-C8
                                                if -1.278249*x +  0.337248*y + -2.299150 > 0:  # Within  0.500000 Angstroms of being left of C8-N7
                                                        if -2.274081*x +  2.148378*y + -4.334190 > 0:  # Within  0.500000 Angstroms of being left of N7-O6
U
        if  1.957545*x + -1.373286*y +  4.928755 > 0:  # Within  0.500000 Angstroms of being left of O4-C5
                if  1.352820*x +  0.018369*y +  2.726141 > 0:  # Within  0.500000 Angstroms of being left of C5-C6
                        if  0.709695*x +  1.177761*y +  2.702815 > 0:  # Within  0.500000 Angstroms of being left of C6-N1
                                if  0.048394*x +  2.292490*y +  4.634094 > 0:  # Within  0.500000 Angstroms of being left of N1-O2
                                        if -2.022286*x + -1.097174*y +  3.411648 > 0:  # Within  0.500000 Angstroms of being left of O2-N3
                                                if -2.046168*x + -1.018160*y +  3.388465 > 0:  # Within  0.500000 Angstroms of being left of N3-O4"""
    #print("Normal Vector: " + str(nt1normal))
    #print("Normal Vector NT2: " + str(nt2normal))
    #angle = angle_between_vectors(nt1normal, nt2normal)
    #print("ANGLE: " + str(angle))

    #Create a screen

    distanceBetweenZs = abs(nt1BaseZcenter-nt2BaseZcenter) #should get distance from the z plane 
    for atom in nt2._atoms:
        x = atom.coordinates()[0]
        y = atom.coordinates()[1]
        if distanceBetweenZs < near_z_cutoff:
            if atom == 'A' or atom == 'DA':
               return 'A'
            elif atom == 'C' or atom == 'DC':
                return 'C'
            elif atom == 'G' or atom == 'DG':
               return 'U'                   
            elif atom == 'U':
                return 'UT'
            elif atom == 'DT':
                return 'DT'

    if distanceBetweenZs < near_z_cutoff:
        s = 2 #Garbage. Should be used after logic for if point is in a plane.
    elif distanceBetweenZs < true_z_cutoff:
        s = 2 #Garbage. Should be used after logic for if point is in a plane.
    else:
        return None

    #print("NT1 CENTER: " + str(nt1BaseZcenter))
    #print("NT2 CENTER: " + str(nt2BaseZcenter)) #Z coordinate for nucleotide bases

    #print("Atoms in nt1 : ")
    #for atom in nt2._atoms:
        #print("coordinates: ")
       # print(atom.coordinates())
       # print(atom.name)

    """ Not using this yet but annotation for stacking on 3 prime face, stacking on 5 prime face, or near stack on 3, near stack on 5
        Logic need added for the other one too
    if true_found:
        if zmin > 0:
            interaction = "s3" 
        else:
            interaction = "s5" 

    elif near_found:
        if zmin > 0:
            interaction = "ns3" 
        else:
            interaction = "ns5"
    """


    return None

def get_basepair_parameters(nt1,nt2,glycosidic_displacement,datapoint):
    """
    Calculate data needed for basepair classification.
    Also, check specific criteria to say that
    they are enough in the same plane to be called coplanar.
    If so, return a number from 0 to 1 to measure the
    degree of coplanarity with 1 being best.
    Criteria for being coplanar or near coplanar:
      Pair.Gap must be < 97th percentile among basepairs (1.5179 Angstroms)
      min_distance must be < 97th percentile among basepairs (2.4589 A)
      Angle between center-center vector and normals must be > 70.2388 degrees
      Angle between normal vectors must be < 39.1315 degrees
    """

    pair_data = {}
    pair_data["coplanar"] = False
    pair_data["coplanar_value"] = -1         # 0 to 1 is coplanar, 1 is the best

    # vector from origin to nt2 when standardized
    displ12 = np.matmul(glycosidic_displacement,nt1.rotation_matrix)

    pair_data["displ12"] = displ12

    # calculate gap and standardized atoms from nt2
    gap12, points2 = calculate_basepair_gap(nt1,nt2)
    pair_data["gap12"] = gap12

    if datapoint:
        datapoint['x'] = displ12[0,0]
        datapoint['y'] = displ12[0,1]
        datapoint['z'] = displ12[0,2]
        datapoint['gap12'] = gap12

    if gap12 >= 1.5179:
        return pair_data, datapoint

    # calculate minimum distance between nt1 base and all atoms of nt2
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

    if datapoint:
        datapoint['min_distance'] = min_distance

    # modified nucleotides don't have hydrogens, so be more flexible with them
    if min_distance >= 3.4589:
        return pair_data, datapoint

    # if working with regular bases, insist on close contact
    if min_distance >= 2.4589 and nt1.sequence in ['A','C','G','U'] and nt2.sequence in ['A','C','G','U']:
        return pair_data, datapoint

    center_displ = np.subtract(nt1.centers["base"],nt2.centers["base"])
    center_displ = center_displ / np.linalg.norm(center_displ) # normalize

    # calculate angle between center_displ and normal vectors to bases
    dot1 = abs(np.matmul(center_displ,nt1.rotation_matrix[:,2]))[0,0]
    if dot1 >= 0.3381:
        return pair_data, datapoint

    dot2 = abs(np.matmul(center_displ,nt2.rotation_matrix[:,2]))[0,0]
    if dot2 >= 0.3381:
        return pair_data, datapoint

    # calculate angle between normal vectors to the bases
    dot3 = abs(np.matmul(nt1.rotation_matrix[:,2].T,nt2.rotation_matrix[:,2]))
    if dot3 <= 0.7757:
        return pair_data, datapoint

    gap21, points1 = calculate_basepair_gap(nt1,nt2,points1)

    if datapoint:
        datapoint['gap21'] = gap21

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

    if datapoint:
        datapoint['coplanar'] = pair_data['coplanar']
        datapoint['coplanar_value'] = pair_data['coplanar_value']

    return pair_data, datapoint

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


def check_basepair_cutoffs(nt1,nt2,pair_data,cutoffs,datapoint):
    """ Given nt1 and nt2 and the dictionary of cutoffs
    for that pair of nucleotides, check cutoffs for each
    basepair interaction type
    """

    displ = pair_data["displ12"]  # vector from origin to nt2 when standardized

    if abs(displ[0,2]) > 3.6:                            # too far out of plane
        return [], datapoint

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
        return [], datapoint

#    return ok_displacement_screen

    rotation_1_to_2 = np.matmul(np.transpose(nt1.rotation_matrix), nt2.rotation_matrix)

    normal_Z = rotation_1_to_2[2,2]   # z component of normal vector to second base

#    print("%s with %s normal_Z is %0.4f" % (nt1.unit_id(),nt2.unit_id(),normal_Z))

    if datapoint:
        datapoint['normal_Z'] = normal_Z

    ok_normal = []

    for interaction,subcategory in ok_displacement_screen:
        cut = cutoffs[interaction][subcategory]
        if normal_Z < cut['normalmin']:
            continue
        if normal_Z > cut['normalmax']:
            continue
        ok_normal.append((interaction,subcategory))

    if len(ok_normal) == 0:
        return [], datapoint

    angle_in_plane = math.atan2(rotation_1_to_2[1,1],rotation_1_to_2[1,0])*57.29577951308232 - 90

    if angle_in_plane <= -90:
        angle_in_plane += 360

    if datapoint:
        datapoint['angle_in_plane'] = angle_in_plane

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
        return [], datapoint

    # (interaction,subcategory) pairs that meet all requirements so far
    ok_gap = []

    for interaction,subcategory in ok_angle_in_plane:
        cut = cutoffs[interaction][subcategory]
        if cut['gapmax'] > 0.1:
            if pair_data["gap12"] > cut['gapmax']:
                continue

        ok_gap.append((interaction,subcategory))

    # deal with the possibility of multiple matching interactions
    if len(ok_gap) == 0:
        return ok_gap, datapoint
    elif len(ok_gap) > 1:
        if len(set([i for i,s in ok_gap])) > 1:
            print("Multiple basepair types for %s, using the first one" % datapoint['url'])
            print(ok_gap)
            print(datapoint)

    if datapoint:
        datapoint['basepair'] = ok_gap[0][0]
        datapoint['basepair_subcategory'] = ok_gap[0][1]

    # return just the first interaction type and subcategory
    return [ok_gap[0][0],ok_gap[0][1]], datapoint


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

def map_PDB_list_to_PDB_IFE_dict(PDB_list):
    """
    map a list of PDB ids or IFEs or URLs to a dictionary whose keys
    are PDB ids and whose values are IFEs in that PDB
    """

    PDB_IFE_Dict = defaultdict(str)   # accumulate PDB-IFE pairs
    for PDB in PDB_list:
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
                    newPDB = fields[1][1:5]   # use only PDB identifier, ignore IFE for now
                    PDB_IFE_Dict[newPDB] += "+" + fields[1].replace('"','')

        elif "+" in PDB:                      # in case multiple chains in an IFE
            newPDB = PDB.split("|")[0]        # in case model, chain is indicated
            PDB_IFE_Dict[newPDB] = PDB

        elif "|" in PDB:                      # in case model, chain is indicated
            newPDB = PDB.split("|")[0]
            PDB_IFE_Dict[newPDB] = PDB

        else:
            PDB_IFE_Dict[PDB] = ""            # indicates to process the whole PDB file

    return PDB_IFE_Dict


#=======================================================================

PDB_list = ['5AJ3']
PDB_list = ['6hiv']
PDB_list = ['3QRQ','5J7L']
PDB_list = ['4V9F','4YBB','4Y4O','6AZ3','4P95']
PDB_list = ['3BT7']
PDB_list = ['5I4A']
PDB_list = ['6A2H']
PDB_list = ['3JB9']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.72/3.0A/csv']
version = "_3.72_3.0"
PDB_list = ['1OCT']
PDB_list = ['4v9fFH.pdb']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/2.5A/csv']
version = "_3.48_2.5"
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.74/4.0A/csv']
version = "_3.74_4.0"
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']
version = "_3.48_3.0"
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/view/NR_4.0_56726.45']
PDB_list = ['5KCR', '4WOI', '6C4I', '5JC9', '5L3P', '5KPW', '3J9Y', '3J9Z', '6BU8', '5WF0', '4V55', '4V54', '4V57', '4V56', '4V50', '4V53', '4V52', '4WF1', '5H5U', '4V5B', '5WFS', '5O2R', '5WFK', '5LZD', '5LZA', '6O9J', '6O9K', '6ORL', '6ORE', '3R8O', '3R8N', '4V85', '5MDV', '5MDW', '4V80', '4U27', '4U26', '4U25', '4U24', '4U20', '5KPS', '6GXM', '5KPX', '4U1U', '3JBU', '4V9P', '3JBV', '6Q9A', '6DNC', '4U1V', '6GXO', '5IQR', '5NWY', '4V9C', '6OSK', '4V9D', '4V9O', '5MGP', '6Q97', '3JCJ', '5J91', '3JCD', '3JCE', '6I7V', '6GXN', '4V64', '5J7L', '5AFI', '6BY1', '6ENU', '4V7V', '4V7U', '4V7T', '4V7S', '3JA1', '6ENF', '6OUO', '6ENJ', '5JU8', '5J8A', '6GWT', '4YBB', '5NP6', '5J88', '5U9G', '5U9F', '4V6D', '4V6E', '4V6C', '5JTE', '6OT3', '5J5B', '4WWW', '6OSQ', '5U4J', '5MDZ', '5U4I', '6NQB', '5UYQ', '5UYP', '5MDY', '5WDT', '6H4N', '5UYK', '4V89', '5UYM', '5UYL', '5UYN', '5WE6', '5WE4', '5KCS', '4V4Q', '4V4H', '5IT8']
PDB_list = ['4V51','4V9K']
PDB_list = ['6WJR']
PDB_list = ['6TPQ']
PDB_list = ['4KTG']
version = "_3.160_2.5"
PDB_list = ['5KCR']
PDB_list = ['7ECF']  # DNA quadruplex
PDB_list = ['4TNA']
PDB_list = ['5J7L']
PDB_list = ['4V9F','5J7L','4ARC']
PDB_list = ['4ARC']
PDB_list = ['4ARC']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.201/all/csv']
PDB_list = ['4V9F','6ZMI','7K00']

PDB_list = ['4V9F']
#PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.216/3.0A/csv']
#PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.216/2.0A/csv']

base_seq_list = ['DA','DT','DC','DG']  # for DNA
base_seq_list = []                     # for all nucleic acids, modified or not
base_seq_list = ['A','U','C','G']      # for RNA

nt_reference_point = "base"
atom_atom_min_distance = 5    # minimum distance between atoms in nts to consider them interacting

# plot one instance of each of the pairwise interactions
PlotPair = False
PlotPair = True
AlreadyPlotted = {}

ShowStructureReadingErrors = False
ShowStructureReadingErrors = True

# this path should be specified in localpath.py
# intended for writing out a .pickle file to be used by the FR3D motif search tool
unit_data_path = "C:/Users/zirbel/Documents/GitHub/fr3d-python/data/units"

# annotate all nucleotides in all chains, even when a representative set is used
annotate_entire_PDB_files = True

if __name__=="__main__":

    timerData = myTimer("start")
    lastwritetime = time()

    allInteractionDictionary = defaultdict(list)

    timerData = myTimer("Making PDB list",timerData)

    PDB_IFE_Dict = map_PDB_list_to_PDB_IFE_dict(PDB_list)

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

        if annotate_entire_PDB_files:

            pair_file = "%s_pairs_%s.pickle" % (PDB,fr3d_classification_version)
            pair_to_data_output_file = outputNAPairwiseInteractions + pair_file

            if not os.path.exists(pair_to_data_output_file):

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

                # Hydrogens are already added when the structure is loaded
                #timerData = myTimer("Inferring hydrogens",timerData)
                #print("  Adding hydrogens")
                #structure.infer_hydrogens()  # add hydrogens to NA bases and amino acids; slow

                # write out data file of nucleotide centers and rotations that can be used by FR3D for searches
                # need to be able to identify each chain that is available
                # write_unit_data_file(PDB,unit_data_path,structure)

                interaction_to_triple_list, pair_to_interaction, pair_to_data, timerData = annotate_nt_nt_in_structure(structure,timerData)

                # turn this off during development and testing
                if False:
                    print("  Annotated these interactions: %s" % interaction_to_triple_list.keys())
                    pickle.dump(interaction_to_triple_list,open(outputDataFilePickle,"wb"),2)

                timerData = myTimer("Recording interactions",timerData)
                print('Writing data file %s' % pair_to_data_output_file)
                pickle.dump(pair_to_data,open(pair_to_data_output_file,"wb"),2)

                myTimer("summary",timerData)


        else:
            # only process individual IFEs
            # this has not been tested recently and may need to be modified
            # This would be most relevant for statistical tallies, finding exemplars, etc.
            # But that could be done by just downloading the annotations, not creating them anew
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
            interaction_to_triple_list, pair_to_interaction, pair_to_data, timerData = annotate_nt_nt_interactions(bases, nt_nt_screen_distance, baseCubeList, baseCubeNeighbors, timerData)

            # used to return Python_pairs, pair_to_data, timerData

            timerData = myTimer("Recording interactions",timerData)

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

            myTimer("summary",timerData)

    # when appropriate, write out HTML files
    """
    if len(PDB_IFE_Dict) > 100:
        print("Writing " + outputDataFile)
        timerData = myTimer("Writing HTML files",timerData)
        pickle.dump((allInteractionDictionary,allAATwoBaseDictionary,PDB_list),open(outputDataFile,"wb"))
        writeInteractionsHTML(allInteractionDictionary,outputHTML,version)
    """

    myTimer("summary",timerData)
