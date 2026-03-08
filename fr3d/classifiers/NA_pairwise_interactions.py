# -*- coding: utf-8 -*-
"""
This program reads one or more CIF/PDB files and produces annotations
of nucleotide-nucleotide interactions.
Basepairs are annotated with Leontis-Westhof annotations like cWW, tHS, ...
Many interactions are annotated as "near" basepairs like ncWW, ntHS; use "near" in the category list
In families like cWW some basepairs are annotated like AA cWw to show which base does what; basepair,lower
A few basepairs have "alternative" geometries like cWWa; use basepair,alternative in the category list
cWB is for Table 13 in Leontis-Stombaugh-Westhof, WC-bifurcated interactions; use basepair,cwb
To get near, lower, alternative, and cwb use category basepair_detail

Usage examples:
python NA_pairwise_interactions.py 4TNA
python NA_pairwise_interactions.py -c basepair,stacking 4TNA
python NA_pairwise_interactions.py -c basepair_detail,sugar_ribose 8B0X

When fr3d_python is changed, reinstall by changing directory to fr3d-python, then:
python -m pip install .
"""

import argparse
from collections import defaultdict
import gzip
import math
import numpy as np
import os
import pickle
import sys
from time import time
import traceback
import urllib

# import the version of urlretrieve appropriate to the Python version
if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
    from urllib import urlopen
    read_mode = 'rb'
    write_mode = 'w'
else:
    from urllib.request import urlretrieve as urlretrieve
    from urllib.request import urlopen
    from urllib import request           # not sure why
    read_mode = 'rt'
    write_mode = 'wt'   # write as text

from fr3d.definitions import NAbaseheavyatoms
from fr3d.definitions import NAbasehydrogens
from fr3d.definitions import NAbaseatoms
from fr3d.definitions import nt_sugar
from fr3d.definitions import nt_phosphate
from fr3d.definitions import aa_fg
from fr3d.definitions import aa_linker
from fr3d.definitions import aa_backbone
from fr3d.definitions import planar_atoms
from fr3d.definitions import NAbaseMassiveAndHydrogens

from fr3d.classifiers.class_limits_2024 import nt_nt_cutoffs   # use the latest cutoffs
from fr3d.classifiers.hydrogen_bonds import load_ideal_basepair_hydrogen_bonds
from fr3d.classifiers.hydrogen_bonds import check_hydrogen_bond

# Modified nucleotide mappings
from fr3d.modified.mapping import modified_base_atom_list,parent_atom_to_modified,modified_atom_to_parent,modified_base_to_parent

# read input and output paths from localpath.py
# note that fr3d.localpath does not synchronize with Git, so you can change it locally to point to your own directory structure
try:
    from fr3d.localpath import outputNAPairwiseInteractions
    from fr3d.localpath import inputPath
except:
    inputPath = ""
    outputNAPairwiseInteractions = ""

oo_distance_center_center_distance_cutoff = 20
base_backbone_center_center_distance_cutoff = 17
standard_center_center_distance_cutoff = 12

all_categories = "basepair,basepair_detail,coplanar,stacking,backbone,so,covalent,sugar_ribose,near,bss,loops,oo_distance"

near_discrepancy_cutoff = 1.0     # maximum discrepancy to report as a near pair
near_discrepancy_cutoff = 2.0     # maximum discrepancy to report as a near pair
near_heavy_distance_cutoff = 4.2  # maximum distance between heavy atoms to be considered a near pair
true_heavy_distance_cutoff = 3.8  # maximum distance between heavy atoms to be considered a true pair

HB_donor_hydrogens = {}
HB_donor_hydrogens['A'] = {"N6":["1H6","2H6"], "C2":["H2"], "C8":["H8"], "O2'":[]}
HB_donor_hydrogens['G'] = {"N1":["H1"], "N2":["2H2","1H2"], "C8":["H8"], "O2'":[]}
HB_donor_hydrogens['C'] = {"N4":["1H4","2H4"], "C5":["H5"], "C6":["H6"], "O2'":[]}
HB_donor_hydrogens['U'] = {"N3":["H3"], "C5":["H5"], "C6":["H6"], "O2'":[]}

standard_bases = ['A','C','G','U','DA','DC','DG','DT']

Leontis_Westhof_basepairs = ['cWW','tWW','cWH','cHW','tWH','tHW','cWS','cSW','tWS','tSW','cHH','tHH','cHS','cSH','tHS','tSH','cSS','tSS']
Leontis_Westhof_basepairs += ['cWWa','tWWa','cHHa','tHSa','tSHa','cWB','cBW']  # alternative cases

Leontis_Westhof_basepairs_lower = [x.lower() for x in Leontis_Westhof_basepairs]

nt_reference_point = "base"
atom_atom_min_distance = 5    # minimum distance between atoms in nts to consider them interacting
base_seq_list = []                     # for all nucleic acids, modified or not

verbose = 2  # also print diagnostic information about basepairs
verbose = 3  # also print information about loops
verbose = 0  # do not print much at all
verbose = 1  # print basic information about input, output, and number of interactions

def print_dictionary(datapoint):
    for key,value in sorted(datapoint.items()):
        if type(value) == dict:
            print(key,' is a dictionary:')
            for k,v in sorted(value.items()):
                print("  %s = %s" % (k,v))
        else:
            print(key,value)


def focus_basepair_cutoffs(basepair_cutoffs,interactions):
    """
    Reduce the dictionary of basepair cutoffs to just the pairs
    that need to be annotated in this run.
    """

    focused_basepair_cutoffs = {}
    desired_families_lower = set([])

    if interactions:
        # reduced set of interactions is passed in, like maybe only cWW
        for interaction in interactions:
            desired_families_lower.add(interaction.lower())
    else:
        # use all available interactions
        for combination in basepair_cutoffs.keys():
            for interaction in basepair_cutoffs[combination]:
                desired_families_lower.add(interaction.lower())

    for combination in basepair_cutoffs.keys():
        focused_basepair_cutoffs[combination] = {}
        focused_basepair_cutoffs[combination][1] = {}   # interactions with positive normal
        focused_basepair_cutoffs[combination][-1] = {}  # interactions with negative normal
        for interaction in basepair_cutoffs[combination].keys():
            # extract just the family, to compare to the input list
            family = interaction.lower().replace("a","").replace("b","").replace("n","")

            if family in desired_families_lower or interaction.lower() in desired_families_lower:
                subcat = list(basepair_cutoffs[combination][interaction].keys())[0]
                if basepair_cutoffs[combination][interaction][subcat]["normalmin"] > 0:
                    focused_basepair_cutoffs[combination][1][interaction] = {}   # interactions with positive normal
                    for subcategory in basepair_cutoffs[combination][interaction]:
                        focused_basepair_cutoffs[combination][1][interaction][subcategory] = basepair_cutoffs[combination][interaction][subcategory]
                else:
                    focused_basepair_cutoffs[combination][-1][interaction] = {}   # interactions with negative normal
                    for subcategory in basepair_cutoffs[combination][interaction]:
                        focused_basepair_cutoffs[combination][-1][interaction][subcategory] = basepair_cutoffs[combination][interaction][subcategory]

    # check how this worked
    # for combination in focused_basepair_cutoffs:
    #     for normal in focused_basepair_cutoffs[combination]:
    #         for interaction in focused_basepair_cutoffs[combination][normal]:
    #             for subcategory in focused_basepair_cutoffs[combination][normal][interaction]:
    #                 print(combination, normal, interaction, subcategory, focused_basepair_cutoffs[combination][normal][interaction][subcategory])
    #                 pass

    return focused_basepair_cutoffs

def myTimer(state,data={}):

    # add elapsed time to the current state of the timer
    if "currentState" in data:
        currentState = data["currentState"]
        data[currentState] += time() - data["lastTime"]

    if state == "summary":
        total = 0.000000000001
        for state in data["allStates"]:
            if not state == "lastTime" and not state == "currentState":
                total += data[state]

        print("Summary of time taken:")
        for state in data["allStates"]:
            if not state == "lastTime" and not state == "currentState" and not state == "start":
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


def load_structure(filename,file_id="",preferred_id=None):
    """
    filename is the full path to a .pdb or .cif file
    file_id could be a 4-character PDB identifier, but could be otherwise
    """

    if not file_id:
        path,file_id = os.path.split(filename)
        file_id = file_id.replace(".cif","").replace(".pdb","").replace(".gz","")

    message = []
    original_filename = filename

    # look for the file, possibly with extensions
    if os.path.exists(filename):
        pass
    elif os.path.exists(filename+".cif.gz"):
        filename = filename + ".cif.gz"
    elif os.path.exists(filename+".cif"):
        filename = filename + ".cif"
    elif os.path.exists(filename+".pdb.gz"):
        filename = filename + ".pdb.gz"
    elif os.path.exists(filename+".pdb"):
        filename = filename + ".pdb"

    if verbose >= 1:
        print("  NA_pairwise_interactions: filename is %s" % filename)

    # if still not available, try to download from PDB and save locally
    # download .gz version when possible for speed and to save disk space
    if not os.path.exists(filename):
        if filename.lower().endswith('.cif.gz'):
            download_id = file_id + '.cif.gz'
        elif filename.lower().endswith('.cif'):
            download_id = file_id + '.cif.gz'
            filename = filename + ".gz"
        elif filename.lower().endswith('.pdb.gz'):
            download_id = file_id + '.pdb'
            filename = filename.rstrip('.gz')  # remove .gz because *.pdb.gz is not availble from PDB
        elif filename.lower().endswith('.pdb'):
            download_id = file_id + '.pdb'
        else:
            download_id = file_id + '.cif.gz'
            filename = filename + '.cif.gz'

        url = "https://files.rcsb.org/download/%s" % download_id

        try:
            urlretrieve(url, filename)
        except:
            message.append("Not able to download %s from %s" % (original_filename,url))
            message.append("Tried filename %s" % filename)
            return None, message

        # TODO: detect when this downloads an error file instead; current code is clumsy
        try:
            with open(filename,read_mode) as f:
                lines = f.read()

            if "404 Not Found" in lines:
                message.append("Not able to download %s from %s" % (download_id,url))
                if os.path.exists(filename):
                    os.remove(filename)
                message.append("Code is not clever enough to find or download %s" % original_filename)
                return None, message
        except:
            message.append("Downloaded %s from %s" % (download_id,url))

    # read the file from the disk
    try:
        rm = read_mode
        if filename.lower().endswith('.cif.gz'):
            with gzip.open(filename, rm) as raw:
                from fr3d.cif.reader import Cif
                cif_access = Cif(raw,preferred_id=preferred_id)
                structure = cif_access.structure()
        elif filename.lower().endswith('.cif'):
            with open(filename, rm) as raw:
                from fr3d.cif.reader import Cif
                structure = Cif(raw,preferred_id=preferred_id).structure()
        elif filename.lower().endswith('.pdb.gz'):
            with gzip.open(filename, rm) as raw:
                from fr3d.pdb.pdb_reader import PDBStructure
                structure = PDBStructure(file_id,raw).structures()
                message.append("No symmetry operators applied to .pdb files")
        elif filename.lower().endswith('.pdb'):
            with open(filename, rm) as raw:
                from fr3d.pdb.pdb_reader import PDBStructure
                structure = PDBStructure(file_id,raw).structures()
                message.append("No symmetry operators applied to .pdb files")

        message.append("Loaded " + filename)
        return structure, message

    except TypeError:
        message.append("TypeError when loading %s, loading a different way" % filename)
        rm = 'r'      # needed on Ubuntu
        if filename.lower().endswith('.cif.gz'):
            with gzip.open(filename, rm) as raw:
                from fr3d.cif.reader import Cif
                structure = Cif(raw,preferred_id=preferred_id).structure()
        elif filename.lower().endswith('.cif'):
            with open(filename, rm) as raw:
                from fr3d.cif.reader import Cif
                structure = Cif(raw,preferred_id=preferred_id).structure()
        elif filename.lower().endswith('.pdb.gz'):
            with gzip.open(filename, rm) as raw:
                from fr3d.pdb.pdb_reader import PDBStructure
                structure = PDBStructure(file_id,raw).structures()
                if verbose >=1:
                    print("  No symmetry operators applied to .pdb files")
        elif filename.lower().endswith('.pdb'):
            with open(filename, rm) as raw:
                from fr3d.pdb.pdb_reader import PDBStructure
                structure = PDBStructure(file_id,raw).structures()
                if verbose >=1:
                    print("  No symmetry operators applied to .pdb files")

        message.append("Loaded " + filename)
        return structure, message

    except Exception as ex:
        message.append("  Could not load %s due to exception %s: %s" % (filename,type(ex).__name__,ex))
        if type(ex).__name__ == "TypeError":
            message.append("  See suggestions in the fr3d-python Readme file")

        traceback_details = traceback.format_exc()
        print("Complete Traceback:")
        print(traceback_details)

        return None, message

    message.append("Could not load %s" % (filename))
    return None, message


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

    return atom_to_part_list


def get_atom_coordinates(nt,atom_names):
    """
    Get coordinates of the specified atoms,
    mapping to a modified nucleotide if necessary.
    """

    coordinates = []
    seq = nt.sequence

    for atom_name in atom_names:
        # check if there is a mapping
        if seq in parent_atom_to_modified:
            # map the atom name
            if atom_name in parent_atom_to_modified[seq]:
                coordinate = nt.centers[parent_atom_to_modified[seq][atom_name]]
            else:
                # hope for the best, or get an empty vector
                coordinate = nt.centers[atom_name]
        else:
            # default
            coordinate = nt.centers[atom_name]

        coordinates.append(coordinate)

    return coordinates


def get_one_atom_coordinates(nt,atom_name):
    """
    Get coordinates of the specified atom,
    mapping to a modified nucleotide if necessary.
    """

    seq = nt.sequence

    # check if there is a mapping
    if seq in ['A','C','G','U','DA','DC','DG','DT']:
        # standard, fast
        coordinates = nt.centers[atom_name]
    elif seq in parent_atom_to_modified:
        # map the atom name
        if atom_name in parent_atom_to_modified[seq]:
            coordinates = nt.centers[parent_atom_to_modified[seq][atom_name]]
        else:
            # hope for the best, or get an empty vector
            coordinates = nt.centers[atom_name]
    else:
        # default
        coordinates = nt.centers[atom_name]

    return coordinates


def make_nt_cubes_full(bases, screen_distance_cutoff, nt_reference="base"):
    """
    Builds cubes with side length screen_distance_cutoff
    using nt_reference as the point for each nucleotide.
    Cubes are named by a rounded value of x,y,z and by model.
    All 26 neighboring cubes are generated.
    """

    # build a set of cubes and record which bases are in which cube
    # also record which other cubes are neighbors of each cube
    baseCubeList = {}
    baseCubeNeighbors = {}

    # build a set of cubes and record which bases are in which cube
    for base in bases:
        center = base.centers[nt_reference]  # chosen reference point
        if len(center) == 3:
            x = math.floor(center[0]/screen_distance_cutoff)
            y = math.floor(center[1]/screen_distance_cutoff)
            z = math.floor(center[2]/screen_distance_cutoff)
            model = base.model
            key = "%d,%d,%d,%s" % (x,y,z,model)
            if key in baseCubeList:
                baseCubeList[key].append(base)
            else:
                baseCubeList[key] = [base]
                baseCubeNeighbors[key] = []
                for a in [-1,0,1]:
                    for b in [-1,0,1]:
                        for c in [-1,0,1]:
                            k = "%d,%d,%d,%s" % (x+a,y+b,z+c,model)
                            baseCubeNeighbors[key].append(k)

    return baseCubeList, baseCubeNeighbors


def make_nt_cubes_half(bases, screen_distance_cutoff, nt_reference="base"):
    """
    Builds cubes with side length screen_distance_cutoff
    using nt_reference as the point for each nucleotide.
    Cubes are named by a rounded value of x,y,z and by model.
    Only 13 neighboring cubes are generated, so each pair
    of nucleotides will only be generated once.
    """

    # build a set of cubes and record which bases are in which cube
    # also record which other cubes are neighbors of each cube
    baseCubeList = {}
    baseCubeNeighbors = {}
    # build a set of cubes and record which bases are in which cube
    for base in bases:
        center = base.centers[nt_reference]  # chosen reference point
        if len(center) == 3:
            x = math.floor(center[0]/screen_distance_cutoff)
            y = math.floor(center[1]/screen_distance_cutoff)
            z = math.floor(center[2]/screen_distance_cutoff)
            model = base.model
            key = "%d,%d,%d,%s" % (x,y,z,model)
            if key in baseCubeList:
                baseCubeList[key].append(base)
            else:
                baseCubeList[key] = [base]
                baseCubeNeighbors[key] = []
                # same cube and 13 neighbors, no two in opposite directions
                cubes = [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [0, 1, -1], [1, 0, 0], [1, 0, 1], [1, 0, -1], [1, 1, 0], [1, 1, 1], [1, 1, -1], [1, -1, 0], [1, -1, 1], [1, -1, -1]]
                for a,b,c in cubes:
                    k = "%d,%d,%d,%s" % (x+a,y+b,z+c,model)
                    baseCubeNeighbors[key].append(k)

    return baseCubeList, baseCubeNeighbors


def reverse_edges(inter):

    if len(inter) <= 2:
        rev = inter
    elif inter == 'N/A':
        rev = inter
    elif len(inter) == 3:
        rev = inter[0] + inter[2] + inter[1]
    elif len(inter) == 4 and inter[0] == 'n':            # like ntSH
        rev = inter[0] + inter[1] + inter[3] + inter[2]
    elif len(inter) == 4 and inter[0] == '!':            # like !tSH
        rev = inter[0] + inter[1] + inter[3] + inter[2]
    elif len(inter) == 4:                                # like tSHa
        rev = inter[0] + inter[2] + inter[1] + inter[3]
    elif len(inter) == 5:
        rev = inter[0] + inter[1] + inter[3] + inter[2] + inter[4]
    elif len(inter.split(" ")) == 2:
        fields = inter.split(" ")
        if fields[0] == 'cur':
            rev = inter[0:5] + inter[6] + inter[5]
        else:
            rev = inter[0] + inter[2] + inter[1] + inter[3:]
    else:
        rev = inter[0:(len(inter)-2)] + inter[len(inter)-1] + inter[len(inter)-2]

    return rev

# def makeListOfNtIndices(baseCubeList, baseCubeNeighbors):
#     """
#     This function returns a sorted list of all the nts indices in ascending order.
#     It was added as a method to be able to extract information about the O3' atom of the previous nucleotide
#     """
#     lastNT = {}
#     for nt1key in baseCubeList:                         # key to first cube
#         for nt2key in baseCubeNeighbors[nt1key]:        # key to each potential neighboring cube, including the first
#             if nt2key in baseCubeList:                  # if this cube was actually made
#                 for nt1 in baseCubeList[nt1key]:
#                     if nt1.index not in lastNT:
#                         lastNT[nt1.index] = nt1
#     return lastNT


def map_unit_id_to_previous_O3(bases):
    """
    Create a dictionary whose key is unit id and whose value
    is the 3d coordinates of the O3' atom of the previous nucleotide,
    if available, otherwise empty vector.
    """

    list_of_nucleotides = []
    for base in bases:
        coordinates = get_one_atom_coordinates(base,"O3'")
        P = get_one_atom_coordinates(base,"P")
        t = (base.model,base.symmetry,base.chain,base.index or -99,base.unit_id(),coordinates,P)
        list_of_nucleotides.append(t)

    list_of_nucleotides.sort()

    previous_O3_coordinates = np.empty([1,3])
    previous_model = None
    previous_symmetry = None
    previous_chain = None
    previous_index = None
    unit_id_to_previous_O3 = {}

    for model,symmetry,chain,index,unit_id,O3_coordinates,P in list_of_nucleotides:

        if not index:
            continue

        if model == previous_model and symmetry == previous_symmetry and chain == previous_chain and index == previous_index + 1:
            unit_id_to_previous_O3[unit_id] = previous_O3_coordinates
        else:
            unit_id_to_previous_O3[unit_id] = np.empty([1,3])

        # save current values for next nucleotide
        previous_model = model
        previous_symmetry = symmetry
        previous_chain = chain
        previous_index = index
        previous_O3_coordinates = O3_coordinates

    return unit_id_to_previous_O3


def check_for_two_interactions_on_same_edge(unit_id_to_basepairs,get_datapoint=False):
    """
    Loop over nucleotides, find those with two or more interactions on the same edge,
    choose the best, remove the others
    """

    # set of tuples of unit ids to leave out of the basepair list
    remove_pairs = set()
    make_near_pairs = set()

    for unit_id, basepairs in unit_id_to_basepairs.items():
        if len(basepairs) > 1:
            # unit_id makes more than one basepair
            for i in range(len(basepairs)-1):
                interaction_1, quality_1, unit_id_1 = basepairs[i]

                if (unit_id,unit_id_1) in remove_pairs:
                    # pair is already set to be removed
                    continue

                e1 = interaction_1.replace("n","")[1].lower()   # base edge
                f1 = unit_id_1.split("|")

                for j in range(i+1,len(basepairs)):
                    interaction_2, quality_2, unit_id_2 = basepairs[j]

                    if (unit_id,unit_id_2) in remove_pairs:
                        # already dealt with this pair
                        continue

                    if (unit_id,unit_id_2) in make_near_pairs:
                        # already dealt with this pair
                        continue

                    e2 = interaction_2.replace("n","")[1].lower()    # base edge

                    if not e1 == e2:
                        # different edges
                        continue

                    if len(f1) >= 7:
                        # possible to have an alternate id
                        f2 = unit_id_2.split("|")
                        if len(f1) == len(f2):
                            # possible to be the same unit with different alternate ids
                            f1[6] = ''
                            f2[6] = ''
                            if f1 == f2:
                                # same unit with different alternate ids
                                continue

                    common_atoms = set(quality_1['atoms1']) & set(quality_2['atoms1'])

                    if len(common_atoms) > 0:

                        if get_datapoint:
                            unit_id_list = unit_id
                            if verbose >= 2:
                                print("  Base %s makes multiple basepairs listed %d and %d below" % (unit_id,i,j))

                                for bp in basepairs:
                                    print(bp)
                                    interaction, quality, u1 = bp
                                    unit_id_list += "," + u1

                                print("  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s" % unit_id_list)
                                print('  Common atoms %s' % common_atoms)

                        atom_names = "".join(common_atoms)
                        if "H" in atom_names or len(common_atoms) > 1:

                            if verbose >= 2:
                                if "H" in atom_names:
                                    print('  Common atoms %s include a hydrogen, checking for conflicts' % common_atoms)
                                else:
                                    print('  Two or more common atoms %s, checking for conflicts' % common_atoms)

                            # conflicting basepairs
                            if interaction_1.startswith("n") and interaction_2.startswith("n"):
                                # both near, remove the worse one if it's pretty bad
                                if quality_1['cutoff_distance'] < quality_2['cutoff_distance']:
                                    if quality_2['cutoff_distance'] > 0.5 * near_discrepancy_cutoff:
                                        remove_pairs.add((unit_id,unit_id_2))
                                        remove_pairs.add((unit_id_2,unit_id))
                                        if verbose >= 2:
                                            print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  removed due to conflicting edge, both near" % (interaction_2,unit_id,unit_id_2,unit_id,unit_id_2))
                                            print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  not removed" % (interaction_1,unit_id,unit_id_1,unit_id,unit_id_1))
                                else:
                                    if quality_1['cutoff_distance'] > 0.5 * near_discrepancy_cutoff:
                                        remove_pairs.add((unit_id,unit_id_1))
                                        remove_pairs.add((unit_id_1,unit_id))
                                        if verbose >= 2:
                                            print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  removed due to conflicting edge, both near" % (interaction_1,unit_id,unit_id_1,unit_id,unit_id_1))
                                            print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  not removed" % (interaction_2,unit_id,unit_id_2,unit_id,unit_id_2))

                            elif not interaction_1.startswith("n") and not interaction_2.startswith("n"):
                                # both true, make one near
                                if quality_1['max_gap'] < quality_2['max_gap']:
                                    make_near_pairs.add((unit_id,unit_id_2))
                                    make_near_pairs.add((unit_id_2,unit_id))
                                    if verbose >= 2:
                                        print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  switched to near due to conflicting edge" % (interaction_2,unit_id,unit_id_2,unit_id,unit_id_2))
                                        print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  not switched" % (interaction_1,unit_id,unit_id_1,unit_id,unit_id_1))
                                else:
                                    make_near_pairs.add((unit_id,unit_id_1))
                                    make_near_pairs.add((unit_id_1,unit_id))
                                    if verbose >= 2:
                                        print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  switched to near due to conflicting edge" % (interaction_1,unit_id,unit_id_1,unit_id,unit_id_1))
                                        print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  not switched" % (interaction_2,unit_id,unit_id_2,unit_id,unit_id_2))

                            else:
                                # one near, one true, remove the near one if it's bad
                                if quality_2['cutoff_distance'] > 0.5 * near_discrepancy_cutoff:
                                    remove_pairs.add((unit_id,unit_id_2))
                                    remove_pairs.add((unit_id_2,unit_id))
                                    if verbose >= 2:
                                        print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  removed due to conflicting edge, cutoffs" % (interaction_2,unit_id,unit_id_2,unit_id,unit_id_2))
                                        print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  not removed" % (interaction_1,unit_id,unit_id_1,unit_id,unit_id_1))
                                elif quality_1['cutoff_distance'] > 0.5 * near_discrepancy_cutoff:
                                    remove_pairs.add((unit_id,unit_id_1))
                                    remove_pairs.add((unit_id_1,unit_id))
                                    if verbose >= 2:
                                        print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  removed due to conflicting edge, cutoffs" % (interaction_1,unit_id,unit_id_1,unit_id,unit_id_1))
                                        print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  not removed" % (interaction_2,unit_id,unit_id_2,unit_id,unit_id_2))

    return remove_pairs, make_near_pairs


def annotate_nt_nt_interactions(bases, center_center_distance_cutoff, baseCubeList, baseCubeNeighbors, categories, focused_basepair_cutoffs, ideal_hydrogen_bonds, timerData, get_datapoint = False):
    """
    loop through nt cubes, loop through neighboring nt cubes,
    then loop through bases in the two cubes,
    screening distances between them, then annotating interactions
    When get_datapoint is True, collect data about each pair to pass back
    """

    count_pair = 0

    interaction_to_pair_list = defaultdict(list) # map interaction to list of pairs
    category_to_interactions = defaultdict(set)  # map category to list of observed interactions

    pair_to_data = defaultdict(dict)             # place to record data for diagnostic purposes

    unit_id_to_basepairs = defaultdict(list) # map unit_id and edge to list of basepairs with their quality

    max_center_center_distance = 0     # record the largest screening distance for which an interaction is found
    bph_center_center_distance = 0     # record the largest screening distance for which an interaction is found
    oo_center_center_distance = 0     # record the largest screening distance for which an interaction is found

    basepair_parent_base_combination_set = set(['A,A','A,C','A,G','A,U','C,C','G,C','C,U','G,G','G,U','U,U','A,DT','C,DT','G,DT','DT,DT'])

    # keep track of overlapping chains
    overlapping_chains = set()

    # For base-backbone interactions, we need to know
    if 'backbone' in categories.keys():
        #ntDict = makeListOfNtIndices(baseCubeList, baseCubeNeighbors)
        unit_id_to_previous_O3 = map_unit_id_to_previous_O3(bases)

    for nt1key in baseCubeList:                         # key to first cube
        for nt2key in baseCubeNeighbors[nt1key]:        # key to each potential neighboring cube, including the first
            if nt2key in baseCubeList:                  # if this cube was actually made
                for nt1 in baseCubeList[nt1key]:        # first nt of a potential pair

                    if len(nt1.centers["base"]) < 3:
                        if verbose >= 2:
                            print("  Missing base center for %s" % nt1.unit_id())
                            print(nt1.centers["base"])
                        continue

                    parent1 = get_parent(nt1.sequence)   # map modified nts to parent nt A, C, G, U, DT
                    if not parent1:
                        if verbose >= 2:
                            print("  No parent for %s" % nt1.unit_id())
                        continue

                    gly1 = get_glycosidic_atom_coordinates(nt1,parent1)

                    if len(gly1) < 3:
                        if verbose >= 2:
                            print("  Missing glycosidic atom for %s" % nt1.unit_id())
                        continue

                    number1 = nt1.number                 # nucleotide number

                    for nt2 in baseCubeList[nt2key]:           # second nt of a potential pair
                        if nt1.unit_id() == nt2.unit_id():
                            continue

                        # only consider each nt1, nt2 pair in one direction
                        # Those in different cubes only occur once
                        # Those from the same cube need a way to select just one pair
                        if nt1key == nt2key:
                            if nt1.chain > nt2.chain:
                                continue
                            elif nt1.chain == nt2.chain:
                                if nt1.index and nt2.index:
                                    if nt1.index > nt2.index:
                                        continue
                                elif nt1.unit_id() > nt2.unit_id():
                                    continue

                        if len(nt2.centers["base"]) < 3:
                            if verbose >= 2:
                                print("  Missing base center for %s" % nt2.unit_id())
                                print(nt2.centers["base"])
                            continue

                        # avoid some strange errors due to overlapping nucleotides
                        # not a complete solution, just a first step
                        if nt1.pdb == '1BVO':
                            # D pairs with symmetry operated E
                            # D and E are on top of each other
                            if nt1.symmetry == nt2.symmetry and nt1.chain != nt2.chain:
                                continue
                        elif nt1.pdb == '1R71' and nt1.chain != nt2.chain:
                            ok_chains = ['EF','FE','GH','HG','IJ','JI','KL','LK']
                            if not nt1.chain+nt2.chain in ok_chains:
                                continue
                        elif nt1.pdb == '3CRX' and nt1.chain != nt2.chain:
                            # D ASM1 is on top of F ASM2
                            ok_chains = ['DE','ED','CF','FC']
                            if not nt1.chain+nt2.chain in ok_chains:
                                continue
                        elif nt1.pdb == '3CZ3':
                            # altid B for chain E can be on top of altid A for chain F
                            if not nt1.alt_id == nt2.alt_id:
                                continue
                        elif nt1.pdb == '4BUL' and nt1.chain != nt2.chain:
                            ok_chains = ['EF','FE','GH','HG']
                            if not nt1.chain+nt2.chain in ok_chains:
                                continue
                        elif nt1.pdb == '4KTG':
                            # both symmetry operators put nucleotides in the same locations
                            if not nt1.symmetry == nt2.symmetry:
                                continue
                        elif nt1.pdb == '4WLS':
                            if nt1.chain in ['U','V'] and nt2.chain in ['X','Y']:
                                continue
                            elif nt1.chain in ['X','Y'] and nt2.chain in ['U','V']:
                                continue
                        elif nt1.pdb == '5A39':
                            ok_chains = ['CG','GC','DH','HD','EF','FE']
                            if not nt1.chain+nt2.chain in ok_chains:
                                continue
                        elif nt1.pdb == '5UA1':
                            ok_chains = ['CD','DC','EF','FE']
                            if not nt1.chain+nt2.chain in ok_chains:
                                continue
                        elif nt1.pdb == '5UA2':
                            if nt1.chain == nt2.chain and nt1.symmetry != nt2.symmetry:
                                continue
                        elif nt1.pdb == '6KHY':
                            if nt1.chain == 'G' and nt2.chain == 'H':
                                continue
                            elif nt1.chain == 'H' and nt2.chain == 'G':
                                continue
                            elif nt1.chain == 'I' and nt2.chain == 'J':
                                continue
                            elif nt1.chain == 'J' and nt2.chain == 'I':
                                continue

                        # vector displacement between base centers
                        displacement = abs(nt2.centers["base"]-nt1.centers["base"]) # center-center

                        # quick screens for base centers being too far apart
                        if displacement[0] > center_center_distance_cutoff or \
                           displacement[1] > center_center_distance_cutoff or \
                           displacement[2] > center_center_distance_cutoff:
                            continue

                        # avoid comparing alternate coordinates of the same nucleotide
                        # check in the order most likely to terminate the fastest
                        # do not check sequence, because sometimes ||A and ||B forms have different sequence
                        if number1 == nt2.number:
                            if nt1.chain == nt2.chain:
                                if nt1.symmetry == nt2.symmetry:
                                    if nt1.insertion_code == nt2.insertion_code:
                                        if nt1.alt_id != nt2.alt_id:
                                            #print("Skipping pair of alternate coordinates", (nt1.unit_id(),nt2.unit_id()))
                                            continue

                        # calculate actual center-center distance, screen
                        center_center_distance = np.linalg.norm(displacement)

                        # base centers are too far apart to interact, screen them out
                        if center_center_distance > center_center_distance_cutoff:
                            continue

                        # too short center_center_distance means overlapping nucleotides
                        if center_center_distance < 1:
                            print("Overlapping nucleotides",nt1.unit_id(),nt2.unit_id())
                            overlapping_chains.add((nt1.symmetry,nt1.chain,nt2.symmetry,nt2.chain))
                            overlapping_chains.add((nt2.symmetry,nt2.chain,nt1.symmetry,nt1.chain))

                        # some structures have overlapping nucleotides, screen those out
                        if center_center_distance < 2:
                            continue

                        unit_id_pair = (nt1.unit_id(),nt2.unit_id())  # tuple for these nucleotides in this order
                        reversed_pair = (nt2.unit_id(),nt1.unit_id())

                        parent2 = get_parent(nt2.sequence)
                        if not parent2:
                            if verbose >= 2:
                                print("  No parent for %s" % nt2.unit_id())
                            continue

                        parent_pair = parent1 + "," + parent2
                        parent_pair_reversed = parent2 + "," + parent1

                        marked_coplanar = False

                        # store data for diagnostics, if requested
                        if get_datapoint:
                            datapoint12 = {}
                            datapoint12['center_center_distance'] = center_center_distance
                            datapoint12['nt1_seq'] = nt1.sequence
                            datapoint12['nt2_seq'] = nt2.sequence
                            datapoint12['nt1_parent'] = parent1
                            datapoint12['nt2_parent'] = parent2
                            datapoint12['url'] = "https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt1.unit_id(),nt2.unit_id())

                            datapoint21 = {}
                            datapoint21['center_center_distance'] = center_center_distance
                            datapoint21['nt1_seq'] = nt2.sequence
                            datapoint21['nt2_seq'] = nt1.sequence
                            datapoint21['nt1_parent'] = parent2
                            datapoint21['nt2_parent'] = parent1
                            datapoint21['url'] = "https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt2.unit_id(),nt1.unit_id())

                        else:
                            datapoint12 = None
                            datapoint21 = None

                        # annotate backbone oxygen-oxygen distances
                        if 'oo_distance' in categories:
                            pair_list = check_oo_distance(nt1,nt2,parent1,parent2)
                            if pair_list:
                                count_pair += 1
                                interaction_to_pair_list['oo_distance'] += pair_list
                                category_to_interactions['oo_distance'].add('oo_distance')
                                oo_center_center_distance = max(oo_center_center_distance,center_center_distance)  # for setting optimally

                        # base centers are too far apart for the remaining interactions
                        if center_center_distance > base_backbone_center_center_distance_cutoff:
                            # store data for diagnostics, if requested
                            if datapoint12:
                                pair_to_data[unit_id_pair] = datapoint12

                            if datapoint21:
                                pair_to_data[reversed_pair] = datapoint21
                            continue

                        # annotate base phosphate and base ribose interactions
                        if 'backbone' in categories:
                            timerData = myTimer("Check backbone interactions", timerData)

                            # get coordinates of O3' of the nucleotide before nt2, part of the phosphate of nt2
                            previousO3 = unit_id_to_previous_O3.get(nt2.unit_id(),np.empty([1,3]))
                            interactionbPh, interactionbR, datapoint12 = check_base_backbone_interactions(nt1, nt2, previousO3, parent1, parent2, datapoint12)

                            if interactionbPh and len(interactionbPh) > 0:
                                count_pair += 1
                                interaction_to_pair_list[interactionbPh].append(unit_id_pair)
                                category_to_interactions['backbone'].add(interactionbPh)
                                bph_center_center_distance = max(bph_center_center_distance,center_center_distance)  # for setting optimally

                            if interactionbR and len(interactionbR) > 0:
                                count_pair += 1
                                interaction_to_pair_list[interactionbR].append(unit_id_pair)
                                category_to_interactions['backbone'].add(interactionbR)
                                bph_center_center_distance = max(bph_center_center_distance,center_center_distance)  # for setting optimally

                            # get coordinates of O3' of the nucleotide before nt1, part of the phosphate of nt1
                            previousO3 = unit_id_to_previous_O3.get(nt1.unit_id(),np.empty([1,3]))
                            interactionbPh, interactionbR, datapoint21 = check_base_backbone_interactions(nt2, nt1, previousO3, parent2, parent1, datapoint21)

                            if interactionbPh and len(interactionbPh) > 0:
                                count_pair += 1
                                interaction_to_pair_list[interactionbPh].append(reversed_pair)
                                category_to_interactions['backbone'].add(interactionbPh)
                                bph_center_center_distance = max(bph_center_center_distance,center_center_distance)  # for setting optimally
                            if interactionbR and len(interactionbR) > 0:
                                count_pair += 1
                                interaction_to_pair_list[interactionbR].append(reversed_pair)
                                category_to_interactions['backbone'].add(interactionbR)
                                bph_center_center_distance = max(bph_center_center_distance,center_center_distance)  # for setting optimally

                        # base centers are too far apart for the remaining interactions
                        if center_center_distance > standard_center_center_distance_cutoff:
                            # store data for diagnostics, if requested
                            if datapoint12:
                                pair_to_data[unit_id_pair] = datapoint12

                            if datapoint21:
                                pair_to_data[reversed_pair] = datapoint21
                            continue

                        # check base to oxygen stack; always base first, oxygen second
                        if 'so' in categories:
                            timerData = myTimer("Check base oxygen stack",timerData)
                            interaction, datapoint12, interaction_reversed = check_base_oxygen_stack_rings(nt1,nt2,parent1,datapoint12)

                            if len(interaction) > 0:
                                count_pair += 1
                                interaction_to_pair_list[interaction].append(unit_id_pair)
                                interaction_to_pair_list[interaction_reversed].append(reversed_pair)
                                max_center_center_distance = max(max_center_center_distance,center_center_distance)  # for setting optimally
                                category_to_interactions['so'].add(interaction)
                                category_to_interactions['so'].add(interaction_reversed)

                            interaction, datapoint21, interaction_reversed = check_base_oxygen_stack_rings(nt2,nt1,parent2,datapoint21)

                            if len(interaction) > 0:
                                count_pair += 1
                                interaction_to_pair_list[interaction].append(reversed_pair)
                                interaction_to_pair_list[interaction_reversed].append(unit_id_pair)
                                max_center_center_distance = max(max_center_center_distance,center_center_distance)  # for setting optimally
                                category_to_interactions['so'].add(interaction)
                                category_to_interactions['so'].add(interaction_reversed)

                        if 'stacking' in categories:
                            timerData = myTimer("Check base base stack", timerData)
                            interaction, datapoint12, interaction_reversed = check_base_base_stacking(nt1, nt2, parent1, parent2, datapoint12)
                            if len(interaction) > 0:
                                count_pair += 1
                                interaction_to_pair_list[interaction].append(unit_id_pair)
                                max_center_center_distance = max(max_center_center_distance,center_center_distance)  # for setting optimally
                                category_to_interactions['stacking'].add(interaction)
                                category_to_interactions['stacking'].add(interaction_reversed)

                        # annotate sugar ribose interactions
                        if 'sugar_ribose' in categories:
                            timerData = myTimer("Check sugar ribose", timerData)
                            if not parent1 in ['DA','DC','DG','DT'] and not parent2 in ['DA','DC','DG','DT']:
                                interaction, datapoint12 = check_sugar_ribose(nt1, nt2, parent1, datapoint12)
                                if len(interaction) > 0:
                                    count_pair += 1
                                    interaction_to_pair_list[interaction].append(unit_id_pair)
                                    category_to_interactions['sugar_ribose'].add(interaction)
                                    interaction_reversed = reverse_edges(interaction)
                                    category_to_interactions['sugar_ribose'].add(interaction_reversed)
                                    max_center_center_distance = max(max_center_center_distance,center_center_distance)  # for setting optimally

                                interaction, datapoint21 = check_sugar_ribose(nt2, nt1, parent2, datapoint21)
                                if len(interaction) > 0:
                                    count_pair += 1
                                    interaction_to_pair_list[interaction].append(reversed_pair)
                                    category_to_interactions['sugar_ribose'].add(interaction)
                                    interaction_reversed = reverse_edges(interaction)
                                    category_to_interactions['sugar_ribose'].add(interaction_reversed)
                                    max_center_center_distance = max(max_center_center_distance,center_center_distance)  # for setting optimally

                        # annotate basepairs
                        if 'basepair' in categories:
                            gly2 = get_glycosidic_atom_coordinates(nt2,parent2)
                            if len(gly2) < 3:
                                if verbose >= 2:
                                    print("  Missing glycosidic atom for %s" % nt2.unit_id())
                                continue

                            # always annotate cWW basepairs to be able to calculate crossing numbers
                            # check coplanar and basepairing for bases in specific orders
                            # AA, CC, GG, UU will be checked in both nucleotide orders, that's important
                            if parent_pair in basepair_parent_base_combination_set:

                                pair_data = {}
                                pair_data["glycosidic_displacement"] = np.subtract(gly2,gly1)
                                # vector from origin to nt2 when standardized
                                pair_data["displ12"] = np.dot(pair_data["glycosidic_displacement"],nt1.rotation_matrix)
                                pair_data["parent1"] = parent1
                                pair_data["parent2"] = parent2

                                if 'coplanar' in categories:
                                    timerData = myTimer("Check coplanar",timerData)
                                    pair_data, datapoint12 = check_coplanar(nt1,nt2,pair_data,datapoint12)

                                    # annotate coplanar relationship when present
                                    if pair_data['coplanar']:
                                        count_pair += 1
                                        interaction_to_pair_list['cp'].append(unit_id_pair)
                                        category_to_interactions['coplanar'].add('cp')
                                        marked_coplanar = True

                                timerData = myTimer("Check basepairing",timerData)

                                cutoffs = focused_basepair_cutoffs[parent1+","+parent2]
                                hydrogen_bonds = ideal_hydrogen_bonds[parent1+","+parent2]
                                interaction12, subcategory12, quality12, datapoint12 = check_basepair_cutoffs(nt1,nt2,pair_data,cutoffs,hydrogen_bonds,datapoint12)

                                interaction12_reversed = reverse_edges(interaction12)

                            else:
                                interaction12 = ""
                                interaction12_reversed = ""

                            # check pair in the other order
                            if parent_pair_reversed in basepair_parent_base_combination_set:

                                pair_data = {}
                                pair_data["glycosidic_displacement"] = np.subtract(gly1,gly2)
                                # vector from origin to nt2 when standardized
                                pair_data["displ12"] = np.dot(pair_data["glycosidic_displacement"],nt2.rotation_matrix)
                                pair_data["parent1"] = parent2
                                pair_data["parent2"] = parent1

                                if 'coplanar' in categories:
                                    timerData = myTimer("Check coplanar",timerData)
                                    pair_data, datapoint21 = check_coplanar(nt2,nt1,pair_data,datapoint21)

                                    # annotate coplanar relationship
                                    if pair_data['coplanar'] and not marked_coplanar:
                                        count_pair += 1
                                        interaction_to_pair_list['cp'].append(unit_id_pair)
                                        category_to_interactions['coplanar'].add('cp')

                                timerData = myTimer("Check basepairing",timerData)
                                cutoffs = focused_basepair_cutoffs[parent2+","+parent1]
                                hydrogen_bonds = ideal_hydrogen_bonds[parent2+","+parent1]
                                interaction21, subcategory21, quality21, datapoint21 = check_basepair_cutoffs(nt2,nt1,pair_data,cutoffs,hydrogen_bonds,datapoint21)

                                interaction21_reversed = reverse_edges(interaction21)

                            else:
                                interaction21 = ""
                                interaction21_reversed = ""

                            # if annotated interaction in both pair orders, choose the better one
                            if len(interaction12) > 0 and len(interaction21) > 0:
                                conflict_message = "%5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  duplicate interaction %-5s in second direction" % (interaction21,nt1.unit_id(),nt2.unit_id(),nt1.unit_id(),nt2.unit_id(),interaction12_reversed)

                                if interaction12_reversed.lower() == interaction21.lower():
                                    # same annotation, just different in order of edges
                                    interaction21 = ""      # ignore this one
                                    conflict_message = ""
                                elif "n" in interaction12 and "n" in interaction21:
                                    if quality12['cutoff_distance'] < quality21['cutoff_distance']:
                                        interaction21 = ""      # knock this one out
                                    else:
                                        interaction12 = ""      # knock this one out
                                elif "n" in interaction21:
                                    interaction21 = ""          # use true instead of near
                                    conflict_message = ""
                                elif "n" in interaction12:
                                    interaction12 = ""          # use true instead of near
                                    conflict_message = ""
                                else:
                                    # both true, but different
                                    conflict_message += "No clear way to decide between them"
                                    interaction21 = ""      # break the tie

                                if len(conflict_message) > 0:
                                    if len(interaction21) > 0:
                                        conflict_message += "  Using %s" % interaction21
                                    else:
                                        conflict_message += "  Using %s" % interaction12_reversed

                                    if verbose >= 2:
                                        print(conflict_message)

                                    # record conflicting interactions if desired
                                    if False and get_datapoint:
                                        with open(os.path.join(outputNAPairwiseInteractions,'conflicting.txt'),'a') as conf:
                                            conf.write(conflict_message+"\n")

                            if len(interaction12) > 0:
                                new_interaction = [interaction12,interaction12_reversed,subcategory12,quality12,nt1.unit_id(),nt2.unit_id()]

                                if datapoint12 and datapoint21:
                                    datapoint21['basepair'] = interaction12_reversed
                                    datapoint21['basepair_subcategory'] = datapoint12['basepair_subcategory']

                            elif len(interaction21) > 0:
                                interaction21_reversed = reverse_edges(interaction21)
                                new_interaction = [interaction21,interaction21_reversed,subcategory21,quality21,nt2.unit_id(),nt1.unit_id()]

                                if datapoint12 and datapoint21:
                                    datapoint12['basepair'] = interaction21_reversed
                                    datapoint12['basepair_subcategory'] = datapoint21['basepair_subcategory']

                            else:
                                new_interaction = []

                            if len(new_interaction) > 0:
                                count_pair += 1
                                max_center_center_distance = max(max_center_center_distance,center_center_distance)

                                # record the basepair interaction in both directions, according to edge
                                interaction, interaction_reversed, subcategory, quality, u1, u2 = new_interaction

                                # remove n and a from interaction, if present
                                interaction_clean = interaction.replace("n","").replace("a","")
                                # interaction_clean_reversed = reverse_edges(interaction_clean)

                                if not u1 in unit_id_to_basepairs:
                                    unit_id_to_basepairs[u1] = []
                                unit_id_to_basepairs[u1].append([interaction,quality,u2])

                                if not u2 in unit_id_to_basepairs:
                                    unit_id_to_basepairs[u2] = []

                                quality_reversed = {}
                                quality_reversed['cutoff_distance'] = quality['cutoff_distance']
                                quality_reversed['max_gap'] = quality['max_gap']
                                quality_reversed['atoms1'] = quality['atoms2']
                                quality_reversed['atoms2'] = quality['atoms1']
                                unit_id_to_basepairs[u2].append([interaction_reversed,quality_reversed,u1])

                        # store data for diagnostics, if requested
                        if datapoint12:
                            pair_to_data[unit_id_pair] = datapoint12

                        if datapoint21:
                            pair_to_data[reversed_pair] = datapoint21

    # check for two basepair interactions on the same edge
    if len(overlapping_chains) == 0:
        remove_pairs, make_near_pairs = check_for_two_interactions_on_same_edge(unit_id_to_basepairs,get_datapoint)
        if get_datapoint:
            # update annotation in datapoint
            for (n1,n2) in remove_pairs:
                pair_to_data[(n1,n2)]["basepair"] = ""
                pair_to_data[(n2,n1)]["basepair"] = ""
            for (n1,n2) in make_near_pairs:
                pair_to_data[(n1,n2)]["basepair"] = "n" + pair_to_data[(n1,n2)]["basepair"]
                pair_to_data[(n2,n1)]["basepair"] = "n" + pair_to_data[(n2,n1)]["basepair"]

    else:
        # when there are overlapping chains, some bases may make two cWW pairs, for example
        remove_pairs = []
        make_near_pairs = []

    # record remaining basepairs, but each one only once
    already_saved = set()
    for unit_id, basepairs in unit_id_to_basepairs.items():
        for interaction, quality, unit_id_2 in basepairs:
            if not (unit_id,unit_id_2) in remove_pairs and not (unit_id,unit_id_2) in already_saved:

                if (unit_id,unit_id_2) in make_near_pairs:
                    interaction = "n" + interaction

                interaction_to_pair_list[interaction].append((unit_id,unit_id_2))

                already_saved.add((unit_id_2,unit_id))
                already_saved.add((unit_id,unit_id_2))

                if not interaction in category_to_interactions['basepair']:
                    interaction_reversed = reverse_edges(interaction)
                    category_to_interactions['basepair'].add(interaction)
                    category_to_interactions['basepair'].add(interaction_reversed)
                    category_to_interactions['basepair_detail'].add(interaction)
                    category_to_interactions['basepair_detail'].add(interaction_reversed)

    file_id = ""
    for nt in bases:
        file_id = nt.pdb
        break

    if verbose >= 1:
        print("  Found %d nucleotide-nucleotide interactions in %s" % (count_pair,file_id))

    if verbose >= 3:
        print("  Maximum screen distance for standard contacts is %8.4f" % max_center_center_distance)
        print("  Maximum screen distance for bph      contacts is %8.4f" % bph_center_center_distance)
        print("  Maximum screen distance for oo       contacts is %8.4f" % oo_center_center_distance)

    if 'basepair' in categories:
        # calculate and save crossing numbers for each annoated interaction
        timerData = myTimer("Calculate crossing",timerData)
        if len(overlapping_chains) == 0:
            interaction_to_list_of_tuples = crossing_bss_loops(bases,interaction_to_pair_list,categories)
        else:
            interaction_to_list_of_tuples = crossing_bss_loops(bases,interaction_to_pair_list,{})

            with open('pdb_with_overlapping_chains','a') as f:
                for s1,c1,s2,c2 in overlapping_chains:
                    f.write("%s\t%s\t%s\t%s\t%s\n" % (file_id,s1,c1,s2,c2))

        if 'bSS' in interaction_to_list_of_tuples:
            category_to_interactions['bss'] = set(['bSS'])
    else:
        # this is probably not going to end well
        interaction_to_list_of_tuples = interaction_to_pair_list

    return interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data


def get_alt_id(unit_id):
    fields = unit_id.split("|")
    if len(fields) >= 7:
        alt_id = fields[6]
    else:
        alt_id = ""
    return alt_id


def select_partner(unitids,already_used):
    """
    If unitids has more than one member, avoid the ones in already_used
    and select the one with blank alt id, then alt id A, then alt id B, etc.
    """

    if len(unitids) == 1:
        return list(unitids)[0]
    elif len(unitids) == 0:
        return None

    remaining = unitids - already_used

    if len(remaining) == 0:
        return None

    unit_id_with_alt_id = []
    for unit_id in remaining:
        alt_id = get_alt_id(unit_id)
        unit_id_with_alt_id.append((unit_id,alt_id))

    unit_id_with_alt_id.sort(key=lambda x: x[1])

    return unit_id_with_alt_id[0][0]


def crossing_bss_loops(bases,interaction_to_pair_list,categories):
    """
    Identify which cWW pairs are nested.
    Then for each interaction, calculate the number of nested cWW pairs it crosses
    Also identify "border single stranded" or "bSS" pairs of nucleotides,
    which start and end a hairpin, or start and end each strand of an internal loop or junction loop.
    If requested, extract hairpin, internal, and multi-helix junction loops.

    The correct unit of study here is a chain, which is specified by a model number,
    chain identifier, and symmetry operator.  Called "MCS".  Should not need to be split often.
    """

    # map unit_id to model, chain, sequence index, base, symmetry
    # index is the hard one to get, the others are just for convenience
    # map model, chain, index to unit id
    # map model, chain to min and max index
    unit_id_to_fields = {}
    MCS_index_to_unit_id = {}  # problematic, with alternate ids having different unit ids
    MCS_to_min_index = defaultdict(lambda: 9999999)
    MCS_to_max_index = defaultdict(lambda: 0)
    for nt in bases:
        unit_id = nt.unit_id()
        fields = unit_id.split("|")

        # file_id = fields[0]
        model = fields[1]
        chain = fields[2]
        base  = fields[3]
        if len(fields) == 9:
            symmetry = fields[8]
        else:
            symmetry = "1_555"

        unit_id_to_fields[unit_id] = (model,chain,nt.index,base,symmetry)

        MCS = (model,chain,symmetry)

        if nt.index:
            MCS_to_min_index[MCS] = min(MCS_to_min_index[MCS],nt.index)
            MCS_to_max_index[MCS] = max(MCS_to_max_index[MCS],nt.index)

            if not MCS in MCS_index_to_unit_id:
                MCS_index_to_unit_id[MCS] = {}

        # in case of alternate ids (like A, B), do not overwrite the first one
        if nt.index and not nt.index in MCS_index_to_unit_id[MCS]:
            MCS_index_to_unit_id[MCS][nt.index] = unit_id

    # make variables to store the chain indices where nested cWW's occur
    MCS_to_nested_cWW_endpoints = {}
    MCS_to_endpoints = defaultdict(set)
    for MCS, index_to_unit_id in sorted(MCS_index_to_unit_id.items()):
        if verbose >= 3:
            print("  MCS %s has %4d nucleotides" % (MCS,len(index_to_unit_id)))

        # at first, each index maps to itself
        nested_cWW_endpoints = []
        for i in range(0,MCS_to_max_index[MCS]+1):
            nested_cWW_endpoints.append(i)
        MCS_to_nested_cWW_endpoints[MCS] = nested_cWW_endpoints

        # set up empty sets of endpoints for use later
        MCS_to_endpoints[MCS] = set()

    # find AU, GC, GU cWW basepairs within each model,chain,symmetry
    MCS_to_canonical_cWW_indices = defaultdict(list)   # separate list for each model and chain
    # two_chain_pairs = []                                 # canonical pairs between chains
    for interaction in interaction_to_pair_list.keys():
        if not interaction.lower() in ['cww','cwwa']:
            continue
        for u1,u2 in interaction_to_pair_list[interaction]:
            model1, chain1, index1, base1, symmetry1 = unit_id_to_fields[u1]
            model2, chain2, index2, base2, symmetry2 = unit_id_to_fields[u2]

            if not model1 == model2:
                # should never happen, but check, for good form
                continue

            # solitary nucleotides do not have an index, do not need this treatment
            if not index1:
                continue
            if not index2:
                continue

            # bases could be PSU or other modified base making WC pair
            parent1 = get_parent_as_RNA(base1,'')
            parent2 = get_parent_as_RNA(base2,'')

            # record AU, GC, GU cWW pairs by index within each chain and symmetry
            if parent1+parent2 in ['AU','UA','CG','GC','GU','UG']:
                if chain1 == chain2 and symmetry1 == symmetry2:
                    # record in increasing order
                    MCS = (model1,chain1,symmetry1)
                    if index1 < index2:
                        MCS_to_canonical_cWW_indices[MCS].append((index1,index2))
                    else:
                        MCS_to_canonical_cWW_indices[MCS].append((index2,index1))
                # else:
                #     # store canonical pairs between chains for processing later
                #     two_chain_pairs.append((model1,chain1,chain2,index1,index2))

    # within each model,chain,symmetry, sort nested cWW by distance between them
    # starting with the shortest-range pairs, record nested cWW pairs
    # by mapping one index to the other in MCS_to_nested_cWW_endpoints
    bss_endpoints = set()
    for MCS, pairs in MCS_to_canonical_cWW_indices.items():
        if verbose >=3:
            print("  Getting nested for model %s chain %s symmetry %s" % MCS)
        cWW_pairs = sorted(pairs, key=lambda p: (p[1]-p[0],p[0]))

        # loop over cWW pairs and if no conflict, record as being nested
        for index1,index2 in cWW_pairs:
            # loop over indices within this pair, see if they map outside this pair
            i = index1+1
            while i < index2 and MCS_to_nested_cWW_endpoints[MCS][i] > index1 and MCS_to_nested_cWW_endpoints[MCS][i] < index2:
                i += 1

            if i == index2:
                # this pair is nested, so record the endpoints
                MCS_to_nested_cWW_endpoints[MCS][index1] = index2
                MCS_to_nested_cWW_endpoints[MCS][index2] = index1

                if 'bss' in categories:
                    # these are the cWW endpoints that may be bss start and stop points
                    # so they can be in different model, chain, symmetry contexts
                    bss_endpoints.add((MCS,index1,MCS,index2))

            else:
                pass

    # loop over all pairs, calculate crossing number, the number of nested pairs crossed
    # record interacting pairs and their crossing number as triples

    # it's really important to do the cww pairs first, otherwise they get skipped and we lose bss
    # move cWW, cWw, cwW, cWWa to the beginning of the list
    interactions = interaction_to_pair_list.keys()
    interactions = sorted(interactions,key=lambda i: i.lower() not in ['cww','cwwa'])

    interaction_to_list_of_tuples = defaultdict(list)
    pairs_to_crossing = {}
    for interaction in interactions:

        if interaction == "":
            continue

        for itpl in interaction_to_pair_list[interaction]:
            u1 = itpl[0]
            u2 = itpl[1]
            if (u1,u2) in pairs_to_crossing:
                # no need to re-compute
                crossing = pairs_to_crossing[(u1,u2)]
            else:
                model1,chain1,index1,base1,symmetry1 = unit_id_to_fields[u1]
                model2,chain2,index2,base2,symmetry2 = unit_id_to_fields[u2]

                if not index1 or not index2:
                    # at least one of these nucleotides is not in a chain, so the
                    # crossing number is 0
                    crossing = 0
                else:

                    MCS1 = (model1,chain1,symmetry1)
                    MCS2 = (model2,chain2,symmetry2)

                    if not model1 == model2:
                        # should never happen, but check for good form
                        continue

                    crossing = 0

                    # interactions within the same chain can have non-zero crossing number
                    # some chains may not have any cWW pairs, then all interactions are nested
                    if chain1 == chain2 and symmetry1 == symmetry2:
                        if MCS1 in MCS_to_nested_cWW_endpoints:
                            # put indices of the current interaction in increasing order
                            index1,index2 = sorted([index1,index2])

                            # count nested cWW that reach outside of [index1,index2]
                            for i in range(index1+1,index2):
                                j = MCS_to_nested_cWW_endpoints[MCS1][i]

                                if not isinstance(j,tuple) and (j < index1 or j > index2):
                                    crossing += 1

                            if verbose >= 3 and crossing > 0:
                                print("%-20s and %-20s make %5s and have crossing number %3d" % (u1,u2,interaction,crossing))
                    else:
                        # different chains or different symmetries
                        # we know model1 == model2 already
                        # count nested pairs in chain1 that cross index1, in chain2 that cross index2
                        # Note: an interaction between chains could cross WC pairs between those
                        # chains, but we are not counting the crossing number for that
                        # In the same way, it would be hard to calculate a crossing number for
                        # crossing WC pairs that go between two chains.
                        for MCS, index in [(MCS1,index1),(MCS2,index2)]:
                            if MCS in MCS_to_nested_cWW_endpoints:
                                m = MCS_to_max_index[MCS]
                                if index < m / 2:
                                    # closer to the beginning of the chain
                                    for i in range(0,index):
                                        j = MCS_to_nested_cWW_endpoints[MCS][i]
                                        if not isinstance(j,tuple) and j > index:
                                            crossing += 1
                                else:
                                    # closer to the end of the chain
                                    for i in range(index+1,m+1):
                                        j = MCS_to_nested_cWW_endpoints[MCS][i]
                                        if not isinstance(j,tuple) and j < index:
                                            crossing += 1

                        # note inter-chain nested cWW basepairs for later bSS calculation
                        if crossing == 0 and interaction.lower() in ['cww','cwwa'] and 'bss' in categories:
                            # record canonical cWW pairs and their endpoints by chain
                            parent1 = get_parent_as_RNA(base1,'')
                            parent2 = get_parent_as_RNA(base2,'')

                            if parent1+parent2 in ['AU','UA','CG','GC','GU','UG']:
                                # here we map MCS and index to a tuple, to indicate it crosses chains
                                # that might be faster, but it does add a layer of complexity
                                MCS_to_nested_cWW_endpoints[MCS1][index1] = (MCS2,index2)
                                MCS_to_nested_cWW_endpoints[MCS2][index2] = (MCS1,index1)

                                bss_endpoints.add((MCS1,index1,MCS2,index2))

                pairs_to_crossing[(u1,u2)] = crossing

            if len(itpl) > 2:
                # change to list to be able to change or add the crossing number
                itpl_list = list(itpl)
                itpl_list[2] = crossing
                new_tuple = tuple(itpl_list)
            else:
                new_tuple = (u1,u2,crossing)
            interaction_to_list_of_tuples[interaction].append(new_tuple)

            # duplicate certain pairs in reversed order; saves time this way
            if interaction in ["s33","s35","s53","s55","ns33","ns35","ns53","ns55"]:
                interaction_to_list_of_tuples[reverse_edges(interaction)].append((u2,u1,crossing))
                pairs_to_crossing[(u2,u1)] = crossing
            elif interaction == 'cp':
                interaction_to_list_of_tuples[interaction].append((u2,u1,crossing))
                pairs_to_crossing[(u2,u1)] = crossing
            elif interaction[0] in ["c","t"]:
                interaction_to_list_of_tuples[reverse_edges(interaction)].append((u2,u1,crossing))
                pairs_to_crossing[(u2,u1)] = crossing
            elif interaction[0:2] in ["nc","nt"]:
                interaction_to_list_of_tuples[reverse_edges(interaction)].append((u2,u1,crossing))
                pairs_to_crossing[(u2,u1)] = crossing

    if 'bss' in categories:
        for MCS1,index1,MCS2,index2 in bss_endpoints:
            MCS_to_endpoints[MCS1].add(index1)
            MCS_to_endpoints[MCS2].add(index2)

        # make it easy to jump to the cww partner
        # set this up here to add the special ncww pairs identified in this section
        # use a list because sometimes there is more than one pairing partner, with ||A and ||B alternate ids
        unitid_to_cww_partner = defaultdict(set)

        # add nested AU, GC, GU ncWW pairs adjacent to canonical cWW pairs so we do not
        # extend a loop past a nested canonical ncWW pair.
        # There may be two in a row, so check until no more additions
        # Unfortunately this adds 4V9F|1|0|C|1558 ncWW 4V9F|1|0|G|1563 and shortens the hairpin
        pair_added = True
        while pair_added:
            pair_added = False
            for interaction in interaction_to_list_of_tuples.keys():
                if interaction.lower() in ['ncww','ncwwa']:
                    for u1,u2,crossing in interaction_to_list_of_tuples[interaction]:
                        if crossing == 0:
                            model1,chain1,index1,base1,symmetry1 = unit_id_to_fields[u1]
                            model2,chain2,index2,base2,symmetry2 = unit_id_to_fields[u2]

                            if not index1:
                                continue
                            if not index2:
                                continue

                            parent1 = get_parent_as_RNA(base1,'')
                            parent2 = get_parent_as_RNA(base2,'')

                            if parent1+parent2 in ['AU','UA','CG','GC','GU','UG']:
                                MCS1 = (model1,chain1,symmetry1)
                                MCS2 = (model2,chain2,symmetry2)

                                if (MCS1,index1+1,MCS2,index2-1) in bss_endpoints \
                                or (MCS1,index1-1,MCS2,index2+1) in bss_endpoints:
                                    if not index1 in MCS_to_endpoints[MCS1] and not index2 in MCS_to_endpoints[MCS2]:
                                        # if (u1,u2,0) in interaction_to_list_of_tuples.get('cp',[]):
                                            if MCS1 == MCS2:
                                                MCS_to_nested_cWW_endpoints[MCS1][index1] = index2
                                                MCS_to_nested_cWW_endpoints[MCS2][index2] = index1
                                            else:
                                                MCS_to_nested_cWW_endpoints[MCS1][index1] = (MCS2,index2)
                                                MCS_to_nested_cWW_endpoints[MCS2][index2] = (MCS1,index1)

                                            MCS_to_endpoints[MCS1].add(index1)
                                            MCS_to_endpoints[MCS2].add(index2)

                                            bss_endpoints.add((MCS1,index1,MCS2,index2))

                                            unitid_to_cww_partner[u1].add(u2)
                                            unitid_to_cww_partner[u2].add(u1)

                                            pair_added = True

                                            # if index1 <= index2:
                                            #     print("Adding %s %s %s to the bSS endpoints" % (u1,interaction,u2))

        # identify nucleotides that border a single-stranded region (bSS relation)
        # this includes strands between nested cWW pairs on one chain
        # it also includes cWW pairs with zero crossing number that go between chains
        # Note that sometimes there are nucleotides between u1 and u2 in the chain
        # that do not have xyz coordinates.  We still allow u1 and u2 to have the bSS relation,
        # because there is no intervening cWW basepair (but with full xyz coordinates there could be).
        unit_id_pair_to_interaction = {}

        # use a list because sometimes there is more than one pairing partner, with ||A and ||B alternate ids
        unitid_to_bss_partner = defaultdict(set)

        bSS_list = []
        bSS_list_ordered = []
        for MCS, endpoints in MCS_to_endpoints.items():
            if verbose >= 3:
                print("  Getting bSS for %s %s %s" % MCS)

            # c is the "lower" index; we increase it in this process
            c = MCS_to_min_index[MCS]

            if not c == 1:
                if verbose >=3:
                    print('  Minimum index is %d' % c)

            # highest index in the chain is also an endpoint
            all_endpoints = sorted(endpoints - set([c])) + [MCS_to_max_index[MCS]]
            pairing_partners = MCS_to_nested_cWW_endpoints[MCS]

            # walk through the chain, looking for single-stranded regions
            # e is the "upper" index; we increase it in this process
            for e in all_endpoints:

                pc = pairing_partners[c]
                MCS2 = None
                if isinstance(pc,tuple):
                    # cWW pair with crossing number 0 between chains
                    MCS2 = pc[0]
                    pc = pc[1]

                pe = pairing_partners[e]
                MCS3 = None
                if isinstance(pe,tuple):
                    # cWW pair with crossing number 0 between chains
                    MCS3 = pe[0]
                    pe = pe[1]

                u1 = MCS_index_to_unit_id[MCS][c]
                u2 = MCS_index_to_unit_id[MCS][e]

                if c < 100 and verbose >= 3:
                    print('  Thinking about bSS between %-20s and %-20s' % (u1,u2), end=" ")
                    print('  c is %d, pc is %d, e is %d, pe is %d' % (c,pc,e,pe))

                if (MCS2 or MCS3) and not MCS2 == MCS3:
                    # c and e are in one chain, but pc and pe are in different chains
                    bSS_list.append((u1,u2,0))
                    bSS_list.append((u2,u1,0))
                    bSS_list_ordered.append((u1,u2))
                    unitid_to_bss_partner[u1].add(u2)

                    if verbose >= 3:
                        if c == MCS_to_min_index[MCS]:
                            print("  %-20s bSS %-20s at start of chain &" % (u1,u2))
                        elif e == MCS_to_max_index[MCS]:
                            print("  %-20s bSS %-20s at end of chain &" % (u1,u2))
                        else:
                            print("  %-20s bSS %-20s between chains" % (u1,u2))
                elif c == MCS_to_min_index[MCS] and pc == c and not MCS2:
                    # c is at start of chain but does not make a cWW pair
                    if verbose >= 3:
                        print("  %-20s bSS %-20s at start of chain *" % (u1,u2))
                        print(c,pc,e,pe,MCS,MCS2,MCS3)
                    bSS_list.append((u1,u2,0))
                    bSS_list.append((u2,u1,0))
                    bSS_list_ordered.append((u1,u2))
                    unitid_to_bss_partner[u1].add(u2)
                elif e == MCS_to_max_index[MCS] and pe == e and not MCS3:
                    # e is at end of chain but does not make a cWW pair
                    if verbose >= 3:
                        print("  %-20s bSS %-20s at end of chain *" % (u1,u2))
                    bSS_list.append((u1,u2,0))
                    bSS_list.append((u2,u1,0))
                    bSS_list_ordered.append((u1,u2))
                    unitid_to_bss_partner[u1].add(u2)
                elif e - c > 1:
                    # space between cWW basepairs
                    if pc-pe == e-c and ((not pc == e and not MCS2) or (MCS2)):
                        # symmetric internal loop, check if complementary and no interactions
                        if MCS2:
                            MCS1 = MCS2
                        else:
                            MCS1 = MCS
                        complementary = True
                        opposite_pairs = []

                        for i in range(c+1,e):
                            v1 = MCS_index_to_unit_id[MCS].get(i,None)
                            if not v1:
                                continue
                            v2 = MCS_index_to_unit_id[MCS1].get(pc-(i-c),None)
                            if not v2:
                                continue
                            parent1 = get_parent_as_RNA(v1.split("|")[3],'')
                            parent2 = get_parent_as_RNA(v2.split("|")[3],'')
                            if parent1+parent2 in ['AU','UA','CG','GC','GU','UG']:
                                opposite_pairs.append((v1,v2))
                            else:
                                complementary = False
                                break

                        if complementary:
                            # check for ncWW interactions at least somewhere
                            found_ncWW = False
                            found_other_interaction = False
                            # build pair to interaction mapping if not already built
                            if not unit_id_pair_to_interaction:
                                for interaction in interaction_to_pair_list.keys():
                                    if interaction.lower().replace("n","").replace("a","") in Leontis_Westhof_basepairs_lower:
                                        for v1,v2 in interaction_to_pair_list[interaction]:
                                            unit_id_pair_to_interaction[(v1,v2)] = interaction

                            for (v1,v2) in opposite_pairs:
                                # check in both orders because pairs are not stored both ways yet
                                if (v1,v2) in unit_id_pair_to_interaction:
                                    interaction = unit_id_pair_to_interaction[(v1,v2)]
                                    if interaction.lower() in ['ncww','ncwwa']:
                                        found_ncWW = True
                                    else:
                                        found_other_interaction = True
                                elif (v2,v1) in unit_id_pair_to_interaction:
                                    interaction = unit_id_pair_to_interaction[(v2,v1)]
                                    if interaction.lower() in ['ncww','ncwwa']:
                                        found_ncWW = True
                                    else:
                                        found_other_interaction = True

                            # if not found_ncWW or found_other_interaction:
                            if found_other_interaction:
                                complementary = False
                                for (v1,v2) in opposite_pairs:
                                    int1 = unit_id_pair_to_interaction.get((v1,v2),None)
                                    int2 = unit_id_pair_to_interaction.get((v2,v1),None)
                                    if verbose >= 3:
                                        print("  Found complementary pair %s and %s making %s or %s" % (v1,v2,int1,int2))
                                if found_other_interaction:
                                    if verbose >= 3:
                                        print('  Found other interaction in this IL')
                                if verbose >= 3:
                                    print('  Recording bSS between %s and %s' % (u1,u2))
                            else:
                                pass
                                # for pair in opposite_pairs:
                                #     print("  Found complementary pair %s and %s" % pair)
                                # print('  Not recording bSS between %s and %s' % (u1,u2))

                        if not complementary:
                            bSS_list.append((u1,u2,0))
                            bSS_list.append((u2,u1,0))
                            bSS_list_ordered.append((u1,u2))
                            unitid_to_bss_partner[u1].add(u2)
                            if verbose >= 3:
                                print('  %-20s bSS %-20s from symmetric IL' % (u1,u2))

                    else:
                        bSS_list.append((u1,u2,0))
                        bSS_list.append((u2,u1,0))
                        bSS_list_ordered.append((u1,u2))
                        unitid_to_bss_partner[u1].add(u2)
                        if verbose >= 3:
                            print("  %-20s bSS %-20s gap between cWW's" % (u1,u2))
                elif abs(pc-pe) > 1:
                    if verbose >= 3:
                        if c == MCS_to_min_index[MCS]:
                            print('  %-20s bSS %-20s at start of chain #' % (u1,u2))
                        else:
                            print('  %-20s bSS %-20s distance 1 apart' % (u1,u2))
                    bSS_list.append((u1,u2,0))
                    bSS_list.append((u2,u1,0))
                    bSS_list_ordered.append((u1,u2))
                    unitid_to_bss_partner[u1].add(u2)

                # move up the "lower" index
                c = e

            if verbose >= 3:
                print("  Last index is %d" % c)
                print("  Max  index is %s" % MCS_to_max_index[MCS])

        if len(bSS_list) > 0:
            interaction_to_list_of_tuples['bSS'] = bSS_list

    if 'loops' in categories:
        # map unit ids making nested cWW to their pairing partners
        # this may pick up cWW pairs that are not GC, AU, or GU but
        # that is OK because we are only going to use the ones that
        # come from the bSS relation as we walk along the chain
        for interaction in sorted(interaction_to_list_of_tuples.keys(),reverse=True):
            if interaction.lower() in ['cww','cwwa']:
                for unitid1, unitid2, crossing in interaction_to_list_of_tuples[interaction]:
                    if crossing == 0:
                        unitid_to_cww_partner[unitid1].add(unitid2)
                        unitid_to_cww_partner[unitid2].add(unitid1)

        # identify the flanking unit ids of hairpin, internal, junction loops
        all_loops = []
        bss_done = set()
        loop_counter = defaultdict(lambda: 0)

        # should sort interaction_to_list_of_tuples['bSS'] by alt id of the first nucleotide
        # also exclude any bSS pairs from alt id A to alt id B
        bSS_tuples = []
        c = 1
        for unitid1, unitid2 in bSS_list_ordered:
            alt_id1 = get_alt_id(unitid1)
            alt_id2 = get_alt_id(unitid2)
            # use a counter to keep the list in the original order as much as possible
            # that way, we get bSS nucleotides in index order on each chain
            bSS_tuples.append((alt_id1,alt_id2,c,unitid1,unitid2))
            c += 1

        bSS_tuples = sorted(bSS_tuples)

        # for alt_id1, alt_id2, c, unitid1, unitid2 in bSS_tuples:
        #     print('bSS tuple',alt_id1,alt_id2,c,unitid1,unitid2)
        # for u,p in unitid_to_bss_partner.items():
        #     print('bss partner',u,p)

        for alt_id1, alt_id2, c, unitid1, unitid2 in bSS_tuples:
            if not unitid1 in unitid_to_cww_partner:
                continue

            if not unitid2 in unitid_to_cww_partner:
                continue

            if (unitid1,unitid2) in bss_done:
                continue

            if (unitid2,unitid1) in bss_done:
                continue

            # print('Starting with bSS pair',unitid1,unitid2)

            bss_done.add((unitid1,unitid2))
            already_used = set()

            # store unit ids in a list of length L called loop
            # L is even, equal to 2 or 4 or 6 or more, with this structure:
            # positions 0 and 1 are bSS partners
            # positions 1 and 2 are cWW partners
            # positions 2 and 3 are bSS partners, and so on
            # positions L-1 and 0 are cWW partners
            if unitid1 in unitid_to_cww_partner[unitid2]:
                # hairpin loop
                loop = [unitid1,unitid2]
            else:
                # walk from unitid1 to unitid2, collecting all nucleotides
                loop = [unitid1,unitid2]
                start = unitid1
                unitid1 = ''
                # follow cWW partner and get its bSS partner, until you get back to the start
                while start != unitid1:
                    # print('  building on ',loop)
                    unitid1 = select_partner(unitid_to_cww_partner[unitid2],already_used)
                    # print('  unitid1',unitid1)
                    if unitid1:
                        unitid2 = select_partner(unitid_to_bss_partner[unitid1],already_used)
                        # print('  unitid2',unitid2)
                        if unitid2:
                            if unitid1 in loop or unitid2 in loop:
                                # bad situation, returned to the loop but not the start
                                # need to avoid infinite while loop
                                print("  Anomalous loop %s" % ",".join(loop))
                                loop = []
                                break
                            loop.append(unitid1)
                            loop.append(unitid2)
                            # print('  Next unitid1 bSS unitid2 are %s and %s' % (unitid1,unitid2))
                            bss_done.add((unitid1,unitid2))
                            unitid1 = select_partner(unitid_to_cww_partner[unitid2],already_used)
                            if unitid1:
                                # print('  Next unitid1 will be %s' % (unitid1))
                                pass
                            else:
                                loop = []
                                break
                        else:
                            loop = []
                            break
                    else:
                        loop = []
                        break

            # print('Final loop:',loop)

            if loop:
                full_loop, loop_counter = fill_in_strands_of_loop(loop,unit_id_to_fields,MCS_index_to_unit_id,loop_counter)
                if full_loop:
                    full_loop['merged_from'] = [len(all_loops)]
                    all_loops.append(full_loop)

                    if verbose >= 3:
                        for i, unitid in enumerate(full_loop['unit_ids']):
                            print('  %s %2s %s %s' % (full_loop['type'],i,full_loop['border_indicators'][i],unitid))
                        print()

        # map flanking cWW pairs to list of loop indices, to find shared cWW pairs
        flanking_cWW_pair_to_loop = defaultdict(list)
        for a, full_loop in enumerate(all_loops):
            loop = full_loop['border_unit_ids']
            for i in range(len(loop)):
                if i % 2 == 1:
                    u1 = loop[i]
                    if i == len(loop)-1:
                        u2 = loop[0]
                    else:
                        u2 = loop[i+1]
                    if u1 < u2:
                        flanking_cWW_pair_to_loop[(u1,u2)].append((a,u1,u2))
                    else:
                        flanking_cWW_pair_to_loop[(u2,u1)].append((a,u1,u2))

        # merge loops that share cWW pairs and have interactions across that pair
        # some loops, like IL_4V9F_008 and IL_4V9F_009, are separated by just
        # one cWW pair and they make interactions across that pair.
        # In a 2d structure, they look like two IL, but in 3D they form a single IL
        # with an embedded cWW basepair.
        # We extract them both ways; thus some loop strands overlap

        # which cWW pairs have more than one loop sharing them?
        cWW_pairs_to_check = []
        for cWW_pair, triple_list in flanking_cWW_pair_to_loop.items():
            if len(triple_list) == 2:
                cWW_pairs_to_check.append(cWW_pair)
            elif len(triple_list) > 2:
                if verbose >= 3:
                    print('  Too many loops between %s and %s' % cWW_pair)
                    for index, v1, v2 in triple_list:
                        print(all_loops[index])

        # if some loops share a cWW pair, check if they should be merged
        if len(cWW_pairs_to_check) > 0:
            # record unit id pairs to interactions
            unit_id_pair_to_interaction = {}
            if verbose >= 3:
                print("Setting up unit_id_pair_to_interaction")
            for interaction in sorted(interaction_to_pair_list.keys()):
                # print("  Processing %s" % interaction)
                # if interaction in ["s33","s35","s53","s55"]+Leontis_Westhof_basepairs:
                # only check basepairs; stacks are not significant enough to merge
                if interaction in Leontis_Westhof_basepairs:
                    for v1,v2 in interaction_to_pair_list[interaction]:
                        unit_id_pair_to_interaction[(v1,v2)] = interaction

            # record stacking partners for each nucleotide face
            unitid_face_to_stacking_partners = defaultdict(set)
            face_to_stacks = {}
            face_to_stacks['3'] = ['s33','s35','ns33','ns35']
            face_to_stacks['5'] = ['s53','s55','ns53','ns55']
            face_to_stacks['3'] = ['s33','s35']
            face_to_stacks['5'] = ['s53','s55']
            for face, face_list in face_to_stacks.items():
                for interaction in face_list:
                    for v1,v2,crossing in interaction_to_list_of_tuples[interaction]:
                        unitid_face_to_stacking_partners[(v1,face)].add(v2)

            # check for shared cWW pairs and merge, until all merges are complete
            while len(cWW_pairs_to_check) > 0:
                # the current cWW pair is shared by two loops; should they be merged?
                (u1,u2) = cWW_pairs_to_check.pop(0)
                index0, a1, a2 = flanking_cWW_pair_to_loop[(u1,u2)][0]
                index1, b1, b2 = flanking_cWW_pair_to_loop[(u1,u2)][1]

                loop0 = all_loops[index0]
                loop1 = all_loops[index1]

                if verbose >= 3:
                    print('  Found %s and %s sharing cWW pair %s and %s' % (loop0['type'],loop1['type'],u1,u2))

                loop0_unit_ids = set(loop0['unit_ids']) - set([u1,u2])
                loop1_unit_ids = set(loop1['unit_ids']) - set([u1,u2])

                merge_loops = False

                # is there a basepair between the two loops?
                for v1 in loop0_unit_ids:
                    for v2 in loop1_unit_ids:
                        if (v1,v2) in unit_id_pair_to_interaction:
                            interaction = unit_id_pair_to_interaction[(v1,v2)]
                            merge_loops = True
                            if verbose >= 3:
                                print('    Found interaction %s between %s and %s' % (interaction,v1,v2))
                        if (v2,v1) in unit_id_pair_to_interaction:
                            interaction = unit_id_pair_to_interaction[(v2,v1)]
                            merge_loops = True
                            if verbose >= 3:
                                print('    Found interaction %s between %s and %s' % (interaction,v2,v1))

                # is there a stacking interaction to the opposite side of the cWW pair?
                # this happens often enough that we need to check for it; e.g., IL_4V9F_008 and IL_4V9F_009
                # we use a1, a2, b1, b2 to keep flanking nucleotides in the correct
                # relationship to the loop they came from
                loop0_unit_ids = loop0_unit_ids - set(loop0['border_unit_ids'])
                loop1_unit_ids = loop1_unit_ids - set(loop1['border_unit_ids'])

                intersect = unitid_face_to_stacking_partners[(a1,'5')] & loop1_unit_ids
                if len(intersect) > 0:
                    merge_loops = True
                    if verbose >= 3:
                        print("    Found 5' face of %s stacking on %s" % (a1,str(intersect)))
                intersect = unitid_face_to_stacking_partners[(a2,'3')] & loop1_unit_ids
                if len(intersect) > 0:
                    merge_loops = True
                    if verbose >= 3:
                        print("    Found 3' face of %s stacking on %s" % (a2,str(intersect)))
                intersect = unitid_face_to_stacking_partners[(b1,'5')] & loop0_unit_ids
                if len(intersect) > 0:
                    merge_loops = True
                    if verbose >= 3:
                        print("    Found 5' face of %s stacking on %s" % (b1,str(intersect)))
                intersect = unitid_face_to_stacking_partners[(b2,'3')] & loop0_unit_ids
                if len(intersect) > 0:
                    merge_loops = True
                    if verbose >= 3:
                        print("    Found 3' face of %s stacking on %s" % (b2,str(intersect)))

                if merge_loops:
                    ids0 = loop0['border_unit_ids']
                    ids1 = loop1['border_unit_ids']
                    new_loop = []

                    # put u1,u2 pair at start and end of loop0 and loop1
                    # ids0 = keep_loop_ids(ids0,u1,u2)
                    # ids1 = keep_loop_ids(ids1,u1,u2)

                    if verbose >= 3:
                        print('    Merging loops')
                        print('    ids0 is %s' % ids0)
                        print('    ids1 is %s' % ids1)
                        print('    cWW pair to merge on is %s to %s' % (u1,u2))
                        print('    ids0 is %s' % ids0)
                        print('    ids1 is %s' % ids1)

                    if len(ids0) == 2:
                        # easy to merge HL into IL or junction
                        new_loop = [id for id in ids1 if not id in [u1,u2]]
                    elif len(ids1) == 2:
                        new_loop = [id for id in ids0 if not id in [u1,u2]]
                    else:
                        # neither loop is a hairpin
                        # accumulate ids from loop0 up to u1,u2
                        i = 0
                        while i < len(ids0) and not ids0[i] in [u1,u2]:
                            new_loop.append(ids0[i])
                            i += 1
                        # go to after ids0[i] in loop1
                        j = 0
                        while j < len(ids1) and not ids1[j] == ids0[i]:
                            j += 1
                        # accumulate ids from loop1
                        j += 1
                        while j < len(ids1) and not ids1[j] in [u1,u2]:
                            new_loop.append(ids1[j])
                            j += 1
                        # accumulate ids from loop1
                        j = 0
                        while not ids1[j] in [u1,u2]:
                            j += 1
                            new_loop.append(ids1[j])
                        # accumulate ids from loop0
                        if i > 0:
                            # u1, u2 are in the middle of loop0
                            i += 2
                        else:
                            # u1, u2 are at the start and end of loop1
                            i = 1
                        while i < len(ids0) and not ids0[i] in [u1,u2]:
                            new_loop.append(ids0[i])
                            i += 1

                    if verbose >= 3:
                        print("    new loop is: %s" % new_loop)

                    full_loop, loop_counter = fill_in_strands_of_loop(new_loop,unit_id_to_fields,MCS_index_to_unit_id,loop_counter)
                    if full_loop:
                        full_loop['merged_from'] = loop0['merged_from'] + loop1['merged_from']
                        full_loop['identifier'] = full_loop['identifier'].replace('single','merged')
                        all_loops.append(full_loop)

                        loop_index = len(all_loops)-1

                        # point to new loop for any additional merges
                        for cWW_pair, triple_list in flanking_cWW_pair_to_loop.items():
                            for i, t in enumerate(triple_list):
                                if t[0] == index0 or t[0] == index1:
                                    flanking_cWW_pair_to_loop[cWW_pair][i] = (loop_index,t[1],t[2])

                        if verbose >= 3:
                            for i, unitid in enumerate(full_loop['unit_ids']):
                                print('  %s %2s %s %s' % (full_loop['type'],i,full_loop['border_indicators'][i],unitid))
                            print()

                            if len(full_loop['merged_from']) > 2:
                                print('  Merged from these loop indices: %s' % full_loop['merged_from'])

        # stuff the loops in here, even though they don't fit the rest of the pattern
        interaction_to_list_of_tuples['loops'] = all_loops

    return interaction_to_list_of_tuples


def fill_in_strands_of_loop(loop,unit_id_to_fields,MCS_index_to_unit_id,loop_counter):
    """
    loop is a list of unit ids of flanking pairs
    loop[0] makes bSS with loop[1]
    loop[1] makes cWW with loop[2]
    loop[2] makes bSS with loop[3], and so on
    loop[-1] makes cWW with loop[0]
    This function fills in the nucleotides between the flanking pairs
    It returns a dictionary with keys 'type', 'unit_ids', 'border_unit_ids', 'border_indicators'
    """
    missing_index = []
    all_unit_ids = []
    all_border_indicators = []

    if len(loop) %2 == 1:
        # do not create a loop when there are an odd number of nucleotides
        # should not happen, but sometimes it does, and there just isn't time
        # to debug every case like that
        return {}, loop_counter

    for i, unitid in enumerate(loop):
        all_unit_ids.append(unitid)
        all_border_indicators.append('1')

        fields = unit_id_to_fields[unitid]
        MCS = (fields[0],fields[1],fields[4])

        if i % 2 == 0:
            index_low = unit_id_to_fields[unitid][2]
            index_high = unit_id_to_fields[loop[i+1]][2]
            for index in range(index_low+1,index_high):
                if index in MCS_index_to_unit_id[MCS]:
                    all_unit_ids.append(MCS_index_to_unit_id[MCS][index])
                    all_border_indicators.append('0')
                else:
                    missing_index.append(index)

    if len(missing_index) > 0:
        # do not create a loop when there are missing nucleotides
        return {}, loop_counter
    else:
        if len(loop) == 2:
            loop_type = 'HL'
        elif len(loop) == 4:
            loop_type = 'IL'
        else:
            j = len(loop) // 2
            loop_type = 'J%d' % j

        loop_counter[loop_type] += 1
        file_id = loop[0].split("|")[0]

        full_loop = {}
        full_loop['type'] = loop_type
        full_loop['border_unit_ids'] = loop
        full_loop['unit_ids'] = all_unit_ids
        full_loop['border_indicators'] = all_border_indicators
        full_loop['identifier'] = '%s_%s_%03d_single' % (loop_type,file_id,loop_counter[loop_type])
    return full_loop, loop_counter

def annotate_covalent_connections(nucleotides, interaction_to_list_of_tuples, category_to_interactions, timerData):
    """
    Loop through bases, sort by model, symmetry, chain, and
    record the distance in the chain between successive
    observed nucleotides.
    """

    nts_to_sort = defaultdict(list)

    # list the nucleotides from chains by model, symmetry, chain, index
    for nt in nucleotides:
        if nt.index:
            nts_to_sort[(nt.model+" "+nt.symmetry+" "+nt.chain,nt.index)].append(nt.unit_id())

    sorted_keys = sorted(nts_to_sort.keys())

    for i in range(0,len(sorted_keys)-1):
        key1 = sorted_keys[i]
        key2 = sorted_keys[i+1]

        # same model, symmetry, chain
        if key1[0] == key2[0]:
            chain_distance = key2[1]-key1[1]
            for u1 in nts_to_sort[key1]:
                for u2 in nts_to_sort[key2]:
                    interaction = "p_" + str(chain_distance)
                    interaction_to_list_of_tuples[interaction].append((u1,u2,0))
                    category_to_interactions["covalent"].add(interaction)

    return interaction_to_list_of_tuples, category_to_interactions, timerData


def annotate_nt_nt_in_structure(structure,categories,focused_basepair_cutoffs={},ideal_hydrogen_bonds={},chains=[],timerData=None,get_datapoint=False):
    """
    This function can be called from the pipeline to annotate a structure
    structure is an output from
    """

    if not focused_basepair_cutoffs and 'basepair' in categories:
        focused_basepair_cutoffs = focus_basepair_cutoffs(nt_nt_cutoffs,categories['basepair'])

    if not ideal_hydrogen_bonds and 'basepair' in categories:
        ideal_hydrogen_bonds = load_ideal_basepair_hydrogen_bonds()

    # structures.py controls what residues are returned, that sometimes needs to be expanded
    # SOLITARY means to check non-NA chains for nucleotides like ATP.  That adds time.
    if chains:
        # bases = structure.residues(chain = chains, type = ["RNA linking","DNA linking"])  # load all RNA/DNA nucleotides
        bases = structure.residues(chain = chains, type = ["RNA","DNA","PNA","SOLITARY"])  # load all RNA/DNA nucleotides from desired chains
    else:
        # bases = structure.residues(type = ["RNA linking","DNA linking"])  # load nice RNA/DNA nucleotides
        bases = structure.residues(type = ["RNA","DNA","PNA","SOLITARY"])  # load all RNA/DNA nucleotides

    # for base in bases:
    #     print(base.unit_id())

    if not timerData:
        timerData = myTimer("start")

    # maximum center-center distance to check for interactions
    nt_nt_screen_distance = standard_center_center_distance_cutoff
    if 'backbone' in categories:
        nt_nt_screen_distance = base_backbone_center_center_distance_cutoff
    if 'oo_distance' in categories:
        nt_nt_screen_distance = oo_distance_center_center_distance_cutoff

    timerData = myTimer("Build cubes for neighbors",timerData)
    baseCubeList, baseCubeNeighbors = make_nt_cubes_half(bases, nt_nt_screen_distance, nt_reference_point)
    # annotate nt-nt interactions
    interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_interactions(bases, nt_nt_screen_distance, baseCubeList, baseCubeNeighbors, categories, focused_basepair_cutoffs, ideal_hydrogen_bonds, timerData, get_datapoint)

    # annotate covalent connections
    interaction_to_list_of_tuples, category_to_interactions, timerData = annotate_covalent_connections(bases, interaction_to_list_of_tuples, category_to_interactions, timerData)

    return interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data


def get_parent(sequence,if_none=None):
    """
    Look up parent sequence for RNA, DNA, and modified nucleotides.
    Return A, C, G, U, DT for cases that treat DT differently than U
    """

    if sequence in ['A','C','G','U','DT']:
        return sequence
    elif sequence in ['DA','DC','DG']:
        return sequence[1]
    elif sequence in modified_base_to_parent:
        parent = modified_base_to_parent[sequence]
        if parent in ['A','C','G','U','DT']:
            return parent
        elif parent in ['DA','DC','DG']:
            return parent[1]
    return if_none


def get_parent_as_RNA(sequence,if_none=None):
    """
    Look up parent sequence for RNA, DNA, and modified nucleotides.
    Return A, C, G, U to make it easier
    """

    if sequence in ['A','C','G','U']:
        return sequence
    elif sequence in ['DA','DC','DG']:
        return sequence[1]
    elif sequence in ['T','DT']:
        return 'U'
    elif sequence in modified_base_to_parent:
        parent = modified_base_to_parent[sequence]
        if parent in ['A','C','G','U']:
            return parent
        elif parent in ['DA','DC','DG']:
            return parent[1]
        elif parent == 'DT':
            return 'U'
    return if_none


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


def check_base_oxygen_stack_rings(nt1,nt2,parent1,datapoint):
    """
    Does one of the backbone oxygens of nt2 stack inside a ring on the base of nt1?
    """

    true_z_cutoff = 3.5
    near_z_cutoff = 3.6

    interaction = ""
    interaction_reversed = ""

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
            interaction_reversed = "s" + oxygenmin + "3"
        else:
            interaction = "s5" + oxygenmin
            interaction_reversed = "s" + oxygenmin + "5"

    elif near_found:  # over a base ring, but z value too large for true
        if zmin > 0:
            interaction = "ns3" + oxygenmin
            interaction_reversed = "ns" + oxygenmin + "3"
        else:
            interaction = "ns5" + oxygenmin
            interaction_reversed = "ns" + oxygenmin + "5"

    else:            # not over a base ring, but maybe close enough

        r2min = 999   # keep track of distance to base center, use the minimum

        for x,y,z,oxygen in oxygen_points:

            nearring5 = False
            nearring6 = False

            # check ellipses only
            if abs(z) < true_z_cutoff:
                if parent1 == 'A' or parent1 == 'DA':
                    if -1.302671*x + -0.512161*y + -0.512114 > 0:  # Left of C4-C5
                        if 1.033454*(x-(-1.138126))**2 + 0.143656*(x-(-1.138126))*(y-(-0.650781)) + (y-(-0.650781))**2 < 2.163590:  # A5 r=0.3
                            nearring5 = True
                    else:
                        if 1.001608*(x-(0.850305))**2 + 0.169100*(x-(0.850305))*(y-(-0.017921)) + (y-(-0.017921))**2 < 2.766745:  # A6 r=0.3
                            nearring6 = True
                elif parent1 == 'C' or parent1 == 'DC':
                    if 0.867183*(x-(-0.298275))**2 + 0.040055*(x-(-0.298275))*(y-(-0.153209)) + (y-(-0.153209))**2 < 2.652492:  # C r=0.3
                        nearring6 = True
                elif parent1 == 'G' or parent1 == 'DG':
                    if -1.306197*x + -0.492373*y + -0.896488 > 0:  # Left of C4-C5
                        if 1.032607*(x-(-1.476126))**2 + 0.129895*(x-(-1.476126))*(y-(-0.541964)) + (y-(-0.541964))**2 < 2.157145:  # G5 r=0.3
                            nearring5 = True
                    else:
                        if 1.082495*(x-(0.521747))**2 + 0.260413*(x-(0.521747))*(y-(0.023305)) + (y-(0.023305))**2 < 2.920747:  # G6 r=0.3
                            nearring6 = True
                elif parent1 == 'DT':
                    if 0.959551*(x-(0.029169))**2 + 0.128151*(x-(0.029169))*(y-(-0.304375)) + (y-(-0.304375))**2 < 2.766276:  # DT r=0.3
                        nearring6 = True
                elif parent1 == 'U':
                    if 0.912164*(x-(-0.302801))**2 + 0.143626*(x-(-0.302801))*(y-(-0.157137)) + (y-(-0.157137))**2 < 2.752991:  # U r=0.3
                        nearring6 = True

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
                        ringmin = "near_ring5"
                    else:
                        ringmin = "near_ring6"

        if near_found:
            if zmin > 0:
                interaction = "ns3" + oxygenmin
                interaction_reversed = "ns" + oxygenmin + "3"
            else:
                interaction = "ns5" + oxygenmin
                interaction_reversed = "ns" + oxygenmin + "5"

    if False and len(interaction) > 0:
        print('%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\t\t=hyperlink("https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s")' % (nt1.unit_id(),nt2.unit_id(),interaction,xmin,ymin,zmin,nt1.unit_id(),nt2.unit_id()))

    if datapoint:
        if len(interaction) > 0:
            datapoint['sOx'] = xmin
            datapoint['sOy'] = ymin
            datapoint['sOz'] = zmin
            datapoint['sOring'] = ringmin
            datapoint['sOoxygen'] = oxygenmin
            datapoint['sOinteraction'] = interaction

    return interaction, datapoint, interaction_reversed

def check_convex_hull_atoms(x,y,z, parent):
    """Method to check and see if an atom that has been translated into standard orientation falls within the
    convex hull of a nucleotide base based on the type of nucleotide (A,C,G,U,DT). Numbers Generated From generate_location_checks.py and
    this method was developed for use in the check_base_base_stacking method.
    Takes in two nucleotides coordinates alongside the nucleotide type.
    Returns True if The points fall within the convex hull of the parent1 and returns False Otherwise."""
    near_z_cutoff = 4.5
    inside = False
    if abs(z) < near_z_cutoff:
        if parent == 'A' or parent == 'DA':
            if -2.327244*x +  4.271447*y +  9.515028 > 0:  # Left of H9'-H2
                if -3.832809*x + -2.350503*y + 10.927316 > 0:  # Left of H2-H61
                    if  0.451014*x + -1.690509*y +  5.259508 > 0:  # Left of H61-H62
                        if  4.252574*x + -2.330898*y + 10.447200 > 0:  # Left of H62-H8
                            if  1.456465*x +  2.100463*y +  7.567280 > 0:  # Left of H8-H9'
                                inside = True
        elif parent == 'C' or parent == 'DC':
            if -0.889476*x +  2.269450*y +  5.323403 > 0:  # Left of H1'-O2
                if -4.532779*x + -1.065616*y +  6.851131 > 0:  # Left of O2-H42
                    if -0.190206*x + -1.731804*y +  5.226294 > 0:  # Left of H42-H41
                        if  1.955107*x + -1.508802*y +  6.480180 > 0:  # Left of H41-H5
                            if  2.523463*x + -0.045961*y +  6.153526 > 0:  # Left of H5-H6
                                if  1.133891*x +  2.082733*y +  5.627625 > 0:  # Left of H6-H1'
                                    inside = True
        elif parent == 'G' or parent == 'DG':
            if -1.107310*x +  4.872516*y + 11.647152 > 0:  # Left of H9'-H21
                if -1.684502*x +  0.422659*y +  6.436199 > 0:  # Left of H21-H22
                    if -1.592264*x + -1.681840*y +  6.230291 > 0:  # Left of H22-H1
                        if -1.019666*x + -2.216349*y +  5.884100 > 0:  # Left of H1-O6
                            if  2.274081*x + -2.148378*y +  5.898397 > 0:  # Left of O6-N7
                                if  1.656548*x + -1.350181*y +  4.208981 > 0:  # Left of N7-H8
                                    if  1.473113*x +  2.101573*y +  7.865111 > 0:  # Left of H8-H9'
                                        inside = True

        elif parent == 'U':
            if -0.960553*x +  2.292490*y +  5.471254 > 0:  # Left of H1'-O2
                if -2.493573*x + -0.200338*y +  4.589448 > 0:  # Left of O2-H3
                    if -1.574881*x + -1.914996*y +  4.563214 > 0:  # Left of H3-O4
                        if  1.403523*x + -2.301733*y +  5.976805 > 0:  # Left of O4-H5
                            if  2.504701*x +  0.041797*y +  6.092950 > 0:  # Left of H5-H6
                                if  1.120783*x +  2.082780*y +  5.621468 > 0:  # Left of H6-H1'
                                    inside = True
        elif parent == 'DT':
            if -1.125648*x +  2.281277*y +  6.199955 > 0:  # Left of C1'-O2
                if -2.368105*x + -0.456021*y +  4.878252 > 0:  # Left of O2-H3
                    if -1.526233*x + -1.897795*y +  4.450270 > 0:  # Left of H3-O4
                        if  1.301401*x + -2.544887*y +  5.949759 > 0:  # Left of O4-C7
                            if  2.031505*x +  1.412190*y +  3.691439 > 0:  # Left of C7-C6
                                if  1.687080*x +  1.205236*y +  3.097805 > 0:  # Left of C6-C1'
                                    inside = True
        else:
            if verbose >= 2:
                print("  Unrecognized parent " + parent + " in function check_convex_hull_atoms. FR3D is currently unable to recognize this modified base.")
            return False
    return inside

def return_overlap(listOfAtoms, nt1, nt2, parent):
    """Function to check if there's overlap between a list of atoms from one nucleotide and the atoms of the base of another
    list of base atoms of nt2, 2 nucleotides, and the parent of nt1 are passed in.
    Checks each atom in the list of atoms. Takes its coordinates and translates them to be in respect to nt1 in standard orientation
    Calls check_convex_hull_atoms to see if there is truly overlap
    Finds the value of z closest to 0.
    If overlap is found:
         a list of the x,y,z coordinates of a point with overlap and the minimum z value are t returned as well as a true flag to show there is overlap
    Otherwise:
        overlap is returned as False, and coordinates are filled with dummy lists filled with -100 (which are not coordinates that would be seen otherwise)"""
    min_z = 1000 # absolute minimum | Used to check the atom of nt2 distance from nt1 after its translated to standard orientation and the same transformation is applied to nt2
    maxz = -1000 #actual maximum | checks to see if a point of an nt is on both sides of the other nt
    minz = 1000 #actual minimum | checks to see if a point of an nt is on both sides of the other nt
    retValue = [-100,-100,-100]
    # min_z shows how close two are together, where as minz shows the actual smallest z value
    inside = False
    overlap = False

    #iteratoe over list of atoms of nt2. Check xyz of each atom and see if projection is found.
    for atom in listOfAtoms:
        point = nt2.centers[atom]
        if len(point) == 3:
            x,y,z = translate_rotate_point(nt1, point) #put nt1 in standard orientation, apply same transformation to nt2, get back coordinates of atom of nt2 from nt1 center
            inside = check_convex_hull_atoms(x,y,z, parent)
            if abs(z) < abs(min_z):
                min_z = z
                retValue = [x,y,z]

            # check to see if a nt has points on both sides of the plane of a nt. See https://rna.bgsu.edu/rna3dhub/display3D/unitid/6ZMI%7C1%7CL5%7CG%7C2605,6ZMI%7C1%7CL5%7CG%7C2668 for an example.
            if z < minz:
                minz = z
            if z > maxz:
                maxz = z

            if inside:
                overlap = True # since we're iterating over the whole list of atoms, inside will be set over and over so a second flag overlap will be set that won't be reset if inside is true at least once
                               # This allows us to check for atoms that may be closer.
    if maxz > 0 and minz < 0:
        return False, [-100, -100, -100] # Don't return true for nts that have atoms on both sides of the other nts
    if overlap:
        return True, retValue
    return False, [-100,-100,-100]


def get_base_atom_names(sequence):
    """
    For standard bases, look up the base heavy and hydrogen atoms.
    For modified bases, map the parent base heavy and hydrogen atoms
    to the corresponding atoms on the modified base.
    Return a set.
    """

    if sequence in NAbaseheavyatoms:
        # standard base
        base_atoms = NAbaseatoms[sequence]

    elif sequence in modified_base_to_parent:
        # modified base
        parent_atoms = NAbaseatoms[modified_base_to_parent[sequence]]

        base_atoms = set()
        for parent_atom in parent_atoms:
            if parent_atom in parent_atom_to_modified[sequence]:
                base_atoms.add(parent_atom_to_modified[sequence][parent_atom])

    else:
        if verbose >= 2:
            print('  Not able to identify base atoms for %s' % sequence)
        base_atoms = set()

    return base_atoms


def check_base_base_stacking(nt1, nt2, parent1, parent2, datapoint):
    """
    Check for nucleotide base stacking.
    Two nucleotides and their parents are passed in.
    Create a list of their outermost atoms.
    Project nucleotides onto one another to find overlap.
    Annotated Near Stacking if the following criteria are met:
        Overlap is found at least one way
        Displacement of z coordinate is less than 4.5 and greater than 1
        Normal line z value is greater than 0.5
    Annotated as True Stacking  if the following criteria are met:
        Overlap is found both ways
        Displacement of z coord is less than 4 and greater than 1
        Normal line z value is greater than 0.6
    No annotation is generated if the criterion for near stacking or true stacking aren't met.
    """

    true_z_cutoff = 4 #angstroms, near stacking of 4.5 angstroms checked in check_convex_hull_atoms function

    interaction = ""
    interaction_reversed = ""
    reverseAnnotation = False
    #Outermost Atoms of NT Bases whose coordinates will be checked to see if they fit in the base of another nt

    # these sets are already defined
    # parent1BaseAtoms = NAbaseatoms[parent1]
    # parent2BaseAtoms = NAbaseatoms[parent2]

    #Create a list in case one of these is nucleotides is a modified nucleotide.
    #This will allow us to project atoms that may not follow the same coordinates as standard
    #nucleotides and see if they will project onto the base of another nt.

    nt1baseAtomsList = get_base_atom_names(nt1.sequence)

    if len(nt1baseAtomsList) == 0:
        if verbose >= 2:
            print("  Can't check base stacking for %s and %s" % (nt1.unit_id(),nt2.unit_id()))
        return "", datapoint, ""

    nt2baseAtomsList = get_base_atom_names(nt2.sequence)

    if len(nt2baseAtomsList) == 0:
        if verbose >= 2:
            print("  Can't check base stacking for %s and %s" % (nt1.unit_id(),nt2.unit_id()))
        return "", datapoint, ""

    #Variables to flag if an atom from nt2 was projected onto nt1 and to check if nt1 atoms project onto nt2
    nt2on1=False
    nt1on2=False

    #Is there overlap?
    #Returns true if an atom is projected inside the atom (overlap). Also returns the x,y,z coordinates of the nt inside and the minimum z value

    nt2on1, coords = return_overlap(nt2baseAtomsList, nt1, nt2, parent1)
    nt1on2, coords2 = return_overlap(nt1baseAtomsList, nt2, nt1, parent2)

    #check near stacking
    if nt2on1 or nt1on2:
        # in a near stacking instance where where nt1 projects onto nt2 but nt2 doesnt project onto nt1 its important to make sure that the normal z and min z are calculated correctly
        # and that you're using the right value.

        # Extract normal z vector and minimum z value depending on how the nts project onto one another
        # projection of nt2 onto nt 1 but not nt1 onto nt2
        if nt1on2 and not nt2on1:
            rotation_2_to_1 = np.dot(np.transpose(nt2.rotation_matrix), nt1.rotation_matrix)
            normal_Z = rotation_2_to_1[2,2]
            min_vertical_distance = coords2[2]
            reverseAnnotation = True #when projection is only found on 1 to 2 and not the other way around the interaction is switched. Use flag to mark this scenario

        # projection of nt1 onto nt2 but not nt2 onto nt1
        elif not nt1on2 and nt2on1:
            rotation_1_to_2 = np.dot(np.transpose(nt1.rotation_matrix), nt2.rotation_matrix)
            normal_Z = rotation_1_to_2[2,2]
            min_vertical_distance = coords[2]

        #projection of both
        elif nt1on2 and nt2on1:
            if coords[2] != -100 and coords2[2] != -100:
                rotation_2_to_1 = np.dot(np.transpose(nt2.rotation_matrix), nt1.rotation_matrix)
                normal_Z = rotation_2_to_1[2,2]
                min_vertical_distance = coords[2]

        if min_vertical_distance > 0:
            if normal_Z > 0:
                interaction = "ns35" # second base above, pointing up
                interaction_reversed = "ns53"
            elif normal_Z < 0:
                interaction = "ns33" #second base above, pointing down
                interaction_reversed = "ns33"

        elif min_vertical_distance < 0:
            if normal_Z > 0:
                interaction = "ns53"  #second base below, pointing up
                interaction_reversed = "ns35"
            elif normal_Z < 0:
                interaction =  "ns55" #second base below, pointing down
                interaction_reversed = "ns55"

        # if projection is only found from nt2 onto nt1, the extracted values are backwards and will lead to a backwards annotation. Swap the values
        if reverseAnnotation:
            interactionPH = interaction_reversed
            interaction_reversed = interaction
            interaction = interactionPH
        if datapoint:
            datapoint['normal_Z'] = normal_Z

        # min_distance = calculate_min_distances(nt1, nt2, None)[1]
        min_distance, heavy_min_distance, base_points, atomname = calculate_base_min_distances(nt1, nt2)

        # check for true stacking. If it meets criteria, strip the n from the annotation
        if abs(min_distance) < true_z_cutoff and abs(min_distance) > 1 and abs(normal_Z) > 0.6 and nt2on1 == True and nt1on2 == True:
            interaction = interaction.replace("n","")
            interaction_reversed = interaction_reversed.replace("n", "")

        #checks the last of the criteria to make sure its near stacking. All others get no annotation
        #Min z must be greater than 1 and the normal z should be greater than 0.5 to be considered near

        if abs(min_distance) < 1 or abs(normal_Z) < 0.5:
            return "", datapoint, ""

    if False and len(interaction) > 0:
        print('%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\t\t=hyperlink("https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s")' % (nt1.unit_id(),nt2.unit_id(),interaction,coords[0],coords[1],coords[2],nt1.unit_id(),nt2.unit_id()))

    if datapoint and len(interaction) > 0:
        datapoint['gap12'], base_points2, atomname2 = calculate_basepair_gap(nt1,nt2)
        try:
            datapoint['angle_in_plane'] = math.atan2(rotation_1_to_2[1,1],rotation_1_to_2[1,0])*57.29577951308232 - 90
        except:
            datapoint['angle_in_plane'] = math.atan2(rotation_2_to_1[1,1],rotation_2_to_1[1,0])*57.29577951308232 - 90
        datapoint['normal_Z'] = normal_Z
        datapoint['sInteraction'] = interaction
        datapoint['xStack'] = coords[0]
        datapoint['yStack'] = coords[1]
        datapoint['zStack'] = coords[2]
        datapoint['nt1on2'] = nt1on2
        datapoint['nt2on1'] = nt2on1
        datapoint['min_distance'] = min_distance
        datapoint['heavy_min_distance'] = heavy_min_distance
        datapoint['normal_Z'] = normal_Z

    return interaction, datapoint, interaction_reversed


def calculate_min_distances(nt1, nt2, base_points2):
    """
    Calculate minimum distances between two nucleotides.
    Return the minimum distance between base 1 and base 2 as base_min_distance
    Return the minimum distance between
    """

    base_min_distance = 1000
    min_distance = 1000
    base_points1 = []

    base1_atoms = get_base_atom_names(nt1.sequence)
    base2_atoms = get_base_atom_names(nt2.sequence)

    if base_points2:
        base_points2_local = [p for p in base_points2]
    else:
        base_points2_local = []

    for atom in nt1.atoms():              # nt1 atoms
        q = [atom.x, atom.y, atom.z]      # nt1 base atoms
        if atom.name in base1_atoms:
            base_points1.append(q)                 # save for later gap21 calculation

            if base_points2:
                for p in base_points2:                 # nt2 atoms
                    d = np.linalg.norm(np.subtract(p,q))
                    if d < min_distance:
                        min_distance = d
                    if d < base_min_distance:
                        base_min_distance = d
            else:
                for atom2 in nt2.atoms():
                    if atom2.name in base2_atoms:
                        p = [atom2.x, atom2.y, atom2.z]
                        base_points2_local.append(p)
                        d = np.linalg.norm(np.subtract(p,q))
                        if d < min_distance:
                            min_distance = d
                        if d < base_min_distance:
                            base_min_distance = d
        else:
            if base_points2:
                for p in base_points2:                 # nt2 atoms
                    d = np.linalg.norm(np.subtract(p,q))
                    if d < min_distance:
                        min_distance = d
            else:
                for atom2 in nt2.atoms():
                    p = [atom2.x, atom2.y, atom2.z]
                    d = np.linalg.norm(np.subtract(p,q))
                    if d < min_distance:
                        min_distance = d

    return min_distance, base_min_distance, base_points1


def calculate_base_min_distances(nt1, nt2, base_points2 = [], atomname2 = []):
    """
    Calculate minimum distances between the bases of two nucleotides
    Return the minimum distance between base 1 and base 2 as base_min_distance
    Also return the points in base 1
    This includes hydrogens, and hydrogen-hydrogen distances, which makes
    these numbers lower than you might think
    The coplanar annotation was developed with that definition of min_distance
    """

    if len(base_points2) == 0:
        base_points2 = []
        atomname2 = []
        base2_atoms = get_base_atom_names(nt2.sequence)
        for atom2 in nt2.atoms():
            if atom2.name in base2_atoms:
                p = [atom2.x, atom2.y, atom2.z]
                base_points2.append(p)
                atomname2.append(atom2.name)

    base_min_distance = 9999
    heavy_min_distance = 9999

    base1_atoms = get_base_atom_names(nt1.sequence)

    base_points1 = []
    atomname1 = []

    for atom in nt1.atoms():              # nt1 atoms
        if atom.name in base1_atoms:
            q = [atom.x, atom.y, atom.z]      # nt1 base atoms including hydrogens
            base_points1.append(q)                 # save for later gap21 calculation
            atomname1.append(atom.name)
            heavy1 = not atom.name.startswith("H")

            for i, p in enumerate(base_points2):                 # nt2 atoms
                d = np.linalg.norm(np.subtract(p,q))
                if d < base_min_distance:
                    base_min_distance = d

                if d < heavy_min_distance and heavy1 and not atomname2[i].startswith("H"):
                    heavy_min_distance = d

    return base_min_distance, heavy_min_distance, base_points1, atomname1


def look_up_atom_coordinates(nt, firstAtoms = [], secondAtoms = []):
    """
    Looping over atom names can be slow, so do that all at once here.
    """

    firstAtomCoordinates = {}
    secondAtomCoordinates = {}

    for atom in nt.atoms():
        if atom.name in firstAtoms:
            firstAtomCoordinates[atom.name] = atom
        if atom.name in secondAtoms:
            secondAtomCoordinates[atom.name] = atom

    return firstAtomCoordinates, secondAtomCoordinates

def base_backbone_modified_nucleotide_dictionary_processing(baseMassiveAndHydrogens,nt1, parent1):
    """
    Method used to add modified nucleotides to a dictionary that is used for processing in function check_base_backbone_interactions.
    This method finds atoms in a modified nucleotide that correspond with the atoms in that modified nucleotides parent.
    Checks to see if atom of modified base has the same name as its parents, if it does the parents relevent information is added to the dictionary
    for the key of the modified bases name.

    Accepts in original dictionary of backbone interactions by base, a nucleotide and its parent. Returns updated dictionary with new key value pairs for the modified nucleotide.

    NOTE: This will miss hydrogen bonds that may be on a heavy atom that isn't normally checked. This will also default to the parent case
    in cases where a methyl group is added onto a heavy atom.
     """
    baseMassiveAndHydrogens[nt1.sequence] = []
    for atoms in nt1.atoms():
        if 'P' not in atoms.name and "'" not in atoms.name: #Eliminate backbone atoms
            for atom in baseMassiveAndHydrogens[parent1]:
                if atoms.name in atom[1]:
                    baseMassiveAndHydrogens[nt1.sequence].append(atom)
    return baseMassiveAndHydrogens


def check_base_backbone_interactions(nt1,nt2,previousO3,parent1,parent2,datapoint):
    """
    Function to check base backbone interactions between the base of nt1 and backbone oxygens of nt2
    """

    # annotations to return
    phosphate = ""
    ribose = ""

    # places to store interactions meeting the requirements
    site_to_phosphate_oxygens = {}
    true_phosphate = []
    near_phosphate = []

    true_ribose = []
    near_ribose = []

    # specify cutoffs for interactions ##########################
    carbonCutoff = 4.0          # maximum massive - oxygen distance
    nCarbonCutoff = 4.5         # near

    nitrogenCutoff = 3.5        # maximum massive - oxygen distance
    nNitrogenCutoff = 4.0       # near

    angleLimit = 130            # angle limit for BPh, BR
    nAngleLimit = 110           # angle limit for near BPh, BR

    #sugarAtoms = ["C1'","C2'","O2'","C3'","O3'","C4'","O4'","C5'","O5'",'P','OP1','OP2','O3 of prev']
    # phosphate oxygens on nt2
    phosphateOxygenNames = [ "O5'", 'OP1', 'OP2']

    # 04-27-2023 For BR, O3' shouldn't be checked, it's considered phosphate, see JAR3D paper for this
    riboseOxygenNames = ["O2'","O4'"]

    # retrieve atom records, mapping to parent atoms if necessary
    phosphateOxygens = get_atom_coordinates(nt2, phosphateOxygenNames)
    riboseOxygens    = get_atom_coordinates(nt2, riboseOxygenNames)

    # use O3' coordinates of nucleotide before nt2 if available
    if previousO3.any():
        phosphateOxygens.append(previousO3)

    # if nt2 has a P atom and it is far from the plane of base 1, don't look for BPh interactions
    Pcoord = get_one_atom_coordinates(nt2, "P")
    if Pcoord.any():
        try:
            p_standard = translate_rotate_point(nt1, Pcoord)
            if not p_standard:
                phosphateOxygens = []
            elif abs(p_standard[2]) > 4.5: # phosphorus far from plane
                phosphateOxygens = []
        except:
            if verbose >= 2:
                print("  Phosphorus calculation failed for %s,%s" % (nt1.unit_id(),nt2.unit_id()))

    # Loop through each donor-hydrogen site on base 1
    for sites in NAbaseMassiveAndHydrogens[parent1]:
        baseHydrogens, baseMassive = get_atom_coordinates(nt1, sites[0:2])

        # Set the cutoff distance depending on which atom is the donor
        if "C" in sites[1]: #atoms[1] is the name of the base massive atom being checked.
            cutoff = carbonCutoff
            nCutoff = nCarbonCutoff
        elif "N" in sites[1]:
            cutoff = nitrogenCutoff
            nCutoff = nNitrogenCutoff

        # Loop through the oxygens in the phosphate backbone to extract info for base phosphate interactions
        for i, oxygen_coordinates in enumerate(phosphateOxygens):
            if oxygen_coordinates.any():
                phosphateAngle = calculate_hb_angle(baseMassive,baseHydrogens,oxygen_coordinates) # angle between base massive, its corresponding hydrogen, and oxygen
                if phosphateAngle:
                    phosphateDistance = distance_between_vectors(baseMassive,oxygen_coordinates) #distance from the oxygen to the base atom
                    if phosphateDistance:
                        if phosphateAngle > angleLimit and phosphateDistance < cutoff:
                            # a rough measure of quality of the bond
                            quality = (cutoff-phosphateDistance) + (phosphateAngle-angleLimit)/20.0
                            true_phosphate.append((-quality,sites[2],oxygen_coordinates,i))
                            if not sites[2] in site_to_phosphate_oxygens:
                                site_to_phosphate_oxygens[sites[2]] = set([i])
                            else:
                                site_to_phosphate_oxygens[sites[2]].add(i)
                        elif phosphateAngle > nAngleLimit and phosphateDistance < nCutoff:
                            # a rough measure of quality of the bond
                            quality = (nCutoff-phosphateDistance) + (phosphateAngle-nAngleLimit)/20.0
                            near_phosphate.append((-quality,"n"+sites[2],oxygen_coordinates))

        # Loop through oxygens in ribose to extract info about angle and distance
        for oxygen_coordinates in riboseOxygens:
            if oxygen_coordinates.any():
                riboseAngle = calculate_hb_angle(baseMassive,baseHydrogens,oxygen_coordinates)
                if riboseAngle:
                    riboseDistance = distance_between_vectors(baseMassive, oxygen_coordinates)
                    if riboseDistance:
                        if riboseAngle > angleLimit and riboseDistance < cutoff:
                            # a rough measure of quality of the bond
                            quality = (cutoff-riboseDistance) + (riboseAngle-angleLimit)/20.0
                            true_ribose.append((-quality,sites[3],oxygen_coordinates))
                        elif riboseAngle > nAngleLimit and riboseDistance < nCutoff:
                            # a rough measure of quality of the bond
                            quality = (nCutoff-riboseDistance) + (riboseAngle-nAngleLimit)/20.0
                            near_ribose.append((-quality,"n"+sites[3],oxygen_coordinates))

        # record the best phosphate interaction
        if len(true_phosphate) == 1:
            phosphate = true_phosphate[0][1]
            phosphate_oxygen = true_phosphate[0][2]
        elif len(site_to_phosphate_oxygens.keys()) > 1:
            # Check for multiple BPh with more than one oxygen
            if '7BPh' in site_to_phosphate_oxygens and '9BPh' in site_to_phosphate_oxygens:
                # make sure there are two different oxygen atoms; union tells if there are distinct ones
                distinct_oxygens = site_to_phosphate_oxygens['7BPh'] | site_to_phosphate_oxygens['9BPh']
                if len(distinct_oxygens) > 1:
                    phosphate = "8BPh" # C N4-1H4 and C5-H5 interact with 2 oxygens of phosphate, called 8BPh
            elif '3BPh' in site_to_phosphate_oxygens and '5BPh' in site_to_phosphate_oxygens:
                # make sure there are two different oxygen atoms; union tells if there are distinct ones
                distinct_oxygens = site_to_phosphate_oxygens['3BPh'] | site_to_phosphate_oxygens['5BPh']
                if len(distinct_oxygens) > 1:
                    phosphate = "4BPh" # G N2-2H2 and N1-H1 interacts with 2 oxygens of phosphate, called 4BPh
            if not phosphate:
                best = sorted(true_phosphate)[0]
                phosphate = best[1]
                phosphate_oxygen = best[2]
        elif len(true_phosphate) > 1:
            best = sorted(true_phosphate)[0]
            phosphate = best[1]
            phosphate_oxygen = best[2]
        elif len(near_phosphate) == 1:
            phosphate = near_phosphate[0][1]
            phosphate_oxygen = near_phosphate[0][2]
        elif len(near_phosphate) > 1:
            best = sorted(near_phosphate)[0]
            phosphate = best[1]
            phosphate_oxygen = best[2]

        # record the best ribose interaction
        if len(true_ribose) == 1:
            ribose = true_ribose[0][1]
            ribose_oxygen = true_ribose[0][2]
        elif len(true_ribose) > 1:
            best = sorted(true_ribose)[0]
            ribose = best[1]
            ribose_oxygen = best[2]
        elif len(near_ribose) == 1:
            ribose = near_ribose[0][1]
            ribose_oxygen = near_ribose[0][2]
        elif len(near_ribose) > 1:
            best = sorted(near_ribose)[0]
            ribose = best[1]
            ribose_oxygen = best[2]

    if datapoint:
        if phosphate:
            datapoint['BPh'] = phosphate
            if phosphate in ['4BPh','8BPh']:
                datapoint['BPh_oxygen'] = []
                for quality,site,oxygen_coordinates,i in true_phosphate:
                    if i in distinct_oxygens:
                        a = translate_rotate_point(nt1, oxygen_coordinates)
                        datapoint['BPh_oxygen'].append(a)
                        #print('%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\t\thttps://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s' % (nt1.unit_id(),nt2.unit_id(),phosphate,a[0],a[1],a[2],nt1.unit_id(),nt2.unit_id()))
            else:
                a = translate_rotate_point(nt1, phosphate_oxygen)
                datapoint['BPh_oxygen'] = [a]
                #print('%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\t\thttps://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s' % (nt1.unit_id(),nt2.unit_id(),phosphate,a[0],a[1],a[2],nt1.unit_id(),nt2.unit_id()))

        if ribose:
            datapoint['BR'] = ribose
            a = translate_rotate_point(nt1, ribose_oxygen)
            datapoint['BR_oxygen'] = [a]
            #print('%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\t\thttps://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s' % (nt1.unit_id(),nt2.unit_id(),ribose,a[0],a[1],a[2],nt1.unit_id(),nt2.unit_id()))

    return phosphate, ribose, datapoint


def add_atom_to_unit_id(unit_id, atom_name):
    """
    Add an atom name to a unit id.
    """

    fields = unit_id.split('|')
    if len(fields) == 5:
        u = unit_id + '|' + atom_name
    else:
        fields[5] = atom_name
        u = '|'.join(fields)

    return u


def check_oo_distance(nt1,nt2,parent1,parent2):
    """
    Find pairs of OP1 and OP2 atoms in different nucleotides that are
    close enough to accommodate ion binding
    """

    pair_list = []

    # identify some cases where nucleotides overlap, exclude those
    p1 = get_one_atom_coordinates(nt1, "P")
    p2 = get_one_atom_coordinates(nt2, "P")
    if p1.any() and p2.any():
        distance = distance_between_vectors(p1, p2)
        if distance < 3.0:
            return pair_list

    oxygens = ["OP1","OP2"]
    for oxygen1 in oxygens:
        oxygen1_coords = get_one_atom_coordinates(nt1, oxygen1)
        if oxygen1_coords.any():
            for oxygen2 in oxygens:
                oxygen2_coords = get_one_atom_coordinates(nt2, oxygen2)
                if oxygen2_coords.any():
                    distance = distance_between_vectors(oxygen1_coords, oxygen2_coords)
                    if distance < 3.5:
                        u1 = nt1.unit_id()
                        u2 = nt2.unit_id()
                        a1 = add_atom_to_unit_id(u1, oxygen1)
                        a2 = add_atom_to_unit_id(u2, oxygen2)

                        # unlike other interactions, store distance and new atom ids
                        pair_list.append((u1, u2, None, a1, a2, distance))
                        # print('OO distance %20s %20s %10.4f' % (a1, a2, distance))

    return pair_list


def check_coplanar(nt1,nt2,pair_data,datapoint):
    """
    Calculate data needed for basepair classification.
    Also, check specific criteria to say that
    the bases are enough in the same plane to be called coplanar.
    If so, return a number from 0 to 1 to measure the
    degree of coplanarity with 1 being best.
    Criteria for being coplanar or near coplanar:
    Pair.Gap must be < 97th percentile among basepairs (1.5179 Angstroms)
    min_distance must be < 97th percentile among basepairs (2.4589 A)
    Angle between center-center vector and normals must be > 70.2388 degrees
    Angle between normal vectors must be < 39.1315 degrees
    """

    pair_data["coplanar"] = False
    pair_data["coplanar_value"] = None         # 0 to 1 is coplanar, 1 is the best

    displ12 = pair_data["displ12"]

    # calculate gap and standardized atoms from nt2
    gap12, base_points2, atomname2 = calculate_basepair_gap(nt1,nt2)
    pair_data["gap12"] = gap12

    if datapoint:
        datapoint['x'] = displ12[0,0]
        datapoint['y'] = displ12[0,1]
        datapoint['z'] = displ12[0,2]
        datapoint['gap12'] = gap12

    if gap12 >= 1.5179:
        return pair_data, datapoint

    min_distance, heavy_min_distance, base_points1, atomname1 = calculate_base_min_distances(nt1, nt2, base_points2, atomname2)

    pair_data["min_distance"] = min_distance
    pair_data["heavy_min_distance"] = heavy_min_distance

    if datapoint:
        datapoint['min_distance'] = min_distance
        datapoint["heavy_min_distance"] = heavy_min_distance

    # modified nucleotides don't always have hydrogens, so be more flexible with them
    if min_distance >= 3.4589:
        return pair_data, datapoint

    # if working with regular bases, insist on close contact
    if min_distance >= 2.4589 and nt1.sequence in ['A','C','G','U'] and nt2.sequence in ['A','C','G','U']:
        return pair_data, datapoint

    center_displ = np.subtract(nt1.centers["base"],nt2.centers["base"])
    center_displ = center_displ / np.linalg.norm(center_displ) # normalize

    # calculate angle between center_displ and normal vectors to bases
    dot1 = abs(np.dot(center_displ,nt1.rotation_matrix[:,2]))[0,0]
    if dot1 >= 0.3381:
        return pair_data, datapoint

    dot2 = abs(np.dot(center_displ,nt2.rotation_matrix[:,2]))[0,0]
    if dot2 >= 0.3381:
        return pair_data, datapoint

    # calculate angle between normal vectors to the bases
    dot3 = abs(np.dot(nt1.rotation_matrix[:,2].T,nt2.rotation_matrix[:,2]))
    if dot3 <= 0.7757:
        return pair_data, datapoint

    gap21, base_points1, atomname1 = calculate_basepair_gap(nt2,nt1,base_points1)

    if datapoint:
        datapoint['gap21'] = gap21

    if gap12 <  0.5062:             # 70th percentile
      gap1val = 1
    elif gap12 <  0.9775:           # 90th percentile
      gap1val = 1+(gap12- 0.5062)*(-1.0609)
    elif gap12 <  1.5179:           # 97th percentile
      gap1val = 0.5+(gap12- 0.9775)*(-0.9252)
    else:
      gap1val = 0

    if gap21 <  0.5062:             # 70th percentile
      gap2val = 1
    elif gap21 <  0.9775:           # 90th percentile
      gap2val = 1+(gap21- 0.5062)*(-1.0609)
    elif gap21 <  1.5179:           # 97th percentile
      gap2val = 0.5+(gap21- 0.9775)*(-0.9252)
    else:
      gap2val = 0

    if dot1 <  0.1139:              # 70th percentile
      dot1val = 1
    elif dot1 <  0.2193:            # 90th percentile
      dot1val = 1+(dot1- 0.1139)*(-4.7408)
    elif dot1 <  0.3381:            # 97th percentile
      dot1val = 0.5+(dot1- 0.2193)*(-4.2103)
    else:
      dot1val = 0

    if dot2 <  0.1139:              # 70th percentile
      dot2val = 1
    elif dot2 <  0.2193:            # 90th percentile
      dot2val = 1+(dot2- 0.1139)*(-4.7408)
    elif dot2 <  0.3381:            # 97th percentile
      dot2val = 0.5+(dot2- 0.2193)*(-4.2103)
    else:
      dot2val = 0

    if -dot3 < -0.9509:             # 70th percentile
      dot3val = 1
    elif -dot3 < -0.8835:           # 90th percentile
      dot3val = 1+(-dot3-(-0.9509))*(-7.4217)
    elif -dot3 < -0.7757:           # 97th percentile
      dot3val = 0.5+(-dot3-(-0.8835))*(-4.6390)
    else:
      dot3val = 0

    if min_distance <  1.8982:      # 70th percentile
      min_dist_val = 1
    elif min_distance <  2.1357:    # 90th percentile
      min_dist_val = 1+(min_distance- 1.8982)*(-2.1050)
    elif min_distance <  2.4859:    # 97th percentile
      min_dist_val = 0.5+(min_distance- 2.1357)*(-1.4280)
    else:
      min_dist_val = 0

    # Pair.Coplanar is 1 if all are within the 70th percentile
    # Pair.Coplanar is 0.5 if all are within the 90th percentile
    # Pair.Coplanar is > 0 if all are within the 97th percentile
    # Between these, it decreases linearly

    pair_data["coplanar"] = True
    pair_data["coplanar_value"] = min([gap1val, gap2val, dot1val, dot2val, dot3val, min_dist_val])

    if datapoint:
        datapoint['coplanar'] = pair_data['coplanar']
        datapoint['coplanar_value'] = pair_data['coplanar_value']

    return pair_data, datapoint


def calculate_basepair_gap(nt1,nt2,base_points2=[],atomname=[]):
    """
    Calculate the vertical distance between nearest edges of two bases,
    from the plane of nt1 to the nearest atom of nt2.
    """

    displacements = []
    distances = []

    if len(base_points2) > 0:
        # calculate distances from base atoms of nt2 to center of nt1 base
        for p in base_points2:
            v = np.subtract(p,nt1.centers["base"])
            d = np.linalg.norm(v)
            displacements.append(v)
            distances.append(d)

    else:
        # look up the base atoms in nt2
        base_points2 = []
        atomname = []

        base_atoms = get_base_atom_names(nt2.sequence)

        for atom in nt2.atoms():
            if atom.name in base_atoms:
                p = [atom.x, atom.y, atom.z]
                base_points2.append(p)
                atomname.append(atom.name)

                v = np.subtract(p,nt1.centers["base"])
                d = np.linalg.norm(v)
                displacements.append(v)
                distances.append(d)

    # sort indices of atoms in nt2 by distance to center of nt1
    indices = np.argsort(distances)

    m = min(3,len(indices))

    gap12 = 100
    for k in range(0,m):              # 3 nearest points
        p = displacements[indices[k]]
        z = abs(np.dot(p,nt1.rotation_matrix[:,2])[0,0])  # distance out of plane of nt1
        if z < gap12:
            gap12 = z                 # gap is smallest z value

    return gap12, base_points2, atomname


def check_sugar_ribose(nt1,nt2,parent1,datapoint):
    """
    Check for O2'-O2' distance being compatible with a hydrogen bond.
    When nt1 is C or U, check O2-O2' distance and O2' being near the plane of base 1
    When nt2 is A or G, check N3-O2' distance and O2' being near the plane of base 1
    """

    nt1_o2p = get_one_atom_coordinates(nt1,"O2'")

    if not len(nt1_o2p) == 3:
        # nt1 has no identified O2' atom
        return "", datapoint

    nt2_o2p = get_one_atom_coordinates(nt2,"O2'")

    if not len(nt2_o2p) == 3:
        # nt2 has no identified O2' atom
        return "", datapoint

    o2p_o2p_displ = np.subtract(nt1_o2p,nt2_o2p)
    o2p_o2p_distance = np.linalg.norm(o2p_o2p_displ)

    if o2p_o2p_distance > 3.8:
        # O2' atoms are too far apart
        return "", datapoint

    if datapoint:
        datapoint['o2p_o2p_distance'] = o2p_o2p_distance

    if parent1 in ['A', 'G']:
        base_point = get_one_atom_coordinates(nt1,'N3')
    elif parent1 in ['C','U']:
        base_point = get_one_atom_coordinates(nt1,'O2')
    else:
        base_point = None

    if not len(base_point) == 3:
        # nt1 does not have the appropriate base atom
        return "", datapoint

    base_o2p_displ = np.subtract(base_point,nt2_o2p)
    base_o2p_distance = np.linalg.norm(base_o2p_displ)

    if base_o2p_distance > 3.8:
        # nt1 base atom is too far from nt2 O2' atom
        return "", datapoint

    if datapoint:
        datapoint['base_o2p_distance'] = base_o2p_distance

    z = abs(np.dot(base_o2p_displ,nt1.rotation_matrix[:,2])[0,0])  # distance of nt2 O2' out of plane of nt1

    if z > 2.0:
        # nt2 O2' atom is too far above or below the plane of base of nt1
        return "", datapoint

    if datapoint:
        datapoint['nt2_o2p_height'] = z

    p0, p1 = get_atom_coordinates(nt1,["C1'","C2'"])
    p2, p3 = get_atom_coordinates(nt2,["C2'","C1'"])

    if len(p0) == 3 and len(p1) == 3 and len(p2) == 3 and len(p3) == 3:
        angle = torsion_angle(p0,p1,p2,p3)
        if abs(angle) > 90:
            # cis case, as in cSS
            annotation = 'cSR'
        else:
            # trans case, base flipped over compared to cSS
            annotation = 'tSR'
    else:
        if verbose >= 2:
            print('  Not able to calculate orientation of SR annotation for %s-%s' % (nt1.unit_id(),nt2.unit_id()))
        return "", datapoint

    if datapoint:
        datapoint['sugar_ribose'] = annotation

    # if annotation in ['cSR','tSR']:
    #     print('%s %8.4f %8.4f %8.4f %8.4f %s-%s https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s https://rna.bgsu.edu/correspondence/variability?id=%s,%s&format=unique' % (annotation,o2p_o2p_distance,base_o2p_distance,z,angle,nt1.sequence,nt2.sequence,nt1.unit_id(),nt2.unit_id(),nt1.unit_id(),nt2.unit_id()))

    return annotation, datapoint


def store_basepair_quality(cutoff_distance,pair_data,LW,hydrogen_bonds):
    quality = {}
    quality['cutoff_distance'] = cutoff_distance
    quality['max_gap'] = max(pair_data["gap12"],pair_data["gap21"])
    quality['atoms1'] = []
    quality['atoms2'] = []
    LW_clean = LW.replace("n","")
    if LW_clean in hydrogen_bonds:
        for donor, hydrogen, acceptor, direction,a,b,c,d in hydrogen_bonds[LW_clean]:
            if direction == "12":
                quality['atoms1'].append(hydrogen)
                quality['atoms2'].append(acceptor)
            else:
                quality['atoms2'].append(hydrogen)
                quality['atoms1'].append(acceptor)

    return quality


def check_basepair_cutoffs(nt1,nt2,pair_data,cutoffs,hydrogen_bonds,datapoint):
    """
    Given nt1 and nt2 and the dictionary of cutoffs for that pair of nucleotides,
    check cutoffs for each basepair interaction type.
    Also compute hydrogen bonds for all basepairs consistent with the normal vector.
    If no cutoffs are fully met but all hydrogen bonds are met, annotate as near.
    datapoint accumulates information about the interaction for diagnostics.
    """

    quality = {}                  # dictionary of quality scores for each interaction type

    displ = pair_data["displ12"]  # vector from origin to nt2 when standardized

    if abs(displ[0,2]) > 3.6:     # too far out of plane for a basepair; don't check further
        return "", "", quality, datapoint

    # check sign of normal vector to cut number of possible families in half
    rotation_1_to_2 = np.dot(np.transpose(nt1.rotation_matrix), nt2.rotation_matrix)
    normal_Z = rotation_1_to_2[2,2]   # z component of normal vector to second base
    if normal_Z > 0:
        normal_sgn = 1
        possible_interactions = list(cutoffs[1].keys())
    else:
        normal_sgn = -1
        possible_interactions = list(cutoffs[-1].keys())

    if datapoint:
        datapoint['normal_Z'] = normal_Z

    if datapoint:
        # calculate and store hydrogen bond data even if not needed
        # in order to save diagnostic information
        check_order = ['hydrogen bonds','cutoffs']
    else:
        # check cutoffs first and exit if no interaction is present
        check_order = ['cutoffs','hydrogen bonds']

    LW_bond_rank = []

    # check possible basepairs two different ways
    for check in check_order:
        if check == 'cutoffs':
            # check cutoffs for each interaction type
            if datapoint:
                cutoff_distance_max = 5.0    # keep checking up to this number to have the data
                cutoff_distance_max = near_discrepancy_cutoff    # keep checking up to this number to have the data
            else:
                cutoff_distance_max = near_discrepancy_cutoff # faster annotation

            ok_normal_displ = []   # interactions with OK normal and displacement
            for interaction in possible_interactions:
                for subcategory in cutoffs[normal_sgn][interaction].keys():
                    cut = cutoffs[normal_sgn][interaction][subcategory]
                    cutoff_distance = 0.0
                    cutoff_distance += max(0,cut['xmin'] - displ[0,0])  # how far below xmin
                    cutoff_distance += max(0,displ[0,0] - cut['xmax'])  # how far above xmax

                    if cutoff_distance >= cutoff_distance_max:
                        continue

                    cutoff_distance += max(0,cut['ymin'] - displ[0,1])  # how far below ymin
                    cutoff_distance += max(0,displ[0,1] - cut['ymax'])  # how far above ymax

                    if cutoff_distance >= cutoff_distance_max:
                        continue

                    if 'radiusmax' in cut:
                        radius = math.sqrt(displ[0,0]**2 + displ[0,1]**2)
                        cutoff_distance += max(0,radius - cut['radiusmax'])  # how far above radiusmax

                        if 'radiusmin' in cut:
                            cutoff_distance += max(0,cut['radiusmin'] - radius)  # how far below radiusmin

                    if cutoff_distance >= cutoff_distance_max:
                        continue

                    cutoff_distance += max(0,cut['zmin'] - displ[0,2])  # how far below zmin
                    cutoff_distance += max(0,displ[0,2] - cut['zmax'])  # how far above zmax

                    if cutoff_distance >= cutoff_distance_max:
                        continue

                    # accentuate wrong normal vector by factor of 3
                    cutoff_distance += 3*max(0,cut['normalmin'] - normal_Z)  # how far below normalmin
                    cutoff_distance += 3*max(0,normal_Z - cut['normalmax'])  # how far above normalmax

                    # if abs(normal_Z) < 0.4:
                    #     # bases are too close to being perpendicular
                    #     # but in some categories, these are important, so comment this out
                    #     cutoff_distance += near_discrepancy_cutoff

                    # cutoffs are met or close enough for now
                    if cutoff_distance < cutoff_distance_max:
                        ok_normal_displ.append((interaction,subcategory,cut,cutoff_distance)) # ("cWW",0), etc.

            # if not close to meeting any cutoffs and we are not collecting data, return now to save time
            if len(ok_normal_displ) == 0 and not datapoint:
                return "", "", quality, datapoint

            # calculation revised to have the right sense to it 2023-07-19 CLZ
            # it was OK for cWW and other families where you see 3 and 5 faces
            # but for cHS and others where both 3' faces point the same direction, the angle is mostly reversed now,
            # but more than just a sign change the more the bases are tilted relative to each other
            angle_in_plane = math.atan2(rotation_1_to_2[1,1],rotation_1_to_2[0,1])*57.29577951308232 - 90

            if angle_in_plane <= -90:
                angle_in_plane += 360

            ok_angle_in_plane = []

            # 3 Angstrom radius rotated by 10 degrees moves 3*10*pi/180 = 0.523 Angstroms
            # Divide angle by 20 to get somewhat equivalent distance in Angstroms
            # but then penalize more because angle in plane is more disruptive than displacement
            for interaction,subcategory,cut,cutoff_distance in ok_normal_displ:
                angle_penalty = 0.1
                if cut['anglemin'] < cut['anglemax']:     # for ranges within -90 to 270 like 50 to 120
                    cutoff_distance += angle_penalty*max(0,cut['anglemin'] - angle_in_plane)  # how far below anglemin
                    cutoff_distance += angle_penalty*max(0,angle_in_plane - cut['anglemax'])  # how far above anglemax
                else:                                     # for ranges straddling 270 like 260 to -75
                    cutoff_distance += angle_penalty*min(max(0,cut["anglemin"]-angle_in_plane),max(0,angle_in_plane-cut["anglemax"]))

                if cutoff_distance < cutoff_distance_max:
                    ok_angle_in_plane.append((interaction,subcategory,cut,cutoff_distance))

            # if not close to meeting any cutoffs and we are not collecting data, return now
            if len(ok_angle_in_plane) == 0 and not datapoint:
                return "", "", quality, datapoint

            if not 'gap12' in pair_data:
                # calculate gap and standardized atoms from nt2
                gap12, base_points2, atomname2 = calculate_basepair_gap(nt1,nt2)
                pair_data["gap12"] = gap12

            if not 'gap21' in pair_data:
                # calculate gap and standardized atoms from nt2
                gap21, base_points1, atomname1 = calculate_basepair_gap(nt2,nt1)
                pair_data["gap21"] = gap21

            if not 'heavy_min_distance' in pair_data:
                min_distance, heavy_min_distance, base_points1, atomname1 = calculate_base_min_distances(nt1, nt2)
                pair_data["min_distance"] = min_distance
                pair_data["heavy_min_distance"] = heavy_min_distance

            if datapoint:
                datapoint['angle_in_plane'] = angle_in_plane
                datapoint['gap12'] = pair_data["gap12"]
                datapoint['gap21'] = pair_data["gap21"]
                datapoint['gapmax'] = max(pair_data["gap21"],pair_data["gap12"])
                datapoint['min_distance'] = pair_data['min_distance']
                datapoint['heavy_min_distance'] = pair_data['heavy_min_distance']

            match = []              # Meets the cutoffs for a category like cWW
            direct_near_match = []  # Like ncWW category
            near_match = []         # Close to a category like cWW
            check_hbonds = []       # Which LW families to check hydrogen bonds for
            for interaction,subcategory,cut,cutoff_distance in ok_angle_in_plane:
                if pair_data["min_distance"] < 0.5:
                    # unrealistically close to one another, cannot be a basepair
                    cutoff_distance += near_discrepancy_cutoff

                if pair_data["heavy_min_distance"] < 1.5:
                    # unrealistically close to one another, cannot be a basepair
                    # distance could be set more carefully
                    cutoff_distance += near_discrepancy_cutoff

                if cut['gapmax'] > 0.1:
                    gap_diff = max(0,max(pair_data["gap12"],pair_data["gap21"])-cut['gapmax'])
                    # accentuate wrong gap
                    cutoff_distance += 3*gap_diff

                    # extra penalty when both are out of plane with each other
                    if min(pair_data["gap12"],pair_data["gap21"]) > cut['gapmax']:
                        cutoff_distance += 3*gap_diff

                # identify cases where there is no base-base hydrogen bond to be examined
                cSS_one_hbond = False
                if interaction == 'cSs' and pair_data['parent2'] in ['C','U']:
                    cSS_one_hbond = True
                elif interaction == 'csS' and pair_data['parent2'] == 'U':
                    cSS_one_hbond = True

                if cutoff_distance > 0:
                    # impose the near discrepancy cutoff now, must be near a true category, not near a near category
                    if cutoff_distance < near_discrepancy_cutoff and not interaction.startswith("n"):
                        # must have minimum distance between bases to be counted as near
                        if pair_data['heavy_min_distance'] < near_heavy_distance_cutoff or cSS_one_hbond:
                            # bases are close enough to list as a near pair
                            near_match.append((interaction,subcategory,cutoff_distance))
                            check_hbonds.append(interaction)
                elif interaction.startswith("n"):
                    # directly classified as near, like for certain single h-bonds
                    # trust it and don't check hydrogen bonds
                    direct_near_match.append([interaction,subcategory,cutoff_distance])
                    if verbose >= 2:
                        print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  directly classified as near" % (interaction,nt1.unit_id(),nt2.unit_id(),nt1.unit_id(),nt2.unit_id()))
                elif pair_data['heavy_min_distance'] > true_heavy_distance_cutoff and not cSS_one_hbond:
                    # matches a true category but the bases are too far apart for a good basepair
                    near_match.append([interaction,subcategory,cutoff_distance])
                    check_hbonds.append(interaction)
                    if verbose >= 2:
                        print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  switched to near due to minimum distance" % (interaction,nt1.unit_id(),nt2.unit_id(),nt1.unit_id(),nt2.unit_id()))
                else:
                    # true pair; could conceivably match more than one category
                    match.append([interaction,subcategory,cutoff_distance])
                    check_hbonds.append(interaction)

            if not datapoint:
                # having just checked cutoffs, we want to check hydrogen bonds
                # only for the interactions in match and near_match
                # Restrict possible interactions to those that meet all cutoffs
                possible_interactions = check_hbonds

        else:
            # check hydrogen bonds
            #print("Checking hydrogen bonds for https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt1.unit_id(),nt2.unit_id()))

            atom_set_to_bond_parameters = {}
            LW_to_atom_sets = {}
            #donor_hydrogen_to_badness = {}

            # check hydrogen bonds for interactions that are possible by the normal vector or cutoffs
            for LW in possible_interactions:
                if LW in hydrogen_bonds:
                    for atom_set in hydrogen_bonds[LW]:

                        if not LW in LW_to_atom_sets:
                            LW_to_atom_sets[LW] = set()
                        LW_to_atom_sets[LW].add(atom_set)

                        if atom_set in atom_set_to_bond_parameters:
                            # this atom_set was already checked for a different LW family
                            result = atom_set_to_bond_parameters[atom_set]

                        else:
                            # check atom set for hydrogen bond

                            # atom sets are listed from donor to acceptor
                            if atom_set[3] == '12':
                                result = check_hydrogen_bond(nt1,nt2,atom_set)
                            else:
                                result = check_hydrogen_bond(nt2,nt1,atom_set)

                            # result is a dictionary with many fields
                            atom_set_to_bond_parameters[atom_set] = result

            # store all information about hydrogen bonds checked
            if datapoint:
                datapoint["LW_to_atom_sets"] = LW_to_atom_sets
                datapoint["atom_set_to_results"] = atom_set_to_bond_parameters

            # count hydrogen bonds for each possible annotation
            LW_bond_counter = []
            # store interactions with enough hydrogen bonds
            hbond_interactions = set()
            # store mapping from LW family to second-shortest bond length
            LW_to_distances = {}
            for LW in possible_interactions:
                LW_to_distances[LW] = []
                if LW in hydrogen_bonds:
                    checked_counter = 0
                    bond_counter = 0
                    badnesses = []
                    for atom_set in hydrogen_bonds[LW]:
                        result = atom_set_to_bond_parameters[atom_set]

                        if result["bond_checked"]:
                            checked_counter += 1
                            # donor     = atom_set[0]
                            # hydrogen  = atom_set[1]
                            # direction = atom_set[3]
                            badnesses.append(result["badness"])

                            # if the badness is not so far from the best that we have seen for this donor
                            # move away from comparing quality of hydrogen bonds, since that works
                            # differently if you check for lots of pairing families versus just one
                            # Results differ depending on order of checking cutoffs versus h-bonds
                            # if True or result["badness"] < donor_hydrogen_to_badness[(donor,hydrogen,direction)] + 0.5:
                                # evaluate the bond in the context of the specific LW family

                            mind = atom_set[4]
                            maxd = atom_set[5]
                            distance = result["donor_acceptor_distance"]

                            LW_to_distances[LW].append(distance)

                            if mind <= distance and distance <= maxd:
                                mina = atom_set[6]
                                maxa = atom_set[7]
                                hb_angle = result["heavy_donor_acceptor_angle"]
                                if mina <= hb_angle and hb_angle <= maxa:
                                    bond_counter += 1

                    # this is some older diagnostics
                    if checked_counter > 0:
                        if len(badnesses) == 1:
                            # single h-bond interactions are compared to each other only
                            second_badness = badnesses[0] + 100.0
                        else:
                            # two or more h-bond interactions are preferred over one
                            second_badness = sorted(badnesses)[1]
                        LW_bond_counter.append((LW,bond_counter,checked_counter,second_badness))

                    # this is where the decision is made
                    if checked_counter == 1 and bond_counter == 1:
                        hbond_interactions.add(LW)
                    elif checked_counter >= 2 and bond_counter >= 2:
                        hbond_interactions.add(LW)
                    elif checked_counter == 0:
                        # can't reject based on hydrogen bonds, if there are none
                        hbond_interactions.add(LW)
                else:
                    # can't reject based on hydrogen bonds, if there are none
                    hbond_interactions.add(LW)

            # best annotation considering only hydrogen bonds
            if len(LW_bond_counter) > 0:
                # sort interactions to find the best hydrogen bonds
                # sort by badness of second worst hydrogen bond
                LW_bond_rank = sorted(LW_bond_counter, key=lambda x : (x[3]))
                LW = LW_bond_rank[0][0]   # best LW category

                if LW_bond_rank[0][3] < 2.0:
                    # if second_badness is not horrible
                    if datapoint:
                        datapoint['hbond_best_pair'] = LW
                    #datapoint['hbond'] = LW_bonds[LW]
                    #datapoint['hbond_messages'] = LW_bond_messages[LW]

    # use hydrogen bond data to update the interactions in match, maybe change from true to near
    still_match = []
    demotion = False
    for m in range(0,len(match)):
        # if not already near, if not enough hydrogen bonds, and if not a basepair subcategory
        if not match[m][0].startswith("n") and not match[m][0] in hbond_interactions:
            if verbose >= 2:
                print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  demoted to near by hydrogen bonds" % (match[m][0],nt1.unit_id(),nt2.unit_id(),nt1.unit_id(),nt2.unit_id()))

            # move to the near match list, don't keep in the match list
            near_match.append(match[m])
            demotion = True
        else:
            still_match.append(match[m])
    match = still_match


    if datapoint and demotion and len(match) == 0:
        # cannot easily and accurately track what the demotion was from, but note it anyway
        datapoint['demoted_hbond'] = True

    if len(match) > 0:
        # at least one perfect match, omit the near matches
        match = sorted(match, key=lambda x: (x[1],len(x[0]))) # sort by subcategory, then interaction name length
    elif len(direct_near_match) > 0:
        match = [direct_near_match[0]]     # use the first one
    elif len(near_match) > 0:
        # check hydrogen bond lengths, then
        # sort near matches by cutoff_distance

        near_matches = []
        for LW, subcategory, cutoff_distance in near_match:
            # check hydrogen bond lengths
            dist2 = 10.0
            if LW in LW_to_distances:
                if len(LW_to_distances[LW]) > 1:
                    # use second shortest
                    dist2 = sorted(LW_to_distances[LW])[1]
                elif len(LW_to_distances[LW]) == 1:
                    # only one hydrogen bond, use that
                    dist2 = LW_to_distances[LW][0]
            if dist2 < 5.0:
                # hydrogen bond length is short enough to be considered a near pair
                # mark the interaction as near
                # be ready to rank both by cutoff_distance and dist2
                near_matches.append(("n"+LW,subcategory,cutoff_distance,dist2))
            else:
                if verbose >= 2:
                    print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  rejected h-bonds, length %8.2f" % (LW,nt1.unit_id(),nt2.unit_id(),nt1.unit_id(),nt2.unit_id(),dist2))
                pass

        if len(near_matches) > 0:
            # sort by product of cutoff distance and second shortest hydrogen bond
            near_matches = sorted(near_matches, key=lambda x: x[2]*x[3])
            match = [near_matches[0][0:3]]     # use the nearest one, call it near
        else:
            # no matches at all, return what we have so far
            return "", "", quality, datapoint
    else:
        # no matches at all, return what we have so far
        return "", "", quality, datapoint

    if len(match) == 1:
        interaction,subcategory,cutoff_distance = match[0]
        if cutoff_distance > 0 and not "n" in interaction:
            LW = "n" + interaction
        else:
            LW = interaction

        quality = store_basepair_quality(cutoff_distance,pair_data,LW,hydrogen_bonds)

        if datapoint:
            datapoint['basepair'] = LW
            datapoint['basepair_subcategory'] = match[0][1]
            datapoint['cut_dist'] = match[0][2]
            #datapoint['hbond'] = LW_bonds[interaction]
            #datapoint['hbond_messages'] = LW_bond_messages[interaction]
        return LW, subcategory, quality, datapoint
    else:
        # multiple matching basepair interactions between these two nucleotides
        LW_remaining = set([i.replace("n","") for i,s,cd in match])
        if len(LW_remaining) == 1:
            # one family, mutiple subcategories, quite OK, they are designed to overlap
            interaction,subcategory,cutoff_distance = match[0]
            if cutoff_distance > 0 and not "n" in interaction:
                LW = "n" + interaction
            else:
                LW = interaction

            quality = store_basepair_quality(cutoff_distance,pair_data,LW,hydrogen_bonds)

            if datapoint:
                datapoint['basepair'] = LW
                datapoint['basepair_subcategory'] = subcategory
                datapoint['cut_dist'] = cutoff_distance
                #datapoint['hbond'] = LW_bonds[LW]
                #datapoint['hbond_messages'] = LW_bond_messages[LW]
            return LW, subcategory, quality, datapoint

        else:
            # loop over hydrogen bond sets from best to worst
            for LW,bond_counter,checked_counter,max_badness in LW_bond_rank:
                for LW2,subcategory,cutoff_distance in match:
                    if LW == LW2:
                        quality = store_basepair_quality(cutoff_distance,pair_data,LW,hydrogen_bonds)
                        if datapoint:
                            datapoint['basepair'] = LW
                            datapoint['basepair_subcategory'] = subcategory
                            datapoint['cut_dist'] = cutoff_distance
                            #datapoint['hbond'] = LW_bonds[LW]
                            #datapoint['hbond_messages'] = LW_bond_messages[LW]
                        if verbose >= 2:
                            print("  %5s %-22s %-22s  https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  multiple annotations meet cutoffs: %s" % (LW,nt1.unit_id(),nt2.unit_id(),nt1.unit_id(),nt2.unit_id(),LW_remaining))
                        return LW, subcategory, quality, datapoint

            if verbose >= 2:
                print("  No match between all cutoffs and all hydrogen bonds, using first match")
            interaction,subcategory,cutoff_distance = match[0]
            if cutoff_distance > 0 and not "n" in interaction:
                LW = "n" + interaction
            else:
                LW = interaction

            quality = store_basepair_quality(cutoff_distance,pair_data,LW,hydrogen_bonds)

            if datapoint:
                datapoint['basepair'] = LW
                datapoint['basepair_subcategory'] = subcategory
                datapoint['cut_dist'] = cutoff_distance
                #datapoint['hbond'] = LW_bonds[LW]
                #datapoint['hbond_messages'] = LW_bond_messages[LW]
            return LW, subcategory, quality, datapoint


def get_glycosidic_atom_coordinates(nt,parent):

    gly = None

    if nt.sequence in ['A','G','DA','DG']:
        gly = nt.centers["N9"]
    elif nt.sequence in ['C','U','DC','DT']:
        gly = nt.centers["N1"]
    elif nt.sequence in parent_atom_to_modified:
        if parent in ['A','G','DA','DG']:
            gly = nt.centers[parent_atom_to_modified[nt.sequence]["N9"]]
        elif parent in ['C','U','DC','DT']:
            gly = nt.centers[parent_atom_to_modified[nt.sequence]["N1"]]

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

    angle = math.acos(np.dot(np.dot(b.T,rotation),b)[0,0] / np.dot(b.T,b)[0,0]).real
    angle = angle * np.sign(np.linalg.det(np.concatenate((b,np.dot(rotation,b),axis),axis=1)))
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


def angle_between_vectors(vec1, vec2):
    # Calculate an angle from 0 to 180 degrees between two vectors
    if len(vec1) == 3 and len(vec2) == 3:
        cosang = np.dot(vec1, vec2)
        sinang = np.linalg.norm(np.cross(vec1, vec2))
        angle = np.arctan2(sinang, cosang)
        return 180*angle/np.pi
    else:
        return None


def angle_between_three_points(P1,P2,P3):
    # Calculate an angle from 0 to 180 degrees between vector P2-P1 and P2-P3
    if len(P1) == 3 and len(P2) == 3 and len(P3) == 3:
        return angle_between_vectors(P1-P2,P3-P2)
    else:
        return None


def distance_between_vectors(vec1, vec2):
    # Calculate the distance between two vectors
    if len(vec1) == 3 and len(vec2) == 3:
        return np.linalg.norm(np.subtract(vec1,vec2))
    else:
        return None


def unit_vector(v):
    return v / np.linalg.norm(v)


def torsion_angle(p0,p1,p2,p3):
    """
    From https://stackoverflow.com/questions/20305272/dihedral-torsion-angle-from-four-points-in-cartesian-coordinates-in-python
    Pass in four vectors.
    """

    b0 = -1.0*(p1 - p0)
    b1 = p2 - p1
    b2 = p3 - p2

    # normalize b1 so that it does not influence magnitude of vector
    # rejections that come next
    b1 /= np.linalg.norm(b1)

    # vector rejections
    # v = projection of b0 onto plane perpendicular to b1
    #   = b0 minus component that aligns with b1
    # w = projection of b2 onto plane perpendicular to b1
    #   = b2 minus component that aligns with b1
    v = b0 - np.dot(b0, b1)*b1
    w = b2 - np.dot(b2, b1)*b1

    # angle between v and w in a plane is the torsion angle
    # v and w may not be normalized but that's fine since tan is y/x
    x = np.dot(v, w)
    y = np.dot(np.cross(b1, v), w)
    return np.degrees(np.arctan2(y, x))

def map_PDB_list_to_PDB_IFE_dict(PDB_list):
    """
    map a list of PDB ids or IFEs or URLs to a dictionary whose keys
    are PDB ids and whose values are representative chains in that PDB.
    If one PDB has a lot of representative chains, several chains will be joined with +.
    """

    PDB_IFE_Dict = defaultdict(str)   # accumulate PDB-IFE pairs
    for PDB in PDB_list:
        try:
            if "nrlist" in PDB and "NR_" in PDB:
                                          # referring to an equivalence class online
                                          # download the entire representative set,
                                          # then find the right line for the equivalence class
                                          # then extract the list
                if sys.version_info[0] < 3:
                    f = urllib.urlopen(PDB)
                    myfile = f.read()
                else:
                    f = urllib.request.urlopen(PDB)
                    myfile = f.read().decode()

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
                if sys.version_info[0] < 3:
                    f = urllib.urlopen(PDB)
                    myfile = f.read()
                else:
                    f = urllib.request.urlopen(PDB)
                    myfile = f.read().decode()
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
        except:
            if verbose >= 0:
                print("  Not able to map %s to its IFEs" % PDB)

    # remove leading + signs
    for PDB in PDB_IFE_Dict:
        if PDB_IFE_Dict[PDB].startswith("+"):
            PDB_IFE_Dict[PDB] = PDB_IFE_Dict[PDB][1:]

    return PDB_IFE_Dict


def write_unit_data_file(PDB,unit_data_path,structure):
    """
    Write out data file(s) of nucleotide centers and rotation matrices,
    primarily for use by the FR3D motif search tool.
    If unit_data_path is empty, no files are written.
    One file for each chain.

    This function may not write out all modified nucleotides
    It will probably miss solitary nucleotides like ATP.
    """

    if len(unit_data_path) > 0:

        nucleotides = structure.residues(type = ["RNA","DNA","PNA"])
        all_nts = {}

        # get the nucleotides in each model and chain, able to sort by symmetry and index
        for nt in nucleotides:
            fields = nt.unit_id().split("|")
            id = "_".join(fields[0:3])   # PDB_model_chain

            if len(fields) == 9:
                symmetry = fields[8]
            else:
                symmetry = ""

            if not id in all_nts:
                all_nts[id] = []

            if nt.index:
                all_nts[id].append((symmetry,nt.index,nt))
            else:
                all_nts[id].append((symmetry,0,nt))

        # loop over models and chains
        for id in all_nts.keys():

            # write out data for each nucleotide
            # note that _NA goes with glycosidic centers, while _RNA would be for base centers
            filename = os.path.join(unit_data_path, "units", id + "_NA.pickle")

            units = []
            order = []
            cntrs = []
            rttns = []

            # sort by symmetry and index
            for symmetry,index,nt in sorted(all_nts[id], key=lambda p: (p[0],p[1])):
                units.append(nt.unit_id())
                order.append(nt.index)
                cntrs.append(nt.centers["glycosidic"])
                rttns.append(nt.rotation_matrix)

            rsset = [units, order, cntrs, rttns]

            with open(filename, 'wb') as fh:
                # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
                pickle.dump(rsset, fh, 2)

            if verbose >= 1:
                print("  Wrote unit data file %s" % filename)


def write_txt_output_file(outputNAPairwiseInteractions,file_id,interaction_to_list_of_tuples,categories,category_to_interactions):
    """
    Write interactions according to category, and within each
    category, write by annotation.
    """

    # loop over types of output files requested
    for category in categories:
        if category in ["near","lower","alternative","cwb","loops"]:
            continue

        filename = os.path.join(outputNAPairwiseInteractions,file_id + "_" + category + ".txt")

        quads_to_write = []
        # loop over all interactions found in this category
        for interaction in sorted(category_to_interactions[category]):
            if category == 'basepair' and not "near" in categories and "n" in interaction:
                continue

            if category == 'basepair' and not 'cwb' in categories and ('cWB' in interaction or 'cBW' in interaction):
                continue

            inter = interaction

            if category == 'basepair' and not "lower" in categories:
                # capitalize base edges to simplify
                inter = interaction.replace("w","W").replace("s","S").replace("h","H")

            if category == 'basepair' and not "alternative" in categories:
                # remove "alternative" designations
                inter = inter.replace("a","")

            # if this category has a restricted list of interactions to output
            if len(categories[category]) == 0 or inter in categories[category]:
                for itpl in interaction_to_list_of_tuples[interaction]:
                    if len(itpl) == 3:
                        a,b,c = itpl
                        quads_to_write.append((a,inter,b,c))
                    else:
                        if interaction == 'oo_distance':
                            # also write a distance and a url to view the interaction
                            # print(itpl)
                            u1, u2, crossing, a1, a2, distance = itpl
                            url = "https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (u1,u2)
                            quads_to_write.append((a1,inter,a2,crossing,"%0.4f" % distance,url))

        # sort quads by model, first chain, first number, first unit id (for alt id, insertion code, symmetry), interaction
        ordered = sorted(quads_to_write, key=lambda x: (int(x[0].split("|")[1]) or 0,x[0].split("|")[2],int(x[0].split("|")[4]),x[0],x[1],x[2]))
        with open(filename,'w') as f:
            for o in ordered:
                if len(o) == 4:
                    f.write("%s\t%s\t%s\t%s\n" % (o))
                elif len(o) == 6:
                    f.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (o))

    if 'loops' in interaction_to_list_of_tuples:
        # follow format used by https://rna.bgsu.edu/rna3dhub/loops/download_with_breaks/8GLP
        filename = os.path.join(outputNAPairwiseInteractions,file_id + "_loops.txt")
        with open(filename,'w') as f:
            for full_loop in interaction_to_list_of_tuples['loops']:
                a = full_loop['identifier']
                b = ",".join(full_loop['unit_ids'])
                c = ",".join(full_loop['border_indicators'])
                f.write('"%s","%s","%s"\n' % (a,b,c))


def write_ebi_json_output_file(outputNAPairwiseInteractions,file_id,interaction_to_list_of_tuples,categories,category_to_interactions,chain,unit_id_to_sequence_position,modified):
    """
    For each chain, write interactions according to category,
    and within each category, write by annotation.
    Other than that, the interactions are listed in no particular order.
    """

    import json

    # loop over types of output files requested
    for category in categories.keys():
        filename = os.path.join(outputNAPairwiseInteractions,file_id + "_" + chain + "_" + category + ".json")

        output = {}
        output["pdb_id"] = file_id
        output["chain_id"] = chain
        output["modified"] = modified

        annotations = []
        for interaction in category_to_interactions[category]:
            inter = interaction
            if "n" in interaction:
                continue
            if "cWB" in interaction or "cBW" in interaction:
                continue
            if category == 'basepair':
                # capitalize base edges
                inter = interaction.replace("w","W").replace("s","S").replace("h","H")
                inter = interaction.replace("a","")
            # if this category has a restricted list of interactions to output
            if len(categories[category]) == 0 or inter in categories[category]:
                for a,b,c in interaction_to_list_of_tuples[interaction]:
                    fields1 = a.split("|")
                    fields2 = b.split("|")
                    if fields1[2] == chain and fields2[2] == chain:
                        if unit_id_to_sequence_position[a] < unit_id_to_sequence_position[b]:
                            ann = {}
                            ann["seq_id1"]  = str(unit_id_to_sequence_position[a])
                            ann["3d_id1"]   = fields1[4]
                            ann["nt1"]      = fields1[3]
                            ann["unit1"]    = fields1[3]
                            ann["bp"]       = inter
                            ann["seq_id2"]  = str(unit_id_to_sequence_position[b])
                            ann["nt2"]      = fields2[3]
                            ann["unit2"]    = fields2[3]
                            ann["3d_id2"]   = fields2[4]
                            ann["crossing"] = str(c)
                            #{"seq_id1":"1","3d_id1":"13","nt1":"C","bp":"cWW","seq_id2":"71","nt2":"G","3d_id2":"83","crossing":"0"}

                            annotations.append(ann)

        output["annotations"] = annotations

        with open(filename,'w') as f:
            f.write(json.dumps(output))


#=======================================================================
def generatePairwiseAnnotation(entry_id, chain_id, inputPath, outputNAPairwiseInteractions, category, output_format):

    if isinstance(entry_id,str):
        entry_id = entry_id.split(",")

    # dictionary to control what specific annotations are output, in a file named for the key
    # empty list means to output all interactions in that category
    # non-empty list specifies which interactions to output in that category
    categories = {}

    if category:
        category_list = category.split(",")
        if len(category_list) == 1 and category_list[0] == "oo_distance":
            categories['oo_distance'] = []
        else:
            for category in category_list:
                c = category.lower()
                if c in all_categories.split(","):
                    categories[category.lower()] = []
                else:
                    print('Category %s is not recognized' % c)
                categories['basepair'] = []         # always basepairs, to get crossing numbers
    else:
        # default is to annotate and write just "true" basepairs
        categories['basepair'] = Leontis_Westhof_basepairs
        categories['basepair'] = []

    if 'loop' in categories or 'loops' in categories:
        categories['loops'] = []

    if 'loops' in categories:
        categories['coplanar'] = []
        categories['bss'] = []
        categories['stacking'] = []
        categories['basepair'] = Leontis_Westhof_basepairs
        categories['basepair'] = []

    if 'bss' in categories:
        categories['basepair_detail'] = []
        categories['basepair'] = Leontis_Westhof_basepairs + ['cWB','cBW']  # bifurcated pairs
        categories['basepair'] = []

    if 'basepair_detail' in categories:
        categories['basepair'] = Leontis_Westhof_basepairs + ['cWB','cBW']  # bifurcated pairs
        categories['basepair'] = []
    elif 'basepair' in categories:
        categories['basepair'] = Leontis_Westhof_basepairs
        categories['basepair'] = []


    # check existence of input path
    if len(inputPath) > 0 and not os.path.exists(inputPath):
        if verbose >= 1:
            print("  Attempting to create input path %s" % inputPath)
        os.mkdir(inputPath)

    # check existence of output path
    if len(outputNAPairwiseInteractions) > 0 and not os.path.exists(outputNAPairwiseInteractions):
        if verbose >= 1:
            print("  Attempting to create output path %s" % outputNAPairwiseInteractions)
        os.mkdir(outputNAPairwiseInteractions)

    # process additional arguments as PDB files
    PDBs = []  # list of (path,filename) entries
    entries = entry_id
    for entry in entries:
        # identify path to the PDB file, if any
        path_split = os.path.split(entry)   # produces a tuple

        if len(path_split[0]) > 0:
            PDBs.append(path_split)
        else:
            PDBs.append((inputPath,entry))

    # annotate each PDB file
    timerData = myTimer("start")
    failed_structures = []
    counter = 0

    if chain_id:
        if len(entry_id) > 1:
            if verbose >= 0:
                print("  Chain argument can only be used with a single PDB file")
            PDBs = []
        else:
            chains = chain_id.split(",")
    else:
        chains = []

    # restrict dictionary of cutoffs to just the basepairs needed here
    if 'basepair' in categories:
        focused_basepair_cutoffs = focus_basepair_cutoffs(nt_nt_cutoffs,categories['basepair'])
        ideal_hydrogen_bonds = load_ideal_basepair_hydrogen_bonds()
    else:
        focused_basepair_cutoffs = {}
        ideal_hydrogen_bonds = {}

    """
    for combination in ideal_hydrogen_bonds:
        for LW in ideal_hydrogen_bonds[combination]:
            print(combination, LW, ideal_hydrogen_bonds[combination][LW])
    """

    for path, PDB in PDBs:
        counter += 1

        # attempt to identify the main file identifier, could be a 4-character pdb id
        file_id = PDB.replace(".cif","").replace(".pdb","").replace(".gz","")

        filename = os.path.join(path,PDB)

        if verbose >= 1:
            print("  Reading file %s, which is number %d out of %d" % (filename, counter, len(PDBs)))
        timerData = myTimer("Read CIF files",timerData)

        # suppress error messages, but report failures at the end
        structure, messages = load_structure(filename,file_id)

        if not structure:
            for message in messages:
                failed_structures.append((file_id,message))
            continue

        interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_in_structure(structure,categories,focused_basepair_cutoffs,ideal_hydrogen_bonds,chains,timerData)
        timerData = myTimer("Record interactions",timerData)
        if verbose >= 1:
            print("  Recording interactions in %s" % outputNAPairwiseInteractions)

        if output_format == 'txt':
            write_txt_output_file(outputNAPairwiseInteractions,file_id,interaction_to_list_of_tuples,categories,category_to_interactions)
        elif output_format == 'ebi_json':
            if chains:
                bases = structure.residues(chain = chains, type = ["RNA linking","DNA linking"])  # load all RNA/DNA nucleotides
            else:
                bases = structure.residues(type = ["RNA linking","DNA linking"])  # load all RNA/DNA nucleotides

            chain_unit_id_to_sequence_position = {}
            chain_modified = {}
            for base in bases:
                chain = base.chain
                if not chain in chain_unit_id_to_sequence_position:
                    chain_unit_id_to_sequence_position[chain] = {}
                    chain_modified[chain] = []
                chain_unit_id_to_sequence_position[chain][base.unit_id()] = base.index

                fields = base.unit_id().split('|')
                if not fields[3] in ['A','C','G','U','DA','DC','DG','DT']:
                    modif = {}
                    modif['seq_id'] = str(base.index)
                    modif['nt1'] = fields[3]
                    modif['unit1'] = fields[3]
                    modif['3d_id'] = fields[4]
                    chain_modified[chain].append(modif)

            for chain in list(chain_unit_id_to_sequence_position.keys()):
                write_ebi_json_output_file(outputNAPairwiseInteractions,file_id,interaction_to_list_of_tuples,categories, category_to_interactions, chain, chain_unit_id_to_sequence_position[chain],chain_modified[chain])

        else:
            if verbose >= 0:
                print('  Output format %s not recognized' % output_format)

    if verbose >= 1:
        myTimer("summary",timerData)
        if len(failed_structures) > 0:
            print("  Error messages:")
            for message in failed_structures:
                print("  %s %s" % message)
        else:
            print("All files read successfully")

if __name__=="__main__":

    # allow user to specify input and output paths
    parser = argparse.ArgumentParser()
    parser.add_argument('PDBfiles', type=str, nargs='+', help='.cif filename(s)')
    parser.add_argument('-o', "--output", help="Output Location of Pairwise Interactions")
    parser.add_argument('-i', "--input", help='Input Path')
    parser.add_argument('-c', "--category", help='Interaction category or categories (%s)' % all_categories)
    parser.add_argument('-f', "--format", help='Output format (txt,ebi_json)')
    parser.add_argument('-v', "--verbose", help='Verbose level (0,1,2,3)')
    parser.add_argument("--chain", help='Chain or chains separated by commas, no spaces; only for one PDB file')

    problem = False
    args = parser.parse_args()

    # Process command line arguments
    if args.input:
        inputPath = args.input
    else:
        if not inputPath:
            inputPath = ""

    if args.output:
        outputNAPairwiseInteractions = args.output     # set output path
    else:
        if not outputNAPairwiseInteractions:
            outputNAPairwiseInteractions = ""

    if args.format:
        outputFormat = args.format
    else:
        outputFormat = 'txt'

    if args.chain:
        chain_id = args.chain
    else:
        chain_id = None

    if args.category:
        category = args.category.replace("-","_")
    else:
        category = 'basepair'

    if "all" in category:
        category = all_categories

    if args.verbose:
        verbose = int(args.verbose)

    entry_id = args.PDBfiles[0].split(",")

    generatePairwiseAnnotation(entry_id, chain_id, inputPath, outputNAPairwiseInteractions, category, outputFormat)

