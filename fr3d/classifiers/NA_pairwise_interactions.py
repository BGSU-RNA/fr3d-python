# -*- coding: utf-8 -*-
"""
This program reads one or more CIF/PDB files and produces annotations
of nucleotide-nucleotide interactions.
Basepairs are annotated with Leontis-Westhof annotations like cWW, tHS.
Basepairs might also be annotated as "near" with ncWW, ntHS; ask for basepair_detail.
A few basepair categories have "alternative" geometries like cWWa.
Alternative geometries are not checked for hydrogen bonds.
cWB is for Table 13 in Leontis-Stombaugh-Westhof, bifurcated interactions.

Usage examples:
python NA_pairwise_interactions.py 4TNA
python NA_pairwise_interactions.py -c basepair_detail 4TNA
python NA_pairwise_interactions.py -c basepair_detail,sugar_ribose 8B0X

When fr3d_python is changed, reinstall with:
change directory to fr3d-python
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

from fr3d.classifiers.class_limits_2023 import nt_nt_cutoffs   # use latest cutoffs
from fr3d.classifiers.hydrogen_bonds import load_ideal_basepair_hydrogen_bonds
from fr3d.classifiers.hydrogen_bonds import check_hydrogen_bond

# Modified nucleotide mappings from atom_mappings_refined.py
from fr3d.data.mapping import modified_base_atom_list,parent_atom_to_modified,modified_atom_to_parent,modified_base_to_parent

# read input and output paths from localpath.py
# note that fr3d.localpath does not synchronize with Git, so you can change it locally to point to your own directory structure
try:
    from fr3d.localpath import outputNAPairwiseInteractions
    from fr3d.localpath import inputPath
except:
    inputPath = ""
    outputNAPairwiseInteractions = ""

nt_nt_screen_distance = 12  # maximum center-center distance to check

near_discrepancy_cutoff = 1.0     # maximum discrepancy to report as a near pair
near_heavy_distance_cutoff = 4.2  # maximum distance between heavy atoms to be considered a near pair
true_heavy_distance_cutoff = 3.8  # maximum distance between heavy atoms to be considered a true pair

HB_donor_hydrogens = {}
HB_donor_hydrogens['A'] = {"N6":["1H6","2H6"], "C2":["H2"], "C8":["H8"], "O2'":[]}
HB_donor_hydrogens['G'] = {"N1":["H1"], "N2":["2H2","1H2"], "C8":["H8"], "O2'":[]}
HB_donor_hydrogens['C'] = {"N4":["1H4","2H4"], "C5":["H5"], "C6":["H6"], "O2'":[]}
HB_donor_hydrogens['U'] = {"N3":["H3"], "C5":["H5"], "C6":["H6"], "O2'":[]}

standard_bases = ['A','C','G','U','DA','DC','DG','DT']

nt_reference_point = "base"
atom_atom_min_distance = 5    # minimum distance between atoms in nts to consider them interacting
base_seq_list = []                     # for all nucleic acids, modified or not


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

    lower_interactions = set([])

    if interactions:
        for interaction in interactions:
            lower_interactions.add(interaction.lower())
    else:
        # use all available interactions
        for combination in basepair_cutoffs.keys():
            for interaction in basepair_cutoffs[combination]:
                lower_interactions.add(interaction.lower())

    for combination in basepair_cutoffs.keys():
        focused_basepair_cutoffs[combination] = {}
        focused_basepair_cutoffs[combination][1] = {}   # interactions with positive normal
        focused_basepair_cutoffs[combination][-1] = {}  # interactions with negative normal
        for interaction in basepair_cutoffs[combination].keys():
            family = interaction.lower()[0:3]   # take off alternative category a, b, etc.

            if family in lower_interactions:
                subcat = list(basepair_cutoffs[combination][interaction].keys())[0]
                if basepair_cutoffs[combination][interaction][subcat]["normalmin"] > 0:
                    focused_basepair_cutoffs[combination][1][interaction] = {}   # interactions with positive normal
                    for subcategory in basepair_cutoffs[combination][interaction]:
                        focused_basepair_cutoffs[combination][1][interaction][subcategory] = basepair_cutoffs[combination][interaction][subcategory]

                else:
                    focused_basepair_cutoffs[combination][-1][interaction] = {}   # interactions with negative normal
                    for subcategory in basepair_cutoffs[combination][interaction]:
                        focused_basepair_cutoffs[combination][-1][interaction][subcategory] = basepair_cutoffs[combination][interaction][subcategory]

    """
    # check how this worked
    for combination in focused_basepair_cutoffs:
        for normal in focused_basepair_cutoffs[combination]:
            for interaction in focused_basepair_cutoffs[combination][normal]:
                for subcategory in focused_basepair_cutoffs[combination][normal][interaction]:
                    print(combination, normal, interaction, subcategory, focused_basepair_cutoffs[combination][normal][interaction][subcategory])
                    pass
    """

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

        url = "http://files.rcsb.org/download/%s" % download_id

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
                structure = Cif(raw,preferred_id=preferred_id).structure()
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
                print("  No symmetry operators applied to .pdb files")
        elif filename.lower().endswith('.pdb'):
            with open(filename, rm) as raw:
                from fr3d.pdb.pdb_reader import PDBStructure
                structure = PDBStructure(file_id,raw).structures()
                print("  No symmetry operators applied to .pdb files")

        message.append("Loaded " + filename)
        return structure, message

    except Exception as ex:
        message.append("  Could not load %s due to exception %s: %s" % (filename,type(ex).__name__,ex))
        if type(ex).__name__ == "TypeError":
            message.append("  See suggestions in the fr3d-python Readme file")
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

#    print(atom_to_part_list)

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
    if seq in parent_atom_to_modified:
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
    else:
        rev = inter[0:(len(inter)-2)] + inter[len(inter)-1] + inter[len(inter)-2]

    return rev

def makeListOfNtIndices(baseCubeList, baseCubeNeighbors):
    """This function returns a sorted list of all the nts indices in ascending order.
    It was added as a method to be able to extract information about the O3' atom of the previous nucleotide"""
    lastNT = {}
    for nt1key in baseCubeList:                         # key to first cube
        for nt2key in baseCubeNeighbors[nt1key]:        # key to each potential neighboring cube, including the first
            if nt2key in baseCubeList:                  # if this cube was actually made
                for nt1 in baseCubeList[nt1key]:
                    if nt1.index not in lastNT:
                        lastNT[nt1.index] = nt1
    return(lastNT)


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
        t = (base.model,base.symmetry,base.chain,base.index,base.unit_id(),coordinates,P)
        list_of_nucleotides.append(t)

    list_of_nucleotides.sort()

    previous_O3_coordinates = np.empty([1,3])
    previous_model = None
    previous_symmetry = None
    previous_chain = None
    previous_index = None
    unit_id_to_previous_O3 = {}

    for model,symmetry,chain,index,unit_id,O3_coordinates,P in list_of_nucleotides:

        #print(model,symmetry,chain,index,unit_id,O3_coordinates)

        if model == previous_model and symmetry == previous_symmetry and chain == previous_chain and index == previous_index + 1:
            unit_id_to_previous_O3[unit_id] = previous_O3_coordinates

            # if P.any():
            #     print("%-20s distance from P to previous O3' is %f" % (unit_id,np.linalg.norm(P-previous_O3_coordinates)))
            # else:
            #     print(P)
        else:
            unit_id_to_previous_O3[unit_id] = np.empty([1,3])
            # print("%-20s has no previous O3' see http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s" % (unit_id,unit_id))

        # save current values for next nucleotide
        previous_model = model
        previous_symmetry = symmetry
        previous_chain = chain
        previous_index = index
        previous_O3_coordinates = O3_coordinates

    return unit_id_to_previous_O3


def check_for_two_interactions_on_same_edge(unit_id_to_basepairs,get_datapoint=False):
    """
    Loop over nucleotides, find those with two or more interactions on the same edge, choose the best, remove the others
    """

    # set of tuples of unit ids to leave out of the basepair list
    remove_pairs = set()
    make_near_pairs = set()

    for unit_id, basepairs in unit_id_to_basepairs.items():
        if len(basepairs) > 1:
            for i in range(len(basepairs)-1):
                basepair = basepairs[i]
                interaction_1, quality_1, unit_id_1 = basepairs[i]

                e1 = interaction_1.replace("n","")[1].lower()   # base edge

                if (unit_id,unit_id_1) in remove_pairs:
                    continue

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

                    common_atoms = set(quality_1['atoms1']) & set(quality_2['atoms1'])

                    if len(common_atoms) > 0:

                        if get_datapoint:
                            unit_id_list = unit_id
                            print("  Base %s makes multiple basepairs listed %d and %d below" % (unit_id,i,j))

                            for bp in basepairs:
                                print(bp)
                                interaction, quality, u1 = bp
                                unit_id_list += "," + u1

                            print("  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s" % unit_id_list)
                            print('  Common atoms %s' % common_atoms)

                        atom_names = "".join(common_atoms)
                        if "H" in atom_names or len(common_atoms) > 1:

                            if get_datapoint:
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
                                        print("  %-5s %-22s %-22s  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  removed due to conflicting edge, max gap" % (interaction_2,unit_id,unit_id_2,unit_id,unit_id_2))
                                else:
                                    if quality_1['cutoff_distance'] > 0.5 * near_discrepancy_cutoff:
                                        remove_pairs.add((unit_id,unit_id_1))
                                        remove_pairs.add((unit_id_1,unit_id))
                                        print("  %-5s %-22s %-22s  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  removed due to conflicting edge, max gap" % (interaction_1,unit_id,unit_id_1,unit_id,unit_id_1))

                            elif not interaction_1.startswith("n") and not interaction_2.startswith("n"):
                                # both true, make one near
                                if quality_1['max_gap'] < quality_2['max_gap']:
                                    make_near_pairs.add((unit_id,unit_id_2))
                                    make_near_pairs.add((unit_id_2,unit_id))
                                    print("  %-5s %-22s %-22s  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  switched to near due to conflicting edge, max gap" % (interaction_2,unit_id,unit_id_2,unit_id,unit_id_2))
                                else:
                                    make_near_pairs.add((unit_id,unit_id_1))
                                    make_near_pairs.add((unit_id_1,unit_id))
                                    print("  %-5s %-22s %-22s  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  switched to near due to conflicting edge, max gap" % (interaction_1,unit_id,unit_id_1,unit_id,unit_id_1))

                            else:
                                # one near, one true, remove the near one if it's bad
                                if quality_2['cutoff_distance'] > 0.5 * near_discrepancy_cutoff:
                                    remove_pairs.add((unit_id,unit_id_2))
                                    remove_pairs.add((unit_id_2,unit_id))
                                    print("  %-5s %-22s %-22s  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  removed due to conflicting edge, cutoffs" % (interaction_2,unit_id,unit_id_2,unit_id,unit_id_2))
                                elif quality_1['cutoff_distance'] > 0.5 * near_discrepancy_cutoff:
                                    remove_pairs.add((unit_id,unit_id_1))
                                    remove_pairs.add((unit_id_1,unit_id))
                                    print("  %-5s %-22s %-22s  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  removed due to conflicting edge, cutoffs" % (interaction_1,unit_id,unit_id_1,unit_id,unit_id_1))

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

    basepair_parent_base_combination_set = set(['A,A','A,C','A,G','A,U','C,C','G,C','C,U','G,G','G,U','U,U','A,DT','C,DT','G,DT','DT,DT'])

    # For base-backbone interactions, we need to know
    if 'backbone' in categories.keys():
        #ntDict = makeListOfNtIndices(baseCubeList, baseCubeNeighbors)
        unit_id_to_previous_O3 = map_unit_id_to_previous_O3(bases)

    for nt1key in baseCubeList:                         # key to first cube
        for nt2key in baseCubeNeighbors[nt1key]:        # key to each potential neighboring cube, including the first
            if nt2key in baseCubeList:                  # if this cube was actually made
                for nt1 in baseCubeList[nt1key]:        # first nt of a potential pair

                    if len(nt1.centers["base"]) < 3:
                        print("  Missing base center for %s" % nt1.unit_id())
                        print(nt1.centers["base"])
                        continue

                    parent1 = get_parent(nt1.sequence)   # map modified nts to parent nt
                    gly1 = get_glycosidic_atom_coordinates(nt1,parent1)

                    if len(gly1) < 3:
                        print("  Missing glycosidic atom for %s" % nt1.unit_id())
                        continue

                    number1 = nt1.number                 # nucleotide number

                    for nt2 in baseCubeList[nt2key]:           # second nt of a potential pair
                        # only consider each nt1, nt2 pair in one direction
                        # Those in different cubes only occur once
                        # Those from the same cube need a way to select just one pair
                        if nt1key == nt2key:
                            if nt1.chain > nt2.chain:
                                continue
                            elif nt1.chain == nt2.chain and nt1.index > nt2.index:
                                continue

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

                        # avoid comparing alternate coordinates of the same nucleotide
                        # check in the order most likely to terminate the fastest
                        if number1 == nt2.number:
                            if nt1.sequence == nt2.sequence:
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

                        # some structures have overlapping nucleotides, screen, those out
                        if center_center_distance < 2:
                            continue

                        unit_id_pair = (nt1.unit_id(),nt2.unit_id())  # tuple for these nucleotides in this order
                        reversed_pair = (nt2.unit_id(),nt1.unit_id())

                        parent2 = get_parent(nt2.sequence)
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
                            datapoint12['url'] = "http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt1.unit_id(),nt2.unit_id())

                            datapoint21 = {}
                            datapoint21['center_center_distance'] = center_center_distance
                            datapoint21['nt1_seq'] = nt2.sequence
                            datapoint21['nt2_seq'] = nt1.sequence
                            datapoint21['nt1_parent'] = parent2
                            datapoint21['nt2_parent'] = parent1
                            datapoint21['url'] = "http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt2.unit_id(),nt1.unit_id())

                        else:
                            datapoint12 = None
                            datapoint21 = None

                        # check base to oxygen stack; always base first, oxygen second
                        if 'sO' in categories.keys():
                            timerData = myTimer("Check base oxygen stack",timerData)
                            interaction, datapoint12, interaction_reversed = check_base_oxygen_stack_rings(nt1,nt2,parent1,datapoint12)

                            if len(interaction) > 0:
                                count_pair += 1
                                interaction_to_pair_list[interaction].append(unit_id_pair)
                                interaction_to_pair_list[interaction_reversed].append(reversed_pair)
                                max_center_center_distance = max(max_center_center_distance,center_center_distance)  # for setting optimally
                                category_to_interactions['sO'].add(interaction)
                                category_to_interactions['sO'].add(interaction_reversed)

                            interaction, datapoint21, interaction_reversed = check_base_oxygen_stack_rings(nt2,nt1,parent2,datapoint21)

                            if len(interaction) > 0:
                                count_pair += 1
                                interaction_to_pair_list[interaction].append(reversed_pair)
                                interaction_to_pair_list[interaction_reversed].append(unit_id_pair)
                                max_center_center_distance = max(max_center_center_distance,center_center_distance)  # for setting optimally
                                category_to_interactions['sO'].add(interaction)
                                category_to_interactions['sO'].add(interaction_reversed)

                        if 'stacking' in categories.keys():
                            timerData = myTimer("Check base base stack", timerData)
                            interaction, datapoint12, interaction_reversed = check_base_base_stacking(nt1, nt2, parent1, parent2, datapoint12)
                            if len(interaction) > 0:
                                count_pair += 1
                                interaction_to_pair_list[interaction].append(unit_id_pair)
                                max_center_center_distance = max(max_center_center_distance,center_center_distance)  # for setting optimally
                                category_to_interactions['stacking'].add(interaction)
                                category_to_interactions['stacking'].add(interaction_reversed)

                        # annotate sugar ribose interactions;
                        if 'sugar_ribose' in categories.keys():
                            timerData = myTimer("Check sugar ribose", timerData)
                            if not parent1 in ['DA','DC','DG','DT'] and not parent2 in ['DA','DC','DG','DT']:
                                interaction, datapoint12 = check_sugar_ribose(nt1, nt2, parent1, datapoint12)
                                if len(interaction) > 0:
                                    count_pair += 1
                                    interaction_to_pair_list[interaction].append(unit_id_pair)
                                    category_to_interactions['sugar_ribose'].add(interaction)

                                interaction, datapoint21 = check_sugar_ribose(nt2, nt1, parent2, datapoint21)
                                if len(interaction) > 0:
                                    count_pair += 1
                                    interaction_to_pair_list[interaction].append(reversed_pair)
                                    category_to_interactions['sugar_ribose'].add(interaction)

                        # annotate base phosphate and base ribose interactions
                        if 'backbone' in categories.keys():
                            timerData = myTimer("Check backbone interactions", timerData)

                            # # you need the O3' atom of the last nucleotide and this dict will help you get that component.
                            # lastNT = None
                            # lastNT2 = None
                            # if nt1.index - 1 > 0 and (nt1.index-1) in ntDict:
                            #     lastNT = ntDict[nt1.index-1]
                            # if nt2.index - 1 > 0 and (nt2.index-1) in ntDict:
                            #     lastNT2 = ntDict[nt2.index-1]

                            # get coordinates of O3' of the nucleotide before nt2, part of the phosphate of nt2
                            previousO3 = unit_id_to_previous_O3.get(nt2.unit_id(),np.empty([1,3]))
                            interactionbPh, interactionbR, datapoint12 = check_base_backbone_interactions(nt1, nt2, previousO3, parent1, parent2, datapoint12)

                            if interactionbPh and len(interactionbPh) > 0:
                                count_pair += 1
                                interaction_to_pair_list[interactionbPh].append(unit_id_pair)
                                category_to_interactions['backbone'].add(interactionbPh)
                            if interactionbR and len(interactionbR) > 0:
                                count_pair += 1
                                interaction_to_pair_list[interactionbR].append(unit_id_pair)
                                category_to_interactions['backbone'].add(interactionbR)

                            #     max_center_center_distance = max(max_center_center_distance,center_center_distance)  # for setting optimally


                        gly2 = get_glycosidic_atom_coordinates(nt2,parent2)
                        if len(gly2) < 3:
                            print("  Missing glycosidic atom for %s" % nt2.unit_id())
                            continue

                        # always annotate cWW basepairs to be able to calculate crossing numbers
                        # check coplanar and basepairing for bases in specific orders
                        # AA, CC, GG, UU will be checked in both nucleotide orders, that's OK
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

                                # annotate coplanar relationship
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

                            # record basepairs made by modified nucleotides
                            if False and len(interaction) > 0 and not (nt1.sequence in standard_bases and nt2.sequence in standard_bases):
                                print('%s\t%s\t%s\t%s\t%s\t%s\t=hyperlink("http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s")' % (nt1.sequence,interaction[0],nt2.sequence,nt1.unit_id(),nt2.unit_id(),interaction,nt1.unit_id(),nt2.unit_id()))
                                try:
                                    with open('C:/Users/zirbel/Documents/FR3D/Modified Nucleotides/list.txt','a') as file:
                                        file.write('%s\t%s\t%s\t%s\t%s\t%s\t=hyperlink("http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s")\n' % (nt1.sequence,interaction[0],nt2.sequence,nt1.unit_id(),nt2.unit_id(),interaction,nt1.unit_id(),nt2.unit_id()))
                                except:
                                    pass

                            if False and len(interaction12) > 0 and 'gap12' in pair_data:
                                print("  Identified parents as %s and %s" % (parent1,parent2))
                                print("  Found %s interaction between %-18s and %-18s" % (interaction12,nt1.unit_id(),nt2.unit_id()))
                                print("  Gap value %0.8f" % pair_data["gap12"])
                                #print("  Coplanar Boolean %s" % pair_data["coplanar"])
                                #if 'coplanar_value' in pair_data and pair_data["coplanar_value"]:
                                #    print("  Coplanar value %0.8f" % pair_data["coplanar_value"])
                                print("  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt1.unit_id(),nt2.unit_id()))
                                print("")

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

                            if 'coplanar' in categories.keys():
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

                            if False and len(interaction21) > 1 and 'gap12' in pair_data:
                                print("  Identified parents as %s and %s" % (parent2,parent1))
                                print("  Found %s interaction between %-18s and %-18s" % (interaction21,nt2.unit_id(),nt1.unit_id()))
                                print("  Gap value %0.8f" % pair_data["gap12"])
                                #print("  Coplanar Boolean %s" % pair_data["coplanar"])
                                #if 'coplanar_value' in pair_data and pair_data["coplanar_value"]:
                                #    print("  Coplanar value %0.8f" % pair_data["coplanar_value"])
                                print("  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt1.unit_id(),nt2.unit_id()))
                                print("")

                        else:
                            interaction21 = ""
                            interaction21_reversed = ""

                        # if annotated interaction in both pair orders, choose the better one
                        if len(interaction12) > 0 and len(interaction21) > 0:
                            conflict_message = "%-5s %-22s %-22s  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  duplicate interaction %-5s in second direction" % (interaction21,nt1.unit_id(),nt2.unit_id(),nt1.unit_id(),nt2.unit_id(),interaction12_reversed)

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
                            interaction_clean_reversed = reverse_edges(interaction_clean)

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
    remove_pairs, make_near_pairs = check_for_two_interactions_on_same_edge(unit_id_to_basepairs,get_datapoint)

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

    print("  Found %d nucleotide-nucleotide interactions" % count_pair)

    if False:
        print("  Maximum screen distance for actual contacts is %8.4f" % max_center_center_distance)

    # calculate and save crossing numbers for each annoated interaction
    timerData = myTimer("Calculate crossing",timerData)
    interaction_to_list_of_tuples = calculate_crossing_numbers(bases,interaction_to_pair_list)

    return interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data

def calculate_crossing_numbers(bases,interaction_to_pair_list):
    """
    Identify which cWW pairs are nested.
    Then for each interaction, calculate the number of nested cWW pairs it crosses
    """

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

    # find AU, GC, GU cWW basepairs within each chain
    for interaction in ['cWW','cWw','cwW','acWW','acWw','acwW']:
        for u1,u2 in interaction_to_pair_list[interaction]:
            chain1, index1 = unit_id_to_index[u1]
            chain2, index2 = unit_id_to_index[u2]

            # record AU, GC, GU cWW pairs by index within each chain
            if chain1 == chain2:
                fields = u1.split('|')
                parent1 = get_parent(fields[3])
                fields = u2.split('|')
                parent2 = get_parent(fields[3])

                if parent1+parent2 in ['AU','UA','CG','GC','GU','UG']:
                    if index1 < index2:
                        chain_to_cWW_pairs[chain1].append((index1,index2))
                    else:
                        chain_to_cWW_pairs[chain1].append((index2,index1))

    chain_nested_cWW_endpoints = {}

    # within each chain, sort nested cWW by distance between them
    # starting with the shortest-range pairs, record nested cWW pairs
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
    interaction_to_list_of_tuples = defaultdict(list)

    # loop over pairs, calculate crossing number
    # record interacting pairs and their crossing number as triples
    for interaction in interaction_to_pair_list.keys():

        if interaction == "":
            continue

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

            interaction_to_list_of_tuples[interaction].append((u1,u2,crossing))

            # duplicate certain pairs in reversed order; saves time this way
            if interaction in ["s33","s35","s53","s55","cp","ns33","ns35","ns53","ns55"]:
                interaction_to_list_of_tuples[reverse_edges(interaction)].append((u2,u1,crossing))
            elif interaction[0] in ["c","t"]:
                interaction_to_list_of_tuples[reverse_edges(interaction)].append((u2,u1,crossing))
            elif interaction[0:2] in ["nc","nt"]:
                interaction_to_list_of_tuples[reverse_edges(interaction)].append((u2,u1,crossing))

    return interaction_to_list_of_tuples

def annotate_covalent_connections(nucleotides, interaction_to_list_of_tuples, category_to_interactions, timerData):
    """
    Loop through bases, sort by model, symmetry, chain, and
    record the distance in the chain between successive
    observed nucleotides.
    """

    nts_to_sort = defaultdict(list)

    # list the nucleotides by model, symmetry, chain, index
    for nt in nucleotides:
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
                    interaction_to_list_of_tuples[interaction].append((u1,u2,None))
                    category_to_interactions["covalent"].add(interaction)
                    #print("%s\t%s\t%s" % (u1,interaction,u2))

    return interaction_to_list_of_tuples, category_to_interactions, timerData


def annotate_nt_nt_in_structure(structure,categories,focused_basepair_cutoffs={},ideal_hydrogen_bonds={},chains=[],timerData=None,get_datapoint=False):
    """
    This function can be called from the pipeline to annotate a structure
    structure is an output from
    """

    if not focused_basepair_cutoffs:
        focused_basepair_cutoffs = focus_basepair_cutoffs(nt_nt_cutoffs,categories['basepair'])

    if not ideal_hydrogen_bonds:
        ideal_hydrogen_bonds = load_ideal_basepair_hydrogen_bonds()


    if chains:
        bases = structure.residues(chain = chains, type = ["RNA linking","DNA linking"])  # load all RNA/DNA nucleotides
    else:
        bases = structure.residues(type = ["RNA linking","DNA linking"])  # load all RNA/DNA nucleotides

    if not timerData:
        timerData = myTimer("start")

    timerData = myTimer("Building cubes",timerData)
    baseCubeList, baseCubeNeighbors = make_nt_cubes_half(bases, nt_nt_screen_distance, nt_reference_point)
    # annotate nt-nt interactions
    timerData = myTimer("Annotating interactions",timerData)
    interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_interactions(bases, nt_nt_screen_distance, baseCubeList, baseCubeNeighbors, categories, focused_basepair_cutoffs, ideal_hydrogen_bonds, timerData, get_datapoint)

    # annotate covalent connections
    interaction_to_list_of_tuples, category_to_interactions, timerData = annotate_covalent_connections(bases, interaction_to_list_of_tuples, category_to_interactions, timerData)

    return interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data


def get_parent(sequence):
    """
    Look up parent sequence for RNA, DNA, and modified nucleotides
    """

    if sequence in ['A','C','G','U']:
        return sequence
    elif sequence in ['DA','DC','DG']:
        return sequence[1]
    elif sequence == 'DT':
        return sequence
    elif sequence in modified_base_to_parent:
        return modified_base_to_parent[sequence]
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


def check_base_oxygen_stack_rings(nt1,nt2,parent1,datapoint):
    '''
    Does one of the backbone oxygens of nt2 stack inside a ring on the base of nt1?
    '''

    true_z_cutoff = 3.5
    near_z_cutoff = 3.6
    outside_z_cutoff = (true_z_cutoff + near_z_cutoff)/2

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
        print('%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\t\t=hyperlink("http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s")' % (nt1.unit_id(),nt2.unit_id(),interaction,xmin,ymin,zmin,nt1.unit_id(),nt2.unit_id()))

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

            # check to see if a nt has points on both sides of the plane of a nt. See http://rna.bgsu.edu/rna3dhub/display3D/unitid/6ZMI%7C1%7CL5%7CG%7C2605,6ZMI%7C1%7CL5%7CG%7C2668 for an example.
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
        print("  Can't check base stacking for %s and %s" % (nt1.unit_id(),nt2.unit_id()))
        return "", datapoint, ""

    nt2baseAtomsList = get_base_atom_names(nt2.sequence)

    if len(nt2baseAtomsList) == 0:
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

    if len(interaction) > 0 and False:
        print('%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\t\t=hyperlink("http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s")' % (nt1.unit_id(),nt2.unit_id(),interaction,coords[0],coords[1],coords[2],nt1.unit_id(),nt2.unit_id()))

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

                    # if "6ERI|1|AA|U|1241" in nt1.unit_id() or "6ERI|1|AA|U|1241" in nt2.unit_id():
                    #     if "6ERI|1|AA|A|558" in nt1.unit_id() or "6ERI|1|AA|A|558" in nt2.unit_id():
                    #         print("Minimum at %s" % (atom.name))

    # investigate a concern with this calculation
    if "6ERI|1|AA|U|1241" in nt1.unit_id() or "6ERI|1|AA|U|1241" in nt2.unit_id():
        if "6ERI|1|AA|A|558" in nt1.unit_id() or "6ERI|1|AA|A|558" in nt2.unit_id():
            print("  %s %s base_min_distance = %f" % (nt1.unit_id(),nt2.unit_id(),base_min_distance))
            print("  %s %s heavy_min_distance = %f" % (nt1.unit_id(),nt2.unit_id(),heavy_min_distance))
            # for atom1 in nt1.atoms():
            #     if atom1.name in base1_atoms:
            #         print("atom1.name = %s" % atom1.name)
            # for atom2 in nt2.atoms():
            #     if atom2.name in base2_atoms:
            #         print("atom2.name = %s" % atom2.name)

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
    """Method used to add modified nucleotides to a dictionary that is used for processing in function check_base_backbone_interactions.
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
    Function to check base backbone interactions
    nt1 base checked for hydrogen bonds with phosphate and ribose of nt2
    """

    # annotations to return
    phosphate = ""
    ribose = ""

    # if the bases are far away from one another, don't check for base backbone interactions
    dis = distance_between_vectors(nt1.centers["base"], nt2.centers["base"])
    if abs(dis) < 16:   # Initial cutoff

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
            #print('Found O3 of nucleotide before %s' % nt2.unit_id())

        # if nt2 has a P atom and it is far from the plane of base 1, don't look for BPh interactions
        Pcoord = get_one_atom_coordinates(nt2, "P")
        if Pcoord.any():
            try:
                p_standard = translate_rotate_point(nt1, Pcoord)
                if abs(p_standard[2]) > 4.5: # phosphorus far from plane
                    phosphateOxygens = []
            except:
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

            #Loop through the oxygens in the phosphate backbone to extract info for base phosphate interactions
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

            #Loop through oxygens in ribose to extract info about angle and distance
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
                #print(true_phosphate)
                #print(near_phosphate)
                datapoint['BPh'] = phosphate
                if phosphate in ['4BPh','8BPh']:
                    datapoint['BPh_oxygen'] = []
                    for quality,site,oxygen_coordinates,i in true_phosphate:
                        if i in distinct_oxygens:
                            a = translate_rotate_point(nt1, oxygen_coordinates)
                            datapoint['BPh_oxygen'].append(a)
                            #print('%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\t\thttp://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s' % (nt1.unit_id(),nt2.unit_id(),phosphate,a[0],a[1],a[2],nt1.unit_id(),nt2.unit_id()))
                else:
                    a = translate_rotate_point(nt1, phosphate_oxygen)
                    datapoint['BPh_oxygen'] = [a]
                    #print('%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\t\thttp://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s' % (nt1.unit_id(),nt2.unit_id(),phosphate,a[0],a[1],a[2],nt1.unit_id(),nt2.unit_id()))

            if ribose:
                #print(true_ribose)
                #print(near_ribose)
                datapoint['BR'] = ribose
                a = translate_rotate_point(nt1, ribose_oxygen)
                datapoint['BR_oxygen'] = [a]
                #print('%s\t%s\t%s\t%0.4f\t%0.4f\t%0.4f\t\thttp://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s' % (nt1.unit_id(),nt2.unit_id(),ribose,a[0],a[1],a[2],nt1.unit_id(),nt2.unit_id()))

            #datapoint['gap12'], base_points2, atomname2 = calculate_basepair_gap(nt1,nt2)
            #datapoint['url'] = "http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt1.unit_id(),nt2.unit_id())

    return phosphate, ribose, datapoint

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

    #print(len(indices))
    #print(displacements)
    #print(nt1.unit_id())
    #print(nt2.unit_id())

    #units = ["7O7Y|1|A2|A|1081","7O7Y|1|A2|PSU|1082"]
    #units = ["8BUU|1|a|A|76","8BUU|1|a|U|77"]

    gap12 = 100
    for k in range(0,m):              # 3 nearest points
        p = displacements[indices[k]]
        z = abs(np.dot(p,nt1.rotation_matrix[:,2])[0,0])  # distance out of plane of nt1
        if z < gap12:
            gap12 = z                 # gap is smallest z value

        #if nt1.unit_id() in units and nt2.unit_id() in units:
        #    print("Atom %s of %s is %8.4f from plane of %s" % (atomname[indices[k]],nt2.unit_id(),z,nt1.unit_id()))

    """
    if nt1.unit_id() in units and nt2.unit_id() in units:
        print(base_points2)
        print(displacements)
        print(distances)
        print(nt1.centers["base"])
        for atom in nt1.atoms():
            print(nt1.unit_id(),atom.name,atom.x,atom.y,atom.z)
        print("Compare %-24s %-24s gap %8.4f nt1_x %8.3f nt2_x %8.3f" % (nt1.unit_id(), nt2.unit_id(), gap12, nt1.centers["C1'"][0], nt2.centers["C1'"][0]))
    """

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
        print('  Not able to calculate orientation of SR annotation for %s-%s' % (nt1.unit_id(),nt2.unit_id()))
        return "", datapoint

    if datapoint:
        datapoint['sugar_ribose'] = annotation

    # if annotation in ['cSR','tSR']:
    #     print('%s %8.4f %8.4f %8.4f %8.4f %s-%s http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s http://rna.bgsu.edu/correspondence/variability?id=%s,%s&format=unique' % (annotation,o2p_o2p_distance,base_o2p_distance,z,angle,nt1.sequence,nt2.sequence,nt1.unit_id(),nt2.unit_id(),nt1.unit_id(),nt2.unit_id()))

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

    # check possible basepairs two different ways
    for check in check_order:
        if check == 'cutoffs':
            # check cutoffs for each interaction type
            if datapoint:
                cutoff_distance_max = 10.0    # keep checking up to this number to have the data
            else:
                cutoff_distance_max = near_discrepancy_cutoff  # faster annotation

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

                    if cutoff_distance >= cutoff_distance_max:
                        continue

                    cutoff_distance += max(0,cut['zmin'] - displ[0,2])  # how far below zmin
                    cutoff_distance += max(0,displ[0,2] - cut['zmax'])  # how far above zmax

                    if cutoff_distance >= cutoff_distance_max:
                        continue

                    # accentuate wrong normal vector by factor of 3
                    cutoff_distance += 3*max(0,cut['normalmin'] - normal_Z)  # how far below normalmin
                    cutoff_distance += 3*max(0,normal_Z - cut['normalmax'])  # how far above normalmax

                    # cutoffs are met or close enough for now
                    if cutoff_distance < cutoff_distance_max:
                        ok_normal_displ.append((interaction,subcategory,cut,cutoff_distance)) # ("cWW",0), etc.

            # if not close to meeting any cutoffs and we are not collecting data, return now to save time
            if len(ok_normal_displ) == 0 and not datapoint:
                return "", "", quality, datapoint

            if False and datapoint and len(possible_interactions) > 0:
                print("\nhttp://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt1.unit_id(),nt2.unit_id()))

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

            for interaction,subcategory,cut,cutoff_distance in ok_normal_displ:
                if cut['anglemin'] < cut['anglemax']:     # for ranges within -90 to 270 like 50 to 120
                    cutoff_distance += 0.05*max(0,cut['anglemin'] - angle_in_plane)  # how far below anglemin
                    cutoff_distance += 0.05*max(0,angle_in_plane - cut['anglemax'])  # how far above anglemax
                else:                                     # for ranges straddling 270 like 260 to -75
                    cutoff_distance += 0.05*min(max(0,cut["anglemin"]-angle_in_plane),max(0,angle_in_plane-cut["anglemax"]))

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
                if cut['gapmax'] > 0.1:
                    # accentuate wrong gap by factor of 4
                    cutoff_distance += 4*max(0,max(pair_data["gap12"],pair_data["gap21"])-cut['gapmax'])  # how far above gapmax

                # identify cases where there is no base-base hydrogen bond
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
                elif pair_data['heavy_min_distance'] > true_heavy_distance_cutoff and not cSS_one_hbond:
                    # matches a true category but the bases are too far apart for a good basepair
                    near_match.append([interaction,subcategory,cutoff_distance])
                    check_hbonds.append(interaction)
                    print("  %-5s %-22s %-22s  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  switched to near due to minimum distance" % (interaction,nt1.unit_id(),nt2.unit_id(),nt1.unit_id(),nt2.unit_id()))
                else:
                    # true pair
                    match.append([interaction,subcategory,cutoff_distance])
                    check_hbonds.append(interaction)

            if not datapoint:
                # having just checked cutoffs, we want to check hydrogen bonds
                # only for the interactions in match and near_match
                # Restrict possible interactions to those that meet all cutoffs
                possible_interactions = check_hbonds

        else:
            # check hydrogen bonds
            #print("Checking hydrogen bonds for http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt1.unit_id(),nt2.unit_id()))

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

                        # old code to avoid choosing a worse acceptor for a given donor-hydrogen pair
                        # store the lowest "badness" for this donor-hydrogen pair
                        # over all acceptors; avoids identifying a worse option
                        # if result["bond_checked"]:
                        #     donor     = atom_set[0]
                        #     hydrogen  = atom_set[1]
                        #     direction = atom_set[3]

                        #     # also calculate a boolean to tell whether the bond meets the
                        #     # specific criteria to form a hydrogen bond for this LW family

                        #     # store according to donor and direction to disqualify worse acceptors
                        #     if (donor,hydrogen,direction) in donor_hydrogen_to_badness:
                        #         if result["badness"] < donor_hydrogen_to_badness[(donor,hydrogen,direction)]:
                        #             donor_hydrogen_to_badness[(donor,hydrogen,direction)] = result["badness"]
                        #     else:
                        #         donor_hydrogen_to_badness[(donor,hydrogen,direction)] = result["badness"]

                #             if datapoint:
                #                 message = '%4s %-4s, %-3s, %-4s, %s bond, distance %6.3f, hydrogen angle %6.1f, badness %6.3f, donor-acceptor %8s, heavy distance %6.3f, angle atoms %12s, heavy angle %6.1f' % (LW,atom_set[0],atom_set[1],atom_set[2],atom_set[3],result["distance"],result["angle"],result["badness"],result["donor_acceptor_atoms"],result["donor_acceptor_distance"],result["heavy_donor_acceptor_atoms"],result["heavy_donor_acceptor_angle"])
                #                 # print(message)
                #                 LW_bonds[LW].append(result)
                #                 LW_bond_messages[LW].append(message)

                # elif datapoint:
                #     message = 'No hydrogen bonds to check for %s %s %s' % (nt1.unit_id(),nt2.unit_id(),LW)
                #     LW_bond_messages[LW].append(message)
                    #print(message)

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
                                    # print('Found a hydrogen bond for %s with %s' % (LW,atom_set))
                                # else:
                                #     print('!!!!!!!!!!! Rejection due to bond angle')

                            # elif datapoint:
                            #     message = 'Rejected %s, %s, %s, %s bond with distance %0.4f, angle %0.4f, badness %0.4f' % (atom_set[0],atom_set[1],atom_set[2],atom_set[3],result["distance"],result["angle"],result["badness"])
                            #     print(message)
                            #     # LW_bond_messages[LW].append(message)

                    # if datapoint:
                    #     message = '%s has %d out of %s checked hydrogen bonds' % (LW,bond_counter,checked_counter)
                    #     #print(message)
                    #     LW_bond_messages[LW] = [message] + LW_bond_messages[LW]

                    # this is some older diagnostics
                    if checked_counter > 0:
                        if len(badnesses) == 1:
                            # single h-bond interactions are compared to each other only
                            second_badness = badnesses[0] + 100.0
                        else:
                            # two or more h-bond interactions are preferred over one
                            second_badness = sorted(badnesses)[1]
                        LW_bond_counter.append((LW,bond_counter,checked_counter,second_badness))

                    # if bond_counter == 0 and checked_counter > 0:
                    #     if datapoint:
                    #         message = 'Rejecting %s basepair since it has no hydrogen bonds' % LW
                    #         LW_bond_messages[LW].append(message)

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

            if datapoint:
                # best annotation considering only hydrogen bonds
                if len(LW_bond_counter) > 0:
                    # sort interactions to find the best hydrogen bonds
                    # sort by badness of second worst hydrogen bond
                    LW_bond_rank = sorted(LW_bond_counter, key=lambda x : (x[3]))
                    LW = LW_bond_rank[0][0]   # best LW category

                    if LW_bond_rank[0][3] < 2.0:
                        # if second_badness is not horrible
                        #print("\nhttp://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt1.unit_id(),nt2.unit_id()))
                        datapoint['hbond_best_pair'] = LW
                        #datapoint['hbond'] = LW_bonds[LW]
                        #datapoint['hbond_messages'] = LW_bond_messages[LW]

    # use hydrogen bond data to update the interactions in match, maybe change from true to near
    still_match = []
    demotion = False
    for m in range(0,len(match)):
        # if not already near, if not enough hydrogen bonds, and if not a basepair subcategory

        # if match[m][0] == 'cWWa':
        #     for key, value in datapoint.items():
        #         if type(value) == dict:
        #             print(key)
        #             for key2, value2 in value.items():
        #                 print("   ",key2,value2)
        #         else:
        #             print(key,value)


        if not match[m][0].startswith("n") and not match[m][0] in hbond_interactions:
            print("  %-5s %-22s %-22s  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s  demoted to near by hydrogen bonds" % (match[m][0],nt1.unit_id(),nt2.unit_id(),nt1.unit_id(),nt2.unit_id()))
            # if len(match) > 1:
            #     print('  All current matches', match)
            # if datapoint:
            #     print_dictionary(datapoint)
            #     print()
            match[m][0] = "n" + match[m][0]
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
            dist2 = 10.0
            if LW in LW_to_distances:
                if len(LW_to_distances[LW]) > 1:
                    dist2 = sorted(LW_to_distances[LW])[1]
                elif len(LW_to_distances[LW]) == 1:
                    dist2 = LW_to_distances[LW][0]
            if dist2 < 5.0:
                # mark the interaction as near
                # be ready to rank both by cutoff_distance and dist2
                near_matches.append(("n"+LW,subcategory,cutoff_distance,dist2))
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

        #print(quality)

        if datapoint:
            datapoint['basepair'] = LW
            datapoint['basepair_subcategory'] = match[0][1]
            #datapoint['hbond'] = LW_bonds[interaction]
            #datapoint['hbond_messages'] = LW_bond_messages[interaction]
        return LW, subcategory, quality, datapoint
    else:
        # multiple matching basepair interactions between these two nucleotides
        LW_remaining = set([i.replace("n","") for i,s,cd in match])
        if len(LW_remaining) == 1:
            # one family, mutiple subcategories, quite OK, they are designed to overlap
            #print("  Family %s, multiple subcategories, using %d" % (match[0][0],match[0][1]))
            interaction,subcategory,cutoff_distance = match[0]
            if cutoff_distance > 0 and not "n" in interaction:
                LW = "n" + interaction
            else:
                LW = interaction

            quality = store_basepair_quality(cutoff_distance,pair_data,LW,hydrogen_bonds)

            if datapoint:
                datapoint['basepair'] = LW
                datapoint['basepair_subcategory'] = subcategory
                #datapoint['hbond'] = LW_bonds[LW]
                #datapoint['hbond_messages'] = LW_bond_messages[LW]
            return LW, subcategory, quality, datapoint

        else:
            print("  http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (nt1.unit_id(),nt2.unit_id()))
            print("  Multiple annotations meet all cutoffs, %s" % LW_remaining)
            # loop over hydrogen bond sets from best to worst
            for LW,bond_counter,checked_counter,max_badness in LW_bond_rank:
                for LW2,subcategory,cutoff_distance in match:
                    if LW == LW2:
                        quality = store_basepair_quality(cutoff_distance,pair_data,LW,hydrogen_bonds)
                        if datapoint:
                            datapoint['basepair'] = LW
                            datapoint['basepair_subcategory'] = subcategory
                            #datapoint['hbond'] = LW_bonds[LW]
                            #datapoint['hbond_messages'] = LW_bond_messages[LW]
                        print("  Using %s\n" % LW)
                        return LW, subcategory, quality, datapoint

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

    print(np.dot(np.dot(b.T,rotation),b)[0,0])
    print(np.dot(b.T,b)[0,0])
    """

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
    """

    if len(unit_data_path) > 0:

        nucleotides = structure.residues(type = ["RNA linking","DNA linking"])
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

            #all_nts[id].append((nt.index,nt.unit_id(),nt.centers["glycosidic"],nt.rotation_matrix))
            all_nts[id].append((symmetry,nt.index,nt))

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

            print("  Wrote unit data file %s" % filename)


def write_txt_output_file(outputNAPairwiseInteractions,file_id,interaction_to_list_of_tuples,categories,category_to_interactions):
    """
    Write interactions according to category, and within each
    category, write by annotation.
    Other than that, the interactions are listed in no particular order.
    """

    if "near" in categories:
        true_near = ["","n"]
    else:
        true_near = [""]

    # loop over types of output files requested
    for category in categories:
        if category == "near":
            continue
        filename = os.path.join(outputNAPairwiseInteractions,file_id + "_" + category + ".txt")
        with open(filename,'w') as f:
            # loop over all interactions found in this category

            for interaction in sorted(category_to_interactions[category]):
                if category == 'basepair':
                    # capitalize base edges to simplify
                    inter = interaction.replace("w","W").replace("s","S").replace("h","H")
                else:
                    inter = interaction

                # if this category has a restricted list of interactions to output
                if len(categories[category]) == 0 or inter in categories[category] or ("near" in categories and "n" in interaction and inter.replace('n','') in categories[category]):
                    for a,b,c in interaction_to_list_of_tuples[interaction]:
                        f.write("%s\t%s\t%s\t%s\n" % (a,inter,b,c))

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
            if category == 'basepair':
                # capitalize base edges
                inter = interaction.replace("w","W").replace("s","S").replace("h","H")
            # if this category has a restricted list of interactions to output
            if len(categories[category]) == 0 or inter in categories[category]:
                for a,b,c in interaction_to_list_of_tuples[tn+interaction]:
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

        #print(json.dumps(output))

        with open(filename,'w') as f:
            f.write(json.dumps(output))


#=======================================================================
def generatePairwiseAnnotation(entry_id, chain_id, inputPath, outputNAPairwiseInteractions, category, output_format):

    if isinstance(entry_id,str):
        entry_id = [entry_id]

    # dictionary to control what specific annotations are output, in a file named for the key
    # empty list means to output all interactions in that category
    # non-empty list specifies which interactions to output in that category
    categories = {}

    Leontis_Westhof_basepairs = ['cWW', 'cSS', 'cHH', 'cHS', 'cHW', 'cSH', 'cSW', 'cWH', 'cWS', 'tSS', 'tHH', 'tHS', 'tHW', 'tSH', 'tSW', 'tWH', 'tWS', 'tWW']

    if category:
        for category in category.split(","):
            categories[category] = []
    else:
        # default is to annotate and write just "true" basepairs
        categories['basepair'] = Leontis_Westhof_basepairs

    if 'basepair_detail' in categories:
        categories['basepair'] = Leontis_Westhof_basepairs + ['cWB','cBW']  # bifurcated pairs
    elif 'basepair' in categories:
        categories['basepair'] = Leontis_Westhof_basepairs
    else:
        categories['basepair'] = ['cWW']  # only annotate cWW, for crossing number

    # check existence of input path
    if len(inputPath) > 0 and not os.path.exists(inputPath):
        print("  Attempting to create input path %s" % inputPath)
        os.mkdir(inputPath)

    # check existence of output path
    if len(outputNAPairwiseInteractions) > 0 and not os.path.exists(outputNAPairwiseInteractions):
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
            print("  Chain argument can only be used with a single PDB file")
            PDBs = []
        else:
            chains = chain_id.split(",")
    else:
        chains = []

    # restrict dictionary of cutoffs to just the basepairs needed here
    focused_basepair_cutoffs = focus_basepair_cutoffs(nt_nt_cutoffs,categories['basepair'])
    ideal_hydrogen_bonds = load_ideal_basepair_hydrogen_bonds()

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

        print("  Reading file %s, which is number %d out of %d" % (filename, counter, len(PDBs)))
        timerData = myTimer("Reading CIF files",timerData)

        # suppress error messages, but report failures at the end
        structure, messages = load_structure(filename,file_id)

        if not structure:
            for message in messages:
                failed_structures.append((file_id,message))
            continue

        interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_in_structure(structure,categories,focused_basepair_cutoffs,ideal_hydrogen_bonds,chains,timerData)
        timerData = myTimer("Recording interactions",timerData)
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
            print('  Output format %s not recognized' % output_format)

    myTimer("summary",timerData)

    if len(failed_structures) > 0:
        print("  Error messages:")
        for message in failed_structures:
            print("  %s %s" % message)
    else:
        print("All files read successfully")

    # note status of stacking annotations
    if 'stacking' in categories:
        print("  Stacking annotations are not yet finalized")

    if 'basepair' in categories:
        print("  Basepair annotations are not yet finalized")

if __name__=="__main__":

    # allow user to specify input and output paths
    parser = argparse.ArgumentParser()
    parser.add_argument('PDBfiles', type=str, nargs='+', help='.cif filename(s)')
    parser.add_argument('-o', "--output", help="Output Location of Pairwise Interactions")
    parser.add_argument('-i', "--input", help='Input Path')
    parser.add_argument('-c', "--category", help='Interaction category or categories (basepair,stacking,sO,backbone,coplanar,basepair_detail,covalent,sugar_ribose,near)')
    parser.add_argument('-f', "--format", help='Output format (txt,ebi_json)')
    parser.add_argument("--chain", help='Chain or chains separated by commas, no spaces; only for one PDB file')

    problem = False
    args = parser.parse_args()

    # Following if statements deal with command line arguments.
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

    entry_id = args.PDBfiles

    generatePairwiseAnnotation(entry_id, chain_id, inputPath, outputNAPairwiseInteractions, category, outputFormat)

