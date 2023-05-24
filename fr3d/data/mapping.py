# Read file called atom_mappings.txt, create dictionary from this where keys are modified 
# nucleotides and values are a triple of parent sequence, parent atom, then corresponding 
# atom from modified nucleotide. 

# 02-15-2023 update atom_mappings.txt to atom_mappings_refined.txt. Improvements on mappings for hydrogens

# NOTES: 
    # How the following dictionaries work

    # modified_atom_to_parent['4EN']['N8'] = 'C8': 4EN is modified, N8 is from the modified 4EN, corresponds to parent A's C8
    # parent_atom_to_modified['4EN']['C8'] = 'N8': 4EN is modified, C8 is from parent A, corresponds to modified 4EN's N8

    # modified_base_to_parent['PSU'] = 'U': PSU key yields parent U
    # modified_base_to_parent['4EN'] = 'A': 4EN key yields parent A
    # modified_base_to_hydrogens['PSU']: list of hydrogens on the base of PSU
    # modified_base_to_hydrogens_coordinates['PSU']['HN1']: triple of coordinates of H5, the hydrogen of U that PSU HN1 is mapped to
    # modified_base_atom_list['PSU']: list of names of all atoms in PSU
import csv
from fr3d import definitions as defs
import os 
import sys

#from fr3d.data.atom_mappings_refined import mapping_text

#from traceback import print_exception
# if sys.version_info[0] < 3:
#     #Deals with opening csv file
#     from io import open as open

if sys.version_info[0] < 3:
    read_mode = 'rb'
else:
    read_mode = 'rt'

def create_modified_nucleotide_to_parent_mappings():
    # Read in mapping file wherever its located on system in python path. Read file line by line and create mapping.
    modified_atom_map = {}

    # the next line caused problems on a user's system
    #path =  os.path.dirname(os.path.abspath(__file__))

    # this works better in the hydrogen_bonds.py program, maybe it will work here:
    current_path,current_program = os.path.split(os.path.abspath(__file__))

    print('mapping.py is being run in path %s' % current_path)

    filename = os.path.join(current_path,"atom_mappings_refined.txt")

    print('mapping.py is trying to open %s' % filename)

    with open(filename, read_mode) as fid:
        lines = fid.readlines()

    for line in lines:
        fields = line.split()
        if len(fields) == 4:
            if not fields[2] in modified_atom_map:
                modified_atom_map[fields[2]] = []
            modified_atom_map[fields[2]].append((fields[0], fields[1], fields[3]))

    # subject = csv.reader(open(os.path.join(path, "atom_mappings_refined.txt"), "r", encoding="utf8"), delimiter="\t")
    # for line in subject:
    #     lastline = ""
    #     if len(line) > 3:
    #         if line[2] not in modified_atom_map.keys():
    #             modified_atom_map[line[2]] = []
    #         modified_atom_map[line[2]].append((line[0], line[1], line[3]))
                                            #parent,    parentAtom, mapped modified atom

    """
    for line in mapping_text.split("\n"):
        fields = line.split()
        if len(fields) == 4:
            if not fields[2] in modified_atom_map:
                modified_atom_map[fields[2]] = []
            modified_atom_map[fields[2]].append((fields[0], fields[1], fields[3]))
    """

    modified_base_to_hydrogens = {}
    modified_atom_to_parent = {}
    parent_atom_to_modified = {}
    modified_base_to_parent = {}
    modified_base_atom_list = {} 
    modified_base_to_hydrogens_coordinates = {}

    for modified_nucleotide in modified_atom_map:
        modified_base_to_hydrogens[modified_nucleotide] = []
        modified_base_to_parent[modified_nucleotide] = {}
        modified_base_to_parent[modified_nucleotide] = modified_atom_map[modified_nucleotide][0][0]
        modified_base_to_hydrogens_coordinates[modified_nucleotide] = {}
        modified_base_atom_list[modified_nucleotide] = []
        modified_atom_to_parent[modified_nucleotide] = {}
        parent_atom_to_modified[modified_nucleotide] = {}

        for atom in modified_atom_map[modified_nucleotide]:
            if len(atom) == 3:
                modified_atom_to_parent[modified_nucleotide][atom[2]] = atom[1]
                parent_atom_to_modified[modified_nucleotide][atom[1]] = atom[2]
                if atom[1] in defs.NAbaseheavyatoms[atom[0]] or atom[1] in defs.NAbasehydrogens[atom[0]]: # The parent mapping is in the base
                    modified_base_atom_list[modified_nucleotide].append(atom[2])
                    modified_base_to_hydrogens[modified_nucleotide].append(atom[2])
                    modified_base_to_hydrogens_coordinates[modified_nucleotide][atom[2]] = (defs.NAbasecoordinates[atom[0]][atom[1]])

    return modified_base_to_hydrogens, modified_atom_to_parent, parent_atom_to_modified, modified_base_to_parent, modified_base_atom_list,  modified_base_to_hydrogens_coordinates


modified_base_to_hydrogens, modified_atom_to_parent, parent_atom_to_modified, modified_base_to_parent, modified_base_atom_list,  modified_base_to_hydrogens_coordinates = create_modified_nucleotide_to_parent_mappings()


try:
    pass
    # print("Modified nucleotide mappings read successfully.")
except Exception as e:
    print("mapping.py is unable to load mappings for modified nucleotides.")
    print('Error message: %s' % str(e))
