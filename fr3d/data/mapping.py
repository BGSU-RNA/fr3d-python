# Read file called atom_mappings.txt, create dictionary from this where keys are modified 
# nucleotides and values are a triple of parent sequence, parent atom, then corresponding 
# atom from modified nucleotide. 
import csv
from traceback import print_exception
from fr3d import definitions as defs
import os 
import sys

modified_atom_map = {}
path =  os.path.dirname(os.path.abspath(__file__))

subject = csv.reader(open(os.path.join(path, "atom_mappings.txt"), "r", encoding="utf8"), delimiter="\t")
for line in subject:
    lastline = ""
    if len(line) > 3:
        if line[2] not in modified_atom_map.keys():
            modified_atom_map[line[2]] = []
        modified_atom_map[line[2]].append((line[0], line[1], line[3]))

# modified_atom_to_parent['PSU']['N1'] = 'C5'
# parent_atom_to_modified['PSU']['C5'] = 'N1'
# modified_base_to_parent['PSU'] = 'U'


modified_hydrogens = {}
modified_atom_to_parent = {}
parent_atom_to_modified = {}
modified_base_to_parent = {}
modified_base_atom_list = {} 

modified_hydrogens_coordinates = {}
for modified_nucleotide in modified_atom_map:
    # print(modified_nucleotide)
    modified_hydrogens[modified_nucleotide] = []
    modified_base_to_parent[modified_nucleotide] = {}
    # print(modified_atom_map[modified_nucleotide])
    modified_base_to_parent[modified_nucleotide] = modified_atom_map[modified_nucleotide][0][0]
    modified_hydrogens_coordinates[modified_nucleotide] = {}
    modified_base_atom_list[modified_nucleotide] = []
    modified_atom_to_parent[modified_nucleotide] = {}
    parent_atom_to_modified[modified_nucleotide] = {}

    for atom in modified_atom_map[modified_nucleotide]:
        if len(atom) == 3:
            modified_atom_to_parent[modified_nucleotide][atom[1]] = atom[2]
            parent_atom_to_modified[modified_nucleotide][atom[2]] = atom[1]
            modified_base_atom_list[modified_nucleotide].append(atom[2])

            if 'H' in atom[2] and "'" not in atom[1] and "P" not in atom[1]: #atom atom of modified but not a backbone atom since modifieds can have unconventional naming conventions

                modified_hydrogens[modified_nucleotide].append(atom[2])
                modified_hydrogens_coordinates[modified_nucleotide][atom[2]] = (defs.NAbasecoordinates[atom[0]][atom[1]])
    print(parent_atom_to_modified[modified_nucleotide])
# print(modified_hydrogens)
# print(modified_hydrogens_coordinates['OMG']['HN1'])
