# Read file called atom_mappings.txt, create dictionary from this where keys are modified
# nucleotides and values are a triple of parent sequence, parent atom, then corresponding
# atom from modified nucleotide.

# NOTES:
    # How the following dictionaries work

    # modified_atom_to_parent['4EN']['N8'] = 'C8': 4EN is modified, N8 is from the modified 4EN, corresponds to parent A's C8
    # parent_atom_to_modified['4EN']['C8'] = 'N8': 4EN is modified, C8 is from parent A, corresponds to modified 4EN's N8

    # modified_base_to_parent['PSU'] = 'U': PSU key yields parent U
    # modified_base_to_parent['4EN'] = 'A': 4EN key yields parent A
    # modified_base_to_hydrogens['PSU']: list of hydrogens on the base of PSU
    # modified_base_to_hydrogen_coordinates['PSU']['HN1']: triple of coordinates of H5, the hydrogen of U that PSU HN1 is mapped to
    # modified_base_atom_list['PSU']: list of names of all atoms in PSU

from fr3d import definitions as defs
import os
import sys

if sys.version_info[0] < 3:
    read_mode = 'rb'
else:
    read_mode = 'rt'

def create_modified_nucleotide_to_parent_mappings():
    # Read in mapping file from the folder where this program is installed
    current_path,current_program = os.path.split(os.path.abspath(__file__))
    filename = os.path.join(current_path,"atom_mappings.txt")

    #print('mapping.py is being run in path %s' % current_path)
    #print('mapping.py is trying to open %s' % filename)

    with open(filename, read_mode) as fid:
        lines = fid.readlines()

    modified_atom_map = {}

    for line in lines:
        fields = line.split()
        if len(fields) == 4:
            if len(fields[1]) > 0 and len(fields[3]) > 0:
                # only process lines that list a parent atom and a modified atom
                if not fields[2] in modified_atom_map:
                    modified_atom_map[fields[2]] = []
                modified_atom_map[fields[2]].append((fields[0], fields[1], fields[3]))

    modified_base_to_parent = {}
    modified_atom_to_parent = {}
    parent_atom_to_modified = {}
    modified_base_atom_list = {}
    modified_base_to_hydrogens = {}
    modified_base_to_hydrogen_coordinates = {}

    for modified_nucleotide in modified_atom_map:
        modified_base_to_parent[modified_nucleotide] = modified_atom_map[modified_nucleotide][0][0]

        modified_base_atom_list[modified_nucleotide] = []
        modified_atom_to_parent[modified_nucleotide] = {}
        parent_atom_to_modified[modified_nucleotide] = {}
        modified_base_to_hydrogens[modified_nucleotide] = []
        modified_base_to_hydrogen_coordinates[modified_nucleotide] = {}

        for fields in modified_atom_map[modified_nucleotide]:
            if len(fields) == 3:
                parent_nucleotide,parent_atom,modified_atom = fields
                modified_atom_to_parent[modified_nucleotide][modified_atom] = parent_atom
                parent_atom_to_modified[modified_nucleotide][parent_atom] = modified_atom
                if parent_atom in defs.NAbaseheavyatoms[parent_nucleotide] or parent_atom in defs.NAbasehydrogens[parent_nucleotide]: # The parent mapping is in the base
                    modified_base_atom_list[modified_nucleotide].append(modified_atom)
                    if modified_atom[0] == 'H':
                        modified_base_to_hydrogens[modified_nucleotide].append(modified_atom)
                        modified_base_to_hydrogen_coordinates[modified_nucleotide][modified_atom] = (defs.NAbasecoordinates[parent_nucleotide][parent_atom])

    return modified_base_to_hydrogens, modified_atom_to_parent, parent_atom_to_modified, modified_base_to_parent, modified_base_atom_list,  modified_base_to_hydrogen_coordinates

try:
    modified_base_to_hydrogens, modified_atom_to_parent, parent_atom_to_modified, modified_base_to_parent, modified_base_atom_list,  modified_base_to_hydrogen_coordinates = create_modified_nucleotide_to_parent_mappings()
    # print("Modified nucleotide mappings read successfully.")

    all_parents = set()
    for modified in modified_base_to_parent.keys():
        parent = modified_base_to_parent[modified]
        all_parents.add(parent)
    # print("All parents: %s" % sorted(all_parents))

    # modified = 'OMG'
    # print(modified_base_to_hydrogens[modified])
    # print(modified_base_to_hydrogen_coordinates[modified])

except Exception as e:
    print("mapping.py is unable to load mappings for modified nucleotides.")
    print("This can happen after installing with 'python setup.py install' with no known fix.")
    print("Instead, from the directory where setup.py is, use 'python -m pip install .'")
    print('Error message: %s' % str(e))
