# Read file called atom_mappings.txt, create dictionary from this where keys are modified
# nucleotides and values are the parent nucleotide

# from fr3d import definitions as defs
import os
import sys

if sys.version_info[0] < 3:
    read_mode = 'rb'
else:
    read_mode = 'rt'

def modified_nucleotide_to_parent_mapping():
    # Read in mapping file from the folder where this program is installed
    current_path,current_program = os.path.split(os.path.abspath(__file__))
    filename = os.path.join(current_path,"atom_mappings.txt")

    #print('mapping.py is being run in path %s' % current_path)
    #print('mapping.py is trying to open %s' % filename)

    with open(filename, read_mode) as fid:
        lines = fid.readlines()

    modified_base_to_parent = {}

    for line in lines:
        fields = line.split()
        if len(fields) == 4:
            if len(fields[1]) > 0 and len(fields[3]) > 0:
                # only process lines that list a parent atom and a modified atom
                if not fields[2] in modified_base_to_parent:
                    modified_base_to_parent[fields[2]] = fields[0]

    return modified_base_to_parent

try:
    modified_base_to_parent = modified_nucleotide_to_parent_mapping()
    # print("Modified nucleotide mappings read successfully.")

    # modified = 'OMG'
    # print(modified_base_to_hydrogens[modified])
    # print(modified_base_to_hydrogen_coordinates[modified])

except Exception as e:
    print("mapping_nucleotide_to_parent.py is unable to load mappings for modified nucleotides.")
    print("This can happen after installing with 'python setup.py install' with no known fix.")
    print("Instead, from the directory where setup.py is, use 'python -m pip install .'")
    print('Error message: %s' % str(e))
