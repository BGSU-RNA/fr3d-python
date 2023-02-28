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
from traceback import print_exception
from fr3d import definitions as defs
import os 
import sys
import pickle
from fr3d.localpath import fr3d_pickle_path

#Deals with opening csv file
if sys.version_info[0] < 3:
    from io import open as open
    
# See if pickle files already exist that contain mapped dictionaries
paths = []
if fr3d_pickle_path: 
    if not os.path.exists(fr3d_pickle_path + '\parent_atom_to_modified.pickle'):
        parent_atom_to_modified_pickle_path = fr3d_pickle_path + '\parent_atom_to_modified.pickle'
        paths.append(parent_atom_to_modified_pickle_path)
    if not os.path.exists(fr3d_pickle_path + '\modified_atom_to_parent.pickle'):
        modified_atom_to_parent_pickle_path = fr3d_pickle_path + '\modified_atom_to_parent.pickle'
        paths.append(modified_atom_to_parent_pickle_path)
    if not os.path.exists(fr3d_pickle_path + '\modified_base_atom_list.pickle'):
        modified_base_atom_list_pickle_path = fr3d_pickle_path + '\modified_base_atom_list.pickle'
        paths.append(modified_base_atom_list_pickle_path)    
    if not os.path.exists(fr3d_pickle_path + '\modified_base_to_parent.pickle'):
        modified_base_to_parent_pickle_path = fr3d_pickle_path + '\modified_base_to_parent.pickle'
        paths.append(modified_base_to_parent_pickle_path)
    if not os.path.exists(fr3d_pickle_path + '\modified_base_to_hydrogens.pickle'):
        modified_base_to_hydrogens_pickle_path = fr3d_pickle_path + '\modified_base_to_hydrogens.pickle'
        paths.append(modified_base_to_hydrogens_pickle_path)       
    if not os.path.exists(fr3d_pickle_path + '\modified_base_to_hydrogens_coordinates.pickle'):
        modified_base_to_hydrogens_coordinates_pickle_path = fr3d_pickle_path + '\modified_base_to_hydrogens_coordinates.pickle'
        paths.append(modified_base_to_hydrogens_coordinates_pickle_path)   

# If no pickle file already exist, then read the text file and create the mapping dictionaries
if len(paths) > 0:  
    # Read in mapping file wherever its located on system in python path. Read file line by line and create mapping.
    modified_atom_map = {}
    path =  os.path.dirname(os.path.abspath(__file__))
    subject = csv.reader(open(os.path.join(path, "atom_mappings_refined.txt"), "r", encoding="utf8"), delimiter="\t")
    for line in subject:
        lastline = ""
        if len(line) > 3:
            if line[2] not in modified_atom_map.keys():
                modified_atom_map[line[2]] = []
            modified_atom_map[line[2]].append((line[0], line[1], line[3]))
                                            #parent,    parentAtom, mapped modified atom
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
    # Dump those dictionaries to pickle files now so we don't have to process the file again
    print("Writing serialized files of modified nucleotide mappings...")
    for file in paths:
        with open(file, 'wb') as handle:
            name = file.replace(fr3d_pickle_path + '\\', '')
            end_name = name.replace('.pickle', '')
            pickle.dump(end_name, handle, protocol=pickle.HIGHEST_PROTOCOL)
        print("Serialized modified mapping file written at: " + file)

# if the pickles exist, load the pickles to fill data
elif len(paths) == 0: 
    print("Serialized files for mappings found. Loading serialized mappings.")
    if fr3d_pickle_path:
        if os.path.exists(fr3d_pickle_path + '\parent_atom_to_modified.pickle'):
            parent_atom_to_modified_pickle_path = fr3d_pickle_path + '\parent_atom_to_modified.pickle'
            with open(parent_atom_to_modified_pickle_path, 'rb') as handle:
                parent_atom_to_modified = pickle.load(handle)
        if os.path.exists(fr3d_pickle_path + '\modified_atom_to_parent.pickle'):
            modified_atom_to_parent_pickle_path = fr3d_pickle_path + '\modified_atom_to_parent.pickle'
            with open(modified_atom_to_parent_pickle_path, 'rb') as handle:
                modified_atom_to_parent = pickle.load(handle)
        if os.path.exists(fr3d_pickle_path + '\modified_base_atom_list.pickle'):
            modified_base_atom_list_pickle_path = fr3d_pickle_path + '\modified_base_atom_list.pickle'
            with open(modified_base_atom_list_pickle_path, 'rb') as handle:
                modified_base_atom_list = pickle.load(handle)
        if os.path.exists(fr3d_pickle_path + '\modified_base_to_parent.pickle'):
            modified_base_to_parent_pickle_path = fr3d_pickle_path + '\modified_base_to_parent.pickle'
            with open(modified_base_to_parent_pickle_path, 'rb') as handle:
                modified_base_to_parent = pickle.load(handle)
        if os.path.exists(fr3d_pickle_path + '\modified_base_to_hydrogens.pickle'):
            modified_base_to_hydrogens_pickle_path = fr3d_pickle_path + '\modified_base_to_hydrogens.pickle'
            with open(modified_base_to_hydrogens_pickle_path, 'rb') as handle:
                modified_base_to_hydrogens = pickle.load(handle)
        if os.path.exists(fr3d_pickle_path + '\modified_base_to_hydrogens_coordinates.pickle'):
            modified_base_to_hydrogens_coordinates_pickle_path = fr3d_pickle_path + '\modified_base_to_hydrogens_coordinates.pickle'
            with open(modified_base_to_hydrogens_coordinates_pickle_path, 'rb') as handle:
                modified_base_to_hydrogens_coordinates = pickle.load(handle)
    print("Serialized mappings for modified nucleotides loaded.")
