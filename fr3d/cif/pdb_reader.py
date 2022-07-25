from Bio.PDB import PDBParser
from Bio.PDB.PDBParser import PDBParser

# import necessary fr3d classes
from fr3d.data import Atom
from fr3d.data import Component
from fr3d.data import Structure

import sys
import operator as op
import itertools as it
import collections as coll
import copy
import argparse
import os


##### LOADING THE FILE #############################################
# This will be done in NA_pairwise_interactions.py in the future ###
# but for testing I'm using it here for now ########################
# add argument to specify pdb file #################################
####################################################################
parser = argparse.ArgumentParser()
parser.add_argument('PDBfiles', type=str, nargs='+', help='.pdb filename(s)')
args = parser.parse_args()
PDBs = []  # list of (path,filename) entries
entries = args.PDBfiles
print(entries)

for entry in entries:
    path_split = os.path.split(entry)

    if len(path_split[0]) > 0:
        PDBs.append(path_split)
    # else:
    #     PDBs.append((inputPath,entry))
print(PDBs)

for entry in entries:
    name = entry.split(".")[0]
print(name)

###################################################################
### Read in the specified PDB file using BioPython's PDB Reader ###
# reading in PDB file #############################################
###################################################################

p = PDBParser(PERMISSIVE=1)
#structure = p.get_structure("Top_Model_11nt_GAAA", "Top Model_11nt_GAAA.pdb")
                            # What the file will be referred to as, file name (must be in root dir)
structure = p.get_structure(name, entries[0])

###################################################################
# structure is the object that contains all relevent information ##
# about the specified 3D structure ################################
###################################################################

###################################################################
# Extract relevant information ####################################
# Header contains info given in the remarks of a pdb file. ########
# IDCode refers to 4 letter PDB Code ##############################
###################################################################
if not structure.header['name']:
    structure.header['name'] = "UNNAMED"
    #name = "UNNAMED"
    #name = "TOP_Model_11nt_GAAA"
else: 
    name = structure.header['idcode']


def __generate_atoms__(structure):
    """Method that loops through all the atoms in a pdb file and then creates
    atom objects defined by fr3d-python's Atom object.

    Uses BioPython Parser's structure to loop through models in a structure, 
    residues in a model, and then atoms in a residue to get all atoms.
    Extracts relevant information about each atom and then appends it to a list.

    Returns a list of atom objects that represent each atom in a structure.  
    """
    atoms = set()     
    for model in structure:
        residues = model.get_residues()
        for residue in residues:
            full_id = residue.get_full_id()
            ins_code = full_id[3][2] 
            pdb=full_id[0]
            this_model = int(full_id[1]) + 1 # BioPython starts at 0 and Fr3d-Python starts at 1. Add 1 to each model so unit ids match 
            this_chain = full_id[2]
            component_number = full_id[3][1]
            if 'H' in full_id[3][0][0]:
                res_group = 'HETATM'
            else:
                res_group = 'ATOM'
            res = residue.get_resname(),
            res=res[0]
            if ins_code == " ":
                ins_code = None

            for atom in residue:
                first = atom.id[0]
                # logic to extract the type of atom from the id
                if 'C' == first: #Carbon
                    atom_type = 'C' 
                elif 'O' == first: #Ox
                    atom_type = 'O'
                elif 'P' == first: #Make sure its Phosphorus and not oxygen
                    atom_type = 'P'
                elif 'N' == first: # nitrogen
                    atom_type = 'N'
                else: #Magnesium, other ions
                    atom_type = atom.id

                x = atom.coord[0]
                y = atom.coord[1]
                z = atom.coord[2]

                alt_id = atom.get_altloc()
                if alt_id == " ":
                    alt_id = None
                atoms.append(Atom(x=x, y=y, z=z,
                    pdb=name,
                    model=this_model,
                    chain=this_chain,
                    component_id=res,
                    component_number=component_number,
                    component_index=component_number,
                    insertion_code=ins_code,
                    alt_id= alt_id,
                    group=res_group,
                    type=atom_type,
                    name=atom.get_name(),
                    symmetry=None,
                    polymeric=False))
    return atoms

                        # group=atom['group_PDB'],
                        # type=atom['type_symbol'],
# Create dictionaries with structure header info
# for model in structure:
#     for chain in model:
#         for residue in chain:
#             for atom in residue:
#                 print(atom.__dict__)



def __residues__(structure):
    key = op.attrgetter(
                'pdb',
                'model',
                'chain',
                'component_id',
                'component_number',
                'insertion_code', #in python 3, the sorted function below cannot accept None values, which sometimes there are in insertion code
                'symmetry'
                )
    mapping = it.groupby(sorted(__generate_atoms__(structure), key=key), key)
    for comp_id, all_atoms in mapping:
        print(comp_id, all_atoms)
        for atoms in __group_alt_atoms__(all_atoms):
            first = atoms[0]
            #atom_type = self._chem.get(first.component_id, {})
            #atom_type = atom_type.get('type', None)
            alt_id = first.alt_id
            if alt_id == '.':
                alt_id = None
            print(first.component_id, atoms, first.component_number)
            yield Component(
                atoms,
                pdb=first.pdb,
                model=first.model,
                type=None, #atom_type,
                alt_id=alt_id,
                chain=first.chain,
                symmetry=first.symmetry,
                sequence=first.component_id,
                number=first.component_number,
                index=first.component_index,
                insertion_code=first.insertion_code,
                polymeric=first.polymeric,
            )

def __group_alt_atoms__(atoms):
        def ordering_key(atoms):
            return atoms[0].alt_id
        alt_ids = coll.defaultdict(list)
        for atom in atoms:
            alt_ids[atom.alt_id].append(atom)

        if len(alt_ids) == 1:
            return list(alt_ids.values())

        if None in alt_ids:
            common = alt_ids.pop(None)
            for alt_id, specific_atoms in list(alt_ids.items()):
                for common_atom in common:
                    copied = copy.deepcopy(common_atom)
                    copied.alt_id = alt_id
                    specific_atoms.append(copied)

        return sorted(list(alt_ids.values()), key=ordering_key)

#__group_alt_atoms__(__generate_atoms__(structure))

print(__residues__(structure))
for iter in __residues__(structure):
    print(iter)
# print(atoms)

                #print(z,y,z,name,models,chains)
                # atoms.append(Atom(x=atom.coord[0], y=atom.coords['y'], z=atoms.coord['z'],
                #     pdb=self.pdb,
                #     model=self.model,
                #     chain=self.chain,
                #     component_id=self.component_id,
                #     component_number=self.component_number,
                #     component_index=self.component_index,
                #     insertion_code=self.insertion_code,
                #     alt_id=self.alt_id,
                #     group=self.group,
                #     type=self.type,
                #     name=self.name,
                #     symmetry=self.symmetry,
                #     polymeric=self.polymeric))


# new = Cif(structure)
new = Component(structure, __generate_atoms__(structure))

def structure(self):
    """Get the structure from the Cif file.
    :returns: The first structure in the cif file.
    """
    if sys.version_info[0] < 3:
        pdb = self.data.getName()
    else: 
        pdb = self.data.name
    residues = __residues__(structure)
    return Structure(list(residues), pdb=name)