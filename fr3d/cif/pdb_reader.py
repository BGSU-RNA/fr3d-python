from importlib.resources import path
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
import logging

###################################################################
# PDBStructure Class ##############################################
# See if I can extend Cif to do this ##############################
###################################################################
class PDBStructure(object):
    """Container for data extracted from PDB Files. Uses PDBParser 
    by BioPython to extract relevant information from PDB files. """
    def __init__(self, filename):
        p = PDBParser(PERMISSIVE=1) # Call BioPython Method to read cif file
        self.structure = p.get_structure(filename, filename) # biopython method (What the file will be referred to as, what the file is named in your local path)
        self.name = self.structure.header['idcode']
        if self.name == " " or not self.name:
            path = filename.split('.')
            names = path[0].split("\\") #slightly sloppy but it works
            self.name = names[-1]
        #structure = p.get_structure("Top_Model_11nt_GAAA", "Top Model_11nt_GAAA.pdb")
        self.logger = logging.getLogger('fr3d.cif.reader.PDBStructure')
        self.residues = self.__residues__(self.name)

    def __generate_atoms__(self, pdb):
        """Method that loops through all the atoms in a pdb file and then creates
        atom objects defined by fr3d-python's Atom object.

        Uses BioPython Parser's structure to loop through models in a structure, 
        residues in a model, and then atoms in a residue to get all atoms.
        Extracts relevant information about each atom and then appends it to a list.

        Returns a list of atom objects that represent each atom in a structure.  
        """
        atoms = [] # Maybe this can be a set 
        for model in self.structure:
            residues = model.get_residues() # Biopython 
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
                        pdb=self.name,
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
                        symmetry='1_555', #I haven't figured out how to extract symmetries from pdb files yet. Resort to identity
                        polymeric=True)) # Need to find a way to parse this from biopython. Important, may be relevent in structures.py
        return atoms
        
    # TODO: See if we can import the same method in the Cif class of reader.py here.
    def __residues__(self, pdb):
        key = op.attrgetter(
                    'pdb',
                    'model',
                    'chain',
                    'component_id',
                    'component_number',
                    'insertion_code', #in python 3, the sorted function below cannot accept None values, which sometimes there are in insertion code
                    'symmetry'
                    ) 
        mapping = it.groupby(sorted(self.__generate_atoms__(self.name), key=key), key)
        for comp_id, all_atoms in mapping:
            for atoms in self.__group_alt_atoms__(all_atoms):
                first = atoms[0]
                #atom_type = self._chem.get(first.component_id, {})
                #atom_type = atom_type.get('type', None)
                alt_id = first.alt_id
                if alt_id == '.':
                    alt_id = None
                yield Component(
                    atoms,
                    pdb=first.pdb,
                    model=first.model,
                    type="RNA linking", #atom_type,
                    alt_id=alt_id,
                    chain=first.chain,
                    symmetry=first.symmetry,
                    sequence=first.component_id,
                    number=first.component_number,
                    index=first.component_index,
                    insertion_code=first.insertion_code,
                    polymeric=first.polymeric,
                )

    def __group_alt_atoms__(self, atoms):
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

    def structures(self):
        """Get the structure from the Cif file.
        :returns: The first structure in the cif file.
        """
        pdb = self.name
        residues = self.__residues__(pdb)
        return Structure(list(residues), pdb=pdb)