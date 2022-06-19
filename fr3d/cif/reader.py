import re
import itertools as it
import collections as coll
import warnings
import logging
import operator as op
import functools as ft
import copy
import sys
import numpy as np

if sys.version_info[0] < 3:
    from itertools import ifilter as filter # old name
    from pdbx.reader.PdbxReader import PdbxReader as Reader #If running python 2.7, use pdbx reader 
else:
    from pdbx import PdbxReader as Reader #if running in python 3, install mmcif package to deal with reading. ([mmcif_pdbx] will need added to setup.py for this)


from fr3d.data import Atom
from fr3d.data import Component
from fr3d.data import Structure

oldStructures = ["6I2N", "4WR6", "6H5Q", "4WRO", "4WZD", "7MSF", "1EO4", "6MSF", "6NUT", "5A79", "5A7A", "5APO", "1VS9", "2I1C", "6QNQ", "5AFI", "1FJF", "5AA0", "5MJV", "5MSF", "5Z9W", "4Z92", "5FN1", "6GV4", "5M74"]

""" The set of symbols that mark an operator expression as complex """
COMPLEX_SYMBOLS = set('()-')


class MissingBlockException(Exception):
    """This class is raised when trying to get a missing block of data.
    """
    pass


class MissingColumn(Exception):
    """This is raised when trying to get a missing column from a table.
    """
    pass


class ComplexOperatorException(Exception):
    """This is raised when we come across complex operators that we cannot
    easily deal with. These tend to show up in viral structures and not things
    we deal with currently.
    """
    pass


class UnusableUnobservedTable(Exception):
    pass


class MissingSymmetry(Exception):
    """This is raised when we cannot determine a symmetry operator for an atom.
    """
    pass


class UnmppedResidueException(Exception):
    """This is raised if we do not map all residues in a chain to the
    experimental sequence.
    """
    pass


class TooManyMappedResidueException(Exception):
    """Raised if too many units are mapped.
    """
    pass


class Cif(object):
    """Top level container for all Cif related data. This assumes that each
    mmCIF file contains a single datablock. This doesn't have to be true but
    makes things easier.
    """

    def __init__(self, handle=None, data=None):
        if data is None:
            reader = Reader(handle)
            self.data = []
            reader.read(self.data)
            self.data = self.data[0]
        else:
            self.data = data

        if handle is None and data is None:
            raise ValueError("Must give either handle or data")
        if sys.version_info[0] < 3:
            self.pdb = self.data.getName()
        else: 
            self.pdb = self.data.name
        self._operators = self.__load_operators__()
        self._assemblies = self.__load_assemblies__()
        self._entities = self.__load_entities__()
        self._chem = self.__load_chem_comp__()
        self.logger = logging.getLogger('fr3d.cif.reader.Cif')

    def __load_operators__(self):
        operators = {}
        #some Cif files don't have a struct_oper_list
        if hasattr(self, 'pdbx_struct_oper_list'):
            for op in self.pdbx_struct_oper_list:
                op['matrix'] = [[None] * 3, [None] * 3, [None] * 3]
                op['vector'] = [None] * 3

                for row in range(3):
                    op['vector'][row] = float(op['vector[%s]' % str(row + 1)])

                    for column in range(3):
                        key = 'matrix[%s][%s]' % (str(row + 1), str(column + 1))
                        op['matrix'][row][column] = float(op[key])

                transform = np.zeros((4, 4))
                transform[0:3, 0:3] = op['matrix']
                transform[0:3, 3] = op['vector']
                transform[3, 3] = 1.0

                op['matrix'] = np.array(op['matrix'])
                op['vector'] = np.array(op['vector'])
                op['transform'] = np.array(transform)

                operators[op['id']] = op

        identity = self.__identity_operator__()
        operators[identity['id']] = identity

        return operators

    def __identity_operator__(self):
        mat = np.identity(3)
        vector = np.array([1, 1, 1])
        trans = np.zeros((4, 4))
        trans[0:3, 0:3] = mat
        trans[0:3, 3] = vector
        trans[3, 3] = 1.0
        return {
            'id': 'I',
            'name': 'I',
            'vector': vector,
            'matrix': mat,
            'transform': trans
        }

    def __load_assemblies__(self):
        listOfNumbers = []
        assemblies = coll.defaultdict(list)
        if hasattr(self, 'pdbx_struct_assembly_gen'): #3% of structures don't have an assembly gen
            for assembly in self.pdbx_struct_assembly_gen:
                oper_expression = assembly['oper_expression']

                # # TODO: Implement computation of complex operators
                # if COMPLEX_SYMBOLS & set(oper_expression):
                #     warnings.warn('Cannot compute symmetries from complex '
                #                   'expressions. Will use a simple identity '
                #                   'transformation if no others possible')
                #     operators = []
                # else:
                operators = oper_expression.split(',')

                for asym_id in assembly['asym_id_list'].split(','):
                    for operator in operators:
                        listOfNumbers = []
                        #Some Cif files such as 5MSF use a whole list of operators, a hyphen should do an adequte job of flagging these situations
                        if '-' in operator and 'X' not in operator and 'P' not in operator:# list of operators, but not crystallographic operators
                            op1 = operator.replace("(", "")
                            op1 = op1.replace(")", "") 
                            startEnd = op1.split('-')
                            if sys.version_info[0] < 3:
                                for number in xrange(int(startEnd[0]), int(startEnd[1])+1):
                                    listOfNumbers.append(number)                                
                            else:
                                for number in range(int(startEnd[0]), int(startEnd[1])+1):
                                    listOfNumbers.append(number)
                            for symmetry in listOfNumbers: #Add all of these symmetry operators and then go to the next operator
                                op = self._operators[str(symmetry)]    
                                previous_id = [x['id'] for x in assemblies[asym_id]] #avoid applying the same operator twice.
                                if not op['id'] in previous_id:
                                    assemblies[asym_id].append(op)
                            continue
                        #Others have a crystal frame transformation
                        elif 'X' in operator:# or 'P' in operator: # P moves coordinates into a "standard" icosahedral point symmettry frame
                            #X# moves x,y,z into crystallographic positions.
                            #op = self._operators['I'] #Just apply Identity
                            pass # I don't think we need to apply anything.
                        elif 'P' in operator:  # P moves coordinates into a "standard" icosahedral point symmettry frame
                            if self.pdb in oldStructures:
                                op = self.operators[operator]  # For our database, unit id needs to be the same as it used to be so this is for backward compatibility with unit ids created for and used by the BGSU database. 
                                                               # This is only applied for the structures in the old structures list.
                        
                        else: #Normal case
                            if '(' in operator or ')' in operator:
                                operator = operator.replace("(","")
                                operator = operator.replace(")","")
                            op = self._operators[operator]

                        previous_id = [x['id'] for x in assemblies[asym_id]] #avoid applying the same operator twice.
                        if not op['id'] in previous_id:
                            assemblies[asym_id].append(op)
        return assemblies

    def __load_entities__(self):
        entities = {}
        if hasattr(self, 'entity'):
            for entity in self.entity:
                entities[entity['id']] = entity
        return entities

    def __load_chem_comp__(self):
        chem = {}
        if hasattr(self, 'chem_comp'):
            for obj in self.chem_comp:
                chem[obj['id']] = obj
        return chem

    def structure(self):
        """Get the structure from the Cif file.
        :returns: The first structure in the cif file.
        """
        if sys.version_info[0] < 3:
            pdb = self.data.getName()
        else: 
            pdb = self.data.name
        residues = self.__residues__(pdb)
        return Structure(list(residues), pdb=pdb)

    def experimental_sequence(self, chain):
        """Get the experimental sequence for a given chain.
        :chain: The chain name to use, should be the pdb_strand_id in the cif
        file.
        :returns: A list of the sequence. The entries in the list may be 1, 2
        or 3 character entries if the chain is RNA, DNA or amino acids
        respectively.
        """

        sequence = []
        for row in self.pdbx_poly_seq_scheme:
            if chain != row['pdb_strand_id']:
                continue
            sequence.append(row['mon_id'])
        return sequence

    def experimental_sequence_mapping(self, chain):
        """Create a mapping between the observed sequences and the experimental
        sequences. This allows the determination of what residues are
        unobserved and which are observed as well as where the gaps in the
        structure are. This will prevent duplicate mappings from being created.
        In some cases, like 4X4N, there are duplicate entries for a single unit
        id like position.
        :chain: Name of the chain to use.
        :returns: An iterable that produces the sequence, the sequence unit id,
        the unit id.
        """

        chain_compare = ft.partial(op.eq, chain)
        if isinstance(chain, (list, tuple, set)):
            chain_compare = ft.partial(op.contains, set(chain))

        pdb = self.data.getName()
        mapping = coll.defaultdict(list)
        for residue in self.__residues__(pdb):
            if chain_compare(residue.chain):
                key = (residue.chain, residue.number, residue.insertion_code)
                mapping[key].append(residue.unit_id())
        mapping = dict(mapping)

        entries = self.pdbx_poly_seq_scheme
        filtered = [r for r in entries if chain_compare(r['pdb_strand_id'])]

        # So in some structures, such as 4X4N, there is more than one entry for
        # the same seq id but with a different sequence, ie, position 29 has
        # two entries one as an A and one as a G. Looking at the PDB page I see
        # that it uses the first entry and A. Without being able to pick which
        # one is 'correct', we will just use the first one. Thus our usage of a
        # prev parameter to allow use to skip producing a mapping if the last
        # unit we have seen is the same as the current unit. This also forces
        # us into having an index variable, and not just using enumerate.
        prev = None
        index = 0
        seen = set()
        for row in filtered:
            current_chain = row['pdb_strand_id']
            insertion_code = row['pdb_ins_code']
            if insertion_code == '.':
                insertion_code = "" #Should this be empty string or none

            number = row['pdb_seq_num']
            if number == '?':
                self.logger.warning("Bad seq number pdbx_poly_seq_scheme "
                                    "entry at %s", row)
                continue

            number = int(number)
            key = (current_chain, number, insertion_code)

            # Here is where we skip if we have a duplicate seq_id entry. We do
            # not skip at the level of unit_id because there may be more than
            # one unit mapping to a seq_id for units with alt ids.
            if prev is not None and key == prev:
                continue

            # It is possible that this is error prone, since I think there is a
            # formal constraint on pdbx_poly_seq_scheme to require that units
            # appear grouped by chain. That said, I've never seen a case where
            # this isn't true.
            if prev is not None and prev[0] != current_chain:
                index = 0

            prev = key
            unit_ids = mapping.get(key, [None])
            seq_data = (pdb, current_chain, row['mon_id'], number)
            seq_id = '%s|Sequence|%s|%s|%s' % seq_data
            if insertion_code:
                seq_id += ('|||%s' % insertion_code)

            if seq_id in seen:
                raise ValueError("Can't map one sequence residue %s twice" %
                                 seq_id)

            seen.add(seq_id)
            for unit_id in unit_ids:
                if unit_id in seen:
                    raise ValueError("Can't map one unit id %s twice" %
                                     unit_id)

                if unit_id is not None:
                    seen.add(unit_id)

                yield {
                    'unit_id': unit_id,
                    'seq_id': seq_id,
                    'seq_unit': row['mon_id'],
                    'index': index,
                    'number': number,
                    'chain': current_chain,
                }
            index += 1

    def __breaks__(self):
        pass

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

    def __residues__(self, pdb):
        if sys.version_info[0] < 3:
            key = op.attrgetter(
                'pdb',
                'model',
                'chain',
                'component_id',
                'component_number',
                'insertion_code',
                'symmetry',
        )
        else:
            key = op.attrgetter(
                'pdb',
                'model',
                'chain',
                'component_id',
                'component_number',
                #'insertion_code', #in python 3, the sorted function below cannot accept None values, which sometimes there are in insertion code
                'symmetry',
            )

        mapping = it.groupby(sorted(self.__atoms__(pdb), key=key), key)


        for comp_id, all_atoms in mapping:
            for atoms in self.__group_alt_atoms__(list(all_atoms)):
                first = atoms[0]
                atom_type = self._chem.get(first.component_id, {})
                atom_type = atom_type.get('type', None)
                alt_id = first.alt_id
                if alt_id == '.':
                    alt_id = None

                yield Component(
                    atoms,
                    pdb=first.pdb,
                    model=first.model,
                    type=atom_type,
                    alt_id=alt_id,
                    chain=first.chain,
                    symmetry=first.symmetry,
                    sequence=first.component_id,
                    number=first.component_number,
                    index=first.component_index,
                    insertion_code=first.insertion_code,
                    polymeric=first.polymeric,
                )

    def __atoms__(self, pdb):
        try:
            max_operators = max(len(op) for op in list(self._assemblies.values()))
        except:
            max_operators=1 #if there aren't any operators, there should be one operator applied and it's the identity

        if not self._assemblies:
            if sys.version_info[0] < 3:
                operators = it.repeat(self.__identity_operator__(), len(self.atom_site))
                pdbs = it.repeat(pdb, len(self.atom_site))
                zipped = it.izip(pdbs, self.atom_site, operators)
                return it.imap(lambda a: self.__atom__(*a), zipped)
            else:
                operators = it.repeat(self.__identity_operator__(), len(self.atom_site)) 
                pdbs = it.repeat(pdb, len(self.atom_site))
                zipped = list(zip(pdbs, self.atom_site, operators))
                return list(map(lambda a: self.__atom__(*a), zipped))

        def operator(entry):
            pdb, atom, number = entry
            operators = self.operators(atom['label_asym_id'])
            if not operators:
                self.logger.warning("No operator found for %s", atom)
                return None
            if number < len(operators):
                return pdb, atom, operators[number]
            return None

        atoms = []
        if sys.version_info[0] < 3:
            for index in xrange(max_operators):
                indexes = it.repeat(index, len(self.atom_site))
                pdbs = it.repeat(pdb, len(self.atom_site))
                zipped = it.izip(pdbs, self.atom_site, indexes)
                with_operators = it.imap(operator, zipped)
                filtered = filter(None, with_operators)
                atoms.append(it.imap(lambda a: self.__atom__(*a), filtered))
        else:
            for index in range(max_operators):
                indexes = it.repeat(index, len(self.atom_site))
                pdbs = it.repeat(pdb, len(self.atom_site))
                zipped = (zip(pdbs, self.atom_site, indexes))
                with_operators = list(map(operator, zipped))
                filtered = filter(None, with_operators)
                atoms.append(map(lambda a: self.__atom__(*a), filtered))
        return it.chain.from_iterable(atoms)

    def __atom__(self, pdb, atom, symmetry):
        x, y, z = self.__apply_symmetry__(atom, symmetry)

        index = atom['label_seq_id'] if 'label_seq_id' in atom else '.'
        if index != '.' and index != '':
            index = int(index)
        else:
            index = None

        symmetry_name = self.__symmetry_name__(symmetry)
        ins_code = atom['pdbx_PDB_ins_code'] if 'pdbx_PDB_ins_code' in atom else '?'
        if ins_code == '?':
            ins_code = None

        alt_id = atom['label_alt_id'] if 'label_alt_id' in atom else '.'
        if alt_id == '.':
            alt_id = None

        model = atom['pdbx_PDB_model_num'] if 'pdbx_PDB_model_num' in atom else 1
        component_id = atom['label_comp_id'] if 'label_comp_id' in atom else atom['auth_comp_id']
        atom_id = atom['label_atom_id'] if 'label_atom_id' in atom else atom['auth_atom_id']

        #Some authors leave letters or non numerical digits in their seq_id. That letter should be found in pdbx_PDB_ins_code if needed
        try:
            atom_auth_seq_id = int(atom['auth_seq_id'])
        except:
            atom_auth_seq_id = int(re.sub('\D','',atom['auth_seq_id']))

        return Atom(pdb=pdb,
                    model=model,
                    chain=atom['auth_asym_id'],
                    component_id=component_id,
                    component_number = atom_auth_seq_id, 
                    #component_number = atom['auth_seq_id'], #Used to be casted to be an int. Notify if changing to be a string causes any issues anywhere. 
                    component_index=index,
                    insertion_code=ins_code,
                    alt_id=alt_id,
                    x=x, y=y, z=z,
                    group=atom['group_PDB'],
                    type=atom['type_symbol'],
                    name=atom_id,
                    symmetry=symmetry_name,
                    polymeric=self.is_polymeric_atom(atom))

    def __apply_symmetry__(self, atom, symmetry):
        coords = [float(atom['Cartn_x']),
                  float(atom['Cartn_y']),
                  float(atom['Cartn_z']),
                  1.0]
        result = np.dot(symmetry['transform'], np.array(coords))
        return result[0:3].T
    
    def __symmetry_name__(self, symmetry):
        symmetry_name = symmetry.get('name')
        # if symmetry.get('type') == 'identity operation': #a small handful of cif files have a missing name for a symmetry that is labelled as an ID matrix. See 5A9Z for an example of this.
        #    return '1_555'
        if self.pdb in oldStructures:
            symmetry_name = 'P_%s' % symmetry['id'] # For our database, unit id needs to be the same as it used to be so this is for backward compatibility with unit ids created for and used by the BGSU database. 
                                                    # This is only applied for the structures in the old symmetry list.
        elif not symmetry_name or symmetry_name == '?': 
            symmetry_name = 'ASM_%s' % symmetry['id'] #we've decided this is the best way to annotate these symmetries going forward as they're not named and this is what Cathy Lawson recommended.

        return symmetry_name
    
    def table(self, name):
        return Table(self, self.__block__(name))

    def has_table(self, name):
        block_name = re.sub('^_', '', name)
        block = self.data.getObj(block_name)
        return bool(block)

    def operators(self, asym_id):
        assemblies = self._assemblies[asym_id]

        if not assemblies:
             self.logger.warning("Asym id %s.%s is not part of any assemblies."
                                 " Defaulting to identity operator",
                                 self.pdb, asym_id)
             assemblies = [self.__identity_operator__()]

        #This logic is now being covered in load_assemblies method
        # seen = set()
        # matching = []
        # for assembly in assemblies:
        #     if assembly['id'] not in seen:
        #         seen.add(assembly['id'])
        #         matching.append(assembly)

        return assemblies

    def is_water(self, entity_id):
        return self._entities[entity_id]['type'] == 'water'

    def is_polymeric(self, entity_id):
        return self._entities[entity_id]['type'] == 'polymer'

    def is_polymeric_atom(self, atom):
        return 'label_entity_id' in atom and self.is_polymeric(atom['label_entity_id'])

    def __block__(self, name):
        block_name = re.sub('^_', '', name)
        if sys.version_info[0] < 3:
            block = self.data.getObj(block_name)
        else: 
            block = self.data.get_object(block_name)
        if not block:
            raise MissingBlockException("Unknown block " + name)
        return block

    def __getattr__(self, name):
        try:
            return self.table(name)
        except MissingBlockException:
            raise AttributeError("Unknown block " + name)


class Table(object):

    """Container for a single table in the data block. This provides some
    useful methods for accessing the data.
    """

    def __init__(self, cif, block, rows=None):
        self._cif = cif
        self.block = block
        self.rows = rows
        
        if sys.version_info[0] < 3:
            self.columns = self.block.getItemNameList()
        else:
            self.columns = self.block.item_name_list
        self.columns = [re.sub('_.+\.', '', name) for name in self.columns]

        if self.rows is None:
            if sys.version_info[0] < 3:
                length = self.block.getRowCount()
                self.rows = [self.__row__(index) for index in xrange(length)]
            else:
                length = self.block.row_count
                self.rows = [self.__row__(index) for index in range(length)]

    def column(self, name):
        """Get a column by name"""
        if name not in self.columns:
            raise MissingColumn("Unknown column")

        values = []
        for row in self.rows:
            values.append(row[name])
        return values

    def size(self):
        """Get a tuple of (rowCount, columnCount).
        """
        return (len(self), len(self.columns))

    def __row__(self, number):
        """Get a row by index. Note that this may or may not be in the same
        order as they appear in the cif file, since cif files are not required
        to be ordered. The row will be a dict of the form { attribute: value }.
        Each attribute will have the name of the block stripped.
        """
        if sys.version_info[0] < 3:
            return dict(list(zip(self.columns, self.block.getRow(number))))
        else:
            return dict(list(zip(self.columns, self.block.row_list[number])))

    def __getattr__(self, name):
        """Get the column with the given name.
        """
        try:
            return self.column(name)
        except MissingColumn:
            raise AttributeError("Unknown column: %s" % name)

    def __getitem__(self, index):
        if isinstance(index, str):
            try:
                return self.column(index)
            except MissingColumn:
                raise KeyError("Unknown column: %s" % index)

        if isinstance(index, int):
            return self.rows[index]

        if isinstance(index, slice):
            return Table(self._cif, self.block, rows=self.rows[index])

        raise TypeError("Unknown key type, should be str, int or slice")

    def __len__(self):
        """Get the number of rows.
        """
        return len(self.rows)
