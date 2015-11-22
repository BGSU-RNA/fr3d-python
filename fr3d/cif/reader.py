import re
import itertools as it
import collections as coll
import warnings
import logging

import numpy as np

from pdbx.reader.PdbxParser import PdbxReader as Reader

from fr3d.data import Atom
from fr3d.data import Component
from fr3d.data import Structure
from fr3d.unit_ids import encode


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

        self.pdb = self.data.getName()
        self._operators = self.__load_operators__()
        self._assemblies = self.__load_assemblies__()
        self._entities = self.__load_entities__()
        self._chem = self.__load_chem_comp__()
        self.logger = logging.getLogger('fr3d.cif.reader.Cif')

    def __load_operators__(self):
        operators = {}
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
        assemblies = coll.defaultdict(list)
        for assembly in self.pdbx_struct_assembly_gen:
            oper_expression = assembly['oper_expression']

            # TODO: Implement computation of complex operators
            if COMPLEX_SYMBOLS & set(oper_expression):
                warnings.warn('Cannot compute symmetries from complex '
                              'expressions. Will use a simple identity '
                              'transformation if no others possible')
                operators = []
            else:
                operators = oper_expression.split(',')

            for asym_id in assembly['asym_id_list'].split(','):
                for operator in operators:
                    op = self._operators[operator]
                    assemblies[asym_id].append(op)

        for asym_id, ops in assemblies.items():
            if not ops:
                self.logger.info("Adding default identity operator for %s",
                                 asym_id)
                assemblies[asym_id].append(self._operators['I'])

        return assemblies

    def __load_entities__(self):
        entities = {}
        for entity in self.entity:
            entities[entity['id']] = entity
        return entities

    def __load_chem_comp__(self):
        chem = {}
        for obj in self.chem_comp:
            chem[obj['id']] = obj
        return chem

    def structure(self):
        """Get the structure from the Cif file.

        :returns: The first structure in the cif file.
        """

        pdb = self.data.getName()
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

        pdb = self.data.getName()
        mapping = coll.defaultdict(list)
        for residue in self.__residues__(pdb):
            if residue.chain == chain:
                key = (residue.number, residue.insertion_code)
                mapping[key].append(residue.unit_id())
        mapping = dict(mapping)

        entries = self.pdbx_poly_seq_scheme
        filtered = it.ifilter(lambda r: r['pdb_strand_id'] == chain, entries)
        # model = self.atom_site[0]['pdbx_PDB_model_num']

        seen = set()
        prev_number = None
        for index, row in enumerate(filtered):
            insertion_code = row['pdb_ins_code']
            if insertion_code == '.':
                insertion_code = None

            number = row['pdb_seq_num']
            if number == '?' or row['auth_seq_num'] == '?':
                unit_ids = [None]
            else:
                if prev_number == number:
                    self.logger.warning("Duplicate pdbx_poly_seq_scheme "
                                        "entry at %s", number)
                    continue

                prev_number = number
                key = (int(number), insertion_code)

                if key not in mapping:
                    raise ValueError("Could not find unit for %s", key)

                unit_ids = mapping[key]

            seq_data = (pdb, chain, row['mon_id'], row['seq_id'])
            seq_id = '%s|Sequence|%s|%s|%s' % seq_data

            if seq_id in seen:
                raise ValueError("Can't map one sequence residue twice")

            seen.add(seq_id)
            for unit_id in unit_ids:
                seen.add(unit_id)
                yield (row['mon_id'], seq_id, unit_id)

    def __breaks__(self):
        pass

    def __residues__(self, pdb):
        mapping = it.groupby(sorted(self.__atoms__(pdb),
                                    key=lambda a: a.component_unit_id()),
                             lambda a: a.component_unit_id())

        for comp_id, atoms in mapping:
            atoms = list(atoms)
            first = atoms[0]
            type = self._chem.get(first.component_id, {})
            type = type.get('type', None)
            alt_id = first.alt_id
            if alt_id == '.':
                alt_id = None
            yield Component(atoms,
                            pdb=first.pdb,
                            model=first.model,
                            type=type,
                            alt_id=alt_id,
                            chain=first.chain,
                            symmetry=first.symmetry,
                            sequence=first.component_id,
                            number=first.component_number,
                            index=first.component_index,
                            insertion_code=first.insertion_code,
                            polymeric=first.polymeric)

    def __atoms__(self, pdb):
        max_operators = max(len(op) for op in self._assemblies.values())

        if not max_operators:
            raise ValueError("Could not find any operators")

        def operator(entry):
            pdb, atom, number = entry
            operators = self.operators(atom['label_asym_id'])
            if number < len(operators):
                return pdb, atom, operators[number]
            return None

        atoms = []
        for index in xrange(max_operators):
            indexes = it.repeat(index, len(self.atom_site))
            pdbs = it.repeat(pdb, len(self.atom_site))
            zipped = it.izip(pdbs, self.atom_site, indexes)
            with_operators = it.imap(operator, zipped)
            filtered = it.ifilter(None, with_operators)
            atoms.append(it.imap(lambda a: self.__atom__(*a), filtered))

        return it.chain.from_iterable(atoms)

    def __atom__(self, pdb, atom, symmetry):
        x, y, z = self.__apply_symmetry__(atom, symmetry)

        index = atom['label_seq_id']
        if index != '.':
            index = int(index)

        symmetry_name = self.__symmetry_name__(symmetry)

        ins_code = atom['pdbx_PDB_ins_code']
        if ins_code == '?':
            ins_code = None

        alt_id = atom['label_alt_id']
        if alt_id == '.':
            alt_id = None

        return Atom(pdb=pdb,
                    model=int(atom['pdbx_PDB_model_num']),
                    chain=atom['auth_asym_id'],
                    component_id=atom['label_comp_id'],
                    component_number=int(atom['auth_seq_id']),
                    component_index=index,
                    insertion_code=ins_code,
                    alt_id=alt_id,
                    x=x, y=y, z=z,
                    group=atom['group_PDB'],
                    type=atom['type_symbol'],
                    name=atom['label_atom_id'],
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
        if not symmetry_name or symmetry_name == '?':
            symmetry_name = 'P_%s' % symmetry['id']
        return symmetry_name

    def table(self, name):
        return Table(self, self.__block__(name))

    def operators(self, asym_id):
        matching = []
        seen = set()
        assemblies = self._assemblies[asym_id]
        for assembly in assemblies:
            if assembly['id'] not in seen:
                seen.add(assembly['id'])
                matching.append(assembly)
        return matching

    def is_water(self, entity_id):
        return self._entities[entity_id]['type'] == 'water'

    def is_polymeric(self, entity_id):
        return self._entities[entity_id]['type'] == 'polymer'

    def is_polymeric_atom(self, atom):
        return self.is_polymeric(atom['label_entity_id'])

    def __block__(self, name):
        block_name = re.sub('^_', '', name)
        block = self.data.getObj(block_name)
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

        self.columns = self.block.getItemNameList()
        self.columns = [re.sub('_.+\.', '', name) for name in self.columns]

        if self.rows is None:
            length = self.block.getRowCount()
            self.rows = [self.__row__(index) for index in xrange(length)]

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
        return dict(zip(self.columns, self.block.getRow(number)))

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
