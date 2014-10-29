import re
import collections as coll

import numpy as np

from pdbx.reader.PdbxParser import PdbxReader as Reader

from fr3d.data import Atom
from fr3d.data import Component
from fr3d.data import Structure


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


class CIF(object):
    """Top level container for all CIF related data. This assumes that each
    mmCIF file contains a single datablock. This doesn't have to be true but
    makes things easier.
    """

    def __init__(self, handle):
        reader = Reader(handle)
        self.data = []
        reader.read(self.data)
        self.data = self.data[0]
        self.pdb = self.data.getName()
        self._assemblies = self.__load_assemblies__()
        self._entities = self.__load_entities__()
        self._chem = self.__load_chem_comp__()

    def __load_assemblies__(self):
        operators = dict((op['id'], op) for op in self.pdbx_struct_oper_list)
        assemblies = coll.defaultdict(list)
        seen = coll.defaultdict(set)
        for assembly in self.pdbx_struct_assembly_gen:
            operator = assembly['oper_expression']
            if operator not in operators:
                raise ComplexOperatorException()
            for asym_ids in assembly['asym_id_list'].split(','):
                for asym_id in asym_ids:
                    if operators[operator]['id'] not in seen[asym_id]:
                        assemblies[asym_id].append(operators[operator])
                        seen[asym_id].add(operators[operator]['id'])
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
        """Get the list of a structures in the CIF file.

        :returns: A list of all structures in the CIF file.
        """

        pdb = self.data.getName()
        residues = self.__residues__(pdb)
        return Structure(residues, pdb=pdb)

    def experimental_sequence(self, chain):
        sequence = []
        for row in self.pdbx_poly_seq_scheme:
            if chain != row['asym_id']:
                continue
            sequence.append(row['mon_id'])
        return sequence

    def experimental_sequence_mapping(self, chain):
        mapping = []
        seen = set()
        for row in self.pdbx_poly_seq_scheme:
            if chain != row['asym_id']:
                continue

            insertion_code = row['pdb_ins_code']
            if insertion_code == '.':
                insertion_code = None

            auth_number = row['auth_seq_num']
            if auth_number == '?':
                unit_id = None
            else:
                unit_id = UIDGenerator({
                    'pdb': self['pdb'],
                    'model': '*',
                    'chain': chain,
                    'residue': row['auth_mon_id'],
                    'number': auth_number,
                    'insertion_code': insertion_code,
                    'symmetry_operator': '*',
                })

            seq_id = '%s|Sequence|%s|%s|%s' % (self['pdb'], self['chain'],
                                               row['mon_id'], row['seq_id'])

            if seq_id in seen:
                raise ValueError("Can't map one sequence residue twice")
            if unit_id and unit_id in seen:
                raise ValueError("Can't map unit %s twice", unit_id)

            seen.add(seq_id)
            seen.add(unit_id)
            mapping.append((row['mon_id'], seq_id, unit_id))

        return mapping

    def __breaks__(self):
        pass

    def __residues__(self, pdb):
        mapping = coll.defaultdict(list)
        for atom in self.__atoms__(pdb):
            mapping[atom.component_unit_id()].append(atom)

        residues = []
        for comp_id, atoms in mapping.items():
            # TODO: Set residue data
            first = atoms[0]
            type = self._chem.get(first.component_id)
            if type:
                type = type['type']
            residues.append(Component(atoms,
                                      pdb=first.pdb,
                                      model=first.model,
                                      type=type,
                                      chain=first.chain,
                                      symmetry=first.symmetry,
                                      sequence=first.component_id,
                                      number=first.component_number,
                                      index=first.component_index,
                                      ins_code=first.ins_code,
                                      polymeric=first.polymeric))
        residues.sort(key=lambda r: r.number)
        return residues

    def __atoms__(self, pdb):
        for atom in self.atom_site:
            for symmetry in self.operators(atom['label_asym_id']):
                if not symmetry:
                    raise MissingSymmetry("Could not find symmetry for %s" %
                                          atom)

                x, y, z = self.__apply_symmetry__(atom, symmetry)
                index = atom['label_seq_id']
                if index != '.':
                    index = int(index)
                ins_code = atom['pdbx_PDB_ins_code']
                if ins_code == '?':
                    ins_code = None

                yield Atom(pdb=pdb,
                           model=int(atom['pdbx_PDB_model_num']),
                           chain=atom['auth_asym_id'],
                           component_id=atom['label_comp_id'],
                           component_number=int(atom['auth_seq_id']),
                           component_index=index,
                           ins_code=ins_code,
                           x=x, y=y, z=z,
                           name=atom['label_atom_id'],
                           symmetry=symmetry['name'],
                           polymeric=self.is_polymeric_atom(atom))

    def __apply_symmetry__(self, atom, symmetry):
        coords = [atom['Cartn_x'], atom['Cartn_y'], atom['Cartn_z']]
        return [float(coord) for coord in coords]

    # def symmetry_operators(self, **kwargs):
    #     atoms = sorted(self.atom_site.rows, key=atom_sorter)
    #     for operator in self.pdbx_struct_oper_list:
    #         fn = lambda a: operator in self.operators(a['label_asym_id'])
    #         yield Symmetry(self, operator, it.ifilter(fn, atoms), **kwargs)

    # def symmetry_operator(self, name, **kwargs):
    #     operator = None
    #     for row in self.pdbx_struct_oper_list:
    #         if row['name'] == name:
    #             operator = row
    #             break

    #     if not operator:
    #         return None

    #     atoms = sorted(self.atom_site.rows, key=atom_sorter)
    #     fn = lambda a: operator in self.operators(a['label_asym_id'])
    #     op = Symmetry(self, operator, it.ifilter(fn, atoms), **kwargs)

    #     if not op:
    #         return None
    #     return op

    def table(self, name):
        return Table(self, self.__block__(name))

    def operators(self, asym_id):
        return self._assemblies[asym_id]

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
