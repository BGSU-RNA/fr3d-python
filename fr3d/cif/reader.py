import re

from collections import defaultdict

import numpy as np

from pdbx.reader.PdbxParser import PdbxReader

from fr3d.data import Atom
from fr3d.data import Component
from fr3d.data import Structure


class MissingOperatorTableError(Exception):
    pass


class InvalidSymmetry(Exception):
    pass


class CIF(object):
    def __init__(self, stream):
        self._data = []
        reader = PdbxReader(stream)
        reader.read(self._data)

    def structures(self):
        return [self.__structure__(entry) for entry in self._data]

    def __structure__(self, entry):
        residues = self.__residues__(entry)
        return Structure({'pdb': entry.getName()}, residues)

    def __residues__(self, entry):
        mapping = defaultdict(list)
        for atom in self.__atoms__(entry):
            mapping[atom.component_unit_id()].append(atom)

        residues = []
        for comp_id, atoms in mapping.items():
            # TODO: Set residue data
            first = atoms[0]
            residues.append(Component({
                'sequence': first.component_id,
                'number': first.component_number,
                'ins_code': first.ins_code,
            }, atoms))
        return residues

    def __atoms__(self, entry):

        for atom in self.__table__(entry, 'atom_site'):
            for symmetry in self.__find_symmetries__(entry, atom):
                if not symmetry:
                    raise InvalidSymmetry

                x, y, z = self.__apply_symmetry__(atom, symmetry)
                # TODO: Set up atom data
                yield Atom({
                    'pdb': entry.getName(),
                    'model': atom['label_entity_id'],
                    'chain': atom['label_asym_id'],
                    'component_id': atom['label_comp_id'],
                    'component_number': atom['label_seq_id'],
                    'ins_code': atom['pdbx_PDB_ins_code'],
                    'x': x,
                    'y': y,
                    'z': z,
                    'name': atom['label_atom_id'],
                    'symmetry': symmetry['name']
                })

    def __apply_symmetry__(self, atom, symmetry):
        # TODO: Apply the symmetries
        return atom['Cartn_x'], atom['Cartn_y'], atom['Cartn_z']
        #coords = np.array([atom['x'], atom['y'], atom['z']])
        #transformed = np.mul(coords, symmetry['matrix'])
        #return transformed[0], transformed[1], transformed[2]

    def __find_symmetries__(self, entry, atom):
        """Compute the symmetry operator for the atom.
        """

        oper_table = entry.getObj('pdbx_struct_oper_list')
        if not oper_table:
            raise MissingOperatorTableError

        assembly_gen = entry.getObj('pdbx_struct_assembly_gen')
        if not assembly_gen:
            return {
                'name': '1_555',
                'translate': 0,
                'rotate': np.array([[1.0, 0.0, 0.0],
                                    [0.0, 1.0, 0.0],
                                    [0.0, 0.0, 1.0]])
            }

        # TODO: Find the symmetries
        return [{
            'name': '1_555',
            'translate': 0,
            'rotate': np.array([[1.0, 0.0, 0.0],
                                [0.0, 1.0, 0.0],
                                [0.0, 0.0, 1.0]])
        }]

    def __table__(self, entry, name):
        block = entry.getObj(name)
        columns = block.getItemNameList()
        columns = [re.sub('_.+\.', '', column) for column in columns]
        for index in xrange(block.getRowCount()):
            yield dict(zip(columns, block.getRow(index)))
