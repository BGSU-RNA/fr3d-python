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
            mapping[atom.component_id()].append(atom)

        residues = []
        for comp_id, atoms in mapping.items():
            # TODO: Set residue data
            residues.append(Component({}, atoms))
        return residues

    def __atoms__(self, entry):
        for atom in entry.getObj('_atom_site'):
            for symmetry in self.__find_symmetries__(entry, atom):
                if not symmetry:
                    raise InvalidSymmetry

                x, y, z = self.__apply_symmetry__(atom, symmetry)
                # TODO: Set up atom data
                yield Atom({
                    'x': x,
                    'y': y,
                    'z': z,
                    'symmetry': symmetry['name']
                })

    def __apply_symmetry__(self, atom, symmetry):
        # TODO: Apply the symmetries
        pass
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

        return None
