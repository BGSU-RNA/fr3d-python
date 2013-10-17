from collections import defaultdict

import numpy as np

import CifFile as cf

from fr3d.data import Atom
from fr3d.data import Component
from fr3d.data import Structure


class MissingOperatorTableError(Exception):
    pass


class InvalidSymmetry(Exception):
    pass


class CIF(object):
    """
    """

    def __init__(self, filename):
        """Create a new CIF reader.

        :filename: The filename to parse.
        """

        self._data = cf.ReadCif(filename)

    def structures(self):
        """Get the list of a structures in the CIF file.

        :returns: A list of all structures in the CIF file.
        """

        return [self.__structure__(entry) for entry in self._data.keys()]

    def __structure__(self, name):
        entry = self._data[name]
        residues = self.__residues__(entry, name)
        return Structure(residues, pdb=name)

    def __residues__(self, entry, pdb):
        mapping = defaultdict(list)
        for atom in self.__atoms__(entry, pdb):
            mapping[atom.component_unit_id()].append(atom)

        residues = []
        for comp_id, atoms in mapping.items():
            # TODO: Set residue data
            first = atoms[0]
            residues.append(Component(atoms,
                                      sequence=first.component_id,
                                      number=first.component_number,
                                      ins_code=first.ins_code))
        return residues

    def __atoms__(self, entry, pdb):

        for atom in self.__table__(entry, 'atom_site'):
            for symmetry in self.__find_symmetries__(entry, atom):
                if not symmetry:
                    raise InvalidSymmetry

                x, y, z = self.__apply_symmetry__(atom, symmetry)
                yield Atom(pdb=pdb,
                           model=atom['label_entity_id'],
                           chain=atom['label_asym_id'],
                           component_id=atom['label_comp_id'],
                           component_number=atom['label_seq_id'],
                           ins_code=atom['pdbx_PDB_ins_code'],
                           x=x, y=y, z=z,
                           name=atom['label_atom_id'],
                           symmetry=symmetry['name'])

    def __apply_symmetry__(self, atom, symmetry):
        coords = [atom['Cartn_x'], atom['Cartn_y'], atom['Cartn_z']]
        return [float(coord) for coord in coords]

    def __find_symmetries__(self, entry, atom):
        """Compute the symmetry operator for the atom.
        """
        # TODO: Find the symmetries

        return [{
            'name': '1_555',
            'translate': 0,
            'rotate': np.array([[1.0, 0.0, 0.0],
                                [0.0, 1.0, 0.0],
                                [0.0, 0.0, 1.0]])
        }]

    def __table__(self, block, name):
        if name[0] != '_':
            name = '_' + name

        for loop in block.loops:
            key = loop.keys()[0]
            if key.startswith(name):
                break

        if loop:
            name = name + '.'
            keys = [key.replace(name, '') for key in loop.GetItemOrder()]
            rows = len(loop[loop.keys()[0]])
            for count in xrange(rows):
                yield dict(zip(keys, loop[count]))
