"""This module contains a single class, the Structure. This class is meant to
represent an entire structure or selected parts of it.
"""

import itertools as it

from fr3d.data.base import EntitySelector
from fr3d.data.base import CoordinateTree
from fr3d.data.pairs import Pairs
from fr3d.unit_ids import encode


class Structure(object):
    """This is a container for a group of residues in a structure. It can serve
    to contain a whole structure with several models and chains or can be
    filtered down to only have selected subsets. It provides useful utility
    methods for grouping and dealing with a collection of components.
    """

    def __init__(self, residues, pdb=None, model=None, chain=None,
                 symmetry=None, breaks=None, operators=None):

        """Create a new Structure.

        :residues: A list of Components that represent the residues in this
        Structure.
        :pdb: The pdb this container is a part of.
        :model: The model number of this Structure.
        :chain: The chain name of this Structure.
        :symmetry: The symmetry operator for this Structure.
        """

        self.pdb = pdb
        self.model = model
        self.chain = chain
        self.symmetry = symmetry
        self._residues = residues
        self._breaks = breaks
        self._operators = operators or {}
        values = self._operators.values()
        self._known_names = set([op['name'] for op in values])
        self._sequence = None

    def residues(self, **kwargs):
        """Get residues from this structure. The keyword arguments work as
        described by EntitySelector.

        :kwargs: Keywords for filtering and ordering
        :returns: The requested residues.
        """
        if 'polymeric' not in kwargs:
            kwargs['polymeric'] = True
        if kwargs.get('polymeric', False) is None:
            kwargs.pop('polymeric')

        return EntitySelector(self._residues, **kwargs)

    def calculate_rotation_matrix(self):
        """ Calculate rotation matrix for bases, modified bases, ...
        """

        for residue in self._residues:
            residue.calculate_rotation_matrix()

    def infer_hydrogens(self):
        """ Infers hydrogen atoms for all residues.
        """

        for residue in self._residues:
            residue.infer_hydrogens()

    def residue(self, unit_id):
        """Get a component by unit id or index. If there is no component at the
        given index, or no component with the given unit id then an IndexError
        is raised.

        :unit_id: The unit id or index of a component.
        :returns: The requested component.
        """

        if isinstance(unit_id, int):
            return self._residues[unit_id]

        filtered = self.residues(unit_id=unit_id)
        found = next(iter(filtered), None)
        if found is None:
            raise IndexError("Unknown residue %s" % unit_id)
        return found

    def select(self, **kwargs):
        """Select a subset of this structure. This can be used to create a new
        structure with only the residues in a specific chain, or chains, with a
        symmetry operator, etc.

        :kwargs: The arguments to use for selecting the residues.
        :returns: A new structure of the selected residues.
        """

        kwargs['pdb'] = self.pdb
        data = {}
        for key in ['pdb', 'model', 'chain', 'symmetry']:
            data[key] = kwargs.get(key)
        data['breaks'] = self._breaks
        return Structure(self.residues(**kwargs), **data)

    def pairs(self, first={}, second={}, distance={}):
        """Create an iterator for all pairs in the structure. The first and
        second arguments are selectors, that is a dictionary that can contain
        the same arguments given to the residues method. The third argument
        distance is a dictionary

        :first: A selector for the first item in the pair.
        :second: A selector for the second item in the pair.
        :distance: A
        :returns: A iterator that will go over all matching pairs in the
        structure.
        """

        pairs = Pairs(self)
        pairs.first(**first)
        pairs.second(**second)
        pairs.distance(**distance)
        return pairs

    def unit_id(self):
        """Get the unit id for this structure
        """

        return encode({
            'pdb': self.pdb,
            'model': self.model,
            'chain': self.chain,
            'symmetry': self.symmetry
        })

    @property
    def sequence(self):
        """Get the sequence for this structure. The sequence is a list of the
        sequences of all residues in this structure.

        :returns: A list of the sequenec of this Structure.
        """

        if self._sequence is None:
            self._sequence = [r.sequence for r in self.residues()]
        return self._sequence

    def distances(self, atoms=None, **kwargs):
        """Create a coordinate tree for the selected residues. Residues are
        selected using the argumetns from kwargs as for the residues method.
        This will use the center for each residue. It will use in order,
        the given atoms, 'base', 'aa_fg', 'modified_nucleotides' or '*'
        (meaning mean of all atoms) if none of the above are defined.

        :param dict atoms: The center to use. If none is given the defaults are
        used.
        :param dict kwargs: The filter to use.
        :returns: A coordinate tree.
        """

        def fn():
            for residue in self.residues(**kwargs):
                if atoms:
                    yield residue, residue.centers[atoms]
                if 'base' in residue.centers:
                    yield residue, residue.centers['base']
                elif 'aa_fg' in residue.centers:
                    yield residue, residue.centers['aa_fg']
                elif 'modified_nucleotides' in residue.centers:
                    yield residue, residue.centers['modified_nucleotides']
                else:
                    yield residue, residue.centers['*']
        return CoordinateTree(fn())

    def atom_distances(self, residues={}, atoms={}):
        """Create a tree for the atom distances. This will filter the residues
        using the residues using the filter given as for .residues() and filter
        the atoms as .atoms() does for components. The coordinate tree will map
        from atom coordinate to residue so it does not provide info on which
        atoms are near, but instead which residue are.

        :param dict atoms: A filter for the atoms to select
        :param dict residues: The filter for the residues to select.
        :returns: A coordinate tree.
        """

        def fn():
            for residue in self.residues(**residues):
                for atom in residue.atoms(**atoms):
                    yield residue, atom.coordinates()
        return CoordinateTree(fn())

    def __len__(self):
        """Compute the length of this Structure. That is the number of residues
        in this structure.

        :returns: The number of atoms.
        """
        return sum(it.imap(lambda _: 1, self._residues))

    def __bool__(self):
        """Check if this structure is true. A structure is true if the list of
        residues is not empty.

        :returns: True if this structure has any residues.
        """
        return bool(self._residues)

    def __repr__(self):
        """Get a simple string representation of this Structure.

        :returns: A string representation of this Structure.
        """
        return '<Structure: %s>' % self.pdb
