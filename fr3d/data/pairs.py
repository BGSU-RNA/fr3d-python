"""This module contains the
"""

import itertools as it
import sys


class Pairs(object):
    """This class provides a way to iterate over pairs in a structure. This
    class is an iterator over pairs of residues in a structure. Without
    specifying anything then this will provide an iterator over all pairs. By
    specifying things using the first, second and distance methods this will
    limit what is iterated over.
    """

    def __init__(self, structure):
        self.structure = structure
        self._first = {}
        self._second = {}
        self._distance = {}

    def first(self, **kwargs):
        """Define the first set of nucleotides.

        :kwargs: Keyword arguments, that will be used to get the residues.
        """
        self._first = dict(kwargs)

    def second(self, **kwargs):
        """Define the second set of nucleotides.
        """
        self._second = dict(kwargs)

    def distance(self, first_atoms=None, second_atoms=None, cutoff=None,
                 use=None):

        """Define the distance cutoffs. This allows the definition of the
        cutoff to use, whether or not to use atom-atom or center-center and the
        sets of atoms to use. If the sets of atoms to use are specified then
        base are the atoms used for all further calculations.

        :first_atoms: A list of atoms to use when comparing distance from the
        first residue. If not given it implies all atoms.
        :second_atoms: A list of atoms to use when comparing distances to the
        second residue. If not given it applies all atoms.
        :cutoff: The distance cutoff.
        :use: The method to use, either center or atoms. Center means things
        will be selected if their center-center distance is within the cutoff,
        while atom means things will only be selected if a pair of atoms is
        within the cutoff.
        """

        if cutoff:
            cutoff = float(cutoff)
            if cutoff <= 0.0:
                raise ValueError("Must give a cutoff greater than 0")
            self._distance['cutoff'] = cutoff

        if first_atoms:
            self._distance['first_atoms'] = first_atoms

        if second_atoms:
            self._distance['second_atoms'] = second_atoms

        if use:
            if use not in set(['atoms', 'center']):
                raise ValueError("Use must be atoms or center")
            self._distance['use'] = use

        if use == 'center':
            if first_atoms is None:
                self._distance['first_atoms'] = '*'

            if second_atoms is None:
                self._distance['second_atoms'] = '*'

    def __iter__(self):
        """Create the iterator.

        When filtering by distances, this will always use a center-center
        method and then if needed the atom-atom method. This is done because
        atom-atom is expensive and many things can be filtered out ahead of
        time using center-center.

        :returns: An iterator over the specified pairs.
        """

        if self._distance:
            if 'cutoff' not in self._distance:
                raise ValueError("Cannot filter by distance without cutoff")

            cutoff = self._distance['cutoff']
            if self._distance.get('use') == 'atoms':
                # Create trees for the first and second residues and then query
                # to for unique residues within the distance cutoff.
                a1 = {}
                if 'first_atoms' in self._distance:
                    a1 = {'name': self._distance['first_atoms']}

                a2 = {}
                if 'second_atoms' in self._distance:
                    a2 = {'name': self._distance['second_atoms']}
                tree1 = self.structure.atom_distances(residues=self._first,
                                                      atoms=a1)
                tree2 = self.structure.atom_distances(residues=self._second,
                                                      atoms=a2)
            else:
                # Create trees for the first and second atoms in the specified
                # residues and then query to for points within the distance
                # cutoff.
                first_atoms = self._distance.get('first_atoms', None)
                tree1 = self.structure.distances(atoms=first_atoms,
                                                 **self._first)

                second_atoms = self._distance.get('second_atoms', None)
                tree2 = self.structure.distances(atoms=second_atoms,
                                                 **self._second)
            pairs = tree1.neighbors(tree2, cutoff, unique=True)
        else:
            # Lazily compute all possible pairs
            pairs = it.product(self.structure.residues(**self._first),
                               self.structure.residues(**self._second))

        # Exclude pairs of 1 component
        if sys.version_info[0] < 3:
            pairs = it.ifilter(lambda (a, b): a != b, pairs)
        else:
            pairs = it.filter(lambda (a, b): a != b, pairs)

        return pairs
