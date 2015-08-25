"""This module contains the
"""

import itertools as it


def by_atom(first_atoms, second_atoms, cutoff):
    """Create a function to filter pairs by atom atom distances.
    """
    def filter(pair):
        print('by_atom', pair[0], pair[1])
        return pair[0].atoms_within(pair[1], using=first_atoms,
                                    to=second_atoms, cutoff=cutoff)
    return filter


def by_center(first_atoms, second_atoms, cutoff):
    """Create a function to filter pairs by center-center distances.
    """
    def filter(p):
        distance = p[0].distance(p[1], using=first_atoms, to=second_atoms)
        return distance <= cutoff
    return filter


class Pairs(object):
    """This class provides a way to iterate over pairs in a structure. This
    class is an iterator over pairs of residues in a structure. Without
    specifing anything then this will provide an iterator over all pairs. By
    specifing things using the first, second and distance methods this will
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
        cutoff to use, weather or not to use atom-atom or center-center and the
        sets of atoms to use. If the sets of atoms to use are specified then
        tese are the atoms used for all further calculations.

        :first_atoms: A list of atoms to use when comparing distance from the
        first residue. If not given it implies all atoms.
        :second_atoms: A list of atoms to use when comparing distances to the
        second residue. If not given it iplies all atoms.
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

        if self._distance and 'cutoff' not in self._distance:
            raise ValueError("Cannot filter by distance without cutoff")

        # Lazily compute all possible pairs
        pairs = it.product(self.structure.residues(**self._first),
                           self.structure.residues(**self._second))

        # Exclude pairs of 1 component
        pairs = it.ifilter(lambda (a, b): a != b, pairs)

        first_atoms = self._distance.get('first_atoms')
        second_atoms = self._distance.get('second_atoms')

        if self._distance:
            fn = by_center(first_atoms, second_atoms, self._distance['cutoff'])
            pairs = it.ifilter(fn, pairs)

        if self._distance.get('use') == 'atoms':
            fn = by_atom(first_atoms, second_atoms, self._distance['cutoff'])
            pairs = it.ifilter(fn, pairs)

        return pairs
