"""This is a package that contains the basic data structures that FR3D works
on, such as Atoms and Components.

"""

import collections as col
import itertools as it
import operator as op

import numpy as np

from scipy import spatial as sp


class EntitySelector(object):
    """This serves as a generic container for entities. We always want to
    provide some methods like getting all entities in a particular order, or
    getting all entities with a particular type. This serves to provide such
    functionality to all classes that contain entities.
    """

    def __init__(self, obj, **kwargs):
        self.obj = obj
        self.options = kwargs

    def __callable_filter__(self, key, func):
        return lambda obj: func(getattr(obj, key, None))

    def __basic_filter__(self, key, value, compare):
        def fn(obj):
            attr = getattr(obj, key, None)
            if callable(attr):
                return compare(attr(), value)
            return compare(attr, value)
        return fn

    def __iter__(self):
        """This method is a way to sort and filter an object easily. The
        keyword arguments, order_by and cmp are used for sorting while
        everything else is used by __checker__ for filtering. The order_by
        keyword may be either a string or a function. If it is a string
        then we create a function which gets that key from the entities
        in obj and uses that as the key function in sorting. All other keywords
        are used as described in __checker__. If no keywords are given then
        the object is simply returned.

        :obj: An iterable to filter and sort
        :kwargs: Keyword arguments for filtering and sorting.
        """

        filtered = iter(self.obj)
        for key, value in self.options.items():
            if key == '_':
                filtered = it.ifilter(value, filtered)
            elif callable(value):
                filtered = it.ifilter(self.__callable_filter__(key, value),
                                      filtered)
            elif isinstance(value, (list, set, tuple)):
                func = self.__basic_filter__(key, set(value),
                                             lambda a, b: op.contains(b, a))
                filtered = it.ifilter(func, filtered)
            else:
                func = self.__basic_filter__(key, value, op.eq)
                filtered = it.ifilter(func, filtered)

        return filtered


class AtomProxy(col.MutableMapping):
    """This class is meant to serve as a way to provide both dictonary like
    access to center data, as well as allow for getting the position of an atom
    as a center.
    """

    def __init__(self, atoms):
        self._atoms = atoms
        self._data = {}
        self._definitions = {}

    def define(self, name, atoms):
        """Define a center to be computed later. This will make it possible to
        use the name to access the center, but unlike simply setting it, the
        center will not be computed until accessed.

        :name: The name of the center.
        :atoms: A list of atoms to use to compute the center.
        """

        self._definitions[name] = atoms

        if isinstance(atoms, basestring):
            self._data[name] = set([atoms])
        else:
            self._data[name] = set(atoms)

    def definition(self, name):
        """Get the definition for the given name. If the name is not defined
        then None is returned.

        :name: The name of the center.
        :returns: The set of atoms, if defined.
        """
        return self._definitions.get(name)

    def definitions(self):
        """Get all defined centers. This will not return the list of atom
        names and is given in no particular order.

        :returns: A list of key names.
        """
        return self._definitions.keys()

    def lookup(self, names, allow_missing=True):
        """Lookup a center but allow for missing atoms. This will attempt to lookup
        all atoms but simply ignore those that are missing when computing the
        center. If no atoms are present an empty numpy array is returned.

        :names: The name(s) to use to compute the center.
        :allow_missing: Boolean for allowing missing atoms, defaults to True.
        """
        return self.__handle_key__(names, allow_missing=allow_missing)

    def __coordinates__(self, names, allow_missing=True):
        coords = []
        if names == set('*'):
            coords = [atom.coordinates() for atom in self._atoms]
        else:
            coords = [a.coordinates() for a in self._atoms if a.name in names]

        if len(coords) < len(names) and not allow_missing:
            raise KeyError("Missing coordinates for: %s" % ', '.join(names))

        if not coords:
            return np.array(coords)

        if len(coords) == 1:
            return coords[0]

        return np.average(coords, axis=0)

    def __handle_key__(self, key, **kwargs):
        if isinstance(key, (list, set, tuple)):
            return self.__coordinates__(set(key), **kwargs)
        elif key not in self._data:
            return self.__coordinates__(set([key]), **kwargs)
        elif isinstance(self._data[key], set):
            self._data[key] = self.__coordinates__(self._data[key], **kwargs)

        if not kwargs.get('allow_missing'):
            return self._data[key]
        return self._data.get(key, np.array([]))

    def __getitem__(self, key):
        return self.__handle_key__(key, allow_missing=True)

    def __setitem__(self, key, value):
        self._data[key] = value

    def __delitem__(self, key):
        del self._data[key]

    def __iter__(self):
        for key in self._data.keys():
            yield key

        for atom in self._atoms:
            yield atom.name

    def __contains__(self, key):
        if key in self._data or key == '*':
            return True

        for atom in self._atoms:
            if atom.name == key:
                return True
        return False

    def __len__(self):
        return len(self._data) + len(self._atoms)

    def __repr__(self):
        return str(self._data)


class CoordinateTree(object):
    """This is a simple wrapper around scipy's KDTree to return components
    instead of just indexes in a list.
    """

    def __init__(self, generator):
        """Create a new CoordinateTree. The given generator should yield 2
        values each time, the residue and a coordinate to use for it. This may
        yield the same residue many times but should not duplicate coordinate
        values.

        :param iterable generator: The generator to use.
        """

        self._residues = []
        coordinates = []
        self.tree = None
        for residue, coordinate in generator:
            if len(coordinate) > 0:
                coordinates.append(coordinate)
                self._residues.append(residue)
        if coordinates:
            self.tree = sp.cKDTree(coordinates)

    def count_neighbors(self, other, r, *p):
        """Return the counts of neighbors in the other tree. Arguments are as
        for cKDTree.count_neighbors, except other is a CoordinateTree. This
        does not uniquify the neighbors before counting.

        :returns: The counts.
        """

        if not self.tree:
            return 0

        return self.tree.count_neighbors(other.tree, r, *p)

    def pairs(self, distance, unique=False, **kwargs):
        """Create a generator over all pairs in this tree which are within the
        given distance cutoff.

        :param float distance: The cutoff.
        :param bool unique: Only get the unique pairs of residues found.
        Uniquess is determined by unit ids.
        :kwargs: Keyword arguments to cKDTree.query_pairs.
        :returns: A generator for the pairs.
        """

        if not self.tree:
            return []

        def fn():
            results = self.tree.query_pairs(distance, **kwargs)
            if results:
                for first, second in results:
                    yield self._residues[first], self._residues[second]

        if unique:
            return self.__as_unique__(fn())
        return fn()

    def neighbors(self, other, distance, unique=False, **kwargs):
        """Create a generator over all points which are within some distance
        cutoff between this tree and another one.

        :param CoordinateTree other: The other tree.
        :param float distance: The cutoff.
        :param bool unique: Only get the unique pairs of found. Uniqueness is
        determined by the unit ids.
        :kwargs: Keyword arguments to cKDTree.query_ball_tree.
        :returns: A generator for the neighbors.
        """

        if not self.tree:
            return []

        def fn():
            results = self.tree.query_ball_tree(other.tree, distance, **kwargs)
            if results:
                for first, r in enumerate(results):
                    for second in r:
                        yield self._residues[first], other._residues[second]

        if unique:
            return self.__as_unique__(fn())
        return fn()

    def __as_unique__(self, generator):
        """Make sure the given generator returns only unique pairs, according
        to the unit_id method.

        :param generator: The generator to make unique.
        :returns: A new generator.
        """

        seen = set()
        for first, second in generator:
            pair = (first.unit_id(), second.unit_id())
            if pair not in seen:
                yield first, second
                seen.add(pair)
