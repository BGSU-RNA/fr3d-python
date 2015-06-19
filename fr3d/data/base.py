"""This is a package that contains the basic data structures that FR3D works
on, such as Atoms and Components.

"""

import collections as col
import itertools as it

import numpy as np


class EntitySelector(object):
    """This serves as a generic container for entities. We always want to
    provide some methods like getting all entities in a particular order, or
    getting all entities with a particular type. This serves to provide such
    functionality to all classes that contain entities.
    """

    def __init__(self, obj, **kwargs):
        self.obj = obj
        self.options = kwargs

    def __collection_filter__(self, key, collection):
        data = set(collection)
        return lambda obj: getattr(obj, key, None) in data

    def __callable_filter__(self, key, func):
        return lambda obj: func(getattr(obj, key, None))

    def __basic_filter__(self, key, value):
        return lambda obj: getattr(obj, key, None) == value

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
                filtered = it.ifilter(self.__collection_filter__(key, value),
                                      filtered)
            else:
                filtered = it.ifilter(self.__basic_filter__(key, value),
                                      filtered)

        return filtered


class AtomProxy(col.MutableMapping):
    """This class is meant to serve as a way to provide both dictonary like
    access to center data, as well as allow for getting the position of an atom
    as a center.
    """

    def __init__(self, atoms):
        self._atoms = atoms
        self._data = {}

    def define(self, name, atoms):
        """Define a center to be computed later. This will make it possible to
        use the name to access the center, but unlike simply setting it, the
        center will not be computed until accessed.

        :name: The name of the center.
        :atoms: A list of atoms to use to compute the center.
        """
        self._data[name] = set(atoms)

    def lookup(self, names, allow_missing=True):
        """Lookup a center but allow for missing atoms. This will attempt to lookup
        all atoms but simply ignore those that are missing when computing the
        center. If no atoms are present an empty list is returned. Unlike a
        simple [] this will not raise a KeyError if an atom is missing.

        :names: The name(s) to use to compute the center.
        :allow_missing: Flag for allowing missing atoms.
        """
        return self.__handle_key__(names, allow_missing=allow_missing)

    def __coordinates__(self, names, allow_missing=False):
        coordinates = []
        for atom in self._atoms:
            if atom.name in names or names == set('*'):
                coordinates.append(atom.coordinates())

        if len(coordinates) < len(names) and not allow_missing:
            raise KeyError("Missing coordinates for: %s" % ', '.join(names))

        if len(coordinates) == 0:
            return []

        if len(coordinates) == 1:
            return coordinates[0]

        return np.average(coordinates, axis=0)

    def __handle_key__(self, key, **kwargs):
        if isinstance(key, (list, set, tuple)):
            return self.__coordinates__(set(key), **kwargs)
        elif key not in self._data:
            return self.__coordinates__(set([key]), **kwargs)
        elif isinstance(self._data[key], set):
            self._data[key] = self.__coordinates__(self._data[key], **kwargs)

        if not kwargs.get('allow_missing'):
            return self._data[key]
        return self._data.get('allow_missing', [])

    def __getitem__(self, key):
        return self.__handle_key__(key, allow_missing=False)

    def __setitem__(self, key, value):
        self._data[key] = value

    def __delitem__(self, key):
        del self._data[key]

    def __iter__(self):
        for key in self._data.keys():
            yield key

        for atom in self._atoms:
            yield atom.name

    def __len__(self):
        return len(self._data) + len(self._atoms)

    def __repr__(self):
        return str(self._data)
