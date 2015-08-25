"""This is a package that contains the basic data structures that FR3D works
on, such as Atoms and Components.

"""

import collections as col
import itertools as it
import operator as op

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

    def lookup(self, names, allow_missing=True):
        """Lookup a center but allow for missing atoms. This will attempt to lookup
        all atoms but simply ignore those that are missing when computing the
        center. If no atoms are present an empty list is returned. Unlike a
        simple [] this will not raise a KeyError if an atom is missing.

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
            return coords

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
        return self._data.get(key, [])

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
