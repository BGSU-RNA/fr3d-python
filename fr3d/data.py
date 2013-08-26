"""This is a package that contains the basic data structures that FR3D works
on, such as Atoms and Components.

"""

from collections import Mapping

from numpy import array

from fr3d.unit_ids import encode


class Entity(Mapping):
    """This class is the base class for things like atoms and other units. It
    is intended to provide a simple dict like access to the data in it as well
    as have a method for generating its unit id. Currently, the data stored in
    this object is not mutable, but this may need to be changed.
    """

    def __init__(self, data):
        self._data = data

    def unit_id(self):
        """Compute the unit id for this Entity.

        :returns: A string of the unit id.
        """

        return encode(self.__rename__())

    def __getitem__(self, item):
        return self._data[item]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)


class EntityContainer(object):
    """This serves as a generic container for entities. We always want to
    provide some methods like getting all entities in a particular order, or
    getting all entities with a particular type. This serves to provide such
    functionality to all classes that contain entities.
    """

    def __getter__(self, obj, **kwargs):
        """This method is a way to sort and filter an object easily. The
        keyword arguments, order_by and cmp are used for sorting while
        everything else is used by __checker__ for filtering. The order_by
        keyword may be either a string or a function. If it is a string
        then we create a function which gets that key from the entities
        in obj and uses that as the key function in sorting. The cmp keyword
        is used as the cmp function in sorting. All other keywords are used as
        described in __checker__. If no keywords are given then the object is
        simply returned.

        :obj: An iterable to filter and sort
        :kwargs: Keyword arguments for filtering and sorting.
        """

        if not kwargs:
            return obj

        orderby = kwargs.pop('order_by', None)
        compare = kwargs.pop('cmp', None)

        checker = self.__checker__(**kwargs)
        raw = [entry for entry in obj if checker(entry)]

        if orderby:
            key = orderby
            if not callable(orderby):
                key = lambda entry: entry.get(orderby, None)
            raw.sort(key=key)

        if compare:
            raw.sort(cmp=compare)

        return raw

    def __make_check__(self, key, value):
        """Generate a function for filtering. If the given value is a callable
        then it is used to filter. The function will be given value of key for
        each object it filters. If value is a list then we check to see if
        the value of key in each object is in the list, otherwise we check to
        see if the value of the key in each object equals the given value.

        :key: The key to filter by.
        :value: The value to use as the filter.
        """

        def check(obj):
            given = obj.get(key, None)

            if callable(value):
                return value(given)

            if isinstance(value, list):
                return given in value

            return value == given

        return check

    def __checker__(self, **kwargs):
        """This generates a function which is used to filter an iterable.
        """

        checkers = []
        for key, value in kwargs.items():
            checkers.append(self.__make_check__(key, value))

        def checker(entity):
            for checker in checkers:
                if not checker(entity):
                    return False
            return True

        return checker


class Atom(Entity):
    """This class represents atoms in a structure. It provides a simple dict
    like access for data as well as a way to get it's coordiantes, unit id
    and the unit id of the component it belongs to.
    """

    def __init__(self, data):
        """Create a new Atom.

        :data: A dictonary of data to provide access to.
        """
        super(Atom, self).__init__(data)

    def component_unit_id(self):
        """Generate the unit id of the component this atom belongs to.

        :returns: A string of the unit id for this atom's component.
        """

        comp_data = self.__rename__()
        if 'atom_name' in comp_data:
            del comp_data['atom_name']
        if 'alt_id' in comp_data:
            del comp_data['alt_id']

        return encode(comp_data)

    def __rename__(self):
        data = dict(self)
        data['atom_name'] = data.pop('name')
        return data

    def coordinates(self):
        """Return a numpy array of the x, y, z coordinates for this atom.

        :returns: A numpy array of the x, y, z coordinates.
        """
        return array([self['x'], self['y'], self['z']])

    def __repr__(self):
        return '<Atom: %s>' % self._data


class Component(Entity, EntityContainer):
    """This represents things like nucleic acids, amino acids, small molecules
    and ligands.
    """

    def __init__(self, data, atoms):
        """Create a new Component.

        :data: The data to provide access to.
        :atoms: The atoms this component is composed of.
        """

        self._atoms = atoms
        super(Component, self).__init__(data)

    def atoms(self, **kwargs):
        """Get, filter and sort the atoms in this component. Access is as
        described by EntityContainer.__getter__.

        :kwargs: The keyword arguments to filter and sort by.
        :returns: A list of the requested atoms.
        """
        return self.__getter__(self._atoms, **kwargs)

    def coordinates(self, **kwargs):
        """Get the coordaintes of all atoms in this component. This will
        filter to the requested atoms, sort and then provide a numpy array
        of all coordinates for the given atoms.

        :kwargs: Arguments to filter and sort by.
        :returns: A numpy array of the coordinates.
        """
        return array([atom.coordinates() for atom in self.atoms(**kwargs)])

    def __rename__(self):
        data = dict(self)
        data['component_id'] = data.pop('sequence')
        data['component_number'] = data.pop('number')
        return data

    def completeness(self):
        # TODO: Compute the completeness of this component
        return False

    def __repr__(self):
        return '<Component %s Atoms: %s>' % (self._data, self._atoms)


class Structure(Entity, EntityContainer):
    """This represents a structure which is composed of components.
    """

    def __init__(self, data, residues):
        self._residues = residues
        super(self, Structure).__init__(data)

    def residues(self, **kwargs):
        """Get residues from this structure. The keyword arguments work as
        described by EntityContainer.__getter__.

        :kwargs: Keywords for filtering and ordering
        :returns: The requested residues.
        """
        return self.__getter__(self._residues, **kwargs)
