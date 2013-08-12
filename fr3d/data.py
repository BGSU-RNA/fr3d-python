"""This is a package that contains the basic data structures that FR3D works
on, such as Atoms and Components.

"""

from collections import Mapping

from numpy import array

from fr3d.unit_ids import encode


class Entity(Mapping):
    def __init__(self, data):
        self._data = data

    def unit_id(self):
        return encode(self.__rename__())

    def __getitem__(self, item):
        return self._data[item]

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)


class EntityContainer(object):
    def __getter__(self, obj, **kwargs):
        orderby = kwargs.pop('order_by', None)
        compare = kwargs.pop('cmp', None)

        checker = self.__checker__(**kwargs)
        raw = [entry for entry in obj if checker(entry)]

        if orderby:
            if not callable(orderby):
                orderby = lambda entry: entry.get(orderby, None)
            raw.sort(key=orderby)

        if compare:
            raw.sort(cmp=compare)

        return raw

    def __make_check__(self, key, value):

        def check(obj):
            given = obj.get(key, None)

            if callable(value):
                return value(given)

            if isinstance(value, list):
                return given in value

            return value == given

        return check

    def __checker__(self, **kwargs):

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
    def __init__(self, data):
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
        """
        return array([self['x'], self['y'], self['z']])

    def __repr__(self):
        return '<Atom: %s>' % self._data


class Component(Entity, EntityContainer):
    def __init__(self, data, atoms):
        self._atoms = atoms
        super(Component, self).__init__(data)

    def atoms(self, **kwargs):
        return self.__getter__(self._atoms, **kwargs)

    def coordinates(self, **kwargs):
        """Get the coordaintes of all atoms in this component.
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
    def __init__(self, data, residues):
        self._residues = residues
        super(self, Structure).__init__(data)

    def residues(self, **kwargs):
        return self.__getter__(self._residues, **kwargs)
