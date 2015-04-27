"""This is a package that contains the basic data structures that FR3D works
on, such as Atoms and Components.

"""

import collections as col
import itertools as it

import numpy as np

from fr3d.unit_ids import encode
from fr3d.definitions import RNAbaseheavyatoms
from fr3d.definitions import RNAbasecoordinates
from fr3d.definitions import RNAbasehydrogens
from fr3d.definitions import aa_sidechain
from fr3d.definitions import aa_backbone
from fr3d.geometry.superpositions import besttransformation


class Entity(object):
    """This class is the base class for things like atoms and other units. It
    is intended to provide a simple dict like access to the data in it as well
    as have a method for generating its unit id. Currently, the data stored in
    this object is not mutable, but this may need to be changed.
    """

    def __init__(self, **kwargs):
        self._data = kwargs
        for key, value in kwargs.items():
            if hasattr(self, key):
                raise ValueError("Can't override properites")
            setattr(self, key, value)

    def unit_id(self):
        """Compute the unit id for this Entity.

        :returns: A string of the unit id.
        """

        return encode(self.__rename__())


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
        in obj and uses that as the key function in sorting. All other keywords
        are used as described in __checker__. If no keywords are given then
        the object is simply returned.

        :obj: An iterable to filter and sort
        :kwargs: Keyword arguments for filtering and sorting.
        """

        if not kwargs:
            return obj

        checker = self.__checker__(**kwargs)
        return list(it.ifilter(checker, obj))

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
            given = getattr(obj, key, None)

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


class AtomProxy(col.MutableMapping):
    """This class is meant to serve as a way to provide both dictonary like
    access to center data, as well as allow for getting the position of an atom
    as a center.
    """

    def __init__(self, atoms):
        self._atoms = atoms
        self._data = {}

    def __coordinates__(self, names):
        coordinates = []
        for atom in self._atoms:
            if atom.name in names:
                coordinates.append(atom.coordinates())

        if len(coordinates) < len(names):
            raise KeyError("Unknown key(s): %s" % names)

        if len(coordinates) == 1:
            return coordinates[0]

        return np.average(coordinates, axis=0)

    def __getitem__(self, key):
        if key not in self._data:
            if isinstance(key, tuple):
                return self.__coordinates__(set(key))
            return self.__coordinates__(set([key]))
        return self._data[key]

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


class Atom(Entity):
    """This class represents atoms in a structure. It provides a simple dict
    like access for data as well as a way to get its coordinates, unit id
    and the unit id of the component it belongs to.
    """

    def __init__(self, **kwargs):
        """Create a new Atom.

        :data: A dictionary of data to provide access to.
        """
        super(Atom, self).__init__(**kwargs)

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
        data = dict(self._data)
        data['atom_name'] = data.pop('name')
        return data

    def transform(self, transform):
        """Create a new atom based of this one, but with transformed
        coordinates.
        """

        coords = [self.x, self.y, self.z, 1.0]
        result = np.dot(transform, np.array(coords))
        x, y, z = result[0:3].T
        data = dict(self._data)
        data['x'] = x
        data['y'] = y
        data['z'] = z

        return Atom(**data)

    def coordinates(self):
        """Return a numpy array of the x, y, z coordinates for this atom.

        :returns: A numpy array of the x, y, z coordinates.
        """
        return np.array([self.x, self.y, self.z])

    def __sub__(self, atom):
        """Compute the distance between this atom and another atom.

        :atom: Another atom.
        :returns: The distance.
        """
        return np.linalg.norm(self.coordinates() - atom.coordinates())

    def __repr__(self):
        return '<Atom: %s>' % self._data


class Component(Entity, EntityContainer):
    """This represents things like nucleic acids, amino acids, small molecules
    and ligands.
    """

    def __init__(self, atoms, **kwargs):
        """Create a new Component.

        :data: The data to provide access to.
        :atoms: The atoms this component is composed of.
        """

        self._atoms = atoms
        super(Component, self).__init__(**kwargs)

        self.centers = AtomProxy(self._atoms)

        if self.sequence in ['A', 'C', 'G', 'U']:
            atoms = RNAbaseheavyatoms[self.sequence]
            self.centers['base'] = self.__compute_center__(atoms)

        if self.sequence in ['ARG', 'LYS', 'HIS', 'GLN','ASN','GLU','ASP','TRP','TYR','PHE','PRO','MET','ILE','LEU','VAL','ALA','SER','CYS','THR']:
            atoms = aa_sidechain[self.sequence]
            self.centers['aa_sidechain'] = self.__compute_center__(atoms)

        if self.sequence in ['ARG', 'LYS', 'HIS', 'GLN','ASN','GLU','ASP','TRP','TYR','PHE','PRO','MET','ILE','LEU','VAL','ALA','GLY','SER','CYS','THR']:
            atoms = aa_backbone[self.sequence]
            self.centers['aa_backbone'] = self.__compute_center__(atoms)

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
        return np.array([atom.coordinates() for atom in self.atoms(**kwargs)])

    def select(self, **kwargs):
        """Select a group of atoms to create a new component out of.

        :kwargs: As for atoms.
        :returns: A new Component
        """
        return Component(self.atoms(**kwargs), **self._data)

    def is_complete(self, names, key='name'):
        """This checks if we can find all atoms in this entity with the given
        names. This assumes that the names for each atom are unique. If you
        wish to use something other than name use the key argument. However, it
        must provide a unique value for each atom, if several atoms with the
        same value are found will cause the function to behave oddly.

        :names: The list of names to check for.
        :key: The key to use for getting atoms. Defaults to name.
        :returns: True if all atoms with the given name are present.
        """
        kwargs = {key: names}
        found = self.atoms(**kwargs)
        return len(found) == len(names)

    def infer_hydrogens(self):
        """Infer the coordinates of the hydrogen atoms of this component.
        Currently, it only works for RNA with .sequence
        """
        if self.sequence not in ['A', 'C', 'G', 'U']:
            return None
        R = []
        S = []
        baseheavy = RNAbaseheavyatoms[self.sequence]

        for atom in self.atoms(name=baseheavy):
            coordinates = atom.coordinates()
            R.append(coordinates)
            S.append(RNAbasecoordinates[self.sequence][atom.name])

        R = np.array(R)
        R = R.astype(np.float)
        S = np.array(S)
        try:
            rotation_matrix, fitted, base_center, rmsd, sse = \
                besttransformation(R, S)
        except:
            return None

        self.rotation_matrix = rotation_matrix

        hydrogens = RNAbasehydrogens[self.sequence]
        coordinates = RNAbasecoordinates[self.sequence]

        for hydrogenatom in hydrogens:
            hydrogencoordinates = coordinates[hydrogenatom]
            newcoordinates = base_center + \
                np.dot(hydrogencoordinates, np.transpose(rotation_matrix))
            self._atoms.append(Atom(name=hydrogenatom,
                                    x=newcoordinates[0, 0],
                                    y=newcoordinates[0, 1],
                                    z=newcoordinates[0, 2]))

    def transform(self, transform):
        """Create a new component from this one by applying a transformation
        matrix.

        :transform: The transformation matrix to apply.
        :returns: A new Component with the same properties by rotated atoms.
        """
        transformed = []
        for atom in self.atoms():
            transformed.append(atom.transform(transform))
        data = dict(self._data)
        return Component(transformed, **data)

    def __rename__(self):
        data = dict(self._data)
        data['component_id'] = data.pop('sequence')
        data['component_number'] = data.pop('number')
        return data

    def __compute_center__(self, atoms):
        """Compute the center position for the given set of atoms. This is done
        through taking the mean position of each atom. If a requested atom
        does not exist it is ignored.

        :atoms: Atoms to use to find the center.
        :returns: The x, y, z coordinates of centers.
        """
        coordinates = [atom.coordinates() for atom in self.atoms(name=atoms)]
        if not coordinates:
            return None
        return np.mean(coordinates, axis=0)

    def __len__(self):
        """Compute the length of this Component. This is the number of atoms in
        this residue.

        :returns: The number of atoms.
        """
        return len(self._atoms)

    def __repr__(self):
        return '<Component %s>' % self.unit_id()


class Model(Entity, EntityContainer):
    def __init__(self, chains, **kwargs):
        self.chains = chains
        super(Model, self).__init__(**kwargs)

    def chain(self, name):
        for chain in self.chains:
            if chain.chain == name:
                return chain
        return None

    def polymers(self):
        for chain in self.chains:
            for polymer in chain.polymers():
                yield polymer

    def residues(self, **kwargs):
        for chain in self.chains:
            for residue in chain.residues(**kwargs):
                yield residue

    def atoms(self, **kwargs):
        for residue in self.residues(**kwargs):
            for atom in residue.atoms():
                yield atom

    def __repr__(self):
        return '<Model: %s|%s>' % (self.pdb, self.model)


class Chain(Entity, EntityContainer):
    def __init__(self, residues, breaks=None, **kwargs):
        self._residues = residues
        self._breaks = breaks
        self._sequence = None
        super(Chain, self).__init__(**kwargs)

    def residues(self, **kwargs):
        if 'polymeric' not in kwargs:
            kwargs['polymeric'] = True
        if kwargs.get('polymeric', False) is None:
            kwargs.pop('polymeric')
        return self.__getter__(self._residues, **kwargs)

    @property
    def sequence(self):
        if self._sequence is None:
            self._sequence = [r['residue'] for r in self.residue_iterator()]
        return self._sequence

    def first(self, **kwargs):
        return self.residues(**kwargs)[0]

    def last(self, **kwargs):
        return self.residues(**kwargs)[-1]

    def endpoints(self):
        return (self.first(), self.last())

    def polymers(self):
        if not self._breaks:
            yield self
            return

        for (start, stop) in self._breaks:
            yield Chain(self._residues[start:stop], **self._data)

    def atoms(self, **kwargs):
        for residue in self.residues(**kwargs):
            for atom in residue.atoms():
                yield atom

    def __getitem__(self, index):
        return self._residues[index]

    def __repr__(self):
        return '<Chain: %s|%s|%s>' % (self.pdb, self.model, self.chain)


class Structure(Entity, EntityContainer):

    def __init__(self, residues, polymers=None, **kwargs):
        """Create a new Structure.

        :residues: A list of Components that represent the residues in this
        Structure.
        :polymers: A list of tuples of the endpoints of polymers in this
        structure.
        :kwargs: Keyword arguments to set as properites of this Structure.
        """
        self._polymers = polymers
        super(Structure, self).__init__(**kwargs)
        self.models = self.__group__(residues)

    def model(self, model):
        """Get a model by number. The number is the same as in the cif file.
        Will raise an IndexException if asking for an unknown model.

        :model: Integer for the model to get.
        :returns: A model.
        """
        # TODO: We are assuming models are sorted correctly
        return self.models[model - 1]

    def chains(self):
        """Get an iterator over all chains in in this Structure.

        :returns: An iterator for all chains.
        """
        for model in self.models:
            for chain in model.chains:
                yield chain

    def chain(self, model_number, chain_id):
        """Get a specific chain.

        :model_number: The model number to use.
        :chain_id: The chain id to get.
        :returns: A Chain or None if no chain is present.
        """
        model = self.model(model_number)
        if not model:
            return None
        return model.chain(chain_id)

    def polymers(self):
        for chain in self.chains():
            for polymer in chain.polymers():
                yield polymer

    def residues(self, **kwargs):
        """Get residues from this structure. The keyword arguments work as
        described by EntityContainer.__getter__.

        :kwargs: Keywords for filtering and ordering
        :returns: The requested residues.
        """
        for chain in self.chains():
            for residue in chain.residues(**kwargs):
                yield residue

    def infer_hydrogens(self):
        """ Infers hydrogen atoms for all bases.
        """
        for residue in self.residues():
            residue.infer_hydrogens()

    def __group__(self, residues):
        """This is a method to group a list of residues into chains and models.
        """
        mapping = col.defaultdict(lambda: col.defaultdict(list))
        for residue in residues:
            chain = residue.chain
            model = residue.model
            mapping[model][chain].append(residue)

        models = []
        for model_id, chains in mapping.items():
            model_chains = []
            for chain_id, residues in chains.items():
                chain = Chain(residues, pdb=self.pdb, model=model_id,
                              chain=chain_id)
                model_chains.append(chain)
            models.append(Model(model_chains, pdb=self.pdb, model=model_id))
        return models

    def atoms(self, **kwargs):
        for residue in self.residues():
            for atom in residue.atoms(**kwargs):
                yield atom

    def __rename__(self):
        return dict(self._data)

    def __len__(self):
        """Compute the length of this Structure. That is the number of residues
        in this structure.

        :returns: The number of atoms.
        """
        return len(self.models)

    def __repr__(self):
        return '<Structure: %s>' % self.pdb
