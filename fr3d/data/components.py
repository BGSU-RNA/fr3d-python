from fr3d.data.base import EntitySelector
from fr3d.data.base import AtomProxy
from fr3d.data.atoms import Atom
from fr3d import definitions as defs
from fr3d.geometry.superpositions import besttransformation

import numpy as np

from fr3d.unit_ids import encode


class Component(EntitySelector):
    """This represents things like nucleic acids, amino acids, small molecules
    and ligands.
    """

    def __init__(self, atoms, pdb=None, model=None, type=None, chain=None,
                 symmetry=None, sequence=None, number=None, index=None,
                 insertion_code=None, polymeric=None, alt_id=None):
        """Create a new Component.

        :atoms: The atoms this component is composed of.
        :pdb: The pdb this is a part of.
        :model: The model number.
        """

        self._atoms = atoms
        self.pdb = pdb
        self.model = model
        self.type = type
        self.chain = chain
        self.symmetry = symmetry
        self.sequence = sequence
        self.number = number
        self.index = index
        self.insertion_code = insertion_code
        self.polymeric = polymeric
        self.alt_id = alt_id
        self.centers = AtomProxy(self._atoms)

        if self.sequence in defs.RNAbaseheavyatoms:
            atoms = defs.RNAbaseheavyatoms[self.sequence]
            self.centers.define('base', atoms)

        if self.sequence in defs.nt_backbone:
            atoms = defs.nt_backbone[self.sequence]
            self.centers.define('nt_backbone', atoms)

        if self.sequence in defs.aa_fg:
            atoms = defs.aa_fg[self.sequence]
            self.centers.define('aa_fg', atoms)

        if self.sequence in defs.aa_backbone:
            atoms = defs.aa_backbone[self.sequence]
            self.centers.define('aa_backbone', atoms)

    def atoms(self, **kwargs):
        """Get, filter and sort the atoms in this component. Access is as
        described by EntitySelector.

        :kwargs: The keyword arguments to filter and sort by.
        :returns: A list of the requested atoms.
        """

        name = kwargs.get('name')
        if isinstance(name, basestring):
            definition = self.centers.definition(name)
            if definition:
                kwargs['name'] = definition

        return EntitySelector(self._atoms, **kwargs)

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
        return Component(list(self.atoms(**kwargs)),
                         pdb=self.pdb,
                         model=self.model,
                         type=self.type,
                         chain=self.chain,
                         symmetry=self.symmetry,
                         sequence=self.sequence,
                         number=self.number,
                         index=self.index,
                         insertion_code=self.insertion_code,
                         alt_id=self.alt_id,
                         polymeric=self.polymeric)

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
        found = list(self.atoms(**kwargs))
        return len(found) == len(names)

    def infer_hydrogens(self):
        """Infer the coordinates of the hydrogen atoms of this component.
        Currently, it only works for RNA with .sequence
        """
        if self.sequence not in defs.RNAbaseheavyatoms:
            return None
        R = []
        S = []
        baseheavy = defs.RNAbaseheavyatoms[self.sequence]

        for atom in self.atoms(name=baseheavy):
            coordinates = atom.coordinates()
            R.append(coordinates)
            S.append(defs.RNAbasecoordinates[self.sequence][atom.name])

        R = np.array(R)
        R = R.astype(np.float)
        S = np.array(S)
        try:
            rotation_matrix, fitted, base_center, rmsd, sse = \
                besttransformation(R, S)
        except:
            return None

        self.rotation_matrix = rotation_matrix

        hydrogens = defs.RNAbasehydrogens[self.sequence]
        coordinates = defs.RNAbasecoordinates[self.sequence]

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
        matrix. This does not keep the rotation matrix if any, but will keep
        any added hydrogens.

        :transform: The transformation matrix to apply.
        :returns: A new Component with the same properties by rotated atoms.
        """

        atoms = [atom.transform(transform) for atom in self.atoms()]
        return Component(atoms, pdb=self.pdb,
                         model=self.model,
                         type=self.type,
                         chain=self.chain,
                         symmetry=self.symmetry,
                         sequence=self.sequence,
                         number=self.number,
                         index=self.index,
                         insertion_code=self.insertion_code,
                         alt_id=self.alt_id,
                         polymeric=self.polymeric)

    def unit_id(self):
        """Compute the unit id of this Component.

        :returns: The unit id.
        """

        return encode({
            'pdb': self.pdb,
            'model': self.model,
            'chain': self.chain,
            'component_id': self.sequence,
            'component_number': self.number,
            'alt_id': self.alt_id,
            'insertion_code': self.insertion_code,
            'symmetry': self.symmetry
        })

    def atoms_within(self, other, using=None, to=None, cutoff=4.0):
        """Determine if there are any atoms from another component within some
        distance.

        :other: Another component to compare agains.
        :using: The atoms from this component to compare with.
        :to: The atoms from the other component to compare against.
        :cutoff: The distances atoms must be within. Default 4.0
        """

        kw1 = {}
        if using:
            kw1['name'] = using

        kw2 = {}
        if to:
            kw2['name'] = to

        for atom1 in self.atoms(**kw1):
            for atom2 in other.atoms(**kw2):
                if atom1.distance(atom2) <= abs(cutoff):
                    return True
        return False

    def distance(self, other, using='*', to='*'):
        """Compute a center center distance between this and another component.

        :other: The other component to get distance to.
        :using: A list of atom names to use for this component. Defaults to '*'
        meaning all atoms.
        :to: A list of atoms names for the second component. Defaults to '*'
        meaning all atoms.
        :returns: The distance between the two centers.
        """
        coordinates = self.centers[using]
        other_coord = other.centers[to]
        distance = np.subtract(coordinates, other_coord)
        return np.linalg.norm(distance)

    def __len__(self):
        """Compute the length of this Component. This is the number of atoms in
        this residue.

        :returns: The number of atoms.
        """
        return len(self._atoms)

    def __eq__(self, other):
        return isinstance(other, Component) and \
            self.pdb == other.pdb and \
            self.model == other.model and \
            self.chain == other.chain and \
            self.symmetry == other.symmetry and \
            self.sequence == other.sequence and \
            self.number == other.number and \
            self.insertion_code == other.insertion_code and \
            self.alt_id == other.alt_id

    def __repr__(self):
        return '<Component %s>' % self.unit_id()
