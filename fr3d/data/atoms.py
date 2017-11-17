"""This module contains classes for representing Atoms from structure files.
"""

from fr3d.unit_ids import encode

import numpy as np


class Atom(object):
    """This class represents atoms in a structure. It provides a simple dict
    like access for data as well as a way to get its coordinates, unit id
    and the unit id of the component it belongs to.
    """

    def __init__(self, pdb=None, model=None, chain=None,
                 component_id=None, component_number=None,
                 component_index=None, insertion_code=None, alt_id=None,
                 x=None, y=None, z=None, group=None, type=None, name=None,
                 symmetry=None, polymeric=None):

        """Create a new Atom.

        :param string pdb: The pdb id this atom is a part of.
        :param int model: The model this atom is a part of.
        :param string chain: The chain this atom is a part of.
        """

        self.pdb = pdb
        self.model = model
        self.chain = chain
        self.component_id = component_id
        self.component_number = component_number
        self.component_index = component_index
        self.insertion_code = insertion_code
        self.alt_id = alt_id
        self.x = x
        self.y = y
        self.z = z
        self.group = group
        self.type = type
        self.name = name
        self.symmetry = symmetry
        self.polymeric = polymeric

    def component_unit_id(self):
        """Generate the unit id of the component this atom belongs to.

        :returns: A string of the unit id for this atom's component.
        """

        return encode({
            'pdb': self.pdb,
            'model': self.model,
            'chain': self.chain,
            'component_id': self.component_id,
            'component_number': self.component_number,
            'alt_id': self.alt_id,
            'insertion_code': self.insertion_code,
            'symmetry': self.symmetry
        })

    def unit_id(self):
        """Create the unit id for this Atom.

        :returns: The unit id string.
        """
        return encode({
            'pdb': self.pdb,
            'model': self.model,
            'chain': self.chain,
            'component_id': self.component_id,
            'component_number': self.component_number,
            'atom_name': self.name,
            'alt_id': self.alt_id,
            'insertion_code': self.insertion_code,
            'symmetry': self.symmetry
        })

    def transform(self, transform):
        """Create a new atom based of this one, but with transformed
        coordinates.

        :transform: A 4x4 numpy array that is the transformation matrix.
        :returns: A new Atom representating as a result of transforming this
        ones coordiantes.
        """

        result = np.dot(transform, np.array([self.x, self.y, self.z, 1.0]))
        x, y, z = result[0:3].T
        return Atom(x=x, y=y, z=z,
                    pdb=self.pdb,
                    model=self.model,
                    chain=self.chain,
                    component_id=self.component_id,
                    component_number=self.component_number,
                    component_index=self.component_index,
                    insertion_code=self.insertion_code,
                    alt_id=self.alt_id,
                    group=self.group,
                    type=self.type,
                    name=self.name,
                    symmetry=self.symmetry,
                    polymeric=self.polymeric)

    def coordinates(self):
        """Return a numpy array of the x, y, z coordinates for this atom.

        :returns: A numpy array of the x, y, z coordinates.
        """
        return np.array([self.x, self.y, self.z])

    def distance(self, atom):
        """Compute the distance between this atom and another atom.

        :atom: Another atom.
        :returns: The distance.
        """
        return np.linalg.norm(self.coordinates() - atom.coordinates())

    def __repr__(self):
        """Creates the string used to represent this Atom when printing it.

        :returns: The string representation.
        """
        return '<Atom: %s>' % self.unit_id()
