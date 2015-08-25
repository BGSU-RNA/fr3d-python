import unittest as ut

from fr3d.data.pairs import Pairs
from fr3d.data import Structure
from fr3d.data import Component
from fr3d.data import Atom


class PairsTest(ut.TestCase):
    def setUp(self):
        self.nt1 = Component(atoms=[
            Atom(x=0, y=0, z=1, name='C1'),
            Atom(x=0, y=0, z=0, name='C2'),
            Atom(x=1, y=0, z=1, name='N'),
        ], sequence='A', model=1, chain='A', number=1, polymeric=True)

        self.nt2 = Component(atoms=[
            Atom(x=5, y=5, z=5, name='C1'),
            Atom(x=5, y=5, z=5, name='C2'),
            Atom(x=5, y=5, z=5, name='N'),
        ], sequence='U', model=1, chain='A', number=2, polymeric=True)

        self.nt3 = Component(atoms=[
            Atom(x=-5, y=-5, z=-5, name='C1'),
            Atom(x=-5, y=-5, z=-5, name='C2'),
            Atom(x=-5, y=-5, z=-5, name='N'),
        ], sequence='C', model=1, chain='A', number=3, polymeric=True)

        self.nt4 = Component(atoms=[
            Atom(x=0, y=0, z=3, name='C1'),
            Atom(x=0, y=0, z=3, name='C2'),
            Atom(x=0, y=0, z=3, name='N'),
        ], sequence='U', model=1, chain='A', number=4, polymeric=True)

        self.nt5 = Component(atoms=[
            Atom(x=0, y=0, z=8, name='C1'),
            Atom(x=0, y=0, z=8, name='C2'),
            Atom(x=0, y=0, z=3, name='N'),
        ], sequence='G', model=1, chain='A', number=5, polymeric=True)

        structure = Structure([self.nt1, self.nt2, self.nt3, self.nt4,
                               self.nt5], pdb='0000', model='1')
        self.pairs = Pairs(structure)

    def test_can_get_all_pairs_by_attributes(self):
        self.pairs.first(sequence='A')
        self.pairs.second(sequence='U')
        val = list(self.pairs)
        ans = [(self.nt1, self.nt2), (self.nt1, self.nt4)]
        self.assertEquals(ans, val)

    def test_can_get_using_center_distance(self):
        self.pairs.distance(use='center', cutoff=3.0)
        val = list(self.pairs)
        ans = [(self.nt1, self.nt4), (self.nt4, self.nt1)]
        self.assertEquals(ans, val)

    def test_can_get_using_atom_atom_distance(self):
        self.pairs.distance(use='atoms', cutoff=3.0)
        val = list(self.pairs)
        ans = [(self.nt1, self.nt4), (self.nt1, self.nt5),
               (self.nt4, self.nt1), (self.nt4, self.nt5),
               (self.nt5, self.nt1), (self.nt5, self.nt4)]
        self.assertEquals(ans, val)

    def test_can_get_using_specific_atoms(self):
        self.pairs.distance(use='atoms', cutoff=3.0, first_atoms=['C1', 'C2'],
                            second_atoms=['C2'])
        val = list(self.pairs)
        ans = [(self.nt1, self.nt4), (self.nt4, self.nt1)]
        self.assertEquals(ans, val)

    def test_can_get_using_atoms_by_definition_name(self):
        self.pairs.distance(use='atoms', cutoff=4.0, first_atoms='base',
                            second_atoms='base')
        val = list(self.pairs)
        ans = [(self.nt1, self.nt4), (self.nt4, self.nt1)]
        self.assertEquals(ans, val)

    def test_can_get_using_definition_name(self):
        self.pairs.distance(use='center', cutoff=3.0, first_atoms='base')
        val = list(self.pairs)
        ans = [(self.nt1, self.nt4), (self.nt4, self.nt1)]
        self.assertEquals(ans, val)
