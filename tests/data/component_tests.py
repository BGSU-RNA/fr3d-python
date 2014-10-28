import unittest as ut

import numpy as np

from fr3d.data import Atom
from fr3d.data import Component


class BasicTest(ut.TestCase):
    def setUp(self):
        self.atoms = [
            Atom(type='C', name='a1', type_name='A', number=3,
                 x=0.0, y=0.0, z=0.0),
            Atom(type='C', name='a2', type_name='B', number=2,
                 x=0.0, y=0.0, z=0.0),
            Atom(type='N', name='b1', type_name='C', number=1,
                 x=0.0, y=0.0, z=0.0),
            Atom(type='N', name='c2', type_name='C', number=0,
                 x=0.0, y=0.0, z=0.0)
        ]
        self.component = Component(self.atoms, type='rna', pdb='1GID', model=1,
                                   chain='A', sequence='C', number=50,
                                   symmetry='6_555')

    def test_has_item_access(self):
        val = self.component.type
        ans = 'rna'
        self.assertEqual(val, ans)

    def test_length_is_number_atoms(self):
        val = len(self.component)
        ans = 4
        self.assertEqual(val, ans)

    def test_computes_unit_id(self):
        val = self.component.unit_id()
        ans = "1GID|1|A|C|50||||6_555"
        self.assertEquals(val, ans)

    def test_can_get_filtered_atoms(self):
        val = self.component.atoms(type='C')
        ans = self.atoms[0:2]
        self.assertEquals(val, ans)

    def test_can_filter_using_a_list(self):
        val = self.component.atoms(type_name=['A', 'C'])
        ans = [self.atoms[0], self.atoms[2], self.atoms[3]]
        self.assertEquals(val, ans)

    def test_can_filter_using_a_function(self):
        val = self.component.atoms(type_name=lambda a: a == 'C')
        ans = self.atoms[2:4]
        self.assertEquals(val, ans)

    def test_can_filter_by_several_attributes(self):
        val = self.component.atoms(type='C', type_name='B')
        ans = [self.atoms[1]]
        self.assertEquals(val, ans)

    def test_can_check_is_complete(self):
        val = self.component.is_complete(['a1', 'a2', 'b1'])
        self.assertTrue(val)

    def test_can_check_is_not_complete(self):
        val = self.component.is_complete(['a1', 'a2', 'g3'])
        self.assertFalse(val)

    def test_can_check_is_complete_using_custom_key(self):
        val = self.component.is_complete([1, 2], key='number')
        self.assertTrue(val)

    def test_can_check_is_not_complete_using_custom_key(self):
        val = self.component.is_complete([1, 2, 10], key='number')
        self.assertFalse(val)


class SubComponentTest(ut.TestCase):
    def setUp(self):
        self.atoms = [
            Atom(type='C', name='a1', type_name='A', number=3,
                 x=0.0, y=0.0, z=0.0),
            Atom(type='C', name='a2', type_name='B', number=2,
                 x=0.0, y=0.0, z=0.0),
            Atom(type='N', name='b1', type_name='C', number=1,
                 x=0.0, y=0.0, z=0.0),
            Atom(type='N', name='c2', type_name='C', number=0,
                 x=0.0, y=0.0, z=0.0)
        ]
        self.component = Component(self.atoms, type='rna', pdb='1GID', model=1,
                                   chain='A', sequence='C', number=50,
                                   symmetry='6_555')
        self.sub = self.component.select(name=['a1', 'a2'])

    def test_can_create_new_component_with_same_data(self):
        val = self.sub.symmetry
        ans = '6_555'
        self.assertEquals(ans, val)

    def test_can_create_new_component_with_corrent_atoms(self):
        val = self.sub.atoms()
        ans = self.atoms[0:2]
        self.assertEquals(ans, val)


class CenterTest(ut.TestCase):
    def setUp(self):
        atoms = [
            Atom(name='N9', x=3.0, y=3.0, z=3.0),
            Atom(name='C4', x=2.0, y=2.0, z=2.0),
            Atom(name='N3', x=1.0, y=1.0, z=1.0),
        ]
        self.residue = Component(atoms, sequence='A')

    def test_computes_center_correctly(self):
        val = self.residue.centers['base']
        ans = np.array([2.0, 2.0, 2.0])
        np.testing.assert_array_equal(val, ans)
