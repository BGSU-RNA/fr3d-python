from unittest import TestCase

import numpy as np

from fr3d.data import Atom
from fr3d.data import Component


class AtomTest(TestCase):
    def setUp(self):
        self.atom = Atom({
            'pdb': '1GID',
            'model': 1,
            'chain': 'A',
            'component_id': 'C',
            'component_number': 50,
            'name': "C1'",
            'symmetry': '6_555',
            'x': -1,
            'y': 0,
            'z': 0,
        })

    def test_has_item_access(self):
        val = self.atom['x']
        ans = -1
        self.assertEqual(val, ans)

    def test_cannot_alter_items(self):
        with self.assertRaises(TypeError):
            self.atom['x'] = 10

    def test_computes_a_unit_id(self):
        val = self.atom.unit_id()
        ans = "1GID|1|A|C|50|C1'|||6_555"
        self.assertEqual(val, ans)

    def test_computes_its_component_id(self):
        val = self.atom.component_unit_id()
        ans = "1GID|1|A|C|50||||6_555"
        self.assertEqual(val, ans)

    def test_can_get_coordinates(self):
        val = self.atom.coordinates()
        ans = np.array([-1, 0, 0])
        np.testing.assert_array_equal(val, ans)


class ComponentTest(TestCase):
    def setUp(self):
        self.atoms = [
            Atom({'type': 'C', 'name': 'A', 'number': 3}),
            Atom({'type': 'C', 'name': 'B', 'number': 2}),
            Atom({'type': 'N', 'name': 'C', 'number': 1}),
            Atom({'type': 'N', 'name': 'C', 'number': 0})
        ]
        self.component = Component({
            'type': 'rna',
            'pdb': '1GID',
            'model': 1,
            'chain': 'A',
            'sequence': 'C',
            'number': 50,
            'symmetry': '6_555',
        }, self.atoms)

    def test_has_item_access(self):
        val = self.component['type']
        ans = 'rna'
        self.assertEqual(val, ans)

    def test_computes_unit_id(self):
        val = self.component.unit_id()
        ans = "1GID|1|A|C|50||||6_555"
        self.assertEqual(val, ans)

    def test_can_get_all_atoms(self):
        val = self.component.atoms()
        ans = self.atoms
        self.assertEquals(val, ans)

    def test_can_get_filtered_atoms(self):
        val = self.component.atoms(type='C')
        ans = self.atoms[0:2]
        self.assertEquals(val, ans)

    def test_can_filter_using_a_list(self):
        val = self.component.atoms(name=['A', 'C'])
        ans = [self.atoms[0], self.atoms[2], self.atoms[3]]
        self.assertEquals(val, ans)

    def test_can_filter_using_a_function(self):
        val = self.component.atoms(name=lambda a: a == 'C')
        ans = self.atoms[2:4]
        self.assertEquals(val, ans)

    def test_can_filter_by_several_attributes(self):
        val = self.component.atoms(type='C', name='B')
        ans = [self.atoms[1]]
        self.assertEquals(val, ans)

    def test_can_get_atoms_in_order_using_order_by(self):
        val = self.component.atoms(order_by='number')
        ans = [self.atoms[3], self.atoms[2], self.atoms[1], self.atoms[0]]
        sorted_atoms = list(self.atoms)
        sorted_atoms.sort(key=lambda a: a['number'])
        self.assertEquals(val, ans)

    def test_can_get_filtered_atoms_in_order(self):
        val = self.component.atoms(type='C', order_by='number')
        ans = [self.atoms[1], self.atoms[0]]
        self.assertEquals(val, ans)


class StructureTest(TestCase):
    def setUp(self):
        pass
