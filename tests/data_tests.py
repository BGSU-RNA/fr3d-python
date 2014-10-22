from unittest import skip
from unittest import TestCase

import numpy as np

from fr3d.data import Atom
from fr3d.data import AtomProxy
from fr3d.data import Component
from fr3d.data import Structure


class AtomProxyTest(TestCase):
    def setUp(self):
        self.atoms = [
            Atom(type='C', name='a1', type_name='A', number=3,
                 x=1.0, y=0.0, z=0.0),
            Atom(type='C', name='a2', type_name='B', number=2,
                 x=2.0, y=0.0, z=0.0),
            Atom(type='N', name='b1', type_name='C', number=1,
                 x=3.0, y=0.0, z=0.0),
            Atom(type='N', name='c2', type_name='C', number=0,
                 x=0.0, y=1.0, z=0.0)
        ]
        self.proxy = AtomProxy(self.atoms)

    def test_can_set_a_value(self):
        ans = 10
        self.proxy['bob'] = ans
        val = self.proxy['bob']
        self.assertEqual(val, ans)

    def test_can_get_atom_value(self):
        ans = np.array([0.0, 1.0, 0.0])
        val = self.proxy['c2']
        np.testing.assert_almost_equal(ans, val)

    def test_knows_is_missing_a_value(self):
        self.assertFalse('bob' in self.proxy)

    def test_can_test_has_value(self):
        self.proxy['bob'] = 1
        self.assertTrue('bob' in self.proxy)

    def test_knows_has_atom(self):
        self.assertTrue('a1' in self.proxy)

    def test_lets_override_proxy_lookup(self):
        ans = 'a'
        self.proxy['a1'] = ans
        val = self.proxy['a1']
        self.assertEqual(val, ans)

    def test_length_counts_atoms(self):
        self.proxy['steve'] = 3
        ans = 5
        val = len(self.proxy)
        self.assertEqual(val, ans)

    def test_can_find_new_atom_positions(self):
        self.atoms.append(Atom(name='s', x=3.0, y=2.0, z=1.0))
        val = self.proxy['s']
        ans = np.array([3.0, 2.0, 1.0])
        np.testing.assert_almost_equal(ans, val)


class AtomTest(TestCase):
    def setUp(self):
        self.atom = Atom(pdb='1GID', model=1, chain='A', component_id='C',
                         component_number=50, name="C1'", symmetry='6_555',
                         x=-1, y=0, z=0)

    def test_has_item_access(self):
        val = self.atom.x
        ans = -1
        self.assertEqual(val, ans)

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

    def test_can_get_atoms_in_order_using_order_by(self):
        val = self.component.atoms(order_by='number')
        ans = [self.atoms[3], self.atoms[2], self.atoms[1], self.atoms[0]]
        sorted_atoms = list(self.atoms)
        sorted_atoms.sort(key=lambda a: a.number)
        self.assertEquals(val, ans)

    def test_can_get_filtered_atoms_in_order(self):
        val = self.component.atoms(type='C', order_by='number')
        ans = [self.atoms[1], self.atoms[0]]
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


class SubComponentTest(TestCase):
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


class InferHydrogenTest(TestCase):
    def setUp(self):
        some_atoms = [
            Atom(ins_code='?', component_id='G', name='P', symmetry='1_555',
                 component_number='118', chain='B', y=54.015, x=47.242,
                 model='1', z=51.393, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name='OP1', symmetry='1_555',
                 component_number='118', chain='B', y=53.619, x=45.943,
                 model='1', z=50.812, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name='OP2', symmetry='1_555',
                 component_number='118', chain='B', y=55.016, x=48.096,
                 model='1', z=50.668, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name="O5'", symmetry='1_555',
                 component_number='118', chain='B', y=54.52, x=47.009,
                 model='1', z=52.887, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name="C5'", symmetry='1_555',
                 component_number='118', chain='B', y=53.848, x=46.12,
                 model='1', z=53.764, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name="C4'", symmetry='1_555',
                 component_number='118', chain='B', y=54.529, x=46.123,
                 model='1', z=55.11, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name="O4'", symmetry='1_555',
                 component_number='118', chain='B', y=54.338, x=47.426,
                 model='1', z=55.73, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name="C3'", symmetry='1_555',
                 component_number='118', chain='B', y=56.037, x=45.94,
                 model='1', z=55.051, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name="O3'", symmetry='1_555',
                 component_number='118', chain='B', y=56.374, x=44.553,
                 model='1', z=54.993, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name="C2'", symmetry='1_555',
                 component_number='118', chain='B', y=56.495, x=46.647,
                 model='1', z=56.326, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name="O2'", symmetry='1_555',
                 component_number='118', chain='B', y=56.349, x=45.878,
                 model='1', z=57.494, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name="C1'", symmetry='1_555',
                 component_number='118', chain='B', y=55.528, x=47.827,
                 model='1', z=56.387, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name='N9', symmetry='1_555',
                 component_number='118', chain='B', y=56.026, x=49.033,
                 model='1', z=55.738, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name='C8', symmetry='1_555',
                 component_number='118', chain='B', y=55.785, x=49.442,
                 model='1', z=54.455, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name='N7', symmetry='1_555',
                 component_number='118', chain='B', y=56.371, x=50.566,
                 model='1', z=54.161, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name='C5', symmetry='1_555',
                 component_number='118', chain='B', y=57.031, x=50.913,
                 model='1', z=55.323, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name='C6', symmetry='1_555',
                 component_number='118', chain='B', y=57.827, x=52.019,
                 model='1', z=55.608, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name='O6', symmetry='1_555',
                 component_number='118', chain='B', y=58.131, x=52.969,
                 model='1', z=54.873, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name='N1', symmetry='1_555',
                 component_number='118', chain='B', y=58.301, x=51.975,
                 model='1', z=56.907, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name='C2', symmetry='1_555',
                 component_number='118', chain='B', y=58.033, x=50.991,
                 model='1', z=57.814, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name='N2', symmetry='1_555',
                 component_number='118', chain='B', y=58.557, x=51.138,
                 model='1', z=59.032, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name='N3', symmetry='1_555',
                 component_number='118', chain='B', y=57.299, x=49.949,
                 model='1', z=57.556, pdb='1GID'),
            Atom(ins_code='?', component_id='G', name='C4', symmetry='1_555',
                 component_number='118', chain='B', y=56.828, x=49.974,
                 model='1', z=56.301, pdb='1GID')
        ]

        self.res = Component(some_atoms, type='rna', pdb='1GID', model=1,
                             chain='A', sequence='G', number=50,
                             symmetry='6_555')

    def test_hydrogen_infers_on_residue(self):
        self.res.infer_hydrogens()
        atoms = list(self.res.atoms(name=['H1', 'H8', 'H9', '1H2', '2H2']))
        self.assertEqual(len(atoms), 5)


class ComponentCenterTest(TestCase):
    def setUp(self):
        atoms = [
            Atom(name='N9', x=3.0, y=3.0, z=3.0),
            Atom(name='C4', x=2.0, y=2.0, z=2.0),
            Atom(name='N3', x=1.0, y=1.0, z=1.0),
        ]
        self.residue = Component(atoms, sequence='A')
        print(self.residue)

    def test_computes_center_correctly(self):
        val = self.residue.centers['base']
        ans = np.array([2.0, 2.0, 2.0])
        np.testing.assert_array_equal(val, ans)


class StructureTest(TestCase):
    def setUp(self):
        self.structure = Structure([], pdb='1S72', model=1)

    @skip
    def test_generates_a_unit_id(self):
        val = self.structure.unit_id()
        ans = '1S72|1'
        self.assertEqual(val, ans)
