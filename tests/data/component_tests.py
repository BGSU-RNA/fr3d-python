import unittest as ut

import pytest
import numpy as np

from fr3d.data import Atom
from fr3d.data import Component


class BasicTest(ut.TestCase):
    def setUp(self):
        self.atoms = [
            Atom(type='C', name='a1', polymeric='A', component_number=3,
                 x=0.0, y=0.0, z=0.0),
            Atom(type='C', name='a2', polymeric='B', component_number=2,
                 x=0.0, y=0.0, z=0.0),
            Atom(type='N', name='b1', polymeric='C', component_number=1,
                 x=0.0, y=0.0, z=0.0),
            Atom(type='N', name='c2', polymeric='C', component_number=0,
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
        val = list(self.component.atoms(type='C'))
        ans = self.atoms[0:2]
        self.assertEquals(val, ans)

    def test_can_filter_using_a_list(self):
        val = list(self.component.atoms(polymeric=['A', 'C']))
        ans = [self.atoms[0], self.atoms[2], self.atoms[3]]
        self.assertEquals(val, ans)

    def test_can_filter_using_a_function(self):
        val = list(self.component.atoms(polymeric=lambda a: a == 'C'))
        ans = self.atoms[2:4]
        self.assertEquals(val, ans)

    def test_can_filter_by_several_attributes(self):
        val = list(self.component.atoms(type='C', polymeric='B'))
        ans = [self.atoms[1]]
        self.assertEquals(val, ans)

    def test_can_check_is_complete(self):
        val = self.component.is_complete(['a1', 'a2', 'b1'])
        self.assertTrue(val)

    def test_can_check_is_not_complete(self):
        val = self.component.is_complete(['a1', 'a2', 'g3'])
        self.assertFalse(val)

    def test_can_check_is_complete_using_custom_key(self):
        val = self.component.is_complete([1, 2], key='component_number')
        self.assertTrue(val)

    def test_can_check_is_not_complete_using_custom_key(self):
        val = self.component.is_complete([1, 2, 10], key='component_number')
        self.assertFalse(val)

    def test_knows_is_equal_to_itself(self):
        self.assertTrue(self.component == self.component)

    def test_knows_is_equal_to_equivelant_component(self):
        self.assertTrue(self.component == self.component.select())


class SubComponentTest(ut.TestCase):
    def setUp(self):
        self.atoms = [
            Atom(type='C', name='a1', polymeric='A',
                 x=0.0, y=0.0, z=0.0),
            Atom(type='C', name='a2', polymeric='B',
                 x=0.0, y=0.0, z=0.0),
            Atom(type='N', name='b1', polymeric='C',
                 x=0.0, y=0.0, z=0.0),
            Atom(type='N', name='c2', polymeric='C',
                 x=0.0, y=0.0, z=0.0)
        ]
        self.component = Component(self.atoms, type='rna', pdb='1GID', model=1,
                                   chain='A', sequence='C',
                                   symmetry='6_555')
        self.sub = self.component.select(name=['a1', 'a2'])

    def test_can_create_new_component_with_same_data(self):
        val = self.sub.symmetry
        ans = '6_555'
        self.assertEquals(ans, val)

    def test_can_create_new_component_with_corrent_atoms(self):
        val = list(self.sub.atoms())
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


class TransformTest(ut.TestCase):
    def setUp(self):
        atoms = [
            Atom(name='N9', x=3.0, y=3.0, z=3.0),
            Atom(name='C4', x=2.0, y=2.0, z=2.0),
            Atom(name='N3', x=1.0, y=1.0, z=1.0),
        ]
        self.residue = Component(atoms, type='rna', pdb='1GID', model=1,
                                 chain='A', sequence='C', number=50,
                                 symmetry='6_555')

    def test_can_transform_atoms(self):
        trans = np.array([[1.0, 0.0, 0.0, 0.0],
                          [0.0, -1.0, 0.0, 97.240],
                          [0.0, 0.0, -1.0, 0.0],
                          [0.0, 0.0, 0.0, 1.0]])
        residue = self.residue.transform(trans)
        val = list(list(residue.atoms())[-1].coordinates())
        ans = [1.0, 96.240, -1.0]
        self.assertEquals(ans, val)

    def test_preserves_unit_id(self):
        trans = np.array([[1.0, 0.0, 0.0, 0.0],
                          [0.0, -1.0, 0.0, 97.240],
                          [0.0, 0.0, -1.0, 0.0],
                          [0.0, 0.0, 0.0, 1.0]])
        residue = self.residue.transform(trans)
        val = residue.unit_id()
        ans = "1GID|1|A|C|50||||6_555"
        self.assertEquals(ans, val)


class InferHydrogenTest(ut.TestCase):
    def setUp(self):
        some_atoms = [
            Atom(insertion_code='?', component_id='G', name='P',
                 symmetry='1_555', component_number='118', chain='B', y=54.015,
                 x=47.242, model='1', z=51.393, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name='OP1',
                 symmetry='1_555', component_number='118', chain='B', y=53.619,
                 x=45.943, model='1', z=50.812, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name='OP2',
                 symmetry='1_555', component_number='118', chain='B', y=55.016,
                 x=48.096, model='1', z=50.668, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name="O5'",
                 symmetry='1_555', component_number='118', chain='B', y=54.52,
                 x=47.009, model='1', z=52.887, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name="C5'",
                 symmetry='1_555', component_number='118', chain='B', y=53.848,
                 x=46.12, model='1', z=53.764, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name="C4'",
                 symmetry='1_555', component_number='118', chain='B', y=54.529,
                 x=46.123, model='1', z=55.11, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name="O4'",
                 symmetry='1_555', component_number='118', chain='B', y=54.338,
                 x=47.426, model='1', z=55.73, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name="C3'",
                 symmetry='1_555', component_number='118', chain='B', y=56.037,
                 x=45.94, model='1', z=55.051, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name="O3'",
                 symmetry='1_555', component_number='118', chain='B', y=56.374,
                 x=44.553, model='1', z=54.993, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name="C2'",
                 symmetry='1_555', component_number='118', chain='B', y=56.495,
                 x=46.647, model='1', z=56.326, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name="O2'",
                 symmetry='1_555', component_number='118', chain='B', y=56.349,
                 x=45.878, model='1', z=57.494, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name="C1'",
                 symmetry='1_555', component_number='118', chain='B', y=55.528,
                 x=47.827, model='1', z=56.387, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name='N9',
                 symmetry='1_555', component_number='118', chain='B', y=56.026,
                 x=49.033, model='1', z=55.738, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name='C8',
                 symmetry='1_555', component_number='118', chain='B', y=55.785,
                 x=49.442, model='1', z=54.455, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name='N7',
                 symmetry='1_555', component_number='118', chain='B', y=56.371,
                 x=50.566, model='1', z=54.161, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name='C5',
                 symmetry='1_555', component_number='118', chain='B', y=57.031,
                 x=50.913, model='1', z=55.323, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name='C6',
                 symmetry='1_555', component_number='118', chain='B', y=57.827,
                 x=52.019, model='1', z=55.608, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name='O6',
                 symmetry='1_555', component_number='118', chain='B', y=58.131,
                 x=52.969, model='1', z=54.873, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name='N1',
                 symmetry='1_555', component_number='118', chain='B', y=58.301,
                 x=51.975, model='1', z=56.907, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name='C2',
                 symmetry='1_555', component_number='118', chain='B', y=58.033,
                 x=50.991, model='1', z=57.814, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name='N2',
                 symmetry='1_555', component_number='118', chain='B', y=58.557,
                 x=51.138, model='1', z=59.032, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name='N3',
                 symmetry='1_555', component_number='118', chain='B', y=57.299,
                 x=49.949, model='1', z=57.556, pdb='1GID'),
            Atom(insertion_code='?', component_id='G', name='C4',
                 symmetry='1_555', component_number='118', chain='B', y=56.828,
                 x=49.974, model='1', z=56.301, pdb='1GID')
        ]

        self.res = Component(some_atoms, type='rna', pdb='1GID', model=1,
                             chain='A', sequence='G', number=50,
                             symmetry='6_555')

    def test_has_no_hydrogens_initially(self):
        atoms = list(self.res.atoms(name=['H1', 'H8', 'H9', '1H2', '2H2']))
        self.assertEqual(len(atoms), 0)

    def test_hydrogen_infers_on_residue(self):
        self.res.infer_hydrogens()
        atoms = list(self.res.atoms(name=['H1', 'H8', 'H9', '1H2', '2H2']))
        self.assertEqual(len(atoms), 5)

    @pytest.mark.skip()
    def test_infers_correct_rotation_matrix_for_normal_base(self):
        pass

    @pytest.mark.skip()
    def test_infers_correct_location(self):
        pass


class AtomsWithin(ut.TestCase):
    def setUp(self):
        self.component1 = Component([
            Atom(name='N9', x=3.0, y=3.0, z=3.0),
            Atom(name='C4', x=2.0, y=2.0, z=2.0),
            Atom(name='N3', x=1.0, y=1.0, z=1.0),
            Atom(name='C3', x=0.0, y=0.0, z=0.0),
        ])
        self.component2 = Component([
            Atom(name='N1', x=0.0, y=-1.0, z=0.0),
            Atom(name='N3', x=-10.0, y=0.0, z=-1.0),
            Atom(name='C2', x=-2.0, y=0.0, z=-3.0)
        ])

    def test_knows_if_an_atom_is_within(self):
        val = self.component1.atoms_within(self.component2, 1.0)
        self.assertTrue(val)

    def test_works_with_negative_cutoff(self):
        val = self.component1.atoms_within(self.component2, -1.0)
        self.assertTrue(val)

    def test_knows_if_no_atoms_within_cutoff(self):
        val = self.component1.atoms_within(self.component2, 0.5)
        self.assertFalse(val)

    def test_knows_if_no_atoms_within_negative_cutoff(self):
        val = self.component1.atoms_within(self.component2, -0.5)
        self.assertFalse(val)

    def test_can_use_given_list_for_first_atoms(self):
        val = self.component1.atoms_within(self.component2, 1.0,
                                           using=['C3', 'N3'])
        self.assertTrue(val)

    def test_knows_if_with_given_list_nothing_is_within(self):
        val = self.component1.atoms_within(self.component2, 1.0, using=['C4', 'N3'])
        self.assertFalse(val)

    def test_can_use_to_list_to_detect_near(self):
        val = self.component1.atoms_within(self.component2, 1.0, to=['N1', 'N3'])
        self.assertTrue(val)

    def test_can_use_to_list_to_detect_not_near(self):
        val = self.component1.atoms_within(self.component2, 1.0, to=['C2', 'N3'])
        self.assertFalse(val)

    def test_can_use_both_to_detect_near(self):
        val = self.component1.atoms_within(self.component2, 1.0, to=['N1', 'N3'],
                                           using=['C3', 'N3'])
        self.assertTrue(val)

    def test_can_use_both_to_detect_not_near(self):
        val = self.component1.atoms_within(self.component2, 1.0, to=['N1', 'N3'],
                                           using=['C4', 'N3'])
        self.assertFalse(val)


class AtomTest(ut.TestCase):
    def setUp(self):
        self.atoms = [
            Atom(name='N9', x=3.0, y=3.0, z=3.0),
            Atom(name='C4', x=2.0, y=2.0, z=2.0),
            Atom(name='N3', x=1.0, y=1.0, z=1.0),
            Atom(name='C3', x=0.0, y=0.0, z=0.0),
        ]
        self.component = Component(self.atoms)

    def test_it_can_get_atoms_by_name(self):
        val = list(self.component.atoms(name=['N3', 'C4']))
        ans = [self.atoms[1], self.atoms[2]]
        self.assertEquals(ans, val)

    def test_it_can_get_atoms_by_defined_names(self):
        self.component.centers.define('example', ['N9', 'C3'])
        val = list(self.component.atoms(name='example'))
        ans = [self.atoms[0], self.atoms[-1]]
        self.assertEquals(ans, val)
