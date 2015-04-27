import unittest as ut

import numpy as np

from fr3d.data import Atom


class AtomTest(ut.TestCase):
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

    def test_can_get_distance_between_atoms(self):
        val = self.atom - Atom(x=1, y=0, z=0)
        ans = 2.0
        self.assertEqual(ans, val)


class AtomTransformationTest(ut.TestCase):

    def setUp(self):
        self.atom = Atom(pdb='1GID', model=1, chain='A', component_id='C',
                         component_number=50, name="C1'", symmetry='6_555',
                         x=-1, y=0, z=2)

    def test_transform_moves_coordiantes(self):
        trans = np.array([[1.0, 0.0, 0.0, 0.0],
                          [0.0, -1.0, 0.0, 97.240],
                          [0.0, 0.0, -1.0, 0.0],
                          [0.0, 0.0, 0.0, 1.0]])
        atom = self.atom.transform(trans)
        val = [atom.x, atom.y, atom.z]
        ans = [-1.0, 97.240, -2.0]
        self.assertEquals(ans, val)

    def test_transform_preserves_unit_id(self):
        trans = np.array([[1.0, 0.0, 0.0, 0.0],
                          [0.0, -1.0, 0.0, 97.240],
                          [0.0, 0.0, -1.0, 0.0],
                          [0.0, 0.0, 0.0, 1.0]])
        atom = self.atom.transform(trans)
        val = atom.unit_id()
        ans = self.atom.unit_id()
        self.assertEquals(ans, val)


class AtomUnitIdTest(ut.TestCase):
    def test_can_build_problematic_id(self):
        atom = Atom(pdb='3V27',
                    model=1,
                    chain='A',
                    component_id='G',
                    component_number=94,
                    component_index=1,
                    insertion_code='A',
                    x=None, y=None, z=None,
                    group='ATOM',
                    type='OP1',
                    name='P',
                    symmetry='1_555',
                    polymeric=False)
        val = atom.component_unit_id()
        ans = '3V27|1|A|G|94|||A'
        self.assertEquals(ans, val)
