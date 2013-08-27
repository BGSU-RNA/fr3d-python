from unittest import TestCase

from numpy import array
from numpy.testing import assert_almost_equal

from fr3d.geometry.superpositions import besttransformation


class TransformationTest(TestCase):

    def setUp(self):
        a = array([[1.0000, 0.0000, 0.0000],
                   [0.0000, 1.0000, 0.0000],
                   [0.0000, 0.0000, 1.0000],
                   [-1.0000, 0.0000, 0.0000],
                   [0.0000, -1.0000, 0.0000],
                   [0.0000, 0.0000, -1.0000]])
        b = array([[0.0000, 1.0000, 0.0000],
                   [-1.0000, 0.0000, 0.0000],
                   [0.0000, 0.0000, 1.0000],
                   [0.0000, -1.0000, 0.0000],
                   [1.0000, 0.0000, 0.0000],
                   [0.0000, 0.0000, -1.0000]])
        self.transformation, _, _ = besttransformation(a, b)

    def test_can_transform_a_3_x_6(self):
        ans = array([[0.0000, -1.0000, 0.0000],
                     [1.0000, 0.0000, 0.0000],
                     [0.0000, 0.0000, 1.0000]])
        assert_almost_equal(ans, self.transformation)

    def test_can_transform_a_3_x_3(self):
        a = array([[1.0000, 0.0000, 0.0000],
                   [0.0000, 1.0000, 0.0000],
                   [0.0000, 0.0000, 1.0000]])
        b = array([[0.0000, 1.0000, 0.0000],
                   [-1.0000, 0.0000, 0.0000],
                   [0.0000, 0.0000, 1.0000]])
        ans = array([[-0.0000, -1.0000, 0.0000],
                     [1.0000, 0.0000, 0.0000],
                     [0.0000, 0.0000, 1.0000]])
        rotation, _, _ = besttransformation(a, b)
        assert_almost_equal(ans, rotation)
