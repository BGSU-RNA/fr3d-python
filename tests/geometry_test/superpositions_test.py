from unittest import TestCase

from numpy import array
from numpy.testing import assert_almost_equal

#from fr3d.geometry.superpositions import bestrotation


class BestRotationTest(TestCase):

    def test_list_almost_equals(self):
        val = array([[1.0, 2.0, 3.0], [4.0, 4.1, 5.0]])
        ans = array([[1.0, 2.0, 3.0],
                     [3.999999999999999999999999999, 4.1, 5.0]])
        assert_almost_equal(val, ans)
