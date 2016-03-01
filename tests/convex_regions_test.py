import unittest
from nose import SkipTest

from fr3d.geometry import convex_regions as cr


class ToTheLeftTest(unittest.TestCase):
    def test_to_left_matches_left(self):
        val = cr.totheleft([0, 1], [1, 0])
        self.assertTrue(val)

    def test_to_left_rejects_right(self):
        val = cr.totheleft([1, 0], [0, 1])
        self.assertFalse(val)

    def test_to_left_matches_on_line(self):
        val = cr.totheleft([0, 1], [0, 1])
        self.assertTrue(val)


class CounterClockWiseConvexTest(unittest.TestCase):
    def test_matches_counter_clockwise(self):
        val = cr.testcounterclockwiseconvex([[0, 0], [2, 0], [3, 1], [2, 2],
                                             [0, 2]])
        self.assertTrue(val)

    def test_rejects_clockwise(self):
        val = cr.testcounterclockwiseconvex([[0, 0], [2, 0], [3, 1], [1, 1],
                                             [2, 2], [0, 2]])
        self.assertFalse(val)

    def test_rejects_another_clockwise(self):
        val = cr.testcounterclockwiseconvex([[0, 0], [-1, -1], [2, 1], [3, 1],
                                             [2, 2], [0, 2]])
        self.assertFalse(val)


class CounterClockWiseInsideTest(unittest.TestCase):
    def test_matches_correctly(self):
        raise SkipTest()
        val = cr.counterclockwiseinside([1, 1], [[0, 0], [2, 0], [3, 1],
                                                 [2, 2], [0, 2]])
        self.assertTrue(val)

    def test_rejects_correctly(self):
        raise SkipTest()
        val = cr.counterclockwiseinside([6, 1], [[0, 0], [2, 0], [3, 1],
                                                 [2, 2], [0, 2]])
        self.assertFalse(val)
