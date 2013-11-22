import unittest

from fr3d.geometry.convex_regions import totheleft
from fr3d.geometry.convex_regions import counterclockwiseinside
from fr3d.geometry.convex_regions import testcounterclockwiseconvex


class ToTheLeftTest(unittest.TestCase):
    def test_to_left_matches_left(self):
        val = totheleft([0, 1], [1, 0])
        self.assertTrue(val)

    def test_to_left_rejects_right(self):
        val = totheleft([1, 0], [0, 1])
        self.assertFalse(val)

    def test_to_left_matches_on_line(self):
        val = totheleft([0, 1], [0, 1])
        self.assertTrue(val)


class CounterClockWiseConvexTest(unittest.TestCase):
    def test_matches_counter_clockwise(self):
        val = testcounterclockwiseconvex([[0, 0], [2, 0], [3, 1], [2, 2],
                                          [0, 2]])
        self.assertTrue(val)

    def test_rejects_clockwise(self):
        val = testcounterclockwiseconvex([[0, 0], [2, 0], [3, 1], [1, 1],
                                          [2, 2], [0, 2]])
        self.assertFalse(val)

    def test_rejects_another_clockwise(self):
        val = testcounterclockwiseconvex([[0, 0], [-1, -1], [2, 0], [3, 1],
                                          [2, 2], [0, 2]])
        self.assertFalse(val)


class CounterClockWiseInsideTest(unittest.TestCase):
    def test_matches_correctly(self):
        val = counterclockwiseinside([1, 1],
                                     [[0, 0], [2, 0], [3, 1], [2, 2], [0, 2]])
        self.assertTrue(val)

    def test_rejects_correctly(self):
        val = counterclockwiseinside([6, 1],
                                     [[0, 0], [2, 0], [3, 1], [2, 2], [0, 2]])
        self.assertFalse(val)
