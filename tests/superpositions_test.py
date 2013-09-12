from unittest import TestCase
import numpy
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
        ans = array([[0.0000, 1.0000, 0.0000],
                     [-1.0000, 0.0000, 0.0000],
                     [0.0000, 0.0000, 1.0000]])
        assert_almost_equal(ans, self.transformation)

    def test_can_transform_a_3_x_3(self):
        a = array([[1.0000, 0.0000, 0.0000],
                   [0.0000, 1.0000, 0.0000],
                   [0.0000, 0.0000, 1.0000]])
        b = array([[0.0000, 1.0000, 0.0000],
                   [-1.0000, 0.0000, 0.0000],
                   [0.0000, 0.0000, 1.0000]])
        ans = array([[0.0000, 1.0000, 0.0000],
                     [-1.0000, 0.0000, 0.0000],
                     [0.0000, 0.0000, 1.0000]])
        rotation, _, _ = besttransformation(a, b)
        assert_almost_equal(ans, rotation)

    def test_can_transform_a_4_x_3_45degreexrotation(self):
        theta=numpy.pi/4        
        ans = array([[1.0000, 0.0000, 0.0000],
                     [0.0000, numpy.cos(theta), -1*numpy.cos(theta)],
                     [0.0000, numpy.sin(theta), numpy.cos(theta)]])
        a = array([[3.0000, 7.0000, 1.0000],
                   [5.0000, 9.0000, 2.0000],
                   [7.0000, 4.0000, 6.0000],
                   [1.0000, 1.0000, 1.0000]])
        J = numpy.ones((4, 4))
        meana=numpy.dot(J,a)/4
        b=numpy.dot(a-meana,ans)
        rotation, _, _ = besttransformation(a, b)
        assert_almost_equal(ans, rotation)
        
    def test_can_transform_a_4_x_3_45degreeyrotation(self):
        theta=numpy.pi/4        
        ans = array([[numpy.cos(theta), 0.0000, numpy.sin(theta)],
                     [0.0000, 1.0000, 0.0000],
                     [-1*numpy.sin(theta), 0.0000, numpy.cos(theta)]])
        a = array([[3.0000, 7.0000, 1.0000],
                   [5.0000, 9.0000, 2.0000],
                   [7.0000, 4.0000, 6.0000],
                   [1.0000, 1.0000, 1.0000]])
        J = numpy.ones((4, 4))
        meana=numpy.dot(J,a)/4
        b=numpy.dot(a-meana,ans)
        rotation, _, _ = besttransformation(a, b)
        assert_almost_equal(ans, rotation)

    def test_can_transform_a_4_x_3_45degreezrotation(self):
        theta=numpy.pi/4        
        ans = array([[numpy.cos(theta), -1*numpy.sin(theta), 0.0000],
                     [numpy.sin(theta), numpy.cos(theta), 0.0000],
                     [0.0000, 0.0000, 1.0000]])
        a = array([[3.0000, 7.0000, 1.0000],
                   [5.0000, 9.0000, 2.0000],
                   [7.0000, 4.0000, 6.0000],
                   [1.0000, 1.0000, 1.0000]])
        J = numpy.ones((4, 4))
        meana=numpy.dot(J,a)/4
        b=numpy.dot(a-meana,ans)
        rotation, _, _ = besttransformation(a, b)
        assert_almost_equal(ans, rotation)

    def test_can_transform_a_4_x_3_generalrotation(self):
        theta1=numpy.pi/4
        theta2=numpy.pi/3
        theta3=numpy.pi/9
        R1 = array([[1.0000, 0.0000, 0.0000],
                     [0.0000, numpy.cos(theta1), -1*numpy.cos(theta1)],
                     [0.0000, numpy.sin(theta1), numpy.cos(theta1)]])
        R2 = array([[numpy.cos(theta2), 0.0000, numpy.sin(theta2)],
                     [0.0000, 1.0000, 0.0000],
                     [-1*numpy.sin(theta2), 0.0000, numpy.cos(theta2)]])
        R3 = array([[numpy.cos(theta3), -1*numpy.sin(theta3), 0.0000],
                     [numpy.sin(theta3), numpy.cos(theta3), 0.0000],
                     [0.0000, 0.0000, 1.0000]])
        ans=numpy.dot(numpy.dot(R1,R2),R3)
        a = array([[3.0000, 7.0000, 1.0000],
                   [5.0000, 9.0000, 2.0000],
                   [7.0000, 4.0000, 6.0000],
                   [1.0000, 1.0000, 1.0000]])
        J = numpy.ones((4, 4))
        meana=numpy.dot(J,a)/4
        b=numpy.dot(a-meana,ans)
        rotation, _, _ = besttransformation(a, b)
        assert_almost_equal(ans, rotation)