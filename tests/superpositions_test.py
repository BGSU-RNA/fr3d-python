from unittest import TestCase

from numpy import array
from numpy.testing import assert_almost_equal

#from fr3d.geometry.superpositions import bestrotation

from superpositions import besttransformation


#TESTS FOR THE BEST ROTATION
a = array([[ 1.0000, 0.0000, 0.0000],[ 0.0000, 1.0000, 0.0000],[ 0.0000, 0.0000, 1.0000],[-1.0000, 0.0000, 0.0000],[ 0.0000,-1.0000, 0.0000],[ 0.0000, 0.0000,-1.0000]])
b = array([[ 0.0000, 1.0000, 0.0000],[-1.0000, 0.0000, 0.0000],[ 0.0000, 0.0000, 1.0000],[ 0.0000,-1.0000, 0.0000],[ 1.0000, 0.0000, 0.0000],[ 0.0000, 0.0000,-1.0000]])
rot = array([[ 0.0000,-1.0000, 0.0000],[ 1.0000, 0.0000, 0.0000],[ 0.0000, 0.0000, 1.0000]])
rotation, newA, newB = besttransformation(a,b)
print rotation
assert_almost_equal(rot,rotation)
print 'Successful Test!'


#TESTS FOR THE BEST TRANSFORMATION
a = array([[ 1.0000, 0.0000, 0.0000],[ 0.0000, 1.0000, 0.0000],[ 0.0000, 0.0000, 1.0000]])
b = array([[ 0.0000, 1.0000, 0.0000],[-1.0000, 0.0000, 0.0000],[ 0.0000, 0.0000, 1.0000]])
rot = array([[-0.0000,-1.0000, 0.0000],[ 1.0000, 0.0000, 0.0000],[ 0.0000, 0.0000, 1.0000]])
rotation, newA,newB = besttransformation(a,b)
assert_almost_equal(rot,rotation)

print rotation
print 'Successful Test!'


#class BestRotationTest(TestCase):
 #   def test_list_almost_equals(self):
  #      a = array([[ 1.0000, 0.0000, 0.0000],[ 0.0000, 1.0000, 0.0000],[ 0.0000, 0.0000, 1.0000],[-1.0000, 0.0000, 0.0000],[ 0.0000,-1.0000, 0.0000],[ 0.0000, 0.0000,-1.0000]])
   #     b = array([[ 0.0000, 1.0000, 0.0000],[-1.0000, 0.0000, 0.0000],[ 0.0000, 0.0000, 1.0000],[ 0.0000,-1.0000, 0.0000],[ 1.0000, 0.0000, 0.0000],[ 0.0000, 0.0000,-1.0000]])
    #    rot = array([[ 0.0000,-1.0000, 0.0000],[ 1.0000, 0.0000, 0.0000],[ 0.0000, 0.0000, 1.0000]])
     #   rotation, shift = besttransformation(a,b)
      #  assert_almost_equal(rot,rotation)


#class BestTransformationTest(TestCase):
 #   def test_list_almost_equals(self):
#		a = array([[ 1.0000, 0.0000, 0.0000],[ 0.0000, 1.0000, 0.0000],[ 0.0000, 0.0000, 1.0000]])
#		b = array([[ 0.0000, 1.0000, 0.0000],[-1.0000, 0.0000, 0.0000],[ 0.0000, 0.0000, 1.0000]])
#		rot = array([[-0.0000,-1.0000, 0.0000],[ 1.0000, 0.0000, 0.0000],[ 0.0000, 0.0000, 1.0000]])
#		rotation, shift = besttransformation(a,b)
 #       assert_almost_equal(rot,rotation)

#class BestRotationTest(TestCase):

 #   def test_list_almost_equals(self):
  #      val = array([[1.0, 2.0, 3.0], [4.0, 4.1, 5.0]])
   #     ans = array([[1.0, 2.0, 3.0],
    #                 [3.999999999999999999999999999, 4.1, 5.0]])
     #   assert_almost_equal(val, ans)
#