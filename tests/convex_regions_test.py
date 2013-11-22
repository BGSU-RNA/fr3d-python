# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 11:12:22 2013

@author: zirbel
"""

from fr3d.geometry.convex_regions import totheleft

v = totheleft([0, 1],[1, 0])
print v
# should be true

v = totheleft([1, 0], [0, 1])
print v
# should be false

v = totheleft([0, 1],[0, 1])
print v
# should be true, being on the line is OK

from fr3d.geometry.convex_regions import testcounterclockwiseconvex

v = testcounterclockwiseconvex([[0,0],[2,0],[3,1],[2,2],[0,2]])
# should be true

v = testcounterclockwiseconvex([[0,0],[2,0],[3,1],[1,1],[2,2],[0,2]])
# should be false

v = testcounterclockwiseconvex([[0,0],[-1,-1],[2,0],[3,1],[2,2],[0,2]])
# should be false

from fr3d.geometry.convex_regions import counterclockwiseinside

v = counterclockwiseinside([1,1], [[0,0],[2,0],[3,1],[2,2],[0,2]])
# should be true

v = counterclockwiseinside([6,1], [[0,0],[2,0],[3,1],[2,2],[0,2]])
# should be false

