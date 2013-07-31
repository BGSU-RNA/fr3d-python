
from numpy import array
from fr3d.geometry.superpositions import bestrotation

a = array([(1,0,1),(0,1,0),(0,0,1)])

rotation = bestrotation(a,a)



