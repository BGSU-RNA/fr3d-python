import numpy as np


def totheleft(a, b):
    """Inputs are two 2-dimensional vectors a and b, thought of as being in the
    plane with their tails at the origin.  Output is true if a points into the 
    half plane that is to the left side of where b points.
    """
    
    a=np.append(a,0)
    b=np.append(b,0)
    cross = np.cross(b,a)

    if cross[2]>=0:
        out = True
    else:
        out = False
    
    return out
    
def counterclockwiseinside(a, b):
    """
    """

    
def testcounterclockwiseconvex(a, b):
    """
    """