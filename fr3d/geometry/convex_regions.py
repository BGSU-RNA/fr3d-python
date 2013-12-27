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


def ptinlefthalf(P1, P2, P3):
    """Points P1, P2, and P3.  Checking if P3 lies in the left half plane 
    defined by the vector from P1 to P2.
    """
    V1 = np.array([P2[0]-P1[0],P2[1]-P1[1]])
    V2 = np.array([P3[0]-P1[0],P3[1]-P1[1]])
    out = totheleft(V2,V1)
    return out
    
    
def testcounterclockwiseconvex(P):
    """Suppose that you have a sequence of points in the plane, 
    P = p1, p2, p3, ... , pN.  These are supposed to go counterclockwise and 
    each new point should be to the left of the line defined by the previous 
    two.  
    I don't want to assume that pN = p1, so the shape will be closed by 
    following p1, p2, p3, ... , pN, p1.  We need a test to make sure that each 
    point lies to the left of the line defined by the previous two.  
    It's important to check all the way to the line defined by pN and p1 and 
    make sure that p2 is to the left of that.
    """
    l = len(P)

    Logic = []
    for i in range(2,l):
        out1 = ptinlefthalf(P[i-2],P[i-1],P[i])
        Logic.append(out1)
    
    out1 = ptinlefthalf(P[l-1],P[0],P[1])
    Logic.append(out1)
    
    if all(Logic)==True:  
        out = True
    else:
        out = False

    return out
    
    
def counterclockwiseinside(a, b):
    """Given a sequence of points that passes the test to be counterclockwise 
    and convex, we need a true/false function that tells whether a given point
    is inside the shape or not.
    """