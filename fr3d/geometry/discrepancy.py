# -*- coding: utf-8 -*-

import numpy as np
from fr3d.geometry.angleofrotation import angle_of_rotation
from fr3d.geometry.superpositions import besttransformation
from fr3d.geometry.superpositions import besttransformation_weighted


class MissingBaseException(Exception):
    """An exception that is raised when a base center that does not exist is
    requested.
    """
    pass

class MissingPhosphateException(Exception):
    """An exception that is raised when the phospate atoms are not present in
    inputted model.
    """
    pass

def discrepancy(ntlist1, ntlist2, centers=['base'], weights=1.0,
                angleweight=1.0):
    """Compute the geometric discrepancy between two lists of components.

    :ntlist1: The first list of components.
    :ntlist2: The second list of components.
    :centers: A list of center names to use, such as
        ['base', 'P', 'C1*', 'ribose']
    :weights: The weights to use. If only one weight is given it is used for
    all centers.
    :angleweight: The weighting for angles to use.
    :returns: The geometric discrepancy.
    """

    assert len(ntlist1) == len(ntlist2)

    # TODO: Should we allow users to pass a tuple too?
    if not isinstance(centers, list):
        centers = [centers]

    if not isinstance(weights, list):
        weights = [weights] * len(centers)
    
    if len(weights) != len(ntlist1):
        weights = np.ones(len(ntlist1))
        
    R = []
    S = []
    W = []
    
    for i in xrange(len(ntlist1)):
        nt1 = ntlist1[i]
        nt2 = ntlist2[i]
        for c in centers:

            if c in nt1.centers:
                R.append(nt1.centers[c])
                S.append(nt2.centers[c])
                W.append(weights[i])
                #The above weights were out of range.  Need to make default weights
                #To have same dimension as the list of nucleotides.
            else:
                if c=='base':
                    raise MissingBaseException(centers)
                if c=='P':
                    if nt1.coordinates(type = 'P')!=[] and nt2.coordinates(type = 'P')!=[]:
                        R.append(nt1.coordinates(type = 'P'))
                        S.append(nt2.coordinates(type = 'P'))
                    else:
                        raise MissingPhosphateException(centers)
    #rotation_matrix, _, _, RMSD = besttransformation(R, S)
    # Superimpose R and S with weights? I need to make changes.
    rotation_matrix, new1, mean1, RMSD, sse = besttransformation_weighted(R, S, W)
    rotation_matrix = np.transpose(rotation_matrix)
    #The rotation_matrix that is outputted from besttransformation is the 
    #transpose of the one you want.

    #print "Sum of squared errors", sse
    
    #print "Rotation matrix"
    #print rotation_matrix

    # loop through the bases and calculate angles between them
    orientationerror = 0
    if 'base' in centers:
        for i in xrange(len(ntlist1)):
            R1 = ntlist1[i].rotation_matrix
            R2 = ntlist2[i].rotation_matrix
            # calculate angle in radians
            angle = angle_of_rotation(np.dot(np.dot(rotation_matrix,R1), np.transpose(R2)))
            orientationerror += np.square(angle)
    discrepancy = np.sqrt(sse + angleweight*orientationerror) / len(ntlist1)
    return discrepancy
    #I must be calculating this part incorrectly, since rotation_matrix is cor