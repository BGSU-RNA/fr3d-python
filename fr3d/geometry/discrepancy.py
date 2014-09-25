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

class LengthofBaseWeightError(Exception):
    """An exception that is raised when the list of base weights is not equal
    in length to the list of nucleotides.
    """
    pass

class LengthofPWeightError(Exception):
    """An exception that is raised when the list of phosphate weights is not
    equal in length to the list of nucleotides.
    """
    pass

class LengthofC1starWeightError(Exception):
    """An exception that is raised when the list of C1* weights is not equal in
    length to the list of nucleotides.
    """
    pass

def discrepancy(ntlist1, ntlist2, centers=['base'], base_weights=1.0,
                P_weights=1.0, C1star_weights = 1.0, angleweight=1.0):
    """Compute the geometric discrepancy between two lists of components.

    :ntlist1: The first list of components.
    :ntlist2: The second list of components.
    :centers: A list of center names to use, such as
        ['base', 'P', 'C1*', 'ribose']
    :base_weights: The base weights to use. If only one weight is given
    it is used for all centers.  Otherwise, provide list of base weights with
    same length as the length of ntlist1
    :P_weigths: The phosphate weights to use. If only one weight is given
    it is used for all.  Otherwise, provide list of phosphate weights
    with same length as the length of ntlist1
    :C1star_weights: The C1* weights to use. If only one weight is given
    it is used for all.  Otherwise, provide list of C1* weights with
    same length as the length of ntlist1
    :angleweight: The weighting for angles to use.
    :returns: The geometric discrepancy.
    """

    assert len(ntlist1) == len(ntlist2)

    # TODO: Should we allow users to pass a tuple too?
    if not isinstance(centers, list):
        centers = [centers]

    if not isinstance(base_weights, list):
        base_weights = [base_weights] * len(ntlist1)

    if len(base_weights)!= len(ntlist1):
        raise LengthofBaseWeightError('Weight length does not match # of nucl.')

    if not isinstance(P_weights, list):
        P_weights = [P_weights] * len(ntlist1)

    if len(P_weights)!= len(ntlist1):
        raise LengthofPWeightError('Weight length does not match # of nucl.')

    if not isinstance(C1star_weights, list):
        C1star_weights = [C1star_weights] * len(ntlist1)

    if len(C1star_weights)!= len(ntlist1):
        raise LengthofC1starWeightError('Weight length does not match # of nucl.')

    R = []
    S = []
    W = []

    for i in xrange(len(ntlist1)):
        nt1 = ntlist1[i]
        nt2 = ntlist2[i]
        for c in centers:
            if c=='base':
                if c in nt1.centers:
                    R.append(nt1.centers[c])
                    S.append(nt2.centers[c])
                    W.append(base_weights[i])
                else:
                    if c=='base':
                        raise MissingBaseException(centers)
            if c=='P':
                if nt1.coordinates(type = 'P')!=[]:
                    R.append(nt1.coordinates(type = 'P')[0])
                    S.append(nt2.coordinates(type = 'P')[0])
                    l=len(nt1.coordinates(type = 'P'))
                    for z in range(0,l):
                        W.append(P_weights[i])
                else:
                    raise MissingPhosphateException(centers)
            if c=='C1*':
                if nt1.coordinates(type = 'C1*')!=[] and nt2.coordinates(type = 'C1*')!=[]:
                    R.append(nt1.coordinates(type = 'C1*')[0])
                    S.append(nt2.coordinates(type = 'C1*')[0])
                    l=len(nt1.coordinates(type = 'C1*'))
                    for q in range(0,l):
                        W.append(C1star_weights[i])
                else:
                    raise MissingPhosphateException(centers)
    #rotation_matrix, _, _, RMSD = besttransformation(R, S)
    # Superimpose R and S with weights? I need to make changes.
    rotation_matrix, new1, mean1, RMSD, sse = besttransformation_weighted(R, S, W)
    rotation_matrix = np.transpose(rotation_matrix)
    #The rotation_matrix that is outputted from besttransformation is the
    #transpose of the one you want.

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
