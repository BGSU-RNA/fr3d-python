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
    assert len(ntlist1) >= 2

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

    if len(ntlist1) == 2:
        centers1  = []
        rotation1 = []
        centers2  = []
        rotation2 = []

        for i in xrange(len(ntlist1)):
            nt1 = ntlist1[i]
            nt2 = ntlist2[i]
            c = 'base'
            if c in nt1.centers:
                centers1.append(nt1.centers[c])
                centers2.append(nt2.centers[c])
                rotation1.append(nt1.rotation_matrix)
                rotation2.append(nt2.rotation_matrix)
            else:
                raise MissingBaseException(centers)

        centers1 = np.array(centers1)
        centers2 = np.array(centers2)
        rotation1 = np.array(rotation1)
        rotation2 = np.array(rotation2)

        discrepancy = matrix_discrepancy(centers1,rotation1,centers1,rotations2)

    else:

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


def matrix_discrepancy(centers1, rotations1, centers2, rotations2,
                       angle_weight=None, center_weight=None):
    """Compute discrepancies given matrices, not components.

    :param list centers1: A list (or list numpy.array) of all centers.
    :param list rotations1: A list of all rotation matrices.
    :param list centers2: A list (or list numpy.array) of all centers.
    :param list rotations2: A list of all rotation matrices.
    :param float angle_weight: The weight to give to the angle component of
    discrepancy.
    :param float center_weight: The weight to give to the center component of
    discrepancy.
    :returns: A float, the discprenacy.
    """

    n = len(centers1)

    assert len(centers2) == n
    assert len(rotations1) == n
    assert len(rotations2) == n
    assert n >= 2

    if not angle_weight:
        angle_weight = 1.0

    if not center_weight:
        center_weight = [1.0] * n

    if n > 2:
        rotation_matrix, new1, mean1, RMSD, sse = \
            besttransformation_weighted(centers1, centers2, center_weight)

        orientation_error = 0
        angles = []
        for r1, r2 in zip(rotations1, rotations2):
            if r1.shape[0] > 0 and r2.shape[0] > 0:
                angle = angle_of_rotation(np.dot(np.dot(rotation_matrix, r2),
                                                 np.transpose(r1)))
                orientation_error += np.square(angle)
        discrepancy = np.sqrt(sse + angle_weight * orientation_error) / n

    else:

        R1 = np.dot(np.transpose(rotations1[1]),rotations1[0])  # rotation from nt 0 to nt1 of 1st motif
        R2 = np.dot(np.transpose(rotations2[0]),rotations2[1])  # rotation from nt 0 to nt1 of 2nd motif

        rot1 = np.dot(R1,R2)                                    #
        ang1 = angle_of_rotation(rot1)

        rot2 = np.dot(np.transpose(R1),np.transpose(R2))
        ang2 = angle_of_rotation(rot2)

        T1 = np.dot(centers1[1] - centers1[0],rotations1[0])
        T2 = np.dot(centers1[0] - centers1[1],rotations1[1])

        S1 = np.dot(centers2[1] - centers2[0],rotations2[0])
        S2 = np.dot(centers2[0] - centers2[1],rotations2[1])

        D1 = T1-S1
        D2 = T2-S2

        discrepancy  = np.sqrt(D1[0]**2 + D1[1]**2 + D1[2]**2 + (angle_weight*ang1)**2)
        discrepancy += np.sqrt(D2[0]**2 + D2[1]**2 + D2[2]**2 + (angle_weight*ang2)**2)

#        factor = 1/(4*np.sqrt(2))    # factor to multiply by discrepancy; faster to precompute?

        discrepancy  = discrepancy * 0.17677669529663687

    return discrepancy

def matrix_discrepancy_cutoff(centers1, rotations1, centers2, rotations2, cutoff,
                       angle_weight=None, center_weight=None):
    """Compute discrepancies given matrices, not components.

    :param list centers1: A list (or list numpy.array) of all centers.
    :param list rotations1: A list of all rotation matrices.
    :param list centers2: A list (or list numpy.array) of all centers.
    :param list rotations2: A list of all rotation matrices.
    :param float cutoff: If discrepancy > cutoff, return none
    :param float angle_weight: The weight to give to the angle component of
    discrepancy.
    :param float center_weight: The weight to give to the center component of
    discrepancy.
    :returns: A float, the discprenacy.
    """

    n = len(centers1)

    assert len(centers2) == n
    assert len(rotations1) == n
    assert len(rotations2) == n
    assert n >= 2

    if not angle_weight:
        angle_weight = [1] * n

    if not center_weight:
        center_weight = [1] * n

    if n > 2:
        # solve np.sqrt(sse + orientation_error) / n > cutoff to give sse + orientation_error > (n*cutoff)**2
        temp_cutoff = (n * cutoff)**2

        rotation_matrix, new1, mean1, RMSD, sse = \
            besttransformation_weighted(centers1, centers2, center_weight)

        total_error = sse
        if total_error > temp_cutoff:
            return None

        angles = []
        for r1, r2 in zip(rotations1, rotations2):
            if r1.shape[0] > 0 and r2.shape[0] > 0:
                angle = angle_of_rotation(np.dot(np.dot(rotation_matrix, r2),
                                                 np.transpose(r1)))
                total_error += np.square(angle)
                if total_error > temp_cutoff:
                    return None

        discrepancy = np.sqrt(total_error) / n

    else:

        R1 = np.dot(np.transpose(rotations1[1]),rotations1[0])  # rotation from nt 0 to nt1 of 1st motif
        R2 = np.dot(np.transpose(rotations2[0]),rotations2[1])  # rotation from nt 0 to nt1 of 2nd motif

        rot1 = np.dot(R1,R2)                                    #
        ang1 = angle_of_rotation(rot1)

        rot2 = np.dot(np.transpose(R1),np.transpose(R2))
        ang2 = angle_of_rotation(rot2)

        T1 = np.dot(centers1[1] - centers1[0],rotations1[0])
        T2 = np.dot(centers1[0] - centers1[1],rotations1[1])

        S1 = np.dot(centers2[1] - centers2[0],rotations2[0])
        S2 = np.dot(centers2[0] - centers2[1],rotations2[1])

        D1 = T1-S1
        D2 = T2-S2

        discrepancy  = np.sqrt(D1[0]**2 + D1[1]**2 + D1[2]**2 + (angle_weight[0]*ang1)**2)
        discrepancy += np.sqrt(D2[0]**2 + D2[1]**2 + D2[2]**2 + (angle_weight[0]*ang2)**2)

#        factor = 1/(4*np.sqrt(2))    # factor to multiply by discrepancy; faster to precompute?

        discrepancy  = discrepancy * 0.17677669529663687

    return discrepancy
