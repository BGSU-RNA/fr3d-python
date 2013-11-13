# -*- coding: utf-8 -*-

import numpy as np

import fr3d.geometry.angleofrotation


class MissingBaseException(Exception):
    """An exception that is raised when a base center that does not exist is
    requested.
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

    assert length(ntlist1) == length(ntlist2)

    # TODO: Should we allow users to pass a tuple too?
    if not isinstance(centers, list):
        centers = [centers]

    if not isinstance(weights, list):
        weights = [weights] * len(centers)

    R = []
    S = []

    for i in xrange(len(ntlist1)):
        nt1 = ntlist1[i]
        nt2 = ntlist2[i]
        for c in centers:

            if c in nt1.centers:
                R.append(nt1.centers[c])
                S.append(nt2.centers[c])
                W.append(weights[i])
            else:
                raise MissingBaseException(center)

    # superimpose R and S
    rotation_matrix, RMSD = superimposeweighted(R, S, W)

    # loop through the bases and calculate angles between them

    orientationerror = 0

    if 'base' in centers:
        for i in xrange(len(ntlist1)):
            R1 = ntlist1[i].rotation_matrix
            R2 = ntlist2[i].rotation_matrix

            # calculate angle in radians, or adjust with angleweight factor
            angle = angleofrotation(R1 * rotation_matrix * R2')

            orientationerror += angle^2

    discrepancy = sqrt(RMSD^2 + angleweight*orientationerror) / length(ntlist1)
