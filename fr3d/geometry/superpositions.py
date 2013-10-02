"""Superpositions contains functions to superpose sets of 3-dimensional
coordinates.
"""

import numpy
from RMSD import RMSD

def besttransformation(set1, set2):
    """This finds the 3x3 rotation matrix which optimally superimposes
    the nx3 matrix of points in set1 onto the nx3 matrix of points set2.
    One reference is this: http://en.wikipedia.org/wiki/Kabsch_algorithm
    Another is a python implementation that goes with pymol, see
    http://www.pymolwiki.org/index.php/Kabsch

    The arguments should be lists or numpy arrays of the form exampled below:
    set1=[[1.0, 2.0, 3.0], [4,5,2], [9,2,4], ....,[3,1,6]]
    sel2=[[4.0, 4.1, 5.0], [2,3,1], [3,1,4], ....,[1, 7, 7]]

    :set1: A list or a numpy array of (n, 3) coordinates.
    :set2: A list or a numpy array of (n, 3) coordinates.
    :returns: The transformation matrix, the new coordinates for the two
    set of coordinates, respectively.
    """

    # Check to make sure same number of (x,y,z) coordinates both sets.
    #If condition is not true this program stops.
    assert len(set1) == len(set2)
    length = len(set1)
    assert length > 0

    # Translation Step is beginning.
    # These add all x_{ij}'s, y_{ij}'s, z_{ij}'s for each element in
    # setj j=1,2
    # i=1,2,3,..length
    # /float(length) divides by length
    # This creates a [mean x_j, mean y_j, mean z_j] for both sets of
    # coordinates setj j=1, or 2, or the centroid of both sets.
    mean1 = numpy.sum(set1, axis=0) / float(length)
    mean2 = numpy.sum(set2, axis=0) / float(length)

    # Next,
    # Subtract x_{ij} by the mean of x_j's. Here i=1,2,..length.
    # Same for y, and z.
    dev1 = set1 - mean1
    dev2 = set2 - mean2
    #  Thus, both sets are translated, so that their centroid coincides with
    # the origin of the coordinate system.
    # Translation Step is now completed.

    # Begin Step to Compute the 3X3 Covariance Matrix, A.
    A = numpy.dot(numpy.transpose(dev2), dev1)
    # Covariance Matrix, A, is now calculated

    # Begin of the Computation of the optimal rotation matrix using
    # Singular Value Decomposition (SVD)
    V, diagS, Wt = numpy.linalg.svd(A)
    # V and Wt are 3x3 orthonormal bases, diagS is the diagonal elements of
    # a 3x3 diagonal matrix, S, in regular SVD.  In SVD, recall that the
    # Covariance matrix, A, is A=V*S*transpose(W) (matrix multiplication).
    #S=numpy.diag(diagS)
    #A=numpy.dot(V,numpy.dot(S,Wt)

    # The next step is to decide whether we need to correct our rotation
    # matrix to ensure a right-handed coordinate system
    # we just need to check for reflections and then produce
    # the rotation.  V and Wt are orthonormal, so their det's
    # are +/-1.
    I = numpy.matrix(numpy.identity(3))
    d = numpy.linalg.det(numpy.dot(numpy.transpose(Wt), numpy.transpose(V)))
    if numpy.isclose(d, -1.0):
        I[2, 2] = d

    # End of the Computation of the optimal rotation matrix

    #The transformation matrix, U, is now V*Wt
    U = numpy.dot(numpy.dot(numpy.transpose(Wt), I), numpy.transpose(V))

    # rotate and translate the molecule
    #sel2 = numpy.dot((set2 - Mean2), U)
    new1 = numpy.dot(dev1, U)
    new2 = dev2
    rmsd = RMSD(new1,new2)
    #Return the transformation matrix, the new coordinates for the two
    #set of coordinates, respectively.
    return U, new1, mean1, rmsd

# Below is an older Matlab implementation, but it would be better to use the
# SVD implementation described above

# % zBestRotation(X,Y) finds the least squares rotation
# % of points X onto points Y
# %
# % X and Y are n by 3 matrices containing the locations of corresponding
# % points
# % What is returned is the best fit to Y = X*R'
# % R is the 3x3 rotation matrix

# % It is assumed that the means have been subtracted from X and Y

# function [R] = zBestRotation(A,B)

# n = length(A(:,1));                 % number of data points

# M = zeros(3,3);                     % initialize the matrix M

# for i=1:n,                          % go through all data points
#   M = M + A(i,:)' * B(i,:);         % "outer product" of two vectors ->matrix
# end

# [v,d] = eig(M'*M);                  % eigenvector matrix, eigenvalue matrix
# [e,i] = sort(-diag(d));             % decreasing sort the of eigenvalues

# a = v(:,i(1));                      % eigenvector with largest eigenvalue
# b = v(:,i(2));
# c = v(:,i(3));

# Ma = M*a/sqrt(-e(1));
# Mb = M*b/sqrt(-e(2));

# if det([a b c]) > 0,
#   g = cross(Ma,Mb);
# else
#   g = -cross(Ma,Mb);
# end

# R = [a b c]*[Ma Mb g]';
