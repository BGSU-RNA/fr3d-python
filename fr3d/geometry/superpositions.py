""" superpositions.py contains methods to superpose sets of 3-dimensional coordinates
"""

# bestrotation(a,b) finds the 3x3 rotation matrix which optimally superimposes
# the nx3 matrix of points a onto the nx3 matrix of points b.
# It is assumed that the coordinates of a and b have mean 0.
# One reference is this: http://en.wikipedia.org/wiki/Kabsch_algorithm
# Another is a python implementation that goes with pymol, see
# http://www.pymolwiki.org/index.php/Kabsch

# Below is an older Matlab implementation, but it would be better to use the
# SVD implementation described above

# % zBestRotation(X,Y) finds the least squares rotation
# % of points X onto points Y
# %
# % X and Y are n by 3 matrices containing the locations of corresponding points
# % What is returned is the best fit to Y = X*R'
# % R is the 3x3 rotation matrix

# % It is assumed that the means have been subtracted from X and Y

# function [R] = zBestRotation(A,B)

# n = length(A(:,1));                 % number of data points

# M = zeros(3,3);                     % initialize the matrix M

# for i=1:n,                          % go through all data points
#   M = M + A(i,:)' * B(i,:);         % "outer product" of two vectors -> matrix
# end

# [v,d] = eig(M'*M);                  % eigenvector matrix, eigenvalue matrix
# [e,i] = sort(-diag(d));             % sort the eigenvalues in decreasing order

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

def bestrotation(a,b):


    return rotation




# % zBestTransformation(X,Y) finds the least squares translation and rotation
# % of points X onto points Y
# %
# % X and Y are n by 3 matrices containing the locations of corresponding points
# % What is returned is the best fit to Y = shift + scale*X*R'
# % R is the 3x3 rotation matrix
# % scale is the scalar scale factor
# % shift is the 3x1 translation between centers of mass

# function [R,scale,shift,sshift] = zBestTransformation(X,Y)

# n = length(X(:,1));                 % number of data points

# mX = mean(X);
# mY = mean(Y);

# A = X-ones(n,1)*mX;                 % subtract the mean vector
# B = Y-ones(n,1)*mY;                 % subtract the mean vector

# R = zBestRotation(A,B);             % find optimal rotation matrix

# shift = (mY - mX*R')';           % the optimal shift with no scaling

def besttransformation(a,b):
    



    return rotation, shift
