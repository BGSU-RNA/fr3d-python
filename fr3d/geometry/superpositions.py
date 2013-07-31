""" superpositions.py contains methods to superpose sets of 3-dimensional coordinates
"""

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

# scale = sqrt(sum(sum(B.*B)) / sum(sum(A.*A)));  % the optimal scale

# shift = (mY - mX*R')';           % the optimal shift with no scaling

# sshift = (mY - scale*mX*R')';    % the optimal shift when rescaling

# % Note: there is no need to calculate and return the optimal scale in our
# % applications, since we don't want to rescale.  There should be two versions
# % of this program, one with scaling, one without.

def besttransformation(a,b):
    



    return rotation, scale, shift, scaledshift
