""" superpositions.py contains methods to superpose sets of 3-dimensional coordinates
"""

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

def besttransformation(X,Y)
	



return R, scale, shift, sshift
