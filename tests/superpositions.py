""" superpositions.py contains methods to superpose sets of 3-dimensional coordinates
"""

# bestrotation(a,b) finds the 3x3 rotation matrix which optimally superimposes
# the nx3 matrix of points a onto the nx3 matrix of points b.
# It is assumed that the coordinates of a and b have mean 0.
# One reference is this: http://en.wikipedia.org/wiki/Kabsch_algorithm
# Another is a python implementation that goes with pymol, see
# http://www.pymolwiki.org/index.php/Kabsch

###############################################################################
#The following is a version that appears to work!
"""
Created on Tue Aug 13 14:33:48 2013
"""
def besttransformation(set1, set2):
 
    #Neat printing.
    import pprint
    #Python math skills.
    import numpy
    
    ##Import a list of [[x,y,z], .....] or as an array.
    #Examples:
    #set1=[[1.0, 2.0, 3.0], [4,5,2], [9,2,4]]
    #sel2=[[4.0, 4.1, 5.0], [2,3,1], [3,1,4]]    
    
    
    # Check to make sure same number of (x,y,z) coordinates in each set. 
    #If condition is not true this program stops.
    
    assert len(set1) == len(set2)
    L = len(set1)
    assert L > 0
    
    #These add all x_{ij}'s, y_{ij}'s, z_{ij}'s for each element in sel_j j=1,2 i=1,2,3,..L
    # /float(L) divides by L
    #This creates a [mean x_j, mean y_j, mean z_j] for each set of coordinates j.
    COM1 = numpy.sum(set1,axis=0) / float(L)
    COM2 = numpy.sum(set2,axis=0) / float(L)
        
    #individual deviation from group means.
    #Specifically,
    #subtract x_{ij} by the mean x_j's. Here i=1,2,..L.  Same for y, and z.

    dev1=set1-COM1
    dev2=set2-COM2
    

    #[sum(x_{i1}^2) + sum(x_{i2}^2) + sum(y_{i1}^2) + sum(y_{i2}^2) + sum(z_{i1}^2)+sum(z_{i2}^2)]
    E0 = numpy.sum( numpy.sum(dev1 * dev1,axis=0),axis=0) + numpy.sum( numpy.sum(dev2*dev2,axis=0),axis=0)

#********* I DO NOT UNDERSTAND STARTING HERE TO......************
    # V and Wt are the orthonormal bases that when multiplied by each other give us
    # the rotation matrix, U.  S, (Sigma, from SVD) provides us with the error.
    V, S, Wt = numpy.linalg.svd( numpy.dot( numpy.transpose(dev2), dev1))
    
    # we already have our solution, in the results from SVD.
    # we just need to check for reflections and then produce
    # the rotation.  V and Wt are orthonormal, so their det's
    # are +/-1.
    
    reflect = float(str(float(numpy.linalg.det(V) * numpy.linalg.det(Wt))))
     
    if reflect == -1.0:
    	S[-1] = -S[-1]
    	V[:,-1] = -V[:,-1]
     
    RMSD = E0 - (2.0 * sum(S))
    RMSD = numpy.sqrt(abs(RMSD / L))

#*********....HERE.************


    #The transformation matrix, U, is now V*Wt
    U = numpy.dot(V, Wt)
     
    # rotate and translate the molecule
    #sel2 = numpy.dot((mol2 - COM2), U)
    new1=numpy.dot(dev1, U)
    new2=numpy.dot(dev2, U)    
    
    #Return the transformation matrix, the new coordinates for the two 
    #set of coordinates, respectively.
    return U, new1, new2
    
    #print 'The transformation matrix'
    #print U
    
    #print 'New Coordinates for Set 1'
    #pp1 = pprint.PrettyPrinter(indent=4)
    #pp1.pprint(set1)

    #print 'New Coordinates for Set 2'
    #pp1 = pprint.PrettyPrinter(indent=4)
    #pp1.pprint(set2)
############################################################################





























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

#def bestrotation(a,b):


 #   return rotation




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

#def besttransformation(a,b):
    



  #  return rotation, shift
