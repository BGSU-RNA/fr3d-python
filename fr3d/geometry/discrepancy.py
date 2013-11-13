# -*- coding: utf-8 -*-
"""
Created on Wed Nov 06 17:19:16 2013

@author: zirbel
"""

import fr3d.geometry.angleofrotation

# centers is a list like ['base','P','C1*','ribose']

def discrepancy(ntlist1,ntlist2,centers=['base'],weights=1.0,angleweight=1.0):
    
    assert(length(ntlist1)) == length(ntlist2))

    if not isinstance(weights, list):
        weights = [weights] * len(centers)

    R = []
    S = []
    
    for i in range(0,length(ntlist1)):
        nt1 = ntlist1[i]
        nt2 = ntlist2[i]
        for c in centers:
            if c in nt1.centers:
                R.append(nt1.centers[c])
                S.append(nt2.centers[c])
                W.append(weights[i])
            elif c in nt1.atoms:
                R.append(nt1._atoms[c])
                S.append(nt2._atoms[c])
                W.append(weights[i])
            else:
                # raise an exception (or just notify the user?)
                
    # superimpose R and S
                
    rotation_matrix, RMSD = superimposeweighted(R,S,W)
    
    # loop through the bases and calculate angles between them

    orientationerror = 0
    
    if 'base' in centers:
        for i in range(0,length(ntlist1)):
            R1 = ntlist1[i].rotation_matrix
            R2 = ntlist2[i].rotation_matrix
            
            # calculate angle in radians, or adjust with angleweight factor
            angle = angleofrotation(R1 * rotation_matrix * R2')

            orientationerror += angle^2
            
    discrepancy = sqrt(RMSD^2 + angleweight*orientationerror) / length(ntlist1)
    