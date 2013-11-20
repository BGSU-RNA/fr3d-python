# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 12:50:43 2013

@author: zirbel
"""


from fr3d.tests import define_1S72_nucleotides
from fr3d.geometry import discrepancy

# Test:  superimpose 7 nucleotides:

d = discrepancy([nt77_9,nt78_9,nt79_9,nt80_9,nt102_9,nt103_9,nt104_9],[nt212_0,nt213_0,nt214_0,nt215_0,nt225_0,nt226_0,nt227_0])
print d

# discrepancy should be about 0.1159
# discrepancy is sqrt(L^2 + A1^2 + A2^2 + ... + A7^2)/7, where L is location error and A1 ... A7 are angles between bases in radians
# angles between bases should be about 0.0700, 0.0287, 0.0547, 0.0754, 0.0262, 0.0580, 0.1083
# or else twice that, not quite sure

# rotation matrix should be about:
#    0.4780    0.4348   -0.7632
#   -0.2844   -0.7455   -0.6028
#   -0.8311    0.5051   -0.2327

# Test:  superimpose 5 nucleotides:

d = discrepancy([nt77_9,nt78_9,nt79_9,nt80_9,nt102_9],[nt212_0,nt213_0,nt214_0,nt215_0,nt225_0])
print d

# discrepancy should be about 0.0748
# discrepancy is sqrt(L^2 + A1^2 + A2^2 + ... + A7^2)/7, where L is location error and A1 ... A7 are angles between bases in radians
# angles between bases should be about     0.0756    0.0266    0.0652    0.0797    0.0282
# or else twice that, not quite sure

# rotation matrix should be about:
#    0.4593    0.4257   -0.7796
#   -0.2964   -0.7539   -0.5863
#   -0.8374    0.5004   -0.2201

