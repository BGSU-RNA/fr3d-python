from fr3d import definitions as defs
from fr3d.classifiers.generic import Classifier as BaseClassifier
import numpy as np


class Classifier(BaseClassifier):
    """A classifier for RNA protein interactions.
    """

    def __init__(self):
        first = {'sequence': defs.RNAbaseheavyatoms.keys()}
        second = {'sequence': defs.aa_backconnect.keys()}
        distance = {'use': 'center', 'cutoff': 4.0}
        super(Classifier, self).__init__(first=first, second=second,
                                         distance=distance)

    def classification(self, first, second):
        squared_xy_dist_list = []
        aa_z =[]
        for aa_atom in second.atoms(name=defs.aa_fg[second.sequence]):
            key = aa_atom.name
            aa_x= second[key][0]
            aa_y= second[key][1]
        
        squared_xy_dist = (aa_x**2) + (aa_y**2)
        squared_xy_dist_list.append(squared_xy_dist)
        
        aa_z.append(second[key][2])
        
        mean_z = np.mean(aa_z)
    
        if min(squared_xy_dist_list) <= 3:
            if second.sequence in set (["TRP", "TYR", "PHE", "HIS", "ARG", "LYS", "ASN", "GLN"]):
                return stacking_angle(first, second)
            else:
                return stacking_tilt(second)
            
        elif 3.1 < min(squared_xy_dist_list)< 35.2 and -2.0 <= mean_z < 2.0:
            if second.sequence in set (["ASP", "GLU", "ASN", "GLN", "HIS", "ARG", "LYS", "SER", "TYR", "TRP", "PHE", "VAL", "LEU", "ILE"]):
                angle= calculate_angle(first, second)
    
            if -1.3 <= angle <= 0.83 or 2.25 <= angle <= 3.14:
                return "pseudopair"
    
    def calculate_angle (base_residue, aa_residue):
        vec1 = vector_calculation(base_residue)
        vec2 = vector_calculation(aa_residue)
                
        angle = angle_between_planes(vec1, vec2)
        return angle

def stacking_angle (base_residue, aa_residue):
    vec1 = vector_calculation(base_residue)
    vec2 = vector_calculation(aa_residue)
                
    angle = angle_between_planes(vec1, vec2)
    
    if aa_residue.sequence in set (["TRP", "TYR", "PHE", "HIS", "ARG"]):
        if angle <=0.79 or 2.35 <= angle <= 3.15:
            return "stacked"
        elif aa_residue.sequence in set (["TYR", "HIS", "ARG", "LYS", "ASN", "GLN"]):
            if 1.32<= angle <=1.64:
                return "perpendicular"

def stacking_tilt(second):
    baa_dist_list = []     
        
    for aa_atom in second.atoms(name=defs.aa_fg[second.sequence]):
        key = aa_atom.name
        aa_z = second[key][2]
        baa_dist_list.append(aa_z)        
    max_baa = max(baa_dist_list)
    min_baa = min(baa_dist_list)
    #print 'max distance: %s' % max_baa + ' min distance: %s' % min_baa
    diff = max_baa - min_baa
    #print aa_residue.unit_id(), diff
    return diff <= defs.tilt_cutoff[second.sequence]
    
def vector_calculation(residue):
    key = residue.sequence
    P1 = residue.centers[defs.Normal_residue[key][0]]
    P2 = residue.centers[defs.Normal_residue[key][1]]
    P3 = residue.centers[defs.Normal_residue[key][2]]
    #print P1, P2, P3
    vector = np.cross((P2 - P1),(P3-P1))
    return vector

def angle_between_planes (vec1, vec2):
    cosang = np.dot(vec1, vec2)
    sinang = np.linalg.norm(np.cross(vec1, vec2))
    angle = np.arctan2(sinang, cosang)
    return angle
