from fr3d import definitions as defs
from fr3d.classifiers.generic import Classifier as BaseClassifier
from fr3d.geometry import angleofrotation as angrot
from fr3d.data import components as compnt
import numpy as np


class Classifier(BaseClassifier):
    """A classifier for RNA protein interactions.
    """

    def __init__(self):
        first = {'sequence': defs.RNAbaseheavyatoms.keys()}
        second = {'sequence': defs.aa_backconnect.keys()}
        distance = {'use': 'center', 'cutoff': 10.0}
        super(Classifier, self).__init__(first=first, second=second,
                                         distance=distance)

    def classification(self, first, second):
        squared_xy_dist_list = []
        aa_z =[]
        """Defines different sets of amino acids"""
        stacked_aa = set (["TRP", "TYR", "PHE", "HIS", "ARG", "LYS", "ASN", "GLN", "LEU", "ILE", "PRO", "THR"])
        pseudopair_aa = set (["ASP", "GLU", "ASN", "GLN", "HIS", "LYS", "ARG", "SER", "TYR", "TRP", "PHE", "VAL", "LEU", "ILE", "MET"])

        for aa_atom in second.atoms(name=defs.aa_fg[second.sequence]):
            key = aa_atom.name
            aa_x= second[key][0]
            aa_y= second[key][1]
        
        squared_xy_dist = (aa_x**2) + (aa_y**2)
        squared_xy_dist_list.append(squared_xy_dist)
        
        aa_z.append(second[key][2])
        
        mean_z = np.mean(aa_z)
    
        if min(squared_xy_dist_list) <= 3:
            if second.sequence in stacked_aa:
                return stacking_angle(first, second)
            else:
                return stacking_tilt(second)
            
        elif 3.1 < min(squared_xy_dist_list)< 35.2 and -2.0 <= mean_z < 2.0:
            if second.sequence in pseudopair_aa:
                angle= compnt.angle_between_normals(first, second)
    
            if -1.3 <= angle <= 0.75 or 2.6 <= angle <= 3.14:
                if second.enough_HBs(first, second):
                    return "pseudopair"
            else:
                return None

    def detect_edge(self, base_residue, aa_residue):
        aa_x = 0
        aa_y = 0
        n = 0
        base_x = 0
        base_y = 0
        aa_coordinates = aa_residue.compnt.transform()
        base_coordinates = base_residue.compnt.transform()
        
        for aa_atom in aa_residue.atoms(name=defs.aa_fg[aa_residue.sequence]):
            key = aa_atom.name
            aa_x+= aa_coordinates[key][0]
            aa_y+= aa_coordinates[key][1]
            n +=1
        aa_center_x = aa_x/n        
        aa_center_y = aa_y/n        
          
        for base_atom in base_residue.atoms(name=defs.RNAbaseheavyatoms[base_residue.sequence]):
            key = base_atom.name
            base_x+= base_coordinates[key][0]
            base_y+= base_coordinates[key][1]
            n +=1
        base_center_x = aa_x/n        
        base_center_y = aa_y/n  
    
        y = aa_center_y - base_center_y
        x = aa_center_x - base_center_x
        angle_aa = np.arctan2(y,x)
    
        if -1 <= angle_aa <= 0:
            return "Sugar"
        elif angle_aa <=1:
            return "WC"
        elif 1.4 <= angle_aa <= 3.2:
            return "Hoogsteen"
    
        
    def enough_HBs(self, base_residue, aa_residue):
        """Calculates atom to atom distance of part "aa_part" of neighboring amino acids
    of type "aa" from each atom of base. Only returns a pair of aa/nt if two
    or more atoms are within the cutoff distance"""
        min_distance = 4
        HB_atoms = set(['N', 'NH1','NH2','NE','NZ','ND1','NE2','O','OD1','OE1','OE2', 'OG', 'OH'])
        n = 0
        base_seq = base_residue.sequence()
        for base_atom in base_residue.atoms(name=defs.RNAbaseheavyatoms[base_seq]):
            for aa_atom in aa_residue.atoms(name=defs.aa_fg[aa_residue.sequence]):

                distance = np.subtract(base_atom.coordinates(), aa_atom.coordinates())
                distance = np.linalg.norm(distance)
            if distance <= min_distance and aa_atom.name in HB_atoms:
                n = n+1
        return n>=2
            

def stacking_angle (base_residue, aa_residue):
    vec1 = compnt.normal_calculation(base_residue)
    vec2 = compnt.normal_calculation(aa_residue)

    stacked_aa = set (["TRP", "TYR", "PHE", "HIS", "ARG", "LYS", "LEU", "ILE", "PRO", "ASN", "GLN"])       
    perpendicular_aa = set (["TYR", "HIS", "ARG", "LYS", "ASN", "GLN", "LEU", "ILE"])
                
    angle = angrot.angle_between_planes(vec1, vec2)
    
    if aa_residue.sequence in stacked_aa:
        if angle <=0.67 or 2.43 <= angle <= 3.15:
            return "stacked"
        elif aa_residue.sequence in perpendicular_aa:
            if 1.2<= angle <=1.64:
                return "perpendicular"

def stacking_tilt(second):
    baa_dist_list = []     
        
    for aa_atom in second.atoms(name=defs.aa_fg[second.sequence]):
        key = aa_atom.name
        aa_z = second[key][2]
        baa_dist_list.append(aa_z)        
    max_baa = max(baa_dist_list)
    min_baa = min(baa_dist_list)
    
    diff = max_baa - min_baa
    return diff <= defs.tilt_cutoff[second.sequence]
    
