import numpy as np

from fr3d import definitions as defs
from fr3d.classifiers.generic import Classifier as BaseClassifier
from fr3d.data import Atom
from fr3d.data import Component


class Classifier(BaseClassifier):
    """A classifier for RNA protein interactions.
    """

    def __init__(self):
        first = {'sequence': defs.RNAbaseheavyatoms.keys(), 'symmetry': '1_555'}
        second = {'sequence': defs.aa_fg.keys()}
        distance = {'use': 'center', 'cutoff': 6.0}
        super(Classifier, self).__init__(first=first, second=second,
                                         distance=distance)
    def classification(self, first, second):
        if not hasattr(first, 'rotation_matrix'):
            return None

        trans_first = first.translate_rotate(first)
        trans_second = first.translate_rotate(second)
        min_xy, mean_z =  self.distance_metrics(trans_first, trans_second)
                
        if min_xy <= 5 and mean_z < 2:
            transformation_matrix = first.standard_transformation()

        if transformation_matrix == None:
            return None
        trans_first = first.transform(transformation_matrix)
        trans_second = second.transform(transformation_matrix)
        min_xy, mean_z =  self.distance_metrics(trans_first, trans_second)

        if min_xy <= 10 and mean_z < 4:

            return self.classify_stacking(trans_first, trans_second)
        elif 5 < min_xy < 20 and mean_z < 2.5:
            return self.classify_pairing(trans_first, trans_second)

    def distance_metrics(self, base_residue, aa_residue):
        squared_xy_dist_list = []
        aa_z_list = []
        base_coord = base_residue.centers["base"]
        for aa_atom in aa_residue.atoms(name=defs.aa_fg[aa_residue.sequence]):
            try:
                aa_x = np.subtract(aa_atom.x, base_coord[0])
                aa_y= np.subtract(aa_atom.y, base_coord[1])
                aa_z = np.subtract(aa_atom.z, base_coord[2])
                squared_xy_dist = (aa_x**2) + (aa_y**2)
                squared_xy_dist_list.append(squared_xy_dist)
                aa_z_list.append(aa_z)
            except:
                print "Incomplete residue"
        
        mean_z = np.mean(aa_z)
        min_xy = min(squared_xy_dist_list)
        return min_xy, mean_z

    def classify_stacking(self, base_residue, aa_residue):
        stacked_planar_aa = set (["TRP", "TYR", "PHE", "HIS", "ARG", "ASN", 
        "GLN", "GLU", "ASP"])
        stacked_aliphatic = set(["LEU", "ILE", "PRO", "THR", "MET", "CYS", "VAL", "ALA", "SER"])
    
        perpendicular_aa = set (["HIS", "ARG", "LYS", "ASN", "GLN", "LEU", "ILE"])
        perpendicular_stack_aa = set(["PHE", "TYR"])
        angle = base_residue.angle_between_normals(aa_residue)
        if aa_residue.sequence in stacked_planar_aa:
            if angle <= 0.67 or 2.45 <= angle <= 3.15:
                face = self.detect_face(base_residue, aa_residue)
                return ("stacked", face)
            elif aa_residue.sequence in stacked_aliphatic:
                return stacking_tilt(aa_residue)
                
            elif 1.2<= angle <=1.64:
                if aa_residue.sequence in perpendicular_stack_aa:
                    return ("perpendicular stacking", None)
                elif aa_residue.sequence in perpendicular_aa:
                    return ("cation-pi", None)

    def classify_pairing(self, first, second):
        pseudopair_aa = set (["ASP", "GLU", "ASN", "GLN", "HIS", "LYS", "ARG",
                              "SER", "TYR", "TRP", "PHE", "VAL", "LEU", "ILE",
                              "MET"])

        if second.sequence in pseudopair_aa:
            angle = first.angle_between_normals(second)
            if 0 <= angle <= 0.75 or 2.6 <= angle <= 3.14:
                edge = self.detect_edge(first, second)
                return ("pseudopair", edge)

    def detect_edge(self, base_residue, aa_residue):
        base_coord = base_residue.centers["base"]
        aa_coord = aa_residue.centers["aa_fg"]

        base_center_x = base_coord[0]
        base_center_y = base_coord[1]

        aa_center_x = aa_coord[0]
        aa_center_y = aa_coord[1]

        y = aa_center_y - base_center_y
        x = aa_center_x - base_center_x
        angle_aa = np.arctan2(y,x)

        if -1 <= angle_aa <= 0:
            return "fgbS"
        elif angle_aa <=1:
            return "fgbWC"
        elif 1.4 <= angle_aa <= 3.2:
            return "fgbH"
