import numpy as np

from fr3d import definitions as defs
from fr3d.classifiers.generic import Classifier as BaseClassifier
from fr3d.geometry import angleofrotation as angrot

stacked_aa = set (["TRP", "TYR", "PHE", "HIS", "ARG", "LYS", "ASN",
                   "GLN", "LEU", "ILE", "PRO", "THR"])

pseudopair_aa = set (["ASP", "GLU", "ASN", "GLN", "HIS", "LYS", "ARG",
                      "SER", "TYR", "TRP", "PHE", "VAL", "LEU", "ILE",
                      "MET"])


class Classifier(BaseClassifier):
    """A classifier for RNA protein interactions.
    """

    def __init__(self):
        first = {'sequence': defs.RNAbaseheavyatoms.keys(), 'symmetry': '1_555'}
        second = {'sequence': defs.aa_fg.keys()}
        distance = {'use': 'center', 'cutoff': 10.0}
        super(Classifier, self).__init__(first=first, second=second,
                                         distance=distance)

    def stacking_tilt(self, second):
        baa_dist_list = []
        for aa_atom in second.atoms(name=defs.aa_fg[second.sequence]):
            key = aa_atom.name
            aa_z = second[key][2]
            baa_dist_list.append(aa_z)
        max_baa = max(baa_dist_list)
        min_baa = min(baa_dist_list)

        diff = max_baa - min_baa
        if diff <= defs.tilt_cutoff[second.sequence]:
            return ("stacked", "N/A")

    def stacking_angle(self, base_residue, aa_residue):
        vec1 = base_residue.normal_calculation()
        vec2 = aa_residue.normal_calculation()

        stacked_aa = set(["TRP", "TYR", "PHE", "HIS", "ARG", "LYS", "LEU", "ILE", "PRO", "ASN", "GLN"])
        perpendicular_aa = set(["TYR", "HIS", "ARG", "LYS", "ASN", "GLN", "LEU", "ILE"])

        angle = angrot.angle_between_planes(vec1, vec2)

        if aa_residue.sequence in stacked_aa:
            if angle <= 0.67 or 2.43 <= angle <= 3.15:
                return ("stacked", "N/A")
            elif aa_residue.sequence in perpendicular_aa:
                if 1.2<= angle <=1.64:
                    return ("perpendicular", "N/A")

    def classify_stacking(first, second):
        pass

    def classification(self, first, second):
        if not hasattr(first, 'rotation_matrix'):
            return None

        transformation_matrix = first.base_transformation_matrix(second)
        if transformation_matrix == None:
            return None

        trans_first = first.transform(transformation_matrix)
        trans_second = second.transform(transformation_matrix)

        min_xy, mean_z =  self.distance_metrics(trans_first, trans_second)
        if min_xy <= 3:
            return self.classify_stacking(first, seecond)
        elif 3 < min_xy < 35.2 and -2.0 <= mean_z < 2.0:
            return self.classify_pairing(first, second)
        return None

    def detect_edge(self, base_residue, aa_residue):
        aa_x = 0
        aa_y = 0
        n = 0
        base_x = 0
        base_y = 0

        #rotation_matrix = base_residue.rotation_matrix
        #base_center = base_residue.centers(defs.RNAbaseheavyatoms[base_residue.sequence])

        for aa_atom in aa_residue.atoms(name=defs.aa_fg[aa_residue.sequence]):

            aa_coordinates = aa_residue.transform(base_residue.base_transformation_matrix())
            key = aa_atom.name
            aa_x+= aa_coordinates[key][0]
            aa_y+= aa_coordinates[key][1]
            n +=1
        aa_center_x = aa_x/n
        aa_center_y = aa_y/n

        for base_atom in base_residue.atoms(name=defs.RNAbaseheavyatoms[base_residue.sequence]):
            base_coordinates = base_residue.transform(base_residue.base_transformation_matrix())
            key = base_atom.name
            base_x+= base_coordinates[key][0]
            base_y+= base_coordinates[key][1]
            n += 1
        base_center_x = aa_x/n
        base_center_y = aa_y/n

        y = aa_center_y - base_center_y
        x = aa_center_x - base_center_x
        angle_aa = np.arctan2(y,x)

        if -1 <= angle_aa <= 0:
            return "fgbS"
        elif angle_aa <=1:
            return "fgbWC"
        elif 1.4 <= angle_aa <= 3.2:
            return "fgbH"
