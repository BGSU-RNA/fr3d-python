import numpy as np

from fr3d import definitions as defs
from fr3d.classifiers.generic import Classifier as BaseClassifier
from fr3d.data import Atom
from fr3d.data import Component
#from fr3d.data import angle_between_normals

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

        transformation_matrix = first.standard_transformation()
        #print first.unit_id(), second.unit_id(),transformation_matrix
        
        """#Test for transformation
        atoms = [
            Atom(name='N9', x=3.0, y=3.0, z=3.0),
            Atom(name='C4', x=2.0, y=2.0, z=2.0),
            Atom(name='N3', x=1.0, y=1.0, z=1.0),
        ]
        self.residue = Component(atoms, type='rna', pdb='1GID', model=1,
                                 chain='A', sequence='C', number=50,
                                 symmetry='6_555')
        trans = np.array([[1.0, 0.0, 0.0, 0.0],
                          [0.0, -1.0, 0.0, 97.240],
                          [0.0, 0.0, -1.0, 0.0],
                          [0.0, 0.0, 0.0, 1.0]])
        
        residue = self.residue.transform(trans)
        val = list(list(residue.atoms())[-1].coordinates())
        
        ans = [1.0, 96.240, -1.0]"""
        
        if transformation_matrix == None:
            return None
        trans_first = first.transform(transformation_matrix)
        trans_second = second.transform(transformation_matrix)
        min_xy, mean_z =  self.distance_metrics(trans_first, trans_second)
                
        if min_xy <= 10 and mean_z < 4:
            return self.classify_stacking(trans_first, trans_second)
        elif 10 < min_xy < 23 and mean_z < 2.5:
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
        """if not squared_xy_dist_list or not aa_z_list:
            print "empty XY squared list"
            min_xy = 0.0
            mean_z = 0
            return min_xy, mean_z"""
            
        min_xy = min(squared_xy_dist_list)
        mean_z = abs(np.mean(aa_z_list))
        if min_xy < 36:
            print "minimum XY dist squared", base_residue.unit_id(), aa_residue.unit_id(), min_xy
            print "mean_z:", mean_z
        return min_xy, mean_z        

    def classify_stacking(self, base_residue, aa_residue):
        stacked_aa = set(["TRP", "TYR", "PHE", "HIS", "ARG", "LYS", "LEU", "ILE",
                          "PRO", "ASN", "GLN"])
        perpendicular_aa = set (["HIS", "ARG", "LYS", "ASN", "GLN", "LEU", "ILE"])
        perpendicular_stack_aa = set(["PHE", "TYR"])
        angle = base_residue.angle_between_normals(aa_residue)
        if aa_residue.sequence in stacked_aa:
            if angle <= 0.67 or 2.45 <= angle <= 3.15:
                edge = self.detect_edge(base_residue, aa_residue)
                return ("stacked", edge)
            elif 1.2<= angle <=1.64:
                if aa_residue.sequence in perpendicular_stack_aa:
                    return "perpendicular stacking"
                elif aa_residue.sequence in perpendicular_aa:
                    return "cation-pi"
                
    
    def classify_pairing(self, first, second):
        pseudopair_aa = set (["ASP", "GLU", "ASN", "GLN", "HIS", "LYS", "ARG",
                      "SER", "TYR", "TRP", "PHE", "VAL", "LEU", "ILE",
                      "MET"])
        if second.sequence in pseudopair_aa:
            
            angle = first.angle_between_normals(second)
            print "pairing first second", first, second, angle
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

          