from fr3d.data.base import EntitySelector
from fr3d.data.base import AtomProxy
from fr3d.data.atoms import Atom
from fr3d import definitions as defs
from fr3d.geometry.superpositions import besttransformation
from fr3d.geometry import angleofrotation as angrot
import numpy as np
import sys
from fr3d.unit_ids import encode

from fr3d.data.mapping import parent_atom_to_modified
from fr3d.data.mapping import modified_atom_to_parent
from fr3d.data.mapping import modified_base_to_parent
from fr3d.data.mapping import modified_base_atom_list

NHBondLength=1

def unit_vector(v):
    return v / np.linalg.norm(v)

# This function calculates an angle from 0 to 180 degrees between two vectors
def angle_between_vectors(vec1, vec2):
    if len(vec1) == 3 and len(vec2) == 3:
        cosang = np.dot(vec1, vec2)
        sinang = np.linalg.norm(np.cross(vec1, vec2))
        angle = np.arctan2(sinang, cosang)
        return 180*angle/np.pi
    else:
        return None

# This function calculates an angle from 0 to 180 degrees between two vectors
def angle_between_three_points(P1,P2,P3):
    if len(P1) == 3 and len(P2) == 3 and len(P3) == 3:
        return angle_between_vectors(P1-P2,P3-P2)
    else:
        return None

# return positions of hydrogens making a tetrahedron with center C and vertices P1 and P2
def pyramidal_hydrogens(P1,C,P2,bondLength=1):

    # infer positions one way
    V1 = P1
    V2 = P2
    # vector from V2 to C
    u = unit_vector(C-V2)
    # construct Rodrigues rotation matrix
    # matrix to rotate 120 degrees around vector u
    W = np.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])
    R = np.identity(3) + (np.sqrt(3)/2)*W + 1.5 * np.dot(W,W)
    # remaining vertices are vector from C to V1 rotated 120 degrees in either direction
    V3 = C + bondLength * unit_vector(np.dot(R,V1-C))
    V4 = C + bondLength * unit_vector(np.dot(np.transpose(R),V1-C))

    # infer positions the other way
    V1 = P2
    V2 = P1
    # vector from V2 to C
    u = unit_vector(C-V2)
    # construct Rodrigues rotation matrix
    # matrix to rotate 120 degrees around vector u
    W = np.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])
    R = np.identity(3) + (np.sqrt(3)/2)*W + 1.5 * np.dot(W,W)
    # remaining vertices are vector from C to V1 rotated 120 degrees in either direction
    VV4 = C + bondLength * unit_vector(np.dot(R,V1-C))
    VV3 = C + bondLength * unit_vector(np.dot(np.transpose(R),V1-C))

    # average the two inferred positions
    P3 = (V3+VV3)/2
    P4 = (V4+VV4)/2

    # diagnostics if desired
    if 1 > 10:
        print("pyramidal_hydrogens V3, VV3 distance",np.linalg.norm(V3-VV3))
        print("pyramidal_hydrogens V4, VV4 distance",np.linalg.norm(V4-VV4))
        print("pyramidal_hydrogens P1-C-P2 angle",angle_between_three_points(P1,C,P2),"as input")
        print("pyramidal_hydrogens P1-C-P3 angle",angle_between_three_points(P1,C,P3),"as inferred")
        print("pyramidal_hydrogens P1-C-P4 angle",angle_between_three_points(P1,C,P4),"as inferred")
        print("pyramidal_hydrogens P2-C-P3 angle",angle_between_three_points(P2,C,P3),"as inferred")
        print("pyramidal_hydrogens P2-C-P4 angle",angle_between_three_points(P2,C,P4),"as inferred")
        print("pyramidal_hydrogens P3-C-P4 angle",angle_between_three_points(P3,C,P4),"as inferred")
    return P3, P4

def planar_hydrogens(P1,P2,P3,bondLength=1):

    A = unit_vector(P2 - P1)
    A1 = P3 + A*bondLength
    B = unit_vector(P3 - P2)
    # Added unit_vector for A->B
    A2 = P3 + unit_vector(B - A)*bondLength

    # diagnostics if desired
    if 1 > 10:
        print("planar_hydrogens P1-P2-P3 angle",angle_between_three_points(P1,P2,P3),"as input")
        print("planar_hydrogens P2-P3-A1 angle",angle_between_three_points(P2,P3,A1),"as inferred")
        print("planar_hydrogens P2-P3-A2 angle",angle_between_three_points(P2,P3,A2),"as inferred")
    return A1, A2

def planar_ring_hydrogen(P1,P2,P3,bondlength=1):

    # vectors P1->P2 and P3->P2
    u=unit_vector(P2-P1)
    v=unit_vector(P2-P3)

    # adding the hydrogens
    w = unit_vector(u + v)
    A1= P2 + bondlength * w

    return A1



class Component(EntitySelector):
    """This represents things like nucleic acids, amino acids, small molecules
    and ligands.
    """

    def __init__(self, atoms, pdb=None, model=None, type=None, chain=None,
                 symmetry=None, sequence=None, number=None, index=None,
                 insertion_code=None, polymeric=None, alt_id=None):
        """Create a new Component.

        :atoms: The atoms this component is composed of.
        :pdb: The pdb this is a part of.
        :model: The model number.
        """

        self._atoms = atoms
        self.pdb = pdb
        self.model = model
        self.type = type
        self.chain = chain
        self.symmetry = symmetry
        self.sequence = sequence
        self.number = number
        self.index = index
        self.insertion_code = insertion_code
        self.polymeric = polymeric
        self.alt_id = alt_id
        self.base_center = None
        self.rotation_matrix = None

        # for bases, calculate and store rotation_matrix
        # calculate and store base_center; especially for modified nt without all heavy atoms
        self.calculate_rotation_matrix()

        # # initialize centers so they can be used to infer hydrogens
        # self.centers = AtomProxy(self._atoms)

        # add hydrogen atoms to standard bases and amino acids
        self.infer_NA_hydrogens()

        # do not routinely add hydrogen atoms to amino acids
        # self.infer_amino_acid_hydrogens()

        # initialize centers again to include hydrogens
        self.centers = AtomProxy(self._atoms)

        # standard and modified bases should have their rotation matrix
        # calculated already, and should have a base center set by that
        # if they don't, there is no sensible way to assign a base center,
        # for example if the structure just has a backbone trace
        if self.base_center is not None:
            self.centers.setcenter('base',self.base_center)

            if self.sequence in ['A','G','DA','DG']:
                self.centers.define('glycosidic',['N9'])
            elif self.sequence in ['C','U','DC','DT']:
                self.centers.define('glycosidic',['N1'])
            elif self.sequence in modified_base_to_parent:
                if modified_base_to_parent[self.sequence] in ['A','G','DA','DG']:
                    self.centers.define('glycosidic',[parent_atom_to_modified[self.sequence]['N9']])
                elif modified_base_to_parent[self.sequence] in ['C','U','DC','DT']:
                    self.centers.define('glycosidic',parent_atom_to_modified[self.sequence]['N1'])

        if self.sequence in defs.nt_sugar:
            atoms = defs.nt_sugar[self.sequence]
            self.centers.define('nt_sugar', atoms)

        if self.sequence in defs.nt_phosphate:
            atoms = defs.nt_phosphate[self.sequence]
            self.centers.define('nt_phosphate', atoms)

        # attempt to add sugar and phosphate centers for all modified nucleotides
        if self.sequence in modified_base_to_parent:
            atoms = defs.nt_sugar['A']
            self.centers.define('nt_sugar', atoms)
            atoms = defs.nt_phosphate['A']
            self.centers.define('nt_phosphate', atoms)

        if self.sequence in defs.aa_fg:
            atoms = defs.aa_fg[self.sequence]
            self.centers.define('aa_fg', atoms)

        if self.sequence in defs.aa_backbone:
            atoms = defs.aa_backbone[self.sequence]
            self.centers.define('aa_backbone', atoms)

    def atoms(self, **kwargs):
        """Get, filter and sort the atoms in this component. Access is as
        described by EntitySelector.

        :kwargs: The keyword arguments to filter and sort by.
        :returns: A list of the requested atoms.
        """

        name = kwargs.get('name')
        
        if sys.version_info[0] < 3:
            if isinstance(name, basestring):
                definition = self.centers.definition(name)
                if definition:
                    kwargs['name'] = definition
        else:
            if isinstance(name, str):
                definition = self.centers.definition(name)
                if definition:
                    kwargs['name'] = definition
            
            
        return EntitySelector(self._atoms, **kwargs)

    def coordinates(self, **kwargs):
        """Get the coordaintes of all atoms in this component. This will
        filter to the requested atoms, sort and then provide a numpy array
        of all coordinates for the given atoms.

        :kwargs: Arguments to filter and sort by.
        :returns: A numpy array of the coordinates.
        """
        return np.array([atom.coordinates() for atom in self.atoms(**kwargs)])

    def select(self, **kwargs):
        """Select a group of atoms to create a new component out of.

        :kwargs: As for atoms.
        :returns: A new Component
        """
        return Component(list(self.atoms(**kwargs)),
                         pdb=self.pdb,
                         model=self.model,
                         type=self.type,
                         chain=self.chain,
                         symmetry=self.symmetry,
                         sequence=self.sequence,
                         number=self.number,
                         index=self.index,
                         insertion_code=self.insertion_code,
                         alt_id=self.alt_id,
                         polymeric=self.polymeric,
                         inferhydrogens=False)

    def is_complete(self, names, key='name'):
        """This checks if we can find all atoms in this entity with the given
        names. This assumes that the names for each atom are unique. If you
        wish to use something other than name use the key argument. However, it
        must provide a unique value for each atom, if several atoms with the
        same value are found will cause the function to behave oddly.

        :names: The list of names to check for.
        :key: The key to use for getting atoms. Defaults to name.
        :returns: True if all atoms with the given name are present.
        """
        kwargs = {key: names}
        found = list(self.atoms(**kwargs))
        return len(found) == len(names)

    def calculate_rotation_matrix(self):
        """Calculate a rotation matrix that will rotate the atoms in an RNA
        base into a standard orientation in the xy plane with the Watson-
        Crick edge in the positive x and y quadrant.
        """

        if self.sequence not in defs.NAbaseheavyatoms and \
                self.sequence not in modified_base_to_parent:
            return None

        R = []   # 3d coordinates of observed base
        S = []   # 3d coordinates of standard base in xy plane

        if self.sequence in defs.NAbaseheavyatoms:
            baseheavy = defs.NAbaseheavyatoms[self.sequence]
            standard_coords = defs.NAbasecoordinates[self.sequence]
            for atom in self.atoms(name=baseheavy):
                R.append(atom.coordinates())
                S.append(standard_coords[atom.name])

        elif self.sequence in modified_base_to_parent:
            # get the mapping from this modified nucleotide to its parent
            #mod_to_parent = defs.modified_nucleotides[self.sequence]
            # print(mod_to_parent)
            # get the standard coordinates for the parent nucleotide
            standard_coords = defs.NAbasecoordinates[modified_base_to_parent[self.sequence]]
            # loop over mapped base atoms in the modified nucleotide
            # parent base atom
            for atom in self.atoms(name=list(modified_base_atom_list[self.sequence])):
                #redundant. Last check should be sufficient. 
                if atom.name in defs.NAbasecoordinates[modified_base_to_parent[self.sequence]]:
                    R.append(atom.coordinates())
                    #print(self.sequence)
                    S.append(standard_coords[parent_atom_to_modified[self.sequence][atom.name]])

            """
            for mod_atom_name, parent_atom_name in mod_to_parent["atoms"].items():
                print(mod_atom_name,parent_atom_name)
                atom = self.atoms(name=[mod_atom_name])
                R.append(atom.coordinates())
                S.append(standard_coords[parent_atom_name])
            """

        R = np.array(R)
        R = R.astype(np.float)
        S = np.array(S)
        S = S.astype(np.float)

        try:
            rotation_matrix, fitted, meanR, rmsd, sse, meanS = \
                besttransformation(R, S)
        except:
            if len(R) != len(S):
                print("%s Rotation matrix calculation failed, sizes %d and %d" % (self.unit_id(),len(R),len(S)))
            elif len(R) < 3:
                print("%s Rotation matrix calculation failed, %d new atoms" % (self.unit_id(),len(R)))
            elif len(S) < 3:
                print("%s Rotation matrix calculation failed, %d standard atoms" % (self.unit_id(),len(S)))
            else:
                print("%s Rotation matrix calculation failed, not sure why" % self.unit_id())

            return None

        self.rotation_matrix = rotation_matrix

        # map the origin out to where the center of the base should be
        if self.sequence in defs.NAbaseheavyatoms:
            # standard bases are designed to have meanS zero; less work
            self.base_center = meanR
        else:
            # some modified bases are missing some heavy atoms, meanS not zero
            # this comes out as a numpy matrix?  different than meanR above
            base_center = np.subtract(meanR,np.dot(rotation_matrix,meanS))
            self.base_center = np.array([base_center[0,0],base_center[0,1],base_center[0,2]])

        """ For the life of me, I could not figure out any other way of
        converting base_center from a 2d array to a 1d array.
        Taking a slice did not work, reshape did not work, etc.
        This is crucially important; a 2d array won't be written to the
        unit_centers table in the database.
        """


    def infer_NA_hydrogens(self):
        """
        Infer the coordinates of the hydrogen atoms of this component.
        RNA and DNA work, and amino acids are being added.
        This code only adds hydrogens with their name, but the Atom entity
        does not have the full unit ID of the atom, unlike the heavy atoms
        taken from the CIF file.
        """

        def get_amino_hydrogen_coords(self, heavy, amino1, amino2):
            """
            Helper function to retrieve the coordinates of amino hydrogens and specified heavy atom.
            Seperate processing for modified nucleotides. 
            If modified nucleotide has a mapping, checks to see if the atom maps to the name of the passed in parent.
            Returns 3 triples of atom coordinates of a (heavy_atom, amino_hydrogen_#1, amino_hydrogen_#2)
            """
            amino1coords = None
            amino2coords = None
            # Standard base
            if self.sequence in defs.NAbasehydrogens:         
                for atom in self._atoms:
                    if atom.name == heavy:
                        heavy = (atom.x, atom.y, atom.z)
                    elif atom.name == amino1:
                        amino1coords = (atom.x, atom.y, atom.z)
                    elif atom.name == amino2:
                        amino2coords = (atom.x, atom.y, atom.z)
            # Mapped modified base, make sure it is in the mappings before going
            elif self.sequence in modified_base_to_parent:
                for atom in self._atoms:
                    if atom.name in modified_base_atom_list[self.sequence]: # Weed out backbone atoms
                        if parent_atom_to_modified[self.sequence][atom.name] == heavy: # parent_atom_to_modified[PSU][C5] would return N1 of parent 
                            heavy = (atom.x, atom.y, atom.z)
                        elif parent_atom_to_modified[self.sequence][atom.name] == amino1:
                            amino1coords = (atom.x, atom.y, atom.z)
                        elif parent_atom_to_modified[self.sequence][atom.name] == amino2:
                            amino2coords = (atom.x, atom.y, atom.z)
            return heavy, amino1coords, amino2coords
        
        amino1coords = None
        amino2coords = None
        heavy = None
        dist1 = 0
        dist2 = 0

        try:
            # going to add or fix hydrogens for standard bases (only on base)
            if self.sequence in defs.NAbasehydrogens:
                hydrogens = set(defs.NAbasehydrogens[self.sequence]) # All hydrogens that should be present on this base
                already = set([atom.name for atom in self._atoms]) # hydrogens that are already observed in the 3D structure
                hydrogens = hydrogens - already 
                if len(already) > 0: #that means high enough resolution that hydrogens don't need to be infered
                    #check distances depending on sequence to make sure labelled correct.
                    if self.sequence == "A" or self.sequence == "DA":
                        heavy = "N7" # H62 is closest to N7
                        amino1 = "H61"
                        amino2 = "H62"
                        heavy, amino1coords, amino2coords = get_amino_hydrogen_coords(self, heavy, amino1, amino2)
                    elif self.sequence == "C" or self.sequence == "DC":
                        heavy = "C5"
                        amino1 = "H41"
                        amino2 = "H42"
                        heavy, amino1coords, amino2coords = get_amino_hydrogen_coords(self, heavy, amino1, amino2)
                    elif self.sequence == "G" or self.sequence == "DG":
                        heavy = "N1"
                        amino1 = "H21"
                        amino2 = "H22"
                        heavy, amino1coords, amino2coords = get_amino_hydrogen_coords(self, heavy, amino1, amino2)
                    elif self.sequence == "DT" or self.sequence == "U":
                        # No amino hydrogens...helps to skip processing
                        heavy = None
                        amino1coords = None
                        amino2coords = None
                    if amino1coords and amino2coords and heavy:
                        if len(amino1coords) == 3 and len(amino2coords) == 3 and len(heavy) == 3:
                            # calculate distances
                            dist1 = (heavy[0]-amino1coords[0])**2 + (heavy[1]-amino1coords[1])**2 + (heavy[2]-amino1coords[2])**2
                            dist2 = (heavy[0]-amino2coords[0])**2 + (heavy[1]-amino2coords[1])**2 + (heavy[2]-amino2coords[2])**2
                            if dist1 < dist2: 
                            #H62 should be closer to N7 than H61, H42 should be closer to C5 than H41
                                for atom in self._atoms:
                                #switch coordinates of amino hydrogens
                                    if atom.name == amino1:
                                        atom.x = amino2coords[0]
                                        atom.y = amino2coords[1]
                                        atom.z = amino2coords[2]
                                    elif atom.name == amino2:
                                        atom.x = amino1coords[0]
                                        atom.y = amino1coords[1] 
                                        atom.z = amino1coords[2]

                # If hydrogens are missing, add missing hydrogens
                if len(hydrogens) > 0: 
                    coordinates = defs.NAbasecoordinates[self.sequence]

                    for hydrogenatom in hydrogens:
                        hydrogencoordinates = coordinates[hydrogenatom] # standard coordinates of hydrogen atoms
                        newcoordinates = self.base_center + \
                            np.dot(self.rotation_matrix, hydrogencoordinates)
                        self._atoms.append(Atom(name=hydrogenatom,
                                                x=newcoordinates[0, 0],
                                                y=newcoordinates[0, 1],
                                                z=newcoordinates[0, 2]))

            # repeat similar logic but for modified nucleotides. Written out twice so that modified nucleotides logic doesn't slow down normal bases as they're much less frequent.
            elif self.sequence in modified_base_to_parent:   
                already = []
                from fr3d.data.mapping import modified_atom_map
                from fr3d.data.mapping import modified_hydrogens
                from fr3d.data.mapping import modified_hydrogens_coordinates
                for atom in self._atoms:
                    if 'H' in atom.name: 
                        already.append(atom)
                # If hydrogens are already observed, make sure they're named correctly according to their mappings and their parents naming conventions
                if len(already) > 0:
                    if modified_base_to_parent[self.sequence] == 'A' or modified_base_to_parent[self.sequence] == 'DA':
                        heavy = "N7" # H62 is closest to N7
                        amino1 = "H61"
                        amino2 = "H62"
                        heavy, amino1coords, amino2coords = get_amino_hydrogen_coords(self, heavy, amino1, amino2)
                    elif modified_base_to_parent[self.sequence] == "C" or modified_base_to_parent[self.sequence] == "DC":
                        heavy = "C5"
                        amino1 = "H41"
                        amino2 = "H42"
                        heavy, amino1coords, amino2coords = get_amino_hydrogen_coords(self, heavy, amino1, amino2)
                    elif modified_base_to_parent[self.sequence] == "G" or modified_base_to_parent[self.sequence] == "DG":
                        heavy = "N1"
                        amino1 = "H21"
                        amino2 = "H22"
                        heavy, amino1coords, amino2coords = get_amino_hydrogen_coords(self, heavy, amino1, amino2)
                    if amino1coords and amino2coords and heavy:
                        if len(amino1coords) == 3 and len(amino2coords) == 3 and len(heavy) == 3:
                            # Calculate distance
                            dist1 = (heavy[0]-amino1coords[0])**2 + (heavy[1]-amino1coords[1])**2 + (heavy[2]-amino1coords[2])**2
                            dist2 = (heavy[0]-amino2coords[0])**2 + (heavy[1]-amino2coords[1])**2 + (heavy[2]-amino2coords[2])**2
                            #H62 should be closer to N7 than H61, H42 should be closer to C5 than H41
                            if dist1 < dist2: 
                                #switch coordinates of amino hydrogens
                                for atom in self._atoms:
                                    if atom.name in modified_base_atom_list[self.sequence]:
                                        if parent_atom_to_modified[self.sequence][atom.name] == amino1:
                                            atom.x = amino2coords[0]
                                            atom.y = amino2coords[1]
                                            atom.z = amino2coords[2]
                                        elif parent_atom_to_modified[self.sequence][atom.name] == amino2:
                                            atom.x = amino1coords[0]
                                            atom.y = amino1coords[1] 
                                            atom.z = amino1coords[2]
                else: # hydrogens aren't observed already. If the heavy atom has a mapping to the parent, infer hydrogen with parent hydrogens coordinates with modified hydrogens name
                    hydrogens = modified_hydrogens[self.sequence]
                    coordinates = modified_hydrogens_coordinates[self.sequence]
                    for hydrogenatom in hydrogens:
                        hydrogencoordinates = coordinates[hydrogenatom]
                        newcoordinates = self.base_center + \
                            np.dot(self.rotation_matrix, hydrogencoordinates)
                        self._atoms.append(Atom(name=hydrogenatom,
                                                x=newcoordinates[0, 0],
                                                y=newcoordinates[0, 1],
                                                z=newcoordinates[0, 2]))

        except:
            print("%s Adding hydrogens failed" % self.unit_id())


    def infer_amino_acid_hydrogens(self):
        """
        Infer the coordinates of the hydrogen atoms of this component.
        RNA and DNA work, and amino acids are being added.
        This code only adds hydrogens with their name, but the Atom entity
        does not have the full unit ID of the atom, unlike the heavy atoms
        taken from the CIF file.
        """
        try:

            if self.sequence == "ALA":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"],NHBondLength)
                self._atoms.append(Atom(name="HB1",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["HB1"])
                self._atoms.append(Atom(name="HB3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB2",x=A2[0],y=A2[1],z=A2[2]))

            elif self.sequence == "ARG":

                A1,A2 = planar_hydrogens(self.centers["NE"],self.centers["CZ"],self.centers["NH1"],NHBondLength)
                self._atoms.append(Atom(name="HH11",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HH12",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["NE"],self.centers["CZ"],self.centers["NH2"])
                self._atoms.append(Atom(name="HH22",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HH21",x=A2[0],y=A2[1],z=A2[2]))

                # the following lines assume that NH2 is closer to CD
                # than NH1 is however, some structures aren't labeled
                # that way, and so either A1 or A2 could be the correct
                # location for HE
                A1,A2 = planar_hydrogens(self.centers["NH1"],self.centers["CZ"],self.centers["NE"])
                self._atoms.append(Atom(name="HE",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CG"],self.centers["CD"],self.centers["NE"])
                self._atoms.append(Atom(name="HD3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HD2",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CB"],self.centers["CG"],self.centers["CD"])
                self._atoms.append(Atom(name="HG2",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HG3",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG"])
                self._atoms.append(Atom(name="HB2",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB3",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

            elif self.sequence == "ASN":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG"])
                self._atoms.append(Atom(name="HB3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB2",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["CB"],self.centers["CG"],self.centers["ND2"],NHBondLength)
                self._atoms.append(Atom(name="HD22",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = planar_hydrogens(self.centers["OD1"],self.centers["CG"],self.centers["ND2"],NHBondLength)
                self._atoms.append(Atom(name="HD21",x=A1[0],y=A1[1],z=A1[2]))

            elif self.sequence == "ASP":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG"])
                self._atoms.append(Atom(name="HB3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB2",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["CB"],self.centers["CG"],self.centers["OD2"],NHBondLength)
                self._atoms.append(Atom(name="HD2",x=A1[0],y=A1[1],z=A1[2]))

            elif self.sequence == "CYS":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["SG"])
                self._atoms.append(Atom(name="HB3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB2",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["SG"],NHBondLength)
                self._atoms.append(Atom(name="HG",x=A1[0],y=A1[1],z=A1[2]))

            elif self.sequence == "GLU":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG"])
                self._atoms.append(Atom(name="HB3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB2",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CB"],self.centers["CG"],self.centers["CD"])
                self._atoms.append(Atom(name="HG3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HG2",x=A2[0],y=A2[1],z=A2[2]))

            elif self.sequence == "GLY":

                A1,A2 = pyramidal_hydrogens(self.centers["N"],self.centers["CA"],self.centers["C"])
                self._atoms.append(Atom(name="HA3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HA2",x=A2[0],y=A2[1],z=A2[2]))

            elif self.sequence == "HIS":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG"])
                self._atoms.append(Atom(name="HB3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB2",x=A2[0],y=A2[1],z=A2[2]))

                A1 = planar_ring_hydrogen(self.centers["CG"],self.centers["ND1"],self.centers["CE1"],NHBondLength)
                self._atoms.append(Atom(name="HD1",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["NE2"],self.centers["CE1"],self.centers["ND1"],NHBondLength)
                self._atoms.append(Atom(name="HE1",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["CE1"],self.centers["NE2"],self.centers["CD2"],NHBondLength)
                self._atoms.append(Atom(name="HE2",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["NE2"],self.centers["CD2"],self.centers["CG"],NHBondLength)
                self._atoms.append(Atom(name="HD2",x=A1[0],y=A1[1],z=A1[2]))

            elif self.sequence == "ILE":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CB"],self.centers["CG1"],self.centers["CD1"])
                self._atoms.append(Atom(name="HG12",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HG13",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["CG1"],self.centers["CB"],self.centers["CG2"],NHBondLength)
                self._atoms.append(Atom(name="HG23",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CB"],self.centers["CG2"],self.centers["HG23"])
                self._atoms.append(Atom(name="HG22",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HG21",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["CB"],self.centers["CG1"],self.centers["CD1"],NHBondLength)
                self._atoms.append(Atom(name="HD11",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CG1"],self.centers["CD1"],self.centers["HD11"])
                self._atoms.append(Atom(name="HD12",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HD13",x=A2[0],y=A2[1],z=A2[2]))

            elif self.sequence == "LEU":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG"])
                self._atoms.append(Atom(name="HB2",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB3",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["CA"],self.centers["N"],self.centers["CG"])
                self._atoms.append(Atom(name="HG",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = planar_hydrogens(self.centers["CB"],self.centers["HB3"],self.centers["CD1"])
                self._atoms.append(Atom(name="HD12",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CG"],self.centers["CD1"],self.centers["HD12"])
                self._atoms.append(Atom(name="HD11",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HD13",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["CB"],self.centers["HB2"],self.centers["CD2"])
                self._atoms.append(Atom(name="HD21",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CG"],self.centers["CD2"],self.centers["HD21"])
                self._atoms.append(Atom(name="HD22",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HD23",x=A2[0],y=A2[1],z=A2[2]))

            elif self.sequence == "LYS":

                A1,A2 = pyramidal_hydrogens(self.centers["CG"],self.centers["CD"],self.centers["CE"])
                self._atoms.append(Atom(name="HD3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HD2",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CB"],self.centers["CG"],self.centers["CD"])
                self._atoms.append(Atom(name="HG3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HG2",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CD"],self.centers["CE"],self.centers["NZ"])
                self._atoms.append(Atom(name="HE3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HE2",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG"])
                self._atoms.append(Atom(name="HB3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB2",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["CD"],self.centers["CE"],self.centers["NZ"],NHBondLength)
                self._atoms.append(Atom(name="HZ3",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CE"],self.centers["NZ"],self.centers["HZ3"])
                self._atoms.append(Atom(name="HZ2",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HZ1",x=A2[0],y=A2[1],z=A2[2]))

            elif self.sequence == "MET":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG"])
                self._atoms.append(Atom(name="HB3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB2",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CB"],self.centers["CG"],self.centers["SD"])
                self._atoms.append(Atom(name="HG3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HG2",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["CG"],self.centers["SD"],self.centers["CE"],NHBondLength)
                self._atoms.append(Atom(name="HE1",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["SD"],self.centers["CE"],self.centers["HE1"])
                self._atoms.append(Atom(name="HE3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HE2",x=A2[0],y=A2[1],z=A2[2]))

            elif self.sequence == "PHE":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG"])
                self._atoms.append(Atom(name="HB3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB2",x=A2[0],y=A2[1],z=A2[2]))

                A1 = planar_ring_hydrogen(self.centers["CG"],self.centers["CD1"],self.centers["CE1"],NHBondLength)
                self._atoms.append(Atom(name="HD1",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["CD1"],self.centers["CE1"],self.centers["CZ"],NHBondLength)
                self._atoms.append(Atom(name="HE1",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["CE1"],self.centers["CZ"],self.centers["CE2"],NHBondLength)
                self._atoms.append(Atom(name="HZ",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["CZ"],self.centers["CE2"],self.centers["CD2"],NHBondLength)
                self._atoms.append(Atom(name="HE2",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["CG"],self.centers["CD2"],self.centers["CE2"],NHBondLength)
                self._atoms.append(Atom(name="HD2",x=A1[0],y=A1[1],z=A1[2]))

            elif self.sequence == "PRO":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["N"],self.centers["CD"])
                self._atoms.append(Atom(name="H",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["N"],self.centers["CD"],self.centers["CG"])
                self._atoms.append(Atom(name="HD2",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HD3",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CD"],self.centers["CG"],self.centers["CB"])
                self._atoms.append(Atom(name="HG2",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HG3",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG"])
                self._atoms.append(Atom(name="HB3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB2",x=A2[0],y=A2[1],z=A2[2]))

            elif self.sequence == "SER":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["OG"])
                self._atoms.append(Atom(name="HB3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB2",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["OG"],NHBondLength)
                self._atoms.append(Atom(name="HG",x=A1[0],y=A1[1],z=A1[2]))

            elif self.sequence == "THR":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG2"])
                self._atoms.append(Atom(name="HB",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = planar_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG2"])
                self._atoms.append(Atom(name="HG21",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CB"],self.centers["CG2"],self.centers["HG21"])
                self._atoms.append(Atom(name="HG23",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HG22",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["CG2"],self.centers["HG23"],self.centers["OG1"])
                self._atoms.append(Atom(name="HG1",x=A1[0],y=A1[1],z=A1[2]))

            elif self.sequence == "TRP":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG"])
                self._atoms.append(Atom(name="HB3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB2",x=A2[0],y=A2[1],z=A2[2]))

                A1 = planar_ring_hydrogen(self.centers["CG"],self.centers["CD1"],self.centers["NE1"],NHBondLength)
                self._atoms.append(Atom(name="HD1",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["CD1"],self.centers["NE1"],self.centers["CE2"],NHBondLength)
                self._atoms.append(Atom(name="HE1",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["CE2"],self.centers["CZ2"],self.centers["CH2"],NHBondLength)
                self._atoms.append(Atom(name="HZ2",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["CZ2"],self.centers["CH2"],self.centers["CZ3"],NHBondLength)
                self._atoms.append(Atom(name="HH2",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["CH2"],self.centers["CZ3"],self.centers["CE3"],NHBondLength)
                self._atoms.append(Atom(name="HZ3",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["CZ3"],self.centers["CE3"],self.centers["CD2"],NHBondLength)
                self._atoms.append(Atom(name="HE3",x=A1[0],y=A1[1],z=A1[2]))

            elif self.sequence == "TYR":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG"])
                self._atoms.append(Atom(name="HB3",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HB2",x=A2[0],y=A2[1],z=A2[2]))

                A1 = planar_ring_hydrogen(self.centers["CG"],self.centers["CD2"],self.centers["CE2"],NHBondLength)
                self._atoms.append(Atom(name="HD2",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["CD2"],self.centers["CE2"],self.centers["CZ"],NHBondLength)
                self._atoms.append(Atom(name="HE2",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["CZ"],self.centers["CE1"],self.centers["CD1"],NHBondLength)
                self._atoms.append(Atom(name="HE1",x=A1[0],y=A1[1],z=A1[2]))

                A1 = planar_ring_hydrogen(self.centers["CG"],self.centers["CD1"],self.centers["CE1"],NHBondLength)
                self._atoms.append(Atom(name="HD1",x=A1[0],y=A1[1],z=A1[2]))

            elif self.sequence == "VAL":

                A1,A2 = pyramidal_hydrogens(self.centers["C"],self.centers["CA"],self.centers["CB"])
                self._atoms.append(Atom(name="HA",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG1"])
                self._atoms.append(Atom(name="HB",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG1"],NHBondLength)
                self._atoms.append(Atom(name="HG11",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CB"],self.centers["CG1"],self.centers["HG11"])
                self._atoms.append(Atom(name="HG12",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HG13",x=A2[0],y=A2[1],z=A2[2]))

                A1,A2 = planar_hydrogens(self.centers["CA"],self.centers["CB"],self.centers["CG2"],NHBondLength)
                self._atoms.append(Atom(name="HG23",x=A1[0],y=A1[1],z=A1[2]))

                A1,A2 = pyramidal_hydrogens(self.centers["CB"],self.centers["CG2"],self.centers["HG23"])
                self._atoms.append(Atom(name="HG21",x=A1[0],y=A1[1],z=A1[2]))
                self._atoms.append(Atom(name="HG22",x=A2[0],y=A2[1],z=A2[2]))

        except:
                print("%s Adding hydrogens failed" % self.unit_id())


    def transform(self, transform_matrix):
        """Create a new component from "self" by applying the 4x4 transformation
        matrix. This does not keep the rotation matrix if any,
        but will keep any added hydrogens.

        :transform_matrix: The 4x4 transformation matrix to apply.
        :returns: A new Component with the same properties except with
        transformed coordinates.
        """

        atoms = [atom.transform(transform_matrix) for atom in self.atoms()]
        comp = Component(atoms, pdb=self.pdb,
                         model=self.model,
                         type=self.type,
                         chain=self.chain,
                         symmetry=self.symmetry,
                         sequence=self.sequence,
                         number=self.number,
                         index=self.index,
                         insertion_code=self.insertion_code,
                         alt_id=self.alt_id,
                         polymeric=self.polymeric)
        return comp

    def translate_rotate(self, residue):
        reference = self.centers["base"]
        rotation = self.rotation_matrix
        for atom in residue.atoms():
            atom_coord = atom.coordinates()
            dist_translate = np.subtract(atom_coord, reference)
            dist_aa_matrix = np.matrix(dist_translate)
            rotated_atom = dist_aa_matrix * rotation
            coord_array = np.array(rotated_atom)
            a = coord_array.flatten()
            transformed_coord = a.tolist()
        return transformed_coord

    def translate_rotate_component(self, component):
        """Translate and rotate the atoms in component according to
        the translation and rotation that will bring self to standard
        position at the origin.
        :param Component residue:  the residue to move
        :returns Component newcomp
        """

        atoms = [self.translate_rotate_atom(atom) for atom in component.atoms()]
        newcomp = Component(atoms, pdb=component.pdb,
                         model=component.model,
                         type=component.type,
                         chain=component.chain,
                         symmetry=component.symmetry,
                         sequence=component.sequence,
                         number=component.number,
                         index=component.index,
                         insertion_code=component.insertion_code,
                         alt_id=component.alt_id,
                         polymeric=component.polymeric)
        # newcomp.infer_hydrogens()  # this is not the time to do this
        return newcomp

    def translate_rotate_atom(self, atom):
        """Translate and rotate atom according to
        the translation and rotation that will bring self to standard
        position at the origin.

        :param Atom atom: The Atom to move.
        :returns Atom: The moved atom.
        """

        atom_coord = atom.coordinates()
        translated_coord = np.subtract(atom_coord, self.base_center)
        translated_coord_matrix = np.matrix(translated_coord)
        rotated_coord = translated_coord_matrix * self.rotation_matrix
        coord_array = np.array(rotated_coord)
        a = coord_array.flatten()
        x, y, z = a.tolist()
        return Atom(x=x, y=y, z=z,
                    pdb=atom.pdb,
                    model=atom.model,
                    chain=atom.chain,
                    component_id=atom.component_id,
                    component_number=atom.component_number,
                    component_index=atom.component_index,
                    insertion_code=atom.insertion_code,
                    alt_id=atom.alt_id,
                    group=atom.group,
                    type=atom.type,
                    name=atom.name,
                    symmetry=atom.symmetry,
                    polymeric=atom.polymeric)

    def standard_transformation(self):
        """Returns a 4X4 transformation matrix which can be used to transform
        any component to the same relative location as the "self" argument in
        its standard location. If this is not an RNA component then this returns
        None.
        :returns: A numpy array suitable for input to self.transform to produce
        a transformed component.
        """

        if 'base' not in self.centers:
            return None
        base_center = self.centers["base"]
        if len(base_center) == 0:
            return None
        seq = self.sequence
        dist_translate = base_center

        rotation_matrix_transpose = self.rotation_matrix.transpose()

        matrix = np.zeros((4, 4))
        matrix[0:3, 0:3] = rotation_matrix_transpose
        matrix[0:3, 3] = -np.dot(rotation_matrix_transpose, dist_translate)
        matrix[3, 3] = 1.0

        return matrix

    def translate(self, aa_residue):
        if 'base' not in self.centers:
            return None
        rotation = self.rotation_matrix
        for atom in aa_residue:
            dist_translate = np.subtract(atom, self.centers["base"])
            rotated_atom = dist_translate*rotation
            coord_array = np.array(rotated_atom)
            a = coord_array.flatten()
            coord = a.tolist()
        return coord


    def unit_id(self):
        """Compute the unit id of this Component.

        :returns: The unit id.
        """

        return encode({
            'pdb': self.pdb,
            'model': self.model,
            'chain': self.chain,
            'component_id': self.sequence,
            'component_number': self.number,
            'alt_id': self.alt_id,
            'insertion_code': self.insertion_code,
            'symmetry': self.symmetry
        })

    def atoms_within(self, other, cutoff, using=None, to=None, min_number=1):
        """Determine if there are any atoms from another component within some
        distance.

        :other: Another component to compare agains.
        :using: The atoms from this component to compare with.
        :to: The atoms from the other component to compare against.
        :cutoff: The distances atoms must be within. Default 4.0
        """

        kw1 = {}
        if using:
            kw1['name'] = using

        kw2 = {}
        if to:
            kw2['name'] = to

        n = 0
        for atom1 in self.atoms(**kw1):
            for atom2 in other.atoms(**kw2):
                if atom1.distance(atom2) <= abs(cutoff):
                    n = n+1

        if n >= min_number:
            return True

    def distance(self, other, using='*', to='*'):
        """Compute a center center distance between this and another component.

        :other: The other component to get distance to.
        :using: A list of atom names to use for this component. Defaults to '*'
        meaning all atoms.
        :to: A list of atoms names for the second component. Defaults to '*'
        meaning all atoms.
        :returns: The distance between the two centers.
        """
        coordinates = self.centers[using]
        other_coord = other.centers[to]
        distance = np.subtract(coordinates, other_coord)
        return np.linalg.norm(distance)

    def __len__(self):
        """Compute the length of this Component. This is the number of atoms in
        this residue.

        :returns: The number of atoms.
        """
        return len(self._atoms)

    def __eq__(self, other):
        return isinstance(other, Component) and \
            self.pdb == other.pdb and \
            self.model == other.model and \
            self.chain == other.chain and \
            self.symmetry == other.symmetry and \
            self.sequence == other.sequence and \
            self.number == other.number and \
            self.insertion_code == other.insertion_code and \
            self.alt_id == other.alt_id

    def __repr__(self):
        return '<Component %s>' % self.unit_id()

    def angle_between_normals(self, aa_residue):
        vec1 = self.normal_calculation()
        vec2 = aa_residue.normal_calculation()
        return angrot.angle_between_planes(vec1, vec2)

    def normal_calculation(self):
        key = self.sequence
        P1 = self.centers[defs.planar_atoms[key][0]]
        P2 = self.centers[defs.planar_atoms[key][1]]
        P3 = self.centers[defs.planar_atoms[key][2]]
        vector = np.cross((P2 - P1), (P3-P1))
        return vector

    def enough_hydrogen_bonds(self, second, min_distance=4, min_bonds=2):
        """Calculates atom to atom distance of part "aa_part" of neighboring
        amino acids of type "aa" from each atom of base. Only returns a pair
        of aa/nt if two or more atoms are within the cutoff distance.
        """

        HB_atoms = set(['N', 'NH1','NH2','NE','NZ','ND1','NE2','O','OD1','OE1','OE2', 'OG', 'OH'])
        n = 0
        base_seq = self.sequence()
        for base_atom in self.atoms(name=defs.NAbaseheavyatoms[base_seq]):
            for aa_atom in second.atoms(name=defs.aa_fg[second.sequence]):
                distance = np.subtract(base_atom.coordinates(), aa_atom.coordinates())
                distance = np.linalg.norm(distance)
                if distance <= min_distance and aa_atom.name in HB_atoms:
                    n = n + 1
                    if n > min_bonds:
                        return True
        return False

    def stacking_tilt(aa_residue, aa_coordinates):
        baa_dist_list = []

        for aa_atom in aa_residue.atoms(name=defs.aa_fg[aa_residue.sequence]):
            key = aa_atom.name
            aa_z = aa_coordinates[key][2]
            baa_dist_list.append(aa_z)
        max_baa = max(baa_dist_list)
        min_baa = min(baa_dist_list)
        diff = max_baa - min_baa
        #print aa_residue.unit_id(), diff
        if diff <= defs.tilt_cutoff[aa_residue.sequence]:
            return "stacked"
