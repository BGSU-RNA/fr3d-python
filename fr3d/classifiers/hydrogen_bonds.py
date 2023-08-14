# This file holds methods for checking hydrogen bonds

import csv
import numpy as np
import os
from fr3d.data.mapping import modified_base_atom_list,parent_atom_to_modified,modified_atom_to_parent,modified_base_to_parent

def load_ideal_basepair_hydrogen_bonds():
    """
    Load a table of ideal hydrogen bonds for each base combination
    in each Leontis-Westhof family, as determined by Jesse Stombaugh.
    For each base combination, store as a mapping from hydrogen bonds
    to basepairs that contain that bond.
    """

    hbond = {}

    current_path,current_program = os.path.split(os.path.abspath(__file__))

    filename = os.path.join(current_path,'H_bonding_Atoms_from_Isostericity_Table.csv')

    data_found = False

    with open(filename) as csvfile:
        bond_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        for row in bond_reader:

            if row and row[0] == 'TYPE':
                data_found = True

            if len(row) == 9 and data_found:
                b1 = row[1]
                b2 = row[2]
                if b1 in ['A','C','G','U'] and b2 in ['A','C','G','U']:
                    combination = b1+','+b2
                    LW = row[0].replace('cis ','c').replace('trans ','t').replace('/','')
                    if combination in ['A,A','A,C','A,G','A,U','C,C','G,C','G,G','G,U','C,U','U,U']:
                        #print(LW,combination)
                        #print(row)

                        a = row[3].replace('-','').replace('*',"'")
                        b = row[4].replace('*',"'")
                        c = row[6].replace('*',"'")
                        d = row[7].replace('-','').replace('*',"'")

                        if not combination in hbond:
                            hbond[combination] = {}

                        if not LW in hbond[combination]:
                            hbond[combination][LW] = []

                        # store as donor, hydrogen, acceptor, and order of nucleotides
                        if len(a) > 0:
                            hbond[combination][LW].append((a,b,c,'12'))
                        else:
                            hbond[combination][LW].append((d,c,b,'21'))

                    elif combination in ['C,A','G,A','U,A','C,G','U,C','U,G']:
                        # reverse order of edges, bases, and atoms
                        #print(LW,combination)
                        #print(row)

                        LW = LW[0]+LW[2]+LW[1]
                        combination = row[2]+','+row[1]

                        #print(LW,combination)

                        a = row[7].replace('-','').replace('*',"'")
                        b = row[6].replace('*',"'")
                        c = row[4].replace('*',"'")
                        d = row[3].replace('-','').replace('*',"'")

                        if not combination in hbond:
                            hbond[combination] = {}

                        if not LW in hbond[combination]:
                            hbond[combination][LW] = []

                        # store as donor, hydrogen, acceptor, and order of nucleotides
                        if len(a) > 0:
                            hbond[combination][LW].append((a,b,c,'12'))
                        else:
                            hbond[combination][LW].append((d,c,b,'21'))

                    # apply RNA hydrogen bonds to DNA interactions
                    # sloppy for now, but better as time permits
                    # will want to remove O2' hbonds from the DNA side

                    RNA_to_DNA = {'A': 'DA', 'C': 'DC', 'G': 'DG', 'U': 'DT'}

                    DNA_pair = RNA_to_DNA[b1] + "," + RNA_to_DNA[b2]
                    if not DNA_pair in hbond:
                        hbond[DNA_pair] = {}
                    hbond[DNA_pair][LW] = hbond[combination][LW]

                    DNA_pair = b1 + "," + RNA_to_DNA[b2]
                    if not DNA_pair in hbond:
                        hbond[DNA_pair] = {}
                    hbond[DNA_pair][LW] = hbond[combination][LW]

                    DNA_pair = RNA_to_DNA[b1] + "," + b2
                    if not DNA_pair in hbond:
                        hbond[DNA_pair] = {}
                    hbond[DNA_pair][LW] = hbond[combination][LW]






    """
    for combination in hbond:
        for LW in hbond[combination]:
            print(combination, LW, hbond[combination][LW])
    """

    return hbond


def calculate_hb_angle(A,B,C):
    if len(A) == 3 and len(B) == 3 and len(C) == 3:
        return angle_between_vectors(np.subtract(A,B),np.subtract(C,B))

# This function calculates an angle from 0 to 180 degrees between two vectors
def angle_between_vectors(vec1, vec2):
    if len(vec1) == 3 and len(vec2) == 3:
        cosang = np.dot(vec1, vec2)
        sinang = np.linalg.norm(np.cross(vec1, vec2))
        angle = np.arctan2(sinang, cosang)
        return 180*angle/np.pi
    else:
        return None

def check_hydrogen_bond(nt1,nt2,atoms):
    """
    Calculate hydrogen bond parameters for hydrogen donor and hydrogen from nt1
    and hydrogen bond acceptor from nt2.
    Return a list telling:
      Whether or not the bond could be checked
      Whether or not the bond meets criteria to form
      Distance between hydrogen and acceptor if available
      Angle in degrees if available
      A measure of the "badness" of the bond
    """

    seq = nt1.sequence
    if seq in ['A','C','G','U','DA','DC','DG','DT']:
        donor    = nt1.centers[atoms[0]]
        hydrogen = nt1.centers[atoms[1]]
    else:
        #print('Checking hydrogen bond for %s, mapping atoms' % nt1.unit_id())

        if seq in parent_atom_to_modified:
            # map the atom name
            if atoms[0] in parent_atom_to_modified[seq]:
                donor    = nt1.centers[parent_atom_to_modified[seq][atoms[0]]]
                #print('Parent atom %s modified atom %s length %d' % (atoms[0],parent_atom_to_modified[seq][atoms[0]],len(donor)))
            else:
                donor    = nt1.centers[atoms[0]]

            if atoms[1] in parent_atom_to_modified[seq]:
                hydrogen = nt1.centers[parent_atom_to_modified[seq][atoms[1]]]
                #print('Parent atom %s modified atom %s length %d' % (atoms[1],parent_atom_to_modified[seq][atoms[1]],len(hydrogen)))
            else:
                hydrogen = nt1.centers[atoms[1]]


        else:
            # uknown modified nucleotide, just try for the named atom
            print('%s is not a known modified nucleotide' % nt1.unit_id())
            donor    = nt1.centers[atoms[0]]
            hydrogen = nt1.centers[atoms[1]]


    seq = nt2.sequence
    if seq in ['A','C','G','U','DA','DC','DG','DT']:
        acceptor = nt2.centers[atoms[2]]
    else:
        #print('Checking hydrogen bond for %s, mapping atoms' % nt2.unit_id())

        if seq in parent_atom_to_modified:
            # map the atom name
            if atoms[2] in parent_atom_to_modified[seq]:
                acceptor = nt2.centers[parent_atom_to_modified[seq][atoms[2]]]
                #print('Parent atom %s modified atom %s length %d' % (atoms[2],parent_atom_to_modified[seq][atoms[2]],len(acceptor)))
            else:
                acceptor = nt2.centers[atoms[2]]

        else:
            # uknown modified nucleotide, just try for the named atom
            print('%s is not a known modified nucleotide' % nt2.unit_id())
            acceptor = nt2.centers[atoms[2]]


    if atoms[1] == "H2'" or len(hydrogen) < 3:
        # hydrogen coordinates are not available, check donor-acceptor distance
        if len(donor) == 3 and len(acceptor) == 3:
            heavy_distance = np.linalg.norm(np.subtract(donor,acceptor))
        else:
            return False, False, float("NaN"), float("NaN"), float("Inf")

        if heavy_distance > 4.5:
            return True, False, heavy_distance, float("NaN"), max(0,heavy_distance-3.0)
        else:
            return True, True, heavy_distance, float("NaN"), max(0,heavy_distance-3.0)


    else:
        if len(hydrogen) == 3 and len(acceptor) == 3:
            distance = np.linalg.norm(np.subtract(hydrogen,acceptor))
        else:
            return False, False, float("NaN"), float("NaN"), float("Inf")

        if len(donor) == 3:
            hb_angle = calculate_hb_angle(donor,hydrogen,acceptor)

            if not hb_angle:
                return False, False, distance, float("NaN"), float("Inf")
            elif hb_angle > 110 and distance < 4.0:
                return True, True, distance, hb_angle, (max(0,distance-2.5)+max(0,150-hb_angle)/20.0)
            else:
                return True, False, distance, hb_angle, (max(0,distance-2.5)+max(0,150-hb_angle)/20.0)

        elif distance < 4.0:
            return True, True, distance, float("NaN"), max(0,distance-2.5)
        else:
            return True, False, distance, float("NaN"), max(0,distance-2.5)


if __name__=="__main__":

    hbond = load_ideal_basepair_hydrogen_bonds()

    angle_sets = set()

    for bc in hbond.keys():
        for interaction in hbond[bc].keys():
            for atom_set in hbond[bc][interaction]:

                print(atom_set)
                a1,h1,a2,direction = atom_set

                base1, base2 = bc.split(",")

                if "D" in base1 or "D" in base2:
                    continue

                if direction == '12':
                    angle_set = "%s,,%s,%s" % (base1,a1,h1)
                else:
                    angle_set = "%s,,%s,%s" % (base2,a1,h1)

                angle_sets.add(angle_set)

    for angle_set in sorted(angle_sets):
        print(angle_set)
