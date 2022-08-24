# This file holds methods for checking hydrogen bonds

import csv
import numpy as np

def load_ideal_basepair_hydrogen_bonds():
    """
    Load a table of ideal hydrogen bonds for each base combination
    in each Leontis-Westhof family, as determined by Jesse Stombaugh.
    For each base combination, store as a mapping from hydrogen bonds
    to basepairs that contain that bond.
    """

    hbond = {}

    with open('H_bonding_Atoms_from_Isostericity_Table.csv', newline='') as csvfile:
        bond_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        for row in bond_reader:
            if len(row) == 9:
                combination = row[1]+','+row[2]
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

                    hbond[combination][LW].append((a,b,c,d))
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

                    hbond[combination][LW].append((a,b,c,d))

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

def check_hbond_donor_acceptor(nt1,nt2,atoms):
    """
    Calculate hydrogen bond parameters for hydrogen donor and hydrogen from nt1
    and hydrogen bond acceptor from nt2.
    """

    donor    = nt1.centers[atoms[0]]
    hydrogen = nt1.centers[atoms[1]]
    acceptor = nt2.centers[atoms[2]]

    if atoms[1] == "H2'":
        if len(donor) == 3 and len(acceptor) == 3:
            distance = np.linalg.norm(np.subtract(donor,acceptor))
        else:
            return False, None, None

        if distance > 4.5:
            return False, distance, None
        else:
            return True, distance, None


    else:
        if len(hydrogen) == 3 and len(acceptor) == 3:
            distance = np.linalg.norm(np.subtract(hydrogen,acceptor))
        else:
            return False, None, None

        if distance > 4:
            return False, distance, None
        else:
            hb_angle = calculate_hb_angle(donor,hydrogen,acceptor)

        if not hb_angle:
            return False, distance, hb_angle
        elif hb_angle > 110:
            return True, distance, hb_angle
        else:
            return False, distance, hb_angle

def check_hydrogen_bond(nt1,nt2,atoms):
    """
    atoms is a list of three atoms in a 4-tuple
    The code first determines which direction is donor and acceptor.
    """

    if len(atoms[0]) > 0:
        hbond = check_hbond_donor_acceptor(nt1,nt2,atoms[0:3])
    else:
        hbond = check_hbond_donor_acceptor(nt2,nt1,atoms[3:0:-1])


    return hbond




