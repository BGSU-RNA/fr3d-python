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

    with open('H_bonding_Atoms_from_Isostericity_Table.csv') as csvfile:
        bond_reader = csv.reader(csvfile, delimiter=',', quotechar='"')
        for row in bond_reader:
            if len(row) == 9:
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

    donor    = nt1.centers[atoms[0]]
    hydrogen = nt1.centers[atoms[1]]
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


