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

    # identify the heavy base atom to use for calculating the heavy-donor-acceptor angle
    base_donor_to_atom_for_angle = {}
    base_donor_to_atom_for_angle["A"] = {}
    base_donor_to_atom_for_angle["A"]["C2"] = "N3"    # atoms N3-C2-H2
    base_donor_to_atom_for_angle["A"]["C8"] = "N7"    # atoms N7-C8-H8
    base_donor_to_atom_for_angle["A"]["N1"] = "C2"    # atoms C2-N1-H1
    base_donor_to_atom_for_angle["A"]["N6"] = "C6"    # atoms C6-N6-H61
    base_donor_to_atom_for_angle["A"]["O2'"] = "C2'"    # atoms C2'-O2'-H2'
    base_donor_to_atom_for_angle["C"] = {}
    base_donor_to_atom_for_angle["C"]["C5"] = "C4"    # atoms C4-C5-H5
    base_donor_to_atom_for_angle["C"]["C6"] = "C5"    # atoms C5-C6-H6
    base_donor_to_atom_for_angle["C"]["N3"] = "C2"    # atoms C2-N3-H3
    base_donor_to_atom_for_angle["C"]["N4"] = "C4"    # atoms C4-N4-H41
    base_donor_to_atom_for_angle["C"]["O2'"] = "C2'"    # atoms C2'-O2'-H2'
    base_donor_to_atom_for_angle["G"] = {}
    base_donor_to_atom_for_angle["G"]["C8"] = "N7"    # atoms N7-C8-H8
    base_donor_to_atom_for_angle["G"]["N1"] = "C2"    # atoms C2-N1-H1
    base_donor_to_atom_for_angle["G"]["N2"] = "C2"    # atoms C2-N2-H21
    base_donor_to_atom_for_angle["G"]["O2'"] = "C2'"    # atoms C2'-O2'-H2'
    base_donor_to_atom_for_angle["U"] = {}
    base_donor_to_atom_for_angle["U"]["C5"] = "C4"    # atoms C4-C5-H5
    base_donor_to_atom_for_angle["U"]["N3"] = "C2"    # atoms C2-N3-H3
    base_donor_to_atom_for_angle["U"]["C6"] = "C5"    # atoms C5-C6-H6
    base_donor_to_atom_for_angle["U"]["O2'"] = "C2'"    # atoms C2'-O2'-H2'
    base_donor_to_atom_for_angle["DA"] = {}
    base_donor_to_atom_for_angle["DA"]["C2"] = "N3"    # atoms N3-C2-H2
    base_donor_to_atom_for_angle["DA"]["C8"] = "N7"    # atoms N7-C8-H8
    base_donor_to_atom_for_angle["DA"]["N1"] = "C2"    # atoms C2-N1-H1
    base_donor_to_atom_for_angle["DA"]["N6"] = "C6"    # atoms C6-N6-H61
    base_donor_to_atom_for_angle["DA"]["O2'"] = "C2'"    # atoms C2'-O2'-H2'
    base_donor_to_atom_for_angle["DC"] = {}
    base_donor_to_atom_for_angle["DC"]["C5"] = "C4"    # atoms C4-C5-H5
    base_donor_to_atom_for_angle["DC"]["C6"] = "C5"    # atoms C5-C6-H6
    base_donor_to_atom_for_angle["DC"]["N3"] = "C2"    # atoms C2-N3-H3
    base_donor_to_atom_for_angle["DC"]["N4"] = "C4"    # atoms C4-N4-H41
    base_donor_to_atom_for_angle["DC"]["O2'"] = "C2'"    # atoms C2'-O2'-H2'
    base_donor_to_atom_for_angle["DG"] = {}
    base_donor_to_atom_for_angle["DG"]["C8"] = "N7"    # atoms N7-C8-H8
    base_donor_to_atom_for_angle["DG"]["N1"] = "C2"    # atoms C2-N1-H1
    base_donor_to_atom_for_angle["DG"]["N2"] = "C2"    # atoms C2-N2-H21
    base_donor_to_atom_for_angle["DG"]["O2'"] = "C2'"    # atoms C2'-O2'-H2'
    base_donor_to_atom_for_angle["DT"] = {}
    base_donor_to_atom_for_angle["DT"]["C5"] = "C4"    # atoms C4-C5-H5
    base_donor_to_atom_for_angle["DT"]["N3"] = "C2"    # atoms C2-N3-H3
    base_donor_to_atom_for_angle["DT"]["C6"] = "C5"    # atoms C5-C6-H6
    base_donor_to_atom_for_angle["DT"]["O2'"] = "C2'"    # atoms C2'-O2'-H2'

    # default return values
    result = {}
    result['bond_checked'] = False
    result['bond_made'] = False
    result['length'] = float("NaN")    # based on distance between hydrogen and acceptor, if available
    result['angle'] = float("NaN")
    result['badness'] = float("NaN")
    result['donor_acceptor_distance'] = float("NaN")
    result['donor_acceptor_atoms'] = ''
    result['heavy_donor_acceptor_angle'] = float("NaN")
    result['heavy_donor_acceptor_atoms'] = ''

    # record atom names according to the standard nucleotide, not modified
    donor_atom = atoms[0]
    acceptor_atom = atoms[2]

    # get atom coordinates and name of heavy atom for angle
    seq = nt1.sequence
    if seq in ['A','C','G','U','DA','DC','DG','DT']:
        donor    = nt1.centers[atoms[0]]
        hydrogen = nt1.centers[atoms[1]]
        atom_for_angle_name = base_donor_to_atom_for_angle[seq][atoms[0]]
        atom_for_angle = nt1.centers[base_donor_to_atom_for_angle[seq][atoms[0]]]  # for heavy-heavy-heavy angle calculation
    else:
        #print('Checking hydrogen bond for %s, mapping atoms' % nt1.unit_id())

        if seq in modified_base_to_parent:
            # map the atom name
            if atoms[0] in parent_atom_to_modified[seq]:
                donor    = nt1.centers[parent_atom_to_modified[seq][atoms[0]]]
                #print('Parent atom %s modified atom %s length %d' % (atoms[0],parent_atom_to_modified[seq][atoms[0]],len(donor)))
            else:
                # unmapped atom, just try the original atom name
                donor    = nt1.centers[atoms[0]]

            if atoms[1] in parent_atom_to_modified[seq]:
                hydrogen = nt1.centers[parent_atom_to_modified[seq][atoms[1]]]
                #print('Parent atom %s modified atom %s length %d' % (atoms[1],parent_atom_to_modified[seq][atoms[1]],len(hydrogen)))
            else:
                hydrogen = nt1.centers[atoms[1]]

            parent = modified_base_to_parent[seq]
            parent_atom_for_angle = base_donor_to_atom_for_angle[parent][atoms[0]]

            if parent_atom_for_angle in parent_atom_to_modified[seq]:
                atom_for_angle_name = parent_atom_to_modified[seq][parent_atom_for_angle]
            else:
                atom_for_angle_name = parent_atom_for_angle

            atom_for_angle = nt1.centers[atom_for_angle_name]

        else:
            # uknown modified nucleotide, just try for the named atom
            print('%s is not a known modified nucleotide' % nt1.unit_id())
            donor    = nt1.centers[atoms[0]]
            hydrogen = nt1.centers[atoms[1]]
            atom_for_angle = np.empty( shape=(0, 0) )    # no way to map it since we don't know the parent nucleotide

    # get acceptor atom coordinates
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

    # record names of heavy atoms from the standard base
    result["donor_acceptor_atoms"] = "%s-%s" % (donor_atom,acceptor_atom)
    result["heavy_donor_acceptor_atoms"] = "%s-%s-%s" % (atom_for_angle_name,donor_atom,acceptor_atom)

    # calculate donor-acceptor distance
    if len(donor) == 3 and len(acceptor) == 3:
        result["donor_acceptor_distance"] = np.linalg.norm(np.subtract(donor,acceptor))

    # calculate heavy-heavy-heavy angle
    if len(donor) == 3 and len(acceptor) == 3 and len(atom_for_angle) == 3:
        result["heavy_donor_acceptor_angle"] = calculate_hb_angle(atom_for_angle,donor,acceptor)

    # calculate distance and angle using hydrogen atom location, if available
    if atoms[1] == "H2'" or len(hydrogen) < 3:
        # hydrogen coordinates are not available, check donor-acceptor distance instead
        if len(donor) == 3 and len(acceptor) == 3:
            heavy_distance = np.linalg.norm(np.subtract(donor,acceptor))
            result["bond_checked"] = True
            result["distance"] = heavy_distance
            result["badness"] = max(0,heavy_distance-3.0)

            if heavy_distance < 4.5:
                result["bond_made"] = True

    else:
        if len(hydrogen) == 3 and len(acceptor) == 3:
            distance = np.linalg.norm(np.subtract(hydrogen,acceptor))
            result["bond_checked"] = True
            result["distance"] = distance
            result["badness"] = max(0,distance-2.5)

            if len(donor) == 3:
                hb_angle = calculate_hb_angle(donor,hydrogen,acceptor)

                if hb_angle:
                    result["badness"] = max(0,distance-2.5)+max(0,150-hb_angle)/20.0
                    result["angle"] = hb_angle

                    if hb_angle > 110 and distance < 4.0:
                        result["bond_made"] = True

            elif distance < 4.0:
                result["bond_made"] = True

    return result


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
