# -*- coding: utf-8 -*-
"""
    This program reads one or more CIF files and produces annotations
    the glycosidic bond orientation.

"""

import numpy as np
import argparse
import math
import os

from fr3d.modified_parent_mapping import modified_nucleotides

from fr3d.cif.reader import Cif
from fr3d.definitions import RNAconnections
from fr3d.definitions import NAbaseheavyatoms
from fr3d.definitions import NAbasehydrogens
from fr3d.definitions import nt_sugar
from fr3d.definitions import nt_phosphate
from fr3d.definitions import Ribophos_connect
from fr3d.definitions import aa_connections
from fr3d.definitions import aa_backconnect
from fr3d.definitions import aa_hydrogen_connections
from fr3d.definitions import aa_fg
from fr3d.definitions import aa_linker
from fr3d.definitions import aa_backbone
from fr3d.definitions import tilt_cutoff
from fr3d.definitions import planar_atoms
from fr3d.definitions import HB_donors
from fr3d.definitions import HB_weak_donors
from fr3d.definitions import HB_acceptors

# read input and output paths from localpath.py
# note that fr3d.localpath does not synchronize with Git, so you can change it locally to point to your own directory structure
try:
    from fr3d.localpath import outputNAPairwiseInteractions
    from fr3d.localpath import inputPath
except:
    inputPath = ""
    outputNAPairwiseInteractions = ""

from NA_pairwise_interactions import load_structure
from NA_pairwise_interactions import get_parent
from NA_pairwise_interactions import myTimer

def annotate_bond_orientation(structure,pdb,pipeline=False):

    bond_annotations = []
    error_message = []

    nts = structure.residues(type = ["RNA linking","DNA linking"])  # load all RNA/DNA nucleotides

    num_nts = 0
    for nt in nts:
        num_nts += 1

        chi = None
        classification = ""

        N1N9 = None
        C2C4 = None
        parent = None

        if nt.sequence in ['A','G','DA','DG']:
            N1N9 = nt.centers["N9"]
            C2C4 = nt.centers["C4"]
            parent = nt.sequence
        elif nt.sequence in ['C','U','DC','DT']:
            N1N9 = nt.centers["N1"]
            C2C4 = nt.centers["C2"]
            parent = nt.sequence
        else:
            parent = get_parent(nt.sequence)
            if parent in ['A','G','DA','DG']:
                N1N9 = nt.centers[modified_nucleotides[nt.sequence]["atoms"]["N9"]]
                C2C4 = nt.centers[modified_nucleotides[nt.sequence]["atoms"]["C4"]]
            elif parent in ['C','U','DC','DT']:
                N1N9 = nt.centers[modified_nucleotides[nt.sequence]["atoms"]["N1"]]
                C2C4 = nt.centers[modified_nucleotides[nt.sequence]["atoms"]["C2"]]
            else:
                if pipeline:
                    error_message.append("%s has no identified parent" % nt.unit_id())
                else:
                    print("%s has no identified parent" % nt.unit_id())

                N1N9 = nt.centers["N9"]      # maybe this is present
                if len(N1N9) == 3:
                    C2C4 = nt.centers["C4"]
                else:
                    N1N9 = nt.centers["N1"]
                    C2C4 = nt.centers["C2"]

        if len(N1N9) == 3 and len(C2C4) == 3:
            C1P = nt.centers["C1'"]
            O4P = nt.centers["O4'"]

            if len(C1P) == 3 and len(O4P) == 3:
                O4P_C1P = C1P - O4P
                N1N9_C1P  = C1P - N1N9
                C2C4_N1N9 = N1N9 - C2C4

                # chi angle definition:
                # O4*_C1*_N1_C2 (for pyrimidines)
                # O4*_C1*_N9_C4 (for purines):

                perp_to_sugar       = np.cross(O4P_C1P,N1N9_C1P)
                norm_perp_to_sugar  = np.linalg.norm(perp_to_sugar)
                if norm_perp_to_sugar != 0:
                    perp_to_sugar   = perp_to_sugar/norm_perp_to_sugar

                perp_to_base        = np.cross(-N1N9_C1P,C2C4_N1N9)
                norm_perp_to_base   = np.linalg.norm(perp_to_base)
                perp_to_base        = perp_to_base/norm_perp_to_base

                cross_cross_chi = np.cross(perp_to_base,perp_to_sugar)

                # Take the dot product of the vectors perp_to_base &
                # perp_to_sugar to get cos(chi).
                # Take norm(cross product) to get sin(chi).

                cos_chi = np.dot(perp_to_sugar,perp_to_base)
                if np.dot(cross_cross_chi,N1N9_C1P) > 0:
                    sin_chi = np.linalg.norm(cross_cross_chi)
                else:
                    sin_chi = -np.linalg.norm(cross_cross_chi)

                # sign of chi_degree matches Bevilacqua 2011 paper on syn and anti
                # this definition matches the IUPAC definition from http://www.chem.qmul.ac.uk/iupac/misc/pnuc2.html#230

                chi  = 180*math.atan2(sin_chi,cos_chi)/math.pi # glycosidic bond angle

                # Giving nomenclature according to chi values: anti (most common), or syn
                if chi > -90 and chi < -45:
                    classification = 'int_syn'
                elif chi >= -45 and chi < 90:
                    classification = 'syn'
                else:
                    classification = 'anti'

    #            print("%-20s has chi = %13.8f and classification %s" % (nt.unit_id(),chi,classification))

        if classification:
            bond_annotations.append({'unit_id'    : nt.unit_id(),
                                    'orientation' : classification,
                                    'chi_degree'  : ("%0.3f" % chi)})

        else:
            if pipeline:
                error_message.append('%s had a calculation error' % nt.unit_id())
            else:
                print('%s had a calculation error' % nt.unit_id())

            bond_annotations.append({'unit_id'    : nt.unit_id(),
                                    'orientation' : 'NA',
                                    'chi_degree'  : None})

    return bond_annotations, error_message

def write_txt_output_file(outputNAPairwiseInteractions,PDBid,bond_annotations):
    """
    Write interactions according to category, and within each
    category, write by annotation.
    Other than that, the interactions are listed in no particular order.
    """

    # loop over types of output files requested
    for category in categories.keys():
        filename = os.path.join(outputNAPairwiseInteractions,PDBid + "_" + category + ".txt")
        with open(filename,'w') as f:
            # loop over all interactions found in this category
            if category == 'glycosidic':
                for d in bond_annotations:
                    a = d['unit_id']
                    if d['orientation'] == 'NA':
                        f.write("%s\t%s\t%s\n" % (d['unit_id'],d['orientation'],d['chi_degree']))
                    else:
                        f.write("%s\t%s\t%0.8f\n" % (d['unit_id'],d['orientation'],d['chi_degree']))

#=======================================================================

ShowStructureReadingErrors = True


if __name__=="__main__":

    # allow user to specify input and output paths
    parser = argparse.ArgumentParser()
    parser.add_argument('PDBfiles', type=str, nargs='+', help='.cif filename(s)')
    parser.add_argument('-o', "--output", help="Output Location of Pairwise Interactions")
    parser.add_argument('-i', "--input", help='Input Path')
    parser.add_argument('-c', "--category", help='Interaction category or categories (glycosidic)')

    # process command line arguments
    args = parser.parse_args()
    if args.input:
        inputPath = args.input
    elif not inputPath:
        inputPath = ""
    if args.output:
        outputNAPairwiseInteractions = args.output     # set output path
    elif not outputNAPairwiseInteractions:
        outputNAPairwiseInteractions = ""

    # dictionary to control what specific annotations are output, in a file named for the key
    # empty list means to output all interactions in that category
    # non-empty list specifies which interactions to output in that category
    categories = {}

    if args.category:
        for category in args.category.split(","):
            categories[category] = []
    else:
        categories['glycosidic'] = []  # default

    # check existence of input path
    if len(inputPath) > 0 and not os.path.exists(inputPath):
        print("Attempting to create input path %s" % inputPath)
        os.mkdir(inputPath)

    # check existence of output path
    if len(outputNAPairwiseInteractions) > 0 and not os.path.exists(outputNAPairwiseInteractions):
        print("Attempting to create output path %s" % outputNAPairwiseInteractions)
        os.mkdir(outputNAPairwiseInteractions)

    # process additional arguments as PDB files
    PDBs = []
    entries = args.PDBfiles
    for entry in entries:
        if '.pdb' in entry.lower():
            x = entry
        else:
            x = entry.replace(".cif","") + ".cif"

        if "/" in entry or "\\" in entry:
            PDBs.append(x)
        else:
            PDBs.append(os.path.join(inputPath,x))

    # annotate each PDB file
    timerData = myTimer("start")
    failed_structures = []
    counter = 0

    for PDB in PDBs:
        counter += 1

        PDBid = PDB[-8:-4]

        print("Reading file " + PDB + ", which is number "+str(counter)+" out of "+str(len(PDBs)))
        timerData = myTimer("Reading CIF files",timerData)

        # suppress error messages, but report failures at the end
        try:
            structure = load_structure(PDB)
        except Exception as ex:
            print("  Could not load structure %s due to exception %s: %s" % (PDB,type(ex).__name__,ex))
            if type(ex).__name__ == "TypeError":
                print("  See suggestions in the fr3d-python Readme file")
            failed_structures.append((PDB,type(ex).__name__,ex))
            continue

        if 'glycosidic' in categories:
            timerData = myTimer("Annotating bond orientation",timerData)
            bond_annotations = annotate_bond_orientation(structure)

            #print(bond_annotations)

            timerData = myTimer("Recording interactions",timerData)
            print("  Recording interactions in %s" % outputNAPairwiseInteractions)
            write_txt_output_file(outputNAPairwiseInteractions,PDBid,bond_annotations)

    myTimer("summary",timerData)

    if len(failed_structures) > 0:
        print("Not able to read these files: %s" % failed_structures)
    else:
        print("All files read successfully")
