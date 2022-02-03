# -*- coding: utf-8 -*-
"""
    This program reads one or more CIF files and produces annotations
    of nucleotide-nucleotide interactions.
    Basepairs are annotated with Leontis-Westhof annotations like cWW, tHS.
    Basepairs might also be annotated as "near" with ncWW, ntHS.
    A few basepair categories have "alternative" geometries like acWW, ctWW.
    Alternative geometries are not checked for hydrogen bonds.
"""

"""
    Developer notes:
        CC acHS needs a tighter gap requirement, avoid http://rna.bgsu.edu/rna3dhub/display3D/unitid/4V9F|1|0|C|2309,4V9F|1|0|C|2281
        AC cWW needs tigher requirements, avoid http://rna.bgsu.edu/rna3dhub/display3D/unitid/4V9F|1|0|A|2465,4V9F|1|0|C|2396
        AG cWS has near that should be true, see big cluster at http://rna.bgsu.edu/webfr3d/Results/60ca1ea55ae61/60ca1ea55ae61.html

    When fr3d is changed, python setup.py install
"""

from fr3d.cif.reader import Cif
from fr3d.definitions import RNAconnections
from fr3d.definitions import NAbasecoordinates
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
from fr3d.modified_parent_mapping import modified_nucleotides

from discrepancy import matrix_discrepancy
import numpy as np
import csv
import urllib
import pickle
import math
import sys

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import defaultdict
# note that fr3d.localpath does not synchronize with Git, so you can change it locally to point to your own directory structure
from fr3d.localpath import outputText
from fr3d.localpath import outputNAPairwiseInteractions
from fr3d.localpath import outputNAPickleInteractions
from fr3d.localpath import contact_list_file
from fr3d.localpath import inputPath
from fr3d.localpath import outputHTML
from fr3d.data.base import EntitySelector

from fr3d.ordering.greedyInsertion import orderWithPathLengthFromDistanceMatrix

#from fr3d.classifiers.base_aafg import distance_metrics
from datetime import datetime
from math import floor
import os
from os import path

from time import time

from class_limits import nt_nt_cutoffs

nt_nt_screen_distance = 12

HB_donor_hydrogens = {}
HB_donor_hydrogens['A'] = {"N6":["1H6","2H6"], "C2":["H2"], "C8":["H8"], "O2'":[]}
HB_donor_hydrogens['G'] = {"N1":["H1"], "N2":["2H2","1H2"], "C8":["H8"], "O2'":[]}
HB_donor_hydrogens['C'] = {"N4":["1H4","2H4"], "C5":["H5"], "C6":["H6"], "O2'":[]}
HB_donor_hydrogens['U'] = {"N3":["H3"], "C5":["H5"], "C6":["H6"], "O2'":[]}


def load_basepair_annotations(filename,all_pair_types):

    with open(filename,'rb') as opener:
        basepairs = pickle.load(opener)

    pair_to_interaction={}

    for bp_type in all_pair_types:
      for (u1,u2,c) in basepairs[bp_type]:
         pair_to_interaction[(u1,u2)]=bp_type

      near_bp_type="n"+bp_type
#      for (u1,u2,c) in basepairs[near_bp_type]:
#         pair_to_interaction[(u1,u2)]=near_bp_type

    return pair_to_interaction


def draw_base(base_seq,dimensions,ax):
    """
    Connects atoms to draw one base in the specified number of dimensions
    """

    NAbasecolor = {}
    NAbasecolor['A'] = [1,0,0]   # red
    NAbasecolor['C'] = [1,214.0/255,0]   # yellowish
    NAbasecolor['G'] = [0.5,1,0]   # green
    NAbasecolor['U'] = [0,1,1]   # cyan
    NAbasecolor['DA'] = [1,0,0]   # red
    NAbasecolor['DC'] = [1,214.0/255,0]   # yellowish
    NAbasecolor['DG'] = [0.5,1,0]   # green
    NAbasecolor['DT'] = [0,0,1]   # blue

    new_base_x = []
    new_base_y = []
    new_base_z = []

    for atomname in RNAconnections[base_seq]:
        coord_base= NAbasecoordinates[base_seq][atomname]
        new_base_x.append(coord_base[0])
        new_base_y.append(coord_base[1])
        new_base_z.append(coord_base[2])

    #print(new_base_x)

    if dimensions == 2:
        ax.plot(new_base_x, new_base_y, color=NAbasecolor[base_seq], linewidth=2.0)
    else:
        ax.plot(new_base_x, new_base_y, new_base_z, color=NAbasecolor[base_seq], linewidth=2.0)

    return


#=======================================================================

base_seq_list = ['DA','DC','DG','DT']  # for DNA
base_seq_list = ['A','C','G','U']      # for RNA

interaction_lists = [["s3O2'","s3O3'","s3O4'","s3O5'","s3OP1","s3OP2","s5O2'","s5O3'","s5O4'","s5O5'","s5OP1","s5OP2"],
                    ["ns3O2'","ns3O3'","ns3O4'","ns3O5'","ns3OP1","ns3OP2","ns5O2'","ns5O3'","ns5O4'","ns5O5'","ns5OP1","ns5OP2"]]

# plot one instance of each of the pairwise interactions
PlotPair = False
PlotPair = True
AlreadyPlotted = {}

unit_data_path = "C:/Users/zirbel/Documents/GitHub/fr3d-python/data/units"

# Just do the annotation, no diagnostics or anything
test_for_pipeline = True

if __name__=="__main__":

    # load all datapoints of interactions after each file in case of a crash
    all_datapoints_output_file = outputNAPairwiseInteractions + "AllDatapoints.pickle"
    all_datapoints = pickle.load(open(all_datapoints_output_file,'rb'))
    print("Loaded %d datapoints from %s" % (len(all_datapoints),all_datapoints_output_file))

    # loop over nt1 bases, then over interactions of interest


    for interaction_list in interaction_lists:

        fig = plt.figure()
        for v, nt1_seq in enumerate(base_seq_list):

            if "n" in interaction_list[0]:
                ax = fig.add_subplot(2, 2, v+1)
            else:
                ax = fig.add_subplot(2, 2, v+1, projection='3d')

            ax.axis("equal")

            datapoints = []

            xvalues = []
            yvalues = []
            zvalues = []
            colors3d  = []
            colors2d  = []

            near_color = [1,0,0]  # red
            true_color = [0,0,0]  # black
            ring5_color = [0,1,1] # cyan
            ring6_color = [0,0,1] # blue

            c = 0
            for datapoint in all_datapoints:
#                if datapoint['nt1_seq'] == nt1_seq and datapoint['interaction'].replace("n","") in interaction_list:
                if datapoint['nt1_seq'] == nt1_seq and datapoint['interaction'] in interaction_list:
                    c += 1
                    if c <= 10000:
                        datapoints.append(datapoint)
                        xvalues.append(datapoint['x'])
                        yvalues.append(datapoint['y'])
                        zvalues.append(datapoint['z'])
                        if datapoint['interaction'].startswith("n"):
                            colors3d.append(near_color)
                        else:
                            colors3d.append(true_color)

                        if datapoint['ring'] == 'ring5':
                            colors2d.append(ring5_color)
                        elif datapoint['ring'] == 'ring6':
                            colors2d.append(ring6_color)
                        else:
                            colors2d.append(near_color)

#                    if datapoint['ring'] and (datapoint['y'] > 0.1 or datapoint['x'] > 0):
#                        print(datapoint)

            if "n" in interaction_list[0]:
                ax.scatter(xvalues,yvalues,color=colors2d,marker=".")
                ax.set_title('Base %s with %s near sO interactions' % (nt1_seq,len(datapoints)))
                draw_base(nt1_seq,2,ax)
            else:
                ax.scatter(xvalues,yvalues,zvalues,color=colors3d,marker=".")
                ax.scatter(xvalues,yvalues,[0 for i in zvalues],color=colors2d,marker=".")
                ax.set_title('Base %s with %s sO interactions' % (nt1_seq,len(datapoints)))
                draw_base(nt1_seq,3,ax)

            print("Plotted %d points" % len(datapoints))

    plt.show()
