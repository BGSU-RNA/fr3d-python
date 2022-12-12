# -*- coding: utf-8 -*-
"""
    plot-phosphate-interactions.py reads a data file and plots points to represent
    base-base phosphate interactions.
"""

from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import pickle
import math
import sys 
import os
from collections import defaultdict
import urllib
import __builtin__

from fr3d.localpath import outputText
from fr3d.localpath import outputNAPairwiseInteractions
from fr3d.localpath import outputNAPickleInteractions
from fr3d.localpath import contact_list_file
from fr3d.localpath import inputPath
from fr3d.localpath import outputHTML
from fr3d.localpath import storeMatlabFR3DPairs
from class_limits import nt_nt_cutoffs

from NA_pairwise_interactions import map_PDB_list_to_PDB_IFE_dict
from draw_residues import draw_base


if sys.version_info[0] < 3:
    from urllib import urlopen
else:
    from urllib.request import urlopen 


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

#=======================================================================
def load_Matlab_FR3D_pairs(PDBID):
    """
    download annotations of RNA basepairs from http://rna.bgsu.edu/pairs/
    Those are triples of (unit_id1,unit_id2,crossingnumber)
    """

    interactionToTriples = defaultdict(list)

    pairsFileName = PDBID + '_RNA_pairs' + '.pickle'
    pathAndFileName = storeMatlabFR3DPairs + pairsFileName

    if not os.path.exists(pathAndFileName):
        print("Downloading %s to %s" % (pairsFileName,pathAndFileName))
        if sys.version_info[0] < 3:
            urllib.urlretrieve("http://rna.bgsu.edu/pairs/"+pairsFileName, pathAndFileName) # testing
        else:
            urllib.request.urlretrieve("http://rna.bgsu.edu/pairs/"+pairsFileName, pathAndFileName) # testing

        # note:  if the file is not present on the server, a text file with a 404 error will be downloaded
        # so it will look like a file was downloaded, but it's not the file you need!

    if os.path.exists(pathAndFileName):
        if sys.version_info[0] < 3:
            try:
                interactionToTriples = pickle.load(open(pathAndFileName,"rb"))
            except:
                print("Could not read "+pairsFileName+", it may not be available")
                try:
                    os.remove(pathAndFileName)
                    print("Removed unsuccessful file %s" % pathAndFileName)
                except:
                    print("Could not remove file %s" % pathAndFileName)
        else:
            try:
                interactionToTriples = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
            except:
                print("Could not read "+pairsFileName+", it may not be available")
                try:
                    os.remove(pathAndFileName)
                    print("Removed unsuccessful file %s" % pathAndFileName)
                except:
                    print("Could not remove file %" % pathAndFileName)

    interactionToPairs = {}
    for interaction in interactionToTriples.keys():
        interactionToPairs[interaction] = [(a,b) for a,b,c in interactionToTriples[interaction]]

    return interactionToPairs
#=======================================================================
def reverse(pair):
    return (pair[1],pair[0])
#=======================================================================
def check_full_data(datapoint):
    return ('gap12' in datapoint and 'angle_in_plane' in datapoint)
#=======================================================================
def print_datapoint(datapoint):

#    if 'url' in datapoint:
#        print("  %s" % datapoint['url'])
    fields = ['x','y','z','gap12','angle_in_plane','normal_Z']

    for f in fields:
        if f in datapoint:
            print("  %20s = %11.6f" % (f,datapoint[f]))

    return None
#=======================================================================
def plot_nt_nt_cutoffs(base_combination,lowercase_list,ax,variables):

    color = ["#BA55D3","#63B8FF","#00EE76","#FF8C00","#CDC9A5","#8A8A8A","#8A8A8A","#8A8A8A","#8A8A8A","#8A8A8A","#8A8A8A"]
    cc = 0

    for bc in nt_nt_cutoffs.keys():
        if bc == base_combination:
            for interaction in nt_nt_cutoffs[bc].keys():
                if interaction.lower().replace("a","") in lowercase_list:
                    for subcategory in nt_nt_cutoffs[bc][interaction].keys():
                        limits = nt_nt_cutoffs[bc][interaction][subcategory]
                        if variables == 1:            # x and y
                            xmin = limits['xmin']
                            xmax = limits['xmax']
                            ymin = limits['ymin']
                            ymax = limits['ymax']
                            ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            ax.text(0,cc*0.7,"%s %d" % (interaction,subcategory),color = color[cc],fontsize=12)
                            cc += 1
                        if variables == 2:            # theta and r
                            xmin = limits['xmin']
                            xmax = limits['xmax']
                            ymin = limits['ymin']
                            ymax = limits['ymax']
                            t1 = math.atan2(ymin,xmax)*180/3.141592654
                            t2 = math.atan2(ymax,xmax)*180/3.141592654
                            t3 = math.atan2(ymax,xmin)*180/3.141592654
                            t4 = math.atan2(ymin,xmin)*180/3.141592654
                            r1 = math.sqrt(xmax**2+ymin**2)
                            r2 = math.sqrt(xmax**2+ymax**2)
                            r3 = math.sqrt(xmin**2+ymax**2)
                            r4 = math.sqrt(xmin**2+ymin**2)

                            ax.plot([t1,t2,t3,t4,t1],[r1,r2,r3,r4,r1],color[cc])
                            ax.text(0,cc*0.7,"%s %d" % (interaction,subcategory),color = color[cc],fontsize=12)
                            cc += 1
                        if variables == 3:            # angle and normal
                            xmin = limits['anglemin']
                            xmax = limits['anglemax']
                            ymin = max(-1.005,limits['normalmin'])  # clamp
                            if ymin < -1:
                                ymin -= cc*0.001                   # don't overlap
                            ymax = min( 1.005,limits['normalmax'])  # clamp
                            if ymax > 1:
                                ymax += cc*0.001                   # don't overlap

                            if xmin < xmax:
                                # angle does not wrap around 270 degrees
                                ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            else:
                                # angle wraps around 270 degrees
                                ax.plot([270,xmin,xmin,270],[ymin,ymin,ymax,ymax],color[cc])
                                ax.plot([-90,xmax,xmax,-90],[ymin,ymin,ymax,ymax],color[cc])

                            cc += 1
                        if variables == 4:            # gap and z
                            xmin = -0.01
                            xmax = limits['gapmax']
                            ymin = limits['zmin'] - cc*0.04
                            ymax = limits['zmax'] + cc*0.04
                            ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            cc += 1
#=======================================================================
def plot_confusion_matrix(confusionMatrix, interaction_list):
    "Tabular construction of confusion matrix created from a dictionary of dictionaries of values."
    row = [['0BPh'],['1BPh'],['2BPh'],['3BPh'],['4BPh'],['5BPh'],['6BPh'],['7BPh'],['8BPh'],['9BPh'],['n0BPh'],['n1BPh'],['n2BPh'],['n3BPh'],['n4BPh'],['n5BPh'],['n6BPh'],['n7BPh'],['n8BPh'],['n9BPh'],['blank']]
    border = ['----','----','----','----','----','----','----','----','----','----','----','----']
    col = 0 
    interaction_list.append('blank')
    for interaction in interaction_list:
        for interaction2 in interaction_list:
                row[col].append(confusionMatrix[interaction][interaction2])
        col+=1

    print("\n\nConfusion Matrix of Annotations. \nColumns Represent Matlab found annotations and Rows Represent Python")
    print_function = getattr(__builtin__, 'print')
    print_function("\t", end = "")
    print_function(*interaction_list, sep='\t')
    print_function(*border, sep = "\t")
    for rows in row:
        print_function(*rows, sep='\t')

def plot_unmatched_pairs(interactionDict, interaction_list, ordering):
    """Creates a 2x9 plot to show which annotations were found by passed in language and not found by other passed in language"""
    if ordering == "python":
        other = "matlab"
    elif ordering == "matlab":
        other = "python"
    if 'blank' in interaction_list:
        interaction_list.remove('blank')

    interaction_list = ["0BPh", "1BPh","2BPh","3BPh","4BPh","5BPh","6BPh","7BPh","8BPh","9BPh",'total']
    print_function = getattr(__builtin__, 'print')
    print("\n\nAnnotations that were found by %s and not by %s" % (ordering, other))
    border = ['----','----','----','----','----','----','----','----', '----','----','----']
    print_function(*interaction_list, sep='\t')
    print_function(*border, sep = "\t")
    for category in interaction_list:
        print_function(interactionDict[category], end = "\t")
    print("\n")

def check_for_matching_pairs(Matlab_pairs, pair_to_data, not_loaded):
    """Loop through all Matlab pairs and all python pairs and see if any aren't found"""
    print("Checking For Unmatched Pairs Between Matlab and Python...")
    found = False
    notFound = []
    for pair in Matlab_pairs:
        found = False
        if(pair[0][0:4] not in not_loaded):
            for pair2, datapoint in pair_to_data.items():
                if pair == pair2:
                    found = True
                    break
            if not found:
                print(str(pair) + " Not Found by python but found by Matlab")
                notFound.append(pair)
    print("Amount of pairs found by Matlab but not by Python: " + str(len(notFound)))

if __name__=="__main__":

    base_seq_list = ['DA','DC','DG','DT']  # for DNA
    base_seq_list = ['A','C','G','U']      # for RNA

    interaction_lists = [["0BPh", "1BPh","2BPh","3BPh","4BPh","5BPh","6BPh","7BPh","8BPh","9BPh"]]

    interaction_lists = [["0BPh", "1BPh","2BPh","3BPh","4BPh","5BPh","6BPh","7BPh","8BPh","9BPh"],
                        ["n0BPh", "n1BPh","n2BPh","n3BPh","n4BPh","n5BPh","n6BPh","n7BPh","n8BPh","n9BPh"]]

    base_combination_list = ['A,A','A,C','A,G','A,U','C,C','G,C','C,U','G,G','G,U','U,U']
    symmetric_base_combination_list = ['A,A','C,C','G,G','U,U']
    
    ###############################################################################################################
    # Data collection for confusion matrix ########################################################################
    # each main key represents annotations found in python. The Keys within represent annotations found in matlab #
    # This dictionary is used to detect instances where python and matlab agree vs where they disagree ############
    ###############################################################################################################
    confusionMatrix = {} 
    confusionMatrix['0BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['1BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['2BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['3BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['4BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['5BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['6BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['7BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['8BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['9BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['n0BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['n1BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['n2BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['n3BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['n4BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['n5BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['n6BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['n7BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['n8BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['n9BPh']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['blank']={'0BPh':0,'1BPh':0,'2BPh':0,'3BPh':0,'4BPh':0,'5BPh':0,'6BPh':0,'7BPh':0,'8BPh':0, '9BPh':0,'n0BPh':0,'n1BPh':0,'n2BPh':0,'n3BPh':0,'n4BPh':0,'n5BPh':0,'n6BPh':0,'n7BPh':0, 'n8BPh':0,'n9BPh':0, 'blank':0}
    confusionMatrix['total'] = 0

    ###############################################################################################################
    # Collection to display totals ################################################################################
    ###############################################################################################################
    incorrectFaces = []
    nearVSTrue = [] 
    matlabNoPythonMatch = []
    pythonNoMatlabMatch = []
    ###############################################################################################################

    # plot one instance of each of the pairwise interactions
    PlotPair = False
    PlotPair = True
    AlreadyPlotted = {}

    unit_data_path = "C:/Users/jimitch/Documents/GitHub/fr3d-python/data/units"

    near_color = [1,0,0]  # red
    true_color = [0,0,0]  # black
    red   = [1,0,0] # red
    black = [0,0,0] # black
    cyan  = [0,1,1] # cyan
    blue  = [0,0,1] # blue

    plot_true_in_3D = False

    PDB_list = ['4V9F','7K00']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.220/1.5A/csv']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.220/2.0A/csv']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.220/2.5A/csv']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.217/3.0A/csv']
    # PDB_list = ['4V9F','6ZMI','7K00']
    PDB_list = ['6ZMI']
#    PDB_list = ['4V9F']
    # PDB_list = ['4TNA']
    PDB_list = ['4V9F']
    PDB_list = ['7K00']
    PDB_list = ['6XRQ']
    PDB_list = ['4V9F']
    PDB_list = ['4MCF']
    PDB_list = ['2C4Q']
    PDB_list = ['7K00']
    #PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.237/2.5A/csv']
    PDB_list = ['4V9F']

    #PDB LIST and Skip Files##################################################
    PDB_IFE_Dict = map_PDB_list_to_PDB_IFE_dict(PDB_list)

    all_PDB_ids = sorted(PDB_IFE_Dict.keys())
    print("Loading NA-pairwise-interactions from %d PDB files" % len(all_PDB_ids))
    PDB_skip_set = set(['1R9F','5NXT','4KTG'])

    pair_to_data = defaultdict(dict)\
    
    #Loading MatLAB Annotations for phosphate##################################
    print('Loading Matlab annotations of these files')
    Matlab_annotation_to_pair = defaultdict(list)
    pair_to_Matlab_annotation = defaultdict(str)
    PDB_IFE_Dict = map_PDB_list_to_PDB_IFE_dict(PDB_list)
    all_PDB_ids = PDB_IFE_Dict.keys()

    for PDB_id in all_PDB_ids:
        interaction_to_pairs = load_Matlab_FR3D_pairs(PDB_id)
        num_pairs = 0
        for interaction in interaction_to_pairs.keys():
            # store all basepairs and only basepairs
            if "BPh" in interaction:
                num_pairs += 1
                Matlab_annotation_to_pair[interaction] += interaction_to_pairs[interaction]
                for pair in interaction_to_pairs[interaction]:
                    pair_to_Matlab_annotation[pair] = interaction
        if num_pairs == 0:
            print("No Matlab-annotated pairs in %s" % PDB_id)

    print("Skipping %d PDB files because they have no Matlab annotation to compare to" % len(PDB_skip_set))
    print("Found these Matlab annotations: %s" % sorted(Matlab_annotation_to_pair.keys()))


    #Loading FR3D Python Annotations########################################
    
    all_PDB_ids = list(set(all_PDB_ids) - PDB_skip_set)
    print("Loading Python annotations from %d PDB files" % len(all_PDB_ids))

    pair_to_data = defaultdict(dict)

    not_loaded = []
        # load output files from NA_pairwise_interactions
    for PDB in all_PDB_ids:
        pair_to_data_file = outputNAPairwiseInteractions + "%s_pairs_v1.pickle" % PDB
        print("Reading %s" % pair_to_data_file)
        try:
            new_dict = pickle.load(open(pair_to_data_file,'rb'))
            pair_to_data.update(new_dict)
        except:
            print("Not able to load annotations for %s" % PDB)
            not_loaded.append(PDB)
            # Not very Efficient way but checking to see if pairs are listed in both python and matlab
            # For whatever reason, this isn't being found accurately where the other processing is.


    for interaction_list in interaction_lists:
        lowercase_list = [i.lower() for i in interaction_list]

        # identify unit id pairs that are annotated as basepairing by Matlab code
        Matlab_pairs = []
        for interaction in interaction_list:
            Matlab_pairs += Matlab_annotation_to_pair[interaction]
        Matlab_pairs = set(Matlab_pairs)

        if False:
            check_for_matching_pairs(Matlab_pairs, pair_to_data, not_loaded)

        for base_combination in base_combination_list:
            nt1_seq, nt2_seq = base_combination.split(",")
         
            # accumulate data specific to this interaction and base combination
            xvalues = []
            yvalues = []
            rvalues = []
            tvalues = []
            zvalues = []
            avalues = []
            gvalues = []
            nvalues = []
            colors2d  = []
            sizes = []

            c = 0

            for pair,datapoint in pair_to_data.items():
                # restrict to the current base combination
                if not datapoint['nt1_seq'] == nt1_seq:
                    continue

                if not datapoint['nt2_seq'] == nt2_seq:
                    continue

                # eliminate Python annotations of alternate locations other than A
                fields1 = pair[0].split("|")

                # some PDB ids are not annotated by Matlab, just skip them here
                if fields1[0] in PDB_skip_set:
                    continue

                if len(fields1) > 5 and not fields1[5] == 'A':
                    continue

                fields2 = pair[1].split("|")
                if len(fields2) > 5 and not fields2[5] == 'A':
                    continue
                
                standardBases = ["A", "C", "G", "U"]
                fields3 = pair[0].split("|")

                if fields3[3] not in standardBases:
                    print(fields3)
                    sys.exit(0)
                have_data_in_other_order = False

                if 'BPh' in datapoint and datapoint['BPh'].lower() in lowercase_list:
                    Python = True
                else:
                    Python = False

                if pair in Matlab_pairs:
                    Matlab = True
                else:
                    Matlab = False

                have_full_data = check_full_data(datapoint)

                if reverse(pair) in pair_to_data:
                    r_datapoint = pair_to_data[reverse(pair)]

                if Matlab and not have_full_data:
                    if reverse(pair) in pair_to_data:
                        r_datapoint = pair_to_data[reverse(pair)]

                        if check_full_data(r_datapoint):
                            have_data_in_other_order = True
                        else:
                            print('Matlab annotation, data in both orders, but no match %s' % datapoint['url'])
                            print(pair)
                            print(reverse(pair))
                    else:
                        print('Matlab annotation but missing data in both orders for %s' % datapoint['url'])
                        print(pair)
                        print(datapoint)
                        print(reverse(pair))

                if (Python or Matlab) and have_full_data:
                    c += 1
                    xvalues.append(datapoint['x'])
                    yvalues.append(datapoint['y'])
                    rvalues.append(math.sqrt(datapoint['x']**2 + datapoint['y']**2))
                    tvalues.append(math.atan2(datapoint['y'],datapoint['x'])*180/3.141592654)
                    zvalues.append(datapoint['z'])
                    gvalues.append(datapoint['gap12'])
                    avalues.append(datapoint['angle_in_plane'])
                    nvalues.append(datapoint['normal_Z'])

                    if 'BPh' in datapoint:
                        phosphate = datapoint['BPh']
                    else:
                        phosphate = "   "
                    print(phosphate)
                    print(pair_to_Matlab_annotation[pair])
                    ###########################################################################################
                    # Add record of annotations by pair #######################################################
                    # confusionMatrix is a dictionary - key is python annotation, value is dictionary where ###
                    # key is matlab annotation - values are numerical total ###################################
                    ###########################################################################################
                    if phosphate != "   " and pair_to_Matlab_annotation[pair] != '':
                        # Python finds an annotation and Matlab also finds an annotation ################
                        confusionMatrix[phosphate][pair_to_Matlab_annotation[pair]] += 1
                        confusionMatrix['total'] += 1

                    elif phosphate != "   " and pair_to_Matlab_annotation[pair] == '':
                        # Python finds an annotation and matlab does not ######################################
                        confusionMatrix[phosphate]['blank'] += 1
                        confusionMatrix['total'] += 1    
                    elif phosphate == "   " and Matlab_annotation_to_pair != "":  
                        # Python does not find an annotation and matlab does ##################################
                        confusionMatrix['blank'][pair_to_Matlab_annotation[pair]] += 1
                        confusionMatrix['total'] += 1

                    ###########################################################################################
                    # Processing for construction of graphs ###################################################
                    ###########################################################################################
                    if Python and Matlab:
                        color = black
                        size = 1

                    if Python and not Matlab:
                        color = blue
                        size = 40
                        print("blue = only Python:  %s Python %4s Matlab %4s %s" % (base_combination,phosphate,pair_to_Matlab_annotation[pair],datapoint['url']))
                        print_datapoint(datapoint)

                    if not Python and Matlab:
                        color = red
                        size = 40
                        print("red = only Matlab: %s Python %4s Matlab %4s %s" % (base_combination,phosphate,pair_to_Matlab_annotation[pair],datapoint['url']))
                        print_datapoint(datapoint)
                    else:
                        color = [0,1,0]
                        size = 100
                    colors2d.append(color)
                    sizes.append(size)

            print("Plotted %5d points for %s %s" % (c,base_combination,",".join(interaction_list)))

            if c > 0:
                fig = plt.figure(figsize=(11.0, 6.0))

                ax = fig.add_subplot(1, 4, 1)
                ax.axis("equal")
                plot_nt_nt_cutoffs(base_combination,lowercase_list,ax,1)
                ax.scatter(xvalues,yvalues,color=colors2d,marker=".",s=sizes)
                ax.set_title('x and y for %d %s %s' % (c,base_combination,interaction_list[0]), rotation=10)

                draw_base(nt1_seq,'default',2,ax)

                ax = fig.add_subplot(1, 4, 2)
                plot_nt_nt_cutoffs(base_combination,lowercase_list,ax,2)
                ax.scatter(tvalues,rvalues,color=colors2d,marker=".",s=sizes)
                ax.set_title('theta and radius', rotation = 10)

                ax = fig.add_subplot(1, 4, 3)
                plot_nt_nt_cutoffs(base_combination,lowercase_list,ax,3)
                ax.scatter(avalues,nvalues,color=colors2d,marker=".",s=sizes)
                ax.set_title('angle and normal, blue=Matlab only',rotation = 10)

                ax = fig.add_subplot(1, 4, 4)
                plot_nt_nt_cutoffs(base_combination,lowercase_list,ax,4)
                ax.scatter(gvalues,zvalues,color=colors2d,marker=".",s=sizes)
                ax.set_title('gap12 and z values, red=Python only',rotation = 10)

                # show all plots for this interaction_list
                figManager = plt.get_current_fig_manager()
                figManager.full_screen_toggle()
                figure_save_file = outputNAPairwiseInteractions + "phosphate_%s_%s_%s.png" % (base_combination,interaction_list[0],len(all_PDB_ids))
                plt.savefig(figure_save_file)
                #plt.show()
                plt.close()

    # # loop over sets of interactions, plotting points on bases
    # for interaction_list in interaction_lists:

    #     fig = plt.figure(figsize=(10.0,8.0))

    #     # loop over individual bases, to make separate plots for each base
    #     for v, nt1_seq in enumerate(base_seq_list):

    #         if "n" in interaction_list[0]:
    #             ax = fig.add_subplot(2, 2, v+1)
    #         elif plot_true_in_3D:
    #             ax = fig.add_subplot(2, 2, v+1, projection='3d')
    #         else:
    #             ax = fig.add_subplot(2, 2, v+1)

    #         ax.axis("equal")

    #         # accumulate data specific to this base
    #         xvalues = []
    #         yvalues = []
    #         zvalues = []
    #         colors3d  = []
    #         colors2d  = []
    #         colorsoxy = []

    #         c = 0

    #         # loop over pairs, finding those with the desired base and interaction
    #         for pair,datapoint in pair_to_data.items():
    #             if not datapoint['nt1_seq'] == nt1_seq:
    #                 continue

    #             if not 'BPh' in datapoint:
    #                 continue

    #             if datapoint['BPh'] in interaction_list:


    #                 inside = False

    #                 #if datapoint['nt1_seq'] == 'G':
    #                 #    print("%3d %5s %s" % (c,datapoint['sOinteraction'],datapoint['url']))

    #                 parent1 = nt1_seq
    #                 x = datapoint['xStack']
    #                 y = datapoint['yStack']

    #                 if parent1 == 'A' or parent1 == 'DA':
    #                     if -2.327244*x +  4.271447*y +  9.515028 > 0:  # Left of H9'-H2
    #                         if -3.832809*x + -2.350503*y + 10.927316 > 0:  # Left of H2-H61
    #                             if  0.451014*x + -1.690509*y +  5.259508 > 0:  # Left of H61-H62
    #                                 if  4.252574*x + -2.330898*y + 10.447200 > 0:  # Left of H62-H8
    #                                     if  1.456465*x +  2.100463*y +  7.567280 > 0:  # Left of H8-H9'
    #                                         inside = True
    #                 elif parent1 == 'C' or parent1 == 'DC':
    #                     if -0.889476*x +  2.269450*y +  5.323403 > 0:  # Left of H1'-O2
    #                         if -4.532779*x + -1.065616*y +  6.851131 > 0:  # Left of O2-H42
    #                             if -0.190206*x + -1.731804*y +  5.226294 > 0:  # Left of H42-H41
    #                                 if  1.955107*x + -1.508802*y +  6.480180 > 0:  # Left of H41-H5
    #                                     if  2.523463*x + -0.045961*y +  6.153526 > 0:  # Left of H5-H6
    #                                         if  1.133891*x +  2.082733*y +  5.627625 > 0:  # Left of H6-H1'
    #                                             inside = True
    #                 elif parent1 == 'G' or parent1 == 'DG':
    #                     if -1.107310*x +  4.872516*y + 11.647152 > 0:  # Left of H9'-H21
    #                         if -1.684502*x +  0.422659*y +  6.436199 > 0:  # Left of H21-H22
    #                             if -1.592264*x + -1.681840*y +  6.230291 > 0:  # Left of H22-H1
    #                                 if -1.019666*x + -2.216349*y +  5.884100 > 0:  # Left of H1-O6
    #                                     if  2.274081*x + -2.148378*y +  5.898397 > 0:  # Left of O6-N7
    #                                         if  1.656548*x + -1.350181*y +  4.208981 > 0:  # Left of N7-H8
    #                                             if  1.473113*x +  2.101573*y +  7.865111 > 0:  # Left of H8-H9'
    #                                                 inside = True
                    
    #                 elif parent1 == 'U':
    #                     if -0.960553*x +  2.292490*y +  5.471254 > 0:  # Left of H1'-O2
    #                         if -2.493573*x + -0.200338*y +  4.589448 > 0:  # Left of O2-H3
    #                             if -1.574881*x + -1.914996*y +  4.563214 > 0:  # Left of H3-O4
    #                                 if  1.403523*x + -2.301733*y +  5.976805 > 0:  # Left of O4-H5
    #                                     if  2.504701*x +  0.041797*y +  6.092950 > 0:  # Left of H5-H6
    #                                         if  1.120783*x +  2.082780*y +  5.621468 > 0:  # Left of H6-H1' 
    #                                             inside = True
    #                 elif parent1 == 'DT':
    #                     if -1.125648*x +  2.281277*y +  6.199955 > 0:  # Left of C1'-O2
    #                         if -2.368105*x + -0.456021*y +  4.878252 > 0:  # Left of O2-H3
    #                             if -1.526233*x + -1.897795*y +  4.450270 > 0:  # Left of H3-O4
    #                                 if  1.301401*x + -2.544887*y +  5.949759 > 0:  # Left of O4-C7
    #                                     if  2.031505*x +  1.412190*y +  3.691439 > 0:  # Left of C7-C6
    #                                         if  1.687080*x +  1.205236*y +  3.097805 > 0:  # Left of C6-C1'
    #                                             inside = True

    #                 if inside:     
    #                     xvalues.append(datapoint['xStack'])
    #                     yvalues.append(datapoint['yStack'])
    #                     zvalues.append(datapoint['zStack'])
    #                     if datapoint['BPh'].startswith("n"):
    #                         colors2d.append(near_color)
    #                     else:
    #                         colors2d.append(true_color)

    #         c = len(xvalues)
    #         if len(xvalues) < 500:
    #             markersize = 4
    #         else:
    #             markersize = 0.5

    #         if "n" in interaction_list[0]:
    #             # 2D plot for near interactions
    #             draw_base(nt1_seq,'CPK',2,ax,zorder=1)
    #             ax.scatter(xvalues,yvalues,color=colors2d,marker=".",s=markersize,zorder=3)
    #             ax.set_title('Base %s with %s near phosphate interactions' % (nt1_seq,c))
    #         elif plot_true_in_3D:
    #             # 3D plot for true interactions
    #             draw_base(nt1_seq,'CPK',3,ax)

    #             ax.scatter(xvalues,yvalues,zvalues,color=colors3d,marker=".",s=markersize)
    #             ax.scatter(xvalues,yvalues,[0 for i in zvalues],color=colors2d,marker=".",s=markersize)
    #             ax.set_title('Base %s with %s phosphate interactions' % (nt1_seq,c))
    #         else:
    #             # 2D plot for true interactions
    #             draw_base(nt1_seq,'CPK',2,ax,zorder=1)
    #             ax.scatter(xvalues,yvalues,color=colors2d,marker=".",s=markersize,zorder=3)
    #             ax.set_title('Base %s with %s phosphate interactions' % (nt1_seq,c))

    #         print("Plotted %d points for %s" % (c,nt1_seq))


    #     figManager = plt.get_current_fig_manager()
    #     figManager.full_screen_toggle()

    #     if "n" in interaction_list[0]:
    #         figure_save_file = outputNAPairwiseInteractions + "phosphate_near_%d" % len(all_PDB_ids)
    #     else:
    #         figure_save_file = outputNAPairwiseInteractions + "phosphate_true_%d" % len(all_PDB_ids)
    #     plt.savefig(figure_save_file+".png")
    #     plt.savefig(figure_save_file+".pdf")
    #     plt.close()


    matlabDictionary = {} 
    # loop over sets of interactions, plotting histograms of z values
    interaction_list = ["0BPh", "1BPh","2BPh","3BPh","4BPh","5BPh","6BPh","7BPh","8BPh","9BPh","n0BPh", "n1BPh","n2BPh","n3BPh","n4BPh","n5BPh","n6BPh","n7BPh","n8BPh","n9BPh"]

    # fig = plt.figure(figsize=(10.0,8.0))



    # # loop over individual bases, to make separate plots for each base
    # for v, nt1_seq in enumerate(base_seq_list):
    #     ax = fig.add_subplot(2, 2, v+1)

    #     # accumulate data specific to this base
    #     zvalues = []
    #     # loop over pairs, finding those with the desired base and interaction
    #     for pair,datapoint in pair_to_data.items():
    #         if not datapoint['nt1_seq'] == nt1_seq:
    #             continue

    #         if not 'BPh' in datapoint:
    #             continue

    #         if abs(datapoint['zStack']) < 2:
    #             continue
    #         if datapoint['BPh'] in interaction_list:
    #             if not "near" in datapoint['BPh']:
    #                 zvalues.append(abs(datapoint['zStack']))

    #     # plot histogram of z values for this base
    #     bins = [2.7+0.05*k for k in range(0,21)]
    #     plt.hist(zvalues,bins=bins)
    #     ax.set_title('Base %s with %s true and near phosphate interactions' % (nt1_seq,len(zvalues)))

    #     print("Plotted %d points for %s" % (len(zvalues),nt1_seq))

    # figManager = plt.get_current_fig_manager()
    # figManager.full_screen_toggle()
    # figure_save_file = outputNAPairwiseInteractions + "phosphate_z_histogram_%d" % len(all_PDB_ids)
    # plt.savefig(figure_save_file+".png")
    # plt.savefig(figure_save_file+".pdf")
    # plt.close()

    # print("Plots are in %s" % outputNAPairwiseInteractions){0:<10}




    ### Construction of confusion matrices ##############################################################################################
    plot_confusion_matrix(confusionMatrix, interaction_list)

    print("Unit ids of pairs that disagree on faces being used. \n ('unit_id_1', 'unit_id_2','PythonAnnotation,MatLabAnnotation) ")
    

    ## Create files in pwd to output pairs stored annotations ##########################################################################
    original_stdout = sys.stdout # Save a reference to the original standard output
    if True: 
        with open('confusionMatrix_phosphate.txt', 'w') as f:
            sys.stdout = f # Change the standard output to the file we created.
            print('Confusion Matrix for Python VS Matlab annotations of base phosphate interactions.')
            sys.stdout.write("{0:<5}".format("\t"))
            for category in interaction_list:
                sys.stdout.write("{0:<5}".format(category))
                sys.stdout.write("{0:<5}".format("\t"))
            sys.stdout.write("\n")
            for row in interaction_list:
                if row != 'total':
                    sys.stdout.write("{0:<5}".format(row))
                    sys.stdout.write("{0:<5}".format("\t"))
                for col in interaction_list:
                    if row != 'total' and col != 'total':
                        if(confusionMatrix[row][col] != 0):
                            sys.stdout.write("{0:<5}".format(str(confusionMatrix[row][col])))
                            sys.stdout.write("{0:<5}".format("\t"))
                        else:
                            sys.stdout.write("{0:<5}".format("0"))
                            sys.stdout.write("{0:<5}".format("\t"))                            
                sys.stdout.write("\n")
            for pair in incorrectFaces:
                print(pair)
            sys.stdout = original_stdout

    # if True: 
    #     with open('nearVsTrue_phosphate.txt', 'w') as f:
    #         sys.stdout = f # Change the standard output to the file we created.
    #         print('Annotations where Python and Matlab phosphate annotations disagree on near and true phosphate')
    #         for pair in nearVSTrue:
    #             print(pair)
    #         sys.stdout = original_stdout
    # if True: 
    #     with open('PythonNoMatlabMatch_phosphate.txt', 'w') as f:
    #         sys.stdout = f # Change the standard output to the file we created.
    #         print('Annotations where Python finds an interaction but Matlab doesnt')
    #         for pair in pythonNoMatlabMatch:
    #             print(pair)
    #         sys.stdout = original_stdout

    # if True: 
    #     with open('MatlabNoPythonMatch_phosphate.txt', 'w') as f:
    #         sys.stdout = f # Change the standard output to the file we created.
    #         print('Annotations where Matlab finds an interaction but Python doesnt')
    #         for pair in matlabNoPythonMatch:
    #             print(pair)
    #         sys.stdout = original_stdout

    #####################################################################################################################################
    
    if False: #Output Near/Total and matching/total percents
        percents = {}
        percents['overallML'] = 0
        percents['overallPy'] = 0 

        pythonTotals['near']=0
        pythonTotals['true']=0
        matlabTotals['near']=0
        matlabTotals['true']=0
        for cat in interaction_list:
            if cat[0] == "n":
                pythonTotals['near'] += pythonTotals[cat]
                matlabTotals['near'] += matlabTotals[cat]
            else: 
                matlabTotals['true'] += matlabTotals[cat]
                pythonTotals['true'] += pythonTotals[cat]
            pythonTotals["total"] += pythonTotals[cat]
            matlabTotals["total"] += matlabTotals[cat]
        print(pythonTotals)
        print(matlabTotals)
        test = pythonTotals['near']/pythonTotals['total']
        print(pythonTotals['total'])

        if pythonTotals > matlabTotals:
            larger = "PYTHON"
        else:
            larger = "MATLAB"
        trueInts=["s33",'s55','s53','s35']
        nearInts=['ns33','ns55','ns53','ns35']
        print("THE NUMBER DIFFERENCE BETWEEN PYTHON AND MATLAB IS: " + str(abs(pythonTotals['total']-matlabTotals['total'])) + " WHERE " + larger + " FINDS MORE ANNOTATIONS")
        print("PYTHON PERCENT NEAR: " + str(float(pythonTotals['near'])/float(pythonTotals['total'])))
        print("MATLAB PERCENT NEAR: " + str(float(matlabTotals['near'])/float(matlabTotals['total'])))
        for interaction in trueInts:
            near = "n" + interaction
            if float(pythonTotals[interaction]+pythonTotals[near] != 0):
                print("PYTHON NEAR/TOTAL FOR " + interaction + " INTERACTIONS: " + str(float(pythonTotals[near])/(float(pythonTotals[interaction]+pythonTotals[near]))))
            if float(matlabTotals[interaction]+matlabTotals[near] != 0):
                print("MATLAB NEAR/TOTAL FOR " + interaction + " INTERACTIONS: " + str(float(matlabTotals[near])/(float(matlabTotals[interaction]+matlabTotals[near]))))
            print("\n")

        matching = 0
        total = 0
        for pyInteraction in interaction_list:
            for mlInteraction in interaction_list:
                if pyInteraction == mlInteraction:
                    matching += confusionMatrix[pyInteraction][mlInteraction]
                    total += confusionMatrix[pyInteraction][mlInteraction] 
                else:
                    total += confusionMatrix[pyInteraction][mlInteraction] 
            total += confusionMatrix[pyInteraction]['blank']

        print("MATCHING: " + str(matching))
        if total != 0: 
            print("Percent Matching: " + str(float(matching)/float(total)))