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

from fr3d.localpath import outputText
from fr3d.localpath import outputNAPairwiseInteractions
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
def make_confusion_matrix(confusionMatrix, interaction_list):
    """
    Tabular construction of confusion matrix created from a dictionary of dictionaries of values.
    """

    output = "      "

    # header line
    for interaction in interaction_list:
        output += '%6s' % interaction

    output += '\n'

    for interaction in interaction_list:
        output += '%6s' % interaction
        for interaction2 in interaction_list:
            output += '%6d' % confusionMatrix[interaction][interaction2]
        output += '\n'

    return output

#=======================================================================
def focused_confusion_matrix(confusionMatrix, interaction_list):
    """
    List the non-zero elements of the confusion matrix.  Easier to focus then.
    """

    output = ""

    for interaction in interaction_list:
        for interaction2 in interaction_list:
            if confusionMatrix[interaction][interaction2] > 0:
                output += 'Matlab %6s is Python %6s %6d times' % (interaction,interaction2,confusionMatrix[interaction][interaction2])
                if not interaction == interaction2:
                    if interaction == '' and interaction2.startswith('n'):
                        output += '  <--- nothing/near'
                    elif interaction2 == '' and interaction.startswith('n'):
                        output += '  <--- near/nothing'
                    elif interaction in interaction2:
                        output += '  <--- true/near'
                    elif interaction2 in interaction:
                        output += '  <--- near/true'
                    elif interaction.startswith('n') and interaction2.startswith('n'):
                        output += '  <--- near/near'
                    elif interaction.startswith('4') and interaction2.startswith('3'):
                        output += '  <--- 4 becomes 3'
                    elif interaction.startswith('4') and interaction2.startswith('5'):
                        output += '  <--- 4 becomes 5'
                    elif interaction.startswith('8') and interaction2.startswith('7'):
                        output += '  <--- 8 becomes 7'
                    elif interaction.startswith('8') and interaction2.startswith('9'):
                        output += '  <--- 8 becomes 9'
                    else:
                        output += '  <--- more significant change'
                output += '\n'

    return output

#=======================================================================
def plot_unmatched_pairs(interactionDict, interaction_list, ordering):
    """Creates a 2x9 plot to show which annotations were found by passed in language and not found by other passed in language"""
    if ordering == "python":
        other = "matlab"
    elif ordering == "matlab":
        other = "python"
    if 'blank' in interaction_list:
        interaction_list.remove('blank')

    print("\n\nAnnotations that were found by %s and not by %s" % (ordering, other))

    # print_function = getattr(__builtin__, 'print')
    # border = ['----','----','----','----','----','----','----','----', '----','----','----']
    # print_function(*interaction_list, sep='\t')
    # print_function(*border, sep = "\t")
    # for category in interaction_list:
    #     print_function(interactionDict[category], end = "\t")

    print(interaction_list)
    print(border)
    for category in interaction_list:
        print(interactionDict[category], end = "\t")

    print("\n")

######
if __name__=="__main__":

    DNA = False
    make_plots = False
    make_plots = True

    PDB_list = ['7K00']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.237/2.5A/csv']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.285/2.0A/csv']
    PDB_list = ['4V9F']

    #PDB LIST and Skip Files####
    PDB_IFE_Dict = map_PDB_list_to_PDB_IFE_dict(PDB_list)

    all_PDB_ids = sorted(PDB_IFE_Dict.keys())
    print("Loading NA-pairwise-interactions from %d PDB files" % len(all_PDB_ids))
    PDB_skip_set = set(['1R9F','5NXT','4KTG'])

    pair_to_data = defaultdict(dict)\

    # Load Matlab annotations for BPh and BR
    print('Loading Matlab annotations of these files')
    pair_to_Matlab_BPh_annotation = defaultdict(str)
    pair_to_Matlab_BR_annotation = defaultdict(str)
    PDB_IFE_Dict = map_PDB_list_to_PDB_IFE_dict(PDB_list)
    all_PDB_ids = PDB_IFE_Dict.keys()

    for PDB_id in all_PDB_ids:
        interaction_to_pairs = load_Matlab_FR3D_pairs(PDB_id)
        num_pairs = 0
        for interaction in interaction_to_pairs.keys():
            # store all basepairs and only basepairs
            if 'BPh' in interaction:
                num_pairs += 1
                for pair in interaction_to_pairs[interaction]:
                    pair_to_Matlab_BPh_annotation[pair] = interaction
            if 'BR' in interaction:
                num_pairs += 1
                for pair in interaction_to_pairs[interaction]:
                    pair_to_Matlab_BR_annotation[pair] = interaction
        if num_pairs == 0:
            print("No Matlab-annotated BPh or BR pairs in %s" % (PDB_id))
            PDB_skip_set.add(PDB_id)

    print("Skipping %d PDB files because they have no Matlab annotation to compare to" % len(PDB_skip_set))

    # Load FR3D Python Annotations##
    # load output files from NA_pairwise_interactions
    all_PDB_ids = list(set(all_PDB_ids) - PDB_skip_set)
    print("Loading Python annotations from %d PDB files" % len(all_PDB_ids))

    pair_to_data = defaultdict(dict)

    not_loaded = []
    for PDB in all_PDB_ids:
        pair_to_data_file = outputNAPairwiseInteractions + "%s_pairs_v1.pickle" % PDB
        print("Reading %s" % pair_to_data_file)
        try:
            new_dict = pickle.load(open(pair_to_data_file,'rb'))
            for pair,datapoint in new_dict.items():
                if 'BPh' in datapoint or 'BR' in datapoint:
                    pair_to_data[pair] = datapoint
        except:
            print("Not able to load annotations for %s" % PDB)
            not_loaded.append(PDB)

    if make_plots:
        near_color = [1,0,0]  # red
        true_color = [0,0,0]  # black
        red   = [1,0,0] # red
        black = [0,0,0] # black
        cyan  = [0,1,1] # cyan
        blue  = [0,0,1] # blue

    if DNA:
        base_seq_list = ['DA','DC','DG','DT']  # for DNA
    else:
        base_seq_list = ['A','C','G','U']      # for RNA

    interaction_list = []
    for bb_type in ['BPh','BR']:
        for i in range(10):
            # true interactions
            interaction_list.append("%d%s" % (i,bb_type))
            # near interactions
            interaction_list.append("n%d%s" % (i,bb_type))

        print('Processing these interactions:')
        print(interaction_list)

    for base in base_seq_list:
        # Data collection for confusion matrix #######
        # each main key represents annotations found in python. The Keys within represent annotations found in matlab #
        # This dictionary is used to detect instances where python and matlab agree vs where they disagree ####

        confusionMatrix = {}
        confusionMatrix['total'] = 0
        for interaction in interaction_list + ['']:
            confusionMatrix[interaction] = {}
            for interaction2 in interaction_list + ['']:
                confusionMatrix[interaction][interaction2] = 0

        # accumulate data specific to this interaction and base combination
        if make_plots:
            xvalues = []
            yvalues = []
            rvalues = []
            tvalues = []
            zvalues = []
            colors2d  = []
            sizes = []

        c = 0
        messages = []

        for pair,datapoint in pair_to_data.items():
            # restrict to the current base combination
            if not datapoint['nt1_seq'] == base:
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

            # this is the place to check for parent of modified nucleotide
            if not fields1[3] == base:
                continue

            python_BPh_annotation = ''
            python_BR_annotation = ''

            if 'BPh' in datapoint and datapoint['BPh']:
                python_BPh_annotation = datapoint['BPh']

            if 'BR' in datapoint and datapoint['BR']:
                python_BR_annotation = datapoint['BR']

            matlab_BPh_annotation = ''
            matlab_BR_annotation = ''

            if pair in pair_to_Matlab_BPh_annotation:
                matlab_BPh_annotation = pair_to_Matlab_BPh_annotation[pair]

            if pair in pair_to_Matlab_BR_annotation:
                matlab_BR_annotation = pair_to_Matlab_BR_annotation[pair]

            if python_BPh_annotation or python_BR_annotation or matlab_BPh_annotation or matlab_BR_annotation:
                c += 1

                if not python_BPh_annotation == matlab_BPh_annotation:
                    message = '%s Matlab %5s becomes Python %5s %s' % (base,matlab_BPh_annotation,python_BPh_annotation,datapoint['url'])
                    messages.append(message)

                if not python_BR_annotation == matlab_BR_annotation:
                    message = '%s Matlab %5s becomes Python %5s %s' % (base,matlab_BR_annotation,python_BR_annotation,datapoint['url'])
                    messages.append(message)

                #######
                # Add record of annotations by pair
                # confusionMatrix is a dictionary - key is Matlab annotation, value is dictionary where
                # key is Python annotation - values are numerical total
                # Thus, the dictionary is like a transition matrix from Matlab annotation to Python
                #######

                if matlab_BPh_annotation or python_BPh_annotation:
                    confusionMatrix[matlab_BPh_annotation][python_BPh_annotation] += 1
                if matlab_BR_annotation or python_BR_annotation:
                    confusionMatrix[matlab_BR_annotation][python_BR_annotation] += 1
                confusionMatrix['total'] += 1

                if make_plots and 'BPh_oxygen' in datapoint:
                    for point in datapoint['BPh_oxygen']:
                        if python_BPh_annotation == matlab_BPh_annotation and not python_BPh_annotation.startswith('n'):
                            color = black
                            size = 5
                        elif python_BPh_annotation == matlab_BPh_annotation:
                            color = cyan
                            size = 5
                            #print("blue = only Python:  %s Python %4s Matlab %4s %s" % (base,python_annotation,matlab_annotation,datapoint['url']))
                            #print_datapoint(datapoint)
                        # elif not python_BPh_annotation and matlab_BPh_annotation:
                        #     color = red
                        #     size = 10
                        #     #print("red = only Matlab: %s Python %4s Matlab %4s %s" % (base,python_annotation,matlab_annotation,datapoint['url']))
                        #     #print_datapoint(datapoint)
                        else:
                            color = red
                            size = 5

                        colors2d.append(color)
                        sizes.append(size)
                        xvalues.append(point[0])
                        yvalues.append(point[1])
                        zvalues.append(point[2])
                        #rvalues.append(math.sqrt(point[0]**2 + point[1]**2))
                        #tvalues.append(math.atan2(point[1],point[0])*180/3.141592654)

                if make_plots and 'BR_oxygen' in datapoint:
                    for point in datapoint['BR_oxygen']:
                        if python_BR_annotation == matlab_BR_annotation and not python_BR_annotation.startswith('n'):
                            color = black
                            size = 5
                        elif python_BR_annotation == matlab_BR_annotation:
                            color = cyan
                            size = 5
                            #print("blue = only Python:  %s Python %4s Matlab %4s %s" % (base,python_BP_annotation,matlab_BR_annotation,datapoint['url']))
                            #print_datapoint(datapoint)
                        # elif not python_BR_annotation and matlab_BR_annotation:
                        #     color = red
                        #     size = 10
                            #print("red = only Matlab: %s Python %4s Matlab %4s %s" % (base,python_BR_annotation,matlab_BR_annotation,datapoint['url']))
                            #print_datapoint(datapoint)
                        else:
                            color = red
                            size = 5
                        colors2d.append(color)
                        sizes.append(size)
                        xvalues.append(point[0])
                        yvalues.append(point[1])
                        zvalues.append(point[2])
                        #rvalues.append(math.sqrt(point[0]**2 + point[1]**2))
                        #tvalues.append(math.atan2(point[1],point[0])*180/3.141592654)

        print("\n")
        print("Found %5d base-backbone interactions made by %s" % (c,base))

        if c > 0 and make_plots:
            fig = plt.figure(figsize=(11.0, 6.0))

            ax = fig.add_subplot(1, 1, 1)
            ax.axis("equal")
            ax.scatter(xvalues,yvalues,color=colors2d,marker=".",s=sizes)
            ax.set_title('x and y for %d %s %s' % (c,base,interaction_list[0]), rotation=0)

            draw_base(base,'default',2,ax)

            # show all plots for this interaction_list
            figManager = plt.get_current_fig_manager()
            figManager.full_screen_toggle()
            figure_save_file = outputNAPairwiseInteractions + "%s_%s_%s_%s.png" % (bb_type,base,interaction_list[0],len(all_PDB_ids))

            print('Saved plot in %s' % figure_save_file)

            plt.savefig(figure_save_file)
            #plt.show()
            plt.close()

        ### Display confusion matrix ##
        cm = make_confusion_matrix(confusionMatrix, interaction_list+[''])
        cm = focused_confusion_matrix(confusionMatrix, interaction_list+[''])
        messages = sorted(messages)
        print(cm)
        print('\n'.join(messages))

        ## Store confusion matrix and messages in a text file
        filename = os.path.join(outputNAPairwiseInteractions, 'confusionMatrix_%s_%s.txt' % (base,bb_type))
        with open(filename, 'w') as f:
            f.writelines(cm)
            for message in messages:
                f.writelines(message + "\n")