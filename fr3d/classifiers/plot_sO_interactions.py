# -*- coding: utf-8 -*-
"""
    plot-sO-interactions.py reads a data file and plots points to represent
    oxygen-base stacking interactions.
"""

from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import pickle
from collections import defaultdict


from fr3d.localpath import outputText
from fr3d.localpath import outputNAPairwiseInteractions
from fr3d.localpath import outputNAPickleInteractions
from fr3d.localpath import contact_list_file
from fr3d.localpath import inputPath
from fr3d.localpath import outputHTML

from NA_pairwise_interactions import map_PDB_list_to_PDB_IFE_dict
from draw_residues import draw_base


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


if __name__=="__main__":

    base_seq_list = ['DA','DC','DG','DT']  # for DNA
    base_seq_list = ['A','C','G','U']      # for RNA

    interaction_lists = [["s3O2'","s3O3'","s3O4'","s3O5'","s3OP1","s3OP2","s5O2'","s5O3'","s5O4'","s5O5'","s5OP1","s5OP2"],
                        ["ns3O2'","ns3O3'","ns3O4'","ns3O5'","ns3OP1","ns3OP2","ns5O2'","ns5O3'","ns5O4'","ns5O5'","ns5OP1","ns5OP2"]]

    # plot one instance of each of the pairwise interactions
    PlotPair = False
    PlotPair = True
    AlreadyPlotted = {}

    unit_data_path = "C:/Users/zirbel/Documents/GitHub/fr3d-python/data/units"

    near_color = [1,0,0]  # red
    true_color = [0,0,0]  # black
    ring5_color = [0,1,1] # cyan
    ring6_color = [0,0,1] # blue
    O4_color = [1,0,0]  # red
    O2_color = [0,1,1]  # cyan
    OP_color = [0,0,0]  # black

    color_by_oxygen = True
    color_by_oxygen = False

    plot_true_in_3D = False

    PDB_list = ['4V9F','7K00']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.217/3.0A/csv']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.217/2.0A/csv']

    PDB_IFE_Dict = map_PDB_list_to_PDB_IFE_dict(PDB_list)

    all_PDB_ids = sorted(PDB_IFE_Dict.keys())
    print("Loading NA-pairwise-interactions from %d PDB files" % len(all_PDB_ids))

    pair_to_data = defaultdict(dict)

    # load output files from NA_pairwise_interactions
    for PDB in all_PDB_ids:
        pair_to_data_file = outputNAPairwiseInteractions + "%s_pairs_v1.pickle" % PDB
        print("Reading %s" % pair_to_data_file)
        try:
            new_dict = pickle.load(open(pair_to_data_file,'rb'))
            pair_to_data.update(new_dict)
        except:
            print("Not able to load annotations for %s" % PDB)

    # loop over sets of interactions, plotting points on bases
    for interaction_list in interaction_lists:

        fig = plt.figure(figsize=(10.0,8.0))

        # loop over individual bases, to make separate plots for each base
        for v, nt1_seq in enumerate(base_seq_list):

            if "n" in interaction_list[0]:
                ax = fig.add_subplot(2, 2, v+1)
            elif plot_true_in_3D:
                ax = fig.add_subplot(2, 2, v+1, projection='3d')
            else:
                ax = fig.add_subplot(2, 2, v+1)

            ax.axis("equal")

            # accumulate data specific to this base
            xvalues = []
            yvalues = []
            zvalues = []
            colors3d  = []
            colors2d  = []
            colorsoxy = []

            c = 0

            # loop over pairs, finding those with the desired base and interaction
            for pair,datapoint in pair_to_data.items():

                if not datapoint['nt1_seq'] == nt1_seq:
                    continue

                if not 'sOinteraction' in datapoint:
                    continue

                if abs(datapoint['sOz']) < 2:
                    continue

                if datapoint['sOinteraction'] in interaction_list:

                    c += 1
                    if c <= 100000:
                        xvalues.append(datapoint['sOx'])
                        yvalues.append(datapoint['sOy'])
                        zvalues.append(datapoint['sOz'])

                        if datapoint['sOinteraction'].startswith("n"):
                            colors3d.append(near_color)
                        else:
                            colors3d.append(true_color)

                        if 'ring5' in datapoint['sOring']:
                            colors2d.append(ring5_color)
                        elif 'ring6' in datapoint['sOring']:
                            colors2d.append(ring6_color)
                        else:
                            colors2d.append(near_color)

                        if datapoint['sOoxygen'] == "O4'":
                            colorsoxy.append(O4_color)
                        elif datapoint['sOoxygen'] in ["O2'","O3'"]:
                            colorsoxy.append(O2_color)
                        else:
                            colorsoxy.append(OP_color)

            # plot accumulated data points for this base
            if color_by_oxygen:
                colors2d = colorsoxy
                colors3d = colorsoxy

            if len(xvalues) < 500:
                markersize = 4
            else:
                markersize = 0.5

            if "n" in interaction_list[0]:
                # 2D plot for near interactions
                draw_base(nt1_seq,2,ax,zorder=1)
                ax.scatter(xvalues,yvalues,color=colors2d,marker=".",s=markersize,zorder=2)
                ax.set_title('Base %s with %s near sO interactions' % (nt1_seq,c))
            elif plot_true_in_3D:
                # 3D plot for true interactions
                draw_base(nt1_seq,3,ax)
                ax.scatter(xvalues,yvalues,zvalues,color=colors3d,marker=".",s=markersize)
                ax.scatter(xvalues,yvalues,[0 for i in zvalues],color=colors2d,marker=".",s=markersize)
                ax.set_title('Base %s with %s sO interactions' % (nt1_seq,c))
            else:
                # 2D plot for true interactions
                draw_base(nt1_seq,2,ax,zorder=1)
                ax.scatter(xvalues,yvalues,color=colors2d,marker=".",s=markersize,zorder=2)
                ax.set_title('Base %s with %s sO interactions' % (nt1_seq,c))

            print("Plotted %d points for %s" % (c,nt1_seq))

        # show all plots for this interaction_list
        if color_by_oxygen:
            print("Red for O4', cyan for O2' and O3', black for O5', O1P, O2P")
        elif "n" in interaction_list[0]:
            print("Cyan for points in or near the 5-sided ring, blue for 6-sided")
        else:
            print("Black for true interactions with z coordinate, cyan for 5-sided ring, blue for 6-sided")

        figManager = plt.get_current_fig_manager()
        figManager.full_screen_toggle()

        if "n" in interaction_list[0]:
            figure_save_file = outputNAPairwiseInteractions + "sO_near_%d" % len(all_PDB_ids)
        else:
            figure_save_file = outputNAPairwiseInteractions + "sO_true_%d" % len(all_PDB_ids)
        plt.savefig(figure_save_file+".png")
        plt.savefig(figure_save_file+".pdf")
        plt.close()



    # loop over sets of interactions, plotting histograms of z values
    interaction_list = ["s3O2'","s3O3'","s3O4'","s3O5'","s3OP1","s3OP2","s5O2'","s5O3'","s5O4'","s5O5'","s5OP1","s5OP2",
                        "ns3O2'","ns3O3'","ns3O4'","ns3O5'","ns3OP1","ns3OP2","ns5O2'","ns5O3'","ns5O4'","ns5O5'","ns5OP1","ns5OP2"]

    fig = plt.figure(figsize=(10.0,8.0))

    # loop over individual bases, to make separate plots for each base
    for v, nt1_seq in enumerate(base_seq_list):

        ax = fig.add_subplot(2, 2, v+1)

        # accumulate data specific to this base
        zvalues = []

        # loop over pairs, finding those with the desired base and interaction
        for pair,datapoint in pair_to_data.items():

            if not datapoint['nt1_seq'] == nt1_seq:
                continue

            if not 'sOinteraction' in datapoint:
                continue

            if abs(datapoint['sOz']) < 2:
                continue

            if datapoint['sOinteraction'] in interaction_list:
                if not "near" in datapoint['sOring']:
                    zvalues.append(abs(datapoint['sOz']))

        # plot histogram of z values for this base
        bins = [2.7+0.05*k for k in range(0,21)]
        plt.hist(zvalues,bins=bins)
        ax.set_title('Base %s with %s true and near sO interactions' % (nt1_seq,len(zvalues)))

        print("Plotted %d points for %s" % (len(zvalues),nt1_seq))

    figManager = plt.get_current_fig_manager()
    figManager.full_screen_toggle()
    figure_save_file = outputNAPairwiseInteractions + "sO_z_histogram_%d" % len(all_PDB_ids)
    plt.savefig(figure_save_file+".png")
    plt.savefig(figure_save_file+".pdf")
    plt.close()

