# -*- coding: utf-8 -*-
"""
    plot-sO-interactions.py reads a data file and plots points to represent
    oxygen-base stacking interactions.
"""

from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import pickle
import math
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

    interaction_lists = [["ns3O2'","ns3O3'","ns3O4'","ns3O5'","ns3OP1","ns3OP2","ns5O2'","ns5O3'","ns5O4'","ns5O5'","ns5OP1","ns5OP2"]]

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
    ring5_color = [0,0,0] # dark gray
    ring6_color = [0,0,0] # black
    O4_color = [1,0,0]  # red
    O2_color = [0,1,1]  # cyan
    OP_color = [0,0,0]  # black

    color_by_oxygen = True
    color_by_oxygen = False

    plot_true_in_3D = False

    PDB_list = ['4V9F','7K00']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.221/1.5A/csv']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.221/2.0A/csv']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.221/2.5A/csv']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.221/3.0A/csv']
    PDB_list = ['4V9F']

    fields = PDB_list[0].split("/")
    if len(fields) > 7:
        release    = PDB_list[0].split("/")[6]
        resolution = PDB_list[0].split("/")[7]
    else:
        release = ""
        resolution = PDB_list[0]

    PDB_IFE_Dict = map_PDB_list_to_PDB_IFE_dict(PDB_list)

    all_PDB_ids = sorted(PDB_IFE_Dict.keys())
    print("Loading NA-pairwise-interactions from %d PDB files" % len(all_PDB_ids))

    pair_to_data = defaultdict(dict)

    fileset = "_".join([release,resolution,str(len(all_PDB_ids))])

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

                if datapoint['sOinteraction'] in interaction_list:

                    """
                    nearring5 = False
                    nearring6 = False

                    #if datapoint['nt1_seq'] == 'G':
                    #    print("%3d %5s %s" % (c,datapoint['sOinteraction'],datapoint['url']))

                    parent1 = nt1_seq
                    x = datapoint['sOx']
                    y = datapoint['sOy']

                    if parent1 == 'A' or parent1 == 'DA':
                        if -1.302671*x + -0.512161*y + -0.512114 > 0:  # Left of C4-C5
                            if 1.033454*(x-(-1.138126))**2 + 0.143656*(x-(-1.138126))*(y-(-0.650781)) + (y-(-0.650781))**2 < 2.163590:  # A5 r=0.3
                                if -0.014382*x + -1.379291*y +  0.934116 > 0:  # Within  0.400000 Angstroms of being left of C5-N7
                                    if  1.286593*x + -0.316949*y +  3.047381 > 0:  # Within  0.400000 Angstroms of being left of N7-C8
                                        if  0.833587*x +  1.089911*y +  3.461823 > 0:  # Within  0.400000 Angstroms of being left of C8-N9
                                            if -0.803127*x +  1.118490*y +  1.698266 > 0:  # Within  0.400000 Angstroms of being left of N9-C4
                                                nearring5 = True
                        else:
                            if 1.001608*(x-(0.850305))**2 + 0.169100*(x-(0.850305))*(y-(-0.017921)) + (y-(-0.017921))**2 < 2.766745:  # A6 r=0.3
                                if  0.363524*x +  1.290539*y +  1.850003 > 0:  # Within  0.400000 Angstroms of being left of C4-N3
                                    if -1.076359*x +  0.793555*y +  3.030628 > 0:  # Within  0.400000 Angstroms of being left of N3-C2
                                        if -1.308429*x + -0.337740*y +  3.174043 > 0:  # Within  0.400000 Angstroms of being left of C2-N1
                                            if -0.319116*x + -1.301200*y +  2.398333 > 0:  # Within  0.400000 Angstroms of being left of N1-C6
                                                if  1.037709*x + -0.957315*y +  1.358356 > 0:  # Within  0.400000 Angstroms of being left of C6-C5
                                                    nearring6 = True
                    elif parent1 == 'C' or parent1 == 'DC':
                        if 0.867183*(x-(-0.298275))**2 + 0.040055*(x-(-0.298275))*(y-(-0.153209)) + (y-(-0.153209))**2 < 2.652492:  # C r=0.3
                            if -0.599253*x +  1.289335*y +  2.254778 > 0:  # Within  0.400000 Angstroms of being left of N1-C2
                                if -1.378522*x +  0.022802*y +  1.824412 > 0:  # Within  0.400000 Angstroms of being left of C2-N3
                                    if -0.676851*x + -1.128767*y +  1.713684 > 0:  # Within  0.400000 Angstroms of being left of N3-C4
                                        if  0.596389*x + -1.312333*y +  2.229696 > 0:  # Within  0.400000 Angstroms of being left of C4-C5
                                            if  1.359882*x + -0.033090*y +  2.615895 > 0:  # Within  0.400000 Angstroms of being left of C5-C6
                                                if  0.698355*x +  1.162053*y +  2.533245 > 0:  # Within  0.400000 Angstroms of being left of C6-N1
                                                    nearring6 = True
                    elif parent1 == 'G' or parent1 == 'DG':
                        if -1.306197*x + -0.492373*y + -0.896488 > 0:  # Left of C4-C5
                            if 1.032607*(x-(-1.476126))**2 + 0.129895*(x-(-1.476126))*(y-(-0.541964)) + (y-(-0.541964))**2 < 2.157145:  # G5
                                if -0.023230*x + -1.376606*y +  1.061418 > 0:  # Within  0.400000 Angstroms of being left of C5-N7
                                    if  1.278249*x + -0.337248*y +  3.488941 > 0:  # Within  0.400000 Angstroms of being left of N7-C8
                                        if  0.841883*x +  1.088640*y +  3.640461 > 0:  # Within  0.400000 Angstroms of being left of C8-N9
                                            if -0.790705*x +  1.117587*y +  1.308988 > 0:  # Within  0.400000 Angstroms of being left of N9-C4
                                                nearring5 = True
                        else:
                            if 1.082495*(x-(0.521747))**2 + 0.260413*(x-(0.521747))*(y-(0.023305)) + (y-(0.023305))**2 < 2.920747:  # G6 r=0.3
                                if  0.449709*x +  1.286231*y +  1.882380 > 0:  # Within  0.400000 Angstroms of being left of C4-N3
                                    if -0.992445*x +  0.855594*y +  2.637045 > 0:  # Within  0.400000 Angstroms of being left of N3-C2
                                        if -1.324604*x + -0.362005*y +  2.800178 > 0:  # Within  0.400000 Angstroms of being left of C2-N1
                                            if -0.533023*x + -1.330285*y +  2.599839 > 0:  # Within  0.400000 Angstroms of being left of N1-C6
                                                if  1.094166*x + -0.941908*y +  1.849906 > 0:  # Within  0.400000 Angstroms of being left of C6-C5
                                                    nearring6 = True
                    elif parent1 == 'DT':
                        if 0.959551*(x-(0.029169))**2 + 0.128151*(x-(0.029169))*(y-(-0.304375)) + (y-(-0.304375))**2 < 2.766276:  # DT r=0.3
                            if -0.675137*x +  1.198579*y +  2.604225 > 0:  # Within  0.400000 Angstroms of being left of N1-C2
                                if -1.365448*x + -0.109817*y +  2.181667 > 0:  # Within  0.400000 Angstroms of being left of C2-N3
                                    if -0.742906*x + -1.165341*y +  1.851614 > 0:  # Within  0.400000 Angstroms of being left of N3-C4
                                        if  0.767749*x + -1.221287*y +  1.936161 > 0:  # Within  0.400000 Angstroms of being left of C4-C5
                                            if  1.338191*x +  0.092630*y +  2.137070 > 0:  # Within  0.400000 Angstroms of being left of C5-C6
                                                if  0.677551*x +  1.205236*y +  2.512771 > 0:  # Within  0.400000 Angstroms of being left of C6-N1
                                                    nearring6 = True
                    elif parent1 == 'U':
                        if 0.912164*(x-(-0.302801))**2 + 0.143626*(x-(-0.302801))*(y-(-0.157137)) + (y-(-0.157137))**2 < 2.752991:  # U r=0.3
                            if -0.589251*x +  1.260286*y +  2.272756 > 0:  # Within  0.400000 Angstroms of being left of N1-C2
                                if -1.384641*x + -0.064970*y +  1.787427 > 0:  # Within  0.400000 Angstroms of being left of C2-N3
                                    if -0.834465*x + -1.135313*y +  1.810304 > 0:  # Within  0.400000 Angstroms of being left of N3-C4
                                        if  0.745842*x + -1.256133*y +  2.408409 > 0:  # Within  0.400000 Angstroms of being left of C4-C5
                                            if  1.352820*x +  0.018369*y +  2.590846 > 0:  # Within  0.400000 Angstroms of being left of C5-C6
                                                if  0.709695*x +  1.177761*y +  2.565310 > 0:  # Within  0.400000 Angstroms of being left of C6-N1
                                                    nearring6 = True

                    #if nearring5 or nearring6:
                    """

                    if True:
                        xvalues.append(datapoint['sOx'])
                        yvalues.append(datapoint['sOy'])
                        zvalues.append(datapoint['sOz'])

                        if datapoint['sOinteraction'].startswith("n"):
                            colors2d.append(near_color)
                            colors3d.append(near_color)
                        else:
                            colors3d.append(true_color)
                            if 'ring5' in datapoint['sOring']:
                                colors2d.append(ring5_color)
                            elif 'ring6' in datapoint['sOring']:
                                colors2d.append(ring6_color)

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

            c = len(xvalues)

            if len(xvalues) < 500:
                markersize = 4
            else:
                markersize = 0.5

            if "n" in interaction_list[0]:
                # 2D plot for near interactions
                draw_base(nt1_seq,'CPK',2,ax,zorder=1)

                """
                # Draw ellipses that coincide with the outer atoms of each ring
                x = []
                y = []
                if nt1_seq == "A":
                    for i in range(0,1000):       # A5
                        t = 6.28318530718*i/1000
                        r = math.sqrt(1.370973/(1.033354*math.cos(t)**2+0.144129*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(-1.138126 + r*math.cos(t))
                        y.append(-0.650781 + r*math.sin(t))
                    ax.plot(x,y,color="black")
                    x = []
                    y = []
                    for i in range(0,1000):       # A6
                        t = 6.28318530718*i/1000
                        r = math.sqrt(1.858876/(1.001635*math.cos(t)**2+0.169630*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(0.850305 + r*math.cos(t))
                        y.append(-0.017921 + r*math.sin(t))
                if nt1_seq == "C":
                    for i in range(0,1000):       # C
                        t = 6.28318530718*i/1000
                        r = math.sqrt(1.764725/(0.864524*math.cos(t)**2+0.036240*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(-0.298275 + r*math.cos(t))
                        y.append(-0.153209 + r*math.sin(t))
                if nt1_seq == "G":
                    for i in range(0,1000):       # G5
                        t = 6.28318530718*i/1000
                        r = math.sqrt(1.365959/(1.033126*math.cos(t)**2+0.130157*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(-1.476126 + r*math.cos(t))
                        y.append(-0.541964 + r*math.sin(t))
                    ax.plot(x,y,color="black")
                    x = []
                    y = []
                    for i in range(0,1000):       # G6
                        t = 6.28318530718*i/1000
                        r = math.sqrt(1.984003/(1.082449*math.cos(t)**2+0.263327*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(0.521747 + r*math.cos(t))
                        y.append(0.023305 + r*math.sin(t))
                if nt1_seq == "U":
                    for i in range(0,1000):       # U
                        t = 6.28318530718*i/1000
                        r = math.sqrt(1.847624/(0.912165*math.cos(t)**2+0.144395*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(-0.302801 + r*math.cos(t))
                        y.append(-0.157137 + r*math.sin(t))
                """

                """
                # draw ellipses that are 0.4 Angstroms outside the ring
                x = []
                y = []
                if nt1_seq == "A":
                    for i in range(0,1000):       # A5
                        t = 6.28318530718*i/1000
                        r = math.sqrt(2.467624/(1.032927*math.cos(t)**2+0.143607*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(-1.138126 + r*math.cos(t))
                        y.append(-0.650781 + r*math.sin(t))

                    ax.plot(x,y,color="black")
                    x = []
                    y = []
                    for i in range(0,1000):       # A6
                        t = 6.28318530718*i/1000
                        r = math.sqrt(3.110185/(1.001835*math.cos(t)**2+0.169297*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(0.850305 + r*math.cos(t))
                        y.append(-0.017921 + r*math.sin(t))

                if nt1_seq == "C":
                    for i in range(0,1000):       # C
                        t = 6.28318530718*i/1000
                        r = math.sqrt(2.987450/(0.864229*math.cos(t)**2+0.035866*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(-0.298275 + r*math.cos(t))
                        y.append(-0.153209 + r*math.sin(t))

                if nt1_seq == "G":
                    for i in range(0,1000):       # G5
                        t = 6.28318530718*i/1000
                        r = math.sqrt(2.461152/(1.032907*math.cos(t)**2+0.129883*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(-1.476126 + r*math.cos(t))
                        y.append(-0.541964 + r*math.sin(t))
                    ax.plot(x,y,color="black")
                    x = []
                    y = []
                    for i in range(0,1000):       # G6
                        t = 6.28318530718*i/1000
                        r = math.sqrt(3.273476/(1.082660*math.cos(t)**2+0.260087*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(0.521747 + r*math.cos(t))
                        y.append(0.023305 + r*math.sin(t))
                if nt1_seq == "U":
                    for i in range(0,1000):       # U
                        t = 6.28318530718*i/1000
                        r = math.sqrt(3.095314/(0.912089*math.cos(t)**2+0.144055*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(-0.302801 + r*math.cos(t))
                        y.append(-0.157137 + r*math.sin(t))
                """

                # draw ellipses that are 0.3 Angstroms outside the ring
                x = []
                y = []

                if nt1_seq == "A":
                    for i in range(0,1000):       # A5
                        t = 6.28318530718*i/1000
                        r = math.sqrt(2.163590/(1.033454*math.cos(t)**2+0.143656*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(-1.138126 + r*math.cos(t))
                        y.append(-0.650781 + r*math.sin(t))
                    ax.plot(x,y,color="green",zorder=2)
                    x = []
                    y = []

                    for i in range(0,1000):       # A6
                        t = 6.28318530718*i/1000
                        r = math.sqrt(2.766745/(1.001608*math.cos(t)**2+0.169100*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(0.850305 + r*math.cos(t))
                        y.append(-0.017921 + r*math.sin(t))
                if nt1_seq == "C":
                    for i in range(0,1000):       # C
                        t = 6.28318530718*i/1000
                        r = math.sqrt(2.652492/(0.867183*math.cos(t)**2+0.040055*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(-0.298275 + r*math.cos(t))
                        y.append(-0.153209 + r*math.sin(t))
                if nt1_seq == "G":
                    for i in range(0,1000):       # G5
                        t = 6.28318530718*i/1000
                        r = math.sqrt(2.157145/(1.032607*math.cos(t)**2+0.129895*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(-1.476126 + r*math.cos(t))
                        y.append(-0.541964 + r*math.sin(t))
                    ax.plot(x,y,color="green",zorder=2)
                    x = []
                    y = []

                    for i in range(0,1000):       # G6
                        t = 6.28318530718*i/1000
                        r = math.sqrt(2.920747/(1.082495*math.cos(t)**2+0.260413*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(0.521747 + r*math.cos(t))
                        y.append(0.023305 + r*math.sin(t))
                if nt1_seq == "U":
                    for i in range(0,1000):       # U
                        t = 6.28318530718*i/1000
                        r = math.sqrt(2.752991/(0.912164*math.cos(t)**2+0.143626*math.cos(t)*math.sin(t)+math.sin(t)**2))
                        x.append(-0.302801 + r*math.cos(t))
                        y.append(-0.157137 + r*math.sin(t))
                ax.plot(x,y,color="green",zorder=2)

                ax.scatter(xvalues,yvalues,color=colors2d,marker=".",s=markersize,zorder=3)
                ax.set_title('Base %s with %s near sO interactions' % (nt1_seq,c))
            elif plot_true_in_3D:
                # 3D plot for true interactions
                draw_base(nt1_seq,'CPK',3,ax)
                ax.scatter(xvalues,yvalues,zvalues,color=colors3d,marker=".",s=markersize)
                ax.scatter(xvalues,yvalues,[0 for i in zvalues],color=colors2d,marker=".",s=markersize)
                ax.set_title('Base %s with %s sO interactions' % (nt1_seq,c))
            else:
                # 2D plot for true interactions
                draw_base(nt1_seq,'CPK',2,ax,zorder=1)
                ax.scatter(xvalues,yvalues,color=colors2d,marker=".",s=markersize,zorder=3)
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
            figure_save_file = outputNAPairwiseInteractions + "sO_%s_near_ellipse" % fileset
        else:
            figure_save_file = outputNAPairwiseInteractions + "sO_%s_true" % fileset
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
    figure_save_file = outputNAPairwiseInteractions + "sO_%s_z_histogram" % fileset
    plt.savefig(figure_save_file+".png")
    plt.savefig(figure_save_file+".pdf")
    plt.close()

    print("Plots are in %s" % outputNAPairwiseInteractions)