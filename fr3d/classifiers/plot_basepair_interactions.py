"""
    plot-basepair-interactions.py reads a data file and plots points to represent
    base pairing interactions.
    Points are colored according to Python and Matlab annotations.
"""

import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from collections import defaultdict

import math
import os
import sys
if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
else:
    from urllib.request import urlretrieve as urlretrieve

from class_limits import basepair_cutoffs
from NA_pairwise_interactions import map_PDB_list_to_PDB_IFE_dict
from draw_residues import draw_base

from fr3d.localpath import outputText
from fr3d.localpath import outputNAPairwiseInteractions
from fr3d.localpath import fr3d_pickle_path
from fr3d.localpath import contact_list_file
from fr3d.localpath import inputPath
from fr3d.localpath import outputHTML
from fr3d.localpath import storeMatlabFR3DPairs

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
        urlretrieve("http://rna.bgsu.edu/pairs/"+pairsFileName, pathAndFileName) # testing

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


""" code for loading basepair annotations

            # compare to previous annotations, may be machine specific

            pathAndFileName = "C:/Users/zirbel/Documents/FR3D/Python FR3D/data/pairs/%s_RNA_pairs.pickle" % PDB

            if not os.path.exists(pathAndFileName):
                pairsFileName = PDB + '_RNA_pairs.pickle'
                print("Downloading "+pairsFileName)
                if sys.version_info[0] < 3:
                    urllib.urlretrieve("http://rna.bgsu.edu/pairs/"+pairsFileName, pathAndFileName)  # python 2
                else:
                    urllib.request.urlretrieve("http://rna.bgsu.edu/pairs/"+pairsFileName, pathAndFileName)  # python 3


            all_pair_types= ['cWW', 'tWW', 'cWH', 'tWH', 'cWS', 'tWS', 'cHH', 'tHH', 'cHS', 'tHS', 'cSS', 'tSS']
            Matlab_pairs = load_basepair_annotations(pathAndFileName,all_pair_types)

            comparison_filename = "C:/Users/zirbel/Documents/FR3D/Python FR3D/comparison/comparison_%s.txt" % PDB
            compare_annotations(Matlab_pairs,Python_pairs,all_pair_types,comparison_filename)

"""



def compare_annotations(pairs_1,pairs_2,all_pair_types,filename):
    """
    Old code to compare annotation systems
    """

    allkeys = set(pairs_1.keys()) | set(pairs_1.keys())

    with open(filename, mode='w') as file:

        file.write("Unit id 1"+"\t"+"Unit id 2"+"\t"+"Base Combination"+"\t"+"Annotation 1"+"\t"+"Annotation 2"+"\t"+"Match Count"+"\t"+"url"+"\n")

        total_annotations  = 0
        count_agreements = 0
        for key in allkeys:
            #print(key)
            x = key[0].split("|")
            y = key[1].split("|")
            base_combination = x[3]+y[3]

            pair_type=""
            if key in pairs_1:
                pair_type = pairs_1[key].replace("n","")

            elif key in pairs_2:
                pair_type = pairs_2[key]

            #print([pair_type])
            if pair_type in ['cHW', 'tHW', 'cSH', 'tSH', 'cSW', 'tSW']:
                continue

            if pair_type in ['cWW', 'tWW','cSS', 'tSS', 'cHH', 'tHH', '0  ']:
                if base_combination in ['UA', 'CG', 'UG', 'GA', 'CA', 'UC']:
                    continue
                if base_combination in ['AA', 'CC', 'GG', 'UU']:
                    if int(x[4])>=int(y[4]):
                        continue

            if key in pairs_1:
                ann1 = pairs_1[key]
            else:
                ann1 = ""
            if key in pairs_2:
                ann2 = pairs_2[key]
            else:
                ann2 = ""

            ## Keep track of how often ann1 and ann2 agree over all pairs
            ## New Python annotations use upper and lowercase to indicate
            ## details of which part of the edge is used with cWw and tHh and such,
            ## so change to lowercase to compare

            if len(ann1) > 0 or len(ann2) > 0:
                total_annotations += 1

            if ann1.lower() == ann2.lower() and len(ann2) > 0:
                count_agreements +=  1

            ## Counting the number of times the annotations match with each other
            ## for this particular pair of neucleotides
            count = 0
            if ann1.lower() == ann2.lower()  and len(ann2) > 0:
                count = 1
            elif ann2.lower() in ann1.lower() and len(ann2) > 0:
                count += 0.6

            url= "http://rna.bgsu.edu/rna3dhub/display3D/unitid/" + key[0] + "," +  key[1]

            if count < 1:
                file.write(key[0]+"\t"+key[1]+"\t"+x[3]+y[3]+"\t"+ann1+"\t"+ann2+"\t"+str(count)+"\t"+url+"\n")

    print("  Number of times ann1 and ann2 agree %d" % count_agreements)
    print("  Total number of ann1 or ann2 annotations %d" % total_annotations)
    print("  Percentage of times ann1 and ann2 agree %0.4f" % (count_agreements*100/total_annotations))


def reverse(pair):
    return (pair[1],pair[0])

def check_full_data(datapoint):

    return ('x' in datapoint and 'y' in datapoint and 'z' in datapoint \
            and 'gap12' in datapoint and 'angle_in_plane' in datapoint and 'normal_Z' in datapoint)


def print_datapoint(datapoint):

#    if 'url' in datapoint:
#        print("  %s" % datapoint['url'])

    fields = ['x','y','z','gap12','angle_in_plane','normal_Z']

    for f in fields:
        if f in datapoint:
            print("  %20s = %11.6f" % (f,datapoint[f]))

    return None


def plot_basepair_cutoffs(base_combination,lowercase_list,ax,variables):

    color = ["#BA55D3","#63B8FF","#00EE76","#FF8C00","#CDC9A5","#8A8A8A","#8A8A8A","#8A8A8A","#8A8A8A","#8A8A8A","#8A8A8A"]
    cc = 0

    for bc in basepair_cutoffs.keys():
        if bc == base_combination:
            for interaction in basepair_cutoffs[bc].keys():
                if interaction.lower().replace("a","") in lowercase_list:
                    for subcategory in basepair_cutoffs[bc][interaction].keys():
                        limits = basepair_cutoffs[bc][interaction][subcategory]
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


if __name__=="__main__":

    base_seq_list = ['DA','DC','DG','DT']  # for DNA
    base_seq_list = ['A','C','G','U']      # for RNA

    symmetric_base_combination_list = ['A,A','C,C','G,G','U,U']
    symmetric_basepair_list = ['cWW','tWW','cHH','tHH','cSS','tSS']

    base_combination_list = ['A,A','A,C','A,G','A,U','C,C','G,C','C,U','G,G','G,U','U,U']

    interaction_lists = [["cWW","cWw","cwW"],["tWW"]]
    interaction_lists = [["cWW"],["tWW"],["cWH"],["tWH"],["cHW"],["tHW"],["cWS"],["tWS"],["cSW"],["tSW"],["cHH"],["tHH"],["cHS"],["tHS"],["cSH"],["tSH"],["cSS"],["tSS"]]

    unit_data_path = "C:/Users/zirbel/Documents/GitHub/fr3d-python/data/units"

    red   = [1,0,0] # red
    black = [0,0,0] # black
    cyan  = [0,1,1] # cyan
    blue  = [0,0,1] # blue

    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.216/2.0A/csv']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.216/3.0A/csv']
    PDB_list = ['4V9F','7K00']

    PDB_IFE_Dict = map_PDB_list_to_PDB_IFE_dict(PDB_list)

    # load all datapoints on pairs of bases, whether annotated as paired or not
    all_PDB_ids = PDB_IFE_Dict.keys()
    print("Loading NA-pairwise-interactions from %d PDB files" % len(all_PDB_ids))

    PDB_skip_set = set(['1R9F','5NXT','4KTG'])

    # load annotations of these PDB files from the BGSU RNA server
    print('Loading Matlab annotations of these files')
    Matlab_annotation_to_pair = defaultdict(list)
    pair_to_Matlab_annotation = defaultdict(str)
    for PDB_id in all_PDB_ids:
        interaction_to_pairs = load_Matlab_FR3D_pairs(PDB_id)
        num_pairs = 0
        for interaction in interaction_to_pairs.keys():
            # store all basepairs and only basepairs
            if "c" in interaction or "t" in interaction:
                num_pairs += 1
                Matlab_annotation_to_pair[interaction] += interaction_to_pairs[interaction]
                for pair in interaction_to_pairs[interaction]:
                    pair_to_Matlab_annotation[pair] = interaction
        if num_pairs == 0:
            print("No Matlab-annotated pairs in %s" % PDB_id)
            PDB_skip_set.add(PDB_id)

    #print("Skipping %d PDB files because they have no Matlab annotation to compare to" % len(PDB_skip_set))
    print("Found these Matlab annotations: %s" % sorted(Matlab_annotation_to_pair.keys()))

    all_PDB_ids = list(set(all_PDB_ids) - PDB_skip_set)
    print("Loading Python annotations from %d PDB files" % len(all_PDB_ids))

    pair_to_data = defaultdict(dict)

    # load output files from NA_pairwise_interactions
    for PDB in all_PDB_ids:
        pair_to_data_file = outputNAPairwiseInteractions + "%s_pairs_v1.pickle" % PDB
        #print("Reading %s" % pair_to_data_file)
        try:
            new_dict = pickle.load(open(pair_to_data_file,'rb'))
            pair_to_data.update(new_dict)
            if len(new_dict.keys()) == 0:
                print("No Python-annotated pairs in %s" % PDB)
        except:
            print("Not able to load annotations for %s" % PDB)

    # loop over specified sets of interactions
    for interaction_list in interaction_lists:

        lowercase_list = [i.lower() for i in interaction_list]

        # identify unit id pairs that are annotated as basepairing by Matlab code
        Matlab_pairs = []
        for interaction in interaction_list:
            Matlab_pairs += Matlab_annotation_to_pair[interaction]
        Matlab_pairs = set(Matlab_pairs)

        # loop over specified base combinations
        for base_combination in base_combination_list:

            # don't show AA cHW because AA cWH will be shown
            if base_combination in symmetric_base_combination_list and interaction in ['cHW','tHW','cSW','tSW','cSH','tSH']:
                continue

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

            # loop over pairs for which we have data, finding those with the interaction and base combination
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

                have_data_in_other_order = False

                # check different annotation schemes to decide how to show this datapoint
                # check Python annotation that generated pair_to_data
                if 'basepair' in datapoint and datapoint['basepair'].lower().replace("a","") in lowercase_list:
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

                if Matlab and not Python and base_combination in symmetric_base_combination_list \
                   and interaction_list[0] in symmetric_basepair_list \
                   and 'basepair' in r_datapoint and r_datapoint['basepair'].lower().replace("a","") in lowercase_list:
                   skip_this_order = True
                else:
                    skip_this_order = False


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

                if (Python or Matlab) and have_full_data and not skip_this_order:

                    c += 1
                    xvalues.append(datapoint['x'])
                    yvalues.append(datapoint['y'])
                    rvalues.append(math.sqrt(datapoint['x']**2 + datapoint['y']**2))
                    tvalues.append(math.atan2(datapoint['y'],datapoint['x'])*180/3.141592654)
                    zvalues.append(datapoint['z'])
                    gvalues.append(datapoint['gap12'])
                    avalues.append(datapoint['angle_in_plane'])
                    nvalues.append(datapoint['normal_Z'])

                    if 'basepair' in datapoint:
                        basepair = datapoint['basepair']
                    else:
                        basepair = "   "

                    if Python and Matlab:
                        color = black
                        size = 1
                    elif Python and not Matlab:
                        color = blue
                        size = 40
                        print("blue = only Python:  %s Python %4s Matlab %4s %s" % (base_combination,basepair,pair_to_Matlab_annotation[pair],datapoint['url']))
                        print_datapoint(datapoint)
                    elif not Python and Matlab:
                        color = red
                        size = 40
                        print("red = only Matlab: %s Python %4s Matlab %4s %s" % (base_combination,basepair,pair_to_Matlab_annotation[pair],datapoint['url']))
                        print_datapoint(datapoint)
                    else:
                        color = [0,1,0]
                        size = 100
                    colors2d.append(color)
                    sizes.append(size)

            print("Plotted %5d points for %s %s" % (c,base_combination,",".join(interaction_list)))

            # make the figure
            if c > 0:
                fig = plt.figure(figsize=(11.0, 4.0))

                ax = fig.add_subplot(1, 4, 1)
                ax.axis("equal")
                plot_basepair_cutoffs(base_combination,lowercase_list,ax,1)
                ax.scatter(xvalues,yvalues,color=colors2d,marker=".",s=sizes)
                ax.set_title('x and y for %d %s %s' % (c,base_combination,interaction_list[0]))
                draw_base(nt1_seq,'default',2,ax)

                ax = fig.add_subplot(1, 4, 2)
                plot_basepair_cutoffs(base_combination,lowercase_list,ax,2)
                ax.scatter(tvalues,rvalues,color=colors2d,marker=".",s=sizes)
                ax.set_title('theta and radius')

                ax = fig.add_subplot(1, 4, 3)
                plot_basepair_cutoffs(base_combination,lowercase_list,ax,3)
                ax.scatter(avalues,nvalues,color=colors2d,marker=".",s=sizes)
                ax.set_title('angle and normal, blue=Matlab only')

                ax = fig.add_subplot(1, 4, 4)
                plot_basepair_cutoffs(base_combination,lowercase_list,ax,4)
                ax.scatter(gvalues,zvalues,color=colors2d,marker=".",s=sizes)
                ax.set_title('gap12 and z values, red=Python only')

                # show all plots for this interaction_list
                figManager = plt.get_current_fig_manager()
                figManager.full_screen_toggle()
                figure_save_file = outputNAPairwiseInteractions + "basepairs_%s_%s_%s.png" % (base_combination,interaction_list[0],len(all_PDB_ids))
                plt.savefig(figure_save_file)
                #plt.show()
                plt.close()

    print('Wrote files to %s' % outputNAPairwiseInteractions)
