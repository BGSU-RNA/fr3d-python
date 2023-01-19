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
import numpy as np
import random

if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
else:
    from urllib.request import urlretrieve as urlretrieve

from class_limits import nt_nt_cutoffs
from NA_pairwise_interactions import map_PDB_list_to_PDB_IFE_dict
from draw_residues import draw_base

from fr3d.localpath import outputText
from fr3d.localpath import outputNAPairwiseInteractions
from fr3d.localpath import fr3d_pickle_path
from fr3d.localpath import contact_list_file
from fr3d.localpath import inputPath
from fr3d.localpath import outputHTML
from fr3d.localpath import storeMatlabFR3DPairs

from orderBySimilarityTemp import treePenalizedPathLength

JS1 = '  <script src="./js/JSmol.min.nojq.js"></script>'
JS2 = '  <script src="./js/jquery.jmolTools.js"></script>'
JS3 = '  <script src="./js/imagehandlinglocal.js"></script>'
JS4 = '<script src="./js/jmolplugin.js" type="text/javascript"></script>'
JS5 = '<script type="text/javascript" src="./js/heatmap.js"></script>'
TEMPLATEPATH = '../search/'
OUTPUTPATH = "C:/Users/zirbel/Documents/FR3D/Python FR3D/output/"

def writeHTMLOutputInefficient(Q,candidates,allvsallmatrix=np.empty( shape=(0, 0) )):
    """
    Deprecated.
    Write the list of candidates in an HTML format that also shows
    the coordinate window and a heat map of all-against-all distances.
    """

    pairTypes = ['pairsStacks']
    pairsToPrint = defaultdict(list)

    pairsToPrint['pairsStacks'] = [(1,2)]

    pagetitle = "%s" % Q['name']

    htmlfilename = Q['name'].replace(" ","_")

    candidatelist = '<table style="white-space:nowrap;">\n'

    numPositions = 2

    sequence_column = 3

    # write header line
    candidatelist += "<tr><th>S.</th><th>Show</th>"

    sequence_column += 1
    for j in range(0,numPositions):
        candidatelist += "<th>Position %d</th>" % (j+1)
        sequence_column += 1
    candidatelist += "<th>Python</th>"
    candidatelist += "<th>Subcat</th>"
    candidatelist += "<th>New Python</th>"
    candidatelist += "<th>Matlab</th>"

    candidatelist += "<th>x</th>"
    candidatelist += "<th>y</th>"
    candidatelist += "<th>z</th>"
    candidatelist += "<th>gap12</th>"
    candidatelist += "<th>angle_in_plane</th>"
    candidatelist += "<th>normal_z</th>"
    candidatelist += "<th>max_distance</th>"
    candidatelist += "<th>min_angle</th>"
    candidatelist += "<th>max_badness</th>"

    candidatelist += "</tr>\n"

    # write one row for each candidate
    for i in range(0,len(candidates)):
        candidate = candidates[i]
        candidatelist += '<tr><td>'+str(i+1)+'.</td><td><label><input type="checkbox" id="'+str(i)+'" class="jmolInline" data-coord="'
        for j in range(0,numPositions):
            candidatelist += candidate[j]
            if j < numPositions-1:
                candidatelist += ','
        candidatelist += '">&nbsp</label></td>'

        PDB_id = candidate[0][0:4]

        """
        if PDB_id in Q["PDB_data_file"]:
            candidatelist += "<td>%s</td>" % format_resolution(Q["PDB_data_file"][PDB_id])
        else:
            candidatelist += "<td>NA</td>"
        """

        # write unit ids
        for j in range(0,numPositions):
            candidatelist += "<td>"+candidate[j]+"</td>"

        # write Python, Matlab interactions
        candidatelist += "<td>%s</td>" % candidate[11]  # python
        candidatelist += "<td>%d</td>" % candidate[12]  # python subcategory
        candidatelist += "<td>%s</td>" % candidate[13]  # new python
        candidatelist += "<td>%s</td>" % candidate[14]  # matlab

        for colnum in range(2,11):
            candidatelist += "<td>%0.2f</td>" % candidate[colnum]

        candidatelist += '</tr>\n'
    candidatelist += '</table>\n'

    discrepancydata = ''

    if np.size(allvsallmatrix) > 0:
        s = allvsallmatrix.shape[0]
        for c in range(0,s):
            ife1 = candidates[c][0]
            for d in range(0,s):
                ife2 = candidates[d][0]

                discrepancydata += '{"discrepancy": ' + "%.4f" % allvsallmatrix[c][d]
                discrepancydata += ', "ife1": "' + ife1 + '", "ife1_index": ' + str(c)
                discrepancydata += ', "ife2": "' + ife2 + '", "ife2_index": ' + str(d) + '}'
                if c < s-1 or d < s-1:
                    discrepancydata += ',\n'

            """
                      {"discrepancy": 0.0, "ife1": "3KIS|1|a", "ife1_index": 0, "ife2": "3KIS|1|a", "ife2_index": 0},
                      {"discrepancy": 0.14582092027413304, "ife1": "3KIS|1|a", "ife1_index": 0, "ife2": "3FIC|1|A", "ife2_index": 1},
                      {"discrepancy": 0.0093556396649322, "ife1": "3KIS|1|a", "ife1_index": 0, "ife2": "3KIQ|1|a", "ife2_index": 2},
                      {"discrepancy": 0.14582092027413304, "ife1": "3FIC|1|A", "ife1_index": 1, "ife2": "3KIS|1|a", "ife2_index": 0},
                      {"discrepancy": 0.0, "ife1": "3FIC|1|A", "ife1_index": 1, "ife2": "3FIC|1|A", "ife2_index": 1},
                      {"discrepancy": 0.14260495908457801, "ife1": "3FIC|1|A", "ife1_index": 1, "ife2": "3KIQ|1|a", "ife2_index": 2},
                      {"discrepancy": 0.0093556396649322, "ife1": "3KIQ|1|a", "ife1_index": 2, "ife2": "3KIS|1|a", "ife2_index": 0},
                      {"discrepancy": 0.14260495908457801, "ife1": "3KIQ|1|a", "ife1_index": 2, "ife2": "3FIC|1|A", "ife2_index": 1},
                      {"discrepancy": 0.0, "ife1": "3KIQ|1|a", "ife1_index": 2, "ife2": "3KIQ|1|a", "ife2_index": 2}
            """

    # read template.html into one string
    with open(TEMPLATEPATH + 'template.html', 'r') as myfile:
        template = myfile.read()

    # replace ###PAGETITLE### with pagetitle
    template = template.replace("###PAGETITLE###",pagetitle)

    sequence_column = 0   # column to set in fixed width font
    template = template.replace("###sequencecolumn###",str(sequence_column))

    queryNote = "Query name: %s.  Found %d candidates from %d of %d files in %0.0f seconds." % (Q['name'].encode('ascii','ignore'),len(candidates),Q["numFilesSearched"],len(Q["searchFiles"]),Q["elapsedCPUTime"])
    queryNote = "Query name: %s.  Found %d candidates from %d of %d files in %0.0f seconds." % (Q['name'],len(candidates),Q["numFilesSearched"],len(Q["searchFiles"]),Q["elapsedCPUTime"])

    if "moreCandidatesThanHeatMap" in Q:
        queryNote += " " + Q["moreCandidatesThanHeatMap"] + "\n"
    else:
        queryNote += "\n"

    template = template.replace("###QUERYNAME###",str(queryNote.encode('ascii','ignore')))

    seeModifyQuery = ''
    template = template.replace("###SEEMODIFYQUERY###",seeModifyQuery)

    template = template.replace("###seeCSVOutput###","")

    # replace ###CANDIDATELIST### with candidatelist
    template = template.replace("###CANDIDATELIST###",candidatelist)

    template = template.replace("###JS1###",JS1)
    template = template.replace("###JS2###",JS2)
    template = template.replace("###JS3###",JS3)
    template = template.replace("###JS4###",JS4)

    refresh = ""
    if "reloadOutputPage" in Q and Q["reloadOutputPage"]:
        refresh = '<meta http-equiv="refresh" content="%d">' % REFRESHTIME
    template = template.replace("###REFRESH###",refresh)

    if np.size(allvsallmatrix) > 0:
        template = template.replace("###JS5###",JS5)    # include heatmap.js code
        discrepancydata = "var data =  [\n" + discrepancydata + "]"
        discrepancydata = '<script type="text/javascript">\n' + discrepancydata + '\n</script>'
        template = template.replace("###DISCREPANCYDATA###",discrepancydata)
    else:
        #template = template.replace("###DISCREPANCYDATA###","")
        template = template.replace("###JS5###","")    # do not display a heat map

    outputfilename = os.path.join(OUTPUTPATH,htmlfilename+".html")

    print("Writing to %s" % outputfilename)

    messages = ""

    messages += "\n<br>"
    if len(Q["userMessage"]) > 0:
        messages += "User messages:<br>\n"
        for line in Q["userMessage"]:
            messages += line + "<br>\n"
    else:
        messages += "No error or warning messages.<br>\n"

    template = template.replace("###MESSAGES###",messages)

    with open(outputfilename, 'w') as myfile:
        myfile.write(template)

def writeHTMLOutput(Q,candidates,allvsallmatrix=np.empty( shape=(0, 0) )):
    """
    Write the list of candidates in an HTML format that also shows
    the coordinate window and a heat map of all-against-all distances.
    Write the heatmap data in an efficient way.
    """

    pairTypes = ['pairsStacks']
    pairsToPrint = defaultdict(list)

    pairsToPrint['pairsStacks'] = [(1,2)]

    pagetitle = "%s" % Q['name']

    htmlfilename = Q['name'].replace(" ","_")

    candidatelist = '<table id="instances" style="white-space:nowrap;">\n'

    numPositions = 2

    sequence_column = 3

    # write header line, make columns sortable, numeric=true, alpha=false
    candidatelist += '<tr><th onclick="sortTable(0,\'instances\',\'numeric\')">S.</th><th onclick="sortTable(1,\'instances\',\'checkbox\')">Show</th>'

    sequence_column += 1
    for j in range(0,numPositions):
        candidatelist += '<th onclick="sortTable(%d,\'instances\',\'alpha\')">Position %d</th>' % (j+2,j+1)
        sequence_column += 1
    candidatelist += '<th onclick="sortTable(4,\'instances\',\'alpha\')">Python</th>'
    candidatelist += '<th onclick="sortTable(5,\'instances\',\'alpha\')">Subcat</th>'
    candidatelist += '<th onclick="sortTable(6,\'instances\',\'alpha\')">New Python</th>'
    candidatelist += '<th onclick="sortTable(7,\'instances\',\'alpha\')">Matlab</th>'
    candidatelist += '<th onclick="sortTable(8,\'instances\',\'numeric\')">x</th>'
    candidatelist += '<th onclick="sortTable(9,\'instances\',\'numeric\')">y</th>'
    candidatelist += '<th onclick="sortTable(10,\'instances\',\'numeric\')">z</th>'
    candidatelist += '<th onclick="sortTable(11,\'instances\',\'numeric\')">gap12</th>'
    candidatelist += '<th onclick="sortTable(12,\'instances\',\'numeric\')">angle_in_plane</th>'
    candidatelist += '<th onclick="sortTable(13,\'instances\',\'numeric\')">normal_z</th>'
    candidatelist += '<th onclick="sortTable(14,\'instances\',\'numeric\')">max_distance</th>'
    candidatelist += '<th onclick="sortTable(15,\'instances\',\'numeric\')">min_angle</th>'
    candidatelist += '<th onclick="sortTable(16,\'instances\',\'numeric\')">max_badness</th>'

    candidatelist += "</tr>\n"

    # write one row for each candidate
    for i in range(0,len(candidates)):
        candidate = candidates[i]
        candidatelist += '<tr><td>'+str(i+1)+'.</td><td><label><input type="checkbox" id="'+str(i)+'" class="jmolInline" data-coord="'
        for j in range(0,numPositions):
            candidatelist += candidate[j]
            if j < numPositions-1:
                candidatelist += ','
        candidatelist += '">&nbsp</label></td>'

        PDB_id = candidate[0][0:4]

        """
        if PDB_id in Q["PDB_data_file"]:
            candidatelist += "<td>%s</td>" % format_resolution(Q["PDB_data_file"][PDB_id])
        else:
            candidatelist += "<td>NA</td>"
        """

        # write unit ids
        for j in range(0,numPositions):
            candidatelist += "<td>"+candidate[j]+"</td>"

        # write Python, Matlab interactions
        candidatelist += "<td>%s</td>" % candidate[11]  # python
        candidatelist += "<td>%d</td>" % candidate[12]  # python subcategory
        candidatelist += "<td>%s</td>" % candidate[13]  # new python
        candidatelist += "<td>%s</td>" % candidate[14]  # matlab

        for colnum in range(2,11):
            candidatelist += "<td>%0.2f</td>" % candidate[colnum]

        candidatelist += '</tr>\n'
    candidatelist += '</table>\n'

    discrepancydata = ''

    if np.size(allvsallmatrix) > 0:
        # write discrepancy data in new 2022 list format
        # first element is a reference to the div in which the heatmap should appear
        discrepancydata = '["#heatmap",['              # start a list, start a matrix

        # second element is a matrix with the numerical values of the discrepancy
        # writing both upper and lower triangles of the matrix
        s = allvsallmatrix.shape[0]
        for c in range(0,s):
            discrepancydata += '['     # start a row of the discrepancy matrix
            ife1 = candidates[c][0]
            for d in range(0,s):
                ife2 = candidates[d][0]
                discrepancydata += "%.4f" % allvsallmatrix[c][d]  # one entry
                if d < s-1:
                    discrepancydata += ','  # commas between entries in a row
                else:
                    discrepancydata += '],\n'  # end a row, newline

        discrepancydata += '],\n'           # end the matrix, continue the list

        # third element is a list of labels of instances
        discrepancydata += '['              # start list of instances
        for c in range(0,s):
            ife1 = candidates[c][0]
            discrepancydata += '"' + ife1 + '"'    # write one instance name in quotes
            if c < s-1:
                discrepancydata += ","  # commas between instances
            else:
                discrepancydata += "]]" # end list of instances, end list of data

    # read template.html into one string
    with open(TEMPLATEPATH + 'template.html', 'r') as myfile:
        template = myfile.read()

    # replace ###PAGETITLE### with pagetitle
    template = template.replace("###PAGETITLE###",pagetitle)

    sequence_column = 0   # column to set in fixed width font
    template = template.replace("###sequencecolumn###",str(sequence_column))

    queryNote = "Query name: %s.  Found %d candidates from %d of %d files in %0.0f seconds." % (Q['name'].encode('ascii','ignore'),len(candidates),Q["numFilesSearched"],len(Q["searchFiles"]),Q["elapsedCPUTime"])
    queryNote = "Query name: %s.  Found %d candidates from %d of %d files in %0.0f seconds." % (Q['name'],len(candidates),Q["numFilesSearched"],len(Q["searchFiles"]),Q["elapsedCPUTime"])

    if "moreCandidatesThanHeatMap" in Q:
        queryNote += " " + Q["moreCandidatesThanHeatMap"] + "\n"
    else:
        queryNote += "\n"

    template = template.replace("###QUERYNAME###",str(queryNote.encode('ascii','ignore')))

    seeModifyQuery = ''
    template = template.replace("###SEEMODIFYQUERY###",seeModifyQuery)

    template = template.replace("###seeCSVOutput###","")

    # replace ###CANDIDATELIST### with candidatelist
    template = template.replace("###CANDIDATELIST###",candidatelist)

    template = template.replace("###JS1###",JS1)
    template = template.replace("###JS2###",JS2)
    template = template.replace("###JS3###",JS3)
    template = template.replace("###JS4###",JS4)

    refresh = ""
    if "reloadOutputPage" in Q and Q["reloadOutputPage"]:
        refresh = '<meta http-equiv="refresh" content="%d">' % REFRESHTIME
    template = template.replace("###REFRESH###",refresh)

    if np.size(allvsallmatrix) > 0:
        template = template.replace("###JS5###",JS5)    # include heatmap.js code
        discrepancydata = "var data =  " + discrepancydata
        discrepancydata = '<script type="text/javascript">\n' + discrepancydata + '\n</script>'
        template = template.replace("###DISCREPANCYDATA###",discrepancydata)
    else:
        #template = template.replace("###DISCREPANCYDATA###","")
        template = template.replace("###JS5###","")    # do not display a heat map

    outputfilename = os.path.join(OUTPUTPATH,htmlfilename+".html")

    print("Writing to %s" % outputfilename)

    messages = ""

    messages += "\n<br>"
    if len(Q["userMessage"]) > 0:
        messages += "User messages:<br>\n"
        for line in Q["userMessage"]:
            messages += line + "<br>\n"
    else:
        messages += "No error or warning messages.<br>\n"

    template = template.replace("###MESSAGES###",messages)

    with open(outputfilename, 'w') as myfile:
        myfile.write(template)


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

    if 'hbond' in datapoint:
        for hbond in datapoint['hbond']:
            print(hbond)

    return None


def plot_basepair_cutoffs(base_combination,lowercase_list,ax,variables):

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


if __name__=="__main__":

    write_html_pages = True

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
    orange = [1,165/255,0] # orange

    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.216/2.0A/csv']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.216/3.0A/csv']
    PDB_list = ['4V9F','7K00']
    PDB_list = ['4TNA']
    PDB_list = ['4V9F','6ZMI','7K00','4TNA']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.261/2.5A/csv']
    PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.262/2.0A/csv']

    PDB_IFE_Dict = map_PDB_list_to_PDB_IFE_dict(PDB_list)

    # load all datapoints on pairs of bases, whether annotated as paired or not
    all_PDB_ids = PDB_IFE_Dict.keys()
    print("Working on %d PDB files" % len(all_PDB_ids))

    PDB_skip_set = set(['1R9F','5NXT','4KTG'])

    # load annotations of these PDB files from the BGSU RNA server
    print('Loading Matlab annotations')
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
    print("Found Matlab annotations in %s files" % len(all_PDB_ids))
    print("Loading Python annotations from %d PDB files" % len(all_PDB_ids))

    pair_to_data = defaultdict(dict)

    # load output files from NA_pairwise_interactions
    for PDB in all_PDB_ids:
        pair_to_data_file = outputNAPairwiseInteractions + "%s_pairs_v1.pickle" % PDB
        try:
            if sys.version_info[0] < 3:
                new_dict = pickle.load(open(pair_to_data_file,'rb'))
            else:
                new_dict = pickle.load(open(pair_to_data_file,'rb'),encoding = 'latin1')
            pair_to_data.update(new_dict)
            if len(new_dict.keys()) == 0:
                print("No Python-annotated pairs in %s" % PDB)
        except:
            print("Not able to load Python annotations for %s from %s" % (PDB,pair_to_data_file))

    # loop over specified sets of interactions
    for interaction_list in interaction_lists:

        lowercase_list = [i.lower() for i in interaction_list]

        # identify unit id pairs that are annotated with this basepairing by Matlab code
        Matlab_pairs = []
        for interaction in interaction_list:
            Matlab_pairs += Matlab_annotation_to_pair[interaction]
        Matlab_pairs = set(Matlab_pairs)

        # loop over specified base combinations
        for base_combination in base_combination_list:

            # don't show AA cHW because AA cWH will be shown, for example
            if base_combination in symmetric_base_combination_list and interaction in ['cHW','tHW','cSW','tSW','cSH','tSH']:
                continue

            nt1_seq, nt2_seq = base_combination.split(",")

            # accumulate data specific to this interaction and base combination
            xvalues = []    # x coordinate of the second base
            yvalues = []    # y coordinate of the second base
            rvalues = []    # radius of second base
            tvalues = []    # angle of second base relative to first
            zvalues = []    # z coordinate of the second base
            avalues = []    # rotation angle of second base
            gvalues = []    # gap
            nvalues = []    # third component of normal vector of second base
            dvalues = []    # worst hydrogen bond distance
            colors2d  = []  # store the color to use
            sizes = []      # store the size of dot to use

            pairs = []      # list of data to print in a table

            # accumulate information about hydrogen bonds; multiple ones per basepair
            hdvalues = []   # hydrogen bond distances
            havalues = []   # hydrogen bond angles
            hbvalues = []   # hydrogen bond badness measures
            hcolors2d = []  # store color to use for hydrogen bond dots
            hsizes    = []  # store size of dot to use for hydrogen bonds

            c = 0           # count points

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

                    if Matlab:
                        matlab_annotation = pair_to_Matlab_annotation[pair]
                    else:
                        matlab_annotation = ""

                    if 'basepair' in datapoint:
                        python_annotation = datapoint['basepair']
                        new_python_annotation = datapoint['basepair']
                        basepair = datapoint['basepair']
                        basepair_subcat = datapoint['basepair_subcategory']
                    else:
                        basepair = "   "
                        basepair_subcat = '0'

                    max_badness = 0
                    max_distance = 0
                    min_angle = 180
                    if 'hbond' in datapoint:
                        for hbond in datapoint['hbond']:
                            if hbond[0][0]:    # able to check the hydrogen bond
                                if hbond[1][4] > max_badness:
                                    max_badness = hbond[1][4]
                                if hbond[1][2] > max_distance:
                                    max_distance = hbond[1][2]
                                if not math.isnan(hbond[1][3]) and hbond[1][3] < min_angle:
                                    min_angle = hbond[1][3]

                    dvalues.append(max_distance)

                    disqualified = False

                    if Python:
                        if max_distance > 4 or min_angle < 100:
                            new_python_annotation = ''
                            disqualified = True

                    pairs.append((pair[0],pair[1],datapoint['x'],datapoint['y'],datapoint['z'],datapoint['gap12'],datapoint['angle_in_plane'],datapoint['normal_Z'],max_distance,min_angle,max_badness-0.000001*random.uniform(0,1),python_annotation,int(basepair_subcat),new_python_annotation,matlab_annotation))

                    # color and print information about bad examples
                    if disqualified:
                        color = orange
                        hcolor = orange
                        size = 40
                        hsize = 40
                    elif Python and not Matlab:
                        color = blue
                        hcolor = blue
                        size = 40
                        hsize = 40
                        #print("blue = only Python:  %s Python: %4s Matlab: %4s\n%s" % (base_combination,basepair,pair_to_Matlab_annotation[pair],datapoint['url']))
                        #print_datapoint(datapoint)
                    elif not Python and Matlab:
                        color = red
                        hcolor = red
                        size = 40
                        hsize = 40
                        #print("red = only Matlab: %s Python: %4s Matlab: %4s\n%s" % (base_combination,basepair,pair_to_Matlab_annotation[pair],datapoint['url']))
                        #print_datapoint(datapoint)
                    elif False and max_badness > 2.5:
                        color = cyan
                        size = 20
                        hcolor = cyan
                        hsize = 20
                        #print("cyan max_badness %8.4f %s Python %4s Matlab %4s\n%s" % (max_badness,base_combination,basepair,pair_to_Matlab_annotation[pair],datapoint['url']))
                        #print_datapoint(datapoint)
                    elif Python and Matlab:
                        color = black
                        size = 1
                        hcolor = color
                        hsize = size
                    else:
                        color = [0,1,0]
                        size = 100
                        hcolor = color
                        hsize = size
                    colors2d.append(color)
                    sizes.append(size)

                    if 'hbond' in datapoint:
                        for hbond in datapoint['hbond']:
                            if hbond[0][0]:    # able to check the hydrogen bond
                                hdvalues.append(hbond[1][2])  # hbond distance
                                havalues.append(hbond[1][3])  # hbond angle
                                hbvalues.append(hbond[1][4])  # hbond badness
                                hcolors2d.append(hcolor)
                                hsizes.append(hsize)

            # make scatterplots for pairwise combinations of data
            if c > 0:
                fig = plt.figure(figsize=(11.0, 7.0))

                ax = fig.add_subplot(2, 3, 1)
                ax.axis("equal")
                plot_basepair_cutoffs(base_combination,lowercase_list,ax,1)
                ax.scatter(xvalues,yvalues,color=colors2d,marker=".",s=sizes)
                ax.set_title('x and y for %d %s %s' % (c,base_combination,interaction_list[0]))
                draw_base(nt1_seq,'default',2,ax)

                ax = fig.add_subplot(2, 3, 2)
                plot_basepair_cutoffs(base_combination,lowercase_list,ax,2)
                ax.scatter(tvalues,rvalues,color=colors2d,marker=".",s=sizes)
                ax.set_title('radius vs theta')

                ax = fig.add_subplot(2, 3, 3)
                plot_basepair_cutoffs(base_combination,lowercase_list,ax,3)
                ax.scatter(avalues,nvalues,color=colors2d,marker=".",s=sizes)
                ax.set_title('normal vs angle blue=P only')

                ax = fig.add_subplot(2, 3, 4)
                plot_basepair_cutoffs(base_combination,lowercase_list,ax,4)
                ax.scatter(gvalues,zvalues,color=colors2d,marker=".",s=sizes)
                ax.set_title('z versus gap12 red=M only')

                ax = fig.add_subplot(2, 3, 5)
                ax.scatter(hdvalues,havalues,color=hcolors2d,marker=".",s=hsizes)
                ax.set_title('h-angle vs h-dist')

                ax = fig.add_subplot(2, 3, 6)
                ax.scatter(hdvalues,hbvalues,color=hcolors2d,marker=".",s=hsizes)
                ax.set_title('h-badness vs h-dist')

                # show all plots for this interaction_list
                figManager = plt.get_current_fig_manager()
                figManager.full_screen_toggle()
                if len(all_PDB_ids) == 1:
                    figure_save_file = outputNAPairwiseInteractions + "basepairs_%s_%s_%s.png" % (interaction_list[0],base_combination,all_PDB_ids[0])
                else:
                    figure_save_file = outputNAPairwiseInteractions + "basepairs_%s_%s_%s.png" % (interaction_list[0],base_combination,len(all_PDB_ids))
                plt.savefig(figure_save_file)
                #plt.show()
                plt.close()

            print("Plotted %5d points for %s %s" % (c,base_combination,",".join(interaction_list)))

            if c > 0 and write_html_pages:
                # mimic how WebFR3D writes result pages
                Q = {}
                Q['name'] = "%s %s %d" % (",".join(interaction_list),base_combination,len(all_PDB_ids))
                Q['numFilesSearched'] = len(all_PDB_ids)
                Q['searchFiles'] = all_PDB_ids
                Q['elapsedCPUTime'] = 0
                Q['userMessage'] = []

                # if too many pairs, show some good and also the worst ones
                if len(xvalues) > 300:
                    pairs = sorted(pairs, key=lambda p: p[9])
                    orderpairs = pairs[0:50] + pairs[-250:]
                    otherpairs = sorted(pairs[51:-251], key=lambda p: (p[11],p[9]))
                else:
                    orderpairs = pairs
                    otherpairs = []

                n = len(orderpairs)  # number of pairs to order by similarity

                dista = np.zeros((n,n))  # matrix to display
                distb = np.zeros((n,n))  # matrix to order by
                maxd = 0
                for i in range(0,n):
                    for j in range(i+1,n):
                        d = 0.0
                        d += (orderpairs[i][2]-orderpairs[j][2])**2
                        d += (orderpairs[i][3]-orderpairs[j][3])**2
                        d += (orderpairs[i][4]-orderpairs[j][4])**2
                        d += (orderpairs[i][5]-orderpairs[j][5])**2
                        d += 0.01*(min(abs(orderpairs[i][6]-orderpairs[j][6]),360-abs(orderpairs[i][6]-orderpairs[j][6])))**2
                        d += (orderpairs[i][7]-orderpairs[j][7])**2
                        dista[i][j] = math.sqrt(d)
                        dista[j][i] = dista[i][j]
                        maxd = max(maxd,dista[i][j])

                        d += 100*(orderpairs[i][12]-orderpairs[j][12])**2
                        distb[i][j] = math.sqrt(d)
                        distb[j][i] = math.sqrt(d)

                # color entries on diagonal according to matching annotations
                for i in range(0,n):
                    # python annotation has changed by new rules
                    if not orderpairs[i][11] == orderpairs[i][13]:
                        dista[i][i] = -2
                    # new python annotation differs from Matlab annotation
                    if not orderpairs[i][13].lower() == orderpairs[i][14].lower():
                        dista[i][i] = -i

                print("Finding order")
                order = treePenalizedPathLength(distb,20)
                print("Found order")
                #print(order)

                reorder_pairs = [orderpairs[o] for o in order] + otherpairs

                reorder_dista = np.zeros((n,n))
                for i in range(0,n):
                    for j in range(0,n):
                        reorder_dista[i][j] = dista[order[i]][order[j]]

                print("Reordered instances and distance matrix")

                writeHTMLOutput(Q,reorder_pairs,reorder_dista)

                print("Wrote efficient HTML file")

                """ Write the list of candidates in an HTML format that also shows
                the coordinate window and a heat map of all-against-all distances.
                """
                # create table entries
                # create heat map



    print("Saved figures in %s" % outputNAPairwiseInteractions)
