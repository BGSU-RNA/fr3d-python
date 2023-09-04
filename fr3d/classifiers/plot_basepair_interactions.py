"""
    plot-basepair-interactions.py reads a data file and plots points to represent
    base pairing interactions.
    Points are colored according to python_fr3d and Matlab annotations.

    HTML output is in C:/Users/zirbel/Documents/FR3D/Python FR3D/output
    PNG output is in C:/Users/zirbel/Documents/FR3D/NAPairwiseInteractions

    cd "C:/Users/zirbel/Documents/FR3D/Python FR3D/output"
    python -m SimpleHTTPServer

"""

from collections import defaultdict
import json
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import pickle
import random
import requests
import sys
import time

# import the version of urlretrieve appropriate to the Python version
if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
    from urllib import urlopen
else:
    from urllib.request import urlretrieve as urlretrieve
    from urllib.request import urlopen

from class_limits_2023 import nt_nt_cutoffs
from NA_pairwise_interactions import map_PDB_list_to_PDB_IFE_dict
from NA_pairwise_interactions import reverse_edges
from draw_residues import draw_base

from fr3d.localpath import outputText
from fr3d.localpath import outputNAPairwiseInteractions
from fr3d.localpath import fr3d_pickle_path
from fr3d.localpath import contact_list_file
from fr3d.localpath import inputPath
from fr3d.localpath import outputHTML
from fr3d.localpath import storeMatlabFR3DPairs

from fr3d.search.file_reading import readPDBDatafile
from fr3d.search.write_output import format_resolution

# very specific directories, for data files in specific formats.  Hard for others to run!
dssr_basepair_path = 'C:/Users/zirbel/Documents/FR3D/Python FR3D/data/pairs_dssr'
rnaview_basepair_path = 'C:/Users/zirbel/Documents/FR3D/Python FR3D/data/pairs_rnaview'
pdb_basepair_path = 'C:/Users/zirbel/Documents/FR3D/Python FR3D/data/pairs_pdb'

OUTPUTPATH = "C:/Users/zirbel/Documents/FR3D/Python FR3D/output/"

#Updated modified nucleotide mappings from atom_mappings_refined.txt
from fr3d.data.mapping import modified_base_atom_list,parent_atom_to_modified,modified_atom_to_parent,modified_base_to_parent

from fr3d.ordering.orderBySimilarity import treePenalizedPathLength
from fr3d.ordering.orderBySimilarity import standardOrder

JS1 = '  <script src="./js/JSmol.min.nojq.js"></script>'
JS2 = '  <script src="./js/jquery.jmolTools.bp.js"></script>'               # special version, superimpose first base
JS3 = '  <script src="./js/imagehandlinglocal.js"></script>'
JS4 = '<script src="./js/jmolplugin.js" type="text/javascript"></script>'
JS5 = '<script type="text/javascript" src="./js/heatmap.js"></script>'
TEMPLATEPATH = '../search/'

def print_dictionary(x):
    for key, value in sorted(x.items()):
        print("  %s: %s" % (key,value))


def writeHTMLOutput(Q,candidates,allvsallmatrix=np.empty( shape=(0, 0) )):
    """
    Write the list of candidates in an HTML format that also shows
    the coordinate window and a heat map of all-against-all distances.
    Write the heatmap data in an efficient way.
    """

    pagetitle = "%s" % Q['name']

    htmlfilename = Q['name'].replace(" ","_")

    # make the table of instances
    candidatelist = '<table id="instances" style="white-space:nowrap;">\n'

    numPositions = 2

    # write header line, make columns sortable, numeric=true, alpha=false
    candidatelist += '<tr><th onclick="sortTable(0,\'instances\',\'numeric\')">S.</th>'
    candidatelist += '<th onclick="sortTable(1,\'instances\',\'checkbox\')">Show</th>'
    candidatelist += '<th onclick="sortTable(2,\'instances\',\'checkbox\')">Note</th>'
    candidatelist += '<th onclick="sortTable(3,\'instances\',\'numeric\')">Res.</th>'

    sequence_column = 3
    for j in range(0,numPositions):
        candidatelist += '<th onclick="sortTable(%d,\'instances\',\'alpha\')">Position %d</th>' % (j+4,j+1)
        sequence_column += 1

    candidatelist += '<th onclick="sortTable(6,\'instances\',\'alpha\')">RNAview</th>'
    candidatelist += '<th onclick="sortTable(7,\'instances\',\'alpha\')">DSSR</th>'
    candidatelist += '<th onclick="sortTable(8,\'instances\',\'alpha\')">pdb</th>'
    candidatelist += '<th onclick="sortTable(9,\'instances\',\'alpha\')">contacts</th>'
    candidatelist += '<th onclick="sortTable(10,\'instances\',\'alpha\')">FR3D</th>'
    candidatelist += '<th onclick="sortTable(11,\'instances\',\'alpha\')">New FR3D</th>'
    candidatelist += '<th onclick="sortTable(12,\'instances\',\'alpha\')">Newest FR3D</th>'
    candidatelist += '<th onclick="sortTable(13,\'instances\',\'alpha\')">Subcat</th>'
    candidatelist += '<th onclick="sortTable(14,\'instances\',\'numeric\')">cut dist</th>'
    candidatelist += '<th onclick="sortTable(15,\'instances\',\'numeric\')">x</th>'
    candidatelist += '<th onclick="sortTable(16,\'instances\',\'numeric\')">y</th>'
    candidatelist += '<th onclick="sortTable(17,\'instances\',\'numeric\')">z</th>'
    candidatelist += '<th onclick="sortTable(18,\'instances\',\'numeric\')">r</th>'
    candidatelist += '<th onclick="sortTable(19,\'instances\',\'numeric\')">maxgap</th>'
    candidatelist += '<th onclick="sortTable(20,\'instances\',\'numeric\')">angle</th>'
    candidatelist += '<th onclick="sortTable(21,\'instances\',\'numeric\')">mindist</th>'
    candidatelist += '<th onclick="sortTable(22,\'instances\',\'numeric\')">normal_z</th>'
    candidatelist += '<th onclick="sortTable(23,\'instances\',\'numeric\')">h dist</th>'
    candidatelist += '<th onclick="sortTable(24,\'instances\',\'numeric\')">h ang</th>'
    candidatelist += '<th onclick="sortTable(25,\'instances\',\'numeric\')">h bad</th>'

    if 'css' in Q['name'].lower():
        candidatelist += '<th onclick="sortTable(26,\'instances\',\'alpha\')">SR12</th>'
        candidatelist += '<th onclick="sortTable(27,\'instances\',\'alpha\')">SR21</th>'

    candidatelist += "</tr>\n"

    # write one row for each candidate
    for i in range(0,len(candidates)):
        candidate = candidates[i]
        candidatelist += '<tr><td>'+str(i+1)+'.</td>'
        candidatelist += '<td><label><input type="checkbox" id="'+str(i)+'" class="jmolInline" data-coord="'
        candidatelist += candidate['unit_id_1'] + ',' + candidate['unit_id_2']
        candidatelist += '">&nbsp</label></td>'

        candidatelist += '<td><label><input type="checkbox" id="note_'+str(i)+'">&nbsp</label></td>'

        PDB_id = candidate['unit_id_1'][0:4]

        if PDB_id in Q["PDB_data_file"]:
            candidatelist += "<td>%s</td>" % format_resolution(Q["PDB_data_file"][PDB_id])
        else:
            candidatelist += "<td>NA</td>"

        # write unit ids
        candidatelist += "<td>"+candidate['unit_id_1']+"</td>"
        candidatelist += "<td>"+candidate['unit_id_2']+"</td>"

        # write python_fr3d, Matlab interactions
        candidatelist += "<td>%s</td>" % candidate['rnaview_annotation']  # rnaview
        candidatelist += "<td>%s</td>" % candidate['dssr_annotation']  # dssr
        candidatelist += "<td>%s</td>" % candidate['pdb_annotation']  # pdb
        candidatelist += "<td>%s</td>" % candidate['contacts_annotation']  # contacts
        candidatelist += "<td>%s</td>" % candidate['matlab_annotation']  # matlab
        candidatelist += "<td>%s</td>" % candidate['python_annotation']  # python
        candidatelist += "<td>%s</td>" % candidate['new_python_annotation']  # new python
        candidatelist += "<td>%d</td>" % candidate['basepair_subcat']  # python subcategory
        candidatelist += "<td>%0.4f</td>" % candidate['best_cutoff_distance']  #
        candidatelist += "<td>%0.2f</td>" % candidate['x']  #
        candidatelist += "<td>%0.2f</td>" % candidate['y']  #
        candidatelist += "<td>%0.2f</td>" % candidate['z']  #
        candidatelist += "<td>%0.2f</td>" % (math.sqrt(candidate['x']**2 + candidate['y']**2))
        candidatelist += "<td>%0.2f</td>" % candidate['maxgap']  #
        candidatelist += "<td>%0.2f</td>" % candidate['angle_in_plane']  #

        if not 'min_distance' in candidate:
            candidate['min_distance'] = -1.0

        candidatelist += "<td>%0.2f</td>" % candidate['min_distance']  #
        candidatelist += "<td>%0.2f</td>" % candidate['normal_Z']  #
        candidatelist += "<td>%0.2f</td>" % candidate['max_distance']  #
        candidatelist += "<td>%0.2f</td>" % candidate['min_angle']  #
        candidatelist += "<td>%0.2f</td>" % candidate['max_badness']  #

        if 'css' in Q['name'].lower():
            if 'sugar_ribose' in candidate:
                candidatelist += "<td>%s</td>" % candidate['sugar_ribose']
            else:
                candidatelist += "<td></td>"

            if 'r_sugar_ribose' in candidate:
                candidatelist += "<td>%s</td>" % candidate['r_sugar_ribose']
            else:
                candidatelist += "<td></td>"

        if 'hbond_messages' in candidate:
            for message in candidate['hbond_messages']:
                # write messages pertaining to the annotated basepair
                candidatelist += "<td>%s</td>" % message   # h-bond message
        else:
            candidatelist += "<td>No hydrogen bond information</td>"

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
            #ife1 = candidates[c]['unit_id_1']
            for d in range(0,s):
                #ife2 = candidates[d]['unit_id_2']
                discrepancydata += "%.4f" % allvsallmatrix[c][d]  # one entry
                if d < s-1:
                    discrepancydata += ','  # commas between entries in a row
                else:
                    discrepancydata += '],\n'  # end a row, newline

        discrepancydata += '],\n'           # end the matrix, continue the list

        # third element is a list of labels of instances
        discrepancydata += '['              # start list of instances
        for c in range(0,s):
            ife1 = candidates[c]['unit_id_1']
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

    sequence_column = 0   # column to set in fixed width font, not important here!
    template = template.replace("###sequencecolumn###",str(sequence_column))

    queryNote = "%s.  Found %d candidates from %d files, showing %d." % (Q['name'],Q['numFound'],len(Q["searchFiles"]),len(candidates))

    template = template.replace("###QUERYNAME###",queryNote)

    seeModifyQuery = ''
    template = template.replace("###SEEMODIFYQUERY### |",seeModifyQuery)

    prefix = ""
    if Q['DNA']:
        prefix += "DNA_"
    if Q['only_modified']:
        prefix += "modified_"

    table, combination_to_LSW_status = generate_LW_family_table(Q['LW'],prefix=prefix)

    template = template.replace("###seeCSVOutput###",table)

    description = "<br>Columns of the table show candidate number in similarity order, checkbox to display coordinates or not, interacting nucleotides, annotations from different programs, variables describing the interaction, and hydrogen bond information<br>"
    # put in the scatterplot
    description += '\n<br><img src="%s"><br>\n' % Q['figure_img_src']

    template = template.replace("###DESCRIPTION###",description)

    # replace ###CANDIDATELIST### with candidatelist
    template = template.replace("###CANDIDATELIST###",candidatelist)

    if Q['DNA']:
        # for when working with DNA
        JS2 = '  <script src="./js/jquery.jmolTools.rnatest.js"></script>'  # special version to get coordinates
    else:
        JS2 = '  <script src="./js/jquery.jmolTools.bp.js"></script>'  # special version for basepairs

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
    if len(Q["userMessage"]) > 0:
        for line in Q["userMessage"]:
            messages += line + "<br>\n"
    else:
        #messages += "No error or warning messages.<br>\n"
        pass

    template = template.replace("###MESSAGES###",messages)

    with open(outputfilename, 'w') as myfile:
        myfile.write(template)


def writeHTMLOutputCompareAnnotators(Q,candidates,allvsallmatrix=np.empty( shape=(0, 0) )):
    """
    Write the list of candidates in an HTML format that also shows
    the coordinate window and a heat map of all-against-all distances.
    Write the heatmap data in an efficient way.
    Do not show x,y,z parameters or scatterplots, but do show annotations from different programs.
    """

    pagetitle = "%s" % Q['name']

    htmlfilename = "compare_v6_" + Q['name'].replace(" ","_")
    htmlfilename = htmlfilename.replace('compare_v6_DNA','DNA_compare_v6')        # awkward, fragile, but works

    print('htmlfilename',htmlfilename)
    print('Qname',Q['name'])

    base_combination = Q['name'].split()[1]

    # make the table of instances
    candidatelist = '<table id="instances" style="white-space:nowrap;">\n'

    numPositions = 2

    # write header line, make columns sortable, numeric=true, alpha=false
    candidatelist += '<tr><th onclick="sortTable(0,\'instances\',\'numeric\')">S.</th><th onclick="sortTable(1,\'instances\',\'checkbox\')">Show</th>'
    candidatelist += '<th onclick="sortTable(2,\'instances\',\'numeric\')">Res.</th>'

    for j in range(0,numPositions):
        candidatelist += '<th onclick="sortTable(%d,\'instances\',\'alpha\')">Position %d</th>' % (j+3,j+1)
    candidatelist += '<th onclick="sortTable(5,\'instances\',\'alpha\')">FR3D</th>'
    candidatelist += '<th onclick="sortTable(6,\'instances\',\'alpha\')">RNAview</th>'
    candidatelist += '<th onclick="sortTable(7,\'instances\',\'alpha\')">DSSR</th>'
    candidatelist += '<th onclick="sortTable(8,\'instances\',\'alpha\')">PDB</th>'
    candidatelist += '<th onclick="sortTable(9,\'instances\',\'alpha\')">contacts</th>'
    candidatelist += '<th onclick="sortTable(10,\'instances\',\'numeric\')">Number</th>'
    candidatelist += '<th onclick="sortTable(11,\'instances\',\'alpha\')">Which</th>'
    candidatelist += '<th onclick="sortTable(12,\'instances\',\'alpha\')">Glyco 1</th>'
    candidatelist += '<th onclick="sortTable(13,\'instances\',\'alpha\')">Glyco 2</th>'

    current_column = 14   # next column after the one list above, keep track manually

    # loop over the hydrogen bonds for this base combination and basepair category
    atom_sets = []
    for candidate in candidates:
        if "atom_sets" in candidate and len(candidate["atom_sets"]) > 0:
            atom_sets = sorted(candidate["atom_sets"])
            #print("Found atom sets %s from %s" % (atom_sets,candidate))
            print("Found atom sets %s" % (atom_sets))

            for atom_set in atom_sets:
                result = candidate["atom_set_to_results"][atom_set]
                atoms = result["donor_acceptor_atoms"]
                candidatelist += '<th onclick="sortTable(%d,\'instances\',\'numeric\')">%s</th>' % (current_column,atoms)
                current_column += 1
                atoms = result["heavy_donor_acceptor_atoms"]
                candidatelist += '<th onclick="sortTable(%d,\'instances\',\'numeric\')">%s</th>' % (current_column,atoms)
                current_column += 1

                print(result)

            candidatelist += '<th onclick="sortTable(%d,\'instances\',\'numeric\')">SBHBB</th>' % (current_column)  # second best hydrogen bond badness
            break

    # pdata["atom_sets"] = sorted(datapoint["LW_to_atom_sets"][interaction])
    # pdata["atom_set_to_results"] = datapoint["atom_set_to_results"]

    if not Q['DNA']:
       candidatelist += '<th>3D structures</th>'
       candidatelist += '<th>Sequence alignments</th>'

    candidatelist += "</tr>\n"

    first_listed_base = 'A'
    for i in range(0,len(candidates)):
        candidate = candidates[i]
        first_listed_base = candidate['unit_id_1'].split("|")[3]
        if first_listed_base in ['A','C','G','U','DA','DC','DG','DT']:
            break

    print('First listed base is %s, base_combination is %s' % (first_listed_base,base_combination))

    # list bases in the order they appear in the filename
    if first_listed_base == base_combination.split(",")[0]:
        id_1 = 'unit_id_1'
        id_2 = 'unit_id_2'
    else:
        id_1 = 'unit_id_2'
        id_2 = 'unit_id_1'

    # write one row for each candidate
    for i in range(0,len(candidates)):
        candidate = candidates[i]
        candidatelist += '<tr><td>'+str(i+1)+'.</td><td><label><input type="checkbox" id="'+str(i)+'" class="jmolInline" data-coord="'

        candidatelist += candidate[id_1] + ',' + candidate[id_2]

        candidatelist += '">&nbsp</label></td>'

        PDB_id = candidate['unit_id_1'][0:4]

        if PDB_id in Q["PDB_data_file"]:
            candidatelist += "<td>%s</td>" % format_resolution(Q["PDB_data_file"][PDB_id])
        else:
            candidatelist += "<td>NA</td>"

        # write unit ids
        candidatelist += "<td>"+candidate[id_1]+"</td>"
        candidatelist += "<td>"+candidate[id_2]+"</td>"

        if id_1 == 'unit_id_1':
            candidatelist += "<td>%s</td>" % candidate['python_annotation'].replace("a","")  # python FR3D, with w, h, s
            candidatelist += "<td>%s</td>" % candidate['rnaview_annotation']  # rnaview
            candidatelist += "<td>%s</td>" % candidate['dssr_annotation']  # dssr
            candidatelist += "<td>%s</td>" % candidate['pdb_annotation']  # pdb
            candidatelist += "<td>%s</td>" % candidate['contacts_annotation']  # contacts
        else:
            candidatelist += "<td>%s</td>" % reverse_edges(candidate['python_annotation'].replace("a",""))  # python FR3D, with w, h, s
            candidatelist += "<td>%s</td>" % reverse_edges(candidate['rnaview_annotation'])  # rnaview
            candidatelist += "<td>%s</td>" % reverse_edges(candidate['dssr_annotation'])  # dssr
            candidatelist += "<td>%s</td>" % candidate['pdb_annotation']  # pdb
            candidatelist += "<td>%s</td>" % reverse_edges(candidate['contacts_annotation'])  # contacts

        candidatelist += "<td>%s</td>" % candidate['annotator_count']  # number that agree
        candidatelist += "<td>%s</td>" % candidate['annotator']  # who annotates

        candidatelist += "<td>%s</td>" % candidate['glycosidic1']  # unit 1 glycosidic bond
        candidatelist += "<td>%s</td>" % candidate['glycosidic2']  # unit 2 glycosidic bond

        badness = []
        for atom_set in atom_sets:
            if "atom_set_to_results" in candidate and atom_set in candidate["atom_set_to_results"]:
                result = candidate["atom_set_to_results"][atom_set]
                candidatelist += "<td>%0.2f</td>" % (result["donor_acceptor_distance"])  # hydrogen bond distance
                candidatelist += "<td>%0.1f</td>" % (result["heavy_donor_acceptor_angle"])  # hydrogen bond angle
                bad = max(0.0, abs(result["donor_acceptor_distance"]-2.9) - 0.1)    # how far from 2.9, penalize when above 0.2
                bad += max(0.0, abs(result["heavy_donor_acceptor_angle"] - 120) - 10)/20.0  # how far from 120, penalize when more than 20 degrees away
                badness.append(bad)
            else:
                candidatelist += "<td></td><td></td>"  # no hydrogen bond information
                badness.append(9.99)

        if len(badness) > 1:
            badness = sorted(badness)
            candidatelist += "<td>%0.2f</td>" % (badness[1])  # second worst badness
        elif len(badness) == 1:
            candidatelist += "<td>%0.2f</td>" % (badness[0])  # the only badness
        else:
            candidatelist += "<td>9.99</td>"  # pretty bad

        if not Q['DNA']:
            candidatelist += '<td><a href="http://rna.bgsu.edu/correspondence/comparison?selection=%s,%s&exp_method=all&resolution=3.0" target="_blank" rel="noopener noreferrer">Compare 3D</a></td>' % (candidate[id_1],candidate[id_2])
            candidatelist += '<td><a href="http://rna.bgsu.edu/correspondence/variability?id=%s,%s&format=unique" target="_blank" rel="noopener noreferrer">Compare sequence</a></td>' % (candidate[id_1],candidate[id_2])

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
            #ife1 = candidates[c]['unit_id_1']
            for d in range(0,s):
                #ife2 = candidates[d]['unit_id_2']
                discrepancydata += "%.4f" % allvsallmatrix[c][d]  # one entry
                if d < s-1:
                    discrepancydata += ','  # commas between entries in a row
                else:
                    discrepancydata += '],\n'  # end a row, newline

        discrepancydata += '],\n'           # end the matrix, continue the list

        # third element is a list of labels of instances
        discrepancydata += '['              # start list of instances
        for c in range(0,s):
            ife1 = candidates[c][id_1]
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

    queryNote = "%s.  Found %d candidates from %d files, showing %d." % (Q['name'],Q['numFound'],len(Q["searchFiles"]),len(candidates))

    template = template.replace("<h2>FR3D output</h2>","<h2>Basepair annotation comparison</h2>")
    template = template.replace("###QUERYNAME###",queryNote)

    seeModifyQuery = ''
    template = template.replace("###SEEMODIFYQUERY### |",seeModifyQuery)

    height = 100+len(candidates)*20
    if height < 300:
        template = template.replace("height:300px","height:%dpx" % (height))

    if Q['DNA']:
        table, combination_to_LSW_status = generate_LW_family_table(Q['LW'],prefix='DNA_compare_v6_')
    else:
        table, combination_to_LSW_status = generate_LW_family_table(Q['LW'],prefix='compare_v6_')

    template = template.replace("###seeCSVOutput###",table)

    description = "<br>Columns of the table show candidate number in similarity order, checkbox to display coordinates or not, interacting nucleotides, annotations from different programs, number of programs which annotate, and a string that tells which programs annotate.  Click on column headings to sort.<br>"

    template = template.replace("###DESCRIPTION###",description)

    # replace ###CANDIDATELIST### with candidatelist
    template = template.replace("###CANDIDATELIST###",candidatelist)

    if Q['DNA']:
        # for when working with DNA
        JS2 = '  <script src="./js/jquery.jmolTools.rnatest.js"></script>'  # special version to get coordinates
    else:
        JS2 = '  <script src="./js/jquery.jmolTools.bp.js"></script>'  # special version for basepairs

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

    print("Compare annotators is writing to %s" % outputfilename)

    messages = ""

    if len(Q["userMessage"]) > 0:
        for line in Q["userMessage"]:
            messages += line + "<br>\n"
    else:
        #messages += "No error or warning messages.<br>\n"
        pass

    template = template.replace("###MESSAGES###",messages)

    with open(outputfilename, 'w') as myfile:
        myfile.write(template)


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
            ## New python_fr3d annotations use upper and lowercase to indicate
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
            and 'gap12' in datapoint and 'gap21' in datapoint and 'angle_in_plane' in datapoint and 'normal_Z' in datapoint)


def make_pretend_datapoint(datapoint):

    if not 'x' in datapoint:
        #print('Adding x to %s' % datapoint)
        datapoint['x'] = 0.000
        datapoint['basepair'] = 'N/A'
        datapoint['basepair_subcategory'] = None
    if not 'y' in datapoint:
        #print('Adding y to %s' % datapoint)
        datapoint['y'] = 0.000
        datapoint['basepair'] = 'N/A'
        datapoint['basepair_subcategory'] = None
    if not 'z' in datapoint:
        #print('Adding z to %s' % datapoint)
        datapoint['z'] = 0.000
        datapoint['basepair'] = 'N/A'
        datapoint['basepair_subcategory'] = None
    if not 'gap12' in datapoint:
        #print('Adding gap12 to %s' % datapoint)
        datapoint['gap12'] = 0.0
    if not 'gap21' in datapoint:
        #print('Adding gap21 to %s' % datapoint)
        datapoint['gap21'] = 0.0
    if not 'normal_Z' in datapoint:
        #print('Adding normal_Z to %s' % datapoint)
        datapoint['normal_Z'] = 0.0
        datapoint['basepair'] = 'N/A'
        datapoint['basepair_subcategory'] = None
    if not 'min_distance' in datapoint:
        datapoint['min_distance'] = -1.0

    # There is no value we can set that does not mess up the scatterplot
    # if not 'angle_in_plane' in datapoint:
    #     #print('Adding angle_in_plane to %s' % datapoint)
    #     datapoint['angle_in_plane'] = None  # otherwise it really messes up the scale

    return datapoint


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


def plot_basepair_cutoffs(base_combination,interaction_list,ax,variables,angle_out_of_order=False):
    """
    Plot the cutoffs in nt_nt_cutoffs as rectangles
    """

    #print('Plotting cutoffs for %s %s' % (base_combination,interaction_list))

    color = ["#BA55D3","#63B8FF","#00EE76","#FF8C00","#CDC9A5","#8A8A8A","#8A8A8A","#8A8A8A","#8A8A8A","#8A8A8A","#8A8A8A"]
    cc = 0

    for bc in nt_nt_cutoffs.keys():
        if bc == base_combination:
            for interaction in nt_nt_cutoffs[bc].keys():

                LW = interaction.replace("n","").replace("a","")

                # if variables == 1:
                #     print('  Checking %s against %s' % (LW,interaction_list))

                if LW in interaction_list:

                    # if variables == 1:
                    #     print('  Success!')

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

                            if t1 < -90:
                                t1 += 360
                            if t2 < -90:
                                t2 += 360
                            if t3 < -90:
                                t3 += 360
                            if t4 < -90:
                                t4 += 360

                            r1 = math.sqrt(xmax**2+ymin**2)
                            r2 = math.sqrt(xmax**2+ymax**2)
                            r3 = math.sqrt(xmin**2+ymax**2)
                            r4 = math.sqrt(xmin**2+ymin**2)

                            ax.plot([t1,t2,t3,t4,t1],[r1,r2,r3,r4,r1],color[cc])

                            if "radiusmax" in limits:
                                ax.plot([min([t1,t2,t3,t4]),max([t1,t2,t3,t4])],[limits["radiusmax"],limits["radiusmax"]],color[cc])

                            #ax.text(0,cc*0.7,"%s %d" % (LW,subcategory),color = color[cc],fontsize=12)
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

                            if angle_out_of_order and xmin < 0:
                                ax.plot([xmin+360,xmax+360,xmax+360,xmin+360,xmin+360],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            elif angle_out_of_order and xmax < 0:
                                ax.plot([xmin+360,xmax,xmax,xmin+360,xmin+360],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            elif xmin < xmax:
                                # angle does not wrap around 270 degrees
                                ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            else:
                                # angle wraps around 270 degrees, plot from maybe 240 to 300
                                # xmin may be 240, xmax may be -50
                                ax.plot([xmin,xmax+360,xmax+360,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                                #ax.plot([270,xmin,xmin,270],[ymin,ymin,ymax,ymax],color[cc])
                                #ax.plot([-90,xmax,xmax,-90],[ymin,ymin,ymax,ymax],color[cc])

                            cc += 1
                        if variables == 4:            # gap and z
                            xmin = -0.01
                            xmax = limits['gapmax']
                            ymin = limits['zmin'] - cc*0.04
                            ymax = limits['zmax'] + cc*0.04
                            ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            cc += 1

                        if variables == 5:            # angle and x
                            xmin = limits['anglemin']
                            xmax = limits['anglemax']
                            ymin = limits['xmin']
                            ymax = limits['xmax']

                            if angle_out_of_order and xmin < 0:
                                ax.plot([xmin+360,xmax+360,xmax+360,xmin+360,xmin+360],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            elif angle_out_of_order and xmax < 0:
                                ax.plot([xmin+360,xmax,xmax,xmin+360,xmin+360],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            elif xmin < xmax:
                                # angle does not wrap around 270 degrees
                                ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            else:
                                # angle wraps around 270 degrees
                                ax.plot([xmin,xmax+360,xmax+360,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                                #ax.plot([270,xmin,xmin,270],[ymin,ymin,ymax,ymax],color[cc])
                                #ax.plot([-90,xmax,xmax,-90],[ymin,ymin,ymax,ymax],color[cc])

                            cc += 1

                        if variables == 6:            # angle and y
                            xmin = limits['anglemin']
                            xmax = limits['anglemax']
                            ymin = limits['ymin']
                            ymax = limits['ymax']

                            if angle_out_of_order and xmin < 0:
                                ax.plot([xmin+360,xmax+360,xmax+360,xmin+360,xmin+360],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            elif angle_out_of_order and xmax < 0:
                                ax.plot([xmin+360,xmax,xmax,xmin+360,xmin+360],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            elif xmin < xmax:
                                # angle does not wrap around 270 degrees
                                ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            else:
                                # angle wraps around 270 degrees
                                ax.plot([xmin,xmax+360,xmax+360,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                                # ax.plot([270,xmin,xmin,270],[ymin,ymin,ymax,ymax],color[cc])
                                # ax.plot([-90,xmax,xmax,-90],[ymin,ymin,ymax,ymax],color[cc])

                            cc += 1


def generate_LW_family_table(LW,prefix=''):
    """
    Make an HTML table of links to pages about all base combinations in a family.
    Also return a dictionary mapping family and base combination to a Boolean telling if that
    one is present in the LSW 2002 paper
    """

    LW_clean = LW.replace("a","").replace("n","")
    LW_lower = LW_clean.lower()

    if LW_lower in ['chw','thw','csw','tsw','csh','tsh']:
        # put edges in the usual order for these families, to make the table
        LW_lower = LW_lower[0] + LW_lower[2] + LW_lower[1]

    LW_upper = LW_lower[0] + LW_lower[1:3].upper()

    output = '<table>\n'
    combination_to_LSW_status = {}

    if DNA:
        bases = ['DA','DC','DG','DT']
        resolutions = ['2.0A']           # number of structures for now, part of the hyperlinks!
        output += "<tr><td>2.0A</td><td></td><td></td><td></td><td></td></tr>\n"
    else:
        bases = ['A','C','G','U']
        resolutions = ['1.5A','2.0A','2.5A','3.0A']
        output += "<tr><td>1.5A</td><td></td><td></td><td></td><td></td><td>2.0A</td><td></td><td></td><td></td><td></td><td>2.5A</td><td></td><td></td><td></td><td></td><td>3.0A</td><td></td><td></td><td></td></tr>\n"

    for b1 in bases:
        output += "<tr>\n"
        for resolution in resolutions:
            for b2 in bases:
                base_combination = b1 + "," + b2

                status = "New"

                if base_combination in nt_nt_cutoffs:
                    for interaction in nt_nt_cutoffs[base_combination].keys():
                        # make base edges uppercase since that's all we get from the pipeline
                        interaction_clean = interaction.replace("n","").replace("a","")[0:3]
                        interaction_lower = interaction_clean.lower()

                        # if there is an interaction with this base combination ...
                        if interaction_lower == LW_lower:
                            status = "LSW 2002"

                # When on C,A cHS, look for A,C cSH because that's how it's stored but print C,A cHS
                else:

                    if LW_lower in ['cww','tww','chh','thh']:
                        status = "Symmetric"
                    else:
                        for interaction in nt_nt_cutoffs[b2 + "," + b1].keys():
                            # reverse interaction since we're looking at the lower part of the family
                            inter = interaction.replace("n","").replace("a","")[0:3]
                            interaction_reverse = inter[0] + inter[2] + inter[1]
                            interaction_lower = interaction_reverse.lower()

                            if interaction_lower == LW_lower:
                                status = "LSW 2002"

                if LW_lower == 'cww' and base_combination == 'G,G':
                    status = "New"

                if LW_lower == 'tss' and b1 in ['C','U']:
                    status = "New"

                if status == "New":
                    base_combination = base_combination.lower()

                if status == "Symmetric":
                    base_combination = "Symmetric"
                    output += '<td>Symmetric</td>\n'
                else:
                    link = prefix+"%s_%s,%s_%s.html" % (LW_upper,b1,b2,resolution)
                    output += '<td><a href="%s">%s %s</a></td>\n' % (link,base_combination,LW_upper)

                combination_to_LSW_status[b1+","+b2] = status

            if not resolution == resolutions[-1]:
                output += "<td>||</td>\n"
        output += "</tr>\n"

    output += '</table><p>\n'

    return output, combination_to_LSW_status


def load_dssr_basepairs(pdb_id):

    pickle_filename = os.path.join(dssr_basepair_path,pdb_id+'_pairs.pickle')

    interaction_to_pair = {}

    if os.path.exists(pickle_filename):
        if sys.version_info[0] < 3:
            interaction_to_pair = pickle.load(open(pickle_filename,"rb"))
        else:
            interaction_to_pair = pickle.load(open(pickle_filename,"rb"), encoding = 'latin1')

    else:
        j = 1
        keep_trying = True
        while keep_trying:

            url = 'http://west.nakb.org/x3dssr/%s_%d.json' % (pdb_id,j)

            time.sleep(0.1)
            try:
                url_data = requests.get(url)
            except:
                print('Pause to re-establish connection for %s' % pdb_id)
                time.sleep(5)
                url_data = requests.get(url)

            if '404 Not Found' in url_data.text:
                keep_trying = False
            else:
                print('Downloading DSSR annotations for %s on try #%d' % (pdb_id,j))
                j += 1
                dssr = json.loads(url_data.text)

                c = 0
                if 'pairs' in dssr:
                    for pair_data in dssr['pairs']:
                        LW = pair_data['LW']

                        f1 = pair_data['nt1'].split(".")   # 1..0.U.12. or 1..B.C.47.H with insertion code
                        f2 = pair_data['nt2'].split(".")

                        if f1[5]:
                            # with insertion code
                            u1 = '%s|%s|%s|%s|%s|||%s' % (pdb_id,f1[0],f1[2],f1[3],f1[4],f1[5])
                            print('DSSR unit id with insertion code %s' % u1)
                        else:
                            u1 = '%s|%s|%s|%s|%s' % (pdb_id,f1[0],f1[2],f1[3],f1[4])

                        if f2[5]:
                            # with insertion code
                            u2 = '%s|%s|%s|%s|%s|||%s' % (pdb_id,f2[0],f2[2],f2[3],f2[4],f2[5])
                            print('DSSR unit id with insertion code %s' % u2)
                        else:
                            u2 = '%s|%s|%s|%s|%s' % (pdb_id,f2[0],f2[2],f2[3],f2[4])

                        if '-' in LW or len(LW) < 3 or (not 'c' in LW and not 't' in LW):
                            #print(u1,u2,LW,'messed up',len(LW))
                            continue

                        rLW = LW[0] + LW[2] + LW[1]

                        if not LW in interaction_to_pair:
                            interaction_to_pair[LW] = []
                            interaction_to_pair[rLW] = []

                        interaction_to_pair[LW].append((u1,u2,0))       # pretend that range is 0; maybe change it some day
                        interaction_to_pair[rLW].append((u2,u1,0))

                        #print(u1,u2,LW)

        with open(pickle_filename, 'wb') as fh:
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(interaction_to_pair, fh, 2)

    #print('Not able to load DSSR basepair annotations for %s' % pdb_id)

    return interaction_to_pair


def load_rnaview_basepairs(pdb_id):

    pickle_filename = os.path.join(rnaview_basepair_path,pdb_id+'_pairs.pickle')

    interaction_to_pair = {}

    if os.path.exists(pickle_filename):
        if sys.version_info[0] < 3:
            interaction_to_pair = pickle.load(open(pickle_filename,"rb"))
        else:
            interaction_to_pair = pickle.load(open(pickle_filename,"rb"), encoding = 'latin1')

    return interaction_to_pair


def load_pdb_basepairs(pdb_id):
    """
    The program PDB_basepair_reader.py downloads and extracts basepairs from .cif files.
    It needs to be run every time we change the representative set.
    """

    pickle_filename = os.path.join(pdb_basepair_path,pdb_id+'_pairs.pickle')

    interaction_to_pair = {}

    if os.path.exists(pickle_filename):
        if sys.version_info[0] < 3:
            interaction_to_pair = pickle.load(open(pickle_filename,"rb"))
        else:
            interaction_to_pair = pickle.load(open(pickle_filename,"rb"), encoding = 'latin1')

    return interaction_to_pair


def load_contacts_basepairs(contacts_filename):

    pair_to_interaction = defaultdict(str)

    if os.path.exists(contacts_filename):
        with open(contacts_filename,'r') as f:
            for line in f:
                fields = line.split()
                u1 = fields[0].replace('||||1_555','')
                interaction = fields[1]
                u2 = fields[2].replace('||||1_555','')

                pair_to_interaction[(u1,u2)] = interaction
                pair_to_interaction[(u2,u1)] = reverse_edges(interaction)

    return pair_to_interaction


def add_pairs_in_order(pair_to_priority,pair_to_data,priority):
    """
    Assign a priority score to each pair.
    Lower priority numbers will be considered earlier.
    This attempts to improve homogeneity in each group.
    """

    for pair, datapoint in pair_to_data.items():
        # no need to consider a pair twice
        if pair in pair_to_priority:
            continue

        # delay pairs with bases in an unusual order
        bc = pair[0].split("|")[3] + pair[1].split("|")[3]
        if bc in ['CA','CG','GA','UA','UC','UG']:
            bc_addition = 0.5
        else:
            bc_addition = 0.0

        if priority == 1 and 'basepair' in datapoint:
            interaction = datapoint['basepair'].replace("n","")
            if interaction and interaction[1].islower():
                # delay cwW until after cWw
                pair_to_priority[pair] = 1.0 + bc_addition + 0.1
            else:
                pair_to_priority[pair] = 1.0 + bc_addition

        else:
            pair_to_priority[pair] = priority + bc_addition
            if bc in ['AA','CC','GG','UU']:
                rev_pair = (pair[1],pair[0])
                if not rev_pair in pair_to_priority:
                    pair_to_priority[rev_pair] = priority

    return pair_to_priority


def evaluate_pair_from_datapoint(datapoint,interaction,nt_nt_cutoffs_bc):
    """
    Using the information in the variable datapoint and the cutoffs passed in,
    find the best-matching category and how close the pair is to that category.
    If data is missing, do your best.
    """

    angle_out_of_order = False

    have_full_data = check_full_data(datapoint)

    if not have_full_data:
        datapoint = make_pretend_datapoint(datapoint)

    datapoint['maxgap'] = max(datapoint['gap21'],datapoint['gap12'])

    if 'basepair' in datapoint:
        python_annotation = datapoint['basepair']
        new_python_annotation = datapoint['basepair']
        basepair = datapoint['basepair']
        basepair_subcat = datapoint['basepair_subcategory']
    else:
        python_annotation = ''
        new_python_annotation = ''
        basepair = "   "
        basepair_subcat = '0'

    # find worst hydrogen bond
    max_badness = -1.0
    max_distance = -1.0
    min_angle = 180
    if 'hbond' in datapoint:
        for hbond in datapoint['hbond']:
            # hbond is like
            # (True, True, 2.391828053713497, 127.44856702606579, 1.1275716486967107)
            if hbond["bond_checked"]:    # able to check the hydrogen bond
                if hbond["badness"] > max_badness:
                    max_badness = hbond["badness"]
                if hbond["distance"] > max_distance:
                    max_distance = hbond["distance"]
                if not math.isnan(hbond["angle"]) and hbond["angle"] < min_angle:
                    min_angle = hbond["angle"]

    # check new cutoffs in every case, call that "Newest FR3D"
    disqualified_hbond = False
    keep_reasons = []
    best_cutoff_distance = 9999
    best_interaction = ''
    new_python_subcat = -1

    near_discrepancy_cutoff = 1.0     # maximum discrepancy to report as a near pair

    # try all variations of this basepair type, including csS and cSs
    possible_interactions = []
    for poss_int in nt_nt_cutoffs_bc.keys():
        if interaction.lower() in poss_int.lower():
            possible_interactions.append(poss_int)

    for bp in possible_interactions:

        subcats = sorted(nt_nt_cutoffs_bc[bp].keys())

        for subcat in subcats:
            reasons = []  # keep track of reasons for losing a classification
            cutoff_distance = 0

            bpc = bp
            bpsc = subcat

            cutoff = nt_nt_cutoffs_bc[bp][subcat]

            if cutoff["anglemin"] > cutoff["anglemax"]:
                angle_out_of_order = True

            if datapoint['x'] < cutoff["xmin"]:
                cutoff_distance += cutoff["xmin"] - datapoint['x']
                reasons.append("xmin")
            elif datapoint['x'] > cutoff["xmax"]:
                cutoff_distance +=  datapoint['x'] - cutoff["xmax"]
                reasons.append("xmax")

            if datapoint['y'] < cutoff["ymin"]:
                cutoff_distance += cutoff["ymin"] - datapoint['y']
                reasons.append("ymin")
            elif datapoint['y'] > cutoff["ymax"]:
                cutoff_distance +=  datapoint['y'] - cutoff["ymax"]
                reasons.append("ymax")

            if "radiusmax" in cutoff:
                radius = math.sqrt(datapoint['x']**2 + datapoint['y']**2)
                if radius > cutoff["radiusmax"]:
                    cutoff_distance +=  radius - cutoff["radiusmax"]
                    reasons.append("radiusmax")

            if datapoint['z'] < cutoff["zmin"]:
                cutoff_distance += cutoff["zmin"] - datapoint['z']
                reasons.append("zmin")
            elif datapoint['z'] > cutoff["zmax"]:
                cutoff_distance +=  datapoint['z'] - cutoff["zmax"]
                reasons.append("zmax")

            if 'normal_Z' in datapoint:
                if np.sign(datapoint['normal_Z']) * np.sign(cutoff["normalmax"]) < 0:
                    # normals point in different directions entirely
                    cutoff_distance += 100
                elif datapoint['normal_Z'] < cutoff["normalmin"]:
                    cutoff_distance += 3*(cutoff["normalmin"] - datapoint['normal_Z'])
                    reasons.append("nmin")
                elif datapoint['normal_Z'] > cutoff["normalmax"]:
                    cutoff_distance += 3*(datapoint['normal_Z'] - cutoff["normalmax"])
                    reasons.append("nmax")
            else:
                cutoff_distance += 1
                reasons.append("nmax")

            if 'angle_in_plane' in datapoint:
                if cutoff["anglemin"] < cutoff["anglemax"]:
                    # usual order where min < max
                    if datapoint['angle_in_plane'] < cutoff["anglemin"]:
                        cutoff_distance += 0.05*(cutoff["anglemin"] - datapoint['angle_in_plane'])
                        reasons.append("angle")
                    if datapoint['angle_in_plane'] > cutoff["anglemax"]:
                        cutoff_distance += 0.05*(datapoint['angle_in_plane'] - cutoff["anglemax"])
                        reasons.append("angle")
                else:
                    # min might be 260 and max might be -60, looking for angles above 260 or below -60
                    if datapoint['angle_in_plane'] < cutoff["anglemin"] and \
                       datapoint['angle_in_plane'] > cutoff["anglemax"]:
                        cutoff_distance += 0.05*min(cutoff["anglemin"]-datapoint["angle_in_plane"],datapoint['angle_in_plane']-cutoff["anglemax"])
                        reasons.append("angle")
            else:
                cutoff_distance += 1
                reasons.append("angle")

            if datapoint['maxgap'] > cutoff["gapmax"]:
                cutoff_distance += 4*(datapoint['maxgap'] - cutoff["gapmax"])  # strong penalty
                reasons.append('gap')

            if cutoff_distance < best_cutoff_distance:
                # met cutoffs better than with previous subcategory
                keep_reasons = reasons
                new_python_annotation = bp   # interaction or reason for rejection
                if cutoff_distance > near_discrepancy_cutoff and not 'n' in bp:
                    new_python_annotation = 'n' + bp

                new_python_subcat = bpsc
                best_interaction = bp
                best_cutoff_distance = cutoff_distance



    if best_cutoff_distance > 0:
        new_python_annotation = ",".join(keep_reasons)

    if max_distance > 4 or min_angle < 100 or max_badness > 3:
        disqualified_hbond = True
        new_python_annotation += ',hbond'

    datapoint['best_interaction'] = best_interaction
    datapoint['best_subcategory'] = new_python_subcat
    datapoint['best_cutoff_distance'] = best_cutoff_distance
    datapoint['keep_reasons'] = keep_reasons
    datapoint['disqualified_hbond'] = disqualified_hbond

    # keep data to put in scatterplots and HTML pages
    # only keep pairs if their parameters are somewhat close to the cutoffs for some category

    # store data about this pair for output
    pdata = datapoint
    pdata['max_distance'] = max_distance
    pdata['min_angle'] = min_angle
    pdata['max_badness'] = max_badness
    pdata['python_annotation'] = python_annotation          # original annotation from NA_pairwise_interactions
    pdata['new_python_annotation'] = new_python_annotation
    pdata['basepair_subcat'] = new_python_subcat
    pdata['best_cutoff_distance'] = best_cutoff_distance

    if 'sugar_ribose' in datapoint:
        pdata['sugar_ribose'] = datapoint['sugar_ribose']
    else:
        pdata['sugar_ribose'] = ''

    if not 'angle_in_plane' in pdata:
        pdata['angle_in_plane'] = 0.0  # so distance can be calculated


    return datapoint, pdata, angle_out_of_order


#=========================================== Main block =============================

if __name__=="__main__":

    # test
    # table, combination_to_LSW_status = generate_LW_family_table('cWW')
    # table, combination_to_LSW_status = generate_LW_family_table('tHS')
    # print(table)
    # print(combination_to_LSW_status)
    # print(crashmenow)

    symmetric_basepair_list = ['cWW','cWw','cwW','tWW','cHH','tHH','cSS','cSs','csS','tSS','tSs','tsS']
    Leontis_Westhof_basepairs = ['cWW','tWW','cWH','tWH','cWS','tWS','cHH','tHH','cHS','tHS','cSS','tSS','cWB']
    LW12 = ['cWW','tWW','cWH','tWH','cWS','tWS','cHH','tHH','cHS','tHS','cSS','tSS']
    LW12b = ['cWw','tWW','cWH','tWH','cWS','tWS','cHh','tHh','cHS','tHS','cSs','tSs']
    LW12c = ['cWw','tWW','cWH','tWH','cWS','tWS','cHh','tHh','cHS','tHS','cSs','tSS']
    LW18 = ['cWW','tWW','cWH','cHW','tWH','tHW','cWS','cSW','tWS','tSW','cHH','tHH','cHS','cSH','tHS','tSH','cSs','csS','tSs','tsS']

    base_combination_list = ['G,C','A,A','A,C','A,G','A,U','C,C','C,U','G,G','G,U','U,U']
    base_combination_list = ['A,A','A,C','A,G','A,U','C,C','G,C','C,U','G,G','G,U','U,U']

    base_combination_to_interaction = {}
    base_combination_to_interaction['A,A'] = ['cWw','tWW','cWH','tWH','cWS','tWS','cHh','tHH','cHS','tHS','cSs','tSs','cWB']
    base_combination_to_interaction['A,C'] = LW18 + ['cWB'] + ['cBW']
    base_combination_to_interaction['A,G'] = LW18
    base_combination_to_interaction['A,U'] = LW18
    base_combination_to_interaction['C,C'] = LW12b + ['cWB']
    base_combination_to_interaction['C,U'] = LW18
    base_combination_to_interaction['G,C'] = LW18
    base_combination_to_interaction['G,G'] = LW12c + ['cWB']
    base_combination_to_interaction['G,U'] = LW18 + ['cWB']
    base_combination_to_interaction['U,U'] = LW12b

    base_combination_to_interaction['DA,DA'] = LW12b
    base_combination_to_interaction['DA,DC'] = LW18
    base_combination_to_interaction['DA,DG'] = LW18
    base_combination_to_interaction['DA,DT'] = LW18
    base_combination_to_interaction['DC,DC'] = LW12b
    base_combination_to_interaction['DC,DT'] = LW18
    base_combination_to_interaction['DG,DC'] = LW18
    base_combination_to_interaction['DG,DG'] = LW12c
    base_combination_to_interaction['DG,DT'] = LW18
    base_combination_to_interaction['DT,DT'] = LW12b

    red   = [1,0,0] # red
    black = [0,0,0] # black
    cyan  = [0,1,1] # cyan
    blue  = [0,0,1] # blue
    orange = [1,165/255,0] # orange
    purple = [73/255,0,146/255] # purple
    purple = [218/255,112/255,214/255] # orchid, shows up better than purple

    write_html_pages = False
    write_html_pages = True

    make_plots = False
    make_plots = True

    all_agree = '11111'            # indicates how many annotators are being compared
    compare_annotators = True      # tally up how many pairs are annotated by FR3D, RNAview, DSSR, PDB, contacts
    compare_annotators = False

    only_modified = False
    only_modified = True           # only show basepairs involving modified nucleotides

    DNA = True
    DNA = False

    # zzz

    resolution_list = ['2.5A','3.0A','1.5A','2.0A']
    resolution_list = ['2.0A','2.5A','3.0A']
    resolution_list = ['1.5A','2.0A','2.5A','3.0A']
    resolution_list = ['2.5A','3.0A']
    resolution_list = ['3.0A']
    resolution_list = ['1.5A']
    resolution_list = ['2.5A']
    resolution_list = ['2.0A']
    resolution_list = ['1.5A','2.0A','2.5A','3.0A']

    if compare_annotators:
        make_plots = False

    for resolution in resolution_list:
        if DNA:
            from DNA_2A_list import PDB_list   # define PDB_list as a list of DNA structures
            PDB_list.remove('5G35')        # has AG pairs that overlap
            PDB_list = sorted(PDB_list)

            base_combination_list = ['DA,DA','DA,DC','DA,DG','DA,DT','DC,DC','DG,DC','DC,DT','DG,DG','DG,DT','DT,DT']
            data_file = []
        else:
            if resolution == '1.5A':
                PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.285/1.5A/csv']
            elif resolution == '2.0A':
                PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.285/2.0A/csv']
            elif resolution == '2.5A':
                PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.285/2.5A/csv','8GLP','8B0X']
            elif resolution == '3.0A':
                PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.285/3.0A/csv','8GLP','8B0X']
            data_file = readPDBDatafile()  # available PDB structures, resolutions, chains

        PDB_IFE_Dict = map_PDB_list_to_PDB_IFE_dict(PDB_list)

        # load all datapoints on pairs of bases, whether annotated as paired or not
        all_PDB_ids = sorted(PDB_IFE_Dict.keys())

        print("Resolution %s, working on %d PDB files" % (resolution,len(all_PDB_ids)))

        representative_chains = set(['8GLP|1|L5','8GLP|1|L8','8GLP|1|S2','8B0X|1|a','8B0X|1|A'])
        for PDB, chains in PDB_IFE_Dict.items():
            for chain in chains.split("+"):
                representative_chains.add(chain)

        PDB_skip_set = set(['1R9F','5NXT','4KTG','7JIL'])

        if '8B0X' in PDB_list:
            PDB_skip_set.add('5J7L')
            PDB_skip_set.add('7K00')
            PDB_skip_set.add('4YBB')

        # load annotations of these PDB files from the BGSU RNA server
        pair_to_interaction_matlab   = defaultdict(str)
        pair_to_interaction_dssr     = defaultdict(str)
        pair_to_interaction_rnaview  = defaultdict(str)
        pair_to_interaction_pdb      = defaultdict(str)
        pair_to_interaction_contacts = defaultdict(str)

        pdb_id_to_annotators = defaultdict(set)   # which programs annotate pairs in each structure?

        pair_to_datapoint = defaultdict(dict)
        unit_id_to_glycosidic = defaultdict(str)

        # load output files from NA_pairwise_interactions
        print("Loading FR3D     annotations from %s" % outputNAPairwiseInteractions)
        for PDB_id in all_PDB_ids:
            pair_to_datapoint_file = outputNAPairwiseInteractions + "%s_pairs_v1.pickle" % PDB_id
            try:
                if sys.version_info[0] < 3:
                    new_dict = pickle.load(open(pair_to_datapoint_file,'rb'))
                else:
                    new_dict = pickle.load(open(pair_to_datapoint_file,'rb'),encoding = 'latin1')
                if len(new_dict.keys()) > 0:
                    # make sure at least one pair was annotated
                    found_basepair = False
                    for pair,datapoint in new_dict.items():
                        if 'basepair' in datapoint and datapoint['basepair']:
                            found_basepair = True

                    if found_basepair:
                        pdb_id_to_annotators[PDB_id].add('FR3D')
                        for pair,datapoint in new_dict.items():
                            # reduce memory usage by only storing pairs we will process
                            u1,u2 = pair

                            # only keep pairs from representative chains
                            fields1 = u1.split("|")
                            chain1 = "|".join(fields1[0:3])
                            if not chain1 in representative_chains and not DNA:
                                continue

                            fields2 = u2.split("|")
                            chain2 = "|".join(fields2[0:3])
                            if not chain2 in representative_chains and not DNA:
                                continue

                            # pair_to_datapoint.update(new_dict)
                            pair_to_datapoint[pair] = datapoint

            except:
                print("Not able to load python_fr3d annotations for %s from %s" % (PDB_id,pair_to_datapoint_file))


            unit_annotation_file = os.path.join(outputNAPairwiseInteractions,"%s_glycosidic.txt" % PDB_id)
            if os.path.exists(unit_annotation_file):
                #print('Reading glycosidic bond conformations for %s' % (PDB_id))
                with open(unit_annotation_file,'r') as f:
                    lines = f.readlines()
                for line in lines:
                    fields = line.split()
                    unit_id = fields[0]
                    glycosidic = fields[1]
                    unit_id_to_glycosidic[unit_id] = glycosidic



        """
        print(pair_to_datapoint[('8GLP|1|L5|U|3914','8GLP|1|L5|U|3915')])
        print(pair_to_datapoint[('8GLP|1|L5|U|3915','8GLP|1|L5|U|3914')])
        """

        # keep track of pairs, prioritizing those with annotations like cWw over those like cWH/cHW
        pair_to_priority = add_pairs_in_order({},pair_to_datapoint,1)

        if not DNA and not compare_annotators:
            print('Loading Matlab   annotations')
            for PDB_id in all_PDB_ids:
                interaction_to_pairs = load_Matlab_FR3D_pairs(PDB_id)
                num_pairs = 0
                for interaction in interaction_to_pairs.keys():
                    # store all basepairs and only basepairs
                    if "c" in interaction or "t" in interaction:
                        num_pairs += 1
                        for pair in interaction_to_pairs[interaction]:
                            pair_to_interaction_matlab[pair] = interaction
                if num_pairs == 0:
                    #print("No Matlab-annotated pairs in %s" % PDB_id)
                    pass
                    #PDB_skip_set.add(PDB_id)

            #print("Skipping %d PDB files because they have no Matlab annotation to compare to" % len(PDB_skip_set))
            #print("Found Matlab annotations in %s files" % (len(all_PDB_ids)-len(PDB_skip_set))

        print('Loading RNAview  annotations from %s' % rnaview_basepair_path)
        for PDB_id in all_PDB_ids:
            interaction_to_triples = load_rnaview_basepairs(PDB_id)
            for interaction in interaction_to_triples.keys():
                # store all basepairs and only basepairs
                if "c" in interaction or "t" in interaction:
                    for u1,u2,crossing in interaction_to_triples[interaction]:
                        pair_to_interaction_rnaview[(u1,u2)] = interaction
                    if len(interaction_to_triples[interaction]) > 0:
                        pdb_id_to_annotators[PDB_id].add('RNAview')

        pair_to_priority = add_pairs_in_order(pair_to_priority,pair_to_interaction_rnaview,2)

        print('Loading DSSR     annotations from %s' % dssr_basepair_path)
        for PDB_id in all_PDB_ids:
            interaction_to_triples = load_dssr_basepairs(PDB_id)
            for interaction in interaction_to_triples.keys():
                # store all basepairs and only basepairs
                if "c" in interaction or "t" in interaction:
                    for u1,u2,crossing in interaction_to_triples[interaction]:
                        pair_to_interaction_dssr[(u1,u2)] = interaction
                    if len(interaction_to_triples[interaction]) > 0:
                        pdb_id_to_annotators[PDB_id].add('DSSR')

        pair_to_priority = add_pairs_in_order(pair_to_priority,pair_to_interaction_dssr,3)

        print('Loading PDB      annotations from %s' % pdb_basepair_path)
        for PDB_id in all_PDB_ids:
            interaction_to_triples = load_pdb_basepairs(PDB_id)
            for interaction in interaction_to_triples.keys():
                # store all basepairs and only basepairs
                if "c" in interaction or "t" in interaction:
                    for u1,u2,crossing in interaction_to_triples[interaction]:
                        pair_to_interaction_pdb[(u1,u2)] = interaction
                    if len(interaction_to_triples[interaction]) > 0:
                        pdb_id_to_annotators[PDB_id].add('PDB')

        pair_to_priority = add_pairs_in_order(pair_to_priority,pair_to_interaction_pdb,5)

        if DNA:
            contacts_filename = os.path.join('C:/Users/zirbel/Documents/RNA/contacts program','pairing_DNA.txt')
        else:
            contacts_filename = os.path.join('C:/Users/zirbel/Documents/RNA/contacts program','pairing_1599_RNA.txt')
        print('Loading contacts annotations from %s' % contacts_filename)
        pair_to_interaction_contacts = load_contacts_basepairs(contacts_filename)
        for u1,u2 in pair_to_interaction_contacts.keys():
            fields = u1.split("|")
            pdb_id = fields[0].upper()
            pdb_id_to_annotators[pdb_id].add('contacts')

        pair_to_priority = add_pairs_in_order(pair_to_priority,pair_to_interaction_pdb,4)

        if compare_annotators:
            all_annotate_counter = 0
            for pdb_id, annotators in sorted(pdb_id_to_annotators.items()):
                if len(annotators) < len(all_agree):
                    if len(annotators) > 1 or not list(annotators)[0] == 'contacts':
                        print('%s only annotated by %s' % (pdb_id,sorted(annotators)))
                else:
                    all_annotate_counter += 1
            print('%d files are annotated by all %d annotators' % (all_annotate_counter,len(all_agree)))


        # remove pairs with symmetry operators, alternate ids, insertion codes when comparing annotators
        complicated_unit_id = set()
        if compare_annotators:
            for u1,u2 in pair_to_priority.keys():
                fields1 = u1.split("|")
                fields2 = u2.split("|")

                if len(fields1) > 5:
                    complicated_unit_id.add(u1)                      # blacklist the full unit id
                    complicated_unit_id.add("|".join(fields1[0:5]))  # blacklist the plain unit id as well
                    complicated_unit_id.add("|".join(fields2[0:5]))  # blacklist the paired plain unit id as well
                    #print('Omitting',u1,"|".join(fields1[0:5]),u2,"|".join(fields2[0:5]))
                if len(fields2) > 5:
                    complicated_unit_id.add(u2)                      # blacklist the full unit id
                    complicated_unit_id.add("|".join(fields2[0:5]))  # blacklist the plain unit id as well
                    complicated_unit_id.add("|".join(fields1[0:5]))  # blacklist the paired plain unit id as well
                    #print('Omitting',u1,"|".join(fields1[0:5]),u2,"|".join(fields2[0:5]))

                # some x-ray structures like 3CMY have a model 2 that RNAview does not mark as model 2,
                # so it gets reported as being in model 1.  Remove every such pair, even from model 1.
                if not fields1[1] == '1' or not fields2[1] == '1':
                    complicated_unit_id.add(u1)                      # blacklist the full unit id
                    complicated_unit_id.add(u2)                      # blacklist the full unit id
                    fields1[1] = '1'
                    complicated_unit_id.add("|".join(fields1[0:5]))  # blacklist any model 1 version of the unit id
                    fields2[1] = '1'
                    complicated_unit_id.add("|".join(fields2[0:5]))  # blacklist any model 1 version of the unit id
                    print('Removing',u1,"|".join(fields1[0:5]),u2,"|".join(fields2[0:5]),' because of model numbers')


        print('Found %d pair_to_priority pairs to work with' % len(pair_to_priority))

        standard_nucleotides = ['A','C','G','U','DA','DC','DG','DT']

        # cull out pairs that we don't want to consider
        pair_to_priority_final = {}
        for pair,priority in pair_to_priority.items():
            u1,u2 = pair
            if u1 in complicated_unit_id:
                continue

            if u2 in complicated_unit_id:
                continue

            fields1 = u1.split("|")

            # only keep pairs from representative chains
            chain1 = "|".join(fields1[0:3])
            if not chain1 in representative_chains and not DNA:
                continue

            # DSSR has hyphens because of symmetry operators being applied
            if '-' in chain1:
                continue

            fields2 = u2.split("|")
            chain2 = "|".join(fields2[0:3])
            if not chain2 in representative_chains and not DNA:
                continue

            # DSSR has hyphens because of symmetry operators being applied
            if '-' in chain2:
                continue

            if only_modified and fields1[3] in standard_nucleotides and fields2[3] in standard_nucleotides:
                continue

            # if both are from a symmetry operator, skip, because we probably already have it
            if len(fields1) == 9 and len(fields2) == 9:
                #print('Skipping %s - %s because both have a symmetry operator' % (u1,u2))
                continue

            pdb_id = fields1[0]

            # only include pairs where all annotators are able to make annotations
            if compare_annotators and not len(pdb_id_to_annotators[pdb_id]) == len(all_agree):
                continue

            # some PDB ids are a problem, just skip them
            if pdb_id in PDB_skip_set:
                continue

            pair_to_priority_final[(u1,u2)] = priority

        # set the order in which to process pairs; doing FR3D first gets more consistent unit orders when symmetric
        pairs_in_order = sorted(pair_to_priority_final.keys(), key=lambda x:pair_to_priority_final[x])

        # try to reduce memory usage
        pair_to_priority_final = {}

        print('Found %d pairs_in_order pairs to work with' % len(pairs_in_order))

        # store pairs by their base combination for faster retrieval
        # store modified nucleotides according to their parent nucleotide
        bc_to_pairs_in_order = {}
        for pair in pairs_in_order:
            b1 = pair[0].split("|")[3]
            if b1 in modified_base_to_parent and not compare_annotators and not b1 in ['DA','DC','DG','DT']:
                b1 = modified_base_to_parent[b1]

            b2 = pair[1].split("|")[3]
            if b2 in modified_base_to_parent and not compare_annotators and not b2 in ['DA','DC','DG','DT']:
                b2 = modified_base_to_parent[b2]

            if not b1 in ['A','C','G','U','DA','DC','DG','DT'] and not compare_annotators:
                print('Unknown nucleotide %s' % pair[0])
            if not b2 in ['A','C','G','U','DA','DC','DG','DT'] and not compare_annotators:
                print('Unknown nucleotide %s' % pair[1])

            bc = b1+","+b2

            if not bc in bc_to_pairs_in_order:
                bc_to_pairs_in_order[bc] = []

            bc_to_pairs_in_order[bc].append(pair)

        # try to reduce memory usage
        pairs_in_order = {}

        count_annotations = {}
        count_annotations_by_group = {}

        big_min_distance_list = []

        # loop over specified base combinations
        for base_combination in base_combination_list:

            interactions_processed = set([])

            # loop over conceivable and distinct interactions for that base combination
            for interaction in base_combination_to_interaction[base_combination]:

                pair_counted = set()  # only display / count each pair once on each HTML page
                pair_shown = set()    # only display each pair once on each HTML page

                python_true_count = 0
                python_near_count = 0

                # faster re-check when working on a specific category; some use lowercase letters
                # if not interaction in ['cWB','cBW']:
                #      continue

                angle_out_of_order = False   # set to True if there is a cutoff like 250 to -60 degrees

                interaction_upper = interaction[0] + interaction[1:3].upper()
                interaction_lower = interaction.lower()

                # standard order for base combination and filename
                # use this for counting and for filenames
                bc_filename = base_combination
                inter_filename = interaction_upper

                # base combination and basepair to use in filenames
                # reverse order of base combination with interactions like tsS because
                # file names for tsS and tSs will be the same on windows
                # also use standard order for the six combinations listed below
                if interaction[1].islower() or interaction_upper in ['cHW','tHW','cSW','tSW','cSH','tSH','cBW']:
                    bases = bc_filename.split(",")
                    bc_filename = bases[1] + "," + bases[0]
                    inter_filename = interaction_upper[0] + interaction_upper[2] + interaction_upper[1]

                print('Working on %s %s' % (base_combination,interaction))
                print('Base combination %s filename %s' % (bc_filename,inter_filename))

                nt1_seq, nt2_seq = base_combination.split(",")

                # accumulate data specific to this interaction and base combination for plottings
                xvalues = []    # x coordinate of the second base
                yvalues = []    # y coordinate of the second base
                rvalues = []    # radius of second base
                tvalues = []    # angle of second base relative to first
                zvalues = []    # z coordinate of the second base
                avalues = []    # rotation angle of second base
                gvalues = []    # maximum gap value between meeting edges
                nvalues = []    # third component of normal vector of second base
                colors2d  = []  # store the color to use
                sizes = []      # store the size of dot to use

                # accumulate information about hydrogen bonds; multiple ones per basepair
                hdvalues = []   # hydrogen bond distances
                havalues = []   # hydrogen bond angles
                hbvalues = []   # hydrogen bond badness measures
                hcolors2d = []  # store color to use for hydrogen bond dots
                hsizes    = []  # store size of dot to use for hydrogen bonds

                # store data to produce an HTML table
                pair_data = []  # list of data dictionaries to print in a table

                # loop over pairs for which we have data, finding those with the interaction and base combination
                for pair in bc_to_pairs_in_order.get(base_combination,[]):

                    if not pair in pair_to_datapoint:
                        print('Pair %s-%s does not have datapoint information' % (pair[0],pair[1]))
                        pair_to_datapoint[pair] = {}
                        pair_to_datapoint[pair]['nt1_seq'] = pair[0].split("|")[3]
                        pair_to_datapoint[pair]['nt2_seq'] = pair[1].split("|")[3]
                        pair_to_datapoint[pair]['basepair'] = 'N/A'
                        pair_to_datapoint[pair]['basepair_subcategory'] = 9

                    datapoint = pair_to_datapoint[pair]

                    if compare_annotators and pair in pair_counted:
                        #print('Already counting and displaying %s,%s' % pair)
                        continue

                    if not compare_annotators and pair in pair_shown:
                        continue

                    fields1 = pair[0].split("|")
                    pdb_id = fields1[0]

                    # check different annotation schemes to decide how to show this datapoint
                    # check python_fr3d annotation that generated pair_to_datapoint
                    if 'basepair' in datapoint and interaction_lower in datapoint['basepair'].lower():
                        python_fr3d = True
                        if "n" in datapoint['basepair']:
                            python_near = True
                            python_true = False
                        else:
                            python_near = False
                            python_true = True
                    else:
                        python_fr3d = False
                        python_near = False
                        python_true = False

                    if 'basepair' in datapoint:
                        python_annotation = datapoint['basepair']
                    else:
                        python_annotation = ""

                    matlab_annotation = pair_to_interaction_matlab[pair]
                    if len(matlab_annotation) > 0 and interaction_lower in matlab_annotation.lower():
                        Matlab = True
                    else:
                        Matlab = False

                    dssr_annotation = pair_to_interaction_dssr[pair]
                    if len(dssr_annotation) > 0 and interaction_lower in dssr_annotation.lower():
                        dssr = True
                    else:
                        dssr = False

                    contacts_annotation = pair_to_interaction_contacts[pair]
                    if len(contacts_annotation) > 0 and interaction_lower in contacts_annotation.lower():
                        contacts = True
                    else:
                        contacts = False

                    rnaview_annotation = pair_to_interaction_rnaview[pair].replace("n","!")
                    if len(rnaview_annotation) > 0 and interaction_lower in rnaview_annotation.lower():
                        rnaview = True
                        if "!" in rnaview_annotation:
                            rnaview_near = True
                            rnaview_true = False
                        else:
                            rnaview_near = False
                            rnaview_true = True
                    else:
                        rnaview = False
                        rnaview_near = False
                        rnaview_true = False

                    # skip cases where RNAview may need to annotate across chains in different bundles
                    if compare_annotators and not rnaview and (python_true or python_near or dssr or contacts) and not chain1 == chain2 and not DNA:
                        print('Skipping %s-%s %s because the chains may be from different bundles' % (pair[0],pair[1],interaction))
                        continue

                    # since PDB annotations only tell the family and since we are only going to list
                    # each pair once, only record a PDB annotation when you can tell what edges are used
                    pdb_annotation = pair_to_interaction_pdb[pair]
                    pdb = False
                    if len(pdb_annotation) > 0 and interaction_lower in pdb_annotation.lower():
                        if (rnaview or python_fr3d or dssr or contacts) or inter_filename in ['cWW','tWW','cHH','tHH']:
                            pdb = True

                    # correctly identify instances as cSs/csS or tSs/tsS; delay csS and tsS to get cSs and tSs
                    if interaction in ['cSs','tSs'] and not nt1_seq == nt2_seq:
                        # you can only delay on cSs or tSs, not on csS or tsS
                        if python_fr3d:
                            if not interaction in datapoint['basepair']:
                                # python_fr3d is csS or tsS, does not match cSs or tSs, so skip until later
                                continue
                        elif rnaview_true or dssr or contacts or pdb:
                            print('Checking %s-%s for %s' % (pair[0],pair[1],interaction))
                            datapoint, pdata, angle_order = evaluate_pair_from_datapoint(datapoint,interaction,nt_nt_cutoffs[base_combination])
                            #print("evaluated datapoint:")
                            #print_dictionary(datapoint)
                            if not interaction in datapoint['best_interaction']:
                                print('Delaying %s-%s with %s in category %s' % (pair[0],pair[1],datapoint['best_interaction'],interaction))
                                continue

                    if interaction in ['csS','tsS'] and not nt1_seq == nt2_seq and not compare_annotators:
                        if python_fr3d:
                            if not interaction in datapoint['basepair']:
                                # python_fr3d is cSs, don't show it here
                                print('Skipping %s %s with %s for %s' % (pair[0],pair[1],datapoint['basepair'],interaction))
                                continue

                    # # why are we missing some pairs that are delayed?
                    # if interaction in ['csS','tsS'] and not nt1_seq == nt2_seq:
                    #     if not python_fr3d and (rnaview_true or dssr or pdb):
                    #         print('Found1 %s-%s on %s' % (pair[0],pair[1],interaction))


                    # if all annotation systems have annotated at least one basepair in this PDB file,
                    # tally how well they agree
                    annotator = ""
                    annotator_count = 0
                    if compare_annotators and len(pdb_id_to_annotators[pdb_id]) == len(all_agree):
                        if python_true or rnaview_true or dssr or contacts or pdb:
                            booleans = (python_true,rnaview_true,dssr,pdb,contacts)
                            for i,truth in enumerate(booleans):
                                if truth:
                                    annotator += "1"
                                    annotator_count += 1
                                else:
                                    annotator += "0"

                            booleans = (python_true,python_near,rnaview_true,rnaview_near,dssr,pdb,contacts)

                            if not inter_filename in count_annotations:
                                count_annotations[inter_filename] = {}
                                count_annotations_by_group[inter_filename] = {}
                            if not bc_filename in count_annotations[inter_filename]:
                                # one for each program, in order
                                count_annotations[inter_filename][bc_filename] = [0,0,0,0,0,0,0]
                                count_annotations_by_group[inter_filename][bc_filename] = defaultdict(int)

                            for i,truth in enumerate(booleans):
                                if truth:
                                    count_annotations[inter_filename][bc_filename][i] += 1

                            count_annotations_by_group[inter_filename][bc_filename][annotator] += 1

                            pair_counted.add(pair)
                            pair_counted.add(reverse(pair))

                    if reverse(pair) in pair_to_datapoint:
                        r_datapoint = pair_to_datapoint[reverse(pair)]
                    else:
                        r_datapoint = {}

                    # do we have all of the parameters that are checked for cutoffs?
                    have_full_data = check_full_data(datapoint)

                    if python_true:
                        python_true_count += 1
                    if python_near:
                        python_near_count += 1

                    if (annotator_count > 0) or (python_true or Matlab or dssr or rnaview_true or pdb or contacts) or (not compare_annotators and python_near):

                        # evaluate the quality of the match to the current pair, for scatterplots and all
                        #print('Evaluating pair %s - %s for interaction %s' % (pair[0],pair[1],interaction))
                        datapoint, pdata, angle_order = evaluate_pair_from_datapoint(datapoint,interaction,nt_nt_cutoffs[base_combination])

                        if angle_order:
                            angle_out_of_order = True

                        # use the reversed pair here?
                        if nt1_seq == nt2_seq and interaction in symmetric_basepair_list and datapoint['best_cutoff_distance'] > 0:
                            #print('Evaluating reversed pair %s - %s' % reverse(pair))
                            r_datapoint, r_pdata, angle_order = evaluate_pair_from_datapoint(r_datapoint,interaction,nt_nt_cutoffs[base_combination])
                            # print_dictionary(r_datapoint)
                            # print()
                            # if new match is better and still matches interaction, reverse the pair, use r_datapoint, etc.
                            if r_datapoint["best_cutoff_distance"] < datapoint["best_cutoff_distance"]:
                                if r_datapoint["python_annotation"].lower() == interaction_lower:
                                    print("Pair %s - %s from %s is better reversed and makes %s" % (pair[0],pair[1],interaction,python_annotation))
                                    pair = reverse(pair)
                                    datapoint = r_datapoint
                                    pdata = r_pdata
                                    pdata['python_annotation'] = reverse_edges(python_annotation)

                            elif nt1_seq == 'G' and interaction == 'cWW' and rnaview_annotation == 'cHW':
                                print("Using reversed pair for %s - %s with cWW and rnaview cHW but I don't remember why" % pair)
                                pair = reverse(pair)
                                datapoint = r_datapoint
                                pdata = r_pdata
                                pdata['python_annotation'] = reverse_edges(python_annotation)

                        best_cutoff_distance = datapoint["best_cutoff_distance"]

                        # save the pair if it is good enough to list in the table
                        if compare_annotators \
                            or (not interaction in nt_nt_cutoffs[base_combination] and (rnaview or pdb or dssr or contacts)) \
                            or best_cutoff_distance < 2 \
                            or ('sugar_ribose' in datapoint and datapoint['sugar_ribose'] == 'cSR' and best_cutoff_distance < 2) \
                            or (rnaview and best_cutoff_distance < 2) \
                            or (pdb and best_cutoff_distance < 2) \
                            or (dssr and best_cutoff_distance < 2) \
                            or (contacts and best_cutoff_distance < 2) \
                            or (matlab_annotation and not "n" in matlab_annotation and best_cutoff_distance < 2):

                            pdata['unit_id_1'] = pair[0]
                            pdata['unit_id_2'] = pair[1]
                            pdata['glycosidic1'] = unit_id_to_glycosidic[pair[0]]
                            pdata['glycosidic2'] = unit_id_to_glycosidic[pair[1]]
                            pdata['matlab_annotation'] = matlab_annotation
                            pdata['dssr_annotation'] = dssr_annotation
                            pdata['rnaview_annotation'] = rnaview_annotation
                            pdata['pdb_annotation'] = pdb_annotation
                            pdata['contacts_annotation'] = contacts_annotation
                            pdata['annotator_count'] = annotator_count  # how many annotate as such
                            pdata['annotator'] = annotator

                            # store information about hydrogen bonds for this base combination and interaction
                            if "LW_to_atom_sets" in datapoint:
                                if interaction in datapoint["LW_to_atom_sets"]:
                                    pdata["atom_sets"] = sorted(datapoint["LW_to_atom_sets"][interaction])
                                    #pdata["atom_set_to_results"] = datapoint["atom_set_to_results"]  # it's already there in most cases
                                else:
                                    # avoid carryover from previous times this datapoint was encountered
                                    pdata["atom_sets"] = []
                            else:
                                pdata["atom_sets"] = []

                            if annotator_count == len(annotator):
                                pdata['annotator_all'] = 1
                            else:
                                pdata['annotator_all'] = 0

                            pair_data.append(pdata)

                            pair_shown.add(pair)
                            pair_shown.add(reverse(pair))

                            if pdata['python_annotation'] and 'min_distance' in datapoint and datapoint['min_distance'] > 3.5:
                                if not (bc_filename[2] in ['C','U'] and inter_filename == 'cSS'):
                                    big_min_distance_list.append((datapoint['min_distance'],bc_filename,inter_filename,pdata['python_annotation']))

                            keep_reasons = datapoint['keep_reasons']

                            # collect data for scatterplots
                            if make_plots and 'angle_in_plane' in datapoint and abs(datapoint['angle_in_plane']) > 0.0001:
                                xvalues.append(datapoint['x'])
                                yvalues.append(datapoint['y'])
                                rvalues.append(math.sqrt(datapoint['x']**2 + datapoint['y']**2))
                                tvalues.append(math.atan2(datapoint['y'],datapoint['x'])*180/3.141592654)
                                zvalues.append(datapoint['z'])
                                if 'gap21' in datapoint:
                                    gvalues.append(max(datapoint['gap12'],datapoint['gap21']))
                                else:
                                    gvalues.append(datapoint['gap12'])
                                    datapoint['gap21'] = 0

                                if angle_out_of_order and datapoint['angle_in_plane'] < 0:
                                    avalues.append(datapoint['angle_in_plane']+360)
                                else:
                                    avalues.append(datapoint['angle_in_plane'])

                                nvalues.append(datapoint['normal_Z'])

                                # colors and sizes for scatterplots
                                if 'gap' in datapoint['keep_reasons']:
                                    color = [0,1,0]  # green
                                    hcolor = [0,1,0]
                                    size = 20
                                    hsize = 20
                                elif len(datapoint['keep_reasons']) > 0:
                                    color = cyan
                                    hcolor = cyan
                                    size = 20
                                    hsize = 20
                                # elif 'disqualified_hbond' in datapoint and datapoint['disqualified_hbond']:
                                #     color = orange
                                #     hcolor = orange
                                #     size = 20
                                #     hsize = 20
                                # elif python_fr3d and not Matlab:
                                #     color = purple
                                #     hcolor = purple
                                #     size = 20
                                #     hsize = 20
                                #     #print_datapoint(datapoint)
                                # elif Matlab and not python_fr3d:
                                #     color = red
                                #     hcolor = red
                                #     if "n" in matlab_annotation:
                                #         size = 5
                                #         hsize = 5
                                #     else:
                                #         size = 20
                                #         hsize = 20
                                #     #print_datapoint(datapoint)
                                # elif python_fr3d and Matlab:
                                #     color = black
                                #     size = 1
                                #     hcolor = color
                                #     hsize = size
                                elif python_true:
                                    color = black
                                    size = 1
                                    hcolor = color
                                    hsize = size
                                else:
                                    # python near or other
                                    color = [0.5,0.5,0.5]  # gray
                                    size = 20       # large
                                    hcolor = color
                                    hsize = size
                                colors2d.append(color)
                                sizes.append(size)

                                # record data about hydrogen bonds
                                hcolors2d.append(hcolor)
                                hsizes.append(hsize)
                                # hdvalues.append(max_distance)  # hbond distance
                                # havalues.append(min_angle)  # hbond angle
                                # hbvalues.append(min(6,max_badness))  # hbond badness

                print('%s %s Python true count %4d Python near count %4d Percentage near %6.2f' % (bc_filename,inter_filename,python_true_count,python_near_count,100.0*python_near_count/(1+python_true_count)))

                # image file name, for writing image and for putting image into HTML file
                if len(all_PDB_ids) <= 10:
                    figure_save_file = os.path.join(OUTPUTPATH,"plots","basepairs_%s_%s_%s.png" % (inter_filename,bc_filename,"_".join(all_PDB_ids)))
                    figure_img_src = "plots/basepairs_%s_%s_%s.png" % (inter_filename,bc_filename,"_".join(all_PDB_ids))
                elif 'nrlist' in PDB_list[0]:
                    figure_save_file = os.path.join(OUTPUTPATH,"plots","basepairs_%s_%s_%s.png" % (inter_filename,bc_filename,resolution))
                    figure_img_src = "plots/basepairs_%s_%s_%s.png" % (inter_filename,bc_filename,resolution)
                elif DNA:
                    figure_save_file = os.path.join(OUTPUTPATH,"plots","DNA_basepairs_%s_%s_%s.png" % (inter_filename,bc_filename,resolution))
                    figure_img_src = "plots/DNA_basepairs_%s_%s_%s.png" % (inter_filename,bc_filename,resolution)
                else:
                    figure_save_file = os.path.join(OUTPUTPATH,"plots","basepairs_%s_%s_%d.png" % (inter_filename,bc_filename,len(all_PDB_ids)))
                    figure_img_src = "plots/basepairs_%s_%s_%d.png" % (inter_filename,bc_filename,len(all_PDB_ids))

                # make scatterplots for pairwise combinations of data
                if len(xvalues) > 0 and make_plots:

                    fig = plt.figure(figsize=(11.0, 7.0))

                    ax = fig.add_subplot(2, 3, 1)
                    ax.axis("equal")
                    plot_basepair_cutoffs(base_combination,[interaction],ax,1)
                    ax.scatter(xvalues,yvalues,color=colors2d,marker=".",s=sizes)
                    ax.set_title('x and y for %d %s %s' % (len(xvalues),base_combination,interaction))
                    draw_base(nt1_seq,'default',2,ax)

                    ax = fig.add_subplot(2, 3, 2)
                    plot_basepair_cutoffs(base_combination,[interaction],ax,2)
                    ax.scatter(tvalues,rvalues,color=colors2d,marker=".",s=sizes)
                    ax.set_title('radius vs theta')

                    ax = fig.add_subplot(2, 3, 3)
                    plot_basepair_cutoffs(base_combination,[interaction],ax,3,angle_out_of_order)
                    ax.scatter(avalues,nvalues,color=colors2d,marker=".",s=sizes)
                    ax.set_title('normal vs angle')

                    ax = fig.add_subplot(2, 3, 4)
                    plot_basepair_cutoffs(base_combination,[interaction],ax,4)
                    ax.scatter(gvalues,zvalues,color=colors2d,marker=".",s=sizes)
                    ax.set_title('z versus maxgap')

                    ax = fig.add_subplot(2, 3, 5)
                    plot_basepair_cutoffs(base_combination,[interaction],ax,5,angle_out_of_order)
                    ax.scatter(avalues,xvalues,color=hcolors2d,marker=".",s=hsizes)
                    ax.set_title('x vs angle')

                    ax = fig.add_subplot(2, 3, 6)
                    plot_basepair_cutoffs(base_combination,[interaction],ax,6,angle_out_of_order)
                    ax.scatter(avalues,yvalues,color=hcolors2d,marker=".",s=hsizes)
                    ax.set_title('y vs angle')

                    """
                    ax = fig.add_subplot(2, 3, 5)
                    ax.scatter(hdvalues,havalues,color=hcolors2d,marker=".",s=hsizes)
                    ax.set_title('h-angle vs h-dist')

                    ax = fig.add_subplot(2, 3, 6)
                    #print('plot_lengths',len(gvalues),len(hbvalues))
                    ax.scatter(gvalues,hbvalues,color=hcolors2d,marker=".",s=hsizes)
                    ax.scatter(0,0,color='white')          # push out the x and y axis
                    ax.scatter(0,3.7,color='white')
                    ax.set_title('h-badness vs maxgap')
                    ax.text(0.1,1.5, 'purple=python_fr3d only')
                    ax.text(0.1,2, 'red=Matlab only')
                    ax.text(0.1,2.5, 'green=bad gap')
                    ax.text(0.1,3, 'orange=bad hbond')
                    ax.text(0.1,3.5, 'cyan=bad cutoffs')
                    """

                    # show all plots for this interaction
                    figManager = plt.get_current_fig_manager()
                    figManager.full_screen_toggle()

                    plt.savefig(figure_save_file)
                    #plt.show()
                    plt.close()

                    print("Plotted %5d points for %s %s" % (len(xvalues),base_combination,interaction))

                # write HTML pages listing instances
                if len(pair_data) >= 0 and write_html_pages:
                    # mimic how WebFR3D writes result pages
                    Q = {}
                    if len(all_PDB_ids) <= 10:
                        Q['name'] = "%s %s %s" % (inter_filename,bc_filename,"_".join(all_PDB_ids))
                    elif 'nrlist' in PDB_list[0]:
                        Q['name'] = "%s %s %s" % (inter_filename,bc_filename,resolution)
                    elif DNA:
                        Q['name'] = "%s %s %s" % (inter_filename,bc_filename,resolution)
                    else:
                        Q['name'] = "%s %s %d" % (inter_filename,bc_filename,len(all_PDB_ids))

                    if DNA:
                        Q['name'] = 'DNA_' + Q['name']

                    if only_modified:
                        Q['name'] = 'modified_' + Q['name']

                    Q['numFilesSearched'] = len(all_PDB_ids)
                    Q['searchFiles'] = all_PDB_ids
                    Q['elapsedCPUTime'] = 0
                    Q['userMessage'] = []
                    Q['figure_img_src'] = figure_img_src
                    Q['DNA'] = DNA
                    Q['only_modified'] = only_modified
                    Q['numFound'] = len(pair_data)
                    Q["PDB_data_file"] = data_file

                    # identify the interaction so a table of HTML links can be made
                    if interaction_upper in Leontis_Westhof_basepairs:
                        Q['LW'] = interaction_upper
                    else:
                        # reverse the edges to get the standard family name
                        Q['LW'] = interaction_upper[0] + interaction_upper[2] + interaction_upper[1]

                    Q['resolution'] = resolution

                    # if too many pairs, show some good and also the worst ones
                    # This ordering is not working!  It often keeps the true ones with bad h-bonds,
                    # puts the true ones beyond 300, so the heat map shows all bad instances.
                    if compare_annotators and len(pair_data) > 300:
                        # sort by full agreement or not, random within each category
                        pair_data = sorted(pair_data, key=lambda p: p['annotator_all']+0.0001*random.uniform(0,1))

                        # print('Sorted instances from less to full agreement')
                        # for pd in pair_data:
                        #     print(pd['annotator_all'],pd['annotator_count'])

                        # keep up to 250 that lack total agreement and up to 50 that agree fully
                        order_pair_data = pair_data[0:250] + pair_data[-50:]

                        # discard all the rest
                        other_pair_data = []

                    elif len(pair_data) > 300:
                        # sort by hydrogen bond badness, putting priority on ones that are not Matlab near
                        pair_data = sorted(pair_data, key=lambda p: p['max_badness']+0.0001*random.uniform(0,1)+10*(len(p['matlab_annotation'])==0)+10*(not "n" in p['matlab_annotation']))

                        # keep 50 of the best pairs, and the 250 worst
                        order_pair_data = pair_data[0:50] + pair_data[-250:]

                        # sort the other pairs by python_fr3d annotation and then by badness
                        other_pair_data = sorted(pair_data[51:-251], key=lambda p: (p['python_annotation'],p['max_badness']+0.0001*random.uniform(0,1)))

                        if len(other_pair_data) > 1000:
                            # keep the worst 1000 of them, otherwise you get 30,000 GC cWW or something
                            other_pair_data = other_pair_data[-1000:]
                    else:
                        order_pair_data = pair_data
                        other_pair_data = []

                    n = len(order_pair_data)  # number of pairs to order by similarity

                    dista = np.zeros((n,n))  # matrix of distances
                    different_normal_penalty = 1600
                    non_flip_distances = []
                    for i in range(0,n):
                        for j in range(i+1,n):
                            d = 0.0
                            d += (order_pair_data[i]['x']-order_pair_data[j]['x'])**2
                            d += (order_pair_data[i]['y']-order_pair_data[j]['y'])**2
                            d += 0.2*(order_pair_data[i]['z']-order_pair_data[j]['z'])**2
                            d += (order_pair_data[i]['maxgap']-order_pair_data[j]['maxgap'])**2
                            d += 0.02*(min(abs(order_pair_data[i]['angle_in_plane']-order_pair_data[j]['angle_in_plane']),360-abs(order_pair_data[i]['angle_in_plane']-order_pair_data[j]['angle_in_plane'])))**2
                            if order_pair_data[i]['normal_Z'] * order_pair_data[j]['normal_Z'] < 0:
                                # base flip with opposite normal vectors, that's a big deal!
                                d += different_normal_penalty
                                dista[i][j] = math.sqrt(d)
                            else:
                                d += 2*(order_pair_data[i]['normal_Z']-order_pair_data[j]['normal_Z'])**2
                                dista[i][j] = math.sqrt(d)
                                non_flip_distances.append(dista[i][j])

                            dista[j][i] = dista[i][j]

                    print("Finding order for %d basepairs" % n)

                    order = treePenalizedPathLength(dista,min(n,40))
                    order = standardOrder(dista,order)

                    if not compare_annotators:
                        # color entries on diagonal according to matching annotations
                        for i in range(0,n):
                            if "gap" in order_pair_data[i]['new_python_annotation']:
                                dista[i][i] = -2  # dark pink for bad gap
                            elif "min" in order_pair_data[i]['new_python_annotation'] or "max" in order_pair_data[i]['new_python_annotation'] or "angle" in order_pair_data[i]['new_python_annotation']:
                                dista[i][i] = -6  # sky blue for other cutoff problem
                            #elif 'hbond' in order_pair_data[i]['new_python_annotation']:
                            #    dista[i][i] = -9  # orange for bad h-bond
                            #elif len(order_pair_data[i]['matlab_annotation']) > 0 and len(order_pair_data[i]['python_annotation']) == 0:
                            #    dista[i][i] = -1  # reddish when matlab annotates but python_fr3d does not
                            #elif len(order_pair_data[i]['python_annotation']) > 0 and not "n" in order_pair_data[i]['python_annotation'] and len(order_pair_data[i]['matlab_annotation']) == 0:
                            #    dista[i][i] = -4  # purple when python_fr3d is true and Matlab is nothing
                            #elif order_pair_data[i]['matlab_annotation'].lower() == "n" + order_pair_data[i]['python_annotation'].lower():
                            #    dista[i][i] = -2  # dark pink when Matlab is near and python_fr3d is true

                    reorder_pairs = [order_pair_data[o] for o in order] + other_pair_data

                    reorder_dista = np.zeros((n,n))
                    if n > 10:
                        # cap dista just below the imact of the different_normal_penalty so that is not part of the value
                        max_dista = np.percentile(non_flip_distances,99)  # cap at 95th percentile to avoid outliers
                        Q["userMessage"].append("Heat map values capped at %0.4f to de-emphasize outliers and flipped bases" % max_dista)
                    else:
                        max_dista = float('inf')

                    for i in range(0,n):
                        for j in range(0,n):
                            reorder_dista[i][j] = min(max_dista,dista[order[i]][order[j]])

                    #print("Reordered instances and distance matrix")

                    """
                    Write the list of candidates in an HTML format that also shows
                    the coordinate window and a heat map of all-against-all distances.
                    """
                    if compare_annotators:
                        writeHTMLOutputCompareAnnotators(Q,reorder_pairs,reorder_dista)
                    else:
                        writeHTMLOutput(Q,reorder_pairs,reorder_dista)

            # print large minimum distances
            big_min_distance_list = sorted(big_min_distance_list)
            for b in big_min_distance_list:
                print("%8.2f %s %5s %5s" % b)


        # write a table of counts of each basepair family and base combination
        if compare_annotators:
            output_text = ""

            print("All four methods annotate at least one basepair in each of these %d structures:" % len(pdb_id_to_annotators.keys()))
            print(",".join(sorted(pdb_id_to_annotators.keys())))

            # show FR3D near
            header = "Family\tLW\tBase 1\tBase 2\tFR3D\tFR3D near\tRNAview\tRNAview near\tDSSR\tPDB\tcontacts\tMax / (1+Min)\tAll\tAll but FR3D\tAll but RNAview\tAll but DSSR\tAll but PDB\tFR3D RNAview\tFR3D DSSR\tFR3D PDB\tRNAview DSSR\tRNAview PDB\tDSSR PDB\tFR3D only\tRNAview only\tDSSR only\tPDB only"

            # leave out FR3D near, add FR3D-All, etc.
            header = "Fam\tLW\tBase 1\tBase 2\tStatus\tFR3D\tRNAview\tRNAview near\tDSSR\tPDB\tcontacts\tMax / (1+Min)\tFR3D-All\tRNAview-All\tDSSR-All\tPDB-All\tcontacts-All\tAll\tAll but FR3D\tAll but RNAview\tAll but DSSR\tAll but PDB\tAll but contacts\t3 agree\t2 agree\tFR3D only\tRNAview only\tDSSR only\tPDB only\tcontacts only\tLink\n"

            all_output_list = [header]
            print(header)

            for index,family in enumerate(LW12):
                output_list = []

                interaction = family

                table, combination_to_LSW_status = generate_LW_family_table(family)

                if interaction in count_annotations:
                    for base_combination in sorted(count_annotations[interaction].keys()):
                        b1,b2 = base_combination.split(",")
                        if combination_to_LSW_status[b1+","+b2] == 'Symmetric':
                            continue

                        output = [str(index+1),family,b1,b2,combination_to_LSW_status[b1+","+b2]]

                        if DNA:
                            filename = "DNA_compare_v6_%s_%s,%s_%s.html" % (interaction,b1,b2,resolution)
                        else:
                            filename = "compare_v6_%s_%s,%s_%s.html" % (interaction,b1,b2,resolution)

                        max_count = 0.0
                        min_count = float('inf')
                        for i,b in enumerate(count_annotations[interaction][base_combination]):
                            if i in [0,2,3,4,5,6]: # skip FR3D near
                                output.append(str(b))
                            if i in [0,2,4,5,6]:  # skip FR3D near, RNAview near
                                max_count = max(max_count,b)
                                min_count = min(min_count,b)

                        if max_count > 0:
                            output.append("%0.2f" % (float(max_count)/(1.0+min_count)))
                            all_count = count_annotations_by_group[interaction][base_combination]["11111"]

                            for i,b in enumerate(count_annotations[interaction][base_combination]):
                                if i in [0,2,4,5,6]:  # skip FR3D near, RNAview near
                                    output.append(str(b-all_count))

                            output.append(str(count_annotations_by_group[interaction][base_combination]["11111"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["01111"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["10111"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["11011"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["11101"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["11110"]))

                            num_agree = defaultdict(int)
                            for mask in count_annotations_by_group[interaction][base_combination].keys():
                                num = mask.count('1')
                                num_agree[num] += 1

                            output.append(str(num_agree[3]))
                            output.append(str(num_agree[2]))

                            output.append(str(count_annotations_by_group[interaction][base_combination]["10000"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["01000"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["00100"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["00010"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["00001"]))
                            output.append("http://rna.bgsu.edu/experiments/annotations/" + filename)
                            output_list.append(output)

                output_list = sorted(output_list)             # sort by base combination

                for output in output_list:
                    output_text = "\t".join(output)
                    all_output_list.append(output_text+"\n")
                    print(output_text)

                if DNA:
                    table_filename = os.path.join(OUTPUTPATH,'DNA_table_v6_%s.txt' % resolution)
                else:
                    table_filename = os.path.join(OUTPUTPATH,'table_v6_%s.txt' % resolution)

                with open(table_filename,'w') as f:
                    f.writelines(all_output_list)

                print('Saved %s table in %s' % (resolution,table_filename))

        print("Saved figures in %s/plots" % outputNAPairwiseInteractions)
        print("Saved HTML files in %s" % OUTPUTPATH)
        print()