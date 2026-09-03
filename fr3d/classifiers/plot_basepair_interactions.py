"""
plot-basepair-interactions.py reads a data file and plots points to represent
base pairing interactions.
Points are colored according to python_fr3d and NAPAIR annotations, maybe

HTML output is in C:/Users/zirbel/Documents/PythonFR3D/output
PNG output is in C:/Users/zirbel/Documents/NAPairwiseInteractions

cd "C:/Users/zirbel/Documents/PythonFR3D/output"
python -m SimpleHTTPServer

"""

from collections import defaultdict
import json
import math
import matplotlib.pyplot as plt
import numpy as np
import os
import pickle
import random
import requests
from scipy.stats import trim_mean
import sys
import time

# import the version of urlretrieve appropriate to the Python version
if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
    from urllib import urlopen
else:
    from urllib.request import urlretrieve as urlretrieve
    from urllib.request import urlopen

from class_limits_2024 import nt_nt_cutoffs
from NA_pairwise_interactions import map_PDB_list_to_PDB_IFE_dict
from NA_pairwise_interactions import reverse_edges
from NA_protein_annotation import make_heatmap_text
from draw_residues import draw_base

from fr3d.localpath import outputNAPairwiseInteractions
from fr3d.localpath import fr3d_pickle_path

from fr3d.search.file_reading import readPDBDatafile
from fr3d.search.file_reading import readUnitAnnotations

# very specific directories, for data files in specific formats.  Hard for others to run!
dssr_basepair_path = 'C:/Users/zirbel/Documents/PythonFR3D/data/pairs_dssr'
rnaview_basepair_path = 'C:/Users/zirbel/Documents/PythonFR3D/data/pairs_rnaview'
mcannotate_basepair_path = 'C:/Users/zirbel/Documents/PythonFR3D/data/pairs_mcannotate'
pdb_basepair_path = 'C:/Users/zirbel/Documents/PythonFR3D/data/pairs_pdb'
napair_basepair_path = 'C:/Users/zirbel/Documents/PythonFR3D/data/pairs_napair'

OUTPUTPATH = "C:/Users/zirbel/Documents/PythonFR3D/output/"

#Updated modified nucleotide mappings from atom_mappings_refined.txt
from fr3d.modified.mapping import modified_base_atom_list,parent_atom_to_modified,modified_atom_to_parent,modified_base_to_parent

from fr3d.ordering.orderBySimilarity import treePenalizedPathLength
from fr3d.ordering.orderBySimilarity import standardOrder

from fr3d.classifiers.NA_pairwise_interactions import get_parent
from fr3d.classifiers.NA_pairwise_interactions import get_parent_as_RNA

JS1 = '  <script src="./js/JSmol.min.nojq.js"></script>'
JS2 = '  <script src="./js/jquery.jmolTools.bp.js"></script>'               # special version, superimpose first base
JS3 = '  <script src="./js/imagehandlinglocal.js"></script>'
JS4 = '<script src="./js/jmolplugin.js" type="text/javascript"></script>'
JS5 = '<script type="text/javascript" src="./js/heatmap.js"></script>'
JS6 = '<script src="./js/sort.table.js" type="text/javascript"></script>'
JS7 = '<script src="./js/shift.click.checkbox.js" type="text/javascript"></script>'

TEMPLATEPATH = '../search/'

near_discrepancy_cutoff = 2.0
unknown_bases = set()

import sys
import gc

def find_variable_names(obj):
    """
    Find variable names in the global scope that reference the given object.
    """
    variable_names = []
    for var_name, var_value in globals().items():
        if var_value is obj:
            variable_names.append(var_name)
    return variable_names

def show_memory_usage():
    """
    Display all variables in memory and their sizes.
    """
    # Get all objects currently in memory
    all_objects = gc.get_objects()

    print(f"{'Type':<20} {'Size (bytes)':<15} {'Variable':<30}")
    print("-" * 65)

    # Iterate through all objects
    for obj in all_objects:
        try:
            # Get the variable name, type, and size
            obj_type = type(obj).__name__
            obj_size = sys.getsizeof(obj)

            if obj_size > 100000:
                # Print the details
                var_names = find_variable_names(obj)

                print("%12s %10d %s" % (obj_type,obj_size,var_names))
        except Exception as e:
            # Handle any errors (e.g., inaccessible objects)
            print(f"Error processing object: {e}")


def print_dictionary(datapoint):
    for key,value in sorted(datapoint.items()):
        if type(value) == dict:
            print(key,' is a dictionary:')
            for k,v in sorted(value.items()):
                print("  %s = %s" % (k,v))
        else:
            print(key,value)


def get_parent_pair(pair):
    """
    Convert DNA and modified nucleotides to A, C, G, U
    Return b1 and b2 as None if there is no match
    """

    if len(pair) < 2:
        print("get_parent_pair not sure what to do with ",pair)
        input("Press Enter")

    if not pair[0] or not pair[1]:
        print("get_parent_pair not sure what to do with ",pair)
        input("Press Enter")

    s1 = pair[0].split("|")[3]
    s2 = pair[1].split("|")[3]

    b1 = get_parent_as_RNA(s1)
    b2 = get_parent_as_RNA(s2)

    if b1 and b2:
        return b1, b2, b1+","+b2
    else:
        return b1, b2, None


def format_resolution(data):

    r = data['resolution']
    try:
        s = "%0.2f" % r
        # remove trailing 0
        if s[-1] == "0":
            s = s[0:(len(s)-1)]
    except:
        if 'NMR' in data['method']:
            s = "NMR"
        else:
            s = "NA"
    return s


def format_distance_atoms(d_atoms,atom_set,nt1_seq,nt2_seq,reverse_bases):

    # print('format_distance_atoms')
    # print('Given d_atoms',d_atoms)
    # print('Given atom_set',atom_set)

    if nt1_seq == nt2_seq:
        s1 = nt1_seq + "1"
        s2 = nt2_seq + "2"
    else:
        s1 = nt1_seq
        s2 = nt2_seq

    f = d_atoms.split("-")

    distance_atoms = ""

    if reverse_bases:
        if atom_set[3] == '12':
            distance_atoms = "%s(%s)..%s(%s)" % (f[1],s2,f[0],s1)
        else:
            distance_atoms = "%s(%s)..%s(%s)" % (f[0],s2,f[1],s1)
    else:
        if atom_set[3] == '12':
            distance_atoms = "%s(%s)..%s(%s)" % (f[0],s1,f[1],s2)
        else:
            distance_atoms = "%s(%s)..%s(%s)" % (f[1],s1,f[0],s2)

    # print('Produced formatted distance_atoms',distance_atoms)

    return distance_atoms


def format_angle_atoms(a_atoms,atom_set,nt1_seq,nt2_seq,reverse_bases):

    # print('format_angle_atoms')
    # print('Given a_atoms',a_atoms)
    # print('Given atom_set',atom_set)

    if nt1_seq == nt2_seq:
        s1 = nt1_seq + "1"
        s2 = nt2_seq + "2"
    else:
        s1 = nt1_seq
        s2 = nt2_seq

    f = a_atoms.split("-")

    angle_atoms = ""

    if reverse_bases:
        if atom_set[3] == '12':
            angle_atoms = "%s(%s)..%s(%s)-%s(%s)" % (f[2],s2,f[1],s1,f[0],s1)
        else:
            angle_atoms = "%s(%s)-%s(%s)..%s(%s)" % (f[0],s2,f[1],s2,f[2],s1)
    else:
        if atom_set[3] == '12':
            angle_atoms = "%s(%s)-%s(%s)..%s(%s)" % (f[0],s1,f[1],s1,f[2],s2)
        else:
            angle_atoms = "%s(%s)..%s(%s)-%s(%s)" % (f[2],s1,f[1],s2,f[0],s2)

    # print('Produced formatted angle_atoms',angle_atoms)

    return angle_atoms


def get_dist2_angle_penalty(interaction_to_atom_sets,candidate,target_distance=None,target_angle=None):
    """
    Loop over hydrogen bonds and find the second-worst distance
    Also calculate angles
    Calculate "badness" relative to a target distance and angle, if available
    """

    # list and score hydrogen bonds
    badness = []
    dist2_list = []           # list of scores for distances too long
    hDistAngle_list = []      # list of scores for distances and/or angles out of specifications
    z_score = 0.0
    z_count = 0
    extra_penalty = 0.0

    candidatelist = ""

    # todo: when more than one interaction for this base combination,
    # only show the hydrogen bonds that match the best
    # this is a problem when there are pairs that are not true FR3D pairs
    # or when there are cWW and cWWa pairs on the same page

    distances_this_interaction = 0
    ok_distances_this_interaction = 0

    for interaction in interaction_to_atom_sets.keys():
        # loop over hydrogen bond atom sets for this interaction
        for atom_set in interaction_to_atom_sets[interaction]:
            if "atom_set_to_results" in candidate and atom_set in candidate["atom_set_to_results"]:
                if len(interaction_to_atom_sets.keys()) > 1 and "python_annotation" in candidate and not "n" in candidate["python_annotation"] and not candidate["python_annotation"] == interaction:
                    candidatelist += "<td></td><td></td>"  # don't show hydrogen bond information
                else:
                    if target_distance:
                        td = target_distance[atom_set]
                    else:
                        td = 2.9

                    if target_angle:
                        ta = target_angle[atom_set]
                    else:
                        ta = 120.0

                    result = candidate["atom_set_to_results"][atom_set]
                    candidatelist += "<td>%0.2f</td>" % (result["donor_acceptor_distance"])  # hydrogen bond distance
                    candidatelist += "<td>%0.1f</td>" % (result["heavy_donor_acceptor_angle"])  # hydrogen bond angle
                    z_score += abs(result["donor_acceptor_distance"] - td) / 0.1
                    z_score += abs(result["heavy_donor_acceptor_angle"] - ta) / 10.0
                    z_count += 2
                    extra_penalty += 3*max(0.0, abs(result["donor_acceptor_distance"]-td) - 1.0)  # extra penalty when more than 1A from target distance

                    bad = max(0.0, abs(result["donor_acceptor_distance"]-2.9) - 0.1)    # how far from 2.9, penalize when above 0.2
                    bad += max(0.0, abs(result["heavy_donor_acceptor_angle"] - 120) - 10)/20.0  # how far from 120, penalize when more than 20 degrees away
                    badness.append(bad)

                    bad = max(0.0, result["donor_acceptor_distance"]-3.5)    # how far above 3.5A
                    dist2_list.append(result["donor_acceptor_distance"])     # distance
                    bad += (max(0.0, abs(result["heavy_donor_acceptor_angle"] - 120) - 20))/20.0  # more than 20 degrees from 120
                    hDistAngle_list.append(bad)

                    distances_this_interaction += 1
                    if candidate["python_annotation"].replace("n","") == interaction:
                        if "C" in result["donor_acceptor_atoms"]:
                            if result["donor_acceptor_distance"] < 4.2:
                                ok_distances_this_interaction += 1
                        elif result["donor_acceptor_distance"] < 3.6:
                            ok_distances_this_interaction += 1

            else:
                candidatelist += "<td></td><td></td>"  # no hydrogen bond information
                z_score = 99.99  # bad
                z_count = 1
                badness.append(9.99)
                dist2_list.append(99.99)
                hDistAngle_list.append(99.99)

    # if len(badness) > 1:
    #     badness = sorted(badness)
    #     candidatelist += "<td>%0.2f</td>" % (badness[1])  # second worst badness
    # elif len(badness) == 1:
    #     candidatelist += "<td>%0.2f</td>" % (badness[0])  # the only badness
    # else:
    #     candidatelist += "<td>9.99</td>"  # pretty bad

    if z_count == 0:
        z_score = 99.99
        z_count = 1

    ok_hbonds = False
    if distances_this_interaction == 1 and ok_distances_this_interaction == 1:
        ok_hbonds = True
    elif distances_this_interaction >= 2 and ok_distances_this_interaction >= 2:
        ok_hbonds = True

    return dist2_list, hDistAngle_list, extra_penalty, candidatelist, z_score, z_count, ok_hbonds


def writeHTMLOutput(Q,candidates,interaction_to_atom_sets,distance_angle_messages,allvsallmatrix=np.empty( shape=(0, 0) )):
    """
    Write the list of candidates in an HTML format that also shows
    the coordinate window and a heat map of all-against-all distances.
    Write the heatmap data in an efficient way.
    Do not show x,y,z parameters or scatterplots, but do show annotations from different programs.
    """

    pagetitle = "%s" % Q['name']

    if show_fr3d_parameters:
        htmlfilename = "fr3d_%s_%s" % (VERSION,Q['name'].replace(" ","_"))
    else:
        htmlfilename = "compare_%s_%s" % (VERSION,Q['name'].replace(" ","_"))
    htmlfilename = htmlfilename.replace('compare_%s_DNA' % VERSION,'DNA_compare_%s' % VERSION)  # awkward, fragile, but works

    # print('htmlfilename',htmlfilename)
    # print('Q name',Q['name'])

    base_combination = Q['name'].split()[1]

    # make the table of instances
    candidatelist = '<table id="instances" style="white-space:nowrap;">\n'

    numPositions = 2

    # write header line, make columns sortable, numeric=true, alpha=false
    candidatelist += '<tr><th onclick="sortTable(0,\'instances\',\'numeric\')">S.</th>'
    candidatelist += '<th onclick="sortTable(1,\'instances\',\'checkbox\')">Show</th>'
    candidatelist += '<th onclick="sortTable(2,\'instances\',\'numeric\')">Res.</th>'

    for j in range(0,numPositions):
        candidatelist += '<th onclick="sortTable(%d,\'instances\',\'alpha\')">Position %d</th>' % (j+3,j+1)
    candidatelist += '<th onclick="sortTable(5,\'instances\',\'alpha\')">FR3D</th>'
    candidatelist += '<th onclick="sortTable(6,\'instances\',\'alpha\')">RNAview</th>'
    candidatelist += '<th onclick="sortTable(7,\'instances\',\'alpha\')">MCA</th>'
    candidatelist += '<th onclick="sortTable(8,\'instances\',\'alpha\')">DSSR</th>'
    candidatelist += '<th onclick="sortTable(9,\'instances\',\'alpha\')">PDB</th>'
    candidatelist += '<th onclick="sortTable(10,\'instances\',\'alpha\')">NAPAIR</th>'
    candidatelist += '<th onclick="sortTable(11,\'instances\',\'numeric\')">#</th>'
    candidatelist += '<th onclick="sortTable(12,\'instances\',\'alpha\')">Which</th>'
    candidatelist += '<th onclick="sortTable(13,\'instances\',\'alpha\')">Glyco 1</th>'
    candidatelist += '<th onclick="sortTable(14,\'instances\',\'alpha\')">Glyco 2</th>'

    if show_fr3d_parameters:
        candidatelist += '<th onclick="sortTable(15,\'instances\',\'alpha\')">Newest</th>'
        candidatelist += '<th onclick="sortTable(16,\'instances\',\'alpha\')">FR3D detail</th>'
        candidatelist += '<th onclick="sortTable(17,\'instances\',\'alpha\')">Subcat</th>'
        candidatelist += '<th onclick="sortTable(18,\'instances\',\'numeric\')">cut dist</th>'
        candidatelist += '<th onclick="sortTable(19,\'instances\',\'numeric\')">x</th>'
        candidatelist += '<th onclick="sortTable(20,\'instances\',\'numeric\')">y</th>'
        candidatelist += '<th onclick="sortTable(21,\'instances\',\'numeric\')">z</th>'
        candidatelist += '<th onclick="sortTable(22,\'instances\',\'numeric\')">r</th>'
        candidatelist += '<th onclick="sortTable(23,\'instances\',\'numeric\')">maxgap</th>'
        candidatelist += '<th onclick="sortTable(24,\'instances\',\'numeric\')">angle</th>'
        candidatelist += '<th onclick="sortTable(25,\'instances\',\'numeric\')">normal_z</th>'
        candidatelist += '<th onclick="sortTable(26,\'instances\',\'numeric\')">hmindist</th>'
        current_column = 27   # next column after the one list above, keep track manually
        if 'css' in Q['name'].lower():
            candidatelist += '<th onclick="sortTable(27,\'instances\',\'alpha\')">SR12</th>'
            candidatelist += '<th onclick="sortTable(28,\'instances\',\'alpha\')">SR21</th>'
            current_column += 2
    else:
        current_column = 15   # next column after the one list above, keep track manually

    # loop over true interaction(s) in this group, like cWW and cWWa
    if len(candidates) > 0:
        for interaction in interaction_to_atom_sets.keys():
            # loop over hydrogen bond atom sets for this interaction
            for atom_set in interaction_to_atom_sets[interaction]:
                # print("Making headers for atom set %s" % (str(atom_set)))

                # look up formatted listings of atoms
                d_atoms = Q['d_atoms'][atom_set]
                a_atoms = Q['a_atoms'][atom_set]

                if len(interaction_to_atom_sets) > 1:
                    # more than one true interaction, so list the interaction
                    d_atoms += " " + interaction
                    a_atoms += " " + interaction

                candidatelist += '<th onclick="sortTable(%d,\'instances\',\'numeric\')">%s</th>' % (current_column,d_atoms)
                current_column += 1
                candidatelist += '<th onclick="sortTable(%d,\'instances\',\'numeric\')">%s</th>' % (current_column,a_atoms)
                current_column += 1

    # make headers for various measures of the quality of the hydrogen bonding
    #candidatelist += '<th onclick="sortTable(%d,\'instances\',\'numeric\')">SBHBB</th>' % (current_column)  # second best hydrogen bond badness
    candidatelist += '<th onclick="sortTable(%d,\'instances\',\'numeric\')">Dist2</th>' % (current_column)  #
    current_column += 1
    if show_fr3d_parameters:
        candidatelist += '<th onclick="sortTable(%d,\'instances\',\'alpha\')">h OK</th>' % (current_column)  #
        current_column += 1

    candidatelist += '<th onclick="sortTable(%d,\'instances\',\'numeric\')">hDistAngle</th>' % (current_column)  #
    current_column += 1
    candidatelist += '<th onclick="sortTable(%d,\'instances\',\'numeric\')">hbond1</th>' % (current_column)  #
    current_column += 1
    candidatelist += '<th onclick="sortTable(%d,\'instances\',\'numeric\')">hbond2</th>' % (current_column)  #
    current_column += 1
    candidatelist += '<th onclick="sortTable(%d,\'instances\',\'numeric\')">hbond3</th>' % (current_column)  #
    current_column += 1

    # pdata["atom_sets"] = sorted(datapoint["LW_to_atom_sets"][interaction])
    # pdata["atom_set_to_results"] = datapoint["atom_set_to_results"]

    if not Q['DNA']:
       candidatelist += '<th>3D structures</th>'
       candidatelist += '<th>Sequence alignments</th>'
       candidatelist += '<th>View</th>'

    candidatelist += "</tr>\n"

    first_listed_base = 'A'
    for i in range(0,len(candidates)):
        candidate = candidates[i]
        first_listed_base = candidate['unit_id_1'].split("|")[3]
        if first_listed_base in ['A','C','G','U','DA','DC','DG','DT']:
            break

    if first_listed_base == base_combination.split(",")[0] or show_fr3d_parameters:
        id_1 = 'unit_id_1'
        id_2 = 'unit_id_2'
    else:
        # switch order of bases to match the order they appear in the filename
        id_1 = 'unit_id_2'
        id_2 = 'unit_id_1'

    # print('First listed base is %s, base_combination is %s' % (first_listed_base,base_combination))

    # write one row for each candidate
    for i in range(0,len(candidates)):

        candidate = candidates[i]

        candidatelist += '<tr><td>'+str(i+1)+'.</td>'
        candidatelist += '<td><label><input type="checkbox" id="'+str(i)+'" class="jmolInline" data-coord="'
        candidatelist += candidate[id_1] + ',' + candidate[id_2]
        candidatelist += '">&nbsp</label></td>'

        pdb_id = candidate['unit_id_1'][0:4]

        if pdb_id in Q["PDB_data_file"]:
            candidatelist += "<td>%s</td>" % format_resolution(Q["PDB_data_file"][pdb_id])
        else:
            candidatelist += "<td>NA</td>"

        # write unit ids
        candidatelist += "<td>"+candidate[id_1]+"</td>"
        candidatelist += "<td>"+candidate[id_2]+"</td>"

        if id_1 == 'unit_id_1':
            candidatelist += "<td>%s</td>" % candidate['python_annotation']  # python FR3D, with w, h, s
            candidatelist += "<td>%s</td>" % candidate['rnaview_annotation']  # rnaview
            candidatelist += "<td>%s</td>" % candidate['mcannotate_annotation']  # mcannotate
            candidatelist += "<td>%s</td>" % candidate['dssr_annotation']  # dssr
            candidatelist += "<td>%s</td>" % candidate['pdb_annotation']  # pdb
            candidatelist += "<td>%s</td>" % candidate['napair_annotation']  # napair
        else:
            # candidatelist += "<td>%s</td>" % reverse_edges(candidate['python_annotation'].replace("a",""))  # python FR3D, with w, h, s
            candidatelist += "<td>%s</td>" % reverse_edges(candidate['python_annotation'])  # python FR3D, with w, h, s
            candidatelist += "<td>%s</td>" % reverse_edges(candidate['rnaview_annotation'])  # rnaview
            candidatelist += "<td>%s</td>" % reverse_edges(candidate['mcannotate_annotation'])  # mcannotate
            candidatelist += "<td>%s</td>" % reverse_edges(candidate['dssr_annotation'])  # dssr
            candidatelist += "<td>%s</td>" % candidate['pdb_annotation']  # pdb
            candidatelist += "<td>%s</td>" % reverse_edges(candidate['napair_annotation'])  # napair

        candidatelist += "<td>%s</td>" % candidate['annotator_count']  # number that agree
        candidatelist += "<td>%s</td>" % candidate['annotator']  # who annotates

        candidatelist += "<td>%s</td>" % candidate['glycosidic1']  # unit 1 glycosidic bond
        candidatelist += "<td>%s</td>" % candidate['glycosidic2']  # unit 2 glycosidic bond

        if show_fr3d_parameters:
            candidatelist += "<td>%s</td>" % candidate['new_python_annotation']  # python
            candidatelist += "<td>%s</td>" % candidate['new_fr3d_detail']  # after checks in this program
            candidatelist += "<td>%d</td>" % candidate['new_basepair_subcat']  # python subcategory
            candidatelist += "<td>%0.4f</td>" % candidate['cut_dist']  #
            candidatelist += "<td>%0.2f</td>" % candidate['x']  #
            candidatelist += "<td>%0.2f</td>" % candidate['y']  #
            candidatelist += "<td>%0.2f</td>" % candidate['z']  #
            candidatelist += "<td>%0.2f</td>" % (math.sqrt(candidate['x']**2 + candidate['y']**2))
            candidatelist += "<td>%0.2f</td>" % candidate['maxgap']  #
            candidatelist += "<td>%0.2f</td>" % candidate['angle_in_plane']  #
            candidatelist += "<td>%0.2f</td>" % candidate['normal_Z']  #

            if not 'heavy_min_distance' in candidate:
                candidate['heavy_min_distance'] = -1.0

            candidatelist += "<td>%0.2f</td>" % candidate['heavy_min_distance']  #

            if 'css' in Q['name'].lower():
                if 'sugar_ribose' in candidate:
                    candidatelist += "<td>%s</td>" % candidate['sugar_ribose']
                else:
                    candidatelist += "<td></td>"

                if 'r_sugar_ribose' in candidate:
                    candidatelist += "<td>%s</td>" % candidate['r_sugar_ribose']
                else:
                    candidatelist += "<td></td>"

        dist2_list, hDistAngle_list, extra_penalty, cl, z_score, z_count, ok_hbonds = get_dist2_angle_penalty(interaction_to_atom_sets,candidate,Q['target_distance'],Q['target_angle'])

        candidatelist += cl

        # hDist is zero when there are two good hydrogen bonds, and non-zero to the extent that the second hydrogen bond is bad
        if len(dist2_list) >= 2:
            dist2 = sorted(dist2_list)[1]   # second shortest distance
            candidatelist += "<td>%0.2f</td>" % (sorted(dist2_list)[1])   # second shortest distance
        elif len(dist2_list) == 1:
            dist2 = dist2_list[0]   # the only distance
            candidatelist += "<td>%0.2f</td>" % (dist2_list[0])   # the only distance
        else:
            dist2 = 99.9
            candidatelist += "<td></td>"

        if show_fr3d_parameters:
            if candidate['ok_hbonds']:
                candidatelist += "<td>OK</td>"
            else:
                candidatelist += "<td>No</td>"

        # hDistAngle is zero when there are two good hydrogen bonds, and non-zero to the extent that the second hydrogen bond is bad
        if len(hDistAngle_list) >= 2:
            candidatelist += "<td>%0.2f</td>" % (sorted(hDistAngle_list)[1])   # second best badness
        elif len(hDistAngle_list) == 1:
            candidatelist += "<td>%0.2f</td>" % (hDistAngle_list[0])   # the only badness
        else:
            candidatelist += "<td></td>"

        # hbond1 is the sum of z-scores
        candidatelist += "<td>%0.2f</td>" % (z_score/z_count)  # total z score

        # hbond2 is the sum of z-scores plus extra penalty for being more than 1A from target distance
        z_score += extra_penalty
        z_count += 1
        candidatelist += "<td>%0.2f</td>" % (z_score/z_count)  # total z score

        # hbond3 has an additional penalty for the maxgap
        z_score += 4*max(0.0,candidate['maxgap']-0.8)
        z_count += 1

        candidatelist += "<td>%0.2f</td>" % (z_score/z_count)  # total z score

        if not Q['DNA']:
            candidatelist += '<td><a href="https://rna.bgsu.edu/correspondence/comparison?selection=%s,%s&exp_method=all&resolution=3.0" target="_blank" rel="noopener noreferrer">Compare 3D</a></td>' % (candidate[id_1],candidate[id_2])
            candidatelist += '<td><a href="https://rna.bgsu.edu/correspondence/variability?id=%s,%s&format=unique" target="_blank" rel="noopener noreferrer">Compare sequence</a></td>' % (candidate[id_1],candidate[id_2])
            candidatelist += '<td><a href="https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" target="_blank" rel="noopener noreferrer">View</a></td>' % (candidate[id_1],candidate[id_2])

        candidatelist += '</tr>\n'
    candidatelist += '</table>\n'

    # make labels for each candidate in the heatmap
    labels = [c[id_1] for c in candidates]
    s = len(allvsallmatrix)
    labels = labels[0:s]
    discrepancy_text = make_heatmap_text(allvsallmatrix,labels)

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

    if VERSION in ['v9']:
        if show_fr3d_parameters:
            prefix = "fr3d_%s_" % VERSION
        else:
            prefix = "compare_%s_" % VERSION
        if Q['DNA']:
            prefix += "DNA_"
        if Q['only_modified']:
            prefix += "modified_"
        table, combination_to_LSW_status = generate_LW_family_table(Q['LW'],prefix=prefix,resolutions=Q['link_blocks'])
    elif show_fr3d_parameters:
        prefix = "fr3d_%s_" % VERSION
        if Q['DNA']:
            prefix += "DNA_"
        if Q['only_modified']:
            prefix += "modified_"
        table, combination_to_LSW_status = generate_LW_family_table(Q['LW'],prefix=prefix)
    else:
        if Q['DNA']:
            table, combination_to_LSW_status = generate_LW_family_table(Q['LW'],prefix='DNA_compare_%s_' % VERSION)
        else:
            table, combination_to_LSW_status = generate_LW_family_table(Q['LW'],prefix='compare_%s_' % VERSION)

    template = template.replace("###seeCSVOutput###",table)

    description = "<br>Columns of the table show candidate number in similarity order, checkbox to display coordinates or not, interacting nucleotides, annotations from different programs, number of programs which annotate, and a string that tells which programs annotate.  Click on column headings to sort.<br>"

    if show_fr3d_parameters:
        # put in the scatterplot
        description += '\n<br><img src="%s"><br>\n' % Q['figure_img_src']

    description += "<br>".join(distance_angle_messages)

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
    template = template.replace("###JS6###",JS6)
    template = template.replace("###JS7###",JS7)

    refresh = ""
    if "reloadOutputPage" in Q and Q["reloadOutputPage"]:
        REFRESHTIME = 2
        refresh = '<meta http-equiv="refresh" content="%d">' % REFRESHTIME
    template = template.replace("###REFRESH###",refresh)

    if np.size(allvsallmatrix) > 0:
        template = template.replace("###JS5###",JS5)    # include heatmap.js code
        template = template.replace("###DISCREPANCYDATA###",discrepancy_text)
    else:
        template = template.replace("###JS5###","")    # do not display a heat map
        template = template.replace("###DISCREPANCYDATA###","")

    outputfilename = os.path.join(OUTPUTPATH,htmlfilename+".html")

    print("plot_basepair_interactions is writing to %s" % outputfilename)

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

            url= "https://rna.bgsu.edu/rna3dhub/display3D/unitid/" + key[0] + "," +  key[1]

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
        datapoint['basepair'] = 'N/A x'
        datapoint['basepair_subcategory'] = None
    if not 'y' in datapoint:
        #print('Adding y to %s' % datapoint)
        datapoint['y'] = 0.000
        datapoint['basepair'] = 'N/A y'
        datapoint['basepair_subcategory'] = None
    if not 'z' in datapoint:
        #print('Adding z to %s' % datapoint)
        datapoint['z'] = 0.000
        datapoint['basepair'] = 'N/A z'
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
        datapoint['basepair'] = 'N/A n'
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

    fields = ['x','y','z','r','gap12','angle_in_plane','normal_Z']

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

                            if "radiusmin" in limits:
                                ax.plot([min([t1,t2,t3,t4]),max([t1,t2,t3,t4])],[limits["radiusmin"],limits["radiusmin"]],color[cc])

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
                                # for example, xmin = -75, xmax = -55
                                # can happen when one subcategory wraps and the other does not
                                ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            elif angle_out_of_order and xmax < 0:
                                # for example, xmin = 260 and xmax = -50
                                ax.plot([xmin-360,xmax,xmax,xmin-360,xmin-360],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            elif angle_out_of_order and xmin > 0:
                                # for example, xmin = 240 and xmax = 270
                                ax.plot([xmin-360,xmax-360,xmax-360,xmin-360,xmin-360],[ymin,ymin,ymax,ymax,ymin],color[cc])
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
                                ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            elif angle_out_of_order and xmax < 0:
                                # for example, xmin = 260 and xmax = -50
                                ax.plot([xmin-360,xmax,xmax,xmin-360,xmin-360],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            elif angle_out_of_order and xmin > 0:
                                # for example, xmin = 240 and xmax = 270
                                ax.plot([xmin-360,xmax-360,xmax-360,xmin-360,xmin-360],[ymin,ymin,ymax,ymax,ymin],color[cc])
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
                                ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            elif angle_out_of_order and xmax < 0:
                                # for example, xmin = 260 and xmax = -50
                                ax.plot([xmin-360,xmax,xmax,xmin-360,xmin-360],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            elif angle_out_of_order and xmin > 0:
                                # for example, xmin = 240 and xmax = 270
                                ax.plot([xmin-360,xmax-360,xmax-360,xmin-360,xmin-360],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            elif xmin < xmax:
                                # angle does not wrap around 270 degrees
                                ax.plot([xmin,xmax,xmax,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                            else:
                                # angle wraps around 270 degrees
                                ax.plot([xmin,xmax+360,xmax+360,xmin,xmin],[ymin,ymin,ymax,ymax,ymin],color[cc])
                                # ax.plot([270,xmin,xmin,270],[ymin,ymin,ymax,ymax],color[cc])
                                # ax.plot([-90,xmax,xmax,-90],[ymin,ymin,ymax,ymax],color[cc])

                            cc += 1


def generate_LW_family_table(LW,prefix='',resolutions=['1.5A','2.0A','2.5A','3.0A'],DNA=False):
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
        output += "<tr><td>2.0A</td><td></td><td></td><td></td><td></td></tr>\n"
    else:
        bases = ['A','C','G','U']
        if len(resolutions) == 4:
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


def map_entry_to_index(fields):
    field_to_index = {}
    for c, field in enumerate(fields):
        field_to_index[field] = c
    return field_to_index


def load_napair_chains(filenames = ['RS_filtered_chains_25','RS_filtered_chains_35']):
    # read RNA and DNA chain lists from NAPAIR, process chains in the union of the two
    chain_set = set()
    pdb_set = set()
    for filename in filenames:
        directory = os.path.join('C:/Users/zirbel/Documents/PythonFR3D/data/pairs_napair',filename)
        # get directory listing
        files = os.listdir(directory)

        # loop over files in the directory
        for file in files:
            # use the filename to get the pdb id and chain
            # 1a9n_Q_filtered.csv

            fields = file.split("_")
            if len(fields) == 3:
                pdb = fields[0].upper()
                pdb_set.add(pdb)
                chain = "%s|1|%s" % (pdb,fields[1])
                chain_set.add(chain)

            if False:
                # read the file to extract every chain that appears, which introduces additional chains
                with open(os.path.join(directory,file),'rt') as f:
                    lines = f.readlines()

                for line in lines:
                    # pdbid,model,family,class,full_subclass,label_alt_id1,label_comp_id1,label_asym_id1,label_entity_id1,label_seq_id1,pdbx_PDB_ins_code1,auth_seq_id1,auth_comp_id1,auth_asym_id1,symmetry_operation1,symmetry_id1,label_alt_id2,label_comp_id2,label_asym_id2,label_entity_id2,label_seq_id2,pdbx_PDB_ins_code2,auth_seq_id2,auth_comp_id2,auth_asym_id2,symmetry_operation2,symmetry_id2,napasco,rmsd,curated_file,knn_metric,coplanarity_angle,coplanarity_shift1,coplanarity_shift2,coplanarity_edge_angle1,coplanarity_edge_angle2,C1_C1_yaw1,C1_C1_pitch1,C1_C1_roll1,C1_C1_yaw2,C1_C1_pitch2,C1_C1_roll2,hb_0_length,hb_0_donor_angle,hb_0_acceptor_angle,hb_0_OOPA1,hb_0_OOPA2,hb_1_length,hb_1_donor_angle,hb_1_acceptor_angle,hb_1_OOPA1,hb_1_OOPA2,hb_2_length,hb_2_donor_angle,hb_2_acceptor_angle,hb_2_OOPA1,hb_2_OOPA2,hb_3_length,hb_3_donor_angle,hb_3_acceptor_angle,hb_3_OOPA1,hb_3_OOPA2
                    # 1b72,1,cWW,cWW_G-C,cWW_G-C_1,,DG,A,1,8,,8,DG,D,,1,,DC,B,2,14,,34,DC,E,,1,88.21301512714346,0.14,6btf_P_DG_4_1_555_T_DC_13_1_555,2.347,168.7993,0.11854984,-0.13904174,10.2320385,10.673146,71.98174,8.523693,170.57051,72.07385,6.348266,168.99513,2.8707104,123.48988,119.513504,5.213443,2.1907403,2.8906677,116.003624,121.427605,-2.7632103,-5.738376,2.852179,122.21896,119.42502,-10.818403,-13.972753,,,,,

                    fields = line.split(",")

                    if fields[0] == 'pdbid':
                        # header line
                        header_to_index = map_entry_to_index(fields)

                        chain1index = header_to_index['auth_asym_id1']
                        chain2index = header_to_index['auth_asym_id2']
                    else:
                        pdb = fields[0].upper()
                        pdb_set.add(pdb)

                        chain = "%s|1|%s" % (pdb,fields[chain1index])

                        # if not chain in chain_set:
                        #     print('NAPAIR annotates',file,pdb,chain)

                        chain_set.add(chain)

                        chain = "%s|1|%s" % (pdb,fields[chain2index])
                        chain_set.add(chain)

        print('Have %4d distinct PDB ids and %4d distinct chain ids after %s' % (len(pdb_set),len(chain_set),filename))

    all_pdb_ids = sorted(pdb_set)
    chain_list = sorted(chain_set)

    print('Have %4d distinct PDB ids and %4d distinct chain ids' % (len(all_pdb_ids),len(chain_list)))

    return all_pdb_ids, chain_list


def load_napair_chains_may_2026(filenames = ['RS35_DNA_c2.csv','RS35_RNA_c2.csv','RS25_DNA_c2.csv','RS25_RNA_c2.csv']):
    # read RNA and DNA chain lists from NAPAIR, process chains in the union of the two
    chain_set = set()
    pdb_set = set()
    for filename in filenames:
        with open(os.path.join('C:/Users/zirbel/Documents/PythonFR3D/data/',filename),'rt') as f:
            lines = f.readlines()

        for line in lines:
            fields = line.split(",")
            pdb = fields[0].upper()
            chain = "%s|1|%s" % (pdb,fields[4])
            pdb_set.add(pdb)
            chain_set.add(chain)

        print('Have %4d distinct PDB ids and %4d distinct chain ids after %s' % (len(pdb_set),len(chain_set),filename))

    all_pdb_ids = sorted(pdb_set)
    chain_list = sorted(chain_set)

    print('Have %4d distinct PDB ids and %4d distinct chain ids' % (len(all_pdb_ids),len(chain_list)))

    return all_pdb_ids, chain_list


def load_dssr_basepairs(pdb_id):
    """
    Load from saved .pickle file, or download from NAKB if not already done
    """

    pickle_filename = os.path.join(dssr_basepair_path,pdb_id+'_pairs.pickle')

    if os.path.exists(pickle_filename):
        if sys.version_info[0] < 3:
            interaction_to_pair = pickle.load(open(pickle_filename,"rb"))
        else:
            interaction_to_pair = pickle.load(open(pickle_filename,"rb"), encoding = 'latin1')

    else:
        interaction_to_pair = {}
        j = 1
        keep_trying = True
        while keep_trying:

            url = 'https://www.nakb.org/x3dssr/%s_%d.json' % (pdb_id,j)

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

    return interaction_to_pair


def load_napair_basepairs_from_csv(filename,resolution):
    # read chains from napair / nabir files and save in .pickle files
    # pdbid,model,family,class,full_subclass,label_alt_id1,label_comp_id1,label_asym_id1,label_entity_id1,label_seq_id1,pdbx_PDB_ins_code1,auth_seq_id1,auth_comp_id1,auth_asym_id1,symmetry_operation1,symmetry_id1,label_alt_id2,label_comp_id2,label_asym_id2,label_entity_id2,label_seq_id2,pdbx_PDB_ins_code2,auth_seq_id2,auth_comp_id2,auth_asym_id2,symmetry_operation2,symmetry_id2,napasco,rmsd,curated_file,knn_metric,coplanarity_angle,coplanarity_shift1,coplanarity_shift2,coplanarity_edge_angle1,coplanarity_edge_angle2,C1_C1_yaw1,C1_C1_pitch1,C1_C1_roll1,C1_C1_yaw2,C1_C1_pitch2,C1_C1_roll2,hb_0_length,hb_0_donor_angle,hb_0_acceptor_angle,hb_0_OOPA1,hb_0_OOPA2,hb_1_length,hb_1_donor_angle,hb_1_acceptor_angle,hb_1_OOPA1,hb_1_OOPA2,hb_2_length,hb_2_donor_angle,hb_2_acceptor_angle,hb_2_OOPA1,hb_2_OOPA2,hb_3_length,hb_3_donor_angle,hb_3_acceptor_angle,hb_3_OOPA1,hb_3_OOPA2
    # 1a1l,1,cWW,cWW_G-C,cWW_G-C_1,"",DG,A,1,2,"",2,DG,B,,1,"",DC,B,2,11,"",61,DC,C,,1,70.59,0.195,4dfk_C_DG_206_1_555_B_DC_111_1_555,1.936,167.90997,-0.16694143,-0.048482712,4.6267805,5.8778987,68.83991,14.175395,-175.51112,69.23006,-9.271735,168.34215,2.9588978,120.70201,119.79926,8.0540495,-2.8031914,2.8780196,119.20008,121.671585,3.3008595,-7.5385513,2.6665695,127.36669,119.720085,-1.037244,-11.929507,NaN,NaN,NaN,NaN,NaN

    # the one below has no percentile and should not be used
    # 5tbw,1,cWS,cWS-U-U,1,2873,U,"", ,,1,2869,U,"", ,,,,,,129.75883,0.12163817,-0.10286882,41.566635,28.494764,30.859066,-31.592983,136.63326,5.0083995,42.776634,147.52533,3.538368,97.25362,137.8698,29.947083,27.315092,3.2720149,128.91217,137.53609,-1.853775,24.402248,,,,,,,,,,

    header = "pdbid,model,family,class,full_subclass,label_alt_id1,label_comp_id1,label_asym_id1,label_entity_id1,label_seq_id1,pdbx_PDB_ins_code1,auth_seq_id1,auth_comp_id1,auth_asym_id1,symmetry_operation1,symmetry_id1,label_alt_id2,label_comp_id2,label_asym_id2,label_entity_id2,label_seq_id2,pdbx_PDB_ins_code2,auth_seq_id2,auth_comp_id2,auth_asym_id2,symmetry_operation2,symmetry_id2,napasco,rmsd,curated_file,knn_metric,coplanarity_angle,coplanarity_shift1,coplanarity_shift2,coplanarity_edge_angle1,coplanarity_edge_angle2,C1_C1_yaw1,C1_C1_pitch1,C1_C1_roll1,C1_C1_yaw2,C1_C1_pitch2,C1_C1_roll2,hb_0_length,hb_0_donor_angle,hb_0_acceptor_angle,hb_0_OOPA1,hb_0_OOPA2,hb_1_length,hb_1_donor_angle,hb_1_acceptor_angle,hb_1_OOPA1,hb_1_OOPA2,hb_2_length,hb_2_donor_angle,hb_2_acceptor_angle,hb_2_OOPA1,hb_2_OOPA2,hb_3_length,hb_3_donor_angle,hb_3_acceptor_angle,hb_3_OOPA1,hb_3_OOPA2"

    interaction_to_triples = defaultdict(list)

    if os.path.exists(filename):
        print('Reading napair file %s' % filename)

        with open(filename,'rt') as f:
            lines = f.readlines()

        header_to_index = map_entry_to_index(lines[0].split(","))

        for line in lines[1:]:
            fields = line.split(",")

            pdb_id = fields[header_to_index['pdbid']].upper()
            model = fields[header_to_index['model']]
            interaction = fields[header_to_index['family']]

            chain1 = fields[header_to_index['auth_asym_id1']]
            base1 = fields[header_to_index['auth_comp_id1']]
            number1 = fields[header_to_index['auth_seq_id1']]
            alt1 = fields[header_to_index['label_alt_id1']]
            ins1 = fields[header_to_index['pdbx_PDB_ins_code1']]
            sym1 = fields[header_to_index['symmetry_operation1']]
            u1 = make_unit_id(pdb_id,model,chain1,base1,number1,"",alt1,ins1,sym1)

            chain2 = fields[header_to_index['auth_asym_id2']]
            base2 = fields[header_to_index['auth_comp_id2']]
            number2 = fields[header_to_index['auth_seq_id2']]
            alt2 = fields[header_to_index['label_alt_id2']]
            ins2 = fields[header_to_index['pdbx_PDB_ins_code2']]
            sym2 = fields[header_to_index['symmetry_operation2']]
            u2 = make_unit_id(pdb_id,model,chain2,base2,number2,"",alt2,ins2,sym2)

            rev = reverse_edges(interaction)

            print('napair',u1,interaction,u2)

            if "." in u1 or "." in u2:
                print('Someone put a . in a unit id! Skipping. %s %s' % (u1,u2))
            else:
                interaction_to_triples[interaction].append((u1,u2,0))       # pretend that range is 0; maybe change it some day
                interaction_to_triples[rev].append((u2,u1,0))

    return interaction_to_triples


def load_pairs_from_pickle_files(all_pdb_ids, method, basepair_path, pdb_id_to_annotators, resolution):
    """
    Load .pickle files of all basepairs from the given method
    """

    pair_to_interaction = defaultdict(str)
    missing = []
    empty = []
    pdb_list = set()

    filename = os.path.join(basepair_path,'all_pairs_%s.pickle' % resolution)
    if os.path.exists(filename):
        pair_to_int, pdb_list, empty, missing = pickle.load(open(filename,'rb'))
        pair_to_interaction = defaultdict(str,pair_to_int)
    elif method == 'napair':
        # load from .csv files according to basepair category
        napair_directory = os.path.join(basepair_path,'RS_filtered_chains_' + resolution.replace('RS',''))
        import glob
        # get file directory and iterate over file names
        csv_files = glob.glob(os.path.join(napair_directory, "*.csv"))

        print('Found %d files in %s' % (len(csv_files),napair_directory))
        # input("Press Enter when ready to read annotations")

        for csv_path_file in csv_files:

            csv_file = os.path.split(csv_path_file)[1]

            print("Reading %s from %s" % (csv_file,csv_path_file))
            # input("Press Enter when ready to read annotations")

            if len(csv_file.split("_")) == 3:
                # example:  cHS_C-A_1.csv

                print('Processing %s for the first time' % csv_path_file)

                interaction_to_triples = load_napair_basepairs_from_csv(csv_path_file,resolution)
                for interaction in interaction_to_triples.keys():
                    # store all basepairs and only basepairs
                    if "c" in interaction or "t" in interaction:
                        for u1,u2,crossing in interaction_to_triples[interaction]:
                            pair_to_interaction[(u1,u2)] = interaction
                        if len(interaction_to_triples[interaction]) > 0:
                            pdb_id_from_u1 = u1.split("|")[0]
                            pdb_list.add(pdb_id_from_u1)

        if len(pdb_list) > 0:
            with open(filename,'wb') as fh:
                pickle.dump((pair_to_interaction, pdb_list, empty, missing), fh, 5)
        else:
            print("Problem with reading NAPAIR annotation files")
            input("Press Enter to continue")
    else:
        for pdb_id in all_pdb_ids:
            # windows tolerates case differences in filenames, that is helpful here

            # pairs from non-filtered data files
            pickle_filename = os.path.join(basepair_path,pdb_id+'_pairs.pickle')

            if not os.path.exists(pickle_filename):
                if method == 'napair':
                        csv_path_file = os.path.join(napair_basepair_path,pdb_id.lower()+".csv")
                        interaction_to_triples = load_napair_basepairs_from_csv(csv_path_file,resolution)
                        with open(pickle_filename,'wb') as fh:
                            pickle.dump(interaction_to_triples, fh, 5)
                elif method == 'DSSR':
                    interaction_to_triples = load_dssr_basepairs(pdb_id)

            if os.path.exists(pickle_filename):
                if sys.version_info[0] < 3:
                    interaction_to_triples = pickle.load(open(pickle_filename,"rb"))
                else:
                    interaction_to_triples = pickle.load(open(pickle_filename,"rb"), encoding = 'latin1')

                    if pdb_id == "4WO2":
                        print(interaction_to_triples)
                        input("Press Enter when ready")

                    if len(interaction_to_triples) == 0:
                        empty.append(pdb_id)
                    else:
                        for interaction in interaction_to_triples.keys():
                            # store all basepairs and only basepairs
                            if "c" in interaction or "t" in interaction:
                                for u1,u2,crossing in interaction_to_triples[interaction]:
                                    pair_to_interaction[(u1,u2)] = interaction
                                if len(interaction_to_triples[interaction]) > 0:
                                    pdb_list.add(pdb_id)
            else:
                missing.append(pdb_id)

        with open(filename,'wb') as fh:
            pickle.dump((pair_to_interaction, pdb_list, empty, missing), fh, 5)

    for pdb_id in pdb_list:
        pdb_id_to_annotators[pdb_id].add(method)

    print('Missing %s annotations from these %d files:\n%s' % (method,len(missing)," ".join(sorted(missing))))
    print('Empty %s annotations from these %d files:\n%s' % (method,len(empty)," ".join(sorted(empty))))

    return pair_to_interaction, pdb_id_to_annotators, empty, missing


def load_curated_pairs():
    # read a list of pairs that represent various clusters of true pairs
    #           0    1 2 3    4 5   6 7 8    9 10
    # cWW_A-U_1,1s72_0_A_1242_1_555_0_U_1122_1_555
    # cSS_G-U_1,4ybb_AA_G_1268_1_555_AA_U_1326_1_555

    curated_pairs_filename = 'curated_files_list.csv'
    curated_pairs_filename = 'curated_list_2026-02-20.csv'

    print("Reading %s" % curated_pairs_filename)

    curated_pairs = set()
    pair_to_interaction_curated = {}
    interactions = set()

    with open(curated_pairs_filename,'rt') as file:
        lines = file.readlines()

    for line in lines:
        if len(line) < 30:
            continue
        f = line.rstrip("\n").split(",")
        interaction = f[0].split("_")[0]
        interactions.add(interaction)
        fields = f[1].split("_")
        pdb = fields[0].upper()
        model = "1"
        chain1 = fields[1]
        base1  = fields[2]
        num1   = fields[3]
        if "." in num1:
            num1,ins1 = num1.split(".")
        else:
            ins1 = ''
        sym1 = fields[4] + "_" + fields[5]
        if sym1 == "1_555":
            sym1 = ""

        chain2 = fields[6]
        base2  = fields[7]
        num2   = fields[8]
        if "." in num2:
            num2,ins2 = num2.split(".")
        else:
            ins2 = ''
        sym2 = fields[9] + "_" + fields[10]
        if sym2 == "1_555":
            sym2 = ""

        u1 = make_unit_id(pdb,model,chain1,base1,num1,'','',ins1,sym1)
        u2 = make_unit_id(pdb,model,chain2,base2,num2,'','',ins2,sym2)
        curated_pairs.add((u1,u2))
        curated_pairs.add((u2,u1))

        pair_to_interaction_curated[(u1,u2)] = interaction
        pair_to_interaction_curated[(u2,u1)] = reverse_edges(interaction)

    print("Found these interactions in NAPAIR curated pairs:")
    print(sorted(interactions))
    print(sorted([reverse_edges(x) for x in interactions]))

    return curated_pairs, pair_to_interaction_curated


def make_unit_id(pdb,model,chain,unit,number,atom,alt,ins,sym):

    if unit == "T":
        unit = "DT"

    alt = alt.replace('"','')
    ins = ins.replace('"','')
    sym = sym.replace('"','').replace('1_555','')

    unit_id = pdb + "|" + model + "|" + chain + "|" + unit + "|" + number + "|" + atom + "|" + alt + "|" + ins + "|" + sym
    while unit_id[-1] == "|":
        unit_id = unit_id[:-1]

    return unit_id


def add_pairs_in_order(pair_to_priority,pair_to_data,priority):
    """
    Assign a priority score to each pair.
    Lower priority numbers will be considered earlier.
    This attempts to improve homogeneity in each group.
    """

    for pair, datapoint in pair_to_data.items():
        # do not override the priority number of a pair
        if pair in pair_to_priority:
            continue

        # delay pairs with bases in an unusual order
        # sometimes b1 or b2 will be None or empty if not recognized
        b1, b2, parent_bc = get_parent_pair(pair)

        if not b1:
            u1, u2 = pair
            s1 = u1.split("|")[3]
            if not s1 in unknown_bases:
                print('Unknown base %s in method with priority %s' % (s1,priority))
                unknown_bases.add(s1)
        if not b2:
            u1, u2 = pair
            s2 = u2.split("|")[3]
            if not s2 in unknown_bases:
                print('Unknown base %s in method with priority %s' % (s2,priority))
                unknown_bases.add(s2)

        (u1,u2) = pair
        if (u1 and "." in u1) or (u2 and "." in u2):
            print('Found a period in a unit id for %s %s' % (u1,u2))

        if b1 and b2:
            bc = b1+b2
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
    datapoint['mingap'] = min(datapoint['gap21'],datapoint['gap12'])

    if 'basepair' in datapoint:
        python_annotation = datapoint['basepair']
    else:
        python_annotation = ''

    new_python_annotation = python_annotation  # in case nothing else is set

    # find worst hydrogen bond
    max_badness = -1.0       # calculated by hydrogen_bonds.py
    # max_distance = -1.0
    # min_angle = 180
    if 'hbond' in datapoint:
        for hbond in datapoint['hbond']:
            # hbond is like
            # (True, True, 2.391828053713497, 127.44856702606579, 1.1275716486967107)
            if hbond["bond_checked"]:    # able to check the hydrogen bond
                if hbond["badness"] > max_badness:
                    max_badness = hbond["badness"]
                # if hbond["distance"] > max_distance:
                #     max_distance = hbond["distance"]
                # if not math.isnan(hbond["angle"]) and hbond["angle"] < min_angle:
                #     min_angle = hbond["angle"]

    # check new cutoffs in every case, call that "FR3D detail"
    keep_reasons = []
    new_cutoff_distance = 9999
    new_python_subcat = -1

    near_discrepancy_cutoff = 2.0     # maximum discrepancy to report as a near pair

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
                cutoff_distance += datapoint['y'] - cutoff["ymax"]
                reasons.append("ymax")

            if "radiusmax" in cutoff:
                radius = math.sqrt(datapoint['x']**2 + datapoint['y']**2)
                if radius > cutoff["radiusmax"]:
                    cutoff_distance += radius - cutoff["radiusmax"]
                    reasons.append("radiusmax")

                if "radiusmin" in cutoff:
                    if radius < cutoff["radiusmin"]:
                        cutoff_distance += cutoff["radiusmin"] - radius
                        reasons.append("radiusmin")

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

                # if abs(datapoint['normal_Z']) < 0.4:
                #     cutoff_distance += near_discrepancy_cutoff
                #     reasons.append("normZ")
            else:
                cutoff_distance += 1
                reasons.append("nmax")

            if "min_distance" in datapoint:
                if abs(datapoint["min_distance"]) < 0.5:
                    # unrealistically close to one another, cannot be a basepair
                    cutoff_distance += near_discrepancy_cutoff
                    reasons.append("mind")

            if 'angle_in_plane' in datapoint:
                angle_penalty = 0.1
                if cutoff["anglemin"] < cutoff["anglemax"]:
                    # usual order where min < max
                    if datapoint['angle_in_plane'] < cutoff["anglemin"]:
                        cutoff_distance += angle_penalty*(cutoff["anglemin"] - datapoint['angle_in_plane'])
                        reasons.append("angle")
                    if datapoint['angle_in_plane'] > cutoff["anglemax"]:
                        cutoff_distance += angle_penalty*(datapoint['angle_in_plane'] - cutoff["anglemax"])
                        reasons.append("angle")
                else:
                    # min might be 260 and max might be -60, looking for angles above 260 or below -60
                    if datapoint['angle_in_plane'] < cutoff["anglemin"] and \
                       datapoint['angle_in_plane'] > cutoff["anglemax"]:
                        cutoff_distance += angle_penalty*min(cutoff["anglemin"]-datapoint["angle_in_plane"],datapoint['angle_in_plane']-cutoff["anglemax"])
                        reasons.append("angle")
            else:
                cutoff_distance += 1
                reasons.append("angle")

            if datapoint['maxgap'] > cutoff["gapmax"]:
                cutoff_distance += 3*(datapoint['maxgap'] - cutoff["gapmax"])  # strong penalty
                reasons.append('gap')

            # extra penalty when both are out of plane with each other
            if datapoint["mingap"] > cutoff['gapmax']:
                cutoff_distance += 3*(datapoint["mingap"] - cutoff['gapmax'])

            if cutoff_distance < new_cutoff_distance:
                # met cutoffs better than with previous subcategory
                keep_reasons = reasons
                if cutoff_distance > near_discrepancy_cutoff:
                    new_python_annotation = ''
                elif cutoff_distance > 0 and not 'n' in bp:
                    new_python_annotation = 'n' + bp
                else:
                    new_python_annotation = bp
                new_python_subcat = bpsc
                new_cutoff_distance = cutoff_distance

    if new_cutoff_distance > 0:
        new_fr3d_detail = ",".join(keep_reasons)
    else:
        new_fr3d_detail = new_python_annotation

    if 'demoted_hbond' in datapoint:
        new_fr3d_detail = 'demoted,' + new_fr3d_detail
        if not 'n' in new_python_annotation:
            new_python_annotation = "n" + new_python_annotation

    # keep data to put in scatterplots and HTML pages
    # only keep pairs if their parameters are somewhat close to the cutoffs for some category

    # store data about this pair for output
    pdata = datapoint
    # pdata['max_distance'] = max_distance
    # pdata['min_angle'] = min_angle
    pdata['max_badness'] = max_badness
    pdata['python_annotation'] = python_annotation          # original annotation from NA_pairwise_interactions
    pdata['new_python_annotation'] = new_python_annotation  # new from checking cutoffs above
    pdata['new_fr3d_detail'] = new_fr3d_detail              # new from checking cutoffs above
    pdata['new_basepair_subcat'] = new_python_subcat
    pdata['new_cutoff_distance'] = new_cutoff_distance
    pdata['keep_reasons'] = keep_reasons

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

    print('Starting')

    symmetric_basepair_list = ['cWW','cWw','cwW','tWW','cHH','tHH','cSS','cSs','csS','tSS','tSs','tsS']
    Leontis_Westhof_basepairs = ['cWW','tWW','cWH','tWH','cWS','tWS','cHH','tHH','cHS','tHS','cSS','tSS','cWB']
    LW12 = ['cWW','tWW','cWH','tWH','cWS','tWS','cHH','tHH','cHS','tHS','cSS','tSS']
    LW12b = ['cWw','tWW','cWH','tWH','cWS','tWS','cHh','tHh','cHS','tHS','cSs','tSs']
    LW12c = ['cWw','tWW','cWH','tWH','cWS','tWS','cHh','tHh','cHS','tHS','cSs','tSS']
    LW18 = ['cWW','tWW','cWH','cHW','tWH','tHW','cWS','cSW','tWS','tSW','cHH','tHH','cHS','cSH','tHS','tSH','cSs','csS','tSs','tsS']

    LW_to_number = {}
    for i,LW in enumerate(Leontis_Westhof_basepairs):
        LW_to_number[LW] = i+1

    base_combination_list = ['G,C','A,A','A,C','A,G','A,U','C,C','C,U','G,G','G,U','U,U']
    base_combination_list = ['A,U','A,A','A,C','A,G','C,C','G,C','C,U','G,G','G,U','U,U']
    base_combination_list = ['C,C','G,C','C,U','G,G','G,U','U,U','A,A','A,C','A,G','A,U']
    base_combination_list = ['C,U','G,G','G,U','U,U','A,A','A,C','A,G','A,U','C,C','G,C']
    base_combination_list = ['A,A','A,C','A,G','A,U','C,C','G,C','C,U','G,G','G,U','U,U']
    base_combination_list = ['G,G','G,U','U,U','A,A','A,C','A,G','A,U','C,C','G,C','C,U']

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

    standard_nucleotides = ['A','C','G','U','DA','DC','DG','DT']

    red   = [1,0,0] # red
    black = [0,0,0] # black
    cyan  = [0,1,1] # cyan
    blue  = [0,0,1] # blue
    orange = [1,165/255,0] # orange
    purple = [73/255,0,146/255] # purple
    purple = [218/255,112/255,214/255] # orchid, shows up better than purple

    write_html_pages = False
    write_html_pages = True

    only_modified = True           # only show basepairs involving modified nucleotides
    only_modified = False

    write_hbond_tables = False     # does not work well
    write_hbond_tables = True      # write a table listing hbond atoms, distances, angles

    show_datapoint_fields = 50           # show this number of datapoint dictionaries
    datapoint_keys = set()

    num_instances_target = 250000    # maximum target number of instances to load and consider

    DNA = True
    DNA = False

    VERSION = ''      # not comparing annotators
    VERSION = 'v6'    # not sure exactly what that was
    VERSION = 'v7'    # only FR3D-annotated basepairs, show hDistAngle, omit C-H.. bonds
    VERSION = 'v8'    # show FR3D or napair annotated basepairs or demoted, show hDist, hDistAngle, include C-H.. bonds
    VERSION = 'v9'    # compare FR3D and napair

    all_agree = '111111'           # indicates how many annotators are being compared
    compare_annotators = False     # write HTML pages for internal evaluation of FR3D annotations
    compare_annotators = True      # write HTML pages to compare FR3D, RNAview, MC-Annotate, DSSR, PDB, napair for basepair group

    show_fr3d_parameters = True
    show_fr3d_parameters = False   # only show 300 candidates, do not show all those details

    # O2_bond_length = 'short'
    # O2_bond_length = 'long'
    O2_bond_length = ''

    # zzz quick way to get to this location

    if VERSION == 'v9':
        resolution_list = ['25','35']
        resolution_list = ['35','25']
    else:
        resolution_list = ['1.5A','2.0A','2.5A','3.0A']

    PDB_skip_set = set()
    unit_id_to_glycosidic = defaultdict(str)

    for resolution in resolution_list:

        pdb_id_to_annotators = defaultdict(set)   # which programs annotate pairs in each structure

        if write_hbond_tables:
            hbond_data_list = []
            dist2_data_list = []   # values to output for each basepair category

        if DNA:
            from DNA_2A_list import PDB_list   # define PDB_list as a list of DNA structures
            PDB_list.remove('5G35')        # has AG pairs that overlap
            PDB_list = sorted(PDB_list)

            base_combination_list = ['DA,DA','DA,DC','DA,DG','DA,DT','DC,DC','DG,DC','DC,DT','DG,DG','DG,DT','DT,DT']
            data_file = []

        elif VERSION in ['v9']:
            print('Loading NAPAIR curated pairs')
            curated_pairs, pair_to_interaction_curated = load_curated_pairs()

            PDB_curated = set([x[0].split("|")[0] for x in curated_pairs])
            curated_chains = set()
            for x in curated_pairs:
                for a in x:
                    fields = a.split("|")
                    curated_chains.add("%s|1|%s" % (fields[0],fields[2]))

            napair_pdbs, napair_chains = load_napair_chains(filenames = ['RS_filtered_chains_' + resolution])

            PDB_list = sorted(napair_pdbs)
            all_pdb_ids = PDB_list

            # write the list of NAPAIR PDB files so that RNAview, MC-Annotate can be run on them
            if resolution == '35':
                with open('reference_list_2026_08.py','wt') as f:
                    f.write('# reference set 35 from August 2026\n')
                    f.write('PDB_lists = [%s]' % PDB_list)

            print('Loading NAPAIR  annotations from %s' % napair_basepair_path)
            pair_to_interaction_napair, pdb_id_to_annotators, empty, missing = load_pairs_from_pickle_files(all_pdb_ids,'napair', napair_basepair_path, pdb_id_to_annotators, resolution)

            print('Found %d NAPAIR chains' % len(napair_chains))

            interaction_pair_to_count = defaultdict(int)
            napair_chain_set = set(napair_chains)
            chains_seen = set()
            for pair, interaction in pair_to_interaction_napair.items():
                b1, b2, parent_bc = get_parent_pair(pair)
                u1, u2 = pair
                # extract PDB and chain, set model to 1 to track chains seen
                c1 = "%s|1|%s" % ((u1.split("|")[0],u1.split("|")[2]))
                c2 = "%s|1|%s" % ((u2.split("|")[0],u2.split("|")[2]))
                if c1 in napair_chain_set or c2 in napair_chain_set:
                    interaction_pair = "%s %s %s" % (interaction,b1,b2)
                    interaction_pair_to_count[interaction_pair] += 1
                    chains_seen.add(c1)
                    chains_seen.add(c2)
                else:
                    print(u1,u2,c1,c2,"seem to be the wrong chains")

            print('NAPAIR interaction counts in resolution %s:' % resolution)
            for interaction_pair, count in sorted(interaction_pair_to_count.items()):
                print(interaction_pair, count)

            print('NAPAIR chains that do not have a basepair in them:')
            print(sorted(napair_chain_set - chains_seen))

            # print('found %5d napair  basepairs' % len(pair_to_interaction_napair.keys()))
            # print(sorted(pair_to_interaction_napair.items())[0:20])

            # print('found %5d curated basepairs' % len(pair_to_interaction_curated.keys()))
            # print(sorted(pair_to_interaction_curated.items())[0:20])

            pair_to_interaction_napair.update(pair_to_interaction_curated)
            representative_chains = napair_chains

            print('Found %5d NAPAIR basepairs' % len(pair_to_interaction_napair.keys()))
            # print(sorted(pair_to_interaction_napair.items())[0:20])

            for pdb_id in all_pdb_ids:
                pdb_id_to_annotators[pdb_id].add('napair')

        else:
            if resolution == '1.5A':
                PDB_list = ['https://rna.bgsu.edu/rna3dhub/nrlist/download/NR/3.349/1.5A/csv']
            elif resolution == '2.0A':
                PDB_list = ['https://rna.bgsu.edu/rna3dhub/nrlist/download/NR/3.349/2.0A/csv']
            elif resolution == '2.5A':
                PDB_list = ['https://rna.bgsu.edu/rna3dhub/nrlist/download/NR/3.349/2.5A/csv','8GLP','8B0X']
            elif resolution == '3.0A':
                PDB_list = ['https://rna.bgsu.edu/rna3dhub/nrlist/download/NR/3.349/3.0A/csv','8GLP','8B0X']

            PDB_IFE_Dict = map_PDB_list_to_PDB_IFE_dict(PDB_list)

            # load all datapoints on pairs of bases, whether annotated as paired or not
            all_pdb_ids = sorted(PDB_IFE_Dict.keys())

            representative_chains = set(['8GLP|1|L5','8GLP|1|L8','8GLP|1|S2','8B0X|1|a','8B0X|1|A'])
            for PDB, chains in PDB_IFE_Dict.items():
                for chain in chains.split("+"):
                    representative_chains.add(chain)

            # PDB_skip_set = set(['1R9F','5NXT','4KTG','7JIL'])

            # if '8B0X' in PDB_list:
            #     PDB_skip_set.add('5J7L')
            #     PDB_skip_set.add('7K00')
            #     PDB_skip_set.add('4YBB')

        Q = readPDBDatafile({"DATAPATHUNITS": os.path.join(fr3d_pickle_path,'units')})  # available PDB structures, resolutions, chains
        data_file = Q["PDB_data_file"]

        print("Resolution %s, working on %d PDB files" % (resolution,len(all_pdb_ids)))

        print('Counting number of pairs of each category')
        category_to_count = defaultdict(int)
        for pair, interaction in pair_to_interaction_napair.items():
            b1, b2, parent_bc = get_parent_pair(pair)
            category = b1 + "_" + b2 + "_" + interaction
            category_to_count[category] += 1


            u1,u2 = pair
            if "." in u1 or "." in u2:
                print('Found . in unit id when counting %s %s' % (u1,u2))
                print(crash)


        pair_to_interaction_dssr       = defaultdict(str)
        pair_to_interaction_rnaview    = defaultdict(str)
        pair_to_interaction_mcannotate = defaultdict(str)
        pair_to_interaction_pdb        = defaultdict(str)

        skip_pair_set = set()
        model_notification_counter = 0

        print('Loading RNAview  annotations from %s' % rnaview_basepair_path)
        pair_to_interaction_rnaview, pdb_id_to_annotators, empty, missing = load_pairs_from_pickle_files(all_pdb_ids,'RNAview', rnaview_basepair_path, pdb_id_to_annotators, resolution)

        print('Loading MC-Annotate annotations from %s' % mcannotate_basepair_path)
        pair_to_interaction_mcannotate, pdb_id_to_annotators, empty, missing = load_pairs_from_pickle_files(all_pdb_ids,'MC-Annotate', mcannotate_basepair_path, pdb_id_to_annotators, resolution)

        print('Loading DSSR     annotations from %s' % dssr_basepair_path)
        pair_to_interaction_dssr, pdb_id_to_annotators, empty, missing = load_pairs_from_pickle_files(all_pdb_ids,'DSSR', dssr_basepair_path, pdb_id_to_annotators, resolution)

        print('Loading PDB      annotations from %s' % pdb_basepair_path)
        pair_to_interaction_pdb, pdb_id_to_annotators, empty, missing = load_pairs_from_pickle_files(all_pdb_ids,'PDB', pdb_basepair_path, pdb_id_to_annotators, resolution)

        if not resolution in ['25','35']:
            print('Loading NAPAIR   annotations from %s' % napair_basepair_path)
            pair_to_interaction_napair, pdb_id_to_annotators, empty, missing = load_pairs_from_pickle_files(all_pdb_ids,'NAPAIR', napair_basepair_path, pdb_id_to_annotators, resolution)

        input("Press Enter to continue")

        count_annotations = {}
        count_annotations_by_group = {}

        # loop over specified base combinations
        for bc_num, base_combination in enumerate(base_combination_list):

            skip_symmetry_pairs = set()
            skip_complicated_unit_id = set()
            skip_over_target = set()
            skip_not_all_annotate = set()

            # show_memory_usage()

            # h-bond data takes up so much space, sometimes it's necessary to work one base combination at a time
            # put '3.0A' or other resolution in the list to load one base combination at a time; memory
            if compare_annotators:
                # will need to keep many additional pairs for comparison
                resolution_memory_challenge_list = ['3.0A']  # RS25 is fine
            else:
                resolution_memory_challenge_list = ['3.0A']

            if bc_num == 0 or resolution in resolution_memory_challenge_list:
                # load each base combination separately to reduce memory usage, though it may take more time
                # load fr3d-python files with datapoint information from NA_pairwise_interactions
                if resolution in resolution_memory_challenge_list:
                    print("Loading FR3D annotations for %s from %s" % (base_combination,outputNAPairwiseInteractions))
                else:
                    print("Loading FR3D annotations for all base combinations from %s" % (outputNAPairwiseInteractions))
                pair_to_datapoint = defaultdict(dict)
                c = 0
                st = time.time()
                for pdb_id in all_pdb_ids:
                    c += 1
                    if c % 50 == 0:
                        if resolution in resolution_memory_challenge_list:
                            print('Reading file %4s # %4d of %d for %s at %s, have %7d pairs, est total time %5.0f seconds, est remaining %8.0f seconds' % (pdb_id,c,len(all_pdb_ids),base_combination,resolution,len(pair_to_datapoint),(time.time()-st)*(len(all_pdb_ids)/c),(time.time()-st)*((len(all_pdb_ids)-c)/c)))
                        else:
                            print('Reading file %4s # %4d of %d for all base combinations at %s, have %7d pairs, est total time %5.0f seconds, est remaining %8.0f seconds' % (pdb_id,c,len(all_pdb_ids),resolution,len(pair_to_datapoint),(time.time()-st)*(len(all_pdb_ids)/c),(time.time()-st)*((len(all_pdb_ids)-c)/c)))

                    pair_to_datapoint_file = outputNAPairwiseInteractions + "%s_datapoint.pickle" % pdb_id

                    if os.path.exists(pair_to_datapoint_file):
                        if sys.version_info[0] < 3:
                            new_dict = pickle.load(open(pair_to_datapoint_file,'rb'))
                        else:
                            new_dict = pickle.load(open(pair_to_datapoint_file,'rb'),encoding = 'latin1')

                        if len(new_dict.keys()) > 0:
                            # make sure at least one pair was annotated before loading all the data
                            found_basepair = False
                            for pair,datapoint in new_dict.items():
                                if 'basepair' in datapoint and datapoint['basepair']:
                                    found_basepair = True
                                    pdb_id_to_annotators[pdb_id].add('FR3D')
                                    break

                            if found_basepair:
                                for pair, datapoint in new_dict.items():
                                    # reduce memory usage by only storing pairs we will process
                                    u1,u2 = pair

                                    # this omits some relevant pairs that are outside of the classification limits
                                    # if not 'basepair' in datapoint:
                                    #     continue

                                    # skip doubly-annotated napair pairs, for one thing
                                    if (u1,u2) in skip_pair_set:
                                        continue

                                    # only keep pairs from representative chains, when that applies
                                    fields1 = u1.split("|")
                                    chain1 = "|".join(fields1[0:3])
                                    fields2 = u2.split("|")
                                    chain2 = "|".join(fields2[0:3])
                                    if VERSION in ['v9'] and representative_chains:
                                        if not chain1 in representative_chains and not chain2 in representative_chains:
                                            # need at least one nucleotide in a representative chain
                                            continue
                                    elif not chain1 in representative_chains or not chain2 in representative_chains:
                                        # need both nucleotides in a representative chain
                                        continue

                                    # if both are from a symmetry operator, skip, because we probably already have it
                                    if len(fields1) == 9 and len(fields2) == 9:
                                        #print('Skipping %s - %s because both have a symmetry operator' % (u1,u2))
                                        skip_symmetry_pairs.add((u1,u2))
                                        continue

                                    # if both have an alternate id, like 7MDL|1|E|U|6	7MDL|1|E|U|67, remove both
                                    if len(fields1) == 7 and len(fields2) == 7:
                                        if fields1[6] == fields2[6]:
                                            interaction = datapoint.get('basepair','')
                                            print('Removing alt ids from %s %s %s' % (u1,interaction,u2))
                                            u1 = "|".join(fields1[0:5])
                                            u2 = "|".join(fields2[0:5])
                                            pair = (u1,u2)
                                            if pair in pair_to_datapoint:
                                                # do not overwrite what was already there
                                                continue
                                            # print('New pair is %s %s %s' % (u1,interaction,u2))

                                    # for some resolutions, focus on one base combination at a time
                                    if not DNA and resolution in resolution_memory_challenge_list:
                                        b1, b2, parent_bc = get_parent_pair(pair)

                                        # keep only pairs with the desired base combination; only keep one copy of each interaction
                                        if not base_combination == b1+","+b2:
                                            continue

                                        # cap the number of instances of very common basepairs at num_instances_target
                                        if 'basepair' in datapoint:
                                            interaction = datapoint['basepair']
                                            category = b1 + "_" + b2 + "_" + interaction
                                            num_instances = category_to_count[category]
                                            if num_instances > num_instances_target and random.random() > num_instances_target/num_instances:
                                                # print('FR3D category %9s, %6d instances, omitting %s %s' % (category,num_instances,u1,u2))
                                                skip_pair_set.add(pair)
                                                skip_pair_set.add((u2,u1))
                                                skip_over_target.add((u1,u2))
                                                continue

                                    # remove unnecessary fields to conserve memory
                                    try:
                                        del datapoint["center_center_distance"]
                                        del datapoint["url"]
                                    except:
                                        pass

                                    if show_datapoint_fields > 0:
                                        datapoint_keys = datapoint_keys | set(datapoint.keys())
                                        show_datapoint_fields -= 1
                                        # if show_datapoint_fields == 0:
                                        #     print('datapoint keys seen so far:')
                                        #     print(sorted(datapoint_keys))

                                    pair_to_datapoint[pair] = datapoint
                    else:
                        print('Cannot open FR3D annotation file %s' % pair_to_datapoint_file)

                if bc_num == 0 and compare_annotators:
                    # show files not annotated by all programs
                    all_annotate_counter = 0
                    for pdb_id, annotators in sorted(pdb_id_to_annotators.items()):
                        if len(annotators) < len(all_agree):
                            if len(annotators) > 1 or not list(annotators)[0] == 'napair':
                                print('%s only annotated by %s' % (pdb_id,sorted(annotators)))
                        else:
                            all_annotate_counter += 1
                    print('%d of %d files are annotated by all %d annotators' % (all_annotate_counter,len(all_pdb_ids),len(all_agree)))

                # keep track of pairs, prioritizing those with annotations like cWw over those like cWH/cHW
                pair_to_priority = add_pairs_in_order({},pair_to_datapoint,1)

                if compare_annotators:
                    # only consider additional pairs if we explicitly want to compare to other annotators
                    pair_to_priority = add_pairs_in_order(pair_to_priority,pair_to_interaction_rnaview,2)
                    pair_to_priority = add_pairs_in_order(pair_to_priority,pair_to_interaction_mcannotate,3)
                    pair_to_priority = add_pairs_in_order(pair_to_priority,pair_to_interaction_dssr,4)
                    pair_to_priority = add_pairs_in_order(pair_to_priority,pair_to_interaction_pdb,5)
                    pair_to_priority = add_pairs_in_order(pair_to_priority,pair_to_interaction_napair,6)

                # focus_pair = ("4V9F|1|0|A|2812","4V9F|1|0|A|2814")
                # print(focus_pair,' 2143 priority ',pair_to_priority.get(focus_pair,'nothing'))

                # remove pairs with symmetry operators, alternate ids, insertion codes when comparing annotators
                if compare_annotators:
                    for u1,u2 in pair_to_priority.keys():
                        fields1 = u1.split("|")
                        fields2 = u2.split("|")

                        if True  and len(fields1) > 5:
                            skip_complicated_unit_id.add(u1)                      # blacklist the full unit id
                            # skip_complicated_unit_id.add("|".join(fields1[0:5]))  # blacklist the plain unit id as well
                            # skip_complicated_unit_id.add("|".join(fields2[0:5]))  # blacklist the paired plain unit id as well
                            #print('Omitting',u1,"|".join(fields1[0:5]),u2,"|".join(fields2[0:5]))
                        if True  and len(fields2) > 5:
                            skip_complicated_unit_id.add(u2)                      # blacklist the full unit id
                            # skip_complicated_unit_id.add("|".join(fields2[0:5]))  # blacklist the plain unit id as well
                            # skip_complicated_unit_id.add("|".join(fields1[0:5]))  # blacklist the paired plain unit id as well
                            #print('Omitting',u1,"|".join(fields1[0:5]),u2,"|".join(fields2[0:5]))

                        # some x-ray structures like 3CMY have a model 2 that RNAview does not mark as model 2,
                        # so it gets reported as being in model 1.  Remove every such pair, even from model 1.
                        if True  and (not fields1[1] == '1' or not fields2[1] == '1'):
                            skip_complicated_unit_id.add(u1)                      # blacklist the full unit id
                            skip_complicated_unit_id.add(u2)                      # blacklist the full unit id
                            fields1[1] = '1'
                            skip_complicated_unit_id.add("|".join(fields1[0:5]))  # blacklist any model 1 version of the unit id
                            fields2[1] = '1'
                            skip_complicated_unit_id.add("|".join(fields2[0:5]))  # blacklist any model 1 version of the unit id
                            if model_notification_counter < 20:
                                print('Removing %24s %24s because of model numbers other than 1' % (u1,u2))
                                model_notification_counter += 1

                print('Found %d pair_to_priority pairs to work with; culling pairs next' % len(pair_to_priority))

                # cull out pairs that we don't want to consider
                # some of these criteria should be used above to save memory
                pair_to_priority_final = {}
                for pair,priority in pair_to_priority.items():

                    if True  and pair in skip_pair_set:
                        continue

                    u1,u2 = pair
                    if True  and u1 in skip_complicated_unit_id:
                        continue

                    if True  and u2 in skip_complicated_unit_id:
                        continue

                    fields1 = u1.split("|")
                    chain1 = "|".join(fields1[0:3])
                    fields2 = u2.split("|")
                    chain2 = "|".join(fields2[0:3])

                    # only keep pairs from representative chains
                    if VERSION in ['v9'] and representative_chains:
                        if not chain1 in representative_chains and not chain2 in representative_chains:
                            continue
                    elif not chain1 in representative_chains or not chain2 in representative_chains:
                        continue

                    # DSSR has hyphens because of symmetry operators being applied
                    if True  and '-' in chain1:
                        continue
                    if True  and '-' in chain2:
                        continue

                    if only_modified and fields1[3] in standard_nucleotides and fields2[3] in standard_nucleotides:
                        continue

                    # if both are from a symmetry operator, skip, because we probably already have it
                    if True  and len(fields1) == 9 and len(fields2) == 9:
                        #print('Skipping %s - %s because both have a symmetry operator' % (u1,u2))
                        skip_symmetry_pairs.add((u1,u2))
                        continue

                    pdb_id = fields1[0]

                    # only include pairs where all annotators are able to make annotations
                    if True  and compare_annotators and not len(pdb_id_to_annotators[pdb_id]) == len(all_agree):
                        skip_not_all_annotate.add((u1,u2))
                        continue

                    # some PDB ids are a problem, just skip them
                    if pdb_id in PDB_skip_set:
                        continue

                    pair_to_priority_final[(u1,u2)] = priority

                print('Found %d pairs_in_order pairs to work with' % len(pair_to_priority_final.keys()))

                # set the order in which to process pairs; doing FR3D first gets more consistent unit orders when symmetric
                pairs_in_order = sorted(pair_to_priority_final.keys(), key=lambda x:pair_to_priority_final[x])

                print('Sorted them')

                # try to reduce memory usage
                pair_to_priority = {}
                pair_to_priority_final = {}

                # store pairs by their base combination for faster retrieval
                # store DNA and modified nucleotides according to their parent nucleotide
                bc_to_pairs_in_order = {}
                for pair in pairs_in_order:
                    b1, b2, parent_bc = get_parent_pair(pair)

                    if b1 and b2:
                        bc = b1+","+b2
                        if not bc in bc_to_pairs_in_order:
                            bc_to_pairs_in_order[bc] = []

                        bc_to_pairs_in_order[bc].append(pair)
                    else:
                        print('base %s became %s, %s became %s' % (pair[0].split("|")[3],b1,pair[1].split("|")[3],b2))

            # try to reduce memory usage
            pairs_in_order = {}
            interactions_processed = set([])

            # loop over conceivable and distinct interactions for that base combination
            for interaction in base_combination_to_interaction[base_combination]:

                pair_counted = set()  # only display / count each pair once on each HTML page
                pair_shown = set()    # only display each pair once on each HTML page

                python_true_count = 0
                python_near_count = 0
                python_demoted_count = 0

                all_dist2_list = []    # dist2 values for each true pair

                angle_out_of_order = False   # set to True if there is a cutoff like 250 to -60 degrees

                interaction_upper = interaction[0] + interaction[1:3].upper()
                interaction_lower = interaction.lower()

                # standard order for base combination and filename
                # use this for counting and for filenames
                bc_filename = base_combination
                inter_filename = interaction_upper
                reverse_bases = False

                # base combination and basepair to use in filenames
                # reverse order of base combination with interactions like tsS because
                # file names for tsS and tSs will be the same on windows
                # also use standard order for the six combinations listed below
                if interaction[1].islower() or interaction_upper in ['cHW','tHW','cSW','tSW','cSH','tSH','cBW']:
                    bases = bc_filename.split(",")
                    bc_filename = bases[1] + "," + bases[0]
                    inter_filename = interaction_upper[0] + interaction_upper[2] + interaction_upper[1]
                    reverse_bases = True

                print('Resolution %s, working on %s %s' % (resolution,base_combination,interaction))
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
                colors  = []    # store the color to use
                sizes = []      # store the size of dot to use

                # accumulate information about hydrogen bonds; multiple ones per basepair
                hdvalues = []   # hydrogen bond distances
                havalues = []   # hydrogen bond angles
                hbvalues = []   # hydrogen bond badness measures

                # store data to produce an HTML table
                pair_data = []  # list of data dictionaries to print in a table

                # identify all atom sets in this family and base combination
                interaction_to_atom_sets = {}
                atom_set_to_d_atoms = {}
                atom_set_to_a_atoms = {}
                bundle_message_count = 0

                # loop over pairs for which we have data, finding those with the interaction and base combination
                for pair in bc_to_pairs_in_order.get(base_combination,[]):

                    if not pair in pair_to_datapoint:
                        u1, u2 = pair
                        if (u2,u1) in pair_to_datapoint:
                            # pick this up next time it is listed
                            continue
                        elif (u2,u1) in pair_to_datapoint and not pair_to_datapoint[(u2,u1)].get('basepair','N/A') == 'N/A':
                            print('Encountered %s %s but other order would be better, with %s' % (u1,u2,pair_to_datapoint[(u2,u1)]['basepair']))
                            if (u2,u1) in bc_to_pairs_in_order:
                                print('Waiting for the next version of it')
                                input("Press Enter")
                                continue
                            else:
                                pair_to_datapoint[pair] = {}
                                pair_to_datapoint[pair]['nt1_seq'] = pair[0].split("|")[3]
                                pair_to_datapoint[pair]['nt2_seq'] = pair[1].split("|")[3]
                                pair_to_datapoint[pair]['basepair'] = reverse_edges(pair_to_datapoint[(u2,u1)]['basepair'])
                                pair_to_datapoint[pair]['basepair_subcategory'] = pair_to_datapoint[(u2,u1)]['basepair']
                                input("Press Enter")

                        elif VERSION in ['v9'] and pair in pair_to_interaction_napair:
                            # print('Pair %s-%s does not have datapoint information https://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s' % (pair[0],pair[1],pair[0],pair[1]))
                            pair_to_datapoint[pair] = {}
                            pair_to_datapoint[pair]['nt1_seq'] = pair[0].split("|")[3]
                            pair_to_datapoint[pair]['nt2_seq'] = pair[1].split("|")[3]
                            pair_to_datapoint[pair]['basepair'] = 'unk'
                            pair_to_datapoint[pair]['basepair_subcategory'] = 9
                        else:
                            continue

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
                    direct_near = False
                    if 'basepair' in datapoint and interaction_lower in datapoint['basepair'].lower():
                        python_fr3d = True
                        if "n" in datapoint['basepair']:
                            python_near = True
                            python_true = False
                            if 'cut_dist' in datapoint and datapoint['cut_dist'] == 0:
                                direct_near = True
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

                    dssr_annotation = pair_to_interaction_dssr[pair]
                    if len(dssr_annotation) > 0 and interaction_lower in dssr_annotation.lower():
                        dssr = True
                    else:
                        dssr = False

                    napair_annotation = pair_to_interaction_napair[pair]
                    # match napair annotation exactly to keep cSs and csS distinct
                    if len(napair_annotation) > 0 and interaction_lower in napair_annotation.lower():
                        napair = True
                    else:
                        napair = False

                    if pair in pair_to_interaction_rnaview:
                        rnaview_annotation = pair_to_interaction_rnaview[pair].replace("n","!")
                    else:
                        # try replacing DA with A, DC with C, DG with G
                        fields = pair[0].split("|")
                        fields[3] = fields[3].replace("DA","A").replace("DC","C").replace("DG","G")
                        p0 = "|".join(fields)
                        fields = pair[1].split("|")
                        fields[3] = fields[3].replace("DA","A").replace("DC","C").replace("DG","G")
                        p1 = "|".join(fields)
                        if (p0,p1) in pair_to_interaction_rnaview:
                            rnaview_annotation = pair_to_interaction_rnaview[(p0,p1)].replace("n","!")
                        else:
                            rnaview_annotation = ''

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

                    if pair in pair_to_interaction_mcannotate:
                        mcannotate_annotation = pair_to_interaction_mcannotate[pair]
                    else:
                        # try replacing DA with A, DC with C, DG with G, DT with T, see if it is annotated
                        fields = pair[0].split("|")
                        fields[3] = fields[3].replace("DA","A").replace("DC","C").replace("DG","G").replace("DT","T")
                        p0 = "|".join(fields)
                        fields = pair[1].split("|")
                        fields[3] = fields[3].replace("DA","A").replace("DC","C").replace("DG","G").replace("DT","T")
                        p1 = "|".join(fields)
                        if (p0,p1) in pair_to_interaction_mcannotate:
                            mcannotate_annotation = pair_to_interaction_mcannotate[(p0,p1)]
                        else:
                            mcannotate_annotation = ''

                    if len(mcannotate_annotation) > 0 and interaction_lower in mcannotate_annotation.lower():
                        mcannotate = True
                    else:
                        mcannotate = False

                    # skip cases where RNAview may need to annotate across chains in different bundles
                    chain1 = pair[0].split("|")[2]
                    chain2 = pair[1].split("|")[2]
                    if True  and compare_annotators and not rnaview and (python_true or python_near or mcannotate or dssr or napair) and not chain1 == chain2 and not DNA:
                        if bundle_message_count < 20:
                            print('Skipping %-24s-%-24s %s because the chains may be from different bundles' % (pair[0],pair[1],interaction))
                            bundle_message_count += 1
                        continue

                    # PDB annotations do not distinguish between edges, so tHS and tSH get the same PDB annotation
                    # only record a PDB annotation when one of the other programs has a matching annotation
                    # and the PDB annotation is consistent with that
                    pdb_annotation = pair_to_interaction_pdb[pair]
                    pdb = False
                    if len(pdb_annotation) > 0 and interaction_lower in pdb_annotation.lower():
                        if (rnaview or mcannotate or python_fr3d or dssr or napair) or inter_filename in ['cWW','tWW','cHH','tHH','cSS','tSS']:
                            pdb = True

                    # special versions for comparing annotators
                    if VERSION in ['v7']:
                        # in v7 we only look at FR3D annotated pairs
                        if not python_true:
                            continue
                    elif VERSION in ['v8']:
                        # in v8 we look at true and near due to demotion due to h-bonds plus napair
                        if python_true:
                            pass
                        elif python_near and 'demoted_hbond' in datapoint:
                            pass
                        elif napair:
                            pass
                        else:
                            continue

                        if not napair and '||' in pair[0] or '||' in pair[1]:
                            # napair does not consistently get all symmetry operators
                            # napair does not have ||A and ||B alternate ids
                            continue

                        if not napair and (not pair[0].split("|")[3] in ['A','C','G','U','DA','DC','DG','DT'] or not pair[1].split("|")[3] in ['A','C','G','U','DA','DC','DG','DT']):
                            # napair does not consistently get all modified bases
                            continue

                    elif VERSION in ['v9']:
                        # look at just fr3d-python true and napair true when you are comparing just them
                        # if python_true:
                        #     pass
                        # elif napair:
                        #     pass
                        # else:
                        #     continue

                        # if any of the six have an an annotation, show that pair
                        if python_true or napair or pdb or mcannotate or rnaview_true or dssr:
                            pass
                        else:
                            continue

                        if True  and not napair and (not pair[0].split("|")[3] in ['A','C','G','U','DA','DC','DG','DT'] or not pair[1].split("|")[3] in ['A','C','G','U','DA','DC','DG','DT']):
                            # skip modified bases because programs cover them inconsistently
                            continue

                    # correctly identify instances as cSs/csS or tSs/tsS; delay csS and tsS to get cSs and tSs
                    if interaction in ['cSs','tSs'] and not nt1_seq == nt2_seq:
                        # you can only delay on cSs or tSs, not on csS or tsS
                        if python_fr3d:
                            if not interaction in datapoint['basepair']:
                                # python_fr3d is csS or tsS, does not match cSs or tSs, so skip until later
                                continue
                        elif compare_annotators and (rnaview_true or mcannotate or dssr or pdb or napair):
                            # print('Checking %s-%s for %s' % (pair[0],pair[1],interaction))
                            datapoint, pdata, angle_order = evaluate_pair_from_datapoint(datapoint,interaction,nt_nt_cutoffs[base_combination])
                            #print("evaluated datapoint:")
                            #print_dictionary(datapoint)
                            if not 'new_python_interaction' in pdata:
                                # not sure how this can happen, because the code always sets this field
                                print('Missing new_python_interaction for %s-%s' % (pair[0],pair[1]))
                                continue
                            if not interaction in pdata['new_python_interaction']:
                                print('Delaying %s-%s with %s in category %s' % (pair[0],pair[1],pdata['new_python_interaction'],interaction))
                                continue

                    if interaction in ['csS','tsS'] and not nt1_seq == nt2_seq and not compare_annotators:
                        if python_fr3d:
                            if not interaction in datapoint['basepair']:
                                # python_fr3d is cSs, don't show it here, wait until later
                                print('Skipping %s %s with %s for %s' % (pair[0],pair[1],datapoint['basepair'],interaction))
                                continue

                    # tally how well the annotators agree on annotating this pair
                    annotator = ""
                    annotator_count = 0
                    annotator_available = 0

                    if python_true or rnaview_true or mcannotate or dssr or pdb or napair:
                        # create a string like 101101 to indicate which of the six methods annotate
                        # use x to indicate a method that did not annotate this PDB file
                        booleans = (python_true,rnaview_true,mcannotate,dssr,pdb,napair)
                        methods = ['FR3D','RNAview','MC-Annotate','DSSR','PDB','napair']

                        for i,truth in enumerate(booleans):
                            if truth:
                                annotator += "1"
                                annotator_count += 1
                                annotator_available += 1
                            elif not methods[i] in pdb_id_to_annotators[pdb_id]:
                                annotator += 'x'
                            else:
                                annotator += "0"
                                annotator_available += 1

                        booleans = (python_true,python_near,rnaview_true,rnaview_near,mcannotate,dssr,pdb,napair)

                        if not inter_filename in count_annotations:
                            count_annotations[inter_filename] = {}
                            count_annotations_by_group[inter_filename] = {}
                        if not bc_filename in count_annotations[inter_filename]:
                            # one for each program, in order
                            count_annotations[inter_filename][bc_filename] = [0,0,0,0,0,0,0,0]
                            count_annotations_by_group[inter_filename][bc_filename] = defaultdict(int)

                        if not "x" in annotator:
                            # omit from the counting any pairs from structures that an annotator count not process
                            # for whatever reason

                            # track how many annotations each method makes
                            for i,truth in enumerate(booleans):
                                if truth:
                                    count_annotations[inter_filename][bc_filename][i] += 1

                            # track how many times we have each annotator string like 101101.
                            count_annotations_by_group[inter_filename][bc_filename][annotator] += 1

                        pair_counted.add(pair)
                        pair_counted.add(reverse(pair))

                    if reverse(pair) in pair_to_datapoint:
                        r_datapoint = pair_to_datapoint[reverse(pair)]
                    else:
                        r_datapoint = {}

                    # do we have all of the parameters that are checked for cutoffs?
                    have_full_data = check_full_data(datapoint)

                    if (python_true) or (annotator_count > 0 and compare_annotators) or \
                        (not compare_annotators and python_near) or (compare_annotators and VERSION == 'v9' and napair):

                        # evaluate the quality of the match to the current pair, for scatterplots and all
                        # print('Evaluating pair %s - %s for interaction %s' % (pair[0],pair[1],interaction))
                        datapoint, pdata, angle_order = evaluate_pair_from_datapoint(datapoint,interaction,nt_nt_cutoffs[base_combination])

                        if angle_order:
                            angle_out_of_order = True

                        if pair[0] == "6E8U|1|B|G|28" or pair[1] == "6E8U|1|B|G|28":
                            print(pair[0],pair[1],interaction,pdata['new_cutoff_distance'],python_true,python_near)
                            # input("Press Enter")

                        # use the reversed pair when bases are the same and basepair is symmetric?
                        if (python_true or python_near) and nt1_seq == nt2_seq and interaction in symmetric_basepair_list and pdata['new_cutoff_distance'] > 0:
                            #print('Evaluating reversed pair %s - %s' % reverse(pair))
                            r_datapoint, r_pdata, angle_order = evaluate_pair_from_datapoint(r_datapoint,interaction,nt_nt_cutoffs[base_combination])
                            # if new match is better and still matches interaction, reverse the pair, use r_datapoint, etc.

                            if pair[0] == "6E8U|1|B|G|28" or pair[1] == "6E8U|1|B|G|28":
                                print(pair[0],pair[1],interaction,pdata['new_cutoff_distance'],r_pdata['new_cutoff_distance'])
                                # input("Press Enter again")


                            if r_pdata['new_cutoff_distance'] < pdata['new_cutoff_distance']:
                                if r_datapoint["python_annotation"].lower() == interaction_lower:
                                    print("Pair %s - %s from %s is better reversed and makes %s" % (pair[0],pair[1],interaction,python_annotation))
                                    pair = reverse(pair)
                                    datapoint = r_datapoint
                                    pdata = r_pdata
                                    # pdata['python_annotation'] = reverse_edges(python_annotation)

                        # elif nt1_seq == nt2_seq and interaction in symmetric_basepair_list:
                        #     if pair in pair_to_datapoint and not reverse(pair) in pair_to_datapoint:
                        #         # keep the original pair order
                        #         pass
                        #     elif False:
                        #         print("Pair %s - %s from %s is better reversed and makes %s" % (pair[0],pair[1],interaction,python_annotation))
                        #         pair = reverse(pair)
                        #         datapoint = r_datapoint
                        #         pdata = r_pdata
                        #         input("Press Enter")

                        # new_cutoff_distance is calculated locally in this program
                        new_cutoff_distance = pdata['new_cutoff_distance']

                        # cut_dist is calculated by NA_pairwise_interactions
                        if not 'cut_dist' in datapoint:
                            datapoint['cut_dist'] = new_cutoff_distance

                        # if interaction.lower() == 'chh':
                        #     print("Testing 2600",interaction,pair)

                        # save the pair if it is good enough to list in the table
                        if python_true \
                            or (VERSION in ['v9'] and napair) \
                            or (not VERSION in ['v9'] and python_near and not compare_annotators and datapoint['cut_dist'] < near_discrepancy_cutoff) \
                            or (not VERSION in ['v9'] and 'sugar_ribose' in datapoint and datapoint['sugar_ribose'] == 'cSR' and datapoint['cut_dist'] < near_discrepancy_cutoff) \
                            or compare_annotators \
                            or (compare_annotators and not interaction in nt_nt_cutoffs[base_combination] and (rnaview or mcannotate or pdb or dssr or napair)) \
                            or (compare_annotators and rnaview and datapoint['cut_dist'] < near_discrepancy_cutoff) \
                            or (compare_annotators and mcannotate and datapoint['cut_dist'] < near_discrepancy_cutoff) \
                            or (compare_annotators and pdb and datapoint['cut_dist'] < near_discrepancy_cutoff) \
                            or (compare_annotators and dssr and datapoint['cut_dist'] < near_discrepancy_cutoff) \
                            or (compare_annotators and napair and datapoint['cut_dist'] < near_discrepancy_cutoff):

                            if python_true:
                                python_true_count += 1
                            if python_near:
                                python_near_count += 1

                            # store for h-bond routine
                            pdata['python_true'] = python_true
                            pdata['napair_true'] = napair

                            pdata['unit_id_1'] = pair[0]
                            pdata['unit_id_2'] = pair[1]
                            pdata['dssr_annotation'] = dssr_annotation
                            pdata['rnaview_annotation'] = rnaview_annotation
                            pdata['mcannotate_annotation'] = mcannotate_annotation
                            pdata['pdb_annotation'] = pdb_annotation
                            # pdata['contacts_annotation'] = contacts_annotation

                            if pair in curated_pairs:
                                pdata['napair_annotation'] = napair_annotation + " cur."
                            else:
                                pdata['napair_annotation'] = napair_annotation

                            pdata['annotator_count'] = annotator_count  # how many annotate as such
                            pdata['annotator'] = annotator

                            # store information about hydrogen bonds for this base combination and interaction
                            # avoid carryover from previous times this datapoint was encountered
                            pdata["atom_sets"] = []

                            # collect atom sets for h-bonds that should be shown in this family and base combination
                            if "LW_to_atom_sets" in datapoint:
                                this_atom_set = []
                                if python_annotation in datapoint["LW_to_atom_sets"]:
                                    # special cases like cWWa
                                    # pdata["atom_sets"] = sorted(datapoint["LW_to_atom_sets"][python_annotation])
                                    this_atom_set = sorted(datapoint["LW_to_atom_sets"][python_annotation])
                                    this_interaction = python_annotation
                                elif interaction in datapoint["LW_to_atom_sets"]:
                                    # pdata["atom_sets"] = sorted(datapoint["LW_to_atom_sets"][interaction])
                                    this_atom_set = sorted(datapoint["LW_to_atom_sets"][interaction])
                                    this_interaction = interaction
                                    #pdata["atom_set_to_results"] = datapoint["atom_set_to_results"]  # it's already there in most cases

                                for atom_set in this_atom_set:
                                    if ("C" in atom_set[0] or "C" in atom_set[2]) and VERSION in ['v7']:
                                        # weak hydrogen bond, ignore in some versions of this output
                                        continue
                                    else:
                                        if not this_interaction in interaction_to_atom_sets:
                                            interaction_to_atom_sets[this_interaction] = []
                                        if not atom_set in interaction_to_atom_sets[this_interaction]:
                                            interaction_to_atom_sets[this_interaction].append(atom_set)

                                            # also save what the header should list for this atom set
                                            result = pdata["atom_set_to_results"][atom_set]

                                            atom_set_to_d_atoms[str(atom_set)] = result["donor_acceptor_atoms"]
                                            atom_set_to_a_atoms[str(atom_set)] = result["heavy_donor_acceptor_atoms"]

                            # if desired, focus on pairs with long or short O2' bond lengths
                            if O2_bond_length:
                                if "atom_set_to_results" in pdata:
                                    atom_set = ('N6', 'H61', "O2'", '12')
                                    if atom_set in pdata["atom_set_to_results"]:
                                        if pdata["atom_set_to_results"][atom_set]["donor_acceptor_distance"] <= 3.5 and O2_bond_length == 'long':
                                            continue
                                        if pdata["atom_set_to_results"][atom_set]["donor_acceptor_distance"] > 3.5 and O2_bond_length == 'short':
                                            continue

                            # set rough position for display; pairs will be sorted from lowest to highest display_position
                            # when there are too many to display
                            if pdata['python_true'] == pdata['napair_true']:
                                # FR3D and NAPAIR agree, whether both true or both false
                                pdata['display_position'] = 2
                            else:
                                # FR3D and NAPAIR disagree, more important to show those
                                pdata['display_position'] = 1

                            # pairs from files that were simply not annotated by some method get
                            pdata['display_position'] -= 0.05 * annotator_available

                            # curated pairs from napair get priority within their display_position
                            u1 = pdata['unit_id_1']
                            u2 = pdata['unit_id_2']
                            if (u1,u2) in curated_pairs:
                                pdata['display_position'] -= 0.4

                            # get second-best hydrogen bond distance and angle
                            # use default target distance and angle
                            dist2_list, hDistAngle_list, extra_penalty, cl, z_score, z_count, ok_hbonds = get_dist2_angle_penalty(interaction_to_atom_sets,pdata,None,None)

                            if not ok_hbonds and not 'dem' in pdata['new_fr3d_detail']:
                                if pdata['cut_dist'] > 0:
                                    # some are directly in a near category, should not be marked
                                    pdata['new_fr3d_detail'] = 'dem,' + pdata['new_fr3d_detail']

                            pdata['ok_hbonds'] = ok_hbonds

                            if len(dist2_list) >= 2:
                                pdata['dist2'] = sorted(dist2_list)[1]
                            elif len(dist2_list) == 1:
                                pdata['dist2'] = dist2_list[0]
                            else:
                                # not clear what to put in this case
                                pdata['dist2'] = 9.99

                            if python_true:
                                all_dist2_list.append(pdata['dist2'])

                            # store this pair
                            pair_data.append(pdata)

                            pair_shown.add(pair)
                            pair_shown.add(reverse(pair))

                            # collect data for scatterplots
                            if show_fr3d_parameters and 'angle_in_plane' in datapoint and abs(datapoint['angle_in_plane']) > 0.0001:
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

                                if angle_out_of_order and datapoint['angle_in_plane'] > 90:
                                    avalues.append(datapoint['angle_in_plane']-360)
                                else:
                                    avalues.append(datapoint['angle_in_plane'])

                                nvalues.append(datapoint['normal_Z'])

                                # colors and sizes for scatterplots
                                if 'demoted_hbond' in datapoint or 'dem' in pdata['new_fr3d_detail']:
                                    color = orange  # orange
                                    size = 5        # medium
                                elif 'gap' in pdata['keep_reasons']:
                                    color = "#ff6db6"  # pink
                                    color = [1,0,0]  # red
                                    size = 5
                                elif len(pdata['keep_reasons']) > 0:
                                    color = cyan
                                    size = 10
                                elif python_true or direct_near:
                                    color = black
                                    size = 1
                                else:
                                    # python near or other or demoted or newly true
                                    color = "#b66dff"  # violet
                                    size = 10       # medium

                                if VERSION in ['v8','v9']:
                                    if python_true and not napair:
                                        color = red
                                        size = 5
                                    elif not python_true and napair:
                                        color = cyan
                                        size = 5
                                    else:
                                        color = black
                                        size = 1

                                colors.append(color)
                                sizes.append(size)

                print('%s %s Python true count %4d Python near count %4d Percentage near %6.2f' % (bc_filename,inter_filename,python_true_count,python_near_count,100.0*python_near_count/(0.00001+python_near_count+python_true_count)))
                print('  pair_data length %3d' % len(pair_data))
                print('  xvalues   length %3d' % len(xvalues))

                # image file name, for writing image and for putting image into HTML file
                if len(all_pdb_ids) <= 10:
                    figure_save_file = os.path.join(OUTPUTPATH,"plots","basepairs_%s_%s_%s.png" % (inter_filename,bc_filename,"_".join(all_pdb_ids)))
                    figure_img_src = "plots/basepairs_%s_%s_%s.png" % (inter_filename,bc_filename,"_".join(all_pdb_ids))
                elif 'nrlist' in PDB_list[0]:
                    figure_save_file = os.path.join(OUTPUTPATH,"plots","basepairs_%s_%s_%s.png" % (inter_filename,bc_filename,resolution))
                    figure_img_src = "plots/basepairs_%s_%s_%s.png" % (inter_filename,bc_filename,resolution)
                elif DNA:
                    figure_save_file = os.path.join(OUTPUTPATH,"plots","DNA_basepairs_%s_%s_%s.png" % (inter_filename,bc_filename,resolution))
                    figure_img_src = "plots/DNA_basepairs_%s_%s_%s.png" % (inter_filename,bc_filename,resolution)
                else:
                    figure_save_file = os.path.join(OUTPUTPATH,"plots","basepairs_%s_%s_%d.png" % (inter_filename,bc_filename,len(all_pdb_ids)))
                    figure_img_src = "plots/basepairs_%s_%s_%d.png" % (inter_filename,bc_filename,len(all_pdb_ids))

                if only_modified:
                    figure_save_file = figure_save_file.replace("basepairs_","modified_basepairs_")
                    figure_img_src = figure_img_src.replace("basepairs_","modified_basepairs_")

                # make scatterplots for pairwise combinations of data
                if len(xvalues) > 0 and show_fr3d_parameters:

                    fig = plt.figure(figsize=(11.0, 7.0))

                    ax = fig.add_subplot(2, 3, 1)
                    ax.axis("equal")
                    plot_basepair_cutoffs(base_combination,[interaction],ax,1)
                    ax.scatter(xvalues,yvalues,color=colors,marker=".",s=sizes)
                    ax.set_title('x, y %d=%d+%d %s %s' % (len(xvalues),python_true_count,python_near_count,base_combination,interaction))
                    draw_base(nt1_seq,'default',2,ax)

                    ax = fig.add_subplot(2, 3, 2)
                    plot_basepair_cutoffs(base_combination,[interaction],ax,2)
                    ax.scatter(tvalues,rvalues,color=colors,marker=".",s=sizes)
                    ax.set_title('radius vs theta')

                    ax = fig.add_subplot(2, 3, 3)
                    plot_basepair_cutoffs(base_combination,[interaction],ax,3,angle_out_of_order)
                    ax.scatter(avalues,nvalues,color=colors,marker=".",s=sizes)
                    ax.set_title('normal vs angle')

                    ax = fig.add_subplot(2, 3, 4)
                    plot_basepair_cutoffs(base_combination,[interaction],ax,4)
                    ax.scatter(gvalues,zvalues,color=colors,marker=".",s=sizes)
                    ax.set_title('z versus maxgap')

                    ax = fig.add_subplot(2, 3, 5)
                    plot_basepair_cutoffs(base_combination,[interaction],ax,5,angle_out_of_order)
                    ax.scatter(avalues,xvalues,color=colors,marker=".",s=sizes)
                    ax.set_title('x vs angle')

                    ax = fig.add_subplot(2, 3, 6)
                    plot_basepair_cutoffs(base_combination,[interaction],ax,6,angle_out_of_order)
                    ax.scatter(avalues,yvalues,color=colors,marker=".",s=sizes)
                    ax.set_title('y vs angle')

                    # show all plots for this interaction
                    figManager = plt.get_current_fig_manager()
                    figManager.full_screen_toggle()

                    plt.savefig(figure_save_file)
                    #plt.show()
                    plt.close()

                    print("Plotted %5d points for %s %s, resolution %s" % (len(xvalues),base_combination,interaction,resolution))

                # a place to store data, similar to a query object in FR3D
                Q = {}

                # write HTML pages listing instances for this base combination and interaction
                # mimic how WebFR3D writes result pages
                if len(all_pdb_ids) <= 10:
                    Q['name'] = "%s %s %s" % (inter_filename,bc_filename,"_".join(all_pdb_ids))
                elif 'nrlist' in PDB_list[0]:
                    Q['name'] = "%s %s %s" % (inter_filename,bc_filename,resolution)
                elif VERSION in ['v9']:
                    Q['name'] = "%s %s %s" % (inter_filename,bc_filename,resolution)
                elif DNA:
                    Q['name'] = "%s %s %s" % (inter_filename,bc_filename,resolution)
                else:
                    Q['name'] = "%s %s %d" % (inter_filename,bc_filename,len(all_pdb_ids))

                if DNA:
                    Q['name'] = 'DNA_' + Q['name']

                if only_modified:
                    Q['name'] = 'modified_' + Q['name']

                Q['resolution'] = resolution
                Q['numFilesSearched'] = len(all_pdb_ids)
                Q['searchFiles'] = all_pdb_ids
                Q['elapsedCPUTime'] = 0
                Q['userMessage'] = []
                Q['figure_img_src'] = figure_img_src
                Q['DNA'] = DNA
                Q['only_modified'] = only_modified
                Q["PDB_data_file"] = data_file

                # for links across the top of the page, tell something about filenames
                if VERSION == 'v9':
                    Q['link_blocks'] = resolution_list
                else:
                    Q['link_blocks'] = resolution_list

                # identify the interaction so a table of HTML links can be made
                if interaction_upper in Leontis_Westhof_basepairs:
                    Q['LW'] = interaction_upper
                else:
                    # reverse the edges to get the standard family name
                    Q['LW'] = interaction_upper[0] + interaction_upper[2] + interaction_upper[1]

                Q['numFound'] = len(pair_data)

                distance_angle_messages = []

                # map atom sets to distance atoms and angle atoms
                Q['d_atoms'] = {}
                Q['a_atoms'] = {}

                Q["target_distance"] = {}
                Q["target_angle"] = {}

                if len(pair_data) > 0:

                    # loop over the hydrogen bonds for this base combination and basepair category
                    # record the atom sets
                    atom_sets = []
                    target_distance = {}
                    target_angle = {}

                    b1,b2 = bc_filename.split(",")

                    if DNA:
                        filename = "DNA_compare_%s_%s_%s,%s_%s.html" % (VERSION,inter_filename,b1,b2,resolution)
                    else:
                        filename = "compare_%s_%s_%s,%s_%s.html" % (VERSION,inter_filename,b1,b2,resolution)

                    family_num = LW_to_number.get(inter_filename,14)   # not sure why this would not work

                    # loop over the hydrogen bonds for this base combination and basepair category
                    for this_interaction in interaction_to_atom_sets.keys():
                        for atom_set in interaction_to_atom_sets[this_interaction]:
                            print('Analyzing interaction %s atom set %s' % (this_interaction,str(atom_set)))

                            d_atoms = atom_set_to_d_atoms[str(atom_set)]
                            a_atoms = atom_set_to_a_atoms[str(atom_set)]

                            # nice way to format the distance atoms
                            distance_atoms = format_distance_atoms(d_atoms,atom_set,nt1_seq,nt2_seq,reverse_bases)
                            angle_atoms = format_angle_atoms(a_atoms,atom_set,nt1_seq,nt2_seq,reverse_bases)

                            # presumably the different interactions have non-overlapping atom sets
                            Q['d_atoms'][atom_set] = distance_atoms
                            Q['a_atoms'][atom_set] = angle_atoms

                            # loop over pdata, accumulate distance and angle from the FR3D-annotated pairs
                            # accumulate a list of bond lengths and angles
                            bond_lengths = []
                            bond_angles = []
                            count = 0

                            # determine mean values among candidates with agreement between annotators
                            num_agree_to_distance = defaultdict(list)
                            num_agree_to_angle = defaultdict(list)
                            for candidate in pair_data:
                                if candidate.get('python_true',False) and candidate['python_annotation'] == this_interaction:
                                    if "atom_set_to_results" in candidate:
                                        num_agree = candidate["annotator_count"]
                                        if atom_set in candidate["atom_set_to_results"]:
                                            num_agree_to_distance[num_agree].append(candidate["atom_set_to_results"][atom_set]["donor_acceptor_distance"])
                                            num_agree_to_angle[num_agree].append(candidate["atom_set_to_results"][atom_set]["heavy_donor_acceptor_angle"])

                            num_agree = 10
                            distances = []
                            angles = []
                            while num_agree > 3 and len(distances) < 5:
                                distances += num_agree_to_distance[num_agree]
                                angles += num_agree_to_angle[num_agree]
                                num_agree -= 1

                            if len(distances) >= 5:
                                target_distance[atom_set] = trim_mean(distances,0.1)
                                target_angle[atom_set] = trim_mean(angles,0.1)
                                distance_angle_messages.append("Target distance for %s bond %10s is %0.2f, target angle for %s is %0.2f, based on %d data points where multiple annotators agree" % (this_interaction,distance_atoms,target_distance[atom_set],angle_atoms,target_angle[atom_set],len(distances)))

                                with open("distance_angle_%s.txt" % Q['resolution'],"a") as outfile:
                                    if "O" in d_atoms:
                                        outfile.write("distance\t%s\t%s\t%0.2f\t%d\tN-O\n" % (distance_atoms,Q["name"],target_distance[atom_set],len(distances)))
                                    else:
                                        outfile.write("distance\t%s\t%s\t%0.2f\t%d\tN-N\n" % (distance_atoms,Q["name"],target_distance[atom_set],len(distances)))
                                    outfile.write("angle\t%s\t%s\t%0.2f\t%d\n" % (angle_atoms,Q["name"],target_angle[atom_set],len(distances)))

                            elif "O" in d_atoms:
                                target_distance[atom_set] = 2.95
                                target_angle[atom_set] = 115.0
                                distance_angle_messages.append("Target distance for %s bond %10s is %0.2f, target angle for %s is %0.2f, based on standard reference" % (this_interaction,distance_atoms,target_distance[atom_set],angle_atoms,target_angle[atom_set]))
                            else:
                                target_distance[atom_set] = 2.8
                                target_angle[atom_set] = 120.0
                                distance_angle_messages.append("Target distance for %s bond %10s is %0.2f, target angle for %s is %0.2f, based on standard reference" % (this_interaction,distance_atoms,target_distance[atom_set],angle_atoms,target_angle[atom_set]))

                            Q["target_distance"] = target_distance
                            Q["target_angle"] = target_angle

                            # record data about hydrogen bonds
                            for candidate in pair_data:
                                if candidate['python_true'] and candidate['python_annotation'] == this_interaction:
                                    if "atom_set_to_results" in candidate and atom_set in candidate["atom_set_to_results"]:
                                        result = candidate["atom_set_to_results"][atom_set]
                                        d_atoms = result["donor_acceptor_atoms"]
                                        a_atoms = result["heavy_donor_acceptor_atoms"]

                                        d = candidate["atom_set_to_results"][atom_set]["donor_acceptor_distance"]
                                        if not math.isnan(d):
                                            count += 1
                                            bond_lengths.append(d)

                                        a = candidate["atom_set_to_results"][atom_set]["heavy_donor_acceptor_angle"]
                                        if not math.isnan(a):
                                            bond_angles.append(a)

                            # write one line for each atom set listing all the data
                            if len(bond_lengths) > 0:
                                d_avg = np.mean(bond_lengths)
                            else:
                                d_avg = 0.0

                            if len(bond_angles) > 0:
                                a_avg = np.mean(bond_angles)
                            else:
                                a_avg = 0.0

                            # check for nan
                            if math.isnan(d_avg):
                                print(bond_lengths)
                                print(crashnow)

                            hbond_data = {}
                            hbond_data['family_num'] = family_num
                            hbond_data['interaction'] = inter_filename
                            hbond_data['b1'] = b1
                            hbond_data['b2'] = b2
                            hbond_data['count'] = count
                            hbond_data['d_atoms'] = distance_atoms
                            hbond_data['d_avg'] = d_avg
                            if len(bond_lengths) > 1:
                                hbond_data['d_std'] = np.std(bond_lengths)
                            else:
                                hbond_data['d_std'] = 0
                            if len(bond_lengths) > 0:
                                hbond_data['percent_over_3.5'] = 100.0*sum([1 for d in bond_lengths if d > 3.5])/len(bond_lengths)
                            else:
                                hbond_data['percent_over_3.5'] = 0.0
                            hbond_data['a_atoms'] = angle_atoms
                            hbond_data['a_avg'] = a_avg
                            hbond_data['a_std'] = np.std(bond_angles)
                            if len(bond_angles) > 0:
                                hbond_data['percent_under_80'] = 100.0*sum([1 for a in bond_angles if a < 80])/len(bond_angles)
                            else:
                                hbond_data['percent_under_80'] = 0.0
                            hbond_data['link'] = "https://rna.bgsu.edu/experiments/annotations/" + filename
                            hbond_data['image_link'] = 'https://www.nakb.org/ndbmodule/bp-catalog//static/img/%s/%s_%s%s_exemplar.png' % (inter_filename,inter_filename,b1,b2)

                            hbond_data_list.append(hbond_data)

                    # record data about second best hydrogen bond distance for this interaction
                    dist2_data = {}
                    dist2_data['family_num'] = family_num
                    dist2_data['interaction'] = inter_filename
                    dist2_data['b1'] = b1
                    dist2_data['b2'] = b2
                    dist2_data['count'] = count
                    dist2_data['python_true_count'] = python_true_count
                    dist2_data['python_near_count'] = python_near_count
                    if len(all_dist2_list) > 0:
                        dist2_data['percent_over_3.7'] = 100.0*sum([1 for d in all_dist2_list if d > 3.7])/len(all_dist2_list)
                        dist2_data['percent_over_3.6'] = 100.0*sum([1 for d in all_dist2_list if d > 3.6])/len(all_dist2_list)
                        dist2_data['percent_over_3.5'] = 100.0*sum([1 for d in all_dist2_list if d > 3.5])/len(all_dist2_list)
                        dist2_data['percent_over_3.4'] = 100.0*sum([1 for d in all_dist2_list if d > 3.4])/len(all_dist2_list)
                    else:
                        dist2_data['percent_over_3.7'] = 0.0
                        dist2_data['percent_over_3.6'] = 0.0
                        dist2_data['percent_over_3.5'] = 0.0
                        dist2_data['percent_over_3.4'] = 0.0
                    dist2_data['link'] = "https://rna.bgsu.edu/experiments/annotations/" + filename
                    dist2_data['image_link'] = 'https://www.nakb.org/ndbmodule/bp-catalog//static/img/%s/%s_%s%s_exemplar.png' % (inter_filename,inter_filename,b1,b2)

                    dist2_data_list.append(dist2_data)

                # write HTML pages listing instances for this base combination and interaction
                if len(pair_data) >= 0 and write_html_pages:

                    # if too many pairs, show some good and also the worst ones
                    if len(pair_data) > 300 and compare_annotators and not show_fr3d_parameters:
                        # sort by full agreement or not, random within each category
                        pair_data = sorted(pair_data, key=lambda p: p['display_position']+0.0001*random.uniform(0,1))

                        # keep 250 with highest priority and 50 of whatever is least
                        order_pair_data = pair_data[0:250]

                        # get the tier 2 pairs by shifting the tier 1 pairs later in the list
                        remaining_pair_data = sorted(pair_data[250:], key=lambda p: p['display_position']+3*(p['display_position']<=1)+0.0001*random.uniform(0,1))

                        # and take 50 of those
                        order_pair_data += remaining_pair_data[0:50]

                        # discard all the rest
                        other_pair_data = []

                    elif len(pair_data) > 300:
                        # sort by true or not, then by second-best hydrogen bond distance
                        pair_data = sorted(pair_data, key=lambda p: (p['python_true'],p['dist2']))

                        # keep 50 of the best pairs, and the 250 worst
                        order_pair_data = pair_data[0:50] + pair_data[-250:]

                        # sort pairs 301 to up to 1300 by python_fr3d annotation and then by badness
                        other_pair_data = sorted(pair_data[51:-251], key=lambda p: (p['python_annotation'],p['dist2']+0.0001*random.uniform(0,1)))

                        # sort pairs 301 to up to 1300 randomly so you get a selection
                        other_pair_data = sorted(pair_data[51:-251], key=lambda p: (random.uniform(0,1),p['python_annotation'],p['dist2']+0.0001*random.uniform(0,1)))

                        if len(other_pair_data) > 1000:
                            # keep 1000 of them, otherwise you get 30,000 GC cWW or something
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

                    print("Finding display order by similarity for %d basepairs" % n)

                    order = treePenalizedPathLength(dista,min(n,100))
                    order = standardOrder(dista,order)

                    if compare_annotators:
                        for i in range(0,n):
                            if VERSION in ['v9']:
                                # indicate curated pairs on the diagonal
                                u1 = order_pair_data[i]['unit_id_1']
                                u2 = order_pair_data[i]['unit_id_2']
                                # print('Checking for %s,%s in curated_pairs' % (u1,u2))
                                if (u1,u2) in curated_pairs:
                                    dista[i][i] = -9
                    else:
                        # color entries on diagonal according to matching annotations
                        for i in range(0,n):
                            if VERSION in ['v8','v9']:
                                if order_pair_data[i]['python_true'] and not order_pair_data[i]['napair_true']:
                                    dista[i][i] = -1  # red for napair maybe being wrong
                                elif not order_pair_data[i]['python_true'] and order_pair_data[i]['napair_true']:
                                    dista[i][i] = -6  # sky blue for python maybe being wrong

                            else:
                                if 'dem' in order_pair_data[i]['new_fr3d_detail']:
                                    dista[i][i] = -9  # orange for bad h-bond
                                elif "gap" in order_pair_data[i]['new_fr3d_detail']:
                                    # dista[i][i] = -2  # dark pink for bad gap
                                    dista[i][i] = -1  # red for bad gap
                                elif "min" in order_pair_data[i]['new_fr3d_detail'] or "max" in order_pair_data[i]['new_fr3d_detail'] or "angle" in order_pair_data[i]['new_fr3d_detail']:
                                    dista[i][i] = -6  # sky blue for other cutoff problem

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

                    # Look up glycosidic bond orientations for the pairs being displayed
                    print('Loading glycosidic bond conformations')
                    QQ = {}
                    QQ['DATAPATHUNITS'] = "C:/Users/zirbel/Documents/PythonFR3D/data/units"
                    for pair in reorder_pairs:
                        u1 = pair['unit_id_1']
                        u2 = pair['unit_id_2']
                        if not u1 in unit_id_to_glycosidic:
                            chain = "|".join(u1.split("|")[0:3])
                            QQ, new_mappings = readUnitAnnotations(QQ,chain)
                            unit_id_to_glycosidic.update(new_mappings)
                        if not u2 in unit_id_to_glycosidic:
                            chain = "|".join(u2.split("|")[0:3])
                            QQ, new_mappings = readUnitAnnotations(QQ,chain)
                            unit_id_to_glycosidic.update(new_mappings)
                        if u1 in unit_id_to_glycosidic and 'orientation' in unit_id_to_glycosidic[u1]:
                            pair['glycosidic1'] = unit_id_to_glycosidic[u1]['orientation']
                        else:
                            pair['glycosidic1'] = ''
                            print('%s has no glycosidic information' % u1)
                        if u2 in unit_id_to_glycosidic and 'orientation' in unit_id_to_glycosidic[u2]:
                            pair['glycosidic2'] = unit_id_to_glycosidic[u2]['orientation']
                        else:
                            pair['glycosidic2'] = ''
                            print('%s has no glycosidic information' % u2)

                    """
                    Write the list of candidates in an HTML format that also shows
                    the coordinate window and a heat map of all-against-all distances.
                    """

                    writeHTMLOutput(Q,reorder_pairs,interaction_to_atom_sets,distance_angle_messages,reorder_dista)

                    show_fr3d_parameters = not show_fr3d_parameters
                    writeHTMLOutput(Q,reorder_pairs,interaction_to_atom_sets,distance_angle_messages,reorder_dista)
                    show_fr3d_parameters = not show_fr3d_parameters


            print('Skipped %5d %s pairs because of symmetry operators' % (len(skip_symmetry_pairs),base_combination))
            print('Skipped %5d %s pairs because of complicated unit ids' % (len(skip_complicated_unit_id),base_combination))
            print('Skipped %5d %s pairs because of being over the target number' % (len(skip_over_target),base_combination))
            print('Skipped %5d %s pairs because not all annotators could process it' % (len(skip_not_all_annotate),base_combination))


        # write a table of counts of each basepair family and base combination
        if compare_annotators:
            output_text = ""

            print("All methods annotate at least one basepair in each of these %d structures:" % len(pdb_id_to_annotators.keys()))
            print(",".join(sorted(pdb_id_to_annotators.keys())))

            # show FR3D near
            header = "Family\tLW\tBase 1\tBase 2\tFR3D\tFR3D near\tRNAview\tRNAview near\tMC-Annotate\tDSSR\tPDB\tNAPAIR\tMax / (1+Min)\tAll\tAll but FR3D\tAll but RNAview\tAll but DSSR\tAll but PDB\tFR3D RNAview\tFR3D DSSR\tFR3D PDB\tRNAview DSSR\tRNAview PDB\tDSSR PDB\tFR3D only\tRNAview only\tDSSR only\tPDB only"

            # leave out FR3D near, add FR3D-All, etc.
            header = "Fam\tLW\tBase 1\tBase 2\tStatus\tFR3D\tRNAview\tRNAview near\tMCA\tDSSR\tPDB\tNAPAIR\tMax / (1+Min)\tIntersection\tFR3D-Int\tRNAview-Int\tMCA-Int\tDSSR-Int\tPDB-Int\tNAPAIR-Int\tAll but FR3D\tAll but RNAview\tAll but MCA\tAll but DSSR\tAll but PDB\tAll but NAPAIR\t4 agree\t3 agree\t2 agree\tFR3D only\tRNAview only\tMCA only\tDSSR only\tPDB only\tNAPAIR only\tInstances\n"

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
                            filename = "DNA_compare_%s_%s_%s,%s_%s.html" % (VERSION,interaction,b1,b2,resolution)
                        else:
                            filename = "compare_%s_%s_%s,%s_%s.html" % (VERSION,interaction,b1,b2,resolution)

                        max_count = 0.0
                        min_count = float('inf')
                        for i,b in enumerate(count_annotations[interaction][base_combination]):
                            if i in [0,2,3,4,5,6,7]: # skip FR3D near
                                output.append(str(b))
                            if i in [0,2,4,5,6,7]:  # skip FR3D near, RNAview near
                                max_count = max(max_count,b)
                                min_count = min(min_count,b)

                        if max_count > 0:
                            output.append("%0.2f" % (float(max_count)/(1.0+min_count)))
                            all_count = count_annotations_by_group[interaction][base_combination]["111111"]

                            output.append(str(all_count)) # All, the number in the intersection

                            for i,b in enumerate(count_annotations[interaction][base_combination]):
                                if i in [0,2,4,5,6,7]:  # skip FR3D near, RNAview near
                                    # FR3D-All, etc.; how much each annotates outside the intersection
                                    output.append(str(b-all_count))

                            output.append(str(count_annotations_by_group[interaction][base_combination]["011111"])) # FR3D-All
                            output.append(str(count_annotations_by_group[interaction][base_combination]["101111"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["110111"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["111011"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["111101"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["111110"]))

                            num_agree = defaultdict(int)
                            for mask, count in count_annotations_by_group[interaction][base_combination].items():
                                num = mask.count('1')
                                num_agree[num] += count

                            output.append(str(num_agree[4]))
                            output.append(str(num_agree[3]))
                            output.append(str(num_agree[2]))

                            output.append(str(count_annotations_by_group[interaction][base_combination]["100000"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["010000"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["001000"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["000100"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["000010"]))
                            output.append(str(count_annotations_by_group[interaction][base_combination]["000001"]))
                            output.append("https://rna.bgsu.edu/experiments/annotations/" + filename)
                            output_list.append(output)

                output_list = sorted(output_list)             # sort by base combination

                for output in output_list:
                    output_text = "\t".join(output)
                    all_output_list.append(output_text+"\n")
                    print(output_text)

                if DNA:
                    table_filename = os.path.join(OUTPUTPATH,'DNA_table_%s_%s.txt' % (VERSION,resolution))
                else:
                    table_filename = os.path.join(OUTPUTPATH,'table_%s_%s.txt' % (VERSION,resolution))

                if only_modified:
                    table_filename = table_filename.replace("table_","modified_table_")

                with open(table_filename,'wt') as f:
                    f.writelines(all_output_list)

                print('Saved %s table in %s' % (resolution,table_filename))

        if write_hbond_tables:
            if DNA:
                hbond_table_file = os.path.join(OUTPUTPATH,"DNA_h-bond_table_%s.txt" % (resolution))
            else:
                hbond_table_file = os.path.join(OUTPUTPATH,"h-bond_table_%s.txt" % (resolution))

            with open(hbond_table_file,"w") as hbond_table:
                hbond_table.write("Family\tLW\tBase 1\tBase 2\tCount\tHbond Atoms\tAvg Distance\tStd Distance\tPercent Over 3.5 Distance\tDonor Angle Atoms\tAvg Angle\tStd Angle\tPercent Under 80 Angle\tEric comment\tOther comment\tImage\tInstances\n")

                # sort by family, base 1, base 2, atoms
                hbond_data_list = sorted(hbond_data_list, key=lambda x: (x['family_num'],x['b1'],x['b2'],x['d_avg']))

                for hbond_data in hbond_data_list:
                    hbond_text = ""
                    hbond_text += "%d\t" % hbond_data['family_num']
                    hbond_text += "%s\t" % hbond_data['interaction']
                    hbond_text += "%s\t" % hbond_data['b1']
                    hbond_text += "%s\t" % hbond_data['b2']
                    hbond_text += "%d\t" % hbond_data['count']
                    hbond_text += "%s\t" % hbond_data['d_atoms']
                    hbond_text += "%0.2f\t" % hbond_data['d_avg']
                    hbond_text += "%0.2f\t" % hbond_data['d_std']
                    hbond_text += "%0.2f\t" % hbond_data['percent_over_3.5']
                    hbond_text += "%s\t" % hbond_data['a_atoms']
                    hbond_text += "%0.2f\t" % hbond_data['a_avg']
                    hbond_text += "%0.2f\t" % hbond_data['a_std']
                    hbond_text += "%0.2f\t" % hbond_data['percent_under_80']
                    hbond_text += "\t\t%s\t" % hbond_data['image_link']
                    hbond_text += "%s" % hbond_data['link']

                    hbond_table.write(hbond_text+"\n")

            print("Saved h-bond table in %s" % hbond_table_file)

        if write_hbond_tables:
            # write table of second best hydrogen bond distances
            if DNA:
                dist2_table_file = os.path.join(OUTPUTPATH,"DNA_dist2_table_%s.txt" % (resolution))
            else:
                dist2_table_file = os.path.join(OUTPUTPATH,"dist2_table_%s.txt" % (resolution))

            with open(dist2_table_file,"w") as dist2_table:
                dist2_table.write("Family\tLW\tBase 1\tBase 2\tTrue\tNear\t% Over 3.4A\t% Over 3.5A\t% Over 3.6A\t% Over 3.7A\tEric comment\tOther comment\tImage\tInstances\n")

                # sort by family, base 1, base 2, atoms
                dist2_data_list = sorted(dist2_data_list, key=lambda x: (x['family_num'],x['b1'],x['b2']))

                for dist2_data in dist2_data_list:
                    dist2_text = ""
                    dist2_text += "%d\t" % dist2_data['family_num']
                    dist2_text += "%s\t" % dist2_data['interaction']
                    dist2_text += "%s\t" % dist2_data['b1']
                    dist2_text += "%s\t" % dist2_data['b2']
                    dist2_text += "%d\t" % dist2_data['python_true_count']
                    dist2_text += "%d\t" % dist2_data['python_near_count']
                    dist2_text += "%0.2f\t" % dist2_data['percent_over_3.4']
                    dist2_text += "%0.2f\t" % dist2_data['percent_over_3.5']
                    dist2_text += "%0.2f\t" % dist2_data['percent_over_3.6']
                    dist2_text += "%0.2f\t" % dist2_data['percent_over_3.7']
                    dist2_text += "\t\t%s\t" % dist2_data['image_link']
                    dist2_text += "%s" % dist2_data['link']

                    dist2_table.write(dist2_text+"\n")

            print("Saved dist2 table in %s" % dist2_table_file)

        if show_fr3d_parameters:
            print("Saved figures in %s/plots" % outputNAPairwiseInteractions)
        if write_html_pages:
            print("Saved HTML files in %s" % OUTPUTPATH)
