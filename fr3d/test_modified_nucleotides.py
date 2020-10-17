# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 12:44:30 2014 @author: Poorna
Name: RNA-protein detection
"""

"""Detect and plot RNA base- amino acid interactions."""
from fr3d.cif.reader import Cif
from fr3d.definitions import RNAconnections
from fr3d.definitions import NAbaseheavyatoms
from fr3d.definitions import NAbasecoordinates
from fr3d.definitions import Ribophos_connect
from fr3d.definitions import aa_connections
from fr3d.definitions import aa_backconnect
from fr3d.definitions import aa_fg
from fr3d.definitions import tilt_cutoff
from fr3d.definitions import planar_atoms
from fr3d.definitions import HB_donors
from fr3d.definitions import HB_acceptors
from fr3d.definitions import modified_nucleotides
import numpy as np
import csv
import urllib
import math

import matplotlib.pyplot as plt
from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D
# note that fr3d.localpath does not synchronize with Git, so you can change it locally to point to your own directory structure
from fr3d.localpath import outputText
from fr3d.localpath import outputBaseAAFG
from fr3d.localpath import inputPath
from fr3d.localpath import outputHTML

from fr3d.ordering.greedyInsertion import orderWithPathLengthFromDistanceMatrix

#from fr3d.classifiers.base_aafg import distance_metrics
from datetime import datetime
from math import floor
import os.path
from os import path

def get_structure(filename):
    if not os.path.exists(filename):
        mmCIFname = filename[-8:]
        print("Downloading "+mmCIFname)
        f = urllib.urlopen("https://files.rcsb.org/download/%s" % mmCIFname)
        myfile = f.read()
        with open(filename, 'w') as outfile:
            outfile.write(myfile)

    with open(filename, 'rb') as raw:
        structure = Cif(raw).structure()
        """All RNA bases are placed in the standard orientation. All Hydrogen
 atoms are inferred. Rotation matrix is calculated for each base."""
#        structure.infer_hydrogens()
        return structure


def atom_dist_basepart(base_residue, aa_residue, base_atoms, c):
    """Calculates atom to atom distance of part "aa_part" of neighboring amino acids
    of type "aa" from each atom of base. Only returns a pair of aa/nt if two
    or more atoms are within the cutoff distance"""
    min_distance = 4
    n = 0
    for base_atom in base_residue.atoms(name=base_atoms):
        for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
            # aa_atom = atom.coordinates()
            distance = np.subtract(base_atom.coordinates(), aa_atom.coordinates())
            distance = np.linalg.norm(distance)
            #print base_residue.unit_id(), aa_residue.unit_id(), distance
            if distance <= min_distance:
                n = n+1
                #print aa_atom.name
    if n>=c:
        #print aa_residue.unit_id()
        return True

def enough_HBs(base_residue, aa_residue, base_atoms):
    """Calculates number of Hydrogen bonds between amino acid part and base_part
    and determines if they are enough to form a pseudopair"""
    min_distance = 4
    base_key = base_residue.sequence
    aa_key = aa_residue.sequence
    n = 0
    base_HB_atoms = base_atoms + ["O2'"]
    base_donors = HB_donors[base_key]
    base_acceptors = HB_acceptors[base_key]
    aa_donors = HB_donors[aa_key]
    aa_acceptors = HB_acceptors[aa_key]

    for base_atom in base_residue.atoms(name=base_HB_atoms):
        for aa_atom in aa_residue.atoms(name=aa_fg[aa_key]):
            distance = np.subtract(base_atom.coordinates(), aa_atom.coordinates())
            distance = np.linalg.norm(distance)
            if distance <= min_distance:
                #print "HB", base_residue.unit_id(), aa_residue.unit_id(), base_atom.name, aa_atom.name, distance
                if base_atom.name in base_donors and aa_atom.name in aa_acceptors:
                    n = n+1
                elif base_atom.name in base_acceptors and aa_atom.name in aa_donors:
                    n = n+1
#    print base_residue.unit_id(), aa_residue.unit_id(), n
    if n>=2:
        return True

def find_neighbors(bases, amino_acids, aa_part, dist_cent_cutoff):
    """Finds all amino acids of type "aa" for which center of "aa_part" is within
    specified distance of center of bases of type "base" """

    # build a set of cubes and record which bases are in which cube
    # also record which other cubes are neighbors of each cube
    baseCubeList = {}
    baseCubeNeighbors = {}
    for base in bases:
        center = base.centers["base"]
        if len(center) == 3:
#            print(base.unit_id() + str(center))
            x = floor(center[0]/dist_cent_cutoff)
            y = floor(center[1]/dist_cent_cutoff)
            z = floor(center[2]/dist_cent_cutoff)
            key = "%d,%d,%d" % (x,y,z)
            if key in baseCubeList:
                baseCubeList[key].append(base)
            else:
                baseCubeList[key] = [base]
                baseCubeNeighbors[key] = []
                for a in [-1,0,1]:
                    for b in [-1,0,1]:
                        for c in [-1,0,1]:
                            k = "%d,%d,%d" % (x+a,y+b,z+c)
                            baseCubeNeighbors[key].append(k)

    # build a set of cubes and record which amino acids are in which cube
    aaCubeList = {}
    for aa in amino_acids:
        center = aa.centers[aa_part]
        if len(center) == 3:
            x = floor(center[0]/dist_cent_cutoff)
            y = floor(center[1]/dist_cent_cutoff)
            z = floor(center[2]/dist_cent_cutoff)
            key = "%d,%d,%d" % (x,y,z)
            if key in aaCubeList:
                aaCubeList[key].append(aa)
            else:
                aaCubeList[key] = [aa]
#        else:
#            print("  Missing center coordinates for " + str(aa))

    return baseCubeList, baseCubeNeighbors, aaCubeList

def annotate_interactions(bases, amino_acids, aa_part, dist_cent_cutoff, baseCubeList, baseCubeNeighbors, aaCubeList):

    # loop through base cubes, loop through neighboring cubes,
    # then loop through bases and amino acids in the two cubes,
    # screening distances between them, then annotating interactions
    """Finds all amino acids of type "aa" for which center of "aa_part" is within
    specified distance of center of bases of type "base" and returns superposed bases"""

    #count_total = 0
    count_pair = 0
    list_aa_coord = []
    list_base_coord = []
    aaList_len = None
    new_aaList_len = None
    list_base_aa = []

    for key in baseCubeList:
        for aakey in baseCubeNeighbors[key]:
            if aakey in aaCubeList:
                for base_residue in baseCubeList[key]:
                    base_seq = base_residue.sequence
                    base_atoms = NAbaseheavyatoms[base_seq]

                    base_center = base_residue.centers["base"]

#                    print("Base center")
#                    print base_center

                    if not base_center.any():
                        continue

                    if not base_residue.centers[planar_atoms[base_seq][0]].any():
                        continue
                    if not base_residue.centers[planar_atoms[base_seq][1]].any():
                        continue
                    if not base_residue.centers[planar_atoms[base_seq][2]].any():
                        continue

                    for aa_residue in aaCubeList[aakey]:
                        aa_center = aa_residue.centers[aa_part]
                        if not aa_center.any():
                            continue

                        if aa_residue.sequence in set (['LYS','SER', 'THR', 'TYR']):
                            c= 1
                        else:
                            c = 2

                        if abs(base_center[0]-aa_center[0]) < dist_cent_cutoff and \
                        abs(base_center[1]-aa_center[1]) < dist_cent_cutoff and \
                        np.linalg.norm(np.subtract(base_center,aa_center)) < dist_cent_cutoff and \
                        atom_dist_basepart(base_residue, aa_residue, base_atoms, c):

                            count_pair = count_pair + 1

                            rotation_matrix = base_residue.rotation_matrix

                            # rotate base atoms into standard orientation
                            base_coordinates = {}
                            centered_base = base_residue.translate_rotate_component(base_residue)
                            for base_atom in centered_base.atoms():
                                base_coordinates[base_atom.name]= base_atom.coordinates()

                            # rotate amino acid atoms into standard orientation
                            aa_coordinates = {}
                            standard_aa = base_residue.translate_rotate_component(aa_residue)
                            for aa_atom in standard_aa.atoms():
                                aa_coordinates[aa_atom.name] = aa_atom.coordinates()

                            standard_aa_center = standard_aa.centers[aa_part]

#                            print aa_residue
#                            print base_residue
                            interaction = type_of_interaction(base_residue, aa_residue, aa_coordinates)

                            base_aa = None
                            if interaction == "pseudopair" and enough_HBs(base_residue, aa_residue, base_atoms):
                                edge = detect_edge(base_residue, base_coordinates,aa_residue, aa_coordinates)
                                base_aa = (base_residue, aa_residue, interaction, edge, standard_aa_center)

                            elif interaction == "SHB":
                                edge = detect_edge(base_residue, base_coordinates,aa_residue, aa_coordinates)
                                base_aa = (base_residue, aa_residue, interaction, edge, standard_aa_center)

                            elif interaction == "perpendicular edge":
                                edge = detect_edge(base_residue, base_coordinates,aa_residue, aa_coordinates)
                                base_aa = (base_residue, aa_residue, interaction, edge, standard_aa_center)

                            elif interaction == "stacked" or interaction == "cation-pi" \
                            or interaction == "perpendicular stacking":
                                edge = detect_face(aa_residue, aa_coordinates)
                                base_aa = (base_residue, aa_residue, interaction, edge, standard_aa_center)

                            if base_aa is not None:
                                list_base_aa.append(base_aa)

                                for base_atom in base_residue.atoms():
                                    list_base_coord.append(base_coordinates)
                                for aa_atom in aa_residue.atoms():
                                    list_aa_coord.append(aa_coordinates)

    return list_base_aa, list_aa_coord, list_base_coord

def type_of_interaction(base_residue, aa_residue, aa_coordinates):
    squared_xy_dist_list = []
    aa_z =[]

    """Defines different sets of amino acids"""
    stacked_planar_aa = set (["TRP", "TYR", "PHE", "HIS", "ARG", "ASN", "GLN", "GLU", "ASP"])
    stacked_aliphatic = set(["LEU", "ILE", "PRO", "THR", "MET", "CYS", "VAL", "ALA", "SER"])
    pseudopair_aa = set (["ASP", "GLU", "ASN", "GLN", "HIS", "ARG", "TYR", "TRP", "PHE", "LYS"])
    shb_aa = set (["SER", "THR", "LYS"])

    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key = aa_atom.name
        aa_x= aa_coordinates[key][0]
        aa_y= aa_coordinates[key][1]

        squared_xy_dist = (aa_x**2) + (aa_y**2)
        squared_xy_dist_list.append(squared_xy_dist)

        aa_z.append(aa_coordinates[key][2])

    mean_z = np.mean(aa_z)

    #print base_residue.unit_id(), aa_residue.unit_id(), min(squared_xy_dist_list), mean_z
    if min(squared_xy_dist_list) <= 5:
        #print base_residue.unit_id(), aa_residue.unit_id(), min(squared_xy_dist_list), mean_z
        if aa_residue.sequence in stacked_planar_aa:
            #print "stacking?", base_residue.unit_id(), aa_residue.unit_id(), min(squared_xy_dist_list), mean_z
            return stacking_angle(base_residue, aa_residue, min(squared_xy_dist_list))

        elif aa_residue.sequence in stacked_aliphatic:
            return stacking_tilt(aa_residue, aa_coordinates)

    elif -1.8 <= mean_z < 1.8 and aa_residue.sequence in pseudopair_aa:
            angle= calculate_angle(base_residue, aa_residue)
            angle = abs(angle)
            #print "pseudopair?", base_residue.unit_id(), aa_residue.unit_id(), angle
            if 0 <= angle <= 0.75 or 2.5 <= angle <= 3.14:
                return "pseudopair"
            elif 0.95 <= angle <=1.64:
                return "perpendicular edge"

    elif -1.8 <= mean_z < 1.8 and aa_residue.sequence in shb_aa:
        base_seq = base_residue.sequence
        base_atoms = NAbaseheavyatoms[base_seq]
        if atom_dist_basepart(base_residue, aa_residue, base_atoms, 1):
            return "SHB"

def calculate_angle (base_residue, aa_residue):
    vec1 = vector_calculation(base_residue)
    vec2 = vector_calculation(aa_residue)

    angle = angle_between_planes(vec1, vec2)
    return angle

def stacking_angle (base_residue, aa_residue, min_dist):
    vec1 = vector_calculation(base_residue)
    vec2 = vector_calculation(aa_residue)
    perpendicular_aa = set (["HIS", "ARG", "LYS", "ASN", "GLN"])
    perpendicular_stack_aa = set(["PHE", "TYR"])
    angle = angle_between_planes(vec1, vec2)
    #print "stacked?"
    #print base_residue.unit_id(), aa_residue.unit_id(), min_dist, angle
    angle = abs(angle)
    if angle <=0.68 or 2.45 <= angle <= 3.15:
        return "stacked"
    elif 1.2<= angle <=1.64:
        if aa_residue.sequence in perpendicular_stack_aa:
            return "perpendicular stacking"
        elif aa_residue.sequence in perpendicular_aa:
            return "cation-pi"

def stacking_tilt(aa_residue, aa_coordinates):
    baa_dist_list = []

    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key = aa_atom.name
        aa_z = aa_coordinates[key][2]
        baa_dist_list.append(aa_z)
    max_baa = max(baa_dist_list)
    min_baa = min(baa_dist_list)
    diff = max_baa - min_baa
    #print aa_residue.unit_id(), diff
    if diff <= tilt_cutoff[aa_residue.sequence]:
        return "stacked"

def vector_calculation(residue):
    key = residue.sequence
    P1 = residue.centers[planar_atoms[key][0]]
    P2 = residue.centers[planar_atoms[key][1]]
    P3 = residue.centers[planar_atoms[key][2]]
#    print key,P1, P2, P3

    vector = np.cross((P2 - P1),(P3-P1))
    return vector

def angle_between_planes (vec1, vec2):
    cosang = np.dot(vec1, vec2)
    sinang = np.linalg.norm(np.cross(vec1, vec2))
    angle = np.arctan2(sinang, cosang)
    return angle

def detect_edge(base_residue, base_coordinates,aa_residue, aa_coordinates):
    aa_x = []
    aa_y = []
    base_x = []
    base_y = []
    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key = aa_atom.name
        aa_x.append(aa_coordinates[key][0])
        aa_y.append(aa_coordinates[key][1])

    aa_center_x = np.mean(aa_x)
    aa_center_y = np.mean(aa_y)

    for base_atom in base_residue.atoms(name=NAbaseheavyatoms[base_residue.sequence]):
        key = base_atom.name
        base_x.append(base_coordinates[key][0])
        base_y.append(base_coordinates[key][1])

    base_center_x = np.mean(base_x)
    base_center_y = np.mean(base_y)

    y = aa_center_y - base_center_y
    x = aa_center_x - base_center_x
    angle_aa = np.arctan2(y,x) #values -pi to pi
    #print base_residue.unit_id(), aa_residue.unit_id(),angle_aa

    purine = set(["A", "G","DA","DG"])
    pyrimidine = set(["C", "U","DC"])
    angle_deg = (180*angle_aa)/3.14159 #values -180 to 180

#    print("  Edge angle in rad and deg" +str(angle_aa) + " " + str(angle_deg))

    if base_residue.sequence in purine:
        if -15 <= angle_deg <= 90:
            return "fgWC"
        elif 90 < angle_deg or angle_deg < -160:
            return "fgH"
        else:
            return "fgS"

    elif base_residue.sequence in pyrimidine:
        if -45 <= angle_deg <= 90:
            return "fgWC"
        elif 90 < angle_deg or angle_deg < -150:
            return "fgH"
        else:
            return "fgS"

def detect_face(aa_residue, aa_coordinates):
    aa_z =[]

    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key = aa_atom.name
        aa_z.append(aa_coordinates[key][2])

    mean_z = np.mean(aa_z)
    if mean_z <= 0:
        return "fgs5"
    else:
        return "fgs3"

def text_output(result_list):
    with open(outputText % PDB, 'wb') as target:
        for result in result_list:
            target.write(str(result))
            target.write("\r\n")
            target.close

def csv_output(result_list):
    with open(outputBaseAAFG % PDB, 'wb') as csvfile:
        fieldnames = ['RNA ID', 'AA ID', 'RNA Chain ID', 'RNA residue','RNA residue number','Protein Chain ID', 'AA residue','AA residue number', 'Interaction', 'Edge']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for base_residue, aa_residue, interaction, edge, standard_aa_center in result_list:
            base = base_residue.unit_id()
            aa = aa_residue.unit_id()
            #print base, aa, interaction
            base_component = str(base).split("|")
            aa_component = str(aa).split("|")
            writer.writerow({'RNA ID': base, 'AA ID': aa, 'RNA Chain ID': base_component[2], 'RNA residue':base_component[3],'RNA residue number': base_component[4],'Protein Chain ID':aa_component[2],'AA residue': aa_component[3],'AA residue number': aa_component[4], 'Interaction': interaction, 'Edge': edge})

        """for base_residue, aa_residue,interaction in result_list:
                    base_component = str(base_residue).split("|")
                    aa_component = str(aa_residue).split("|")
                    writer.writerow({'RNA Chain ID': base_component[2], 'RNA residue':base_component[3],'RNA residue number': base_component[4],'Protein Chain ID':ChainNames[PDB][aa_component[2]],'AA residue': aa_component[3],'AA residue number': aa_component[4], 'Interaction': interaction})"""



def draw_base(base_seq, ax):
    """Connects atoms to draw neighboring bases and amino acids for 3D plots"""
     #creates lists of rotated base coordinates
    for basecoord_list in list_base_coord:
        new_base_x = []
        new_base_y = []
        new_base_z = []

        back_base_x = []
        back_base_y = []
        back_base_z = []


        try:
            for atomname in RNAconnections[base_seq]:
                coord_base = []
                coord_base= basecoord_list[atomname]
                new_base_x.append(coord_base[0])
                new_base_y.append(coord_base[1])
                new_base_z.append(coord_base[2])
            base_lines= ax.plot(new_base_x, new_base_y, new_base_z, label= 'Base')
            #ax.scatter(basecenter[0], basecenter[1], basecenter[2], zdir='y', color='b', marker='o')
            #ax.scatter(x = 0, y= 0, z= 0, color='b', marker='o')
            plt.setp(base_lines, 'color', 'b', 'linewidth', 1.0)

            for atomname in Ribophos_connect[base_seq]:
                back_base=[]
                back_base= basecoord_list[atomname]
                back_base_x.append(back_base[0])
                back_base_y.append(back_base[1])
                back_base_z.append(back_base[2])
            base_lines= ax.plot(back_base_x, back_base_y, back_base_z, label= 'Base')
            plt.setp(base_lines, 'color', 'r', 'linewidth', 1.0)
            #ax.text(9, 1, 1, base_residue)
        except:
            print "Missing residues"
            continue

def draw_aa(aa, ax):
    #Connects atoms to draw neighboring bases and amino acids for 3D plots
    for aacoord_list in list_aa_coord:
        new_aa_x=[]
        new_aa_y=[]
        new_aa_z=[]

        back_aa_x=[]
        back_aa_y=[]
        back_aa_z=[]

        try:
            for atomname in aa_connections[aa]:
                coord_aa=[]
                coord_aa= aacoord_list[atomname]
                new_aa_x.append(coord_aa[0])
                new_aa_y.append(coord_aa[1])
                new_aa_z.append(coord_aa[2])
            aa_lines= ax.plot(new_aa_x, new_aa_y, new_aa_z, label= 'Amino acid')
            plt.setp(aa_lines, 'color', 'r', 'linewidth', 1.0)

            for atomname in aa_backconnect[aa]:
                back_aa=[]
                back_aa= aacoord_list[atomname]
                back_aa_x.append(back_aa[0])
                back_aa_y.append(back_aa[1])
                back_aa_z.append(back_aa[2])
            aa_lines= ax.plot(back_aa_x, back_aa_y, back_aa_z, label= 'Amino acid')
            plt.setp(aa_lines, 'color', 'y', 'linewidth', 1.0)
        except:
            print "Missing residues"
            continue

def draw_aa_cent(aa, aa_part, ax):
    #Connects atoms to draw neighboring bases and amino acids for 3D plots
    for aacoord_list in list_aa_coord:
        new_aa_x=[]
        new_aa_y=[]
        new_aa_z=[]

        aa_center_x = 0
        aa_center_y = 0
        aa_center_z = 0
        n = 0

        if aa_part == 'aa_fg':
            connections = aa_connections
        elif aa_part == 'aa_backbone':
            connections = aa_backconnect
        try:
            for atomname in connections[aa]:
                coord_aa=[]
                coord_aa= aacoord_list[atomname]
                new_aa_x.append(coord_aa[0])
                new_aa_y.append(coord_aa[1])
                new_aa_z.append(coord_aa[2])

                aa_center_x = aa_center_x + coord_aa[0]
                aa_center_y = aa_center_y + coord_aa[1]
                aa_center_z = aa_center_z + coord_aa[2]
                n = n + 1
            ax.scatter(aa_center_x/n, aa_center_y/n, aa_center_z/n, c= 'r', marker = 'o')
        except:
            print "Missing residues"
            continue

def writeInteractionsHTML(allInteractionDictionary,outputHTML):

    for key in allInteractionDictionary:
        pagetitle = key.replace(" ","-")
        htmlfilename = key.replace(" ","-")

#        print("Writing HTML file for "+key+" in "+htmlfilename+".html, found "+ str(len(allInteractionDictionary[key])) + " instances")

        # limit the number of instances shown, to be able to compute and display discrepancy
        numForDiscrepancy = min(300,len(allInteractionDictionary[key]))

        fields = key.split("_")
        print(fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t\t"+str(numForDiscrepancy)+"\t\t"+"http://rna.bgsu.edu/RNAprotein/"+key.replace(" ","-")+".html")

        # calculate discrepancies between all instances, up to 300
        discrepancy = np.zeros((numForDiscrepancy,numForDiscrepancy))
        for i in range(0,numForDiscrepancy):
            instance_1 = allInteractionDictionary[key][i]
            standard_aa_center_1 = instance_1[4]
            for j in range(i+1,numForDiscrepancy):
                instance_2 = allInteractionDictionary[key][j]
                standard_aa_center_2 = instance_2[4]
                discrepancy[i][j] = np.linalg.norm(standard_aa_center_1 - standard_aa_center_2)
                discrepancy[j][i] = np.linalg.norm(standard_aa_center_1 - standard_aa_center_2)

        # base_aa = (base_residue, aa_residue, interaction, edge, standard_aa_center)

        # use greedy insertion 100 times to find a decent ordering of the instances
        newOrder, bestPathLength, distances = orderWithPathLengthFromDistanceMatrix(discrepancy,100)

        # rewrite the list of instances, limiting it to numForDiscrepancy
        newList = []
        for i in range(0,len(newOrder)):
            newList.append(allInteractionDictionary[key][newOrder[i]])
        allInteractionDictionary[key] = newList

        # write out text for radio boxes to display each individual interaction
        instancelist = "<h2>"+pagetitle+"</h2>\n<ol>"
        i = 1
        for base_id, aa_id, interaction, edge, standard_aa_center in allInteractionDictionary[key]:
            instancelist += '<li><label><input type="checkbox" id="'+str(i)+'" class="jmolInline" data-coord="'
            instancelist += base_id +","+ aa_id
            instancelist += '">&nbsp'
            instancelist += base_id +" "+ aa_id +" "+ interaction +" "+ edge
            instancelist += '</label></li>\n'
            i += 1
        instancelist += '</ol>\n'

        # write out text to tell what values to put in the heat map
        discrepancyText = ''
        for c in range(0,numForDiscrepancy):
            instance1 = allInteractionDictionary[key][c][0]  # id of base
            for d in range(0,numForDiscrepancy):
                instance2 = allInteractionDictionary[key][d][0]  # id of base

                discrepancyText += '{"discrepancy": ' + str(discrepancy[newOrder[c]][newOrder[d]])
                discrepancyText += ', "ife1": "' + str(c+1) + "-" + instance1 + '", "ife1_index": ' + str(c)
                discrepancyText += ', "ife2": "' + str(d+1) + "-" + instance2 + '", "ife2_index": ' + str(d) + '}'
                if c < numForDiscrepancy-1 or d < numForDiscrepancy-1:
                    discrepancyText += ',\n'

        # read template.html into one string
#        with open(outputHTML+'/localtemplate.html', 'r') as myfile:
        with open(outputHTML+'/onlinetemplate.html', 'r') as myfile:
            template = myfile.read()

        # replace ###PAGETITLE### with pagetitle
        template = template.replace("###PAGETITLE###",pagetitle)

        # replace ###CANDIDATELIST### with instancelist
        template = template.replace("###CANDIDATELIST###",instancelist)

        # replace ###DISCREPANCYDATA### with discrepancyText
        template = template.replace("###DISCREPANCYDATA###",discrepancyText)

        # write htmlfilename
        with open(outputHTML+'/'+htmlfilename+'.html', 'w') as myfile:
            myfile.write(template)

        # upload the files to /var/www/html/RNAprotein

def PlotGivenBase(GivenBase,ax):

    for base_atom_1 in GivenBase.atoms():
        if (not "'" in base_atom_1.name and not "P" in base_atom_1.name) or "C1'" in base_atom_1.name:
            coord1 = base_atom_1.coordinates()
            for base_atom_2 in GivenBase.atoms():
                if (not "'" in base_atom_2.name and not "P" in base_atom_2.name) or "C1'" in base_atom_2.name:
                    coord2 = base_atom_2.coordinates()

                    xvalues = [coord1[0], coord2[0]]
                    yvalues = [coord1[1], coord2[1]]
                    zvalues = [coord1[2], coord2[2]]
                    distance = math.sqrt((xvalues[0]-xvalues[1])**2 + (yvalues[0]-yvalues[1])**2 + (zvalues[0]-zvalues[1])**2)
                    if 0 < distance and distance < 1.7:
                        #print("Given red atom names:" + base_atom_1.name + " " + base_atom_2.name + " distance %8.4f" % distance)
                        base_lines= ax.plot(xvalues, yvalues, zvalues, label= 'Base')
                        #ax.scatter(basecenter[0], basecenter[1], basecenter[2], zdir='y', color='b', marker='o')
                        #ax.scatter(x = 0, y= 0, z= 0, color='b', marker='o')
                        plt.setp(base_lines, 'color', 'b', 'linewidth', 1.0)

def PlotStandardBase(sequence,ax):

    for base_atom_1 in NAbasecoordinates[sequence].keys():
        if (not "'" in base_atom_1 and not "P" in base_atom_1) or "C1'" in base_atom_1:
            coord1 = NAbasecoordinates[sequence][base_atom_1]
            for base_atom_2 in NAbasecoordinates[sequence].keys():
                if (not "'" in base_atom_2 and not "P" in base_atom_2) or "C1'" in base_atom_2:
                    coord2 = NAbasecoordinates[sequence][base_atom_2]

                    xvalues = [coord1[0], coord2[0]]
                    yvalues = [coord1[1], coord2[1]]
                    zvalues = [coord1[2], coord2[2]]
                    distance = math.sqrt((xvalues[0]-xvalues[1])**2 + (yvalues[0]-yvalues[1])**2 + (zvalues[0]-zvalues[1])**2)
                    if 0 < distance and distance < 1.7:
                        #print("Standard blue atom names:" + base_atom_1 + " " + base_atom_2 + " distance %8.4f" % distance)
                        base_lines= ax.plot(xvalues, yvalues, zvalues, label= 'Base')
                        #ax.scatter(basecenter[0], basecenter[1], basecenter[2], zdir='y', color='b', marker='o')
                        #ax.scatter(x = 0, y= 0, z= 0, color='b', marker='o')
                        plt.setp(base_lines, 'color', 'r', 'linewidth', 1.0)



#=======================================================================

"""Inputs a list of PDBs of interest to generate super-imposed plots"""
PDB_List = ['5AJ3']
PDB_List = ['6hiv']
PDB_List = ['3QRQ','5J7L']
PDB_List = ['5J7L']
PDB_List = ['4V9F','4YBB','4Y4O','6AZ3','4P95']
PDB_List = ['3BT7']
PDB_List = ['5I4A']
PDB_List = ['4V9F']
PDB_List = ['6A2H','6A2I','6IDE','6N60','6N61','6N62','5YUU','6DKS','6E8C','6GVQ','6GVT','6GVT','6GVU','5YUX','6BHX','6FWR','6FWS','6MG2','6MG3','6BQU','6BSE','6BSF','6FBQ','6FBR','6IG1','5YTC','5YTD','5YTE','5YTF','5YTG','5YTH','5Z3N','5ZVA','5ZVB','5YBD','5YUZ','5YV0','5YUS','5YUW','5YV3','5ZLV','6A47','6A4B','6E93','6E94','5WCU','5Z6Z','5ZFW','5ZFY','5ZFZ']
PDB_List = ['6BHJ','6AI6','6CVO','6DOI','6DPB','6DOF','6DP4','6DO8','6DOX','6DOQ','6DPM','6DPF','6DOJ','6DP8','6DOC','6DP1','6DOU','6DON','6DPJ','6DPC','6DOG','6DP5','6DO9','6DOY','6DOR','6DPN','6DPG','6DOK','6DP9','6DOD','6DP2','6DMN','6DOV','6DOO','6DPK','6DPD','6DOH','6DP6','6DOA','6DOZ','6DOS','6DPO','6DOL','6DPH','6DPA','6DOE','6DP3','6DMV','6DOW','6DOP','6DPL','6DPE','6DP7','6DOB','6DP0','6DPP','6DOT','6DOM','6DPI','5ZRF','6D95','5ZQF','6D8P','6D9L','6D92','6D8A','6D8F','6D9K','6DT8','6DTA','6CVT','6CVP','6CVQ','6CVR','5W4U','5W51','5XPA','5USB','5XPG','5USN','5USO','6BM2','6BM4','6BLO','6BQF','6BLP','5XN0','5XN2','5XMA','6BSH','6BSI','6BSJ','6BSG','6AR3','6AR1','5OT2','5WTI','5OLA','5XOW','5VAJ','5KW1','5XOG','5O6U','5XUT','5XUU','5WJR','5XUZ','5XUS','5VI5','5X21','5VZH','5VZ8','5X22','5VZB','5TWS','5VZE','5W7O','5W7N','5MGA','5N9G','5XH6','5NFV','5XH7','5VO8','5UX0','5UHC','5UH6','5UH8','5UH9','5UH5','5X2G','5X2H','5U30','5U31','5U33','5SWM','5KK5','5AWH','5I2D','5FW2','5FW3','5B43','5FW1','4Z7K','5CR2','5IPL','5IPM','5IPN','5B2S','5B2T','5FQ5','5F0Q','5B2R','5F0S','5E18','5E17','5B2P','5B2Q','5B2O','5EV3','5EV4','5EV1','5EV2','5H9E','5H9F','5F9R','4XLN','4XLR','5CZZ','5AXW','5C4X','5C44','5C4A','5C4J','5C3E','4Y7N','4Y52','3X1L','4WB3','4PGY','4YLN','4YLO','4YLP','4WQS','4X67','4X6A','4S20','4Q5V','4QCL','4QYZ','4PY5','4UN3','4UN4','4UN5','4Q5S','4PQU','4PUO','4Q0B','4PWD','4PUQ','4O9M','4OL8','4OO8','4NDF','4NDG','4NDI','4KHW','4KHY','4KHS','4KI4','4KHU','4KI6','4H8K','4BOC','4BWM','4BXX','4BY1','4BY7','4K4Y','4K4V','4FXD','4FYD','4HKQ','4GZY','4GZZ','4HHT','4B3Q','4B3O','4B3P','3ULD','4BBS','4G7O','4GG4','3TWH','4DB4','4FO6','4DQS','3UQ2','3UQ0','4E7A','4A93','4A3G','4A3M','4A3D','4A3J','4A3E','4A3K','4A3B','4A3F','4A3L','4A3C','3S2H','3S1Q','3RZO','3S17','3S1R','3S14','3S1M','3S2D','3S15','3S1N','3RZD','3S16','3QRQ','3PO3','3PO2','3PGW','3OLA','3O3H','3AOH','3AOI','3O3F','3O3G','3M3Y','3M4O','3KYL','3IIN','3I55','3KJO','3HXM','3HK2','3HM9','3HVR','3HO1','3HJF','3I4M','3I4N','3HOU','3HOY','3HOV','3HOZ','3HOW','3HOX','3H3V','3ER8','3GTP','3GTL','3GTQ','3GTG','3GTM','3GTJ','3GTO','3GTK','3F73','3E2E','3E3J','2VUM','2R7Y','3BO4','3BSU','2R7Z','2QK9','2QKB','2QKK','2PPB','2Q2T','2Q2U','2O5J','2O5I','2YU9','2JA7','2JA8','2JA5','2JA6','2NVX','2E2I','2NVZ','2E2H','2NVQ','2E2J','2NVT','2HVS','2HVR','2G8F','2G8U','2G8H','2G8V','2G8I','2G8W','2G8K','1ZBI','1ZBL','1Y77','1Y1W','1Y1Y','1R9T','1R9S','1SI2','1S76','1S77','1SFO','1S0V','1Q7Y','1NH3','1H38','1MSW','1I6H','1HYS','1QLN','1D9D','1D9F']
PDB_List = ['4NLF']
PDB_List = ['1SK5','5XUT','4NLF','4X1A','1EM0','1DPL']
PDB_List = ['4RKV'] # 0.9 Angstrom RNA structure, has A and B locations that confuse fr3d-python
PDB_List = ['4TNA'] # RNA structure with several modified modified_nucleotides
PDB_List = ['5XUT'] # 5XUT|1|A|ARG|1076 has no sidechain
PDB_List = ['6UKF'] # DNA structure
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/2.5A/csv']

base_seq_list = ['A','U','C','G']      # for RNA
base_seq_list = ['DT']  # for DNA DT only
base_seq_list = ['DA','DC','DG','DT']  # for DNA
base_seq_list = modified_nucleotides.keys()  # to study known modified bases

plot_nts = False
plot_nts = True

aa_list = ['ALA','VAL','ILE','LEU','ARG','LYS','HIS','ASP','GLU','ASN','GLN','THR','SER','TYR','TRP','PHE','PRO','CYS','MET']
#aa_list = ['HIS']


"""Inputs base, amino acid, aa_part of interest and cut-off distance for subsequent functions"""
if __name__=="__main__":

    aa_part = 'aa_fg'
    base_part = 'base'

    MinMaxZ = 50       # largest z value for base atoms in one base, but smallest value seen so far

    allInteractionDictionary = defaultdict(list)
    result_nt_aa = []               # for accumulating a complete list over all PDB files

    PDB_File_List = []              # accumulate all PDB file names
    for PDB in PDB_List:
        if "nrlist" in PDB:           # referring to a list online
            f = urllib.urlopen(PDB)
            myfile = f.read()
            alllines = myfile.split("\n")
            for line in alllines:
                fields = line.split(",")

                if len(fields) > 1 and len(fields[1]) > 4:
                    newPDB = fields[1][1:5]
                    PDB_File_List.append(newPDB)
        else:
            PDB_File_List.append(PDB)

    PDB_File_List = sorted(list(set(PDB_File_List)))

    heavyatoms = ['N1','C2','O2','N3','C4','O4','C6','C5','C7']  # not a complete list ...

    bases_plotted = []

    counter = 0
    for PDB in PDB_File_List:
        counter += 1

        print("Reading file " + PDB + ", which is number "+str(counter)+" out of "+str(len(PDB_File_List)))

        start = datetime.now()
        structure = get_structure(inputPath % PDB)  # this is a local function
        bases = structure.residues(sequence=base_seq_list)   # show all standard RNA bases
#        bases = structure.residues()   # show all RNA bases, maybe other units as well

        amino_acids = structure.residues(sequence=aa_list)
        print("Time required to load " + PDB + " " + str(datetime.now() - start))

        numBases = 0
        for base in bases:
            numBases += 1
            print("Analyzing %s",base.unit_id())
            fields = base.unit_id().split("|")
            if base.base_center is None or base.rotation_matrix is None:
                print("No center or rotation matrix found")
            else:

                print("center " + str(base.centers["base"]))
                print("rotation ")
                print(str(base.rotation_matrix))

                print("Rotated into standard position")
                centered_base = base.translate_rotate_component(base)
                for base_atom in centered_base.atoms():
                    if not "'" in base_atom.name and not "P" in base_atom.name:
                        print("['%s']['%s'] = [%9.6f,%9.6f,%9.6f]" % (base.sequence, base_atom.name, base_atom.coordinates()[0], base_atom.coordinates()[1], base_atom.coordinates()[2]))

                # is the observed base more planar than previous instances?
                maxZ = 0
                for base_atom in centered_base.atoms():
                    if not "'" in base_atom.name and not "P" in base_atom.name:
                        maxZ = max(maxZ,abs(base_atom.coordinates()[2]))

                if maxZ < MinMaxZ:
                    MinMaxZ = maxZ
                    print("This base is more planar than previous ones, maxZ is %0.6f" % maxZ)
                    heavy_atoms = []
                    for base_atom in centered_base.atoms():
                        if not "'" in base_atom.name and not "P" in base_atom.name and not "H" in base_atom.name:
                            print("NAbasecoordinates['%s']['%s'] = [%9.6f,%9.6f,%9.6f]" % (base.sequence, base_atom.name, base_atom.coordinates()[0], base_atom.coordinates()[1], base_atom.coordinates()[2]))
                            heavy_atoms.append([base_atom.coordinates()[0], base_atom.coordinates()[1], base_atom.coordinates()[2]])

                    # special code to center and rotate DT to give an exemplar
                    if base.sequence == "DT":
                        heavy_mean = np.average(heavy_atoms,axis=0)
                        heavy_atoms = heavy_atoms - heavy_mean
                        segment0 = centered_base.centers["N1"][0]-centered_base.centers["C1'"][0]
                        segment1 = centered_base.centers["N1"][1]-centered_base.centers["C1'"][1]
                        theta = np.arctan2(segment0,segment1)
                        c = np.cos(theta)
                        s = np.sin(theta)
                        Rotate = np.array(((c, -s), (s, c)))

                        for base_atom in centered_base.atoms():
                            xypoint = [base_atom.coordinates()[0], base_atom.coordinates()[1]]
                            rotatedpoint = np.dot(Rotate,xypoint)
                            if "H7" in base_atom.name:
                                print("NAbasecoordinates['%s']['%s'] = [%9.6f,%9.6f,%9.6f]" % (base.sequence, base_atom.name, rotatedpoint[0], rotatedpoint[1], base_atom.coordinates()[2]))
                            elif not "'" in base_atom.name:
                                print("NAbasecoordinates['%s']['%s'] = [%9.6f,%9.6f,%9.6f]" % (base.sequence, base_atom.name, rotatedpoint[0], rotatedpoint[1], 0))

                if plot_nts and not base.sequence in bases_plotted:

                    bases_plotted.append(base.sequence)

                    fig = plt.figure()
                    ax = fig.add_subplot(111, projection='3d')

                    PlotGivenBase(centered_base,ax)   # plot the actual base in blue

                    if fields[3] in modified_nucleotides:
                        PlotStandardBase(modified_nucleotides[fields[3]]["standard"],ax)  # plot the standard parent base in red
                    elif fields[3] in NAbasecoordinates:
                        PlotStandardBase(fields[3],ax)  # plot the standard base in red


                    ax.set_xlabel('X Axis')
                    ax.set_ylabel('Y Axis')
                    ax.set_zlabel('Z Axis')
                    ax.set_xlim3d(5,-5)
                    ax.set_ylim3d(5,-5)
                    ax.set_zlim3d(5,-5)
                    ax.view_init(elev=29.0, azim=70.0)
                    plt.title("Base " +str(len(bases_plotted)) + " " + base.unit_id() + " in blue in standard position")
                    plt.show()
                    plt.close()


    print(sorted(bases_plotted))
