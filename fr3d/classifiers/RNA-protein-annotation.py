# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 12:44:30 2014 @author: Poorna Roy
Name: RNA-protein detection
"""

"""Detect and plot RNA base- amino acid interactions."""
from fr3d.cif.reader import Cif
from fr3d.definitions import RNAconnections
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

from discrepancy import matrix_discrepancy
import numpy as np
import csv
import urllib
import pickle
import math
import sys

import matplotlib.pyplot as plt
from collections import defaultdict
from mpl_toolkits.mplot3d import Axes3D
# note that fr3d.localpath does not synchronize with Git, so you can change it locally to point to your own directory structure
from fr3d.localpath import outputText
from fr3d.localpath import outputBaseAAFG
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

HB_donor_hydrogens = {}
HB_donor_hydrogens['A'] = {"N6":["1H6","2H6"], "C2":["H2"], "C8":["H8"], "O2'":[]}
HB_donor_hydrogens['G'] = {"N1":["H1"], "N2":["2H2","1H2"], "C8":["H8"], "O2'":[]}
HB_donor_hydrogens['C'] = {"N4":["1H4","2H4"], "C5":["H5"], "C6":["H6"], "O2'":[]}
HB_donor_hydrogens['U'] = {"N3":["H3"], "C5":["H5"], "C6":["H6"], "O2'":[]}

distance_for_base_atom_type_vdw = {}

distance_for_base_atom_type_vdw['CHn'] = 1.92
distance_for_base_atom_type_vdw['CH'] = 1.82
distance_for_base_atom_type_vdw['C'] = 1.74
distance_for_base_atom_type_vdw['S'] = 1.92
distance_for_base_atom_type_vdw['N'] = 1.66
distance_for_base_atom_type_vdw['O'] = 1.51
distance_for_base_atom_type_vdw['H2O'] = 1.68

distance_for_atom_type_vdw = {}

distance_for_atom_type_vdw['CA'] = 1.9
distance_for_atom_type_vdw['C'] = 1.75
distance_for_atom_type_vdw['CH'] = 2.01
distance_for_atom_type_vdw['CH2'] = 1.92
distance_for_atom_type_vdw['CH2b'] = 1.91
distance_for_atom_type_vdw['CH2ch'] = 1.88
distance_for_atom_type_vdw['CH3'] = 1.92
distance_for_atom_type_vdw['CHar'] = 1.82
distance_for_atom_type_vdw['Car'] = 1.74
distance_for_atom_type_vdw['CHim'] = 1.74
distance_for_atom_type_vdw['Cco'] = 1.81
distance_for_atom_type_vdw['Ccoo'] = 1.76
distance_for_atom_type_vdw['SH'] = 1.88
distance_for_atom_type_vdw['S'] = 1.94
distance_for_atom_type_vdw['N'] = 1.71
distance_for_atom_type_vdw['NH'] = 1.66
distance_for_atom_type_vdw['NH+'] = 1.65
distance_for_atom_type_vdw['NH2'] = 1.62
distance_for_atom_type_vdw['NH2+'] = 1.67
distance_for_atom_type_vdw['NH3+'] = 1.67
distance_for_atom_type_vdw['O'] = 1.49
distance_for_atom_type_vdw['Oco'] = 1.52
distance_for_atom_type_vdw['Ocoo'] = 1.49
distance_for_atom_type_vdw['OH'] = 1.54
distance_for_atom_type_vdw['H2O'] = 1.68

distance_for_base_atom_type_coulombic = {}

distance_for_base_atom_type_coulombic['N'] = 1.47
distance_for_base_atom_type_coulombic['O'] = 1.38
distance_for_base_atom_type_coulombic['H2O'] = 1.37

distance_for_atom_type_coulombic = {}

distance_for_atom_type_coulombic['N'] = 1.49
distance_for_atom_type_coulombic['NH'] = 1.55
distance_for_atom_type_coulombic['NH+'] = 1.41
distance_for_atom_type_coulombic['NH2'] = 1.49
distance_for_atom_type_coulombic['NH2+'] = 1.49
distance_for_atom_type_coulombic['NH3+'] = 1.4
distance_for_atom_type_coulombic['O'] = 1.41
distance_for_atom_type_coulombic['Oco'] = 1.41
distance_for_atom_type_coulombic['Ocoo'] = 1.35
distance_for_atom_type_coulombic['OH'] = 1.35
distance_for_atom_type_coulombic['H2O'] = 1.37

distance_limit_vdw = defaultdict(dict)

for base_name in ["A","C","G","U"]:
    distance_limit_vdw[base_name]["C1'"] = distance_for_base_atom_type_vdw['CHn']
    distance_limit_vdw[base_name]["C2'"] = distance_for_base_atom_type_vdw['CHn']
    distance_limit_vdw[base_name]["C3'"] = distance_for_base_atom_type_vdw['CHn']
    distance_limit_vdw[base_name]["C4'"] = distance_for_base_atom_type_vdw['CHn']
    distance_limit_vdw[base_name]["C5'"] = distance_for_base_atom_type_vdw['CHn']
    distance_limit_vdw[base_name]["O2'"] = distance_for_base_atom_type_vdw['O']
    distance_limit_vdw[base_name]["O3'"] = distance_for_base_atom_type_vdw['O']
    distance_limit_vdw[base_name]["O4'"] = distance_for_base_atom_type_vdw['O']
    distance_limit_vdw[base_name]["O5'"] = distance_for_base_atom_type_vdw['O']
    distance_limit_vdw[base_name]["PO1"] = distance_for_base_atom_type_vdw['O']
    distance_limit_vdw[base_name]["PO2"] = distance_for_base_atom_type_vdw['O']

distance_limit_vdw["A"]["N1"] = distance_for_base_atom_type_vdw['N']
distance_limit_vdw["A"]["C2"] = distance_for_base_atom_type_vdw['CH']
distance_limit_vdw["A"]["N3"] = distance_for_base_atom_type_vdw['N']
distance_limit_vdw["A"]["C4"] = distance_for_base_atom_type_vdw['C']
distance_limit_vdw["A"]["C5"] = distance_for_base_atom_type_vdw['C']
distance_limit_vdw["A"]["C6"] = distance_for_base_atom_type_vdw['C']
distance_limit_vdw["A"]["N6"] = distance_for_base_atom_type_vdw['N']
distance_limit_vdw["A"]["N7"] = distance_for_base_atom_type_vdw['N']
distance_limit_vdw["A"]["C8"] = distance_for_base_atom_type_vdw['CH']
distance_limit_vdw["A"]["N9"] = distance_for_base_atom_type_vdw['N']

distance_limit_vdw["G"]["N1"] = distance_for_base_atom_type_vdw['N']
distance_limit_vdw["G"]["C2"] = distance_for_base_atom_type_vdw['C']
distance_limit_vdw["G"]["N2"] = distance_for_base_atom_type_vdw['N']
distance_limit_vdw["G"]["N3"] = distance_for_base_atom_type_vdw['N']
distance_limit_vdw["G"]["C4"] = distance_for_base_atom_type_vdw['C']
distance_limit_vdw["G"]["C5"] = distance_for_base_atom_type_vdw['C']
distance_limit_vdw["G"]["C6"] = distance_for_base_atom_type_vdw['C']
distance_limit_vdw["G"]["O6"] = distance_for_base_atom_type_vdw['O']
distance_limit_vdw["G"]["N7"] = distance_for_base_atom_type_vdw['N']
distance_limit_vdw["G"]["C8"] = distance_for_base_atom_type_vdw['CH']
distance_limit_vdw["G"]["N9"] = distance_for_base_atom_type_vdw['N']

distance_limit_vdw["C"]["N1"] = distance_for_base_atom_type_vdw['N']
distance_limit_vdw["C"]["C2"] = distance_for_base_atom_type_vdw['C']
distance_limit_vdw["C"]["O2"] = distance_for_base_atom_type_vdw['O']
distance_limit_vdw["C"]["N3"] = distance_for_base_atom_type_vdw['N']
distance_limit_vdw["C"]["C4"] = distance_for_base_atom_type_vdw['C']
distance_limit_vdw["C"]["N4"] = distance_for_base_atom_type_vdw['N']
distance_limit_vdw["C"]["C5"] = distance_for_base_atom_type_vdw['CH']
distance_limit_vdw["C"]["C6"] = distance_for_base_atom_type_vdw['CH']

distance_limit_vdw["U"]["N1"] = distance_for_base_atom_type_vdw['N']
distance_limit_vdw["U"]["C2"] = distance_for_base_atom_type_vdw['C']
distance_limit_vdw["U"]["O2"] = distance_for_base_atom_type_vdw['O']
distance_limit_vdw["U"]["N3"] = distance_for_base_atom_type_vdw['N']
distance_limit_vdw["U"]["C4"] = distance_for_base_atom_type_vdw['C']
distance_limit_vdw["U"]["O4"] = distance_for_base_atom_type_vdw['O']
distance_limit_vdw["U"]["C5"] = distance_for_base_atom_type_vdw['CH']
distance_limit_vdw["U"]["C6"] = distance_for_base_atom_type_vdw['CH']

all_aa_list = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]

for aa_name in all_aa_list:
    distance_limit_vdw[aa_name]['CA'] = distance_for_atom_type_vdw['CA']
    distance_limit_vdw[aa_name]['C'] = distance_for_atom_type_vdw['C']
    distance_limit_vdw[aa_name]['O'] = distance_for_atom_type_vdw['O']
    distance_limit_vdw[aa_name]['N'] = distance_for_atom_type_vdw['N']

distance_limit_vdw['ALA']['CB'] = distance_for_atom_type_vdw['CH3']

distance_limit_vdw['CYS']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['CYS']['SG'] = distance_for_atom_type_vdw['SH']

distance_limit_vdw['SER']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['SER']['OG'] = distance_for_atom_type_vdw['OH']

distance_limit_vdw['THR']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['THR']['CG2'] = distance_for_atom_type_vdw['CH']
distance_limit_vdw['THR']['OG1'] = distance_for_atom_type_vdw['OH']

distance_limit_vdw['ASN']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['ASN']['CG'] = distance_for_atom_type_vdw['Cco']
distance_limit_vdw['ASN']['OD1'] = distance_for_atom_type_vdw['Oco']
distance_limit_vdw['ASN']['ND2'] = distance_for_atom_type_vdw['NH2']

distance_limit_vdw['GLN']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['GLN']['CG'] = distance_for_atom_type_vdw['CH2']
distance_limit_vdw['GLN']['CD'] = distance_for_atom_type_vdw['Cco']
distance_limit_vdw['GLN']['OE1'] = distance_for_atom_type_vdw['Oco']
distance_limit_vdw['GLN']['NE2'] = distance_for_atom_type_vdw['NH2']

distance_limit_vdw['ASP']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['ASP']['CG'] = distance_for_atom_type_vdw['Ccoo']
distance_limit_vdw['ASP']['OD1'] = distance_for_atom_type_vdw['Ocoo']
distance_limit_vdw['ASP']['OD2'] = distance_for_atom_type_vdw['Ocoo']

distance_limit_vdw['GLU']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['GLU']['CG'] = distance_for_atom_type_vdw['CH2ch']
distance_limit_vdw['GLU']['CD'] = distance_for_atom_type_vdw['Ccoo']
distance_limit_vdw['GLU']['OE1'] = distance_for_atom_type_vdw['Ocoo']
distance_limit_vdw['GLU']['OE2'] = distance_for_atom_type_vdw['Ocoo']

distance_limit_vdw['PHE']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['PHE']['CG'] = distance_for_atom_type_vdw['Car']
distance_limit_vdw['PHE']['CD1'] = distance_for_atom_type_vdw['CHar']
distance_limit_vdw['PHE']['CD2'] = distance_for_atom_type_vdw['CHar']
distance_limit_vdw['PHE']['CE1'] = distance_for_atom_type_vdw['CHar']
distance_limit_vdw['PHE']['CE2'] = distance_for_atom_type_vdw['CHar']
distance_limit_vdw['PHE']['CZ'] = distance_for_atom_type_vdw['CHar']

distance_limit_vdw['TYR']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['TYR']['CG'] = distance_for_atom_type_vdw['Car']
distance_limit_vdw['TYR']['CD1'] = distance_for_atom_type_vdw['CHar']
distance_limit_vdw['TYR']['CD2'] = distance_for_atom_type_vdw['CHar']
distance_limit_vdw['TYR']['CE1'] = distance_for_atom_type_vdw['CHar']
distance_limit_vdw['TYR']['CE2'] = distance_for_atom_type_vdw['CHar']
distance_limit_vdw['TYR']['CZ'] = distance_for_atom_type_vdw['Car']
distance_limit_vdw['TYR']['OH'] = distance_for_atom_type_vdw['OH']

distance_limit_vdw['HIS']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['HIS']['CG'] = distance_for_atom_type_vdw['Car']
distance_limit_vdw['HIS']['ND1'] = distance_for_atom_type_vdw['NH+']
distance_limit_vdw['HIS']['CD2'] = distance_for_atom_type_vdw['CHim']
distance_limit_vdw['HIS']['CE1'] = distance_for_atom_type_vdw['CHim']
distance_limit_vdw['HIS']['NE2'] = distance_for_atom_type_vdw['NH+']

distance_limit_vdw['LEU']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['LEU']['CG'] = distance_for_atom_type_vdw['CH']
distance_limit_vdw['LEU']['CD1'] = distance_for_atom_type_vdw['CH3']
distance_limit_vdw['LEU']['CD2'] = distance_for_atom_type_vdw['CH3']

distance_limit_vdw['ILE']['CB'] = distance_for_atom_type_vdw['CH']
distance_limit_vdw['ILE']['CG1'] = distance_for_atom_type_vdw['CH2']
distance_limit_vdw['ILE']['CG2'] = distance_for_atom_type_vdw['CH3']
distance_limit_vdw['ILE']['CD'] = distance_for_atom_type_vdw['CH3']

distance_limit_vdw['VAL']['CB'] = distance_for_atom_type_vdw['CH']
distance_limit_vdw['VAL']['CG1'] = distance_for_atom_type_vdw['CH3']
distance_limit_vdw['VAL']['CG2'] = distance_for_atom_type_vdw['CH3']

distance_limit_vdw['MET']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['MET']['CG'] = distance_for_atom_type_vdw['CH2']
distance_limit_vdw['MET']['SD'] = distance_for_atom_type_vdw['S']
distance_limit_vdw['MET']['CE'] = distance_for_atom_type_vdw['CH3']

distance_limit_vdw['PRO']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['PRO']['CG'] = distance_for_atom_type_vdw['CH2']
distance_limit_vdw['PRO']['CD'] = distance_for_atom_type_vdw['CH2']

distance_limit_vdw['LYS']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['LYS']['CG'] = distance_for_atom_type_vdw['CH2']
distance_limit_vdw['LYS']['CD'] = distance_for_atom_type_vdw['CH2']
distance_limit_vdw['LYS']['CE'] = distance_for_atom_type_vdw['CH2ch']
distance_limit_vdw['LYS']['NZ'] = distance_for_atom_type_vdw['NH3+']

distance_limit_vdw['ARG']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['ARG']['CG'] = distance_for_atom_type_vdw['CH2']
distance_limit_vdw['ARG']['CD'] = distance_for_atom_type_vdw['CH2ch']
distance_limit_vdw['ARG']['NE'] = distance_for_atom_type_vdw['NH2+']
distance_limit_vdw['ARG']['CZ'] = distance_for_atom_type_vdw['Car']
distance_limit_vdw['ARG']['NH1'] = distance_for_atom_type_vdw['NH2+']
distance_limit_vdw['ARG']['NH2'] = distance_for_atom_type_vdw['NH2+']

distance_limit_vdw['TRP']['CB'] = distance_for_atom_type_vdw['CH2b']
distance_limit_vdw['TRP']['CG'] = distance_for_atom_type_vdw['Car']
distance_limit_vdw['TRP']['CD1'] = distance_for_atom_type_vdw['CHar']
distance_limit_vdw['TRP']['CD2'] = distance_for_atom_type_vdw['Car']
distance_limit_vdw['TRP']['NE1'] = distance_for_atom_type_vdw['NH']

distance_limit_vdw['TRP']['CE2'] = distance_for_atom_type_vdw['Car']
distance_limit_vdw['TRP']['CZ2'] = distance_for_atom_type_vdw['CHar']
distance_limit_vdw['TRP']['CH2'] = distance_for_atom_type_vdw['CHar']
distance_limit_vdw['TRP']['CE3'] = distance_for_atom_type_vdw['CHar']
distance_limit_vdw['TRP']['CZ3'] = distance_for_atom_type_vdw['CHar']

distance_limit_vdw['GLY']['CA'] = distance_for_atom_type_vdw['CH2b']

distance_limit_coulombic = defaultdict(dict)

for base_name in ["A","C","G","U"]:
    distance_limit_coulombic[base_name]["O2'"] = distance_for_base_atom_type_coulombic['O']
    distance_limit_coulombic[base_name]["O3'"] = distance_for_base_atom_type_coulombic['O']
    distance_limit_coulombic[base_name]["O4'"] = distance_for_base_atom_type_coulombic['O']
    distance_limit_coulombic[base_name]["O5'"] = distance_for_base_atom_type_coulombic['O']
    distance_limit_coulombic[base_name]["PO1"] = distance_for_base_atom_type_coulombic['O']
    distance_limit_coulombic[base_name]["PO2"] = distance_for_base_atom_type_coulombic['O']

distance_limit_coulombic['A']['N1'] = distance_for_base_atom_type_coulombic['N']
distance_limit_coulombic['A']['N3'] = distance_for_base_atom_type_coulombic['N']
distance_limit_coulombic['A']['N6'] = distance_for_base_atom_type_coulombic['N']
distance_limit_coulombic['A']['N7'] = distance_for_base_atom_type_coulombic['N']
distance_limit_coulombic['A']['N9'] = distance_for_base_atom_type_coulombic['N']

distance_limit_coulombic['G']['N1'] = distance_for_base_atom_type_coulombic['N']
distance_limit_coulombic['G']['N2'] = distance_for_base_atom_type_coulombic['N']
distance_limit_coulombic['G']['N3'] = distance_for_base_atom_type_coulombic['N']
distance_limit_coulombic['G']['O6'] = distance_for_base_atom_type_coulombic['O']
distance_limit_coulombic['G']['N7'] = distance_for_base_atom_type_coulombic['N']
distance_limit_coulombic['G']['N9'] = distance_for_base_atom_type_coulombic['N']

distance_limit_coulombic['C']['N1'] = distance_for_base_atom_type_coulombic['N']
distance_limit_coulombic['C']['O2'] = distance_for_base_atom_type_coulombic['O']
distance_limit_coulombic['C']['N3'] = distance_for_base_atom_type_coulombic['N']
distance_limit_coulombic['C']['N4'] = distance_for_base_atom_type_coulombic['N']

distance_limit_coulombic['U']['N1'] = distance_for_base_atom_type_coulombic['N']
distance_limit_coulombic['U']['O2'] = distance_for_base_atom_type_coulombic['O']
distance_limit_coulombic['U']['N3'] = distance_for_base_atom_type_coulombic['N']
distance_limit_coulombic['U']['O4'] = distance_for_base_atom_type_coulombic['O']


for aa_name in all_aa_list:
    distance_limit_coulombic[aa_name]['O'] = distance_for_atom_type_coulombic['O']
    distance_limit_coulombic[aa_name]['N'] = distance_for_atom_type_coulombic['N']

distance_limit_coulombic['SER']['OG'] = distance_for_atom_type_coulombic['OH']

distance_limit_coulombic['THR']['OG1'] = distance_for_atom_type_coulombic['OH']

distance_limit_coulombic['ASN']['OD1'] = distance_for_atom_type_coulombic['Oco']
distance_limit_coulombic['ASN']['ND2'] = distance_for_atom_type_coulombic['NH2']

distance_limit_coulombic['GLN']['OE1'] = distance_for_atom_type_coulombic['Oco']
distance_limit_coulombic['GLN']['NE2'] = distance_for_atom_type_coulombic['NH2']

distance_limit_coulombic['ASP']['OD1'] = distance_for_atom_type_coulombic['Ocoo']
distance_limit_coulombic['ASP']['OD2'] = distance_for_atom_type_coulombic['Ocoo']

distance_limit_coulombic['GLU']['OE1'] = distance_for_atom_type_coulombic['Ocoo']
distance_limit_coulombic['GLU']['OE2'] = distance_for_atom_type_coulombic['Ocoo']

distance_limit_coulombic['TYR']['OH'] = distance_for_atom_type_coulombic['OH']

distance_limit_coulombic['HIS']['ND1'] = distance_for_atom_type_coulombic['NH+']
distance_limit_coulombic['HIS']['NE2'] = distance_for_atom_type_coulombic['NH+']

distance_limit_coulombic['LYS']['NZ'] = distance_for_atom_type_coulombic['NH3+']

distance_limit_coulombic['ARG']['NE'] = distance_for_atom_type_coulombic['NH2+']
distance_limit_coulombic['ARG']['NH1'] = distance_for_atom_type_coulombic['NH2+']
distance_limit_coulombic['ARG']['NH2'] = distance_for_atom_type_coulombic['NH2+']

distance_limit_coulombic['TRP']['NE1'] = distance_for_atom_type_coulombic['NH']

def myTimer(state,data={}):

    # add elapsed time to the current state of the timer
    if "currentState" in data:
        currentState = data["currentState"]
        data[currentState] += time() - data["lastTime"]

    if state == "summary":
        print("Summary of time taken:")
        for state in data["allStates"]:
            if not state == "lastTime" and not state == "currentState":
                print("%-31s: %10.3f seconds %10.3f minutes" % (state,data[state],data[state]/60))
    elif not state in data:
        data[state] = 0
        # keep track of states and the order in which they were seen
        if "allStates" in data:
            data["allStates"].append(state)
        else:
            data["allStates"] = [state]

    # change to the state just starting now
    data["currentState"] = state
    data["lastTime"] = time()

    return data

def get_structure(filename):

    if ".pdb" in filename:
        filename = filename.replace(".cif","")

    if os.path.exists(filename+".pickle"):
        print("  Loading " + filename + ".pickle")
        structure = pickle.load(open(filename+".pickle","rb"))
        return structure
    print(filename)
    if not os.path.exists(filename):
        mmCIFname = filename[-8:]               # last 8 characters ... awkward
        print("  Downloading "+mmCIFname)
        print("https://files.rcsb.org/download/%s" % mmCIFname)
        if sys.version_info[0] < 3:
            urllib.urlretrieve("http://files.rcsb.org/download/%s" % mmCIFname, filename)  # python 2
        else:
            urllib.request.urlretrieve("http://files.rcsb.org/download/%s" % mmCIFname, filename)  # python 3
        #f = urllib.urlopen("https://files.rcsb.org/download/%s" % mmCIFname)
        #myfile = f.read()
        #print(myfile[0:1000])
        #with open(filename, 'w') as outfile:
            #outfile.write(myfile)

    # uncomment the following line to focus on downloading CIF files; sometimes it hangs and you restart
#    raise Exception("Skipping CIF reading for now")

    with open(filename, 'rb') as raw:
        print("  Loading " + filename)
        structure = Cif(raw).structure()
        """All RNA bases are placed in the standard orientation.
        Rotation matrix is calculated for each base."""

        structure.infer_hydrogens()  # add hydrogens to RNA bases and amino acids

#        pickle.dump(structure,open(filename+".pickle","wb"))  # larger file sizes than .cif ... not sure why

        return structure

def build_atom_to_unit_part_list():

    atom_to_part_list = defaultdict(lambda: "unknown")

    for base in NAbaseheavyatoms.keys():
        for atom in NAbaseheavyatoms[base]:
            atom_to_part_list[(base,atom)] = "base"
        for atom in NAbasehydrogens[base]:
            atom_to_part_list[(base,atom)] = "base"
        for atom in nt_phosphate[base]:
            atom_to_part_list[(base,atom)] = "nt_phosphate"
        for atom in nt_sugar[base]:
            atom_to_part_list[(base,atom)] = "nt_sugar"

    for aa in aa_backbone.keys():
        for atom in aa_backbone[aa]:
            atom_to_part_list[(aa,atom)] = "aa_backbone"
        for atom in aa_linker[aa]:
            atom_to_part_list[(aa,atom)] = "aa_linker"
        for atom in aa_fg[aa]:
            atom_to_part_list[(aa,atom)] = "aa_fg"

#    print(atom_to_part_list)

    return atom_to_part_list


def find_atom_atom_contacts(bases,amino_acids,atom_atom_min_distance):
    """Find all atoms within atom_atom_min_distance of each other,
    one from a nt, one from an aa, and report these"""
    # build a set of cubes and record which bases are in which cube
    # also record which other cubes are neighbors of each cube

    atom_atom_min_distance_squared = atom_atom_min_distance**2
    contact_list = []

    atom_to_part_list = build_atom_to_unit_part_list()

    # build a set of cubes and record which nt atoms are in which cube
    ntAtomCubeList = {}
    ntAtomCubeNeighbors = {}
    for base in bases:
        for atom in base.atoms():
            center = atom.coordinates()
            if len(center) == 3:
                x = floor(center[0]/atom_atom_min_distance)
                y = floor(center[1]/atom_atom_min_distance)
                z = floor(center[2]/atom_atom_min_distance)
                key = "%d,%d,%d" % (x,y,z)
                entry = (center[0],center[1],center[2],atom.name,base.unit_id(),atom_to_part_list[(base.sequence,atom.name)])
                if key in ntAtomCubeList:
                    ntAtomCubeList[key].append(entry)
                else:
                    ntAtomCubeList[key] = [entry]
                    ntAtomCubeNeighbors[key] = []
                    for a in [-1,0,1]:
                        for b in [-1,0,1]:
                            for c in [-1,0,1]:
                                k = "%d,%d,%d" % (x+a,y+b,z+c)
                                ntAtomCubeNeighbors[key].append(k)

    # build a set of cubes and record which amino acids are in which cube
    aaAtomCubeList = {}
    for aa in amino_acids:
        for atom in aa.atoms():
            center = atom.coordinates()
            if len(center) == 3:
                x = floor(center[0]/atom_atom_min_distance)
                y = floor(center[1]/atom_atom_min_distance)
                z = floor(center[2]/atom_atom_min_distance)
                key = "%d,%d,%d" % (x,y,z)
                entry = (center[0],center[1],center[2],atom.name,aa.unit_id(),atom_to_part_list[(aa.sequence,atom.name)])
                if key in aaAtomCubeList:
                    aaAtomCubeList[key].append(entry)
                else:
                    aaAtomCubeList[key] = [entry]

    # loop over nt atoms and aa atoms and find those within atom_atom_min_distance
    for key in ntAtomCubeList:                           # one nt atom cube
        for aakey in ntAtomCubeNeighbors[key]:           # neighboring cube
            if aakey in aaAtomCubeList:                  # if an aa atom lies in the neighbor cube
                for ntAtom in ntAtomCubeList[key]:       # loop over nt atoms in first cube
                    for aaAtom in aaAtomCubeList[aakey]: # and over aa atoms in neighbor cube
                        x = abs(ntAtom[0] - aaAtom[0])   # coordinates are stored in 0, 1, 2 of entry
                        y = abs(ntAtom[1] - aaAtom[1])   # coordinates are stored in 0, 1, 2 of entry
                        z = abs(ntAtom[2] - aaAtom[2])   # coordinates are stored in 0, 1, 2 of entry

                        if x > atom_atom_min_distance or \
                           y > atom_atom_min_distance or \
                           x*x + y*y + z*z > atom_atom_min_distance_squared:
                            continue

                        distance = math.sqrt(x*x + y*y + z*z)

                        contact_list.append("%s\t%s\t%s\t%s\t%s\t%s\t%8.4f\n" % (ntAtom[4],ntAtom[5],ntAtom[3],aaAtom[4],aaAtom[5],aaAtom[3],distance))

    return contact_list

def find_neighbors(bases, amino_acids, screen_distance_cutoff, IFE, nt_reference="C1'", aa_reference="aa_fg"):
    """Finds all amino acids for which center of "aa_part" is within
    specified distance of center of bases
    For annotating files in a representative set, it also screens for the
    base being in the given IFE. """

    # build a set of cubes and record which bases are in which cube
    # also record which other cubes are neighbors of each cube
    baseCubeList = {}
    baseCubeNeighbors = {}

    # build a set of cubes and record which bases are in which cube
    for base in bases:

        center = base.centers[nt_reference]               # a reasonably central atom
        if len(center) == 3:
            x = floor(center[0]/screen_distance_cutoff)
            y = floor(center[1]/screen_distance_cutoff)
            z = floor(center[2]/screen_distance_cutoff)
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
        center = aa.centers[aa_reference]     # as of 9/1/2019, some aa lack C atom because of alt ids A or B
        if len(center) == 3:
            x = floor(center[0]/screen_distance_cutoff)
            y = floor(center[1]/screen_distance_cutoff)
            z = floor(center[2]/screen_distance_cutoff)
            key = "%d,%d,%d" % (x,y,z)
            if key in aaCubeList:
                aaCubeList[key].append(aa)
            else:
                aaCubeList[key] = [aa]
        else:
            print("  Missing center coordinates for " + str(aa))

    return baseCubeList, baseCubeNeighbors, aaCubeList

# produce a list of interacting atoms marked by their group
def get_interacting_atoms(nt_residue,aa_residue):

    interacting_atoms = defaultdict(list)

    atom_atom_min_distance_squared = atom_atom_min_distance**2
    nt_parts = ['base','nt_sugar','nt_phosphate']
    aa_parts = ['aa_fg','aa_backbone','aa_linker']
    nt_parts = ['base']
    aa_parts = ['aa_fg']
    nt_atoms = []
    aa_atoms = []

    for nt_part in nt_parts:
        if nt_part == "base" and nt_residue.sequence in NAbaseheavyatoms:
            # important: start with an empty list and add this list,
            # rather than making nt_atoms equal to NAbaseheavyatoms[...]
            # and then appending, because that will append to the list
            # in the dictionary, with unpredictable results.
            nt_atoms += NAbaseheavyatoms[nt_residue.sequence]
            if nt_residue.sequence in NAbasehydrogens:
                nt_atoms += NAbasehydrogens[nt_residue.sequence]
        elif nt_part == "nt_sugar" and nt_residue.sequence in nt_sugar:
            nt_atoms += nt_sugar[nt_residue.sequence]
        elif nt_residue.sequence in nt_phosphate:
            nt_atoms += nt_phosphate[nt_residue.sequence]
        else:
            continue

        for aa_part in aa_parts:
            if aa_part == "aa_fg" and aa_residue.sequence in aa_fg:
                aa_atoms += aa_fg[aa_residue.sequence]
            elif aa_part == "aa_backbone" and aa_residue.sequence in aa_backbone:
                aa_atoms += aa_backbone[aa_residue.sequence]
            elif aa_residue.sequence in aa_linker:
                aa_atoms += aa_linker[aa_residue.sequence]
            else:
                continue

            for nt_atom in nt_residue.atoms(name=nt_atoms):
                nt = nt_atom.coordinates()
                for aa_atom in aa_residue.atoms(name=aa_atoms):
                    aa = aa_atom.coordinates()
                    x = abs(nt[0]-aa[0])
                    y = abs(nt[1]-aa[1])
                    z = abs(nt[2]-aa[2])

                    if x <= atom_atom_min_distance and \
                       y <= atom_atom_min_distance and \
                       x**2 + y**2 + z**2 <= atom_atom_min_distance_squared:
                        distance = math.sqrt(x**2 + y**2 + z**2)
                        interacting_atoms[(nt_part,aa_part)].append((nt_atom,aa_atom,distance))

#    if len(interacting_atoms) > 0:
#        print(interacting_atoms)

    return interacting_atoms

def annotate_interactions(bases, amino_acids, screen_distance_cutoff, baseCubeList, baseCubeNeighbors, aaCubeList):

    # loop through base cubes, loop through neighboring amino acid cubes,
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
    hbond_aa_dict = {}
    contact_list = []
    output = []

    max_screen_distance = 0

    for key in baseCubeList:
        for aakey in baseCubeNeighbors[key]:
            if aakey in aaCubeList:
                for base_residue in baseCubeList[key]:
                    if len(base_residue.centers["base"]) < 3:
                        print("Missing base center for %s" % base_residue.unit_id())
                        print(base_residue.centers["base"])
                        continue
                    for aa_residue in aaCubeList[aakey]:
                        if len(aa_residue.centers["aa_fg"]) < 3:
                            print("Missing aa_fg center for %s" % aa_residue.unit_id())
                            continue
                        displacement = abs(base_residue.centers["base"]-aa_residue.centers["aa_fg"])
                        #print(("  Missing center atom for ") + base_residue.unit_id() + " " + aa_residue.unit_id())

                        if displacement[0] > screen_distance_cutoff or \
                           displacement[1] > screen_distance_cutoff or \
                           np.linalg.norm(displacement) > screen_distance_cutoff:
                            continue

                        screen_distance = np.linalg.norm(displacement)
                        interacting_atoms = get_interacting_atoms(base_residue,aa_residue)

                        if len(interacting_atoms) > 0:
                            if screen_distance > max_screen_distance:
                                max_screen_distance = screen_distance

                        # annotate interactions between base and functional group in detail
                        if ("base","aa_fg") in interacting_atoms:

                            if aa_residue.sequence in set(['LYS', 'SER', 'THR', 'TYR']):
                                if len(interacting_atoms[("base","aa_fg")]) < 1:
                                    continue
                            else:
                                if len(interacting_atoms[("base","aa_fg")]) < 2:
                                    continue

                            base_center = base_residue.centers["base"]

                            if not base_center.any():
                                continue

                            base_seq = base_residue.sequence
                            base_atoms = NAbaseheavyatoms[base_seq]

                            # check sets of three atoms that lie in the plane of the base, for normal calculations
                            if not base_residue.centers[planar_atoms[base_seq][0]].any():
                                continue
                            if not base_residue.centers[planar_atoms[base_seq][1]].any():
                                continue
                            if not base_residue.centers[planar_atoms[base_seq][2]].any():
                                continue

                            aa_center = aa_residue.centers[aa_part]
                            if not aa_center.any():
                                continue

                            count_pair = count_pair + 1

                            rotation_matrix = base_residue.rotation_matrix

                            # note:  translate_rotate_component is in components.py and calls infer_hydrogens

                            # rotate base atoms into standard orientation
                            base_coordinates = {}
                            standard_base = base_residue.translate_rotate_component(base_residue)

                            for base_atom in standard_base.atoms():
                                base_coordinates[base_atom.name] = base_atom.coordinates()
                                # heavy atoms are close to ideal, hydrogens are almost exactly right
                                # print("Standard orientation " + base_atom.name + " coordinates")
                                # print(base_atom.coordinates())

                            # rotate amino acid atoms into standard orientation
                            aa_coordinates = {}
                            standard_aa = base_residue.translate_rotate_component(aa_residue)
                            for aa_atom in standard_aa.atoms():
                                aa_coordinates[aa_atom.name] = aa_atom.coordinates()

                            standard_aa_center = standard_aa.centers[aa_part]

                            # get a preliminary annotation of the interaction
#                            (interaction,interaction_parameters) = type_of_interaction(base_residue, aa_residue, aa_coordinates, standard_aa_center, base_atoms)
                            (interaction,interaction_parameters) = type_of_interaction(standard_base, standard_aa, aa_coordinates, standard_aa_center, base_atoms)

                            base_aa = None
                            edge = None
                            face = None
                            if interaction in ["pseudopair","SHB","perpendicular-edge","other-edge"]:
                                (edge,angle) = detect_base_edge(base_residue, base_coordinates,aa_residue, aa_coordinates)
                                interaction_parameters["angle-in-plane"] = angle
                                base_aa = (base_residue, aa_residue, interaction, edge, standard_aa, interaction_parameters)
                                (face,height) = detect_face(aa_residue, aa_coordinates)

                            elif interaction in ["stacked","pi-pi-stacking","cation-pi","perpendicular-stacking","other-stack"]:
                                (face,height) = detect_face(aa_residue, aa_coordinates)
                                base_aa = (base_residue, aa_residue, interaction, face, standard_aa, interaction_parameters)
                                (edge,angle) = detect_base_edge(base_residue, base_coordinates,aa_residue, aa_coordinates)

                            else:
                                (face,height) = detect_face(aa_residue, aa_coordinates)
                                base_aa = (base_residue, aa_residue, interaction, face, standard_aa, interaction_parameters)
                                (edge,angle) = detect_base_edge(base_residue, base_coordinates,aa_residue, aa_coordinates)

                            if base_aa is not None:
                                list_base_aa.append(base_aa)

                                for base_atom in base_residue.atoms():
                                    list_base_coord.append(base_coordinates)
                                for aa_atom in aa_residue.atoms():
                                    list_aa_coord.append(aa_coordinates)

                            # detect one amino acid interacting with two bases
                            if interaction in ["pseudopair","SHB"]:
                                if aa_residue.unit_id() in hbond_aa_dict:
                                    hbond_aa_dict[aa_residue.unit_id()].append((base_residue, aa_residue, interaction, edge, standard_aa, interaction_parameters))
#                                    print("  %s makes hbonds with %d bases" % (aa_residue.unit_id(),len(hbond_aa_dict[aa_residue.unit_id()])))
                                else:
                                    hbond_aa_dict[aa_residue.unit_id()] = [(base_residue, aa_residue, interaction, edge, standard_aa, interaction_parameters)]

                            #try:
                                #output.append((base_residue.unit_id(),aa_residue.unit_id(),base_residue.sequence,aa_residue.sequence,standard_aa_center[0],standard_aa_center[1],standard_aa_center[2],aa_coordinates['CA'][0],aa_coordinates['CA'][1],aa_coordinates['CA'][2],interaction,edge,face))
                            #except:
                                #print("Missing CA")


    save_path = '/Users/katelandsipe/Documents/Research/FR3D/nt_aa_interactions'
    #output_file = "nt_aa_coordinates_3A"
    #protein_aa_interactions = os.path.join(save_path, output_file+".csv")
    #file = open(protein_aa_interactions, 'a+')

    #writing the data into the file
    #with file:
        #write = csv.writer(file)
        #write.writerows(output)

    #file.close()

    print("  Found %d nucleotide-amino acid pairs" % count_pair)
    print("  Recorded %d nucleotide-amino acid pairs" % len(list_base_aa))
    print("  Maximum screen distance for actual contacts is %8.4f" % max_screen_distance)

    return list_base_aa, list_aa_coord, list_base_coord, hbond_aa_dict

def type_of_interaction(base_residue, aa_residue, aa_coordinates, standard_aa_center, base_atoms):
    """ This function works with base and aa in standard position """
    squared_xy_dist_list = []

    """Defines different sets of amino acids"""
    planar_aa = set (["ARG", "ASN", "ASP", "GLU", "GLN", "HIS", "PHE", "TRP", "TYR"])
    stacked_aliphatic = set(["ALA", "CYS", "ILE", "LEU", "MET", "PRO", "SER", "THR", "VAL"])
    # Note:  LYS and GLY are not in the previous lists

    edge_to_edge_aa = set (["ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "HIS", "LYS", "PHE", "SER", "THR", "TYR", "TRP"])
    shb_aa = set (["ARG", "ASN", "ASP", "GLU", "GLN", "HIS", "LYS", "SER", "THR", "TYR"])

    # calculate distances from aa atoms to base center
    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key  = aa_atom.name
        aa_x = aa_coordinates[key][0]
        aa_y = aa_coordinates[key][1]

        squared_xy_dist = (aa_x**2) + (aa_y**2)
        squared_xy_dist_list.append(squared_xy_dist)

    #print base_residue.unit_id(), aa_residue.unit_id(), min(squared_xy_dist_list), mean_z

    # for a stacking interaction, the x,y coordinate of at least one atom of the amino acid group needs to be
    # within sqrt(5) = 2.236 of the base center 0,0
    min_dist = np.sqrt(min(squared_xy_dist_list))

    if min_dist <= 2.236:
        #print base_residue.unit_id(), aa_residue.unit_id(), min(squared_xy_dist_list), mean_z
        if aa_residue.sequence in planar_aa:
            return stacking_planar_annotation(base_residue, aa_residue, min_dist)

        elif aa_residue.sequence in stacked_aliphatic:
            return stacking_non_planar_annotation(aa_residue, aa_coordinates, min_dist)

        else:
            return ("other-stack",{"dist-xy-from-center":min_dist})

    # check for interactions in the plane of the base
    mean_z = standard_aa_center[2]
    (num_hydrogen_bonds,hydrogen_bond_list) = count_hydrogen_bonds(base_residue, aa_residue, base_atoms)

    if -1.8 <= mean_z < 1.8:
        if aa_residue.sequence in edge_to_edge_aa:
            angle = calculate_angle_between_planar_residues(base_residue, aa_residue)
            if angle:
                if 0 <= angle <= 45 and num_hydrogen_bonds >= 2:
                    return ("pseudopair",{"hydrogen-bonds":hydrogen_bond_list,"angle-between-planes":angle})
                elif 45 <= angle:
                    return ("perpendicular-edge",{"hydrogen-bonds":hydrogen_bond_list,"angle-between-planes":angle})

        if aa_residue.sequence in shb_aa:
#            base_seq = base_residue.sequence
#            base_atoms = NAbaseheavyatoms[base_seq]
            if num_hydrogen_bonds >= 1:
                return ("SHB",{"hydrogen-bonds":hydrogen_bond_list})

        return ("other-edge",{"hydrogen-bonds":hydrogen_bond_list})

    return ("other",{"dist-xy-from-center":min_dist,"hydrogen-bonds":hydrogen_bond_list})

def count_hydrogen_bonds(base_residue, aa_residue, base_atoms):
    """Calculates number of Hydrogen bonds between amino acid part and base_part
    and returns the number and a list of pairs of atoms from (base,aa)
    """

    n = 0                            # number of independent hydrogen bonds being counted toward pseudopairs
    hydrogen_bond_list = []

    aa_key = aa_residue.sequence

    if not (aa_key in HB_donors or aa_key in HB_acceptors or aa_key in HB_weak_donors):
        return (n,hydrogen_bond_list)

    n0 = 0
    hb0 = []

    hb_angle_cutoff = 100
    screen_distance = 4.3                 # around 85% are higher than 4; 4 is a good first screen

    used_base_atoms = []
    used_aa_atoms = []
    carboxylate_donor_used = False   # track the two oxygens on ASP and GLU; only one can be a donor
    HIS_acceptor_used = False        # track the two nitrogens on HIS; only one can be an acceptor

    base_key = base_residue.sequence
    base_donors = HB_donor_hydrogens[base_key].keys()
    base_acceptors = HB_acceptors[base_key]
    base_HB_atoms = list(set(base_donors + base_acceptors))  # don't list atoms twice

    # these amino acids can be mis-modeled with the functional group flipped 180 degrees
    if aa_key in ["ASN","GLN"]:
        num_flip_states = 2
    else:
        num_flip_states = 1

    # check some amino acids in original and flipped orientations
    for flip in range(0,num_flip_states):
        n = 0
        hydrogen_bond_list = []
        if flip == 0:
            aa_donors = HB_donors[aa_key]
            if aa_key in HB_weak_donors:
                base_donors += HB_weak_donors[aa_key]
            aa_acceptors = HB_acceptors[aa_key]
            flip_name = ""
        else:
            # pretend that these three amino acids are flipped; they are often mis-modeled
            if aa_key == "ASN":
                aa_donors = ["OD1"]
                aa_acceptors = ["ND2"]
            elif aa_key == "GLN":
                aa_donors = ["OE1"]
                aa_acceptors = ["NE2"]
            # for now, don't try to recognize flipped HIS, just check for donors/acceptors
#                elif aa_key == "HIS":
#                    aa_donors = ['CD2', 'CE1']
#                    aa_acceptors = ['ND1', 'NE2']
            flip_name = "f"

        for base_atom in base_residue.atoms(name=base_HB_atoms):
            for aa_atom in aa_residue.atoms(name=aa_fg[aa_key]):
                difference = np.subtract(base_atom.coordinates(), aa_atom.coordinates())
                distance = np.linalg.norm(difference)

                if distance > screen_distance:
                    continue

                # check the distance between heavy atoms, depending on the two interacting atoms and the residues
                if ("O" in base_atom.name or "N" in base_atom.name) and ("O" in aa_atom.name or "N" in aa_atom.name):
                    h_bond_ideal_distance  = distance_limit_coulombic[base_residue.sequence][base_atom.name]
                    h_bond_ideal_distance += distance_limit_coulombic[aa_residue.sequence][aa_atom.name]
                    h_bond_over_distance = distance - h_bond_ideal_distance

                    hbondtext = base_residue.sequence +"\t"+ base_atom.name +"\t"+ str(distance_limit_coulombic[base_residue.sequence][base_atom.name]) +"\t"
                    hbondtext += aa_residue.sequence +"\t"+ aa_atom.name +"\t"+ str(distance_limit_coulombic[aa_residue.sequence][aa_atom.name])
                    hbondtext += "\t"+ str(h_bond_ideal_distance) +"\t"+ str(distance) +"\t"+ str(h_bond_over_distance)

                else:
                    h_bond_ideal_distance  = distance_limit_vdw[base_residue.sequence][base_atom.name]
                    h_bond_ideal_distance += distance_limit_vdw[aa_residue.sequence][aa_atom.name]
                    h_bond_over_distance = distance - h_bond_ideal_distance

                    hbondtext = base_residue.sequence +"\t"+ base_atom.name +"\t"+ str(distance_limit_vdw[base_residue.sequence][base_atom.name]) +"\t"
                    hbondtext += aa_residue.sequence +"\t"+ aa_atom.name +"\t"+ str(distance_limit_vdw[aa_residue.sequence][aa_atom.name])
                    hbondtext += "\t"+ str(h_bond_ideal_distance) +"\t"+ str(distance) +"\t"+ str(h_bond_over_distance)

#                print("hbond\t"+hbondtext)

                if distance > h_bond_ideal_distance + 0.4:
                    continue
                    h_bond_over_distance = distance - h_bond_ideal_distance

                if base_atom.name in base_donors and aa_atom.name in aa_acceptors:
                    # loop through hydrogens whose locations are known
                    for hydrogen_atom in base_residue.atoms(name=HB_donor_hydrogens[base_key][base_atom.name]):
                        hb_angle = calculate_hb_angle(base_atom.coordinates(),hydrogen_atom.coordinates(),aa_atom.coordinates())
#                            print("hb_angle %10.8f" % hb_angle)
                        if hb_angle > hb_angle_cutoff:
                            if not base_atom.name in used_base_atoms and not aa_atom.name in used_aa_atoms:
                                if aa_key == "HIS":
                                    if not HIS_acceptor_used:
                                        n = n + 1
                                    HIS_acceptor_used = True
                                else:
                                    n = n + 1
                            # print(hydrogen_atom.name, hydrogen_atom.coordinates(), hb_angle)
                            hydrogen_bond_list.append((base_atom.name,hydrogen_atom.name,aa_atom.name+flip_name,h_bond_ideal_distance,distance,hb_angle))
                            used_base_atoms.append(base_atom.name)
                            used_aa_atoms.append(aa_atom.name)

                    # for O2', coordinates of H are not known, so check angles separately
                    if base_atom.name == "O2'":
                        hb_angle = calculate_hb_angle(base_residue.centers["C2'"],base_atom.coordinates(),aa_atom.coordinates())
#                        print("C2'-O2'-A angle %10.8f" % hb_angle)
                        if hb_angle > 80:
                            if not base_atom.name in used_base_atoms and not aa_atom.name in used_aa_atoms:
                                if aa_key == "HIS":
                                    if not HIS_acceptor_used:
                                        n = n + 1
                                    HIS_acceptor_used = True
                                else:
                                    n = n + 1
                            hydrogen_bond_list.append((base_atom.name,"O2'H",aa_atom.name+flip_name,h_bond_ideal_distance,distance,hb_angle))
                            used_base_atoms.append(base_atom.name)
                            used_aa_atoms.append(aa_atom.name)


                elif base_atom.name in base_acceptors and aa_atom.name in aa_donors:
                    # check OH group on certain amino acids; coordinates of H are not known
                    if (aa_key == "SER" and aa_atom.name == "OG") or \
                       (aa_key == "THR" and aa_atom.name == "OG1") or \
                       (aa_key == "TYR" and aa_atom.name == "OH"):

                        if aa_key == "TYR":
                            carbon = "CZ"
                        else:
                            carbon = "CB"

                        hb_angle = calculate_hb_angle(aa_residue.centers[carbon],aa_atom.coordinates(),base_atom.coordinates())
#                        if hb_angle:
#                            print("%s-OH-%s angle %10.8f ------------ http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (carbon,base_atom.name,hb_angle,base_residue.unit_id(),aa_residue.unit_id()))

                        if hb_angle > 80:
                            if not base_atom.name in used_base_atoms and not aa_atom.name in used_aa_atoms:
                                n = n + 1
                            hydrogen_bond_list.append((base_atom.name,"OHH",aa_atom.name+flip_name,h_bond_ideal_distance,distance,hb_angle))
                            used_base_atoms.append(base_atom.name)
                            used_aa_atoms.append(aa_atom.name)

                    if not base_atom.name in used_base_atoms and not aa_atom.name in used_aa_atoms:
                        # aa with carboxylate, only count one oxygen as a donor
                        if aa_key == "ASP" or aa_key == "GLU":
                            if not carboxylate_donor_used:
                                n = n + 1
                            carboxylate_donor_used = True
                        else:
                            n = n+1
                    hydrogen_bond_list.append((base_atom.name,"H?",aa_atom.name+flip_name,h_bond_ideal_distance,distance,0))
                    used_base_atoms.append(base_atom.name)
                    used_aa_atoms.append(aa_atom.name)
                if flip == 0 and num_flip_states == 2:
                    n0 = n
                    hb0 = hydrogen_bond_list

    if num_flip_states == 2:
        if n0 >= n:       # second flip is not strictly better
            n = n0        # use the first set of hydrogen bonds
            hydrogen_bond_list = hb0
        else:
            print("  Found a flipped amino acid "+aa_residue.unit_id()+" "+base_key+" "+aa_key+" "+str(n)+" $$$$$$$$$$$$$$$$$")
#    print(aa_residue.unit_id(),hydrogen_bond_list[])

    if len(hydrogen_bond_list) > 0:
        for hbond in hydrogen_bond_list:
            print(hbond)
        print("http://rna.bgsu.edu/rna3dhub/display3D/unitid/%s,%s" % (base_residue.unit_id(),aa_residue.unit_id()))
        print("")
    return (n,hydrogen_bond_list)

def stacking_planar_annotation (base_residue, aa_residue, min_dist):
    """ For planar amino acids, determine the stacking classification
    according to the angle between the plane of the amino acid and
    the plane of the base """

    angle = calculate_angle_between_planar_residues(base_residue, aa_residue)

    # cation is about the type of amino acid.  List them ... HIS is positive sometimes.

    perpendicular_stack_aa = set(["HIS", "PHE", "TRP", "TYR"])
    perpendicular_aa = set (["HIS", "ARG", "LYS", "ASN", "GLN"])  # why is HIS in both lists?

    if angle:
        if angle <= 45:
            return ("pi-pi-stacking",{"angle-between-planes":angle})
        elif angle > 45:
            if aa_residue.sequence in perpendicular_stack_aa:
                return ("perpendicular-stacking",{"angle-between-planes":angle,"dist-xy-from-center":min_dist})
            elif aa_residue.sequence in perpendicular_aa:
                return ("cation-pi",{"angle-between-planes":angle,"dist-xy-from-center":min_dist})

    return ("other-stack",{"angle-between-planes":None,"dist-xy-from-center":min_dist})

def stacking_non_planar_annotation(aa_residue, aa_coordinates, min_dist):
    """ For non-planar amino acids, determine the stacking type
    by looking at how spread out the atoms are above/below the
    plane of the base """
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
        return ("stacked",{"stacking-diff":diff,"dist-xy-from-center":min_dist})

    return ("other-stack",{"stacking-diff":diff,"dist-xy-from-center":min_dist})

def calculate_angle_between_planar_residues (base_residue, aa_residue):
    vec1 = normal_vector_calculation(base_residue)
    vec2 = normal_vector_calculation(aa_residue)

    if len(vec1) == 3 and len(vec2) == 3:
        angle = smaller_angle_between_vectors(vec1, vec2)
        return angle
    else:
        print("Missing a normal vector %s %s" % (base_residue.unit_id(),aa_residue.unit_id()))
        return None

def normal_vector_calculation(residue):
    key = residue.sequence
    P1 = residue.centers[planar_atoms[key][0]]
    P2 = residue.centers[planar_atoms[key][1]]
    P3 = residue.centers[planar_atoms[key][2]]
#    print key, residue.unit_id(), P1, P2, P3

    if len(P1) == 3 and len(P2) == 3 and len(P3) == 3:
        normal_vector = np.cross((P2 - P1),(P3 - P1))
        return normal_vector
    else:
        return []

# this function calculates the angle made from A to B to C from 0 to 180 degrees
def calculate_hb_angle(A,B,C):
    if len(A) == 3 and len(B) == 3 and len(C) == 3:
        return angle_between_vectors(np.subtract(A,B),np.subtract(C,B))

# This function calculates an angle from 0 to 90 degrees between two vectors
def smaller_angle_between_vectors(vec1, vec2):
    if len(vec1) == 3 and len(vec2) == 3:
        # the following line sometimes causes "RuntimeWarning: invalid value encountered in double_scalars" on 5JTE
        cosang = abs(np.dot(vec1, vec2) / (np.linalg.norm(vec1) * np.linalg.norm(vec2)))
        angle = np.arccos(cosang)
        return 180*abs(angle)/np.pi
    else:
        return None

# This function calculates an angle from 0 to 180 degrees between two vectors
def angle_between_vectors(vec1, vec2):
    if len(vec1) == 3 and len(vec2) == 3:
        cosang = np.dot(vec1, vec2)
        sinang = np.linalg.norm(np.cross(vec1, vec2))
        angle = np.arctan2(sinang, cosang)
        return 180*angle/np.pi
    else:
        return None

# This function calculates an angle from 0 to 180 degrees between two vectors
def angle_between_three_points(P1,P2,P3):
    if len(P1) == 3 and len(P2) == 3 and len(P3) == 3:
        return angle_between_vectors(P1-P2,P3-P2)
    else:
        return None

# This function calculates an angle from 0 to 180 degrees between two vectors
def distance_between_vectors(vec1, vec2):
    if len(vec1) == 3 and len(vec2) == 3:
        return np.linalg.norm(np.subtract(vec1,vec2))
    else:
        return None

def detect_base_edge(base_residue, base_coordinates, aa_residue, aa_coordinates):
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

    NAsixsidedringatoms = ['N1','C2','N3','C4','C5','C6']

#    for base_atom in base_residue.atoms(name=NAbaseheavyatoms[base_residue.sequence]):
    for base_atom in base_residue.atoms(name=NAsixsidedringatoms):
        key = base_atom.name
        base_x.append(base_coordinates[key][0])
        base_y.append(base_coordinates[key][1])

    base_center_x = np.mean(base_x)
    base_center_y = np.mean(base_y)

    """
    reference_atom = {'A':'C2', 'C':'O2', 'G':'N2', 'U':'O2'}    # for WC versus Sugar edge
    reference_atom = {'A':'N6', 'C':'N4', 'G':'O6', 'U':'O4'}    # for WC versus Hoogsteen edge
    x = base_coordinates[reference_atom[base_residue.sequence]][0] - base_center_x
    y = base_coordinates[reference_atom[base_residue.sequence]][1] - base_center_y
    angle_aa = np.arctan2(y,x)         # values -pi to pi
    angle_deg = (180*angle_aa)/np.pi # values -180 to 180
    print("Base %s has S/WC reference angle %10.8f" % (base_residue.sequence,angle_deg))

    reference_atom = {'A':'N9', 'C':'N1', 'G':'N9', 'U':'N1'}    # for WC versus Hoogsteen edge
    print("Base %s has N1/N9 x value %10.8f" % (base_residue.sequence,base_coordinates[reference_atom[base_residue.sequence]][0]))
    """

    x = aa_center_x - base_center_x
    y = aa_center_y - base_center_y
    angle_aa = np.arctan2(y,x)         # values -pi to pi
    angle_deg = (180*angle_aa)/np.pi # values -180 to 180

    purine = set(["A", "G", "DA", "DG"])
    pyrimidine = set(["C", "U", "DC", "DT"])

    if base_residue.sequence in purine:
        if -14 <= angle_deg <= 104:
            return ("fgWC",angle_deg)
        elif 104 < angle_deg or aa_center_x < -1.3:
            return ("fgH",angle_deg)
        else:
            return ("fgS",angle_deg)

    elif base_residue.sequence in pyrimidine:
        if -32 <= angle_deg <= 86:
            return ("fgWC",angle_deg)
        elif 86 < angle_deg or aa_center_x < -0.37:
            return ("fgH",angle_deg)
        else:
            return ("fgS",angle_deg)

def detect_face(aa_residue, aa_coordinates):
    aa_z =[]

    for aa_atom in aa_residue.atoms(name=aa_fg[aa_residue.sequence]):
        key = aa_atom.name
        aa_z.append(aa_coordinates[key][2])

    mean_z = np.mean(aa_z)
    if mean_z <= 0:
        return ("fgs5",mean_z)
    else:
        return ("fgs3",mean_z)

def unit_vector(v):
    return v / np.linalg.norm(v)

def text_output(result_list):
    with open(outputText % PDB, 'wb') as target:
        for result in result_list:
            target.write(str(result))
            target.write("\r\n")
            target.close

def csv_output(result_list):
    with open(outputBaseAAFG % PDB, 'wb') as csvfile:
        fieldnames = ['RNA ID', 'AA ID', 'RNA Chain ID', 'RNA residue','RNA residue number','Protein Chain ID', 'AA residue','AA residue number', 'Interaction', 'Edge', 'Param']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for base_residue, aa_residue, interaction, edge, standard_aa_center, param in result_list:
            base = base_residue.unit_id()
            aa = aa_residue.unit_id()
            #print base, aa, interaction
            base_component = str(base).split("|")
            aa_component = str(aa).split("|")
            writer.writerow({'RNA ID': base, 'AA ID': aa, 'RNA Chain ID': base_component[2], \
                'RNA residue':base_component[3],'RNA residue number': base_component[4],\
                'Protein Chain ID':aa_component[2],'AA residue': aa_component[3],\
                'AA residue number': aa_component[4], 'Interaction': interaction, 'Edge': edge, 'Param': param})

        """for base_residue, aa_residue,interaction in result_list:
                    base_component = str(base_residue).split("|")
                    aa_component = str(aa_residue).split("|")
                    writer.writerow({'RNA Chain ID': base_component[2], 'RNA residue':base_component[3],\
                    'RNA residue number': base_component[4],'Protein Chain ID':ChainNames[PDB][aa_component[2]],\
                    'AA residue': aa_component[3],'AA residue number': aa_component[4], 'Interaction': interaction})"""



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
            plt.setp(base_lines, 'color', 'g', 'linewidth', 1.0)
            #ax.text(9, 1, 1, base_residue)
        except:
            print("Missing residues")
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

def draw_one_aa(aa, ax):
    #Connects atoms to draw neighboring bases and amino acids for 3D plots
    new_aa_x=[]
    new_aa_y=[]
    new_aa_z=[]

    back_aa_x=[]
    back_aa_y=[]
    back_aa_z=[]

    for atom in aa.atoms():
        ax.text(atom.x,atom.y,atom.z,atom.name)

    # This code goes from atom to atom to atom, should go pair by pair
    for atomname in aa_connections[aa.sequence]:
        coord_aa=[]
        coord_aa= aa.centers[atomname]
        new_aa_x.append(coord_aa[0])
        new_aa_y.append(coord_aa[1])
        new_aa_z.append(coord_aa[2])
    aa_lines= ax.plot(new_aa_x, new_aa_y, new_aa_z, label= 'Amino acid')
    plt.setp(aa_lines, 'color', 'r', 'linewidth', 1.0)

    for hpair in aa_hydrogen_connections[aa.sequence]:
        if hpair[0] in aa.centers and hpair[1] in aa.centers:
            first  = aa.centers[hpair[0]]
            second = aa.centers[hpair[1]]

            aa_lines = ax.plot([first[0],second[0]],[first[1],second[1]],[first[2],second[2]], label= 'Amino acid')
            plt.setp(aa_lines, 'color', 'k', 'linewidth', 1.0)

    for atomname in aa_backconnect[aa.sequence]:
        back_aa=[]
        back_aa= aa.centers[atomname]
        back_aa_x.append(back_aa[0])
        back_aa_y.append(back_aa[1])
        back_aa_z.append(back_aa[2])
    aa_lines= ax.plot(back_aa_x, back_aa_y, back_aa_z, label= 'Amino acid')
    plt.setp(aa_lines, 'color', 'y', 'linewidth', 1.0)

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

def PlotAndAnalyzeAAHydrogens(aa):
    if aa.sequence == "ARG":

        # this output checks angles
        A = unit_vector(aa.centers["NE"] - aa.centers["CZ"])
        B = unit_vector(aa.centers["NH1"] - aa.centers["CZ"])
        C = unit_vector(aa.centers["NH2"] - aa.centers["CZ"])
        print("Fitted hydrogen atoms for ARG")
        print("B-A length",np.linalg.norm(B-A))
        print("C-A length",np.linalg.norm(C-A))
        print("B-A angle",angle_between_vectors(B,A))
        print("C-B angle",angle_between_vectors(C,B))
        print("A-C angle",angle_between_vectors(A,C))
        print("HH12-NH1-HH11 angle",angle_between_three_points(aa.centers["HH12"],aa.centers["NH1"],aa.centers["HH11"]))
        print("HH22-NH2-HH21 angle",angle_between_three_points(aa.centers["HH22"],aa.centers["NH2"],aa.centers["HH21"]))
        print("CG-CD-NE  angle",angle_between_three_points(aa.centers["CG"],aa.centers["CD"],aa.centers["NE"]))
        print("CG-CD-HD2 angle",angle_between_three_points(aa.centers["CG"],aa.centers["CD"],aa.centers["HD2"]))
        print("CG-CD-HD3 angle",angle_between_three_points(aa.centers["CG"],aa.centers["CD"],aa.centers["HD3"]))
        print("NE-CD-HD2 angle",angle_between_three_points(aa.centers["NE"],aa.centers["CD"],aa.centers["HD2"]))
        print("NE-CD-HD3 angle",angle_between_three_points(aa.centers["NE"],aa.centers["CD"],aa.centers["HD3"]))

        print("Hydrogen center HH12")
        print(aa.centers["HH12"])

        print(aa.centers["HD2"])
        print(aa.centers["HD3"])
    if aa.sequence in aa_hydrogen_connections:
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.axis("equal")
        draw_one_aa(aa, ax)
        plt.title('Hydrogens added to '+aa.sequence)
        plt.show()

def WriteProteinUnits(PDB,amino_acids):
    all_aa_list = []
    for aa in amino_acids:
        fields = aa.unit_id().split("|")
        my_tuple = (aa.unit_id(),fields[2],aa.index,aa.centers['aa_fg'])
        all_aa_list.append(my_tuple)

    new_all_aa_list = sorted(all_aa_list, key=lambda x: (x[1],x[2]))

    ids = []
    chainPositions = []
    centers = []
    rotations = []
    for my_tuple in new_all_aa_list:
        if len(my_tuple[3]) == 3:
#            print(my_tuple)
            ids.append(my_tuple[0])
            chainPositions.append(my_tuple[2])
            centers.append(my_tuple[3])
            rotations.append(np.zeros((0,3,3)))
        else:
            print("missing center",my_tuple)

    dataPathUnits = "C:/Users/zirbel/Dropbox/2018 FR3D intersecting pairs/data/units/"

    with open(dataPathUnits + PDB + '_protein.pickle', 'wb') as handle:
        pickle.dump((ids,chainPositions,centers,rotations), handle, protocol = 2)  # protocol 2 is safe for Python 2.7


def writeInteractionsHTML(allInteractionDictionary,outputHTML,version,aa_part):

    SERVER = False
    SERVER = True
    if not SERVER:
        JS1 = '<script src="./js/JSmol.min.nojq.js"></script>'
        JS2 = '<script src="./js/jquery.jmolToolsRNAProtein.js"></script>'
        JS3 = '<script src="./js/imagehandlinglocal.js"></script>'
        JS4 = '<script src="./js/jmolplugin.js" type="text/javascript"></script>'
        JS5 = '<script type="text/javascript" src="./js/heatmap.js"></script>'
        JS6 = '<script src="./js/scroll_table_addon.js"></script>'   # scrollable table
        JS7 = '<script src="./js/table_js.js"></script>'   # scrollable table
        CSS1 = '<link rel="stylesheet" type="text/css" href="./css/table_style.css">'
        CSS2 = '<link rel="stylesheet" type="text/css" href="./css/scroll_table_addon.css"/>'

    else:
        JS1 = '<script src="http://rna.bgsu.edu/rna3dhub/js/jsmol/JSmol.min.nojq.js"></script>'
        JS2 = '<script src="http://rna.bgsu.edu/rna3dhub/js/jquery.jmolTools.js"></script>'
        JS2 = '<script src="./js/jquery.jmolToolsRNAProtein.js"></script>'   # special code to superimpose bases
        JS3 = '<script src="http://rna.bgsu.edu/webfr3d/js/imagehandling.js"></script>'
        JS4 = '<script src="http://rna.bgsu.edu/webfr3d/js/jmolplugin.js" type="text/javascript"></script>'
        JS5 = '<script type="text/javascript" src="http://rna.bgsu.edu/webfr3d/js/heatmap.js"></script>'
        JS6 = '<script src="./js/scroll_table_addon.js"></script>'   # scrollable table
        JS7 = '<script src="./js/table_js.js"></script>'   # scrollable table
        CSS1 = '<link rel="stylesheet" type="text/css" href="./css/table_style.css">'
        CSS2 = '<link rel="stylesheet" type="text/css" href="./css/scroll_table_addon.css"/>'

    # not working yet, so omit:
    JS6 = ""
    JS7 = ""

    count_pair = 0

    for key in allInteractionDictionary:
        pagetitle = key.replace(" ","-")
        htmlfilename = key.replace(" ","-") + version

#        print("Writing HTML file for "+key+" in "+htmlfilename+".html, found "+ str(len(allInteractionDictionary[key])) + " instances")

        fields = key.split("_")
        print(fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t"+str(len(allInteractionDictionary[key])))
        count_pair += len(allInteractionDictionary[key])

        # limit the number of instances shown, to be able to compute and display discrepancy
        numForDiscrepancy = min(300,len(allInteractionDictionary[key]))

        # calculate discrepancies between all instances, up to 300
        discrepancy = np.zeros((numForDiscrepancy,numForDiscrepancy))
        for i in range(0,numForDiscrepancy):
            instance_1 = allInteractionDictionary[key][i]
            aa_1 = instance_1[4]
            for j in range(i+1,numForDiscrepancy):
                instance_2 = allInteractionDictionary[key][j]
                aa_2 = instance_2[4]
                s = 0
                for atom_name in aa_fg[aa_1.sequence]:
                    d = distance_between_vectors(aa_1.centers[atom_name],aa_2.centers[atom_name])
                    if d:
                        s += d**2

                discrepancy[i][j] = np.sqrt(s)/len(aa_fg[aa_1.sequence])  # average distance between corresponding atoms
                discrepancy[j][i] = discrepancy[i][j]
                # discrepancy[j][i] = np.linalg.norm(standard_aa_center_1 - standard_aa_center_2)


        # base_aa = (base_residue, aa_residue, interaction, edge, standard_aa_center, param)

        # use greedy insertion 100 times to find a decent ordering of the instances
        newOrder, bestPathLength, distances = orderWithPathLengthFromDistanceMatrix(discrepancy,10)

        # rewrite the list of instances, limiting it to numForDiscrepancy
        newList = []
        for i in range(0,len(newOrder)):
            newList.append(allInteractionDictionary[key][newOrder[i]])
        allInteractionDictionary[key] = newList

        # write out text for radio boxes to display each individual interaction
        i = 1
        queryNote = "<h2>"+pagetitle+"</h2>\n"

        candidatelist = '<table style="white-space:nowrap;" id="table">\n'
        candidatelist += '<thead><tr><th><span onclick="sortTable(1)">Number</span></th>'
        candidatelist += '<th>View</th>'
        candidatelist += '<th><span onclick="sortTable(3)">Nucleotide</span></th>'
        candidatelist += '<th><span onclick="sortTable(4)">Amino acid</span></th>'
        candidatelist += '<th><span onclick="sortTable(5)">Interaction</span></th>'
        candidatelist += '<th><span onclick="sortTable(6)">Edge</span></th>'
        candidatelist += '<th><span onclick="sortTable(7)">a.a. x</span></th>'
        candidatelist += '<th><span onclick="sortTable(8)">a.a. y</span></th>'
        candidatelist += '<th><span onclick="sortTable(9)">a.a. z</span></th>'
        param = allInteractionDictionary[key][0][5]
        param_list = sorted(list(param.keys()))
        col = 9
        for header in param_list:
            col += 1
            candidatelist += '<th><span onclick="sortTable(%d)">%s</span></th>' % (col,header)
        candidatelist += "</tr></thead>\n"
        candidatelist += '<tbody id="table_rows">\n'
        for base_id, aa_id, interaction, edge, standard_aa, param in allInteractionDictionary[key]:
            candidatelist += '<tr><td>'+str(i)+'.</td><td><label><input type="checkbox" id="'+str(i-1)+'" class="jmolInline" data-coord="'
            candidatelist += base_id +","+ aa_id
            candidatelist += '">&nbsp</td>'
            candidatelist += '<td>%s</td>' % base_id
            candidatelist += '<td>%s</td>' % aa_id
            candidatelist += '<td>%s</td>' % interaction
            candidatelist += '<td>%s</td>' % edge
            candidatelist += '<td>%0.4f</td>' % standard_aa.centers[aa_part][0]
            candidatelist += '<td>%0.4f</td>' % standard_aa.centers[aa_part][1]
            candidatelist += '<td>%0.4f</td>' % standard_aa.centers[aa_part][2]
            for header in param_list:
                if isinstance(param[header],float):
                    candidatelist += '<td>%0.4f</td>' % param[header]
                elif isinstance(param[header],list):
                    for hbond in param[header]:
                        candidatelist += '<td>%s-%s-%s,%0.2fA,%0.2fA,%0.2fd</td>' % hbond
                else:
                    candidatelist += '<td>Error</td>'

            candidatelist += '</tr>\n'
            i += 1
        candidatelist += '</tbody></table>\n'
        candidatelist = candidatelist

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
        with open('template.html', 'r') as myfile:
            template = myfile.read()

        # replace ###PAGETITLE### with pagetitle
        template = template.replace("###PAGETITLE###",pagetitle)

        # replace ###CANDIDATELIST### with candidatelist
        template = template.replace("###CANDIDATELIST###",candidatelist)

        # replace ###DISCREPANCYDATA### with discrepancyText
        discrepancyText = "var data =  [\n" + discrepancyText + "]"
        discrepancyText = '<script type="text/javascript">\n' + discrepancyText + '\n</script>'
        template = template.replace("###DISCREPANCYDATA###",discrepancyText)

        template = template.replace("###QUERYNAME###",queryNote)
        template = template.replace("###JS1###",JS1)
        template = template.replace("###JS2###",JS2)
        template = template.replace("###JS3###",JS3)
        template = template.replace("###JS4###",JS4)
        template = template.replace("###REFRESH###","")
        template = template.replace("###JS5###",JS5)    # include heatmap.js code
        template = template.replace("###JS6###",JS6)
        template = template.replace("###JS7###",JS7)
        template = template.replace("###CSS1###",CSS1)
        template = template.replace("###CSS2###",CSS2)

        # write htmlfilename
        with open(outputHTML+'/'+htmlfilename+'.html', 'w') as myfile:
            myfile.write(template)

    print("Wrote out %d pairwise interactions" % count_pair)

        # upload the files to /var/www/html/RNAprotein

def writeAATwoBaseHTML(allAATwoBaseDictionary,outputHTML,version,aa_part):

    SERVER = False
    SERVER = True
    if not SERVER:
        JS1 = '<script src="./js/JSmol.min.nojq.js"></script>'
        JS2 = '<script src="./js/jquery.jmolToolsAATwoBase.js"></script>'
        JS3 = '<script src="./js/imagehandlinglocal.js"></script>'
        JS4 = '<script src="./js/jmolplugin.js" type="text/javascript"></script>'
        JS5 = '<script type="text/javascript" src="./js/heatmap.js"></script>'
        JS6 = '<script src="./js/scroll_table_addon.js"></script>'   # scrollable table
        JS7 = '<script src="./js/table_js.js"></script>'   # scrollable table
        CSS1 = '<link rel="stylesheet" type="text/css" href="./css/table_style.css">'
        CSS2 = '<link rel="stylesheet" type="text/css" href="./css/scroll_table_addon.css"/>'

    else:
        JS1 = '<script src="http://rna.bgsu.edu/rna3dhub/js/jsmol/JSmol.min.nojq.js"></script>'
        JS2 = '<script src="http://rna.bgsu.edu/rna3dhub/js/jquery.jmolTools.js"></script>'
        JS2 = '<script src="./js/jquery.jmolToolsAATwoBase.js"></script>'   # special code to superimpose bases
        JS3 = '<script src="http://rna.bgsu.edu/webfr3d/js/imagehandling.js"></script>'
        JS4 = '<script src="http://rna.bgsu.edu/webfr3d/js/jmolplugin.js" type="text/javascript"></script>'
        JS5 = '<script type="text/javascript" src="http://rna.bgsu.edu/webfr3d/js/heatmap.js"></script>'
        JS6 = '<script src="./js/scroll_table_addon.js"></script>'   # scrollable table
        JS7 = '<script src="./js/table_js.js"></script>'   # scrollable table
        CSS1 = '<link rel="stylesheet" type="text/css" href="./css/table_style.css">'
        CSS2 = '<link rel="stylesheet" type="text/css" href="./css/scroll_table_addon.css"/>'

    # not working yet, so omit:
    JS6 = ""
    JS7 = ""

    count_pair = 0

    for key in allAATwoBaseDictionary:
        pagetitle = key.replace(" ","-")
        htmlfilename = key.replace(" ","-") + version

#        print("Writing HTML file for "+key+" in "+htmlfilename+".html, found "+ str(len(allAATwoBaseDictionary[key])) + " instances")

        fields = key.split("_")
        print(key+"\t"+str(len(allAATwoBaseDictionary[key])))
        count_pair += len(allAATwoBaseDictionary[key])

        # limit the number of instances shown, to be able to compute and display discrepancy
        numForDiscrepancy = min(300,len(allAATwoBaseDictionary[key]))

        # calculate discrepancies between all instances, up to 300
        discrepancy = np.zeros((numForDiscrepancy,numForDiscrepancy))
        for i in range(0,numForDiscrepancy):
            group1 = allAATwoBaseDictionary[key][i]
            centers1 = [group1[0][1].centers[aa_part]]
            centers1.append(group1[0][0].centers["base"])
            centers1.append(group1[1][0].centers["base"])
            rotations1 = [np.empty((0,0))]
            rotations1.append(group1[0][0].rotation_matrix)
            rotations1.append(group1[1][0].rotation_matrix)

            for j in range(i+1,numForDiscrepancy):
                group2 = allAATwoBaseDictionary[key][j]
                centers2 = [group2[0][1].centers[aa_part]]
                centers2.append(group2[0][0].centers["base"])
                centers2.append(group2[1][0].centers["base"])
                rotations2 = [np.empty((0,0))]
                rotations2.append(group2[0][0].rotation_matrix)
                rotations2.append(group2[1][0].rotation_matrix)
                s = 0
                discrepancy[i][j] = matrix_discrepancy(centers1,rotations1,centers2,rotations2)
                discrepancy[j][i] = discrepancy[i][j]

        # use greedy insertion 100 times to find a decent ordering of the instances
        newOrder, bestPathLength, distances = orderWithPathLengthFromDistanceMatrix(discrepancy,100)

        # rewrite the list of instances, limiting it to numForDiscrepancy
        newList = []
        for i in range(0,len(newOrder)):
            newList.append(allAATwoBaseDictionary[key][newOrder[i]])
        allAATwoBaseDictionary[key] = newList

        # write out text for radio boxes to display each individual interaction
        i = 1
        queryNote = "<h2>"+pagetitle+"</h2>\n"

        candidatelist = '<table style="white-space:nowrap;" id="table">\n'
        candidatelist += '<thead><tr><th><span onclick="sortTable(1)">Number</span></th>'
        candidatelist += '<th>View</th>'
        candidatelist += '<th><span onclick="sortTable(3)">Amino acid</span></th>'
        candidatelist += '<th><span onclick="sortTable(4)">Nucleotide1</span></th>'
        candidatelist += '<th><span onclick="sortTable(5)">Nucleotide2</span></th>'
        candidatelist += '<th><span onclick="sortTable(6)">Edge</span></th>'
        candidatelist += '<th><span onclick="sortTable(7)">a.a. x</span></th>'
        candidatelist += '<th><span onclick="sortTable(8)">a.a. y</span></th>'
        candidatelist += '<th><span onclick="sortTable(9)">a.a. z</span></th>'
        col = 9
        param_list = []
        for header in param_list:
            col += 1
            candidatelist += '<th><span onclick="sortTable(%d)">%s</span></th>' % (col,header)
        candidatelist += "</tr></thead>\n"
        candidatelist += '<tbody id="table_rows">\n'
        for group in allAATwoBaseDictionary[key]:
            aa = group[0][1]
            base1 = group[0][0]
            base2 = group[1][0]
            candidatelist += '<tr><td>'+str(i)+'.</td><td><label><input type="checkbox" id="'+str(i-1)+'" class="jmolInline" data-coord="'
            candidatelist += aa.unit_id() +","+ base1.unit_id() +","+ base2.unit_id()
            candidatelist += '">&nbsp</td>'
            candidatelist += '<td>%s</td>' % aa.unit_id()
            candidatelist += '<td>%s</td>' % base1.unit_id()
            candidatelist += '<td>%s</td>' % base2.unit_id()
#            candidatelist += '<td>%s</td>' % edge
#            candidatelist += '<td>%0.4f</td>' % standard_aa.centers[aa_part][0]
#            candidatelist += '<td>%0.4f</td>' % standard_aa.centers[aa_part][1]
#            candidatelist += '<td>%0.4f</td>' % standard_aa.centers[aa_part][2]

            candidatelist += '</tr>\n'
            i += 1
        candidatelist += '</tbody></table>\n'
        candidatelist = candidatelist

        # write out text to tell what values to put in the heat map
        discrepancyText = ''
        for c in range(0,numForDiscrepancy):
            instance1 = allAATwoBaseDictionary[key][c][0][1].unit_id()  # id of aa
            for d in range(0,numForDiscrepancy):
                instance2 = allAATwoBaseDictionary[key][d][0][1].unit_id()  # id of aa

                discrepancyText += '{"discrepancy": ' + str(discrepancy[newOrder[c]][newOrder[d]])
                discrepancyText += ', "ife1": "' + str(c+1) + "-" + instance1 + '", "ife1_index": ' + str(c)
                discrepancyText += ', "ife2": "' + str(d+1) + "-" + instance2 + '", "ife2_index": ' + str(d) + '}'
                if c < numForDiscrepancy-1 or d < numForDiscrepancy-1:
                    discrepancyText += ',\n'

        # read template.html into one string
#        with open(outputHTML+'/localtemplate.html', 'r') as myfile:
        with open('template.html', 'r') as myfile:
            template = myfile.read()

        # replace ###PAGETITLE### with pagetitle
        template = template.replace("###PAGETITLE###",pagetitle)

        # replace ###CANDIDATELIST### with candidatelist
        template = template.replace("###CANDIDATELIST###",candidatelist)

        # replace ###DISCREPANCYDATA### with discrepancyText
        discrepancyText = "var data =  [\n" + discrepancyText + "]"
        discrepancyText = '<script type="text/javascript">\n' + discrepancyText + '\n</script>'
        template = template.replace("###DISCREPANCYDATA###",discrepancyText)

        template = template.replace("###QUERYNAME###",queryNote)
        template = template.replace("###JS1###",JS1)
        template = template.replace("###JS2###",JS2)
        template = template.replace("###JS3###",JS3)
        template = template.replace("###JS4###",JS4)
        template = template.replace("###REFRESH###","")
        template = template.replace("###JS5###",JS5)    # include heatmap.js code
        template = template.replace("###JS6###",JS6)
        template = template.replace("###JS7###",JS7)
        template = template.replace("###CSS1###",CSS1)
        template = template.replace("###CSS2###",CSS2)

        # write htmlfilename
        with open(outputHTML+'/'+htmlfilename+'.html', 'w') as myfile:
            myfile.write(template)

    print("Wrote out %d AA-two-base interactions" % count_pair)

        # upload the files to /var/www/html/RNAprotein

#=======================================================================

"""Inputs a list of PDBs of interest to generate super-imposed plots"""
PDB_List = ['5AJ3']
PDB_List = ['6hiv']
PDB_List = ['3QRQ','5J7L']
PDB_List = ['5J7L']
PDB_List = ['4V9F','4YBB','4Y4O','6AZ3','4P95']
PDB_List = ['3BT7']
PDB_List = ['5I4A']
PDB_List = ['6A2H']
PDB_List = ['4V9F']
PDB_List = ['3JB9']
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.72/3.0A/csv']
version = "_3.72_3.0"
PDB_List = ['1OCT']
PDB_List = ['4v9fFH.pdb']
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/2.5A/csv']
version = "_3.48_2.5"
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.74/4.0A/csv']
version = "_3.74_4.0"
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']
version = "_3.48_3.0"
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/view/NR_4.0_56726.45']
PDB_List = ['4V9F|1|9']
PDB_List = ['5KCR', '4WOI', '6C4I', '5JC9', '5L3P', '5KPW', '3J9Y', '3J9Z', '6BU8', '5WF0', '4V55', '4V54', '4V57', '4V56', '4V50', '4V53', '4V52', '4WF1', '5H5U', '4V5B', '5WFS', '5O2R', '5WFK', '5LZD', '5LZA', '6O9J', '6O9K', '6ORL', '6ORE', '3R8O', '3R8N', '4V85', '5MDV', '5MDW', '4V80', '4U27', '4U26', '4U25', '4U24', '4U20', '5KPS', '6GXM', '5KPX', '4U1U', '3JBU', '4V9P', '3JBV', '6Q9A', '6DNC', '4U1V', '6GXO', '5IQR', '5NWY', '4V9C', '6OSK', '4V9D', '4V9O', '5MGP', '6Q97', '3JCJ', '5J91', '3JCD', '3JCE', '6I7V', '6GXN', '4V64', '5J7L', '5AFI', '6BY1', '6ENU', '4V7V', '4V7U', '4V7T', '4V7S', '3JA1', '6ENF', '6OUO', '6ENJ', '5JU8', '5J8A', '6GWT', '4YBB', '5NP6', '5J88', '5U9G', '5U9F', '4V6D', '4V6E', '4V6C', '5JTE', '6OT3', '5J5B', '4WWW', '6OSQ', '5U4J', '5MDZ', '5U4I', '6NQB', '5UYQ', '5UYP', '5MDY', '5WDT', '6H4N', '5UYK', '4V89', '5UYM', '5UYL', '5UYN', '5WE6', '5WE4', '5KCS', '4V4Q', '4V4H', '5IT8']
PDB_List = ['4V51','4V9K']
PDB_List = ['6WJR']
PDB_List = ['4v9f']
PDB_List = ['6TPQ']
PDB_List = ['4KTG']
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.160/2.5A/csv']
version = "_3.160_2.5"
PDB_List = ['4KTG']
PDB_List = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.167/3.0A/csv']
version = "_3.167_3.0"
PDB_List = ['4V9F']


ReadPickleFile = True                  # when true, just read the .pickle file from a previous run
ReadPickleFile = False                 # when true, just read the .pickle file from a previous run

base_seq_list = ['A']
base_seq_list = ['DA','DT','DC','DG']  # for DNA
base_seq_list = ['A','U','C','G']      # for RNA

aa_list = ['HIS']
aa_list = ['ILE']
aa_list = ['ALA','VAL','ILE','LEU','ARG','LYS','HIS','ASP','GLU','ASN','GLN','THR','SER','TYR','TRP','PHE','PRO','CYS','MET']
aa_list = ['THR']
aa_list = ['ALA','VAL','ILE','LEU','ARG','LYS','HIS','ASP','GLU','ASN','GLN','SER','TYR','TRP','PHE','PRO','CYS','MET']
aa_list = ['HIS']
aa_list = ['ALA','VAL','ILE','LEU','ARG','LYS','HIS','ASP','GLU','ASN','GLN','SER','TYR','TRP','THR','PHE','PRO','CYS','MET','GLY']

atom_atom_min_distance = 4.5    # minimum atom-atom distance to note an interaction
base_aa_screen_distance = 18    #
nt_reference = "C1'"
aa_reference = "aa_fg"

hasNoProteinFilename = inputPath % "hasNoProtein.pickle"
hasNoProteinFilename = hasNoProteinFilename.replace(".cif","")
try:
    hasNoProtein = pickle.load(open(hasNoProteinFilename,"rb"))
except:
    hasNoProtein = []

# plot one instance of each of the amino acids, showing the hydrogen atoms added
PlotAA = False
PlotAA = True
PlotAA = False
AlreadyPlotted = {}

# The following lines use this program to write out centers by unit ids, for use in FR3D.
# That's not really what this program is for, though.
WriteProteinUnitsFile = True
WriteProteinInteractionFile = True
WriteProteinUnitsFile = False
WriteProteinInteractionFile = False

"""Inputs base, amino acid, aa_part of interest and cut-off distance for subsequent functions"""
if __name__=="__main__":

    timerData = myTimer("start")

    aa_part = 'aa_fg'               # other choices would be ... aa_linker and aa_backbone
    base_part = 'base'

    allInteractionDictionary = defaultdict(list)
    allAATwoBaseDictionary = defaultdict(list)

    result_nt_aa = []               # for accumulating a complete list over all PDB files

    outputDataFile = outputBaseAAFG % ""
    outputDataFile = outputDataFile.replace("aa-fg_base_.csv","RNAProtein"+version+".pickle")

    timerData = myTimer("Making PDB list",timerData)

    if ReadPickleFile:
        print("Reading " + outputDataFile)
        timerData = myTimer("Reading pickle file",timerData)
        allInteractionDictionary,allAATwoBaseDictionary,PDB_List = pickle.load(open(outputDataFile,'rb'))
        writeAATwoBaseHTML(allAATwoBaseDictionary,outputHTML,version,aa_part)
        writeInteractionsHTML(allInteractionDictionary,outputHTML,version,aa_part)
    else:
        PDB_IFE_Dict = defaultdict(str)      # accumulate PDB-IFE pairs
        for PDB in PDB_List:
            if "nrlist" in PDB and "NR_" in PDB:            # referring to an equivalence class online
                                          # download the entire representative set,
                                          # then find the right line for the equivalence class
                                          # then extract the list
                f = urllib.urlopen(PDB)
                myfile = f.read()
                alltbody = myfile.split("tbody")
                alllines = alltbody[1].split("<a class='pdb'>")
                del alllines[0]
                for line in alllines:
                    fields = line.split("</a>")
                    print(fields[0])
                    if len(fields[0]) > 1:
                        newIFE = fields[0].replace(" ","")   # remove spaces
                        newPDB = newIFE[0:4]
                        PDB_IFE_Dict[newPDB] += "+" + newIFE

            elif "nrlist" in PDB:           # referring to a representative set online
                f = urllib.urlopen(PDB)
                myfile = f.read()
                alllines = myfile.split("\n")
                for line in alllines:
                    fields = line.split(",")

                    if len(fields) > 1 and len(fields[1]) > 4:
                        newPDB = fields[1][1:5]       # use only PDB identifier, ignore IFE for now
                        PDB_IFE_Dict[newPDB] += "+" + fields[1].replace('"','')

            elif "+" in PDB:
                newPDB = PDB[0:4]
                PDB_IFE_Dict[newPDB] = PDB

            elif "|" in PDB:
                newPDB = PDB[0:4]
                PDB_IFE_Dict[newPDB] = PDB

            else:
                PDB_IFE_Dict[PDB] = ""            # indicates to process the whole PDB file

        print(PDB_IFE_Dict)


        counter = 0
        count_pair = 0

        # loop through 3D structures and annotate interactions
        PDBs = PDB_IFE_Dict.keys()
        #PDBs = PDBs[::-1]  #reverse the order of the list

        for PDB in PDBs:
            counter += 1

            if 0 > 1 and PDB in hasNoProtein:
                print("Skipping file " + PDB + " since it has no amino acids")
                continue

            print("Reading file " + PDB + ", which is number "+str(counter)+" out of "+str(len(PDB_IFE_Dict)))
            timerData = myTimer("Reading CIF files",timerData)

            structure = get_structure(inputPath % PDB)

            try:
                #structure = get_structure(inputPath % PDB)
                aaa=1
            except:
                print("Could not load structure")
                continue

            IFE = PDB_IFE_Dict[PDB]
            if len(IFE) == 0:
                bases = structure.residues(sequence= base_seq_list)  # load just the types of bases in base_seq_list
            else:
                chain_ids = []
                print(IFE)
                chains = IFE.split("+")
                for chain in chains[1:]:            #skip element zero, leading +
                    fields = chain.split("|")
                    chain_ids.append(fields[2])
                bases = structure.residues(chain = chain_ids, sequence= base_seq_list)

            numBases = 0
            for base in bases:
                numBases += 1

            amino_acids = structure.residues(sequence=aa_list)

            # if desired, write files for FR3D
            if WriteProteinUnitsFile:
                timerData = myTimer("Writing units file",timerData)
                WriteProteinUnits(PDB,amino_acids)

            numAA = 0
            for aa in amino_acids:
                numAA += 1

                if PlotAA and not aa.sequence in AlreadyPlotted:
                    timerData = myTimer("Plotting amino acids",timerData)
                    PlotAndAnalyzeAAHydrogens(aa)
                    AlreadyPlotted[aa.sequence] = 1

            print("  Found " + str(numBases) + " bases and " + str(numAA) + " amino acids in " + PDB)

            if numAA == 0:
                hasNoProtein.append(PDB)
                pickle.dump(hasNoProtein,open(hasNoProteinFilename,"wb"))

            # find atom-atom contacts
            timerData = myTimer("Finding atom-atom contacts",timerData)
            contact_list = find_atom_atom_contacts(bases,amino_acids,atom_atom_min_distance)

            # find neighbors in order to annotate base-aa_fg interactions
            timerData = myTimer("Finding neighbors",timerData)
            print("  Finding neighbors in " + PDB)
            nt_reference = "base"
            aa_reference = "aa_fg"
            base_aa_screen_distance = 10
            baseCubeList, baseCubeNeighbors, aaCubeList = find_neighbors(bases, amino_acids, base_aa_screen_distance, PDB_IFE_Dict[PDB], nt_reference, aa_reference)

            # annotate base-aa_fg interactions
            timerData = myTimer("Annotating interactions",timerData)
            list_base_aa, list_aa_coord, list_base_coord, hbond_aa_dict = annotate_interactions(bases, amino_acids, base_aa_screen_distance, baseCubeList, baseCubeNeighbors, aaCubeList)

            timerData = myTimer("Recording interactions",timerData)

            # accumulate list of interacting units by base, amino acid, interaction type, and edges
            for base_residue, aa_residue, interaction, edge, standard_aa, param in list_base_aa:
                base = base_residue.unit_id()
                # skip symmetry operated instances; generally these are just duplicates anyway
                if not "||||" in str(base):
                    aa = aa_residue.unit_id()
                    base_component = str(base).split("|")
                    aa_component = str(aa).split("|")
                    key = base_component[3]+"_"+aa_component[3]+"_"+interaction+"_"+edge
                    count_pair += 1
                    allInteractionDictionary[key].append((base,aa,interaction,edge,standard_aa,param))  # store tuples

            # accumulate a list of amino acids hydrogen bonding with more than one base
            for aa_unit_id in hbond_aa_dict:
                if len(hbond_aa_dict[aa_unit_id]) > 1:
                    # skip symmetry operated instances; generally these are just duplicates anyway
                    if not "||||" in str(aa_unit_id):
                        group = hbond_aa_dict[aa_unit_id]
                        key = group[0][1].sequence + "-" + group[0][0].sequence + "-" + group[1][0].sequence
                        key = group[0][1].sequence
                        allAATwoBaseDictionary[key].append(group)

            """ 3D plots of base-aa interactions
            for base, aa, interaction in list_base_aa:
                base_seq = base.sequence
                aa= aa.sequence

                draw_base(base_seq, ax)
                draw_aa(aa, ax)
                #draw_aa_cent(aa, aa_part, ax)

                ax.set_xlabel('X Axis')
                ax.set_ylabel('Y Axis')
                ax.set_zlabel('Z Axis')
                ax.set_xlim3d(10, -15)
                ax.set_ylim3d(10, -15)
                ax.set_zlim3d(10, -15)
                plt.title('%s with ' % base_seq +'%s' % aa + ' %s' % aa_part)
                plt.show()
                          """
            #accumulate a full list of resultant RNA-aa pairs
    #        result_nt_aa.extend(list_base_aa)

            #writing out output files
            #text_output(result_nt_aa)

            csv_output(list_base_aa)
            print("  Wrote output to " + outputBaseAAFG % PDB)

            with open(contact_list_file % PDB, 'w') as clf:
                clf.writelines(contact_list)

            myTimer("summary",timerData)

        print("Recorded %d pairwise interactions" % count_pair)

        # when appropriate, write out HTML files
        if len(PDB_IFE_Dict) > 100:
            print("Writing " + outputDataFile)
            timerData = myTimer("Writing HTML files",timerData)
            pickle.dump((allInteractionDictionary,allAATwoBaseDictionary,PDB_List),open(outputDataFile,"wb"))
            writeAATwoBaseHTML(allAATwoBaseDictionary,outputHTML,version,aa_part)
            writeInteractionsHTML(allInteractionDictionary,outputHTML,version,aa_part)

    myTimer("summary",timerData)
