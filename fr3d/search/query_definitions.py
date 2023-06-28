from collections import defaultdict
from query_processing import emptyInteractionList
from query_processing import emptyInteractionMatrix

import numpy as np
import math
import urllib
import json
import sys
import pickle
import string
import os


# defineuserquery has a set of predefined queries
def defineUserQuery(name):

    Q = defaultdict(dict)

    # if name is not already defined, use the last of the following
    if len(name) == 0:
        name = "AU cWW"
        name = "SR triple"
        name = "SR core"
        name = "stacked cWW"
        name = "sarcin3geometric"
        name = "cWWpairWithAminoAcids"
        name = "GU 2 into junction"
        name = "bSS and cWW"
        name = "bSS"

        name = "GU Case 1 UG cWW"
        name = "continuity test"
        name = "unary test"

        name = "UGGU tandem"
        name = "GUUG tandem"
        name = "GGUC_GGUC tandem"
        name = "GU tandem at end of helix"
        name = "GU tandem at end of helix adjacent above"
        name = "GU tandem at end of helix adjacent below"
        name = "GU tandem at end of helix one unpaired"
        name = "GU tandem at end of helix two unpaired"

        name = "UG tandem at end of helix"
        name = "UG tandem at end of helix adjacent above"
        name = "UG tandem at end of helix adjacent below"
        name = "UG tandem at end of helix one unpaired"
        name = "UG tandem at end of helix two unpaired"

        name = ""
        name = "AG tHS"
        name = "tSS LR"
        name = "bSS distance 5,6"
        name = "bPh"

        name = "stacked bases with aa"
        name = "base surrounded by bases"
        name = "base surrounded by amino acids"
        name = "base almost surrounded by amino acids"
        name = "junction 10"
        name = "internal loop"
        name = "Hairpin flanking pair"
        name = "Hairpin interacts with single strand"
        name = "borderSS"
        name = "SR triple"
        name = "cWW with amino acid in minor groove"

        name = "UU Case 0"
        name = "UU Case 2 UU interior"
        name = "UU Case 2 UU cWW"
        name = ["UU Case 0","UU Case 2 UU interior","UU Case 2 UU cWW"]

        name = "GC Case 1b CG cWW"
        name = "GC Case 1c CG cWW"
        name = "GC Case 1 CG interior"
        name = "GC Case 2 GC cWW"
        name = "GC Case 2b GC cWW"
        name = "GC Case 2c GC cWW"
        name = "GC Case 2 GC interior"
        name = "GC Case 0"
        name = "GC Case 1 CG cWW"
        name = ["GC Case 0","GC Case 1 CG cWW","GC Case 1b CG cWW","GC Case 1c CG cWW","GC Case 1 CG interior","GC Case 2 GC cWW","GC Case 2b GC cWW","GC Case 2c GC cWW","GC Case 2 GC interior"]

        name = "AU Case 2 AU cWW"
        name = "AU Case 2b AU cWW"
        name = "AU Case 2c AU cWW"
        name = "AU Case 2 AU interior"
        name = "AU Case 0"
        name = "AU Case 1 UA cWW"
        name = "AU Case 1b UA cWW"
        name = "AU Case 1c UA cWW"
        name = "AU Case 1 UA interior"
        name = ["AU Case 0","AU Case 1 UA cWW","AU Case 1b UA cWW","AU Case 1c UA cWW","AU Case 1 UA interior","AU Case 2 AU cWW","AU Case 2b AU cWW","AU Case 2c AU cWW","AU Case 2 AU interior"]

        name = "GU Case 0"
        name = "GU Case 1 UG cWW"
        name = "GU Case 1b UG cWW"
        name = "GU Case 1c UG cWW"
        name = "GU Case 1 UG interior"
        name = "GU Case 2 GU cWW"
        name = "GU Case 2b GU cWW"
        name = "GU Case 2c GU cWW"
        name = "GU Case 2 GU interior"
        name = ["GU Case 0","GU Case 1 UG cWW","GU Case 1b UG cWW","GU Case 1c UG cWW","GU Case 1 UG interior","GU Case 2 GU cWW","GU Case 2b GU cWW","GU Case 2c GU cWW","GU Case 2 GU interior"]

        name = "AU hairpin Case 1 UA cWW"
        name = "AU hairpin Case 1 UA after cWW"
        name = "GU Case 2b GU cWW"

        name = "unary no pair"

        name = "CC cWW almost"
        name = 'tHH BPh'
        name = "kink turn 65553"

        name = "sarcin5geometric"
        name = "AG tHS"
        name = "Decoding loop"
        name = "Z step n+"
        name = "cWW in NMR"
        name = "cWW with amino acid in minor groove"
        name = "RNA-protein-cWW-minor-groove"
        name = "GNRA hairpin symbolic"
        name = "sO4'5"
        name = "sO4'3"
        name = "Z step"
        name = "IL_4Y4O_235"
        name = "syn pair"
        name = "syn stack"
        name = "NUNNGN tSW"
        name = "next"
        name = "chi_angle"
        name = "RNA-protein4"
        name = 'modified'
        name = 'cWW with modified'
        name = "sO"
        name = "Python modified bases"
        name = "LR cWW"
        name = "chain length"
        name = "count BP"
        name = "count NT"
        name = 'Compare Matlab Python'
        name = "AA cWW"
        name = "DNA"
        name = "sarcin5geometric"
        name = "sarcin5mixed"
        name = "sarcin13mixed"
        name = "GoU 6S0Z"
        name = "GoG 5NJT"
        name = "cWW with modified"
        name = "Stacked cWW GC"
        name = 'cWW in DNA'
        name = 'cWW stack in DNA'
        name = "sarcin5mixed2"

    if name == 'Python modified bases':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = 'RNA'
        Q["interactionMatrix"][0][1] = "cWW_exp"
        Q["interactionMatrix"][1][1] = 'RNA modified'
        #Q["interactionMatrix"][1][0] = '>'
        Q["searchFiles"] = ['4TNA']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.237/2.5A/csv']   # list of IFEs to search

    elif name == 'Compare Matlab Python':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = 'N'
        Q["interactionMatrix"][0][1] = "tWW tWW_exp"
        Q["interactionMatrix"][1][1] = 'N'
        Q["interactionMatrix"][1][0] = '>'
        Q["searchFiles"] = ['4TNA']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.237/2.5A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['7K00']   # set of IFEs to search

    elif name == 'DNA':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = 'DNA'
        Q["interactionMatrix"][0][1] = "pair_exp"
        Q["interactionMatrix"][1][1] = 'RNA DNA'
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.237/2.5A/csv']   # list of IFEs to search

    elif name == 'chi_angle':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = "chi_100_120"
        Q["interactionMatrix"][0][0] = "chi_170_-170"
        Q["interactionMatrix"][0][0] = "chi(80:85)"
        Q["interactionMatrix"][1][0] = "next"
        Q["searchFiles"] = ['7K00','4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['6AZ3']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.233/3.0/csv']   # list of IFEs to search

    elif name == 'count BP':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = "chainlength_100_700"
        Q["interactionMatrix"][1][1] = "chainlength_100_700"
        Q["interactionMatrix"][0][1] = "cWW AU UA GC CG GU UG"
        Q["searchFiles"] = ['4ARC']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F','5J7L','6ZMI']   # set of IFEs to search
        Q["searchFiles"] = ['7JQQ']
        Q["searchFiles"] = ['4TNA']
        Q["searchFiles"] = ['4tna_local']
        Q["searchFiles"] = ['7azs-local']
        Q["searchFiles"] = ['4tna-local']
        Q["searchFiles"] = ['6K0A-local_now']
        Q["searchFiles"] = ['4V9F']
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.242/2.5A/csv']   # list of IFEs to search

    elif name == 'count NT':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 1
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = "chainlength_2230_inf"
        Q["interactionMatrix"][0][0] = "chainlength_730_2229"
        Q["interactionMatrix"][0][0] = "chainlength_100_729"
        Q["interactionMatrix"][0][0] = "chainlength_140_729"
        Q["interactionMatrix"][0][0] = "chainlength_1_139"
        Q["searchFiles"] = ['4V9F']
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.243/3.0A/csv']   # list of IFEs to search

    elif name == 'AG tHS':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2

        directlySpecify = True    # not working as of 10-24-2021, not sure why
        directlySpecify = False
        if directlySpecify:
            Q["requiredInteractions"] = emptyInteractionMatrix(Q["numpositions"])
            Q["requiredInteractions"][0][1] = ["tHS","and"]
            Q['activeInteractions'] = ['tHS']
            Q["errorMessage"] = []
            Q["requiredMoleculeType"] = [['RNA'], ['RNA']]
            Q["combinationConstraint"] = emptyInteractionMatrix(Q["numpositions"])
            Q["combinationConstraint"][0][1] = [('A', 'G')]
        else:
            Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
            Q["interactionMatrix"][0][1] = "AG tHS"
            Q["interactionMatrix"][1][0] = "=4,-4"

        Q["searchFiles"] = ['4V9F|1|9']   # set of IFEs to search
        Q["searchFiles"] = ['1A1T|1|B','1A34|1|B','1A34|1|C','17RA|7|A'];  # for testing new format 3/12/19
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['1S72|1|9']
        Q["searchFiles"] = ['4V9F','5J7L','6ZMI']  # archeal, bacterial, and eukaryotic ribosomes
        Q["searchFiles"] = ['4V9F|1|9+4V9F|1|0', '5J7L|1|CA+5J7L|1|DA+5J7L|1|DB+5J7L|1|AA+5J7L|1|BA+5J7L|1|CB']

    elif name == 'next':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][1][0] = "next"
        Q["searchFiles"] = ['4ARC']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F','5J7L','6ZMI']   # set of IFEs to search
        Q["searchFiles"] = ['3AM1|1|B']   # set of IFEs to search

    elif name == 'syn pair':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "pair"
        Q["interactionMatrix"][0][0] = "syn"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search

    elif name == 'syn stack':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "stack"
        Q["interactionMatrix"][0][0] = "syn"
        Q["interactionMatrix"][0][0] = "~anti"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['2N1Q']   # set of IFEs to search

    elif name == 'Z step':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "s33 and sO4'3"
        Q["interactionMatrix"][1][0] = "next"
        Q["interactionMatrix"][1][1] = "glyco"
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.201/all/csv']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['3DIL', '4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['7K00', '5OB3']   # set of IFEs to search

    elif name == 'NUNNGN tSW':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][1][4] = "UG tSW"
        Q["interactionMatrix"][1][0] = "next"
        Q["interactionMatrix"][2][1] = "next"
        Q["interactionMatrix"][3][2] = "next"
        Q["interactionMatrix"][4][3] = "next"
        Q["interactionMatrix"][5][4] = "next"
        Q["interactionMatrix"][4][4] = "glyco"
        Q["searchFiles"] = ['7K00','5Y85','4ARC','5Y85','5T5H','6DCB','5OB3','3U4M','7KKV','7EOG']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.230/2.0A/csv']   # list of IFEs to search

    elif name == "Stacked cWW GC":
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 4
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][3] = "cWW CG GC"
        Q["interactionMatrix"][1][2] = "cWW CG GC"
        Q["interactionMatrix"][1][0] = "next"
        Q["interactionMatrix"][3][2] = "next"
        Q["interactionMatrix"][2][1] = ">"
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.230/2.0A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['7K00']   # set of IFEs to search

    elif name == "sO4'5":
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "sO4'5 nsO4'5"
        Q["interactionMatrix"][0][1] = "sO4'3"
        Q["searchFiles"] = ['3BNS']   # set of IFEs to search
        Q["searchFiles"] = ['7K00']   # set of IFEs to search

    elif name == "sO":
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = "glyco"
        Q["interactionMatrix"][0][1] = "s3O s5O"
        Q["searchFiles"] = ['3BNS']   # set of IFEs to search
        Q["searchFiles"] = ['1EG0']   # set of IFEs to search
        Q["searchFiles"] = ['4TNA']   # set of IFEs to search
        Q["searchFiles"] = ['7OTC']   # set of IFEs to search
        Q["searchFiles"] = ['7K00']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.221/3.0A/csv']   # list of IFEs to search

    elif name == "sO4'3":
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = "glyco"
        Q["interactionMatrix"][0][1] = "sOP3"
        Q["interactionMatrix"][1][0] = "next"
        Q["searchFiles"] = ['7K00|1|a','7K00|1|A']   # set of IFEs to search

    elif name == 'Z step n+':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][1][0] = "n+sO =1 >"
        Q["interactionMatrix"][1][1] = "syn anti is"
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.201/all/csv']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['3DIL', '4V9F', '1S72']   # set of IFEs to search

    elif name == 'tHH BPh':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "tHH BPh"
        Q["searchFiles"] = ['5J7L']   # set of IFEs to search

    elif name == 'AA cWW':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "AA cWW"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.224/2.5A/csv']   # list of IFEs to search


    elif name == 'continuity test':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = "G"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][0][1] = "cWW"
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F|1|9']   # set of IFEs to search

    elif name == 'tSS LR':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "tSS LR"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search

    elif name == 'borderSS':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "bSS"
        Q["interactionMatrix"][1][0] = ">"  # break the symmetry
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search

    elif name == 'bSS distance 5,6':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "bSS"
        Q["interactionMatrix"][1][0] = "=5,6"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search

    elif name == 'bPh':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "BPh"
        Q["interactionMatrix"][1][0] = "BPh"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search

    elif name == 'unary test':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = "G"
        Q["interactionMatrix"][1][1] = "C"
        Q["interactionMatrix"][0][1] = "cWW"
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F|1|9']   # set of IFEs to search

    elif name == 'AU cWW':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["requiredInteractions"] = emptyInteractionMatrix(Q["numpositions"])
        Q["requiredInteractions"][0][1] = ["cWW","and"]
        Q["activeInteractions"] = ["cWW"]
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['1A1T|1|B','1A34|1|B','1A34|1|C','17RA|7|A'];  # for testing new format 3/12/19
        Q["searchFiles"] = ['4V9F|1|9']   # set of IFEs to search

    elif name == 'LR cWW':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "LR cWW"
        Q["interactionMatrix"][0][1] = "crossing_3_3 cWW"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][0][0] = "N"
        Q["interactionMatrix"][1][1] = "N"
        Q["searchFiles"] = ['1A1T|1|B','1A34|1|B','1A34|1|C','17RA|7|A'];  # for testing new format 3/12/19
        Q["searchFiles"] = ['4V9F|1|0']   # set of IFEs to search

    elif name == 'chain length':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = "N chainlength_2000_inf"
        Q["interactionMatrix"][1][1] = "N chainlength_0_500"
        Q["interactionMatrix"][0][1] = "cWW"
        Q["searchFiles"] = ['6ZMI']   # set of IFEs to search

    elif name == 'SR triple':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 3
        Q["searchFiles"] = ['4V9F|1|9']   # set of IFEs to search
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "cSH"
        Q["interactionMatrix"][1][2] = "tWH"
        Q["searchFiles"] = ['4V9F']   # list of IFEs to search

    elif name == 'cWW in NMR':
        Q["name"] = name
        Q["type"] = "geometric"
        Q["unitID"] = ["2N1Q|2|A|G|234","2N1Q|2|A|C|284"]
        Q["numpositions"] = len(Q["unitID"])
        Q["discrepancy"] = 0.8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "GC CG"
        Q["searchFiles"] = ['2N1Q']   # list of IFEs to search

    elif name == 'cWW in DNA':
        Q["name"] = name
        Q["type"] = "geometric"
        Q["unitID"] = ["2N1Q|2|A|G|234","2N1Q|2|A|C|284"]
        Q["numpositions"] = len(Q["unitID"])
        Q["discrepancy"] = 0.4
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = "DA DC DG DT"
        Q["interactionMatrix"][1][1] = "DA DC DG DT"
        Q["interactionMatrix"][0][0] = "DNA"
        Q["interactionMatrix"][1][1] = "DNA"
        Q["interactionMatrix"][0][1] = "DA,DT DC,DG"

        Q["searchFiles"] = ['7JQQ']   # list of IFEs to search

    elif name == 'cWW stack in DNA':
        Q["name"] = name
        Q["type"] = "geometric"
        Q["unitID"] = ["2N1Q|2|A|G|234","2N1Q|2|A|C|284","2N1Q|1|A|C|233","2N1Q|1|A|G|285"]
        Q["numpositions"] = len(Q["unitID"])
        Q["discrepancy"] = 0.5
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = "DA DC DG DT"
        Q["interactionMatrix"][1][1] = "DA DC DG DT"
        Q["interactionMatrix"][0][0] = "DNA"
        Q["interactionMatrix"][1][1] = "DNA"
        Q["interactionMatrix"][2][2] = "DNA"
        Q["interactionMatrix"][3][3] = "DNA"

        Q["searchFiles"] = ['7JQQ']   # list of IFEs to search

    elif name == "GoU 6S0Z":
        Q["name"] = name
        Q["type"] = "geometric"
        Q["unitID"] = ["6S0Z|1|A|G|2331","6S0Z|1|A|U|2339"]
        Q["numpositions"] = len(Q["unitID"])
        Q["discrepancy"] = 1.0
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "GU"
        Q["searchFiles"] = ['8b0x|1|a','7K00|1|a','5E81','4Y4O','5NJT','5NGM','6S0Z']   # list of IFEs to search

    elif name == "GoG 5NJT":
        Q["name"] = name
        Q["type"] = "geometric"
        Q["unitID"] = ["5NJT|1|A|G|673","5NJT|1|A|G|750"]
        Q["numpositions"] = len(Q["unitID"])
        Q["discrepancy"] = 1.0
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "GG"
        Q["searchFiles"] = ['8b0x','7K00','5E81','4Y4O','5NJT','5NGM','6S0Z']   # list of IFEs to search

    elif name == 'cWW with modified':
        Q["name"] = name
        Q["type"] = "geometric"
        Q["unitID"] = ["2N1Q|2|A|G|234","2N1Q|2|A|C|284"]
        Q["numpositions"] = len(Q["unitID"])
        Q["discrepancy"] = 0.4
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = "RNA modified"
        Q["searchFiles"] = ['2N1Q']   # list of IFEs to search
        Q["searchFiles"] = ['4TNA']
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.267/2.5A/csv']   # list of IFEs to search

    elif name == 'modified':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = "RNA modified"
        Q["interactionMatrix"][1][1] = "RNA"
        Q["interactionMatrix"][1][0] = "next"
        Q["searchFiles"] = ['2N1Q']   # list of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.233/2.5A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['4TNA']

    elif name == 'CC cWW almost':
        Q["name"] = name
        Q["type"] = "geometric"
        Q["unitID"] = ["4V9F|1|0|C|356","4V9F|1|0|C|295"]
        Q["numpositions"] = len(Q["unitID"])
        Q["discrepancy"] = 0.7
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "CC"
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.167/2.5A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F|1|0','4V9F|1|9']   # list of IFEs to search

    elif name == 'sarcin5geometric':
        Q["name"] = name
        Q["type"] = "geometric"
        Q["unitID"] = ["4V9F|1|0|G|2692","4V9F|1|0|U|2693","4V9F|1|0|A|2694","4V9F|1|0|G|2701","4V9F|1|0|A|2702"]
        Q["numpositions"] = len(Q["unitID"])
        Q["discrepancy"] = 0.3
        Q["searchFiles"] = ['4V9F|1|0','4V9F|1|9']   # list of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F']   # list of IFEs to search

    elif name == 'sarcin3geometric':
        Q["name"] = name
        Q["type"] = "geometric"
        Q["numpositions"] = 3
        Q["unitID"] = [None] * Q["numpositions"]
        Q["unitID"][0] = "4V9F|1|0|G|2692"
        Q["unitID"][1] = "4V9F|1|0|U|2693"
        Q["unitID"][2] = "4V9F|1|0|A|2702"
        Q["unittype"] = ["RNA"] * Q["numpositions"]
        Q["discrepancy"] = 0.3
        Q["queryMoleculeType"][0] = "RNA"
        Q["queryMoleculeType"][1] = "RNA"
        Q["queryMoleculeType"][2] = "RNA"
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["requiredMoleculeType"][2] = ["RNA"]
        Q["searchFiles"] = ['4V9F|1|0']   # list of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F|1|0','4V9F|1|9']   # list of IFEs to search

    elif name == 'sarcin5mixed':
        Q["name"] = name
        Q["type"] = "mixed"
        Q["numpositions"] = 5
        Q["unitID"] = [None] * Q["numpositions"]
        Q["unitID"][0] = "4V9f|1|0|G|2692"
        Q["unitID"][1] = "4v9F|1|0|U|2693"
        Q["unitID"][2] = "4V9F|1|0|A|2694"
        Q["unitID"][3] = "4v9f|1|0|G|2701"
        Q["unitID"][4] = "4V9F|1|0|A|2702"
        Q["unittype"] = ["RNA"] * Q["numpositions"]
        Q["queryMoleculeType"][0] = "RNA"   # help FR3D identify where to get the query data
        Q["queryMoleculeType"][1] = "RNA"
        Q["queryMoleculeType"][2] = "RNA"
        Q["queryMoleculeType"][3] = "RNA"
        Q["queryMoleculeType"][4] = "RNA"

        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["requiredMoleculeType"][2] = ["RNA"]
        Q["requiredMoleculeType"][3] = ["RNA"]
        Q["requiredMoleculeType"][4] = ["RNA"]

        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "cSH"
        Q["interactionMatrix"][1][4] = "tWH"

        Q["discrepancy"] = 0.5
        Q["searchFiles"] = ['4V9F|1|0']   # list of IFEs to search
        Q["searchFiles"] = ['7k00']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F','7K00','4Y4O','6ZMI']   # list of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.259/2.5A/csv']   # list of IFEs to search

    elif name == 'sarcin13mixed':
        Q["name"] = name
        Q["type"] = "mixed"
        Q["numpositions"] = 13
        Q["unitID"] = "4V9F|1|0|A|2689,4V9F|1|0|U|2690,4V9F|1|0|A|2691,4V9F|1|0|G|2692,4V9F|1|0|U|2693,4V9F|1|0|A|2694,4V9F|1|0|C|2695,4V9F|1|0|G|2700,4V9F|1|0|G|2701,4V9F|1|0|A|2702,4V9F|1|0|A|2703,4V9F|1|0|C|2704,4V9F|1|0|U|2705".split(",")
        Q["numpositions"] = len(Q["unitID"])
        Q["unittype"] = ["RNA"] * Q["numpositions"]

        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][12] = "cWW ncWW"
        Q["interactionMatrix"][6][7] = "cWW ncWW"
        Q["interactionMatrix"][0][6] = "bSS"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][3][2] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = ">"
        Q["interactionMatrix"][6][5] = ">"

        Q["interactionMatrix"][1][0] = "next"
        Q["interactionMatrix"][2][1] = "next"
        Q["interactionMatrix"][3][2] = "next"
        Q["interactionMatrix"][4][3] = "next"
        Q["interactionMatrix"][5][4] = "next"
        Q["interactionMatrix"][6][5] = "next"

        Q["discrepancy"] = 0.5
        Q["searchFiles"] = ['7k00']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F','7K00','4Y4O','6ZMI']   # list of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.259/2.5A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F']   # list of IFEs to search

    elif name == 'IL_4Y4O_235':
        Q["name"] = name
        Q["type"] = "mixed"
        Q["numpositions"] = 5
        Q["unitID"] = [None] * Q["numpositions"]
# {0: '4Y4O|1|2A|C|1450|||A', 1: '4Y4O|1|2A|C|1451', 2: '4Y4O|1|2A|G|1459', 3: '4Y4O|1|2A|G|1461', 4: '4Y4O|1|2A|A|1460'}
        Q["unitID"][0] = "4Y4O|1|2A|C|1450|||A"
        Q["unitID"][1] = "4Y4O|1|2A|C|1451"
        Q["unitID"][2] = "4Y4O|1|2A|G|1459"
        Q["unitID"][3] = "4Y4O|1|2A|A|1460"
        Q["unitID"][4] = "4Y4O|1|2A|G|1461"
        Q["unittype"] = ["RNA"] * Q["numpositions"]
        Q["queryMoleculeType"][0] = "RNA"   # help FR3D identify where to get the query data
        Q["queryMoleculeType"][1] = "RNA"
        Q["queryMoleculeType"][2] = "RNA"
        Q["queryMoleculeType"][3] = "RNA"
        Q["queryMoleculeType"][4] = "RNA"

        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["requiredMoleculeType"][2] = ["RNA"]
        Q["requiredMoleculeType"][3] = ["RNA"]
        Q["requiredMoleculeType"][4] = ["RNA"]

        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][4] = "cWW"
        Q["interactionMatrix"][1][2] = "cWW"

        Q["discrepancy"] = 0.1
        Q["searchFiles"] = ['4Y4O']   # list of IFEs to search

    elif name == 'sarcin5mixed2':
        Q["name"] = name
        Q["type"] = "mixed"
        Q["numpositions"] = 5
        Q["unitID"] = [None] * Q["numpositions"]
        Q["unitID"][0] = "4V9F|1|0|G|2692"
        Q["unitID"][1] = "4V9F|1|0|U|2693"
        Q["unitID"][2] = "4V9F|1|0|A|2694"
        Q["unitID"][3] = "4V9F|1|0|G|2701"
        Q["unitID"][4] = "4V9F|1|0|A|2702"
        Q["unittype"] = ["RNA"] * Q["numpositions"]
        Q["queryMoleculeType"][0] = "RNA"
        Q["queryMoleculeType"][1] = "RNA"
        Q["queryMoleculeType"][2] = "RNA"
        Q["queryMoleculeType"][3] = "RNA"
        Q["queryMoleculeType"][4] = "RNA"

        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["requiredMoleculeType"][2] = ["RNA"]
        Q["requiredMoleculeType"][3] = ["RNA"]
        Q["requiredMoleculeType"][4] = ["RNA"]

        Q["requiredInteractions"] = emptyInteractionList(Q["numpositions"])
#        Q["requiredInteractions"][0][1] = ["cSH"]
        Q["requiredInteractions"][1][4] = ["tWH"]
        Q["discrepancy"] = 0.5
        Q["searchFiles"] = ['4V9F|1|0']   # list of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F']   # list of IFEs to search

    elif name == 'Decoding loop':
        Q["name"] = name
        Q["type"] = "mixed"
        Q["numpositions"] = 11
        Q["unitID"] = [None] * Q["numpositions"]
        Q["unitID"][0]  = "4V50|1|AA|G|1405"
        Q["unitID"][1]  = "4V50|1|AA|U|1406"
        Q["unitID"][2]  = "4V50|1|AA|C|1407"
        Q["unitID"][3]  = "4V50|1|AA|A|1408"
        Q["unitID"][4]  = "4V50|1|AA|C|1409"
        Q["unitID"][5]  = "4V50|1|AA|G|1491"
        Q["unitID"][6]  = "4V50|1|AA|A|1492"
        Q["unitID"][7]  = "4V50|1|AA|A|1493"
        Q["unitID"][8]  = "4V50|1|AA|G|1494"
        Q["unitID"][9]  = "4V50|1|AA|U|1495"
        Q["unitID"][10] = "4V50|1|AA|C|1496"
        Q["unittype"] = ["RNA"] * Q["numpositions"]
        Q["queryMoleculeType"] = ["RNA"] * Q["numpositions"]
        Q["requiredMoleculeType"] = ["RNA"] * Q["numpositions"]

        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][10] = "cWW"
        Q["discrepancy"] = 0.1
        Q["searchFiles"] = ['4V50']   # list of IFEs to search

    elif name == 'RNA-protein4':
        Q["name"] = name
        Q["type"] = "geometric"
        Q["numpositions"] = 4
        Q["unitID"] = [None] * Q["numpositions"]
        Q["unitID"][0] = "4V9F|1|0|U|19"
        Q["unitID"][1] = "4V9F|1|0|A|524"
        Q["unitID"][2] = "4V9F|1|R|SER|5"
        Q["unitID"][3] = "4V9F|1|R|VAL|6"
        Q["unittype"] = [None] * Q["numpositions"]
        # the following values should be automatically calculated from field 3 of the unit ID
        # this tells FR3D where to find the 3D coordinates of these units
        Q["queryMoleculeType"][0] = "RNA"
        Q["queryMoleculeType"][1] = "RNA"
        Q["queryMoleculeType"][2] = "protein"
        Q["queryMoleculeType"][3] = "protein"
        # the user should be able to decide what type of molecule can appear in a candidate
        # you might want to allow an amino acid to be replaced by a nucleotide; we see that happen sometimes
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["requiredMoleculeType"][2] = ["protein"]
        Q["requiredMoleculeType"][3] = ["protein"]

        Q["discrepancy"] = 0.3
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F|1|0']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F|1|9']   # set of IFEs to search
        Q["searchFiles"] = ['6K0A-local_now']
        Q["searchFiles"] = ['5DC3']   # set of IFEs to search



    elif name == 'stacked bases with aa':
        Q["name"] = name
        Q["type"] = "geometric"
        Q["numpositions"] = 3
        Q["unitID"] = [None] * Q["numpositions"]
        Q["unitID"][0] = "4V9F|1|0|G|1302"
        Q["unitID"][1] = "4V9F|1|0|G|1354"
        Q["unitID"][2] = "4V9F|1|L|LYS|5"
        Q["unittype"] = [None] * Q["numpositions"]
        # the following values should be automatically calculated from field 3 of the unit ID
        # this tells FR3D where to find the 3D coordinates of these units
        Q["queryMoleculeType"][0] = "RNA"
        Q["queryMoleculeType"][1] = "RNA"
        Q["queryMoleculeType"][2] = "protein"
        # the user should be able to decide what type of molecule can appear in a candidate
        # you might want to allow an amino acid to be replaced by a nucleotide; we see that happen
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["requiredMoleculeType"][2] = ["protein"]

        Q["discrepancy"] = 0.4
        Q["searchFiles"] = ['4V9F|1|9']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == "cWW with amino acid in minor groove":
        Q["name"] = name
        Q["type"] = "geometric"
        Q["numpositions"] = 3
        Q["unitID"] = ["4V9F|1|0|A|2089","4V9F|1|0|U|2655","4V9F|1|B|GLN|254"]
        Q["discrepancy"] = 0.3
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][0] = "N"
        Q["interactionMatrix"][1][1] = "N"
        Q["interactionMatrix"][2][2] = "X"
#        Q["interactionMatrix"][0][2] = "A,ARG C,GLU"
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search

    elif name == 'RNA-protein-cWW-minor-groove':
        Q["name"] = name
        Q["type"] = "geometric"
        Q["numpositions"] = 4
        Q["unitID"] = [None] * Q["numpositions"]
        Q["unitID"][0] = "4V9F|1|0|C|2035"
        Q["unitID"][1] = "4V9F|1|0|G|1744"
        Q["unitID"][2] = "4V9F|1|K|HIS|44"
        Q["unitID"][3] = "4V9F|1|K|GLU|13"
        Q["discrepancy"] = 0.5
        # the following values should be automatically calculated from field 3 of the unit ID
        # this tells FR3D where to find the 3D coordinates of these units
        Q["queryMoleculeType"][0] = "RNA"
        Q["queryMoleculeType"][1] = "RNA"
        Q["queryMoleculeType"][2] = "protein"
        Q["queryMoleculeType"][3] = "protein"
        # the user should be able to decide what type of molecule can appear in a candidate
        # you might want to allow an amino acid to be replaced by a nucleotide; we see that happen
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["requiredMoleculeType"][2] = ["protein"]
        Q["requiredMoleculeType"][3] = ["protein"]
        # the user can decide what the base or side-chain should be.  These are generically called UnitType
        Q["requiredUnitType"][0] = []
        Q["requiredUnitType"][1] = []
        Q["requiredUnitType"][2] = ["HIS"]
        Q["requiredUnitType"][3] = []

        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "cWW"

        Q["searchFiles"] = ['4V9F|1|0']   # set of IFEs to search
    #   Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'cWWpairWithAminoAcids':
        Q["name"] = name
        Q["type"] = "mixed"
        Q["numpositions"] = 4
        Q["unitID"] = [None] * Q["numpositions"]
        Q["unitID"][0] = "4V9F|1|0|C|2035"
        Q["unitID"][1] = "4V9F|1|0|G|1744"
        Q["unitID"][2] = "4V9F|1|K|HIS|44"
        Q["unitID"][3] = "4V9F|1|K|GLU|13"
        Q["unittype"] = [None] * Q["numpositions"]
        # the following values should be automatically calculated from field 3 of the unit ID
        # this tells FR3D where to find the 3D coordinates of these units
        Q["queryMoleculeType"][0] = "RNA"
        Q["queryMoleculeType"][1] = "RNA"
        Q["queryMoleculeType"][2] = "protein"
        Q["queryMoleculeType"][3] = "protein"
        # the user should be able to decide what type of molecule can appear in a candidate
        # you might want to allow an amino acid to be replaced by a nucleotide; we see that happen
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["requiredMoleculeType"][2] = ["protein"]
        Q["requiredMoleculeType"][3] = ["protein"]

        Q["discrepancy"] = 0.5
        Q["searchFiles"] = ['4V9F|1|0']   # set of IFEs to search

    elif name == 'SR core':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 5
        Q["searchFiles"] = ['4V9F|1|0']   # set of IFEs to search
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["requiredMoleculeType"][2] = ["RNA"]
        Q["requiredMoleculeType"][3] = ["RNA"]
        Q["requiredMoleculeType"][4] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][1][2] = "cSH"
        Q["interactionMatrix"][2][3] = "tWH"
        Q["interactionMatrix"][0][4] = "tHH"
        Q["interactionMatrix"][0][1] = "s55"

    elif name == 'stacked cWW':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 4
        Q["searchFiles"] = ['4V9F|1|0']   # set of IFEs to search
        Q["requiredMoleculeType"][0] = ["RNA"]
        Q["requiredMoleculeType"][1] = ["RNA"]
        Q["requiredMoleculeType"][2] = ["RNA"]
        Q["requiredMoleculeType"][3] = ["RNA"]
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "s35"
        Q["interactionMatrix"][1][2] = "cWW"
        Q["interactionMatrix"][2][3] = "s35"
        Q["interactionMatrix"][0][3] = "cWW"

    elif name == 'internal loop':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 4
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "bSS"
        Q["interactionMatrix"][1][2] = "cWW"
        Q["interactionMatrix"][2][3] = "bSS"
        Q["interactionMatrix"][0][3] = "cWW"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][3][2] = ">"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search

    elif name == 'GNRA hairpin symbolic':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW"
        Q["interactionMatrix"][1][4] = "tSH GA"
        Q["interactionMatrix"][2][3] = "stack"
        Q["interactionMatrix"][3][4] = "stack"
        Q["interactionMatrix"][1][0] = "=1 after"
        Q["interactionMatrix"][2][1] = "<4 after"
        Q["interactionMatrix"][3][2] = "<4 after"
        Q["interactionMatrix"][4][3] = "<4 after"
        Q["interactionMatrix"][5][4] = "=1 after"
        Q["searchFiles"] = ['4V9F','5J7L']   # set of IFEs to search

    elif name == 'Hairpin flanking pair':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "bSS and cWW"
        Q["interactionMatrix"][1][0] = ">"
        Q["searchFiles"] = ['4V9F|1|9']   # set of IFEs to search

    elif name == 'Hairpin interacts with single strand':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][2] = "cWW and bSS"
        Q["interactionMatrix"][3][5] = "bSS"
        Q["interactionMatrix"][1][4] = "pair stack BPh BR"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = ">"
        Q["searchFiles"] = ['5J7L|1|AA']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F|1|0']   # set of IFEs to search

    elif name == 'base surrounded by bases':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "cSW tSW cSH tSH cSS tSS"
        Q["interactionMatrix"][0][2] = "cWW tWW cWH tWH cWS tWS"
        Q["interactionMatrix"][0][3] = "cHW tHW cHH tHH cHS tHS"
        Q["interactionMatrix"][0][4] = "s33 s35"
        Q["interactionMatrix"][0][5] = "s53 s55"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search

    elif name == 'base surrounded by amino acids':
        Q["name"] = name
        Q["type"] = "geometric"
        Q["discrepancy"] = 0.4
        Q["numpositions"] = 6
        Q["unitID"] = ["4V9F|1|0|G|2471", "4V9F|1|0|A|2633", "4V9F|1|0|C|2114", "4V9F|1|0|U|2278", "4V9F|1|0|C|2472", "4V9F|1|0|A|2470"]
        Q["requiredMoleculeType"] = ["RNA","protein","protein","protein","protein","protein"]
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'base almost surrounded by amino acids':
        Q["name"] = name
        Q["type"] = "geometric"
        Q["discrepancy"] = 0.4
        Q["numpositions"] = 5
        Q["unitID"] = ["4V9F|1|0|G|2471", "4V9F|1|0|C|2114", "4V9F|1|0|U|2278", "4V9F|1|0|C|2472", "4V9F|1|0|A|2470"]
        Q["requiredMoleculeType"] = ["RNA","protein","protein","protein","protein"]
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'junction 10':
        Q["name"] = name
        Q["type"] = "geometric"
        Q["discrepancy"] = 0.15
        Q["unitID"] = ["5WE4|1|a|A|935","5WE4|1|a|C|936","5WE4|1|a|A|937","5WE4|1|a|A|938","5WE4|1|a|G|939","5WE4|1|a|G|1343","5WE4|1|a|C|1344","5WE4|1|a|U|1345","5WE4|1|a|A|1346","5WE4|1|a|G|1347","5WE4|1|a|U|1348","5WE4|1|a|A|1349","5WE4|1|a|A|1350","5WE4|1|a|U|1351","5WE4|1|a|G|1371","5WE4|1|a|U|1372","5WE4|1|a|G|1373","5WE4|1|a|A|1374","5WE4|1|a|A|1375","5WE4|1|a|U|1376","5WE4|1|a|A|1377","5WE4|1|a|C|1378","5WE4|1|a|G|1379","5WE4|1|a|U|1380"]
        Q["numpositions"] = len(Q["unitID"])
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["searchFiles"] = ['5WE4|1|a']   # set of IFEs to search
        Q["searchFiles"] = ['5MDY|1|2']   # set of IFEs to search

    elif name == "kink turn 65553":
        Q["name"] = name
        Q["type"] = "geometric"
        Q["unitID"] = ["4V9F|1|0|U|1026","4V9F|1|0|A|1032","4V9F|1|0|G|1034","4V9F|1|0|C|936","4V9F|1|0|G|940"]
        Q["numpositions"] = len(Q["unitID"])
        Q["discrepancy"] = 0.6                   # maximum discrepancy cutoff
        # add symbolic constraints to greatly speed up the query
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][2] = "bSS"   # border single strand
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][2][3] = "cWW"   #
        Q["interactionMatrix"][0][4] = "cWW"
        Q["interactionMatrix"][1][0] = ">"     # chain direction constraint
        Q["interactionMatrix"][2][1] = ">"     # O(n^2) constraint; implies [2][0] is >
        Q["interactionMatrix"][4][3] = ">"     # no [3][2] constraint, different strands
        Q["searchFiles"] = ['4V9F']   # list of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'UGGU tandem':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][7] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][6] = "UG"
        Q["interactionMatrix"][2][5] = "GU"
        Q["interactionMatrix"][1][0] = "<=1 >"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][3][2] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["interactionMatrix"][6][5] = "<=1 >"
        Q["interactionMatrix"][7][6] = "<=1 >"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GUUG tandem':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][7] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][6] = "GU"
        Q["interactionMatrix"][2][5] = "UG"
        Q["interactionMatrix"][1][0] = "<=1 >"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][3][2] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["interactionMatrix"][6][5] = "<=1 >"
        Q["interactionMatrix"][7][6] = "<=1 >"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GGUC_GGUC tandem':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][7] = "cWW GC"
        Q["interactionMatrix"][3][4] = "cWW CG"
        Q["interactionMatrix"][1][6] = "GU"
        Q["interactionMatrix"][2][5] = "UG"
        Q["interactionMatrix"][1][0] = "<=1 >"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][3][2] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["interactionMatrix"][6][5] = "<=1 >"
        Q["interactionMatrix"][7][6] = "<=1 >"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU tandem at end of helix':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][7] = "cWW AU UA CG GC"
        Q["interactionMatrix"][3][4] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][6] = "UG"
        Q["interactionMatrix"][2][5] = "GU"
        Q["interactionMatrix"][0][1] = "bSS"
        Q["interactionMatrix"][6][7] = "bSS"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][3][2] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["interactionMatrix"][6][5] = "<=1 >"
        Q["interactionMatrix"][7][6] = ">"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU tandem at end of helix adjacent above':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][7] = "cWW AU UA CG GC"
        Q["interactionMatrix"][3][4] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][6] = "UG"
        Q["interactionMatrix"][2][5] = "GU"
        Q["interactionMatrix"][0][1] = "~bSS"
        Q["interactionMatrix"][6][7] = "bSS"
        Q["interactionMatrix"][1][0] = "> =1"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][3][2] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["interactionMatrix"][6][5] = "<=1 >"
        Q["interactionMatrix"][7][6] = ">"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU tandem at end of helix adjacent below':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][7] = "cWW AU UA CG GC"
        Q["interactionMatrix"][3][4] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][6] = "UG"
        Q["interactionMatrix"][2][5] = "GU"
        Q["interactionMatrix"][0][1] = "bSS"
        Q["interactionMatrix"][6][7] = "~bSS"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][3][2] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["interactionMatrix"][6][5] = "<=1 >"
        Q["interactionMatrix"][7][6] = "<=1 >"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU tandem at end of helix one unpaired':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][7] = "cWW AU UA CG GC"
        Q["interactionMatrix"][3][4] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][6] = "UG"
        Q["interactionMatrix"][2][5] = "GU"
        Q["interactionMatrix"][0][2] = "bSS"
        Q["interactionMatrix"][5][7] = "bSS"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][3][2] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["interactionMatrix"][6][5] = "<=1 >"
        Q["interactionMatrix"][7][6] = ">"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU tandem at end of helix two unpaired':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][7] = "cWW AU UA CG GC"
        Q["interactionMatrix"][3][4] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][6] = "UG"
        Q["interactionMatrix"][2][5] = "GU"
        Q["interactionMatrix"][0][3] = "bSS"
        Q["interactionMatrix"][4][7] = "bSS"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][3][2] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["interactionMatrix"][6][5] = "<=1 >"
        Q["interactionMatrix"][7][6] = ">"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search


    elif name == 'UG tandem at end of helix':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][7] = "cWW AU UA CG GC"
        Q["interactionMatrix"][3][4] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][6] = "GU"
        Q["interactionMatrix"][2][5] = "UG"
        Q["interactionMatrix"][0][1] = "bSS"
        Q["interactionMatrix"][6][7] = "bSS"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][3][2] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["interactionMatrix"][6][5] = "<=1 >"
        Q["interactionMatrix"][7][6] = ">"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'UG tandem at end of helix adjacent above':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][7] = "cWW AU UA CG GC"
        Q["interactionMatrix"][3][4] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][6] = "GU"
        Q["interactionMatrix"][2][5] = "UG"
        Q["interactionMatrix"][0][1] = "~bSS"
        Q["interactionMatrix"][6][7] = "bSS"
        Q["interactionMatrix"][1][0] = "> =1"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][3][2] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["interactionMatrix"][6][5] = "<=1 >"
        Q["interactionMatrix"][7][6] = ">"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'UG tandem at end of helix adjacent below':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][7] = "cWW AU UA CG GC"
        Q["interactionMatrix"][3][4] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][6] = "GU"
        Q["interactionMatrix"][2][5] = "UG"
        Q["interactionMatrix"][0][1] = "bSS"
        Q["interactionMatrix"][6][7] = "~bSS"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][3][2] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["interactionMatrix"][6][5] = "<=1 >"
        Q["interactionMatrix"][7][6] = "<=1 >"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'UG tandem at end of helix one unpaired':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][7] = "cWW AU UA CG GC"
        Q["interactionMatrix"][3][4] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][6] = "GU"
        Q["interactionMatrix"][2][5] = "UG"
        Q["interactionMatrix"][0][2] = "bSS"
        Q["interactionMatrix"][5][7] = "bSS"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][3][2] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["interactionMatrix"][6][5] = "<=1 >"
        Q["interactionMatrix"][7][6] = ">"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'UG tandem at end of helix two unpaired':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 8
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][7] = "cWW AU UA CG GC"
        Q["interactionMatrix"][3][4] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][6] = "GU"
        Q["interactionMatrix"][2][5] = "UG"
        Q["interactionMatrix"][0][3] = "bSS"
        Q["interactionMatrix"][4][7] = "bSS"
        Q["interactionMatrix"][1][0] = ">"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][3][2] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["interactionMatrix"][6][5] = "<=1 >"
        Q["interactionMatrix"][7][6] = ">"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU Case 0':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
#        Q["interactionMatrix"][1][4] = "GU"
        Q["interactionMatrix"][1][1] = "G"
        Q["interactionMatrix"][4][4] = "U"
        Q["interactionMatrix"][1][0] = "<=1 >"
        Q["interactionMatrix"][2][1] = "<=1 >"
        Q["interactionMatrix"][4][3] = "<=1 >"
        Q["interactionMatrix"][5][4] = "<=1 >"
        Q["searchFiles"] = ['4Y4O|1|1a']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU Case 1 UG cWW':   # in some structures bSS does not allow for positions 1,2 or 3,4 to be adjacent
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][1][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
#        Q["interactionMatrix"][1][4] = "UG"
        Q["interactionMatrix"][1][1] = "U"
        Q["interactionMatrix"][4][4] = "G"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU Case 1b UG cWW':   # make positions 1,2 adjacent and specifically not bSS
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][2] = "~bSS"
        Q["interactionMatrix"][1][4] = "UG"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = "=1 >"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU Case 1c UG cWW':   # make positions 3,4 adjacent and specifically not bSS
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][2] = "bSS"
        Q["interactionMatrix"][1][4] = "UG"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "~bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = "=1 >"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU Case 1 UG interior':   # does not allow for positions 1,2 or 3,4 to be adjacent in some structures?; need separate queries
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "UG"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][5] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU Case 2 GU cWW':   # in some structures bSS does not allow for positions 1,2 or 3,4 to be adjacent
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][1][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "GU"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU Case 2b GU cWW':   # make positions 1,2 adjacent and specifically not bSS
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][2] = "~bSS"
        Q["interactionMatrix"][1][4] = "GU"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = "=1 >"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU Case 2c GU cWW':   # make positions 3,4 adjacent and specifically not bSS
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][2] = "bSS"
        Q["interactionMatrix"][1][4] = "GU"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "~bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = "=1 >"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU Case 2 GU interior':   # does not allow for positions 1,2 or 3,4 to be adjacent in some structures?; need separate queries
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "GU"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][5] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GC Case 0':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "GC"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = "=1 >"
        Q["interactionMatrix"][4][3] = "=1 >"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GC Case 1 CG cWW':   # in some structures bSS does not allow for positions 1,2 or 3,4 to be adjacent
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][1][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "CG"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GC Case 1b CG cWW':   # make positions 1,2 adjacent and specifically not bSS
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][2] = "~bSS"
        Q["interactionMatrix"][1][4] = "CG"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = "=1 >"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GC Case 1c CG cWW':   # make positions 3,4 adjacent and specifically not bSS
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][2] = "bSS"
        Q["interactionMatrix"][1][4] = "CG"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "~bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = "=1 >"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GC Case 1 CG interior':   # does not allow for positions 1,2 or 3,4 to be adjacent in some structures?; need separate queries
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "CG"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][5] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GC Case 2 GC cWW':   # in some structures bSS does not allow for positions 1,2 or 3,4 to be adjacent
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][1][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "GC"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GC Case 2b GC cWW':   # make positions 1,2 adjacent and specifically not bSS
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][2] = "~bSS"
        Q["interactionMatrix"][1][4] = "GC"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = "=1 >"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GC Case 2c GC cWW':   # make positions 3,4 adjacent and specifically not bSS
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][2] = "bSS"
        Q["interactionMatrix"][1][4] = "GC"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "~bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = "=1 >"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GC Case 2 GC interior':   # does not allow for positions 1,2 or 3,4 to be adjacent in some structures?; need separate queries
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "GC"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][5] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4Y1J']   # set of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'AU Case 0':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "AU"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = "=1 >"
        Q["interactionMatrix"][4][3] = "=1 >"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'AU Case 1 UA cWW':   # in some structures bSS does not allow for positions 1,2 or 3,4 to be adjacent
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][1][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "UA"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'AU Case 1b UA cWW':   # make positions 1,2 adjacent and specifically not bSS
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][2] = "~bSS"
        Q["interactionMatrix"][1][4] = "UA"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = "=1 >"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'AU Case 1c UA cWW':   # make positions 3,4 adjacent and specifically not bSS
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][2] = "bSS"
        Q["interactionMatrix"][1][4] = "UA"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "~bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = "=1 >"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'AU Case 1 UA interior':   # does not allow for positions 1,2 or 3,4 to be adjacent in some structures?; need separate queries
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "UA"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][5] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'AU Case 2 AU cWW':   # in some structures bSS does not allow for positions 1,2 or 3,4 to be adjacent
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][1][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "AU"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'AU Case 2b AU cWW':   # make positions 1,2 adjacent and specifically not bSS
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][2] = "~bSS"
        Q["interactionMatrix"][1][4] = "AU"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = "=1 >"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'AU Case 2c AU cWW':   # make positions 3,4 adjacent and specifically not bSS
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][2] = "bSS"
        Q["interactionMatrix"][1][4] = "AU"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "~bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = "=1 >"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'AU Case 2 AU interior':   # does not allow for positions 1,2 or 3,4 to be adjacent in some structures?; need separate queries
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "AU"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][5] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'UU Case 0':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "UU"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = "=1 >"
        Q["interactionMatrix"][4][1] = ">"          # break the UU symmetry
        Q["interactionMatrix"][4][3] = "=1 >"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'UU Case 2 UU cWW':   # in some structures bSS does not allow for positions 1,2 or 3,4 to be adjacent
        # might not find any candidates because UU is not considered to be the end of a helix
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][1][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "UU"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][4] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'UU Case 2 UU interior':   # does not allow for positions 1,2 or 3,4 to be adjacent in some structures?; need separate queries
        # since UU is the same as UU, there is no Case 1 and Case 2, it's all the same
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][2] = "bSS"
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][4] = "UU"
        Q["interactionMatrix"][2][3] = "cWW AU UA GC CG"
        Q["interactionMatrix"][3][5] = "bSS"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][4][3] = ">"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'GU 2 into junction':   # does not allow for positions 1,2 or 3,4 to be adjacent; need separate queries
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 6
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][5] = "cWW AU UA GC CG"  # last cWW before a junction
        Q["interactionMatrix"][1][4] = "GU"
        Q["interactionMatrix"][2][3] = "~cWW"
        Q["interactionMatrix"][0][2] = "bSS"              # single strand on this side
        Q["interactionMatrix"][3][5] = "bSS"              # single strand on this side
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = "=1 >"
        Q["interactionMatrix"][4][3] = "=1 >"
        Q["interactionMatrix"][5][4] = "=1 >"
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search

    elif name == 'AU hairpin Case 1 UA cWW':   # in some structures bSS does not allow for positions 1,2 or 3,4 to be adjacent
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "bSS cWW UA"
        Q["interactionMatrix"][1][0] = ">"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'AU hairpin Case 1 UA after cWW':   # in some structures bSS does not allow for positions 1,2 or 3,4 to be adjacent
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 4
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][3] = "bSS cWW"
        Q["interactionMatrix"][1][2] = "UA"
        Q["interactionMatrix"][1][0] = "=1 >"
        Q["interactionMatrix"][2][1] = ">"
        Q["interactionMatrix"][3][2] = "=1 >"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.48/3.0A/csv']   # list of IFEs to search

    elif name == 'bSS and cWW':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 3
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "cWW AU UA GC CG"
        Q["interactionMatrix"][1][2] = "bSS"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search

    elif name == 'unary no pair':
        Q["name"] = name
        Q["type"] = "symbolic"
        Q["numpositions"] = 2
        Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
        Q["interactionMatrix"][0][1] = "stack"
        Q["interactionMatrix"][1][0] = "> =1"
        Q["interactionMatrix"][0][0] = "~pair N"
        Q["interactionMatrix"][1][1] = "cWW"
        Q["searchFiles"] = ['4V9F']   # set of IFEs to search
        Q["searchFiles"] = ['6ZMI']   # set of IFEs to search

    else:
        print("Query name %s not recognized by query_definitions.py" % name)

#    with open('../pythoncode/jsonfiles/'+ name+'.json','w') as outfile:
#        json.dump(Q,outfile)

    return Q
