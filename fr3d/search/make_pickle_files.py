import numpy as np
import os.path
from collections import defaultdict
from time import time
import pickle
import urllib
import sys
from fr3d_configuration import CIFPATH
from fr3d_configuration import DATAPATH
from fr3d_configuration import SERVER

# import the version of urlretrieve appropriate to the Python version
if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
else:
    from urllib.request import urlretrieve as urlretrieve



def readPDBDatafile():
    """
    Read .pickle file containing data about each nucleic-
    acid-containing PDB file.
    """

    datafile = {}

    filename = "NA_datafile.pickle"
    pathAndFileName = os.path.join(DATAPATH,'units',filename)

    #if not os.path.exists(pathAndFileName) and not SERVER:
    try:
        if not SERVER:
            urlretrieve("http://rna.bgsu.edu/units/" + filename, pathAndFileName)
            print("Downloaded %s because it changes every week" % filename)
    except:
        print("Unable to download %s" % filename)

    if os.path.exists(pathAndFileName):
        try:
            if sys.version_info[0] < 3:
                datafile = pickle.load(open(pathAndFileName,"rb"))
            else:
                datafile = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
        except:
            print("Could not read "+filename)
            Q["userMessage"].append("Could not retrieve PDB datafile "+filename)

    return datafile

def readNAPositionsFile(Q, chainString, starting_index):
    """
    Read .pickle file of RNA base center and rotation matrix;
    download if necessary
    """

    ids = []
    chainIndices = []
    centers = []
    rotations = []
    line_num = starting_index

    centers = np.zeros((0,3))
    rotations = np.zeros((0,3,3))
    id_to_index = defaultdict()
    index_to_id = defaultdict()

    filename = chainString + "_NA_base_rotation.pickle"
    filename = chainString + "_RNA" + '.pickle'
    pathAndFileName = os.path.join(DATAPATH,'units',filename)

    if not os.path.exists(pathAndFileName) and not SERVER:
        try:
            urlretrieve("http://rna.bgsu.edu/units/" + filename, pathAndFileName) # testing
            print("Downloaded "+filename)
        except:
            print("Unable to download %s" % filename)

    if os.path.exists(pathAndFileName):
        try:
            if sys.version_info[0] < 3:
                ids, chainIndices, centers, rotations = pickle.load(open(pathAndFileName,"rb"))
            else:
                ids, chainIndices, centers, rotations = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
        except:
            print("Could not read "+filename)
            Q["errorMessage"].append("Could not retrieve RNA unit file "+filename)

        for i in range(0,len(ids)):
            id_to_index[ids[i]] = line_num
            index_to_id[line_num] = ids[i]
            line_num += 1

    else:
        print("Could not find "+filename)
        Q["errorMessage"].append("Could not retrieve RNA unit file "+filename)

    centers = np.asarray(centers)
    rotations = np.asarray(rotations)

    return Q, centers, rotations, ids, id_to_index, index_to_id, chainIndices

def readNABackboneFile(Q, chainString, starting_index):
    """
    Read .pickle file of RNA phosphate and sugar centers; download if necessary
    """

    ids = []
    chainIndices = []
    centers = []
    rotations = []
    line_num = starting_index

    phosphate = np.zeros((0,3))
    sugar = np.zeros((0,3))
    id_to_index = defaultdict()
    index_to_id = defaultdict()

    filename = chainString + "_NA_phosphate_sugar" + '.pickle'
    filename = chainString + "_protein" + '.pickle'
    pathAndFileName = os.path.join(DATAPATH,'units',filename)

    if not os.path.exists(pathAndFileName) and not SERVER:
        urlretrieve("http://rna.bgsu.edu/units/"+filename, pathAndFileName)
        print("Downloaded "+filename)

    if os.path.exists(pathAndFileName):
        if sys.version_info[0] < 3:
            try:
                ids, chainIndices, phosphate, sugar = pickle.load(open(pathAndFileName,"rb"))
            except:
                print("Could not read "+filename)
                Q["userMessage"].append("Could not retrieve RNA unit file "+filename)
        else:
            try:
                ids, chainIndices, phosphate, sugar = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
            except:
                print("Could not read "+filename)
                Q["userMessage"].append("Could not retrieve RNA unit file "+filename)

        for i in range(0,len(ids)):
            id_to_index[ids[i]] = line_num
            index_to_id[line_num] = ids[i]
            line_num += 1

    else:
        print("Could not find "+filename)
        Q["userMessage"].append("Could not retrieve RNA unit file "+filename)

    phosphate = np.asarray(phosphate)
    sugar     = np.asarray(sugar)

    return Q, phosphate, sugar, ids, id_to_index, index_to_id, chainIndices

def readProteinPositionsFile(Q, PDBID, starting_index):

    ids = []
    chainIndices = []
    centers = []
    rotations = []
    line_num = starting_index

    centers = np.zeros((0,3))  # of lines/nucleotides in file
    rotations = np.zeros((0,3,3))
    id_to_index = defaultdict()
    index_to_id = defaultdict()

    filename = PDBID + "_protein" + '.pickle'
    pathAndFileName = os.path.join(DATAPATH,'units',filename)

    if not os.path.exists(pathAndFileName) and not SERVER:
        urlretrieve("http://rna.bgsu.edu/units/"+filename, pathAndFileName)
        print("Downloaded "+filename)

    if os.path.exists(pathAndFileName):
        if sys.version_info[0] < 3:
            try:
                ids, chainIndices, centers, rotations = pickle.load(open(pathAndFileName,"rb"))
            except:
                print("Could not read "+filename+" C=================================")
                Q["userMessage"].append("Could not retrieve protein unit file "+filename)
        else:
            try:
                ids, chainIndices, centers, rotations = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
            except:
                print("Could not read "+filename+" C*********************************")
                Q["userMessage"].append("Could not retrieve protein unit file "+filename)

        for i in range(0,len(ids)):
            id_to_index[ids[i]] = line_num
            index_to_id[line_num] = ids[i]
            line_num += 1

    else:
        print("Could not find "+filename)
        Q["userMessage"].append("Could not retrieve protein unit file "+filename)

    centers = np.asarray(centers)
    rotations = np.asarray(rotations)

    return Q, centers, rotations, ids, id_to_index, index_to_id, chainIndices

def readNAPairsFile(Q, PDBID, id_to_index, alternate = ""):
    """
    We need to read the pairwise interactions whether or not there are
    constraints, to be able to display the interactions in the output.
    alternate is for new and/or alternate lists of interacting pairs.
    """

    interactionToPairs = defaultdict(list)
    interactionToIndexPairs = {}
    pairToInteractions = defaultdict(list)
    pairToCrossingNumber = {}

    pairsFileName = PDBID + '_RNA_pairs' + alternate + '.pickle'
    pathAndFileName = os.path.join(DATAPATH,'pairs',pairsFileName)

    if not os.path.exists(pathAndFileName) and not SERVER:
        urlretrieve("http://rna.bgsu.edu/pairs/"+pairsFileName, pathAndFileName) # testing
        print("Downloaded "+pairsFileName)

    # note:  if the file is not present on the server, a text file with a 404 error will be downloaded
    # so it will look like a file was downloaded, but it's not the file you need!

    if os.path.exists(pathAndFileName):
        if sys.version_info[0] < 3:
            try:
                interactionToPairs = pickle.load(open(pathAndFileName,"rb"))
            except:
                print("Could not read "+pairsFileName+", it may not be available D=================================")
                Q["userMessage"].append("Could not retrieve RNA pairs file "+pairsFileName)
        else:
            try:
                interactionToPairs = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
            except:
                print("Could not read "+pairsFileName+", it may not be available D*********************************")
                Q["userMessage"].append("Could not retrieve RNA pairs file "+pairsFileName)

    # interactionToPairs has this structure:
    # interactionToPairs['cWW'] is a list of triples, each triple being (unit_id_1,unit_id_2,range)
    # However, what is passed back is interactionToIndexPairs, which has a different structure
    # interactionToIndexPairs['cWW'][0] is the list of pairs
    # interactionToIndexPairs['cWW'][1] is the list of crossing numbers

    # convert pairs of unit ids to pairs of indices for use in FR3D.py
    for interaction in interactionToPairs:
        if interaction in Q["activeInteractions"]:
            newListOfPairs = []
            crossingNumber = []
            if "BPh" in interaction or "BR" in interaction:
                newListOfPairsSelf = []                      # 0BPh and maybe others can make self interactions
                crossingNumberSelf = []
                for pair in interactionToPairs[interaction]:
                    if pair[0] in id_to_index and pair[1] in id_to_index:
                        indexPair = (id_to_index[pair[0]],id_to_index[pair[1]])
                        pairToInteractions[indexPair].append(interaction)
                        pairToCrossingNumber[indexPair] = pair[2]
                        if pair[0] == pair[1]:
                            newListOfPairsSelf.append(indexPair)
                            crossingNumberSelf.append(pair[2])
                        else:
                            newListOfPairs.append(indexPair)
                            crossingNumber.append(pair[2])
                interactionToIndexPairs[interaction] = (newListOfPairs,crossingNumber)
                interactionToIndexPairs[interaction+"_self"] = (newListOfPairsSelf,crossingNumberSelf)
            else:
                for pair in interactionToPairs[interaction]:
                    if pair[0] in id_to_index and pair[1] in id_to_index:
                        indexPair = (id_to_index[pair[0]],id_to_index[pair[1]])
                        pairToInteractions[indexPair].append(interaction)
                        pairToCrossingNumber[indexPair] = pair[2]
                        newListOfPairs.append(indexPair)
                        crossingNumber.append(pair[2])
                interactionToIndexPairs[interaction] = (newListOfPairs,crossingNumber)
        else:
            for pair in interactionToPairs[interaction]:
                if pair[0] in id_to_index and pair[1] in id_to_index:
                    indexPair = (id_to_index[pair[0]],id_to_index[pair[1]])
                    pairToInteractions[indexPair].append(interaction)
                    pairToCrossingNumber[indexPair] = pair[2]

    return Q, interactionToIndexPairs, pairToInteractions, pairToCrossingNumber


def readUnitAnnotations(Q, ifename):
    # split up ifename and read unit annotations for each individual chain
    # unit annotations are a dictionary mapping unit id to a tuple t
    # t[0] is the chi angle in degrees
    # t[1] is anti/syn/intermediate syn
    # When more annotations are added, they can be added to the tuple

    chains =  ifename.split('+')

    unit_id_to_annotation = {}

    for chain in chains:
        filename        = "%s_NA_unit_annotations.pickle" % chain.replace("|","-")
        pathAndFileName = os.path.join(DATAPATH,'units',filename)

        if not os.path.exists(pathAndFileName) and not SERVER:
            urlretrieve("http://rna.bgsu.edu/units/" + filename, pathAndFileName)
            print("Downloaded "+filename)

        # note:  if the file is not present on the server, a text file with a 404 error will be downloaded
        # so it will look like a file was downloaded, but it's not the file you need!

        if os.path.exists(pathAndFileName):
            if sys.version_info[0] < 3:
                try:
                    new_mapping = pickle.load(open(pathAndFileName,"rb"))
                    unit_id_to_annotation.update(new_mapping)
                except:
                    print("Could not read "+filename+", it may not be available D=================================")
                    Q["userMessage"].append("Could not retrieve RNA pairs file "+filename)
                    try:
                        os.remove(pathAndFileName)
                    except:
                        Q["userMessage"].append("Could not remove RNA pairs file "+filename)
            else:
                try:
                    new_mapping = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
                    unit_id_to_annotation.update(new_mapping)
                except:
                    print("Could not read "+filename+", it may not be available D*********************************")
                    Q["userMessage"].append("Could not retrieve RNA pairs file "+filename)
                    try:
                        os.remove(pathAndFileName)
                    except:
                        Q["userMessage"].append("Could not remove RNA pairs file "+filename)


    return Q, unit_id_to_annotation
