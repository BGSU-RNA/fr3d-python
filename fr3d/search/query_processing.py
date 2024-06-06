"""
Interpret the constraints in the query.
For geometric and mixed searches, set up limits on pairwise distances.
"""

from collections import defaultdict
import json
import math
import numpy as np
import os
import pickle
import sys
import urllib

from file_reading import readNAPositionsFile
from file_reading import readProteinPositionsFile
from file_reading import readPDBDatafile
from file_reading import get_DATAPATHUNITS
from file_reading import get_CIFPATH

from fr3d_configuration import SERVER

from fr3d.data.mapping import modified_base_atom_list,parent_atom_to_modified,modified_atom_to_parent,modified_base_to_parent

# identify codes that go with each type of molecule
RNA_unit_types = set(["A","C","G","U"])
DNA_unit_types = set(["DA","DC","DG","DT"])

RNA_modified_list = set()
DNA_modified_list = set()

for modified, parent in modified_base_to_parent.items():
    if parent in RNA_unit_types:
        RNA_modified_list.add(modified)
    else:
        DNA_modified_list.add(modified)

protein_unit_types = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","PYL","SER","SEC","THR","TRP","TYR","VAL","ASX","GLX","XAA","XLE"]
# it would be good to have a list of modified amino acids


def getMoleculeType(unitType):
    """
    Determine molecule type from unit ID by identifying unit type specifiers
    """

    if "|" in unitType:
        fields = unitType.split('|')
        unitType = fields[3]

    unitType = unitType.upper()

    # if sys.version_info[0] < 3:
    #     unitType = string.upper(unitType)    # convert to uppercase
    # else:
    #     unitType = unitType.upper()

    if unitType in RNA_unit_types:
        return "RNA"
    elif unitType in DNA_unit_types:
        return "DNA"
    elif unitType in protein_unit_types:
        return "protein"
    elif unitType in RNA_modified_list:
        return "RNA"
    elif unitType in RNA_modified_list:
        return "DNA"
    else:
        return ""

# synonyms and abbreviations for pairwise constraints
synonym = {}
synonym["pair"] = ['cWW','acWW','cWWa','tWW','tWWa','cHH','tHH','cSS','tSS','cWH','cHW','tWH','tHW','cWS','cSW','tWS','tSW','cHS','cSH','tHS','tSH']
synonym["cis"] = ['cWW','acWW','cWWa','cHH','cSS','cWH','cHW','cWS','cSW','cHS','cSH']
synonym["trans"] = ['tWW','tHH','tSS','tWH','tHW','tWS','tSW','tHS','tSH']
synonym["stack"] = ['s35','s53','s33','s55']
synonym["BPh"] = ['0BPh','1BPh','2BPh','3BPh','4BPh','5BPh','6BPh','7BPh','8BPh','9BPh']
synonym["BR"] = ['0BR','1BR','2BR','3BR','4BR','5BR','6BR','7BR','8BR','9BR']
synonym["coplanar"] = ['cp']
synonym["sugar_ribose"] = ['cSR','cRS','tSR','tRS']

synonym["s3O"] = ["s3O2'","s3O3'","s3O4'","s3O5'","s3OP1","s3OP2"]
synonym["s5O"] = ["s5O2'","s5O3'","s5O4'","s5O5'","s5OP1","s5OP2"]
synonym["sO3"] = ["sO2'3","sO3'3","sO4'3","sO5'3","sOP13","sOP23"]
synonym["sO5"] = ["sO2'5","sO3'5","sO4'5","sO5'5","sOP15","sOP25"]
synonym["sO2'"] = ["sO2'3","sO2'5"]
synonym["sO3'"] = ["sO3'3","sO3'5"]
synonym["sO4'"] = ["sO4'3","sO4'5"]
synonym["sO5'"] = ["sO5'3","sO5'5"]
synonym["sOP1"] = ["sOP13","sOP15"]
synonym["sOP2"] = ["sOP23","sOP25"]
synonym["sOP"] = ["sOP13","sOP23","sOP15","sOP25"]
synonym["sO"]  = synonym["sO3"] + synonym["sO5"]

# list of interactions, to facilitate reversing them
sOF_interactions = ["s3O","s5O","s3O2'","s3O3'","s3O4'","s3O5'","s3OP1","s3OP2","s5O2'","s5O3'","s5O4'","s5O5'","s5OP1","s5OP2"]
sFO_interactions = ["sO3","sO5","sO2'3","sO3'3","sO4'3","sO5'3","sOP13","sOP23","sO2'5","sO3'5","sO4'5","sO5'5","sOP15","sOP25"]

map_sFO_to_sOF = {}
for i in range(0,len(sOF_interactions)):
    map_sFO_to_sOF[sFO_interactions[i]] = sOF_interactions[i]
    map_sFO_to_sOF["n"+sFO_interactions[i]] = "n"+sOF_interactions[i]

# the ones above can also start with "n" for "near"
for k in list(synonym.keys()):
    synonym["n"+k] = ["n"+s for s in synonym[k]]

synonym["borderSS"] = ["bSS"]
synonym["flankSS"] = ["bSS"]
synonym["flank"] = ["bSS"]

synonym["IS"] = ["int_syn"]   # anti, syn, intermediate syn are possible bond annotations

synonym["AND"] = ["and"]
synonym["&"] = ["and"]
synonym["&&"] = ["and"]

synonym["before"] = "<"
synonym["after"] = ">"

# synonyms for crossing number constraints
crossing = {}
crossing["nested"] = (0,0)
crossing["local"] = (0,2)
crossing["LR"] = (3,float("inf"))

# collect together all known interaction constraints to be able to control what is interpreted as a constraint
# basically, in order to be listed here, there needs to be a synonym for it ... strange organization!
allInteractionConstraints = []
for key in synonym:
    allInteractionConstraints.extend(synonym[key])
allInteractionConstraints = set(allInteractionConstraints)


def emptyInteractionList(n):
    """
    set up an empty list of required interactions,
    then the user only needs to set the actual constraints needed
    """

    emptyList = defaultdict(dict)
    for i in range(0,n):
        for j in range(0,n):
            emptyList[i][j] = []
    return emptyList

def emptyInteractionMatrix(n):
    """
    set up an empty list of required interactions,
    then the user only needs to set the actual constraints needed
    """

    emptyList = defaultdict(dict)
    for i in range(0,n):
        for j in range(0,n):
            emptyList[i][j] = ""
    return emptyList


def readUnitFileNames(Q,PDB_data_file):
    """
    Read all filenames in DATAPATHUNITS to know which chains are available,
    that are not downloaded from PDB.
    A typical filename is 4V9F-1-0_NA.pickle or user-data_file-1-A_NA.pickle
    """

    import glob

    Q, DATAPATHUNITS = get_DATAPATHUNITS(Q)

    units_path = os.path.join(DATAPATHUNITS,'*')
    allfiles = glob.glob(units_path)
    print("  query_processing: There are %d files in %s" % (len(allfiles),units_path))

    # mimic the structure of PDB_data_file; map file_id to list of chains
    unit_data_file = {}

    for file in allfiles:
        # remove path
        chain_filename = os.path.split(file)[1]

        # split by '_' character to identify molecule type
        fields = chain_filename.split('_')

        if len(fields) < 2:
            continue
        else:
            if fields[-1] == 'RNA.pickle':
                molecule_type = 'NA'
                chain_id = chain_filename.replace('_RNA.pickle','')
            elif fields[-1] == 'DNA.pickle':
                molecule_type = 'NA'
                chain_id = chain_filename.replace('_NA.pickle','')
            elif fields[-1] == 'NA.pickle':
                molecule_type = 'NA'
                chain_id = chain_filename.replace('_NA.pickle','')
            elif fields[-1] == 'protein.pickle':
                molecule_type = 'protein'
                chain_id = chain_filename.replace('_protein.pickle','')
            else:
                continue

            chain_fields = chain_id.split('-')

            if len(chain_fields) < 3:
                continue
            else:
                model = chain_fields[-2]
                chain = chain_fields[-1]
                file_id = '-'.join(chain_fields[0:-2])

                # print(chain_filename)
                # print(chain_id)
                # print(model)
                # print(chain)
                # print(file_id)

                # skip files that are already in PDB_data_file
                if not file_id in PDB_data_file:
                    if not file_id in unit_data_file:
                        unit_data_file[file_id] = {}
                        unit_data_file[file_id]['chains'] = {}
                        unit_data_file[file_id]['chains']['NA'] = []
                        unit_data_file[file_id]['chains']['protein'] = []
                        unit_data_file[file_id]['model'] = set()

                    unit_data_file[file_id]['chains'][molecule_type].append(chain)
                    unit_data_file[file_id]['model'].add(model)

                    # print('  query_processing: readUnitFileNames: chain_filename %s has file_id %s and molecule_type %s' % (chain_filename,file_id,molecule_type))

    return unit_data_file


def readQueryFromJSON(JSONfilename):
    """
    Read query specification in JSON file and interpret as a query specification
    """

    JSONfilename = JSONfilename.replace('"','')
    JSONfilename = JSONfilename.replace("'","")

    filename = None

    if os.path.exists(JSONfilename):
        # when given a direct reference to the JSON file, use that
        filename = JSONfilename

    if not filename:
        # look in the usual places for the file
        try:
            from fr3d_configuration import JSONPATH
            pathAndFileName = os.path.join(JSONPATH,JSONfilename)
            if os.path.exists(pathAndFileName):
                filename = pathAndFileName
        except:
            pass

    if not filename:
        pathAndFileName = os.path.join('queries',JSONfilename)
        if os.path.exists(pathAndFileName):
            filename = pathAndFileName

    if not filename:
        # name could be from a WebFR3D search, so try to download it
        if JSONfilename.startswith("Query_"):

            try:
                if not os.path.exists(JSONPATH):
                    os.makedirs(JSONPATH)
            except:
                print("Error: Could not create directory " + JSONPATH + " to store query files")
                Q = {}
                Q["errorMessage"] = []
                Q["errorMessage"].append("Error: Could not create directory " + JSONPATH + " to store query files")
                return Q

            try:
                id = JSONfilename.replace("Query_","").replace(".json","")

                queryURL = "http://rna.bgsu.edu/webfr3d/Results/" + id + "/" + JSONfilename
                print("Downloading %s from %s to %s" % (JSONfilename,queryURL,JSONPATH))

                pathAndFileName = os.path.join(JSONPATH,JSONfilename)

                if sys.version_info[0] < 3:
                    urllib.urlretrieve(queryURL, pathAndFileName)  # python 2
                else:
                    urllib.request.urlretrieve(queryURL, pathAndFileName)  # python 3

                filename = JSONPATH + JSONfilename
            except:
                print("Error: Could not find or download query file " + JSONfilename)

    # read the json file
    if filename:
        try:
            with open(filename) as json_file:
                lines = json_file.readlines()

            # strip out comments from each line, starting at # and going to the end of the line
            cleaned_lines = []
            for line in lines:
                fields = line.split("#")   # the only allowed comment marker
                if len(fields) > 1:
                    cleaned_lines.append(fields[0])
                else:
                    cleaned_lines.append(line)

            Q = json.loads("\n".join(cleaned_lines))

        except:
            print("Error: Could not read query file " + JSONfilename)
            Q = {}
            Q["errorMessage"] = []
            Q["errorMessage"].append("Error: Could not read query file " + JSONfilename)
            Q["errorStatus"] = "write and exit"
            Q["numpositions"] = 0
            Q["type"] = "symbolic"
            Q["searchFiles"] = []
            Q["numFilesSearched"] = 0
            Q["elapsedClockTime"] = 0
            Q["userMessage"] = ["Error: Could not read query file " + JSONfilename]

    else:
        print("Error: Could not find query file " + JSONfilename)
        Q = {}
        Q["errorMessage"] = []
        Q["errorMessage"].append("Error: Could not find query file " + JSONfilename)
        Q["errorStatus"] = "write and exit"
        Q["numpositions"] = 0
        Q["type"] = "symbolic"
        Q["searchFiles"] = []
        Q["numFilesSearched"] = 0
        Q["elapsedClockTime"] = 0
        Q["userMessage"] = ["Error: Could not find query file " + JSONfilename]

    # use the filename if no name is specified
    if not "name" in Q:
        Q["name"] = os.path.basename(JSONfilename).replace(".json","")

    if not Q["name"]:
        Q["name"] = os.path.basename(JSONfilename).replace(".json","")

    if not filename:
        return Q

    if "queryMoleculeType" in Q:
        # make sure all keys are numbers and not strings
        if type(Q["queryMoleculeType"]) is dict:
            for key in list(Q["queryMoleculeType"].keys()):
                if type(key) is str:
                    Q["queryMoleculeType"][int(key)] = Q["queryMoleculeType"][key]
                    del Q["queryMoleculeType"][key]

    if "requiredMoleculeType" in Q:
        # make sure all keys are numbers and not strings
        # if Q["requiredMoleculeType"] is a dictionary, do this:
        if type(Q["requiredMoleculeType"]) is dict:
            for key in list(Q["requiredMoleculeType"].keys()):
                if type(key) is str:
                    Q["requiredMoleculeType"][int(key)] = Q["requiredMoleculeType"][key]
                    del Q["requiredMoleculeType"][key]

    if not "numpositions" in Q:
        if "unitID" in Q:
            Q["numpositions"] = len(Q["unitID"])
        elif "unittype" in Q:
            Q["numpositions"] = len(Q["unittype"])
        elif "interactionMatrix" in Q:
            m = 0
            for key1 in Q["interactionMatrix"].keys():
                m = max(m,int(key1))
                for key2 in Q["interactionMatrix"][key1].keys():
                    m = max(m,int(key1))
            Q["numpositions"] = m + 1
        elif "queryMoleculeType" in Q:
            Q["numpositions"] = max(Q["queryMoleculeType"].keys()) + 1
        elif "requiredMoleculeType" in Q:
            Q["numpositions"] = max(Q["requiredMoleculeType"].keys()) + 1
        else:
            np = 0
            for key in Q.keys():
                fields = key.split(",")
                if len(fields) == 3 and fields[0] == "constraintMatrix":
                    np = max(np,int(fields[1]),int(fields[2]))
            if np > 0:
                Q["numpositions"] = np

    if not "numpositions" in Q:
        Q["errorMessage"].append("Error: Could not determine number of positions in query")
        return Q

    if "interactionMatrix" in Q:
        # make sure all keys are present in interactionMatrix
        # make sure all keys are numbers and not strings
        iM = emptyInteractionMatrix(Q["numpositions"])

        for i in range(0,Q["numpositions"]):
            for j in range(0,Q["numpositions"]):
                if str(i) in Q["interactionMatrix"]:
                    if str(j) in Q["interactionMatrix"][str(i)]:
                        iM[i][j] = Q["interactionMatrix"][str(i)][str(j)]

        Q["interactionMatrix"] = iM

    return Q


def retrieveQueryInformation(Q):
    """
    retrieveQueryInformation reads data files to get necessary information for a query
    """

    # check for constraint matrix, move entries found there to interactionMatrix
    # constraintMatrix counts positions from 1
    # example:  Q["constraintMatrix,1,2"] = "cWW"
    for key in list(Q.keys()):
        fields = key.split(",")
        if len(fields) == 3 and fields[0] == "constraintMatrix":
            # convert position numbers to python indices
            i = int(fields[1])-1
            j = int(fields[2])-1
            if not "interactionMatrix" in Q:
                Q["interactionMatrix"] = emptyInteractionMatrix(Q["numpositions"])
            Q["interactionMatrix"][i][j] += " " + Q[key]

    # infer query type from other fields being present
    if not "type" in Q:
        if "unitID" in Q:
            if "interactionMatrix" in Q:
                Q["type"] = "mixed"
            else:
                Q["type"] = "geometric"
        else:
            Q["type"] = "symbolic"

    Q['centers'] = []
    Q['rotations'] = []

    moleculeTypes = ["RNA","DNA","protein"]

    # look up centers and rotation matrices for each unit id in the query
    if Q["type"] == "geometric" or Q["type"] == "mixed":
        chainData = {}

        if "unitID" in Q:
            if len(Q["unitID"]) > 0:
                Q["numpositions"] = len(Q["unitID"])
            else:
                print("Error:  Need to specify unit IDs for geometric or mixed query")
                Q["errorMessage"].append("Need to specify unit IDs for geometric or mixed query")
        else:
            print("Error:  Need to specify unit IDs for geometric or mixed query")
            Q["errorMessage"].append("Need to specify unit IDs for geometric or mixed query")

        for i in range(0,len(Q["unitID"])):
            Q["unitID"][i] = Q["unitID"][i].replace(" ","")  # strip spaces
            unitID = Q["unitID"][i]
            originalUnitID = Q["unitID"][i]

            fields = unitID.split('|')

            file_id = fields[0]
            if len(file_id) == 4:
                # read information about PDB files, if not already done
                if not "PDB_data_file" in Q:
                    Q["PDB_data_file"] = readPDBDatafile()  # available PDB structures, resolutions, chains

                if file_id.upper() in Q["PDB_data_file"]:
                    # use uppercase for PDB identifiers
                    file_id = file_id.upper()
                    fields[0] = fields[0].upper()

            chainString = fields[0] + '|' + fields[1] + '|' + fields[2]

            fields[3] = fields[3].upper()                # force unittype to be uppercase, A, PHE, DT, etc.
            moleculeType = getMoleculeType(fields[3])

            unitID = "|".join(fields)

            # if unitType is not specified, try finding it as NA, then as protein
            foundID = False
            for mt in moleculeTypes:
                if moleculeType == '' or moleculeType == mt:

                    print("  query_processing: Retrieving data about %s in chain %s" % (unitID,chainString))

                    if not chainString in chainData:         # only load the chain once
                        if mt == "RNA" or mt == "DNA":
                            Q, centers, rotations, ids, id_to_index, index_to_id, chainIndices = readNAPositionsFile(Q,chainString,0)
                        elif mt == "protein":
                            Q, centers, ids, id_to_index, index_to_id, chainIndices = readProteinPositionsFile(Q,fields[0],0)

                        chainData[chainString] = [centers,rotations,id_to_index,index_to_id]

                    # this code trusts that centers all have 3 components; might not always be true
                    centers = chainData[chainString][0]
                    rotations = chainData[chainString][1]
                    id_to_index = chainData[chainString][2]
                    index_to_id = chainData[chainString][3]

                    if not unitID in id_to_index:            # covers the case that unitType (sequence) is missing
                        for uid in list(id_to_index.keys()):
                            f = uid.split("|")
                            f[3] = ""                        # blank out the unitType
                            newID = "|".join(f)              # collapse list into string
                            id_to_index[newID] = id_to_index[uid] # add new unit ids without unit type to dictionary

                        if unitID in id_to_index:
                            newID = index_to_id[id_to_index[unitID]]

                    if unitID in id_to_index:
                        foundID = True
                        moleculeType = getMoleculeType(index_to_id[id_to_index[unitID]])
                        if originalUnitID != unitID:
                            print("Given " + originalUnitID + ", converted to " + unitID + ", using " + index_to_id[id_to_index[unitID]])

                        Q['centers'].append(centers[id_to_index[unitID]])
                        if moleculeType == "RNA":
                            Q['rotations'].append(rotations[id_to_index[unitID]])
                        elif moleculeType == "protein":
                            Q['rotations'].append(np.empty( shape=(0, 0) ))
                        else:
                            print("Unknown molecule type " + moleculeType + "**********")

            if not foundID:
                print("Not able to find coordinates for unit ID " + unitID + " ************** ")
                Q["errorMessage"].append("Not able to find coordinates for unit ID " + unitID)
                Q["errorStatus"] = "write and exit"
    return Q


def vectorDistance(x,y):
    return math.sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)


def calculateQueryConstraints(Q):
    """
    Parse the constraints in interactionMatrix and convert to the form
    that the search code needs.
    For geometric and mixed searches, find the minimum and maximum pairwise distances
    to maintain the specified geometric discrepancy cutoff.
    """

    # If interactionMatrix is part of Q, then parse it and overwrite various other fields of Q
    if "interactionMatrix" in Q and not Q["interactionMatrix"]:
        del Q["interactionMatrix"]

    if "interactionMatrix" in Q:

        # make sure the matrix has all of the entries expected
        while len(Q["interactionMatrix"]) < Q["numpositions"]:
            Q["interactionMatrix"].append([])

        for i in range(0,Q["numpositions"]):
            while len(Q["interactionMatrix"][i]) < Q["numpositions"]:
                Q["interactionMatrix"][i].append('')

        # parse continuity constraints
        Q["continuityConstraint"] = defaultdict(dict)
        foundContinuityConstraint = False
        for i in range(0,Q["numpositions"]):
            for j in range(i):                      # look at entries below the diagonal

                iM = Q["interactionMatrix"][i][j]
                iM = iM.replace(";"," ")            # allow entries to be separated by semicolon
                iM = iM.replace("  "," ")           # replace multiple spaces
                iM = iM.replace("  "," ")           # replace multiple spaces
                iM = iM.replace("  "," ")           # replace multiple spaces
                iM = iM.replace(",<"," <")          # be tolerant of some commas
                iM = iM.replace(",>"," >")          # be tolerant of some commas
                iM = iM.replace(",="," =")          # be tolerant of some commas

                continuityConstraint = None         # default should be that there is no constraint

                if iM != None and len(iM) > 0:
                    continuityConstraint = ["", -1000000000, 1000000000]  # default, very open
                    iM = iM.replace("next","=1 >").replace("Next","=1 >")
                    iM = iM.replace("previous","=1 <").replace("Previous","=1 <")
                    continuityInfo = iM.split(" ")

                    for k in range(0,len(continuityInfo)):
                        continuityInfo[k] = continuityInfo[k].replace("before","<")
                        continuityInfo[k] = continuityInfo[k].replace("after",">")
                        continuityInfo[k] = continuityInfo[k].replace("=<","<=")
                        continuityInfo[k] = continuityInfo[k].replace("=>",">=")

                    for k in range(0,len(continuityInfo)):
                        constraint = continuityInfo[k]

                        if constraint == "=1":
                            continuityInfo[k] = "<=1"      # code as a between constraint; faster
                        if constraint == "=1,2":
                            continuityInfo[k] = "<=2"      # code as a between constraint; faster
                        if constraint == "=1,2,3":
                            continuityInfo[k] = "<=3"      # code as a between constraint; faster

                        if constraint == ">":
                            foundContinuityConstraint = True
                            continuityConstraint[0] = "between" #NTi > NTj
                            continuityConstraint[1] = 0
                            continuityInfo[k] = " "    # remove this constraint, already done
                        elif constraint == "<":
                            foundContinuityConstraint = True
                            continuityConstraint[0] = "between" #NTi < NTj
                            continuityConstraint[2] = 0
                            continuityInfo[k] = " "    # remove this constraint, already done

                    for k in range(0,len(continuityInfo)):
                        constraint = continuityInfo[k]

                        if constraint.startswith("<"):
                            val = constraint.replace("<","").replace("=","")
                            if not val.isdigit():
                                Q["errorMessage"].append("Unrecognized digit, skipping continuity constraint %s" % constraint)
                                continue
                            else:
                                val = int(val)

                            foundContinuityConstraint = True

                            if constraint.startswith("<="):
                                dist = val+1         # maximum sequential distance
                            else:
                                dist = val
                            continuityConstraint[0] = "between"
                            continuityConstraint[1] = max(continuityConstraint[1],-dist)
                            continuityConstraint[2] = min(continuityConstraint[2], dist)
                            continuityInfo[k] = " "    # remove this constraint, already done

                    for k in range(0,len(continuityInfo)):
                        constraint = continuityInfo[k]
                        if constraint.startswith("="):
                            distanceText = constraint.replace("=","").split(",")

                            try:
                                distances = [int(d) for d in distanceText if d.isdigit() or d.replace("-","").isdigit]
                            except:
                                Q["errorMessage"].append("Could not interpret %d-%d continuity constraint %s, skipping it." % (i+1,j+1,constraint))
                                continue

                            if len(distances) == 0:
                                Q["errorMessage"].append("No digit found, skipping %d-%d continuity constraint %s" % (i+1,j+1,constraint))
                                continue

                            continuityConstraint[0] = "equal"
                            foundContinuityConstraint = True

                            newList = []
                            for d in distances:
                                if d > continuityConstraint[1] and d < continuityConstraint[2]:
                                    newList.append(d)
                                if -d > continuityConstraint[1] and -d < continuityConstraint[2]:
                                    newList.append(-d)
                            continuityConstraint[1] = min(newList)-1   # minimum, like a between search
                            continuityConstraint[2] = max(newList)+1   # maximum, like a between search
                            continuityConstraint.append(set(newList))
                            continuityInfo[k] = " "    # remove this constraint, already done

                    for k in range(0,len(continuityInfo)):
                        constraint = continuityInfo[k]

                        if constraint.startswith(">"):
                            val = constraint.replace(">","").replace("=","")
                            if not val.isdigit():
                                Q["errorMessage"].append("Unrecognized digit, skipping continuity constraint %s" % constraint)
                                continue
                            else:
                                val = int(val)

                            foundContinuityConstraint = True

                            if constraint.startswith(">="):
                                dist = val-1         # maximum sequential distance
                            else:
                                dist = val

                            if continuityConstraint[0] == "":                # no previous constraint
                                continuityConstraint[0] = "outside"
                                continuityConstraint[1] = -dist
                                continuityConstraint[2] = dist
                            elif continuityConstraint[0] == "between":
                                if continuityConstraint[1] >= dist:    # outside constraint already met
                                    continue
                                elif continuityConstraint[2] <= -dist: # outside constraint already met
                                    continue
                                elif continuityConstraint[1] > -dist: # no < possibilities remain
                                    continuityConstraint[1] = max(continuityConstraint[1],dist)
                                    continuityConstraint[2] = max(continuityConstraint[2],dist)
                                elif continuityConstraint[2] < dist:  # no > possibilities remain
                                    continuityConstraint[1] = min(continuityConstraint[1],-dist)
                                    continuityConstraint[2] = min(continuityConstraint[2],-dist)
                                else:
                                    leftList = range(continuityConstraint[1]+1,-dist)
                                    rightList = range(dist+1,continuityConstraint[2])
                                    newList = leftList + rightList
                                    continuityConstraint[0] = "equal"
                                    continuityConstraint[1] = min(newList)-1   # minimum, like a between search
                                    continuityConstraint[2] = max(newList)+1   # maximum, like a between search
                                    continuityConstraint.append(newList)
                            elif continuityConstraint[0] == "equal":
                                newList = []
                                for d in continuityConstraint[3]:
                                    if d < -dist or d > dist:
                                        newList.append(d)
                                continuityConstraint[1] = min(newList)-1   # minimum, like a between search
                                continuityConstraint[2] = max(newList)+1   # maximum, like a between search
                                continuityConstraint[3] = set(newList)

                    Q["continuityConstraint"][j][i] = continuityConstraint

        if not foundContinuityConstraint:
            del Q["continuityConstraint"]          # remove this field so we don't take time to check

        # Parse interaction matrix to look for interaction constraints

        RNACombinationConstraints = ['AA','AC','AG','AU','CA','CC','CG','CU','GA','GC','GG','GU','UA','UC','UG','UU']

        Q["requiredInteractions"] = emptyInteractionList(Q["numpositions"])
        Q["prohibitedInteractions"] = emptyInteractionList(Q["numpositions"])
        Q["crossingNumber"] = emptyInteractionList(Q["numpositions"])
        Q["combinationConstraint"] = emptyInteractionList(Q["numpositions"])

        foundRequiredInteraction = False
        foundProhibitedInteraction = False
        foundCrossingNumber = False
        foundCombinationConstraint = False
        Q["alternateInteractions"] = set([])       # for _exp and possibly others
        foundAlternateInteractions = False

        for i in range(Q["numpositions"]):
            for j in range(Q["numpositions"]):

                #print("interactions",i,j)

                requiredInteractions = []
                prohibitedInteractions = []
                combinationConstraints = []

                iM = Q["interactionMatrix"][i][j]
                if iM != None and len(iM) > 0:

                    rawInteractionConstraints = iM.split(" ")

                    # allow the abbreviation "n+interaction" for "interaction ninteraction"
                    interactionConstraints = []
                    for constraint in rawInteractionConstraints:
                        if "n+" in constraint:
                            interactionConstraints.append(constraint.replace("n+",""))  # true interaction
                            interactionConstraints.append(constraint.replace("n+","n")) # near interaction
                        else:
                            interactionConstraints.append(constraint)

                    for constraint in interactionConstraints:
                        #print("Constraint: " + constraint)

                        # enable reading alternate annotations from different annotation files
                        # these are distinguished by "_" in the annotation and search terms
                        constraintType = ""
                        alternate = ""
                        if "_" in constraint and not "chi" in constraint.lower() and not "crossing" in constraint.lower() and not "chainlength" in constraint.lower():
                            constraintType = constraint.split("_")[0].replace("~","")
                            alternate = "_" + constraint.split("_")[1]
                            Q["alternateInteractions"].add(alternate)
                            foundAlternateInteractions = True

                        # special treatment for certain pairs typed into blue boxes below diagonal
                        asymmetricPair = False
                        if "BPh" in constraint or "BR" in constraint:
                            asymmetricPair = True
                        elif "sO" in constraint or "s3O" in constraint or "s5O" in constraint or "sO3" in constraint or "sO5" in constraint:
                            asymmetricPair = True

                        if j >= i or asymmetricPair:  # yellow box or special interaction
                            if len(constraint) > 0 and constraint[0] == "~":
                                constraint = constraint[1:]
                                if constraint in synonym:
                                    constraints = synonym[constraint]
                                    prohibitedInteractions.extend(constraints)
                                    foundProhibitedInteraction = True
                                elif constraintType in synonym:
                                    constraints = [c+alternate for c in synonym[constraintType]]
                                    prohibitedInteractions.extend(constraints)
                                    foundProhibitedInteraction = True
                                elif constraint in allInteractionConstraints:
                                    prohibitedInteractions.append(constraint)
                                    foundProhibitedInteraction = True
                                elif constraintType in allInteractionConstraints:
                                    prohibitedInteractions.append(constraint)
                                    foundProhibitedInteraction = True
                            else:
                                if constraint in synonym:
                                    constraints = synonym[constraint]
                                    requiredInteractions.extend(constraints)
                                    foundRequiredInteraction = True
                                elif constraintType in synonym:
                                    constraints = [c+alternate for c in synonym[constraintType]]
                                    requiredInteractions.extend(constraints)
                                    foundRequiredInteraction = True
                                elif constraint in crossing:
                                    Q["crossingNumber"][i][j] = crossing[constraint]
                                    foundCrossingNumber = True
                                elif "crossing" in constraint:
                                    limits = constraint.replace("crossing","").split("_")
                                    print("Crossing number limits " + str(limits))
                                    Q["crossingNumber"][i][j] = [int(limits[1]),int(limits[2])]
                                    print("Crossing number limits " + str(Q["crossingNumber"][i][j]))
                                    foundCrossingNumber = True
                                elif j>i and "," in constraint:
                                    units = constraint.split(",")
                                    combinationConstraints.append((units[0],units[1]))
                                    foundCombinationConstraint = True
                                elif j>i and constraint in RNACombinationConstraints:  # for backward compatibility
                                    combinationConstraints.append((constraint[0],constraint[1]))
                                    foundCombinationConstraint = True
                                elif constraint in allInteractionConstraints:
                                    requiredInteractions.append(constraint)
                                    foundRequiredInteraction = True
                                elif constraintType in allInteractionConstraints:
                                    requiredInteractions.append(constraint)
                                    foundRequiredInteraction = True

                    if len(requiredInteractions) > 0:
                        if not requiredInteractions[-1] == "and":
                            requiredInteractions.append("and")

                        while len(requiredInteractions) > 0 and requiredInteractions[0] == "and":
                            requiredInteractions = requiredInteractions[1:]


                Q["requiredInteractions"][i][j] = requiredInteractions
                Q["prohibitedInteractions"][i][j] = prohibitedInteractions
                Q["combinationConstraint"][i][j] = combinationConstraints

        if not foundRequiredInteraction:
            del Q["requiredInteractions"]          # remove this field so we don't take time to check

        if not foundProhibitedInteraction:
            del Q["prohibitedInteractions"]          # remove this field so we don't take time to check

        if not foundAlternateInteractions:
            del Q["alternateInteractions"]          # remove this field so we don't take time to check

        if not foundCrossingNumber:
            del Q["crossingNumber"]

        if not foundCombinationConstraint:
            del Q["combinationConstraint"]

        # Parse Unary Constraints

        letterToConstraint ={ 'A': ['A'],
        'C': ['C'],
        'G': ['G'],
        'U': ['U'],
        'M' :[ 'A', 'C' ],
        'R' : ['A', 'G'],
        'W' : ['A' , 'U'],
        'S' : ['C' , 'G'],
        'Y' : ['C' , 'U'],
        'K' : ['G' , 'U'],
        'V' : ['A', 'C', 'G'],
        'H' : ['A', 'C', 'U'],
        'D' : ['A', 'G', 'U'],
        'B' : ['C', 'G', 'U'],
        'N' : ['A', 'C', 'G', 'U']
        }

        Q["requiredUnitType"] = [None] * Q["numpositions"]
        Q["requiredMoleculeType"] = [None] * Q["numpositions"]

        Q["glycosidicBondOrientation"] = [None] * Q["numpositions"]
        foundGlycosidicBondOrientation = False
        Q["chiAngle"] = [None] * Q["numpositions"]
        Q["chainLength"] = [None] * Q["numpositions"]
        foundChiAngle = False
        foundChainLength = False

        # process unary constraints
        for i in range(Q["numpositions"]):
            Q["requiredUnitType"][i] = []
            Q["requiredMoleculeType"][i] = []
            Q["glycosidicBondOrientation"][i] = []
            Q["chiAngle"][i] = []
            Q["chainLength"][i] = []

            if Q["interactionMatrix"][i][i] == None or len(Q["interactionMatrix"][i][i]) == 0:
                if 'RNA' not in Q["requiredMoleculeType"][i]:
                    Q["requiredMoleculeType"][i].append('RNA')

            else:
                iMtext = Q["interactionMatrix"][i][i].replace(","," ")

                iMarray = iMtext.split(" ")
                for iM in iMarray:
                    if len(iM) == 0:
                        continue
                    if iM == "x" or iM == "X":
                        iM = "protein"

                    if iM.lower() == 'modified':
                        if 'rna' in iMtext.lower():
                            Q["requiredUnitType"][i] += RNA_modified_list
                            if not 'RNA' in Q["requiredMoleculeType"][i]:
                                Q["requiredMoleculeType"][i].append('RNA')
                        elif 'dna' in iMtext.lower():
                            Q["requiredUnitType"][i] += DNA_modified_list
                            if not 'DNA' in Q["requiredMoleculeType"][i]:
                                Q["requiredMoleculeType"][i].append('DNA')
                        else:
                            Q["requiredUnitType"][i] += RNA_modified_list
                            Q["requiredUnitType"][i] += DNA_modified_list
                            if not 'RNA' in Q["requiredMoleculeType"][i]:
                                Q["requiredMoleculeType"][i].append('RNA')
                            if not 'DNA' in Q["requiredMoleculeType"][i]:
                                Q["requiredMoleculeType"][i].append('DNA')


                    elif iM in RNA_unit_types:
                        Q["requiredUnitType"][i].append(iM)
                        if not 'RNA' in Q["requiredMoleculeType"][i]:
                            Q["requiredMoleculeType"][i].append('RNA')

                    elif iM in letterToConstraint:
                        Q["requiredUnitType"][i] += letterToConstraint[iM]
                        if not 'RNA' in Q["requiredMoleculeType"][i]:
                            Q["requiredMoleculeType"][i].append('RNA')

                    elif iM in DNA_unit_types:
                        Q["requiredUnitType"][i].append(iM)
                        if not 'DNA' in Q["requiredMoleculeType"][i]:
                            Q["requiredMoleculeType"][i].append('DNA')

                    elif iM in protein_unit_types:
                        Q["requiredUnitType"][i].append(iM)
                        if 'protein' not in Q["requiredMoleculeType"][i]:
                            Q["requiredMoleculeType"][i].append('protein')

                    elif iM == "RNA":
                        if not 'RNA' in Q["requiredMoleculeType"][i]:
                            Q["requiredMoleculeType"][i].append('RNA')

                    elif iM == "DNA":
                        if not 'DNA' in Q["requiredMoleculeType"][i]:
                            Q["requiredMoleculeType"][i].append('DNA')

                    elif iM == "protein":
                        if not 'protein' in Q["requiredMoleculeType"][i]:
                            Q["requiredMoleculeType"][i].append('protein')

                    elif iM.lower() == "syn":
                        Q["glycosidicBondOrientation"][i].append("syn")
                        foundGlycosidicBondOrientation = True

                    elif iM.lower() == "anti":
                        Q["glycosidicBondOrientation"][i].append("anti")
                        foundGlycosidicBondOrientation = True

                    elif iM.lower() == "is" or iM.lower() == "int syn" or iM.lower() == "int_syn" or iM.lower() == "int-syn" or iM.lower() == "intermediate syn":
                        Q["glycosidicBondOrientation"][i].append("int_syn")
                        foundGlycosidicBondOrientation = True

                    elif iM.lower() == "~syn":
                        Q["glycosidicBondOrientation"][i].append("anti")
                        Q["glycosidicBondOrientation"][i].append("int_syn")
                        foundGlycosidicBondOrientation = True

                    elif iM.lower() == "~anti":
                        Q["glycosidicBondOrientation"][i].append("syn")
                        Q["glycosidicBondOrientation"][i].append("int_syn")
                        foundGlycosidicBondOrientation = True

                    elif iM.lower() == "~is" or iM.lower() == "~int syn" or iM.lower() == "~int_syn" or iM.lower() == "~int-syn" or iM.lower() == "~intermediate syn":
                        Q["glycosidicBondOrientation"][i].append("anti")
                        Q["glycosidicBondOrientation"][i].append("syn")
                        foundGlycosidicBondOrientation = True

                    elif iM.lower() == "glyco" or iM.lower() == "glycosidic":
                        Q["showGlycosidicBondOrientation"] = True
                        foundGlycosidicBondOrientation = True

                    elif "chi" in iM.lower():
                        fields = iM.lower().replace("(","_").replace(")","")
                        fields = fields.replace("[","_").replace("]","")
                        fields = fields.replace(":","_")
                        fields = fields.split('_')
                        if len(fields) == 3:
                            try:
                                a = float(fields[1])
                                b = float(fields[2])
                                foundChiAngle = True
                                if a < b:
                                    Q["chiAngle"][i] = ["between",a,b]
                                else:
                                    Q["chiAngle"][i] = ["outside",a,b]
                            except:
                                Q["errorMessage"].append('Could not parse chi constraint %d' % i)
                        else:
                            Q["errorMessage"].append('Could not parse chi constraint %d' % i)

                    elif "chainlength" in iM.lower():
                        fields = iM.lower().replace("(","_").replace(")","")
                        fields = fields.replace("[","_").replace("]","")
                        fields = fields.replace(":","_")
                        fields = fields.split('_')
                        if len(fields) == 3:
                            try:
                                if fields[1].lower() == 'inf':
                                    a = float('inf')
                                else:
                                    a = int(fields[1])
                                if fields[2].lower() == 'inf':
                                    b = float('inf')
                                else:
                                    b = int(fields[2])
                                foundChainLength = True
                                if a < b:
                                    Q["chainLength"][i] = ["between",a,b]
                                else:
                                    Q["chainLength"][i] = ["outside",a,b]
                            except:
                                Q["errorMessage"].append('Could not parse chainlength constraint %d' % i)
                        else:
                            Q["errorMessage"].append('Could not parse chainlength constraint %d' % i)

                    else:
                        # for example an unrecognized modified nucleotide
                        Q["requiredUnitType"][i].append(iM)
                        Q["errorMessage"].append('Did not recognize %s' % iM)


            if len(Q["requiredMoleculeType"][i]) == 0:
                Q["requiredMoleculeType"][i] = ['RNA']               # default constraint

        if not foundGlycosidicBondOrientation:
            del Q["glycosidicBondOrientation"]
        if not foundChiAngle:
            del Q["chiAngle"]
        if not foundChainLength:
            del Q["chainLength"]

        # for unary constraints, we need to consider the possibility of things separated by spaces or commas
        # could be "AG" for RNA A or G
        # could be "SER VAL" for protein, and wanting just these two amino acids
        # could be "DG" which is "DNA, G"
        # could be "DG DA" which is DNA, G or A
        # could be "RNA DNA" which allows either
        # could be "RNA protein"
        # could be "protein"

        # WORRY ABOUT EXCLUSION LATER


    Q["activeInteractions"] = []

    for i in range(0,Q["numpositions"]):
        for j in range(0,Q["numpositions"]):
            if 'requiredInteractions' in Q:
                Q["activeInteractions"].extend(Q["requiredInteractions"][i][j])
            if 'prohibitedInteractions' in Q:
                Q["activeInteractions"].extend(Q["prohibitedInteractions"][i][j])

    Q["activeInteractions"] = list(set(Q["activeInteractions"]))

    if "and" in Q["activeInteractions"]:
        Q["activeInteractions"].remove("and")

    if not "requiredMoleculeType" in Q or not Q["requiredMoleculeType"]:
        Q["requiredMoleculeType"] = [["RNA"]] * Q["numpositions"]

    if Q["type"] == "mixed" or Q["type"] == "geometric":
        # determine the type of each unit to facilitate retrieving its information
        # this must be fragile with modified nucleotides and amino acids
        if not "queryMoleculeType" in Q:
            Q["queryMoleculeType"] = [None] * Q["numpositions"]
            for i in range(Q["numpositions"]):
                Q["queryMoleculeType"][i] = getMoleculeType(Q["unitID"][i])

    if not "locationWeight" in Q:
        Q["locationWeight"] = [1] * Q["numpositions"]

    for i in range(0,Q["numpositions"]):
        s = sum(Q["locationWeight"])
        Q["locationWeight"][i] = Q["numpositions"] * Q["locationWeight"][i]/s

    # incorporate query units into the search constraints
    if Q["type"] == "geometric" or Q["type"] == "mixed":
        if not "discrepancy" in Q:
            Q["discrepancy"] = 0.3
            print("No discrepancy specified, using 0.3")
            Q["userMessage"].append("No discrepancy specified, using 0.3")
        Q["requireddistanceminimum"] = defaultdict(dict)
        Q["requireddistancemaximum"] = defaultdict(dict)
        Q["SSCutoff"] = [None] * Q["numpositions"] # SSCutoff = sum of squares cutoff
        Q["largestMaxRange"] = 0
        Q["distance"] = np.zeros((Q["numpositions"],Q["numpositions"]))
        s = 0

        for i in range(0,Q["numpositions"]):
            s = s + Q["locationWeight"][i]

            Q["SSCutoff"][i] = (Q["numpositions"]**2)*(Q["discrepancy"]**2)*s;
            for j in range(i+1,Q["numpositions"]):

                Q["distance"][i][j] = np.linalg.norm(Q["centers"][i] -Q["centers"][j])
                Q["distance"][j][i] = Q["distance"][i][j]

                if Q["numpositions"] > 2:
                    wi = Q["locationWeight"][i]
                    wj = Q["locationWeight"][j]
                    delta = math.sqrt((wi + wj) / (wi * wj)) * Q["numpositions"] * Q["discrepancy"]
                else:
                    delta = Q["numpositions"] * Q["discrepancy"]

                Q["requireddistanceminimum"][i][j] = Q["distance"][i][j] - delta
                Q["requireddistancemaximum"][i][j] = Q["distance"][i][j] + delta
                Q["requireddistanceminimum"][j][i] = Q["distance"][i][j] - delta
                Q["requireddistancemaximum"][j][i] = Q["distance"][i][j] + delta

                Q["largestMaxRange"] = max(Q["largestMaxRange"],Q["distance"][i][j] + delta)


    # note what molecule type we are searching for, to load the right chains
    searchingRNA = False
    for i in range(0,Q["numpositions"]):
        if 'RNA' in Q["requiredMoleculeType"][i]:
            searchingRNA = True
    searchingDNA = False
    for i in range(0,Q["numpositions"]):
        if 'DNA' in Q["requiredMoleculeType"][i]:
            searchingDNA = True
    searchingProtein = False
    for i in range(0,Q["numpositions"]):
        if 'protein' in Q["requiredMoleculeType"][i]:
            searchingProtein = True

    # process the given list of files to search and resolve into chains and IFEs
    IFEList = []            # accumulate all chain and PDB names
    unit_data_file = []     # list of files stored in DATAPATHUNITS
    # unit_file_id_to_chains = defaultdict(list)   # map from file_id of unit_data_file to list of chains

    search_file_list = []

    # process any wildcards
    for search_file in Q["searchFiles"]:
        if "*" in search_file:
            print('  query_processing: Found wildcard in %s' % search_file)

            search_file = search_file.strip()

            (path_to_file,search_filename) = os.path.split(search_file)

            if os.path.exists(path_to_file):
                print('  query_processing: Found path %s' % path_to_file)
                file_list = os.listdir(path_to_file)

                import fnmatch

                print('  query_processing: Found %d files to search' % len(file_list))

                # Loop over the filenames
                for filename in file_list:
                    if fnmatch.fnmatch(filename, "*.pdb"):
                        search_file_list.append(os.path.join(path_to_file,filename))
        else:
            search_file_list.append(search_file)

    # loop over the specified search files
    for search_file in search_file_list:

        search_file = search_file.strip()

        (path_to_file,search_filename) = os.path.split(search_file)

        # not sure exactly how this works when search_file already has a path; probably OK
        search_file_id = search_filename.replace('.gz','').replace('.cif','').replace('.pdb','')

        if len(search_file) == 0:
            continue

        # look up representative sets, if requested, and replace with IFE names
        if "nrlist" in search_file:           # referring to lists that are posted online
            listLoaded = False
            Q, DATAPATHUNITS = get_DATAPATHUNITS(Q)
            pathAndFileName = os.path.join(DATAPATHUNITS,'Representative_sets.pickle')

            # check in local file of representative sets first, in case already loaded
            if os.path.exists(pathAndFileName):
                with open(pathAndFileName, 'rb') as fh:
                    representativeSets = pickle.load(fh)
                if search_file in representativeSets and len(representativeSets[search_file]) > 0:
                    IFEList.extend(representativeSets[search_file])
                    listLoaded = True

            if not listLoaded:
                newList = []

                # requesting NMR structures only
                # this is a bit of a hack
                if 'NMR/csv' in search_file:
                    url = search_file.replace('NMR/csv','all/csv')
                    NMRonly = True
                else:
                    url = search_file
                    NMRonly = False

                try:
                    # download representative set list from URL given in search_file variable
                    if sys.version_info[0] < 3:
                        f = urllib.urlopen(url)          # python 2
                    else:
                        f = urllib.request.urlopen(url)  # python 3
                    myfile = f.read()
                except Exception as e:
                    Q["errorMessage"].append("Not able to download representative set %s" % search_file)
                    print("Not able to download representative set; check the internet connection.")
                    print(e)
                    myfile = []

                if "Arial" in str(myfile):
                    print("Problem with downloading representative set, got:")
                    print(myfile)
                elif len(myfile) > 0:
                    allLines = myfile.split(b'\n')
                    for line in allLines:
                        # CSV file, split by comma
                        fields = line.split(b",")
                        if len(fields) > 1:
                            # second column has representative IFE; remove ""
                            IFE = fields[1].replace(b'"',b'')
                            if NMRonly:
                                file_id = IFE.split("|")[0]
                                if len(file_id) == 4:
                                    if not "PDB_data_file" in Q:
                                        Q["PDB_data_file"] = readPDBDatafile()  # available PDB structures, resolutions, chains

                                    file_id = IFE[0:4]
                                    if file_id in list(Q["PDB_data_file"]):
                                        if 'method' in list(Q["PDB_data_file"][file_id]):
                                            if 'NMR' in Q["PDB_data_file"][file_id]['method']:
                                                newList.append(IFE.decode("ascii"))
                                                # newList.append(IFE)
                            elif len(IFE) > 1:
                                newList.append(IFE.decode("ascii"))
                                # newList.append(IFE)

                    IFEList.extend(newList)

                    if os.path.exists(pathAndFileName):
                        with open(pathAndFileName, 'rb') as fh:
                            representativeSets = pickle.load(fh)
                    else:
                        representativeSets = {}

                    representativeSets[search_file] = newList
                    pickle.dump(representativeSets, open(pathAndFileName, "wb" ), 2)

            IFEList = sorted(IFEList)

            continue

        # print('  query_processing: search_file_id is %s' % search_file_id)

        if "|" in search_file:
            # a chain or IFE is specified
            chains = search_file.split("+")
            for k in range(0,len(chains)):
                fields = chains[k].split("|")
                # convert PDB ids to uppercase to avoid case problems
                fields[0] = fields[0].upper()
                chains[k] = "|".join(fields)
            IFEList.append("+".join(chains))

            continue

        if len(search_file_id) == 4:
            # process as a PDB id
            # if RNA or DNA is being searched for, convert any 4-letter PDB IDs to strings with chain strings separated by + signs
            # which is how RNA chains are stored in .pickle files for speed
            # not needed for protein because they are not stored by chain

            # load mapping from PDB id to chain and other information, if not already done
            if not "PDB_data_file" in Q:
                Q, DATAPATHUNITS = get_DATAPATHUNITS(Q)
                Q["PDB_data_file"] = readPDBDatafile(DATAPATHUNITS)  # available PDB structures, resolutions, chains

            search_file_id_upper = search_file_id.upper()

            if search_file_id_upper in Q["PDB_data_file"].keys():
                chains = []

                #print(Q["PDB_data_file"][search_file])

                if searchingRNA and 'RNA' in Q["PDB_data_file"][search_file_id_upper]['chains']:
                    chains += ["|".join([search_file_id_upper,'1',x]) for x in Q["PDB_data_file"][search_file_id_upper]['chains']['RNA']]
                if searchingDNA and 'DNA' in Q["PDB_data_file"][search_file_id_upper]['chains']:
                    chains += ["|".join([search_file_id_upper,'1',x]) for x in Q["PDB_data_file"][search_file_id_upper]['chains']['DNA']]
                IFEList.append("+".join(chains))

                print("  query_processing: Found chains %s in PDB file %s" % (chains,search_file_id_upper))

                continue

        # on the server, only PDB files will be searched
        # if we got this far, no PDB file will be found
        if SERVER:
            print('  Unknown file %s' % search_file)
            Q["errorMessage"].append("Did not recognize %s as a file from the Protein Data Bank" % search_file)
            continue

        # look for .pickle files for search_file in the units directory
        Q, DATAPATHUNITS = get_DATAPATHUNITS(Q)
        units_path = os.path.join(DATAPATHUNITS)
        if os.path.exists(units_path):

            # read the units folder to see what chains are available and list those
            # read directory listing for units directory, get filenames that include search_file in the name
            if not unit_data_file:
                print('  query_processing: Running readUnitFilenames to find previously processed unit files')
                if "PDB_data_file" in Q:
                    unit_data_file = readUnitFileNames(Q,Q["PDB_data_file"])
                else:
                    unit_data_file = readUnitFileNames(Q,set())

            if search_file_id in unit_data_file:

                # print("  query_processing: found the chains %s in the units directory" % unit_data_file[search_file_id]['chains'])

                # found the list of chains in the /units directory
                chains = []
                if (searchingRNA or searchingDNA) and 'NA' in unit_data_file[search_file_id]['chains']:
                    chains += ["|".join([search_file_id,'1',x]) for x in unit_data_file[search_file_id]['chains']['NA']]
                IFEList.append("+".join(chains))

                continue

            # other code that tries to do the same thing but not as well!
            # get directory listing

            # print('  query_processing: Looking for unit file for %s in %s' % (search_file_id,units_path))

            # if not unit_file_id_to_chains:
            #     file_list = os.listdir(units_path)
            #     # Loop over the filenames
            #     for filename in file_list:
            #         fields = filename.split("_")
            #         # use all but last three fields as the file_id
            #         file_id = "_".join(fields[0:-3])
            #         # print(fields)
            #         # print(file_id)
            #         # print("_".join(fields[0:-1]))
            #         # store file_id, model, and chain
            #         unit_file_id_to_chains[file_id].append("_".join(fields[0:-1]))


            # if search_file_id in unit_file_id_to_chains:
            #     chains = unit_file_id_to_chains[search_file_id]
            #     IFEList.append("+".join(chains))
            #     print("  query_processing: Found chains %s in units directory for %s" % (chains,search_file_id_upper))
            #     continue

        # if search_file already has a path and the file exists
        if os.path.exists(search_file) or os.path.exists(search_file + ".gz"):
            from file_reading import processPDBFile

            print('  query_processing: About to run processPDBFile')
            chains, file_id, messages = processPDBFile(Q,search_file,search_file_id)

            print('  query_processing: Ran processPDBFile')
            # print(chains)
            # print(file_id)
            # print(messages)

            IFEList.append("+".join(chains))
            continue

        # look for file in the cif or pdb directory
        # if no .pickle files are available, process the .cif or .pdb file
        Q, CIFPATH = get_CIFPATH(Q)

        if CIFPATH and os.path.exists(CIFPATH):
            print('  query_processing: Looking for %s in %s' % (search_file_id,CIFPATH))

            CIFPATH_search_file = os.path.join(CIFPATH,search_file)

            if os.path.exists(CIFPATH_search_file):

                from file_reading import processPDBFile

                print('  query_processing: about to run processPDBFile')
                chains, file_id, messages = processPDBFile(Q,CIFPATH_search_file)

                print('  query_processing: Ran processPDBFile')
                print(chains)
                print(file_id)
                print(messages)

                IFEList.append("+".join(chains))
                continue

        print('Unable to find search file %s' % search_file)
        Q["errorMessage"].append("Unable to find search file %s" % search_file)

    # finalize the list of IFEs to search
    Q["searchFiles"] = [x for x in IFEList if len(x) > 0]

    #set cutoffs for geometric or mixed searches
    if Q["type"] == "geometric" or Q["type"] == "mixed":
        cutoff = Q['SSCutoff'] #rename sscutoffs
        for i in range(len(cutoff)):
            if cutoff[i] != 'Inf':
                cutoff[i] = float(cutoff[i])
            else:
                cutoff[i] = float("inf")
        Q["cutoff"] = cutoff

    return Q
