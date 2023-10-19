from file_reading import readNAPositionsFile
from file_reading import readProteinPositionsFile
from file_reading import readPDBDatafile
from collections import defaultdict
from fr3d_configuration import DATAPATH
from fr3d_configuration import JSONPATH
from fr3d_configuration import SERVER

from fr3d.modified_parent_mapping import modified_nucleotides

import numpy as np
import math
import urllib
import json
import sys
import pickle
import string
import os

# identify codes that go with each type of molecule
RNAUnitTypes = ["A","C","G","U"] + list(modified_nucleotides.keys())
DNAUnitTypes = ["DA","DC","DG","DT"]
proteinUnitTypes = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","PYL","SER","SEC","THR","TRP","TYR","VAL","ASX","GLX","XAA","XLE"]

RNA_modified_list = []
DNA_modified_list = []

for x in modified_nucleotides.keys():
    if modified_nucleotides[x]['standard'] in ['A','C','G','U']:
        RNA_modified_list.append(x)
    else:
        DNA_modified_list.append(x)
RNA_modified_list = sorted(RNA_modified_list)
DNA_modified_list = sorted(DNA_modified_list)


# synonyms and abbreviations for pairwise constraints
synonym = {}
synonym["pair"] = ['acWW','cWW','tWW','cHH','tHH','cSS','tSS','cWH','cHW','tWH','tHW','cWS','cSW','tWS','tSW','cHS','cSH','tHS','tSH']
synonym["cis"] = ['acWW','cWW','cHH','cSS','cWH','cHW','cWS','cSW','cHS','cSH']
synonym["trans"] = ['tWW','tHH','tSS','tWH','tHW','tWS','tSW','tHS','tSH']
synonym["stack"] = ['s35','s53','s33','s55']
synonym["BPh"] = ['0BPh','1BPh','2BPh','3BPh','4BPh','5BPh','6BPh','7BPh','8BPh','9BPh']
synonym["BR"] = ['0BR','1BR','2BR','3BR','4BR','5BR','6BR','7BR','8BR','9BR']
synonym["coplanar"] = ['cp']

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

# Determine molecule type from unit ID by identifying unit type specifiers
def getMoleculeType(unitType):
    if "|" in unitType:
        fields = unitType.split('|')
        unitType = fields[3]

    if sys.version_info[0] < 3:
        unitType = string.upper(unitType)    # convert to uppercase
    else:
        unitType = unitType.upper()

    if unitType in RNAUnitTypes:
        return "RNA"
    elif unitType in proteinUnitTypes:
        return "protein"
    elif unitType in DNAUnitTypes:
        return "DNA"
    elif unitType in modified_nucleotides.keys():
        return modified_nucleotides[unitType]["standard"]
    else:
        return ""

# set up an empty list of required interactions,
# then the user only needs to set the actual constraints needed
def emptyInteractionList(n):
    emptyList = defaultdict(dict)
    for i in range(0,n):
        for j in range(0,n):
            emptyList[i][j] = []
    return emptyList

# set up an empty list of required interactions,
# then the user only needs to set the actual constraints needed
def emptyInteractionMatrix(n):
    emptyList = defaultdict(dict)
    for i in range(0,n):
        for j in range(0,n):
            emptyList[i][j] = ""
    return emptyList


def readUnitFileNames(PDB_data_file):
    """
    Read all filenames in DATAPATH + '/units' to know which chains are available,
    that are not downloaded from PDB
    """

    import glob

    units_path = os.path.join(DATAPATH,'units','*')
    allfiles = glob.glob(units_path)
    print("  query_processing: There are %d files in %s" % (len(allfiles),units_path))

    # mimic the structure of PDB_data_file
    unit_data_file = {}

    for file in allfiles:
        # remove path
        chain_filename = os.path.split(file)[1]

        # split by '_' character to identify molecule type
        type_field = chain_filename.split('_')
        if len(type_field) < 3:
            continue
        else:
            # move toward all the NA files being _NA.pickle
            if type_field[-1] == 'RNA.pickle':
                molecule_type = 'NA'
            elif type_field[-1] == 'DNA.pickle':
                molecule_type = 'NA'
            elif type_field[-1] == 'NA.pickle':
                molecule_type = 'NA'
            elif type_field[-1] == 'protein.pickle':
                molecule_type = 'protein'
            else:
                continue

            # split out chain to see if it is already known
            fields = chain_filename.split('-')

            if len(fields) < 3:
                continue
            else:
                chain = fields[-1].split('_')[0]
                model = fields[-2]
                filename = '-'.join(fields[0:-2])

                # skip files that are already in PDB_data_file
                if not filename in PDB_data_file.keys():
                    if not filename in unit_data_file:
                        unit_data_file[filename] = {}
                        unit_data_file[filename]['chains'] = {}
                        unit_data_file[filename]['chains']['NA'] = []
                        unit_data_file[filename]['chains']['protein'] = []
                        unit_data_file[filename]['model'] = []

                    unit_data_file[filename]['chains'][molecule_type].append(chain)
                    unit_data_file[filename]['model'].append(model)

                    #print('  query_processing: Found units for chain %s|%s|%s' % (filename,model,chain))

    # store only unique model names
    for filename in unit_data_file.keys():
        unit_data_file[filename]['model'] = list(set(unit_data_file[filename]['model']))

    return unit_data_file

def readQueryFromJSON(JSONfilename):
    # Read query specification in JSON file and interpret as a query specification

    directory = JSONPATH
    try:
        os.stat(directory)
    except:
        os.mkdir(directory)
        print("Made " + directory + " directory")

    pathAndFileName = JSONPATH + JSONfilename

    id = JSONfilename.replace("Query_","").replace(".json","")

    if not os.path.exists(pathAndFileName):
        queryURL = "http://rna.bgsu.edu/webfr3d/Results/" + id + "/" + JSONfilename
        print("Downloading "+JSONfilename+" from "+queryURL)
        if sys.version_info[0] < 3:
            urllib.urlretrieve(queryURL, pathAndFileName)  # python 2
        else:
            urllib.request.urlretrieve(queryURL, pathAndFileName)  # python 3

    with open(JSONPATH + JSONfilename) as json_file:
        Q = json.load(json_file)

    if not "name" in Q:
        Q["name"] = id

    if not Q["name"]:
        Q["name"] = id

    return Q

def retrieveQueryInformation(Q):
    """
    retrieveQueryInformation reads data files to get necessary information for a query
    """

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

        for i in range(0,Q["numpositions"]):
            Q["unitID"][i] = Q["unitID"][i].replace(" ","")  # strip spaces
            unitID = Q["unitID"][i]
            originalUnitID = Q["unitID"][i]
            foundID = False

            fields = unitID.split('|')

            fields[0] = fields[0].upper()                # force PDB ID to be uppercase for unix
            fields[3] = fields[3].upper()                # force unittype to be uppercase, A, PHE, DT, etc.

            unitID = "|".join(fields)

            chainString = fields[0] + '-' + fields[1] + '-' + fields[2]
            moleculeType = getMoleculeType(fields[3])
            # if unitType is not specified, try reading it as RNA, then as protein, ...
            for mt in moleculeTypes:
                if moleculeType == '' or moleculeType == mt:
                    if not chainString in chainData:         # only load once
                        if mt == "RNA" or mt == "DNA":
                            Q, centers, rotations, ids, id_to_index, index_to_id, chainIndices = readNAPositionsFile(Q,chainString,0)
                        elif mt == "protein":
                            Q, centers, ids, id_to_index, index_to_id, chainIndices = readProteinPositionsFile(Q,fields[0],0)

                        chainData[chainString] = [centers,rotations,id_to_index,index_to_id]

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

def synonymAlernateConstraint(synonym,alternate,constraint):

    return constraints

def calculateQueryConstraints(Q):
    """
    Parse the constraints in interactionMatrix and convert to the form
    that the search code needs.
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
                if(iM != None and len(iM) > 0):

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

        #RNAUnitTypes = ["A","C","G","U","M","R","W","S","Y","K","V","H","D","B","N"]
        #DNAUnitTypes = ["DA","DC","DG","DT"]
        #proteinUnitTypes = ["ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE","LEU","LYS","MET","PHE","PRO","PYL","SER","SEC","THR","TRP","TYR","VAL","ASX","GLX","XAA","XLE"]

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

            if(Q["interactionMatrix"][i][i] == None or len(Q["interactionMatrix"][i][i]) == 0):
                if('RNA' not in Q["requiredMoleculeType"][i]):
                    Q["requiredMoleculeType"][i].append('RNA')

            else:
                iMtext = Q["interactionMatrix"][i][i].replace(","," ")

                iMarray = iMtext.split(" ")
                for iM in iMarray:
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


                    elif iM in RNAUnitTypes:
                        Q["requiredUnitType"][i].append(iM)
                        if('RNA' not in Q["requiredMoleculeType"][i]):
                            Q["requiredMoleculeType"][i].append('RNA')

                    elif iM in letterToConstraint:
                        Q["requiredUnitType"][i] += letterToConstraint[iM]
                        if('RNA' not in Q["requiredMoleculeType"][i]):
                            Q["requiredMoleculeType"][i].append('RNA')

                    elif(iM in DNAUnitTypes):
                        Q["requiredUnitType"][i].append(iM)
                        if('DNA' not in Q["requiredMoleculeType"][i]):
                            Q["requiredMoleculeType"][i].append('DNA')

                    elif(iM in proteinUnitTypes):
                        Q["requiredUnitType"][i].append(iM)
                        if('protein' not in Q["requiredMoleculeType"][i]):
                            Q["requiredMoleculeType"][i].append('protein')

                    elif(iM == "RNA"):
                        if('RNA' not in Q["requiredMoleculeType"][i]):
                            Q["requiredMoleculeType"][i].append('RNA')

                    elif(iM == "DNA"):
                        if('DNA' not in Q["requiredMoleculeType"][i]):
                            Q["requiredMoleculeType"][i].append('DNA')

                    elif(iM == "protein"):
                        if('protein' not in Q["requiredMoleculeType"][i]):
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
    unit_data_file = []     # list of files stored in DATAPATH + '/units'
    for search_file in Q["searchFiles"]:

        search_file = search_file.strip()

        if len(search_file) == 0:
            continue

        # look up representative sets, if requested, and replace with IFE names
        if "nrlist" in search_file:           # referring to lists that are posted online
            listLoaded = False
            pathAndFileName = os.path.join(DATAPATH,'units','Representative_sets.pickle')

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
                                if len(IFE) >= 4:
                                    if not "PDB_data_file" in Q:
                                        Q["PDB_data_file"] = readPDBDatafile()  # available PDB structures, resolutions, chains

                                    pdb_id = IFE[0:4]
                                    if pdb_id in list(Q["PDB_data_file"].keys()):
                                        if 'method' in list(Q["PDB_data_file"][pdb_id].keys()):
                                            if 'NMR' in Q["PDB_data_file"][pdb_id]['method']:
                                                newList.append(IFE.decode("ascii"))
                            elif len(IFE) > 1:
                                newList.append(IFE.decode("ascii"))

                    IFEList.extend(newList)

                    if os.path.exists(pathAndFileName):
                        with open(pathAndFileName, 'rb') as fh:
                            representativeSets = pickle.load(fh)
                    else:
                        representativeSets = {}

                    representativeSets[search_file] = newList
                    pickle.dump(representativeSets, open(pathAndFileName, "wb" ), 2)

        elif "|" in search_file:
            # a chain or IFE is specified
            chains = search_file.split("+")
            for k in range(0,len(chains)):
                fields = chains[k].split("|")
                # convert PDB ids to uppercase to avoid case problems
                fields[0] = fields[0].upper()
                chains[k] = "|".join(fields)
            IFEList.append("+".join(chains))

        elif len(search_file) == 4 and not "|" in search_file and not "." in search_file:
            # process as a PDB id
            # if RNA or DNA is being searched for, convert any 4-letter PDB IDs to strings with chain strings separated by + signs
            # which is how RNA chains are stored in .pickle files for speed
            # not needed for protein because they are not stored by chain

            # load mapping from PDB id to chain and other information, if not already done
            if not "PDB_data_file" in Q:
                Q["PDB_data_file"] = readPDBDatafile()  # available PDB structures, resolutions, chains

            search_file = search_file.upper()

            chains = []

            if search_file in Q["PDB_data_file"].keys():

                #print(Q["PDB_data_file"][search_file])

                if searchingRNA and 'RNA' in Q["PDB_data_file"][search_file]['chains']:
                    chains += ["|".join([search_file,'1',x]) for x in Q["PDB_data_file"][search_file]['chains']['RNA']]
                if searchingDNA and 'DNA' in Q["PDB_data_file"][search_file]['chains']:
                    chains += ["|".join([search_file,'1',x]) for x in Q["PDB_data_file"][search_file]['chains']['DNA']]
                IFEList.append("+".join(chains))
                # print("Found these chains %s in %s" % (chains,file))
            else:
                print("Could not find %s in units/NA_datafile.pickle" % search_file)

                # read directory listing for units directory, get filenames that include search_file in the name
                units_path = os.path.join(DATAPATH,'units')
                if os.path.exists(units_path):
                    # get directory listing
                    file_list = os.listdir(units_path)
                    # Loop over the filenames
                    for filename in file_list:
                        # Check if the search_file string is in the filename
                        print(filename)
                        if search_file in filename:
                            fields = filename.split("_")
                            print(fields)
                            chains += fields[2]

                    if len(chains) == 0:
                        # no pre-computed pairwise annotations found
                        # attempt to annotate pairs from the .cif file


                else:
                    Q["errorMessage"].append("Data directory %s does not exist" % units_path)




        elif SERVER:
            print('  Unknown file %s' % search_file)
            Q["errorMessage"].append("Did not recognize %s as a file from the Protein Data Bank" % search_file)

        else:
            # process as a user-defined file ... but that may not work

            print('  query_processing: Processing %s as a local file' % search_file)

            # read the units folder to see what chains are available and list those
            if not unit_data_file:
                unit_data_file = readUnitFileNames(Q["PDB_data_file"])

            #print("  query_processing: list of local chain files found:")
            #print(unit_data_file)

            if search_file in unit_data_file.keys():

                print("  query_processing: found the chains in %s in the units directory" % search_file)

                # found the list of chains in the /units directory
                chains = []
                if searchingRNA and 'NA' in unit_data_file[search_file]['chains']:
                    chains += ["|".join([search_file,'1',x]) for x in unit_data_file[search_file]['chains']['NA']]
                if searchingDNA and 'NA' in Q["PDB_data_file"][search_file]['chains']:
                    chains += ["|".join([search_file,'1',x]) for x in unit_data_file[search_file]['chains']['NA']]
                IFEList.append("+".join(chains))

            else:
                # if none are available yet, load the .cif or .pdb file
                from file_reading import processPDBFile
                chains, messages = processPDBFile(search_file)
                IFEList.append("+".join(chains))

    Q["searchFiles"] = [x for x in IFEList if len(x) > 0]


    #set cutoffs for geometric or mixed searches
    if(Q["type"] == "geometric" or Q["type"] == "mixed"):
        cutoff = Q['SSCutoff'] #rename sscutoffs
        for i in range(len(cutoff)):
            if cutoff[i] != 'Inf':
                cutoff[i] = float(cutoff[i])
            else:
                cutoff[i] = float("inf")
        Q["cutoff"] = cutoff

    return Q
