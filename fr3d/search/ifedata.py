"""
Read centers and rotations and pairwise interactions for an IFE or entire 3D structure file
"""

import numpy as np

from file_reading import readNAPairsFile
from file_reading import readNAPositionsFile
from file_reading import readProteinPositionsFile
from file_reading import readUnitAnnotations

from query_processing import getMoleculeType


def combine_dicts(x,y):
    z = x.copy()
    z.update(y)
    return z

def readPositionsAndInteractions(Q, ifename, alternate=""):
    """

    """

    fields =  ifename.split('|')

    # field 0 is generally PDB identifier, but could be user-defined instead
    file_id = fields[0]

    # lists should start empty, append with RNA if necessary, with DNA if necessary, with protein if necessary, etc.
    starting_index = 0
    ifedata = {}
    ifedata['index_to_id'] = {}
    ifedata['id_to_index'] = {}
    ifedata['ids'] = []
    ifedata['units'] = []
    ifedata['centers'] = np.empty((0, 3))
    ifedata['models'] = []

    # accumulate the set of all requiredMoleculeType lists
    requiredMoleculeTypes = []
    for index in range(len(Q["requiredMoleculeType"])):
        requiredMoleculeTypes += Q["requiredMoleculeType"][index]
    requiredMoleculeTypes = list(set(requiredMoleculeTypes))

    # check to see if RNA or DNA is a required unit type, and if so, read NA data
    if "RNA" in requiredMoleculeTypes or "DNA" in requiredMoleculeTypes:

        # split chains in an IFE, which are separated by the + character
        for chainString in ifename.split("+"):

            Q, centers, rotations, ids, id_to_index, index_to_id, chainIndices = readNAPositionsFile(Q,chainString,starting_index)

            if len(centers) > 0:
                # the following lines set up separate but parallel lists
                # if you want to insert, or re-order, you need to do that for all lists simultaneously
                # should we instead make a list of units, and each unit is a dictionary with these fields?
                ifedata['index_to_id'] =  combine_dicts(ifedata['index_to_id'], index_to_id)
                ifedata['id_to_index'] =  combine_dicts(ifedata['id_to_index'], id_to_index)
                ifedata['centers'] = np.append(ifedata['centers'], centers, axis = 0) # this is new
                ifedata['ids'].extend(ids)

                for unitID in ids:
                    unit_information = {}
                    unit_information["centers"] = centers[id_to_index[unitID]-starting_index]
                    unit_information["rotations"] = rotations[id_to_index[unitID]-starting_index]
                    data = unitID.split("|")
                    ifedata['models'].append(data[1])
                    unit_information["unitType"] = data[3] # extract out base from unitID

                    # use component type to infer RNA or DNA
                    unit_information["moleculeType"] = getMoleculeType(data[3])

                    unit_information["chainindex"] = chainIndices[id_to_index[unitID]-starting_index]

                    ifedata["units"].append(unit_information) # append center and rotation information for each unit ID

                starting_index += len(centers)

        if "glycosidicBondOrientation" in Q or "chiAngle" in Q or "showGlycosidicBondOrientation" in Q:
            Q, unit_id_to_annotation = readUnitAnnotations(Q,ifename)

            for index in range(0,len(ifedata["units"])):
                unitID = ifedata['index_to_id'][index]
                if unitID in unit_id_to_annotation:
                    gly = unit_id_to_annotation[unitID]['orientation']
                    chi = unit_id_to_annotation[unitID]['chi_degree']
                else:
                    gly = ""
                    chi = None
                    fields = unitID.split("|")
                    while len(fields) < 7:
                        fields.append('')
                    for alternate_id in ['A','B','C','O']:
                        fields[6] = alternate_id
                        newID = "|".join(fields)
                        if not 'server' in Q:
                            print("Trying unit id %s" % newID)
                        if newID in unit_id_to_annotation:
                            gly = unit_id_to_annotation[newID]['orientation']
                            chi = unit_id_to_annotation[newID]['chi_degree']
                            break
                if gly:
                    ifedata["units"][index]["glycosidicBondOrientation"] = gly
                    try:
                        chi = float(chi)
                        chi = round(chi)  # round to nearest integer to match the display
                        ifedata["units"][index]["chiDegree"] = chi
                    except:
                        if not 'server' in Q:
                            Q["errorMessage"].append("No chi angle for %s" % unitID)
                        ifedata["units"][index]["chiDegree"] = None
                else:
                    Q["errorMessage"].append("No glycosidic bond orientation for %s" % unitID)
                    ifedata["units"][index]["glycosidicBondOrientation"] = None
                    ifedata["units"][index]["chiDegree"] = None

    # load NA pair data; even if not part of the query, it's part of the results
    if "RNA" in requiredMoleculeTypes or "DNA" in requiredMoleculeTypes:
        if len(ifedata['units']) > 1:
            Q, interactionToPairs, pairToInteractions, pairToCrossingNumber = readNAPairsFile(Q, file_id, ifedata["id_to_index"], alternate)
            ifedata['interactionToPairs'] = interactionToPairs
            ifedata['pairToInteractions'] = pairToInteractions
            ifedata['pairToCrossingNumber'] = pairToCrossingNumber
        else:
            ifedata['interactionToPairs'] = {}
            ifedata['pairToInteractions'] = {}
            ifedata['pairToCrossingNumber'] = {}

    # also read protein position file if necessary and append to centers, rotations, ids, ...
    if "protein" in requiredMoleculeTypes:
        Q, centers, ids, id_to_index, index_to_id, chainIndices = readProteinPositionsFile(Q, file_id, starting_index)

        if len(ids) > 0:
            ifedata['index_to_id'] =  combine_dicts(ifedata['index_to_id'], index_to_id)
            ifedata['id_to_index'] =  combine_dicts(ifedata['id_to_index'], id_to_index)
            ifedata['centers'] = np.append(ifedata['centers'], centers, axis = 0)
            ifedata['ids'].extend(ids)

            for unitID in ids:
                unit_information = {}
                unit_information["centers"] = centers[id_to_index[unitID]-starting_index]
                unit_information["rotations"] = None  # amino acids do not have a rotation matrix
                data = unitID.split("|")
                ifedata['models'].append(data[1])
                unit_information["unitType"] = data[3] #extract out base from unitID
                unit_information["moleculeType"] = "protein"
                unit_information["chainindex"] = chainIndices[id_to_index[unitID]-starting_index]

                ifedata["units"].append(unit_information) #append center and rotation information for each unit ID
            starting_index += len(centers)

    return Q, ifedata
