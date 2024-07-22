# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 3:00:00 2021

"""

from collections import defaultdict
from copy import copy
import numpy as np
import sys
from time import time

from discrepancy import matrix_discrepancy_cutoff
from myTimer import myTimer
from query_processing import synonym
from pair_processing import get_pairlist

if sys.version_info[0] < 3:
    from time import clock as cputime  # true cpu time
else:
    from time import process_time as cputime   # placeholder until we figure out how in 3.10


def getPairTypes(interactions):

    pairTypes = []
    for interaction in interactions:
        interactionType = interaction.split("_")[0]  # in case of _exp or similar
        if interactionType in synonym['stack']:
            pairTypes.append('pairsStacks')
        elif interactionType in synonym['pair']:
            pairTypes.append('pairsStacks')
        elif interactionType in synonym['nstack']:
            pairTypes.append('pairsStacks')
        elif interactionType in synonym['npair']:
            pairTypes.append('pairsStacks')
        elif interactionType in synonym['coplanar']:
            pairTypes.append('coplanar')
        elif interactionType in synonym['BPh']:
            pairTypes.append('BPh')
        elif interactionType in synonym['nBPh']:
            pairTypes.append('BPh')
        elif interactionType in synonym['BR']:
            pairTypes.append('BR')
        elif interactionType in synonym['nBR']:
            pairTypes.append('BR')
        elif interactionType in synonym['sO3'] or interactionType in synonym['sO5']:
            pairTypes.append('sO')
        elif interactionType in synonym['nsO3'] or interactionType in synonym['nsO5']:
            pairTypes.append('sO')
        else:
            pairTypes.append('misc')

    return interactions, pairTypes


def lookUpInteractions(Q, indices, pairToInteractions, pairToCrossingNumber, units):
    """
    store the interactions in each candidate according to type,
    so the different types can be output in order
    """

    interactions = defaultdict(list)

    for a in range(0, len(indices)):             # first position
        for b in range(0, len(indices)):         # second position
            pair = (indices[a], indices[b])      # pair of indices in the file

            if pair in pairToInteractions:      # there is an interaction between these units
                # it would be nice to calculate crossing number for all pairs, not just interacting pairs
                inter, pairTypes = getPairTypes(pairToInteractions[pair])

                crossingAppended = False

                for i in range(0,len(inter)):
                    if not inter[i].startswith('cp') or 'cp' in Q["activeInteractions"] or 'cp_exp' in Q["activeInteractions"]:
                        if not crossingAppended:
                            interactions[(a, b, "crossingNumber")].append(str(pairToCrossingNumber[pair]))
                            crossingAppended = True
                        interactions[(a, b, pairTypes[i])].append(inter[i])
                    # print("Missing interaction positions %d,%d crossing %d" % (a,b,pairToCrossingNumber[pair]))
    for a in range(0, len(indices)):             # first position
        if "glycosidicBondOrientation" in units[indices[a]]:
            gly = units[indices[a]]["glycosidicBondOrientation"]
            chi = units[indices[a]]["chiDegree"]
            if len(gly) > 0:
                interactions[(a, a, "glycosidicBondOrientation")].append(gly)
                interactions[(a, a, "chiDegree")].append(str(chi))

    return interactions


def myIntersect(a, b):
    """
    intersect two lists, allowing that one or both could be "full",
    which means that there is no known restriction on the list,
    aside from the elements coming from the current universe of
    possible values
    """

    if("full" in a):
        return b
    elif("full" in b):
        return a
    elif isinstance(a, set):
        return a.intersection(b)
    else:
        return list(set(a).intersection(set(b)))


def mySetIntersect(a, b):
    """
    intersect two sets, allowing that one or both could be "full",
    which means that there is no known restriction on the list,
    aside from the elements coming from the current universe of
    possible values
    """

    if a == "full":
        return b
    elif b == "full":
        return a
    else:
        return a & b


def getDistanceError(Q, units, perm, i, j, a, b):
    """
    calculate how far off the distance between units a and b is,
    compared to distance between i and j in the query
    """

    pair_distance = np.linalg.norm(units[a]["centers"] - units[b]["centers"])
    return (Q["permutedDistance"][i][j] - pair_distance)**2


def makeFullList(universe1, universe2):
    """
    Pair up the two universes, but avoid pairs like (a,a)
    because that would put the same unit into two positions.
    """

    newList = []
    for a in universe1:
        for b in universe2:
            if not a == b:
                newList.append((a, b))

    return newList


def printListLengths(Q, numpositions, universe, listOfPairs, text=""):

    print("Universe sizes followed by list lengths, upper triangle. %s" % text)
    for i in range(0, numpositions):
        line = "%6d" % len(universe[i])
        for j in range(0, numpositions):
            if j < i+1:
                line  = "        " + line
            elif listOfPairs[i][j] == "full":
                line += "    full"
            else:
                line += "%8d" % len(listOfPairs[i][j])
        print(line)


def sameAlternateId(Q, ifedata, possibilities):
    """
    screen to make sure that no possibility has units with different
    alternate id (typically A or B)
    """

    for i in reversed(range(0, len(possibilities))):
        possibility = possibilities[i]
        alternateIdsFound = set()
        for index in possibility:
            if index in ifedata['index_to_id']:
                fields = ifedata['index_to_id'][index].split("|")
                if len(fields) >= 7:
                    if len(fields[6]) > 0:
                        alternateIdsFound.add(fields[6])

        if len(alternateIdsFound) > 1:
            ids = []
            for index in possibility:
                ids.append(ifedata['index_to_id'][index])
            print("Found multiple alternate ids", ids)
            del possibilities[i]

    return possibilities

def symmetryDifference(u1,u2):
    """
    If the units are the same except for the symmetry operator, return the symmetry
    operator of the second nucleotide
    """

    # It's possible that they are the same unit id, then return no difference
    if u1 == u2:
        return ""

    fields2 = u2.split("|")

    # second id must have a symmetry operator, otherwise no match
    if len(fields2) == 9:
        fields1 = u1.split("|")

        # quickly check that they are the same unit
        if fields1[4] == fields2[4] and fields1[3] == fields2[3] and fields1[2] == fields2[2]:
            # should already be the same model but might as well check
            if fields1[1] == fields2[1]:
                if len(fields1) == 5 and fields2[6] == "" and fields2[7] == "":
                    # only difference is the symmetry operator, return that
                    return fields2[8]
                elif len(fields1) == 7 and fields1[6] == fields2[6] and fields2[7] == "":
                    return fields2[8]
                elif len(fields1) == 8 and fields1[6] == fields2[6] and fields1[7] == fields2[7]:
                    return fields2[8]
                elif len(fields1) == 9 and fields1[6] == fields2[6] and fields1[7] == fields2[7] and not fields1[8] == fields2[8]:
                    return fields2[8]

    return ""


def oneSymmetryVersion(Q, ifedata, possibilities, index_to_id):
    """
    Avoid repeating the same candidate with just a different symmetry operator.
    2022-12-06 code takes about 1 minute on 13,000 possibilities, that's too much.

    """

    if len(possibilities) > 0:
        # sort the possibilities according to unit id of the first nucleotide
        possibilities = sorted(possibilities, key=lambda x: index_to_id[x[0]])
        toRemove = []
        n = len(possibilities[0])

        for i in range(0,len(possibilities)-1):
            j = i+1
            if i in toRemove:
                continue
            while j < len(possibilities):
                sd = symmetryDifference(index_to_id[possibilities[i][0]],index_to_id[possibilities[j][0]])
                if len(sd) > 0:
                    print("Matching %s inside %s" % (index_to_id[possibilities[i][0]],index_to_id[possibilities[j][0]]))
                    k = 0
                    sdk = sd
                    while k+1 < n and sdk == sd:
                        k = k + 1
                        sdk = symmetryDifference(index_to_id[possibilities[i][k]],index_to_id[possibilities[j][k]])
                    if k == n-1 and sdk == sd:
                        toRemove.append(j)
                        print("Removing %s" % sd)
                j = j + 1

        if len(toRemove) > 0:
            toKeep = set(list(range(0,len(possibilities)))) - set(toRemove)
            possibilities = [possibilities[i] for i in toKeep]

    return possibilities


def extendFragment(Q, ifedata, perm, currentFragment, secondElementList, possibilityArray,
    numpositions, previousDistanceError=0):
    """
    Starting with a list of m-unit matches which meet the first m pairwise
    constraints, return a list of (m+1)-unit matches that meet the first
    m+1 pairwise constraints.
    """

    numpositions = Q["numpositions"]

    # if the current fragment of a possibility is the full length, return it,
    # nothing more to be added
    if(len(currentFragment) == numpositions):
        return Q, [currentFragment]

    possibilities = []
    n = len(currentFragment) - 1

    # using the last position of currentFragment, reduce the possible lists
    # for the later positions in the motif
    newPossibilityArray = []
    for i in range(n + 1, numpositions):
        if currentFragment[-1] in secondElementList[n][i]:
            choices = mySetIntersect(possibilityArray[i - n],
            # choices for i vertex given current possibilities
            secondElementList[n][i][currentFragment[-1]])
        else:
            choices = set([])

        # if the current list is empty, there is no way to complete the current fragment
        if len(choices) == 0:
            return Q, []

        newPossibilityArray.append(choices)

    # if the next set to be considered is full, there will be way too many matches
    if newPossibilityArray[0] == "full":
        Q["errorMessage"].append("Query is underspecified, halting search. Add more constraints.")
        Q["halt"] = True
        print("Problem: next position to be added is a full list," +
            " which suggests that the query is underspecified")
        return Q, []

    # loop through the choices for the next position, that can be added to the current fragment
    if (Q["type"] == "symbolic"):
        for extension in newPossibilityArray[0]:
            if not extension in currentFragment:       # enforce that no nucleotide can be repeated
                Q, new_poss = extendFragment(Q, ifedata, perm, currentFragment + (extension,),
                    secondElementList, newPossibilityArray, numpositions, 0)
                possibilities += new_poss
    else:
        for extension in newPossibilityArray[0]:
            totalDistanceError = previousDistanceError
            # add to the sum of distance errors, and compare to maximum values for that
            for j in range(0, n):
                totalDistanceError  +=  getDistanceError(Q, ifedata["units"], perm, j, n + 1,
                    currentFragment[j], extension)
            if (totalDistanceError <= Q["cutoff"][n+1]):
                Q, new_poss = extendFragment(Q, ifedata, perm, currentFragment + (extension,),
                    secondElementList, newPossibilityArray, numpositions, totalDistanceError)
                possibilities += new_poss

    return Q, possibilities


def getPossibilities(Q, ifedata, perm, secondElementList, numpositions, listOfPairs):
    """
    intersect lists in listOfPairs to get possibilities which satisfy all pairwise constraints.
    possibilities is a list of
    """

    # An underspecified query with two positions needs special treatment
    if Q['numpositions'] == 2 and listOfPairs[0][1] == 'full':
        Q["errorMessage"].append("Query is underspecified, halting search. Add more constraints.")
        Q["halt"] = True
        print("Problem: next position to be added is a full list," +
            " which suggests that the query is underspecified")
        return Q, []

    possibilities = []
    possibilityArray = []

    # start with all pairs corresponding to positions 0 and 1
    # for each of these, build complete possibilities
    # possibilityArray is a list of lists of still possible positions later in the fragment being built

    CPUEndTime = cputime() + Q['MAXTIME'] - Q["CPUTimeUsed"]

    for firstPair in listOfPairs[0][1]:
        possibilityArray = []
        emptyArray = False
        for i in range(1,numpositions):
            if firstPair[0] in secondElementList[0][i]:
                possibilityArray.append(secondElementList[0][i][firstPair[0]])
            else:
                emptyArray = True
        if Q["type"] == "geometric" or Q["type"] == "mixed":
            distanceError = getDistanceError(Q, ifedata['units'], perm, 0, 1, firstPair[0], firstPair[1])
        else:
            distanceError = 0

        if not emptyArray:
            Q, new_poss = extendFragment(Q, ifedata, perm, firstPair, secondElementList,
                possibilityArray, numpositions,  distanceError)
            possibilities += new_poss

        # truncate searches that are taking too long
        if cputime() > CPUEndTime:
            print('Maximum time exceeded, %d possibilities found' % len(possibilities))
            return Q, possibilities

    return Q, possibilities


def buildSecondElementList(Q, listOfPairs, universe, ifedata):
    """
    turn lists of pairs into lists of second elements of pairs
    """
    # does this need Q/ifedata?

    numpositions = len(listOfPairs) + 1
    secondElementList = [0] * (numpositions)

    for i in range(0, numpositions - 1):
        secondElementList[i] = [0] * (numpositions)
        for j in range(i + 1, numpositions):
            if listOfPairs[i][j] == "full":
                secondElementList[i][j] = defaultdict()
                for c in universe[i]:
                    secondElementList[i][j][c] = "full"
            else:
                secondElementList[i][j] = defaultdict(set)
                for c in listOfPairs[i][j]:
                    secondElementList[i][j][c[0]].add(c[1])

    return secondElementList


def permuteListOfPairs(listOfPairs, perm):
    """return new listOfPairs given a permutation"""

    newListOfPairs = {}
    indexOf = {}
    for i in range(0, len(perm)):
        indexOf[perm[i]] = i

    for i in listOfPairs:
        for j in listOfPairs[i]:
            if indexOf[i] > indexOf[j]:
                if listOfPairs[i][j] == "full":
                    pairs = "full"
                else:
                    pairs = []
                    for pair in listOfPairs[i][j]:
                        pairs.append((pair[1],pair[0]))
                if indexOf[j] in newListOfPairs:
                    newListOfPairs[indexOf[j]][indexOf[i]] = pairs
                else:
                    newListOfPairs[indexOf[j]] = {indexOf[i]:pairs}
            else:
                if indexOf[i] in newListOfPairs:
                    newListOfPairs[indexOf[i]][indexOf[j]] = copy(listOfPairs[i][j])
                else:
                    newListOfPairs[indexOf[i]] = {indexOf[j]:copy(listOfPairs[i][j])}

    return newListOfPairs


def reorderPositions(listOfPairs, universe): #reorder positions
    """
    Reorder positions so that intersecting lists is more efficient
    """

    numpositions = len(listOfPairs) + 1
    listOfPairs_sizes = []
    lengths = np.zeros((numpositions,numpositions))
    for i in listOfPairs:
        for j in listOfPairs[i]:
            if listOfPairs[i][j] == "full":
                L = 1000000000                        # needs to be large but not infinite
                lengths[i][j] = L
                lengths[j][i] = L
            else:
                L = len(listOfPairs[i][j])
                lengths[i][j] = L
                lengths[j][i] = L
            listOfPairs_sizes.append(([i,j], L))
    listOfPairs_sizes = sorted(listOfPairs_sizes, key=lambda x: x[1])
    perm = listOfPairs_sizes[0][0]
    del listOfPairs_sizes
    objects = set(range(0, numpositions))
    objects.remove(perm[0])
    objects.remove(perm[1])
    while len(perm) < numpositions:
        bestLen = float("inf")
        bestCon = -1
        for obj in objects:
            myLen = 1
            for p in perm:
                myLen += lengths[p][obj]
            if myLen < bestLen:
                bestLen = myLen
                bestCon = obj

        perm.append(bestCon)
        objects.remove(bestCon)

    newListOfPairs = permuteListOfPairs(listOfPairs, perm)
    newUniverse = [universe[i] for i in perm]

    return newListOfPairs, newUniverse, perm


def pruneUniversesWithPairs(universe, listOfPairs, positions_and_counts = None):
    """
    loop through the i,j,length triples in positions_and_counts
    from shortest list to longest, identify elements corresponding to i,j,
    and keep only those in universe[i] and universe[j]
    """

    if not positions_and_counts:
        positions_and_counts = []
        numpositions = len(universe)

        for i in range(0, numpositions - 1):
            for j in range(i + 1,numpositions):
                positions_and_counts.append((i,j,len(listOfPairs[i][j])))

    # sort from shortest list to longest to impose the strongest requirement first
    triples = sorted(positions_and_counts, key = lambda x: x[2])

    pairs_seen = set([])

    for i,j,l in triples:
        if not listOfPairs[i][j] == 'full' and not (i,j) in pairs_seen:
            pairs_seen.add((i,j))
            universe_i = set([])
            universe_j = set([])
            for (a,b) in listOfPairs[i][j]:
                universe_i.add(a)
                universe_j.add(b)
            universe[i] = universe[i] & universe_i
            universe[j] = universe[j] & universe_j

    # if a universe is empty, the search has no candidates and we can stop
    emptyUniverse = False
    for i in range(len(universe.keys())):
        if len(universe[i]) == 0:
            emptyUniverse = True

    return universe, emptyUniverse


def prunePairsWithUniverses(universe, listOfPairs, positions_and_counts = None):
    """
    Loop through the i,j,length triples in positions_and_counts to identify
    the positions to focus on; not all positions need to be considered yet.
    Get the easy universes down, then loop over the large listOfPairs entries once.
    Loop over i and j and if one of those if in the list of positions to focus on,
    then reduce listOfPairs[i][j] to only include pairs from universe[i] and universe[j].
    """

    numpositions = len(universe)

    if positions_and_counts:
        focus = set([])
        for i,j,l in positions_and_counts:
            focus.add(i)
            focus.add(j)
    else:
        focus = set(list(range(numpositions)))

    emptyUniverse = False

    # collect lengths of pair lists to focus on
    triples = []
    before_counter = 0
    for i in range(0, numpositions - 1):
        for j in range(i + 1,numpositions):
            if i in focus or j in focus:
                if listOfPairs[i][j] != "full":
                    triples.append((i,j,len(listOfPairs[i][j])))
                    before_counter += len(listOfPairs[i][j])

    # sort from shortest list to longest
    triples = sorted(triples, key = lambda x : x[2])

    # remove pairs with at least one element not in a universe
    # also track and reduce universes at the same time
    after_counter = 0
    for i,j,l in triples:
        newListOfPairs = []
        universe_i = set([])
        universe_j = set([])
        for (a, b) in listOfPairs[i][j]:
            if a in universe[i] and b in universe[j]:
                newListOfPairs.append((a, b))
                universe_i.add(a)
                universe_j.add(b)
        listOfPairs[i][j] = newListOfPairs
        after_counter += len(listOfPairs[i][j])
        universe[i] = universe[i] & universe_i
        universe[j] = universe[j] & universe_j

        if len(universe[i]) == 0 or len(universe[j]) == 0:
            emptyUniverse = True
            return universe, listOfPairs, emptyUniverse, 1.0

    if before_counter > 0:
        reduction = (before_counter-after_counter)/before_counter
    else:
        reduction = 0.0

    return universe, listOfPairs, emptyUniverse, reduction


def collectModelChainSymmetry(ifedata):
    # look up chain and symmetry once for each unit
    modelChainSymmetry = []        # continuity constraints implicitly require same model, chain, symmetry
    chains = []

    for a in range(0, len(ifedata["index_to_id"])):
        fields = ifedata['index_to_id'][a].split("|")
        if len(fields) == 9:
            MCS = fields[1] + "_" + fields[2] + "_" + fields[8]
        else:
            MCS = fields[1] + "_" + fields[2]
        modelChainSymmetry.append(MCS)
        chains.append(fields[2])

    return modelChainSymmetry, chains


def FR3D_search(Q, ifedata, ifename, timerData):

    IFEStartTime = time()
    CPUStartTime = cputime()

    numpositions = Q['numpositions']

    interactionToPairs = ifedata['interactionToPairs']
    pairToInteractions = ifedata['pairToInteractions']
    pairToCrossingNumber = ifedata['pairToCrossingNumber']
    index_to_id = ifedata['index_to_id']
    units = ifedata["units"]
    modelChainSymmetry = []
    chains = []

    # find lists due to required pairwise constraints
    # do this early because an empty list prevents calculating distances in mixed searches
    if "requiredInteractions" in Q and numpositions > 1:
        timerData = myTimer("Required constraints")
        emptyList = False

        requiredListOfPairs = defaultdict(dict)

        for i in range(0, numpositions):
            for j in range(i + 1, numpositions):
                requiredListOfPairs[i][j] = "full"

                # interactions above the diagonal
                if len(Q["requiredInteractions"][i][j]) > 0:
                    newListOfPairs = "full"
                    tempList = []

                    for interaction in Q["requiredInteractions"][i][j]:
                        # apply previous tempList restrictions
                        if interaction == "and":
                            newListOfPairs = myIntersect(newListOfPairs, tempList)
                            tempList = []
                        elif interaction in interactionToPairs and len(interactionToPairs[interaction]) > 0:
                            if not interaction == "bSS" and "crossingNumber" in Q and Q["crossingNumber"][i][j]:
                                for a in range(0,len(interactionToPairs[interaction][1])):
                                    cn = interactionToPairs[interaction][1][a]  # current crossing number
                                    if not cn == None and Q["crossingNumber"][i][j][0] <= cn and (
                                        cn <= Q["crossingNumber"][i][j][1]):
                                        tempList.append(interactionToPairs[interaction][0][a])
                            else:
                                tempList.extend(interactionToPairs[interaction][0])

                    requiredListOfPairs[i][j] = myIntersect(requiredListOfPairs[i][j], newListOfPairs)
                    if len(requiredListOfPairs[i][j]) == 0:
                        emptyList = True

                # below the diagonal, for asymmetric pairs like BPh, BR, sO
                if j in Q["requiredInteractions"] and (
                i in Q["requiredInteractions"][j]) and len(Q["requiredInteractions"][j][i]) > 0:
                    newListOfPairs = []
                    for interaction in Q["requiredInteractions"][j][i]:
                        if interaction in interactionToPairs and len(interactionToPairs[interaction]) > 0:
                            if not interaction == "bSS" and "crossingNumber" in Q and Q["crossingNumber"][i][j]:
                                for a in range(0,len(interactionToPairs[interaction][1])):
                                    cn = interactionToPairs[interaction][1][a]  # current crossing number
                                    if Q["crossingNumber"][i][j][0] <= cn and cn <= Q["crossingNumber"][i][j][1]:
                                        u = interactionToPairs[interaction][0][a][0]
                                        v = interactionToPairs[interaction][0][a][1]
                                        newListOfPairs.append((v, u))
                            else:
                                for a in range(0, len(interactionToPairs[interaction][1])):
                                    u = interactionToPairs[interaction][0][a][0]
                                    v = interactionToPairs[interaction][0][a][1]
                                    newListOfPairs.append((v, u))

                    requiredListOfPairs[i][j] = myIntersect(requiredListOfPairs[i][j], newListOfPairs)
                    if len(requiredListOfPairs[i][j]) == 0:
                        emptyList = True

        if emptyList:
            if Q.get('printRequiredNoMatches',False):
                print("A required constraint has no matches, returning an empty list of candidates")
            Q["CPUTimeUsed"] += cputime()-CPUStartTime
            return Q, [], timerData

    # Impose pairwise distance constraints for geometric and mixed searches
    # This can take 20% or more of the runtime in geometric and mixed searches
    timerData = myTimer("Calculating pairwise distances")
    listOfPairs = get_pairlist(Q, ifedata["models"], ifedata["centers"])

    # define the initial universe for each position in the query; they are sets
    universe = {}
    for i in range(0, numpositions):
        universe[i] = set(ifedata['index_to_id'].keys())

    if Q.get('printListLengths', False):
        printListLengths(Q, numpositions, universe, listOfPairs, "After setting up universes and imposing distance constraints.")

    timerData = myTimer("Unary constraints")

    # use unary constraints to reduce each universe
    if "requiredUnitType" in Q:
        for i in range(0, numpositions):
            if(len(Q["requiredUnitType"][i]) > 0): # nonempty unit type constraint
                temp_universe = set([])
                for index in universe[i]:
                    if(units[index]['unitType'] in Q["requiredUnitType"][i]):
                        temp_universe.add(index)
                universe[i] = universe[i] & temp_universe

    if "requiredMoleculeType" in Q:
        for i in range(0, numpositions):
            if len(Q["requiredMoleculeType"][i]) > 0: #nonempty molecule type constraint
                temp_universe = set([])
                for index in universe[i]:
                    if(units[index]["moleculeType"] in Q["requiredMoleculeType"][i]):
                        temp_universe.add(index)
                universe[i] = universe[i] & temp_universe

    if "requiredInteractions" in Q:
        for i in range(0, numpositions):
            if len(Q["requiredInteractions"][i][i]) > 0: # required interaction constraint on diagonal
                temp_universe = set([])
                for interaction in Q["requiredInteractions"][i][i]:
                    if interaction == "and":
                        pass
                    elif interaction in interactionToPairs:
                        for a,b in interactionToPairs[interaction][0]:
                            temp_universe.add(a)
                universe[i] = universe[i] & temp_universe

    if "prohibitedInteractions" in Q:
        for i in range(0, numpositions):
            if len(Q["prohibitedInteractions"][i][i]) > 0: # nonempty constraint
                for interaction in Q["prohibitedInteractions"][i][i]:
                    if interaction in interactionToPairs:
                        # the following line might be faster with a for loop and .add
                        indices = [a for (a, b) in interactionToPairs[interaction][0]]
                        universe[i] = universe[i] - set(indices)

    if "glycosidicBondOrientation" in Q:
        for i in range(0, numpositions):
            if(len(Q["glycosidicBondOrientation"][i]) > 0): # nonempty orientation constraint
                temp_universe = set([])
                for index in universe[i]:
                    if units[index]["glycosidicBondOrientation"] in Q["glycosidicBondOrientation"][i]:
                        temp_universe.add(index)
                universe[i] = universe[i] & temp_universe

    if "chiAngle" in Q:
        for i in range(0, numpositions):
            if(len(Q["chiAngle"][i]) > 0): # nonempty chi angle constraint
                temp_universe = set([])
                a = Q["chiAngle"][i][1]
                b = Q["chiAngle"][i][2]
                if Q["chiAngle"][i][0] == 'between':
                    for index in universe[i]:
                        if units[index]["chiDegree"]:
                            if a <= units[index]["chiDegree"] and units[index]["chiDegree"] <= b:
                                temp_universe.add(index)
                else:
                    for index in universe[i]:
                        if units[index]["chiDegree"]:
                            if a <= units[index]["chiDegree"] or units[index]["chiDegree"] <= b:
                                temp_universe.add(index)

                universe[i] = universe[i] & temp_universe

    if "chainLength" in Q:
        if len(modelChainSymmetry) == 0:
            modelChainSymmetry, chains = collectModelChainSymmetry(ifedata)
        chain_counter = defaultdict(int)
        for chain in chains:
            chain_counter[chain] += 1
        for i in range(0, numpositions):
            if len(Q["chainLength"][i]) > 0:
                temp_universe = set([])
                a = Q["chainLength"][i][1]
                b = Q["chainLength"][i][2]
                if Q["chainLength"][i][0] == 'between':
                    for index in universe[i]:
                        if a <= chain_counter[chains[index]] and chain_counter[chains[index]] <= b:
                            temp_universe.add(index)
                else:
                    for index in universe[i]:
                        if a <= chain_counter[chains[index]] or chain_counter[chains[index]] <= b:
                            temp_universe.add(index)

                universe[i] = universe[i] & temp_universe

    if Q.get('printListLengths', False):
        printListLengths(Q, numpositions, universe, listOfPairs, "After unary constraints.")

    # if one of the universes is empty, no candidates will be found
    for i in range(len(universe)):
        if len(universe[i]) == 0:
            Q["CPUTimeUsed"] += cputime()-CPUStartTime
            return Q, [], timerData
    emptyUniverse = False

    # reduce lists of pairs according to required pairwise constraints
    if "requiredInteractions" in Q and numpositions > 1:
        timerData = myTimer("Required constraints again")
        positions_and_counts = []        # keep track of where interactions were imposed

        for i in range(0, numpositions):
            for j in range(i + 1, numpositions):
                if not requiredListOfPairs[i][j] == "full":
                    listOfPairs[i][j] = myIntersect(listOfPairs[i][j], requiredListOfPairs[i][j])
                    positions_and_counts.append((i,j,len(listOfPairs[i][j])))

        # propagate each shortened list of pairs to the affected universes, very quick
        timerData = myTimer("Reduce universes after required")
        universe, emptyUniverse = pruneUniversesWithPairs(universe, listOfPairs, positions_and_counts)
        if Q.get('printListLengths', False):
            printListLengths(Q, numpositions, universe, listOfPairs, "After required constraints and reducing their universes.")

        if not emptyUniverse:
            # reduce the other pair lists using the reduced universes from above, slow but effective
            timerData = myTimer("Reduce pairs after required")
            universe, listOfPairs, emptyUniverse, reduction = prunePairsWithUniverses(universe, listOfPairs, positions_and_counts)
            if Q.get('printListLengths', False):
                printListLengths(Q, numpositions, universe, listOfPairs, "After required constraints and reducing pairs from universes.")

    # reduce list of pairs according to continuity constraints
    if not emptyUniverse and "continuityConstraint" in Q and numpositions > 1:
        timerData = myTimer("Continuity constraints")
        positions_and_counts = []        # keep track of where interactions were imposed

        if len(modelChainSymmetry) == 0:
            modelChainSymmetry, chains = collectModelChainSymmetry(ifedata)

        # make sorted universes to be able to walk along chains
        sorted_universe = {}
        for i in range(0, numpositions):
            sorted_universe[i] = sorted(universe[i])

        for i in range(0, numpositions):
            for j in range(i + 1, numpositions):
                if i in Q["continuityConstraint"] and (
                j in Q["continuityConstraint"][i]) and Q["continuityConstraint"][i][j]:
                    constraint = Q["continuityConstraint"][i][j]

                    if constraint[0] == "between":
                        newList = []
                        if listOfPairs[i][j] == "full":
                            for m in range(0, len(universe[i])):
                                a = sorted_universe[i][m]
                                p = ifedata['units'][a]["chainindex"]  # sequence position
                                n = 0
                                b = sorted_universe[j][n]
                                q = ifedata['units'][b]["chainindex"]  # sequence position

                                # probe for the first match
                                while n < len(universe[j]) - 1 and q - p <= constraint[1]:
                                    n += 1
                                    b = sorted_universe[j][n]
                                    q = ifedata['units'][b]["chainindex"]  # sequence position

                                # accumulate matches
                                while n < len(universe[j]):
                                    if q - p > constraint[1] and q - p < constraint[2]:
                                        if a != b and modelChainSymmetry[a] == modelChainSymmetry[b]:
                                            newList.append((a, b))
                                    n += 1
                                    if n < len(universe[j]):
                                        b = sorted_universe[j][n]
                                        q = ifedata['units'][b]["chainindex"]  # sequence position
                        else:
                            for (a, b) in listOfPairs[i][j]:
                                p = ifedata['units'][a]["chainindex"]  # sequence position
                                q = ifedata['units'][b]["chainindex"]  # sequence position

                                if modelChainSymmetry[a] == modelChainSymmetry[b] and (
                                q - p > constraint[1]) and q - p < constraint[2]:
                                    newList.append((a, b))
                        listOfPairs[i][j] = newList
                        positions_and_counts.append((i,j,len(listOfPairs[i][j])))

                    elif constraint[0] == "equal":
                        newList = []
                        if listOfPairs[i][j] == "full":
                            starting_n = 0
                            for m in range(0, len(universe[i])):
                                foundOne = False
                                a = sorted_universe[i][m]
                                p = ifedata['units'][a]["chainindex"]  # sequence position
                                n = starting_n
                                b = sorted_universe[j][n]
                                q = ifedata['units'][b]["chainindex"]  # sequence position
                                # probe for the first match
                                while n < len(universe[j]) - 1 and (
                                modelChainSymmetry[a] != modelChainSymmetry[b] or q - p <= constraint[1]):
                                    n += 1
                                    b = sorted_universe[j][n]
                                    q = ifedata['units'][b]["chainindex"]  # sequence position
                                # accumulate matches
                                while n < len(universe[j]) and (
                                modelChainSymmetry[a] == modelChainSymmetry[b] and q - p < constraint[2]):
                                    if a != b and q-p in constraint[3]:
                                        newList.append((a, b))
                                    if not foundOne:
                                        starting_n = n
                                        foundOne = True
                                    n += 1
                                    if n < len(universe[j]):
                                        b = sorted_universe[j][n]
                                        q = ifedata['units'][b]["chainindex"]  # sequence position
                        else:
                            for (a, b) in listOfPairs[i][j]:
                                p = ifedata['units'][a]["chainindex"]  # sequence position
                                q = ifedata['units'][b]["chainindex"]  # sequence position

                                if modelChainSymmetry[a] == modelChainSymmetry[b] and q - p in constraint[3]:
                                    newList.append((a, b))
                        listOfPairs[i][j] = newList
                        positions_and_counts.append((i,j,len(listOfPairs[i][j])))

                    elif constraint[0] == "outside":
                        newList = []
                        if listOfPairs[i][j] == "full":
                            starting_n = 0
                            for m in range(0, len(universe[i])):
                                foundOne = False
                                a = sorted_universe[i][m]
                                p = ifedata['units'][a]["chainindex"]  # sequence position
                                n = starting_n
                                b = sorted_universe[j][n]
                                q = ifedata['units'][b]["chainindex"]  # sequence position
                                # probe for the first match
                                while n < len(universe[j])-1 and modelChainSymmetry[a] != modelChainSymmetry[b]:
                                    n += 1
                                    b = sorted_universe[j][n]
                                    q = ifedata['units'][b]["chainindex"]  # sequence position
                                # accumulate matches
                                while n < len(universe[j]) and modelChainSymmetry[a] == modelChainSymmetry[b]:
                                    if a != b and (q-p < constraint[1] or q-p > constraint[2]):
                                        newList.append((a,b))
                                    if not foundOne:
                                        starting_n = n
                                        foundOne = True
                                    n += 1
                                    if n < len(universe[j]):
                                        b = sorted_universe[j][n]
                                        q = ifedata['units'][b]["chainindex"]  # sequence position
                        else:
                            for (a, b) in listOfPairs[i][j]:
                                p = ifedata['units'][a]["chainindex"]  # sequence position
                                q = ifedata['units'][b]["chainindex"]  # sequence position

                                if modelChainSymmetry[a] == modelChainSymmetry[b] and (
                                q-p < constraint[1] or q-p > constraint[2]):
                                    newList.append((a,b))
                        listOfPairs[i][j] = newList
                        positions_and_counts.append((i,j,len(listOfPairs[i][j])))

        # propagate each shortened list of pairs to the affected universes, very quick
        timerData = myTimer("Reduce universes after cont'y")
        universe, emptyUniverse = pruneUniversesWithPairs(universe, listOfPairs, positions_and_counts)
        if Q.get('printListLengths', False):
            printListLengths(Q, numpositions, universe, listOfPairs, "After continuity constraints and reducing their universes.")

        if not emptyUniverse:
            # reduce the other pair lists using the reduced universes from above, slow but effective
            timerData = myTimer("Reduce pairs after continuity")
            universe, listOfPairs, emptyUniverse, reduction = prunePairsWithUniverses(universe, listOfPairs, positions_and_counts)
            if Q.get('printListLengths', False):
                printListLengths(Q, numpositions, universe, listOfPairs, "After continuity constraints and reducing pairs from universes.")

    # reduce list of pairs subject to a unit type combination constraint, like CG GC
    if not emptyUniverse and "combinationConstraint" in Q and numpositions > 1:
        timerData = myTimer("Combination constraints")
        positions_and_counts = []        # keep track of where interactions were imposed

        for i in range(0, numpositions):
            for j in range(i + 1, numpositions):
                if len(Q["combinationConstraint"][i][j]) > 0:
                    temp_pair_list = []

                    if listOfPairs[i][j] == "full":
                        listOfPairs[i][j] = makeFullList(universe[i], universe[j])

                    for pair in listOfPairs[i][j]:
                        if((ifedata['units'][pair[0]]['unitType'],
                            ifedata['units'][pair[1]]['unitType']) in Q["combinationConstraint"][i][j]):
                            temp_pair_list.append(pair)

                    listOfPairs[i][j] = temp_pair_list

        # propagate each shortened list of pairs to the affected universes, very quick
        timerData = myTimer("Reduce universes after comb")
        universe, emptyUniverse = pruneUniversesWithPairs(universe, listOfPairs, positions_and_counts)
        if Q.get('printListLengths', False):
            printListLengths(Q, numpositions, universe, listOfPairs, "After combination constraints and reducing universes.")

        if not emptyUniverse:
            # reduce the other pair lists using the reduced universes from above, slow but effective
            timerData = myTimer("Reduce pairs after combination")
            universe, listOfPairs, emptyUniverse, reduction = prunePairsWithUniverses(universe, listOfPairs, positions_and_counts)
            if Q.get('printListLengths', False):
                printListLengths(Q, numpositions, universe, listOfPairs, "After combination constraints and reducing pairs.")


    # reduce list of pairs according to prohibited pairwise constraints
    # generally will not reduce lists by much, so don't prune afterward
    if not emptyUniverse and "prohibitedInteractions" in Q and numpositions > 1:
        timerData = myTimer("Prohibited constraints")
        for i in range(0, numpositions):
            for j in range(i + 1, numpositions):
                # above the diagonal
                if len(Q["prohibitedInteractions"][i][j]) > 0:
                    if listOfPairs[i][j] == "full":
                        listOfPairs[i][j] = makeFullList(universe[i], universe[j])

                    for interaction in Q["prohibitedInteractions"][i][j]:
                        if interaction in interactionToPairs and len(interactionToPairs[interaction]) > 0:
                            listOfPairs[i][j] = list(set(listOfPairs[i][j]) -
                                set(interactionToPairs[interaction][0]))

                # below the diagonal, for asymmetric pairs like BPh, BR, sO
                if (j in Q["prohibitedInteractions"] and (
                    i in Q["prohibitedInteractions"][j]) and (
                    len(Q["prohibitedInteractions"][j][i]))) > 0:
                    if listOfPairs[i][j] == "full":
                        listOfPairs[i][j] = makeFullList(universe[i],universe[j])

                    for interaction in Q["prohibitedInteractions"][j][i]:
                        if interaction in interactionToPairs and len(interactionToPairs[interaction]) > 0:
                            prohibitPairs = [(b,a) for (a,b) in interactionToPairs[interaction][0]]
                            listOfPairs[i][j] = list(set(listOfPairs[i][j]) - set(prohibitPairs))


    if not emptyUniverse:
        # propagate each shortened list of pairs to the affected universes, very quick
        timerData = myTimer("Full universe pruning")
        universe, emptyUniverse = pruneUniversesWithPairs(universe, listOfPairs)
        if Q.get('printListLengths', False):
            printListLengths(Q, numpositions, universe, listOfPairs, "Full universe pruning")

    if not emptyUniverse:
        reduction = 0.5
        i = 1
        while reduction > 0.02 and i < 10:
            # reduce the other pair lists using the reduced universes from above, slow but effective
            timerData = myTimer("Full pair pruning %d" % i)
            universe, listOfPairs, emptyUniverse, reduction = prunePairsWithUniverses(universe, listOfPairs)
            if Q.get('printListLengths', False):
                printListLengths(Q, numpositions, universe, listOfPairs, "Full pair pruning #%d, reduction fraction %0.4f" % (i,reduction))
            i += 1

    # No candidates
    if emptyUniverse:
        Q["CPUTimeUsed"] += cputime()-CPUStartTime
        return Q, [], timerData

    if numpositions == 1:
        # just one position, just one universe, those are the possibilities
        possibilities = [[p] for p in universe[0]]
        inverseperm = [0]

    else:
        # intersect pair lists to find complete candidates

        # reorder positions for more efficient intersections
        listOfPairs, universe, perm = reorderPositions(listOfPairs, universe)
        inverseperm = [perm.index(i) for i in range(len(perm))]

        # compute permuted Distance Matrix
        if (Q["type"] == "geometric" or Q["type"] == "mixed"):
            Q["permutedDistance"] = np.zeros((numpositions, numpositions))
            for i in range(numpositions):
                for j in range(numpositions):
                    Q["permutedDistance"][i][j] = Q["distance"][perm[i]][perm[j]]

            if Q.get('printListLengths', False):
                printListLengths(Q, numpositions, universe, listOfPairs, "After reordering positions")

        timerData = myTimer("Intersecting pair lists")

        # turn listOfPairs into lists of possible second indices
        secondElementList = buildSecondElementList(Q, listOfPairs, universe, ifedata)

        # intersect lists of pairs; this can take a very long time
        Q, possibilities = getPossibilities(Q, ifedata, perm, secondElementList, numpositions, listOfPairs)

    # screen to make sure that no possibility has units with different alternate id (typically A or B)
    timerData = myTimer("Same alternate id")
    possibilities = sameAlternateId(Q, ifedata, possibilities)

    if len(possibilities) > 0 and Q.get('printFoundPossibilities',False):
        print("Found %5d possibilities from %s in %10.4f seconds" % (len(possibilities),ifename,(time() - IFEStartTime)))

    candidates = []

    # for geometric and mixed searches, compute discrepancy of possibilities to query motif
    # for purely symbolic searches, every possibility is a candidate
    if((Q["type"] == "geometric" or Q["type"] == "mixed")):
        querycenters = [Q["centers"][i] for i in perm]
        queryrotations = [Q["rotations"][i] for i in perm]
        timerData = myTimer("Discrepancy from query")

        possibility_to_discrepancy = {}
        for possibility in possibilities:
            possibilitycenters = []
            for i in range(0, numpositions):
                possibilitycenters.append(units[possibility[i]]["centers"])
            possibilityrotations = []
            for i in range(0, numpositions):
                possibilityrotations.append(units[possibility[i]]["rotations"])

            d = matrix_discrepancy_cutoff(querycenters, queryrotations, possibilitycenters,
                possibilityrotations, Q["discrepancy"])

            if d is not None and d < Q["discrepancy"]:
                possibility_to_discrepancy[possibility] = d

        # turns out, it's faster to calculate discrepancies than to compare symmetry operators
        # since the comparison is like O(n^2) and calculating discripancies is like O(n)
        timerData = myTimer("One symmetry version")
        possibilities = oneSymmetryVersion(Q, ifedata, possibility_to_discrepancy.keys(), index_to_id)

        for possibility in possibilities:
            d = possibility_to_discrepancy[possibility]
            newcandidate = {}
            indices = [possibility[inverseperm[i]] for i in range(numpositions)]
            newcandidate['indices'] = indices
            newcandidate['unitids'] = [index_to_id[index] for index in indices]
            newcandidate['chainindices'] = [units[index]["chainindex"] for index in indices]
            newcandidate['centers'] = [units[index]["centers"] for index in indices]
            newcandidate['rotations'] = [units[index]["rotations"] for index in indices]
            newcandidate['discrepancy'] = d
            newcandidate['interactions'] = lookUpInteractions(Q,indices,
                pairToInteractions, pairToCrossingNumber, units)
            candidates.append(newcandidate)

    else:

        # screen to keep only one version of symmetry operated candidates
        timerData = myTimer("One symmetry version")
        possibilities = oneSymmetryVersion(Q, ifedata, possibilities, index_to_id)

        for possibility in possibilities:
            newcandidate = {}
            indices = [possibility[inverseperm[i]] for i in range(numpositions)]
            newcandidate['indices'] = indices
            newcandidate['unitids'] = [index_to_id[index] for index in indices]
            newcandidate['chainindices'] = [units[index]["chainindex"] for index in indices]
            newcandidate['centers'] = [units[index]["centers"] for index in indices]
            newcandidate['rotations'] = [units[index]["rotations"] for index in indices]
            newcandidate['interactions'] = lookUpInteractions(Q,indices,
                pairToInteractions, pairToCrossingNumber, units)
            candidates.append(newcandidate)

    Q["CPUTimeUsed"] += cputime()-CPUStartTime

    return Q, candidates, timerData