#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 21:18:59 2018

@author:william xi, Craig Zirbel
"""

import numpy as np
import os
import sys
import datetime
from time import time
from myTimer import myTimer
from write_output import writeHTMLOutput
from write_output import writeCSVOutput

from discrepancy import matrix_discrepancy  # should replace with line below
#from fr3d.geometry.discrepancy import matrix_discrepancy
from orderBySimilarity import treePenalizedPathLength
from orderBySimilarity import reorderSymmetricMatrix

from query_processing import retrieveQueryInformation
from query_processing import calculateQueryConstraints
from query_processing import readQueryFromJSON
from ifedata import readPositionsAndInteractions
from file_reading import readNAPairsFile
from file_reading import checkDirectories

from fr3d_configuration import SERVER
from fr3d_configuration import MAXTIME
from fr3d_configuration import MAXCANDIDATES
from fr3d_configuration import MAXCANDIDATESHEATMAP
from fr3d_configuration import REFRESHTIME

from search import FR3D_search

if sys.version_info[0] < 3:
    from time import clock as cputime  # true cpu time
else:
    from time import process_time as cputime   # placeholder until we figure out how in 3.10

# provide support to send email
import smtplib
from email.mime.text import MIMEText


def send_email_to_user(Q):

    try:
        # Create a text/plain message
        T = "Thank you for requesting a WebFR3D query.  "
        T += "The results can be found at http://rna.bgsu.edu/webfr3d/Results/"
        T += "%s/%s.html" % (queryID,queryID)
        T += "\n\nBest regards,"
        T += "\n\nThe BGSU RNA group"

        # me == the sender's email address
        # you == the recipient's email address
        if "name" in Q and Q["name"]:
            QName = Q["name"]
            QName = QName[0:80]
            QName = QName.replace("|", "_")
            QName = QName.replace("&", "_")
            QName = QName.replace("^", "_")
            QName = QName.replace("%", "_")
            QName = QName.replace("#", "_")
            QName = QName.replace("@", "_")
            QName = QName.replace("!", "_")
            QName = QName.replace("*", "_")
            QName = QName.replace("(", "_")
            QName = QName.replace(")", "_")
            QName = QName.replace("?", "_")
        else:
            QName = ""

        sender = u'rnabgsu@gmail.com'
        recipient = Q["email"]
        recipient = recipient.replace("\t", ",")
        recipient = recipient.replace(";", ",")
        recipient = recipient.replace(" ", ",")
        recipient = recipient.replace(",,", ",")
        recipient = recipient.replace(",,", ",")
        recipient = recipient.replace(",,", ",")
        recipient = recipient.replace(",,", ",")
        recipient = recipient.replace(",,", ",")
        recipient = recipient.replace(",,", ",")

        print("Email recipient:")
        print(recipient)

        msg = MIMEText(T)
        msg['Subject'] = "WebFR3D query: " + QName
        msg['From'] = sender
        msg['To'] = recipient
        msg['Cc'] = sender

        recipientList = recipient.split(",")
        recipientList.append(sender)

        print("Email recipient list:")
        print(recipientList)

        # Send the message via our own SMTP server.
        s = smtplib.SMTP('localhost')
        s.sendmail(sender, recipientList, msg.as_string())
        s.quit()
    except:
        Q["errorMessage"].append("Could not send email")


def readAlternateInteractions(Q, ifedata, ifename):
    """
    Add additional annotations of interactions ending with _exp for example
    """
    interactionToPairs = ifedata['interactionToPairs']
    pairToInteractions = ifedata['pairToInteractions']
    pairToCrossingNumber = ifedata['pairToCrossingNumber']

    for alternate in Q["alternateInteractions"]:
        PDBID = ifename.split("|")[0]
        (Q, altInteractionToPairs, altPairToInteractions,
            altPairToCrossingNumber) = readNAPairsFile(Q, PDBID, ifedata["id_to_index"], alternate)

        # loop over new pairs of indices
        for pair in altPairToInteractions.keys():
            # extend previous list of interactions or start a new one
            # crossing number is just an integer, not a list
            # don't replace previous value, even if that was None

            pairToInteractions[pair].extend(altPairToInteractions[pair])

        for pair in altPairToCrossingNumber.keys():
            if not pair in pairToCrossingNumber:
                pairToCrossingNumber[pair] = altPairToCrossingNumber[pair]

        # loop over mapping of interactions to pairs
        for interaction in altInteractionToPairs.keys():
            # only store data for interactions needed in the query
            # self interactions are treated separately
            if interaction.replace("_self", "") in Q["activeInteractions"]:
                if interaction in interactionToPairs:
                    listOfPairs, crossingNumber = interactionToPairs[interaction]
                else:
                    listOfPairs = []
                    crossingNumber = []
                newListOfPairs = listOfPairs + altInteractionToPairs[interaction][0]
                newCrossingNumber = crossingNumber + altInteractionToPairs[interaction][1]
                interactionToPairs[interaction] = (newListOfPairs, newCrossingNumber)

    ifedata['interactionToPairs'] = interactionToPairs
    ifedata['pairToInteractions'] = pairToInteractions
    ifedata['pairToCrossingNumber'] = pairToCrossingNumber

    return Q, ifedata


def main(argv):

    if type(argv) is str:
        argv = [argv]

    print("FR3D is starting at " + str(datetime.datetime.now()) + " in " + os.getcwd())

    lastWriteTime = cputime()  # CPU time, for pacing the writing of HTML files
    timerData = myTimer("start")

    if len(argv) == 0:
        # load a user-defined query and calculate additional constraints from it
        print("When running FR3D.py from the command line, specify a query using the name of .json file")

        # set to blank to use user-defined query from query_definitions
        queryNames = []

    else:
        if SERVER:
            queryID = argv[0].replace("Query_", "").replace(".json", "").replace(".m", "")
            queryNames = [queryID + "/" + "Query_" + queryID + ".json"]
        else:
            # split input on commas to allow multiple queries to be run
            queryNames = argv[0].split(",")

    # run one or more queries
    for queryName in queryNames:

        timerData = myTimer("Loading query file")
        print("Loading query file: %s" % queryName)
        Q = readQueryFromJSON(queryName)

        if "errorStatus" in Q:
            Q["numFilesSearched"] = 0
            Q["elapsedClockTime"] = 0
            writeHTMLOutput(Q, [])
            continue    # go on to the next query in the loop

        # print the name of the query being run
        if "name" in Q and Q["name"]:
            print("Running query: %s" % Q["name"])
        elif "queryID" in Q and Q["queryID"]:
            print("Running query: %s" % Q["queryID"])

        if Q.get("printInputQuery", False):
            print("Input query:")
            for key in sorted(Q.keys()):
                if key == "PDB_data_file":
                    continue
                elif key == "searchFiles":
                    L = len(Q[key])
                    if L <= 10:
                        print("  %s:%s" % (key,Q[key]))
                    else:
                        print("  %s:%s ..." % (key,Q[key][0:10]))
                else:
                    print("  %s:%s" % (key,Q[key]))

        # pass along information to be able to terminate long, slow searches
        if not "MAXTIME" in Q:
            Q["MAXTIME"] = MAXTIME

        if not "MAXCANDIDATES" in Q:
            Q["MAXCANDIDATES"] = MAXCANDIDATES

        if not "MAXCANDIDATESHEATMAP" in Q:
            Q["MAXCANDIDATESHEATMAP"] = MAXCANDIDATESHEATMAP

        Q["FR3Dstarttime"] = time()           # clock time when this job started
        Q["CPUTimeUsed"] = 0                  # charge users for searching, not file loading
        Q["userMessage"] = []
        Q["errorMessage"] = []

        if SERVER:
            Q['server'] = True

            # on the server, send email now so the user can look up the results later
            if Q["email"]:
                Q = send_email_to_user(Q)

        # check directories needed for data files and output
        Q = checkDirectories(Q)

        # retrieve information about query nucleotides, if any
        Q = retrieveQueryInformation(Q)
        if "errorStatus" in Q:
            Q["numFilesSearched"] = 0
            Q["elapsedClockTime"] = 0
            writeHTMLOutput(Q, [])
            continue    # go on to the next query in the loop

        # interpret the constraints given in interactionMatrix,
        # set cutoffs for geometric or mixed searches
        Q = calculateQueryConstraints(Q)

        if Q.get("printProcessedQuery", False):
            print("Processed version of the query:")
            for key in sorted(Q.keys()):
                if key == "PDB_data_file":
                    continue
                elif key == "searchFiles":
                    L = len(Q[key])
                    if L <= 10:
                        print("  %s:%s" % (key,Q[key]))
                    else:
                        print("  %s:%s ..." % (key,Q[key][0:10]))
                else:
                    print("  %s:%s" % (key,Q[key]))

        # prepare a place to store candidates resulting from the search
        candidates = []

        overallStartTime = time()


        if Q.get("reverseSearchOrder", False):
            ife_search_list = Q["searchFiles"][::-1]
        else:
            ife_search_list = Q["searchFiles"]

        if Q.get("printSearchFiles", 0) > 0:
            numToShow = min(Q["printSearchFiles"],len(ife_search_list))
            print("Searching %d files or chains; first %d are:" % (len(ife_search_list),numToShow))
            print(ife_search_list[0:numToShow])

        numFilesSearched = 0
        # loop over the 3D structure files or IFEs to be searched
        for ifename in ife_search_list:
            if Q.get("printCumulativeCandidates", False) and numFilesSearched > 0:
                if len(candidates) == 1:
                    print("Found %d candidate from %d of %d files in %0.0f seconds so far." % (len(candidates),numFilesSearched,len(Q["searchFiles"]),time() - overallStartTime))
                else:
                    print("Found %d candidates from %d of %d files in %0.0f seconds so far." % (len(candidates),numFilesSearched,len(Q["searchFiles"]),time() - overallStartTime))

            numFilesSearched += 1

            if len(ifename) == 0:
                print("IFE name %s has length 0" % ifename)
                continue

            if Q.get("printSearchingFile", False):
                print("Searching %s, file %d of %d" % (ifename, numFilesSearched, len(Q["searchFiles"])))

            timerData = myTimer("File reading")
            IFEStartTime = time()

            # read RNA, DNA, protein centers, rotations, and pairwise interactions
            Q, ifedata = readPositionsAndInteractions(Q, ifename)

            # if there are not enough units to search, skip the rest of the processing for this IFE
            if len(ifedata["units"]) < Q["numpositions"]:
                if Q.get("printNotEnoughUnits", False):
                    print("%s has only %d units which is not enough for this search" % (ifename,len(ifedata["units"])))
                continue

            # if the query has a constraint such as cWW_exp, for experimental,
            # load the file of those interactions in a different way
            if "alternateInteractions" in Q:
                Q, ifedata = readAlternateInteractions(Q, ifedata, ifename)

            # search for candidates that meet all requirements of the query
            Q, newCandidates, timerData = FR3D_search(Q, ifedata, ifename, timerData)
            candidates += newCandidates

            if Q.get("printFoundCandidates", False):
                if len(newCandidates) == 1:
                    print("Found %d candidate from %s in %0.2f seconds." % (len(newCandidates),ifename,time() - IFEStartTime))
                else:
                    print("Found %d candidates from %s in %0.2f seconds." % (len(newCandidates),ifename,time() - IFEStartTime))

            # write out a provisional list of candidates if enough time has elapsed,
            # using clock time not CPU time
            if cputime() - lastWriteTime > REFRESHTIME:
                # for geometric or mixed searches, sort candidates by discrepancy from query
                if((Q["type"] == "geometric" or Q["type"] == "mixed")):
                    candidates.sort(key = lambda candidate: candidate["discrepancy"])

                # limit the number of candidates to output
                if len(candidates) > MAXCANDIDATES:
                    print("Found %d candidates but MAXCANDIDATES is %d" % (len(candidates), MAXCANDIDATES))
                    candidates = candidates[:MAXCANDIDATES]
                    Q["hitMaxCandidates"] = True

                # write output with just the candidates, no heat map, for default ordering
                Q["reloadOutputPage"] = True
                Q["numFilesSearched"] = numFilesSearched
                Q["elapsedClockTime"] = time() - Q["FR3Dstarttime"]

                writeHTMLOutput(Q, candidates)
                lastWriteTime = cputime()

            if len(candidates) > MAXCANDIDATES or "hitMaxCandidates" in Q:
                print("Found %d candidates but MAXCANDIDATES is %d" % (len(candidates), MAXCANDIDATES))
                Q["hitMaxCandidates"] = True
                Q["userMessage"].append("Found %d candidates; maximum number to output is %d" %
                    (len(candidates), MAXCANDIDATES))
                break

            if Q["CPUTimeUsed"] > Q["MAXTIME"]:
                print("Used %0.2f CPU seconds for searching but maximum allowed time is %0.0f seconds.  Add symbolic constraints and/or reduce the discrepancy cutoff." % (Q["CPUTimeUsed"],Q["MAXTIME"]))
                Q["hitMaxTime"] = True
                Q["userMessage"].append(
                    "Used %0.2f CPU seconds for searching but maximum time allowed is %0.0f seconds.  Add symbolic constraints and/or reduce the discrepancy cutoff." %
                    (Q["CPUTimeUsed"],Q["MAXTIME"]))
                break

            if 'halt' in Q:
                break

            # end of the loop over all IFEs

        Q["reloadOutputPage"] = False
        Q["numFilesSearched"] = numFilesSearched

        if ("hitMaxTime" in Q and Q["hitMaxTime"]) or ("hitMaxCandidates" in Q and Q["hitMaxCandidates"]):
            Q["userMessage"].append("Structures searched:<br>\n")
            ifeList = ""
            for i in range(0, numFilesSearched):
                ifeList += Q["searchFiles"][i] + ","
            Q["userMessage"].append(ifeList + "<br>\n")
            Q["userMessage"].append("Structures not searched:<br>\n")
            ifeList = ""
            for i in range(numFilesSearched, len(Q["searchFiles"])):
                ifeList += Q["searchFiles"][i] + ","
            Q["userMessage"].append(ifeList + "<br>\n")

        # process the list of candidates
        if len(candidates) == 1:
            print("Found %d candidate from %d files in %0.2f seconds." % (len(candidates),Q["numFilesSearched"],time() - overallStartTime))
        else:
            print("Found %d candidates from %d files in %0.2f seconds." % (len(candidates),Q["numFilesSearched"],time() - overallStartTime))

        # for geometric or mixed searches, sort candidates by discrepancy from query
        if((Q["type"] == "geometric" or Q["type"] == "mixed")):
            candidates.sort(key = lambda candidate: candidate["discrepancy"])

        # limit the number of candidates to output
        if len(candidates) > MAXCANDIDATES:
            candidates = candidates[:MAXCANDIDATES]

        timerData = myTimer("Calculate all vs all matrix")

        if Q["numpositions"] > 1 and len(candidates) > 1:
            # compute all against all discrepancies, up to a certain limit;
            # more than 600 is too large even locally, and too hard to see in a heat map
            matrix_dim = min(MAXCANDIDATESHEATMAP, len(candidates))
            allvsallmatrix = np.zeros((matrix_dim, matrix_dim))

            # TODO: future plan for faster all against all comparisons:
            # for i in range(matrix_dim):
            # subtract mean from candidate1centers and store in a new variable
            # see besttransformation for how to subtract them
            # make a new besttransformation_basic function to replace it
            for i in range(matrix_dim):
                candidate1centers = candidates[i]["centers"]
                candidate1rotations = candidates[i]["rotations"]
                for j in range(i + 1, matrix_dim):
                    candidate2centers = candidates[j]["centers"]
                    candidate2rotations = candidates[j]["rotations"]
                    d = matrix_discrepancy(candidate1centers, candidate1rotations,
                        candidate2centers, candidate2rotations)
                    allvsallmatrix[i][j] = d
                    allvsallmatrix[j][i] = d
        else:
            allvsallmatrix = np.zeros((0, 0))

        # reorder the candidates according to ordering by similarity
        timerData = myTimer("Ordering by similarity")
        if allvsallmatrix.shape[0] > 1:
            newOrder = treePenalizedPathLength(allvsallmatrix, 100, 59)
            allvsallmatrix = reorderSymmetricMatrix(allvsallmatrix, newOrder)
            newCandidates = []
            for i in range(0, matrix_dim):
                newCandidates.append(candidates[newOrder[i]])

            # there may be more candidates to view than we have a heat map for
            if len(candidates) > matrix_dim:
                Q["moreCandidatesThanHeatMap"] = "First %d candidates listed in similarity order and shown in heat map" % matrix_dim
                Q["userMessage"].append("First %d candidates listed in similarity order and shown in heat map" % matrix_dim)
                for i in range(matrix_dim, len(candidates)):
                    newCandidates.append(candidates[i])
        else:
            newCandidates = candidates

        timerData = myTimer("Writing output")
        Q["elapsedClockTime"] = time() - Q["FR3Dstarttime"]
        Q["userMessage"].append("FR3D completed successfully")

        writeHTMLOutput(Q, newCandidates, allvsallmatrix)
        writeCSVOutput(Q, newCandidates)

        if len(Q["errorMessage"]) > 0:
            print("Error message:")
            print(Q["errorMessage"])

        if Q.get("printTimer", False):
            print(myTimer("summary"))

if __name__== "__main__":
  main(sys.argv[1:])         # pass in arguments besides the .py filename
