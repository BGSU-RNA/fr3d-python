#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Oct 21 21:18:59 2018

@author:william xi
"""

import inspect #for debugging purposes ONLY
local_vars = {}

import numpy as np
import os
import sys
import datetime
from time import time
from myTimer import myTimer
from write_output import writeHTMLOutput
from write_output import writeCSVOutput

from discrepancy import matrix_discrepancy
#from fr3d.geometry.discrepancy import matrix_discrepancy
from orderBySimilarity import optimalLeafOrder
from orderBySimilarity import treePenalizedPathLength
from orderBySimilarity import reorderSymmetricMatrix

from query_definitions import defineUserQuery
from query_processing import retrieveQueryInformation
from query_processing import calculateQueryConstraints
from query_processing import readQueryFromJSON
from ifedata import readPositionsAndInteractions
from file_reading import readUnitAnnotations # potentially no longer needed here
from file_reading import readNAPairsFile
from file_reading import readPDBDatafile

from fr3d_configuration import DATAPATH
from fr3d_configuration import OUTPUTPATH
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


def main(argv):

    global local_vars #FOR DEBUGGING PURPOSES ONLY

    print("FR3D is starting at " + str(datetime.datetime.now()) + " in " + os.getcwd())

    print("SERVER is " + str(SERVER))

    lastWriteTime = cputime()  # CPU time, for pacing the writing of HTML files

    timerData = myTimer("start")

    # ======================================================================
    # check for required directories for data and output

    if not SERVER:
        directory = DATAPATH
        try:
            os.stat(directory)
        except:
            os.mkdir(directory)
            print("Made " + directory + " directory")

        directory = os.path.join(DATAPATH,'units')
        try:
            os.stat(directory)
        except:
            os.mkdir(directory)
            print("Made " + directory + " directory")

        directory = os.path.join(DATAPATH,'pairs')
        try:
            os.stat(directory)
        except:
            os.mkdir(directory)
            print("Made " + directory + " directory")

        directory = OUTPUTPATH
        try:
            os.stat(directory)
        except:
            os.mkdir(directory)
            print("Made " + directory + " directory")

    # ======================================================================

    if len(argv) == 0:
        # load a user-defined query and calculate additional constraints from it

        # use a query defined in a .json file
        queryName = "sarcin3geometric.json"
        queryName = "Query_627ed1213fa26.json"  # 6EK0 has no glycosidic annotations
        queryName = "Query_6280f9d31d16b.json"  # one chi angle was None
        queryName = "Query_6341a36090c96.json"
        queryName = "Query_638fce3674003.json"
        queryName = "Query_62fc83d4094d8.json"

        # set to blank to use user-defined query from query_definitions
        queryName = ""

    else:
        if SERVER:
            queryID = argv[0].replace("Query_", "").replace(".json", "").replace(".m", "")
            queryName = queryID + "/" + "Query_" + queryID + ".json"
        else:
            queryName = argv[0]

    if not type(queryName) is list:
        queryName = [queryName]

    queryNames = queryName

    # run one or more queries
    for queryName in queryNames:

        timerData = myTimer("Loading query")

        if ".json" in queryName:
            Q = readQueryFromJSON(queryName)
        else:
            Q = defineUserQuery(queryName)

        Q["userMessage"] = []
        Q["errorMessage"] = []

        # pass along information to be able to terminate long, slow searches
        Q["MAXTIME"] = MAXTIME
        Q["FR3Dstarttime"] = time()           # clock time when this job started
        Q["CPUTimeUsed"] = 0                  # charge users for searching, not file loading

        if SERVER:
            Q['server'] = True

        # send email so the user can look up the results later
        if SERVER and Q["email"]:

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

        # log the name of the query being run
        if "name" in Q and Q["name"]:
            print("Running query: %s" % Q["name"].encode('utf-8'))
        elif "queryID" in Q and Q["queryID"]:
            print("Running query: %s" % Q["queryID"])
        else:
            print("Running a query with no name or queryID")

        print("Input query:")
        print(Q)

        Q = retrieveQueryInformation(Q)
        if "errorStatus" in Q:
            Q["numFilesSearched"] = 0
            Q["elapsedClockTime"] = 0
            writeHTMLOutput(Q, [], {})
            continue    # go on to the next query in the loop

        Q = calculateQueryConstraints(Q)

        print("Processed version of the query:")
        for key in sorted(Q.keys()):
            if not key == "PDB_data_file" and not key == "searchFiles":
                print("  %s:%s" % (key,Q[key]))

        # prepare a place to store candidates resulting from the search
        candidates = []

        overallStartTime = time()

        print("Searching %d files or chains; first 10 are:" % len(Q["searchFiles"]))
        print(Q["searchFiles"][0:min(10,len(Q["searchFiles"]))])

        numFilesSearched = 0
        # Load the required PDB files
        # for loop to run over the desired search files
        for ifeNum, ifename in enumerate(Q["searchFiles"]):
            numFilesSearched += 1

            timerData = myTimer("File reading")
            IFEStartTime = time()

            if len(ifename) == 0:
                print("IFE name %s has length 0" % ifename)
                continue

            if not SERVER:
                print("Loading %s, file %d of %d" % (ifename, ifeNum + 1, len(Q["searchFiles"])))

            # read RNA, protein, depending on unittype field
            Q, ifedata = readPositionsAndInteractions(Q, ifename)

            # if there is not enough to search, skip the rest of the processing for this IFE
            if len(ifedata["units"]) < Q["numpositions"]:

                print(ifedata["units"])

                if not 'server' in Q:
                    print("%s has only %d units which is not enough for this search" % (ifename,len(ifedata["units"])))
                continue

            # if the query has a constraint such as cWW_exp, for experimental,
            # load the file of those interactions in a different way
            if "alternateInteractions" in Q:
                interactionToPairs = ifedata['interactionToPairs']
                pairToInteractions = ifedata['pairToInteractions']
                pairToCrossingNumber = ifedata['pairToCrossingNumber']

                for alternate in Q["alternateInteractions"]:
                    PDBID = ifename.split("|")[0]
                    (Q, altInteractionToPairs, altPairToInteractions,
                        altPairToCrossingNumber) = readNAPairsFile(Q, PDBID, ifedata["id_to_index"], alternate)

#                    for key in altPairToInteractions.keys():
#                        print("Interaction between %s and %s is %s" % (index_to_id[key[0]],index_to_id[key[1]],altPairToInteractions[key]))

#                    print(pairToInteractions)
#                    print(pairToCrossingNumber)

                    # loop over new pairs of indices
                    for pair in altPairToInteractions.keys():
                        # extend previous list of interactions or start a new one
                        # crossing number is just an integer, not a list
                        # don't replace previous value, even if that was None

                        pairToInteractions[pair].extend(altPairToInteractions[pair])

                        """
                        if pair in pairToInteractions:
                        else:
                            pairToInteractions[pair] = altPairToInteractions[pair]
                        """

                    for pair in altPairToCrossingNumber.keys():
                        if not pair in pairToCrossingNumber:
                            pairToCrossingNumber[pair] = altPairToCrossingNumber[pair]

#                    print(pairToInteractions)
#                    print(pairToCrossingNumber)

#                    print(Q["activeInteractions"])
#                    print("original",interactionToPairs)

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


            # search for candidates that meet all requirements of the query
            Q, newCandidates, CPUtimeelapsed = FR3D_search(Q, ifedata, ifename, timerData)
            candidates += newCandidates
            Q["CPUTimeUsed"] += CPUtimeelapsed

            if not SERVER or len(newCandidates) > 0:
                # process the list of candidates
                if len(newCandidates) == 1:
                    print("Found %d candidate from %s in %0.2f seconds." % (len(newCandidates),ifename,time() - IFEStartTime))
                else:
                    print("Found %d candidates from %s in %0.2f seconds." % (len(newCandidates),ifename,time() - IFEStartTime))

                if not "PDB_data_file" in Q:
                    Q["PDB_data_file"] = readPDBDatafile()  # available PDB structures, resolutions, chains

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
            for i in range(0, ifeNum + 1):
                ifeList += Q["searchFiles"][i] + ","
            Q["userMessage"].append(ifeList + "<br>\n")
            Q["userMessage"].append("Structures not searched:<br>\n")
            ifeList = ""
            for i in range(ifeNum + 1, len(Q["searchFiles"])):
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
            # more than 600 is too large even locally
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
            """
            if not SERVER:
                st = time()
                newOrder = optimalLeafOrder(allvsallmatrix)
                allvsallmatrix1 = reorderSymmetricMatrix(allvsallmatrix, newOrder)
#                print("OLO time        ",time()-st)
                newCandidates = []
                for i in range(0, matrix_dim):
                    newCandidates.append(candidates[newOrder[i]])
                Q["numFilesSearched"] = numFilesSearched
                Q["elapsedClockTime"] = time() - Q["FR3Dstarttime"]
                writeHTMLOutput(Q, newCandidates, allvsallmatrix1, "_OLO")
            """

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

        print(myTimer("summary"))

        local_vars = inspect.currentframe().f_locals

if __name__== "__main__":
  main(sys.argv[1:])         # pass in arguments besides the .py filename
