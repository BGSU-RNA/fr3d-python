# write the list of candidate ids and all against all discrepancies to an .html file

import numpy as np
import os
from collections import defaultdict

def format_resolution(data):

    r = data['resolution']
    try:
        s = "%0.2f" % r
        # remove trailing 0
        if s[-1] == "0":
            s = s[0:(len(s)-1)]
    except:
        if 'NMR' in data['method']:
            s = "NMR"
        else:
            s = "NA"
    return s

def writeHTMLOutput(Q,candidates,allvsallmatrix=np.empty( shape=(0, 0) )):
    """
    Write the list of candidates in an HTML format that also shows
    the coordinate window and a heat map of all-against-all distances.
    """

    pairTypes = ['glycosidicBondOrientation','chiDegree','pairsStacks','coplanar','BPh','BR','sO','crossingNumber']
    pairsToPrint = defaultdict(list)

    # record which of the many possible pairwise interaction columns contain data
    for pairType in pairTypes:
        for candidate in candidates:
            interactions = candidate["interactions"]
            for a in range(0,len(candidate['indices'])):        # first position
                for b in range(0,len(candidate['indices'])):    # second position
                    if (a,b,pairType) in interactions:
                        if not pairType == 'crossingNumber':
                            pairsToPrint[pairType].append((a,b))
                        elif len(interactions[(a,b,pairType)]) > 0 \
                         and interactions[(a,b,pairType)][0] != '0' \
                         and interactions[(a,b,pairType)][0] != 'None':
                            pairsToPrint[pairType].append((a,b))

    if Q.get('name',None):
        pagetitle = "FR3D %s" % Q['name']
    elif Q.get('description',None):
        pagetitle = "FR3D %s" % Q['description']
    else:
        pagetitle = "FR3D search results"

    candidatelist = '<table id="instances"; style="white-space:nowrap;">\n'

    numPositions = Q["numPositions"]

    # write header line, with instructions about how to sort each column
    candidatelist += "<tr><th onclick=\"sortTable(0,\'instances\',\'numeric\')\">S.</th><th onclick=\"sortTable(1,\'instances\',\'checkbox\')\">Show</th>"
    header_column = 2
    if Q["type"] == "geometric" or Q["type"] == "mixed":
        candidatelist += "<th onclick=\"sortTable(2,\'instances\',\'numeric\')\">Discrepancy</th>"
        header_column = 3
    candidatelist += "<th onclick=\"sortTable(%d,\'instances\',\'numeric\')\">Res. &#8491;</th>" % header_column # or try  &#x212B;
    header_column += 1
    for j in range(0,numPositions):
        candidatelist += "<th onclick=\"sortTable(%d,\'instances\',\'alpha\')\">Position %d</th>" % (header_column,j+1)
        header_column += 1
    candidatelist += "<th onclick=\"sortTable(%d,\'instances\',\'alpha\')\">Sequence</th>" % header_column
    header_column += 1
    sequence_column = header_column
    for pairType in pairTypes:
        for c in sorted(list(set(pairsToPrint[pairType]))):
            if c[0] < c[1] or (c[0] != c[1] and pairType in ["BPh","BR","sO"]) or (c[0]==c[1] and pairType in ["glycosidicBondOrientation","chiDegree"]):
                if c[0] == c[1] and pairType == 'glycosidicBondOrientation':
                    candidatelist += "<th onclick=\"sortTable(%d,\'instances\',\'alpha\')\">Conf. %d</th>" % (header_column,c[0]+1)
                    header_column += 1
                elif c[0] == c[1] and pairType == 'chiDegree':
                    candidatelist += "<th onclick=\"sortTable(%d,\'instances\',\'numeric\')\">Chi %d</th>" % (header_column,c[0]+1)
                    header_column += 1
                elif pairType == 'crossingNumber':
                    candidatelist += "<th onclick=\"sortTable(%d,\'instances\',\'alpha\')\">Crossing %d--%d</th>" % (header_column,c[0]+1,c[1]+1)
                    header_column += 1
                else:
                    candidatelist += "<th onclick=\"sortTable(%d,\'instances\',\'alpha\')\">%d--%d</td>" % (header_column,c[0]+1,c[1]+1)
                    header_column += 1
    candidatelist += "</tr>\n"

    # write one row for each candidate
    for i in range(0,len(candidates)):
        candidate = candidates[i]
        candidatelist += '<tr><td>'+str(i+1)+'.</td><td><label><input type="checkbox" id="'+str(i)+'" class="jmolInline" data-coord="'
        for j in range(0,numPositions):
            candidatelist += candidate["unitids"][j]
            if j < numPositions-1:
                candidatelist += ','
        candidatelist += '">&nbsp</label></td>'
        if(Q["type"] == "geometric" or Q["type"] == "mixed"):
            candidatelist += "<td>%0.4f" % candidate["discrepancy"] + "</td>"

        file_id = candidate["unitids"][0].split("|")[0]
        if "PDB_data_file" in Q and file_id in Q["PDB_data_file"]:
            candidatelist += "<td>%s</td>" % format_resolution(Q["PDB_data_file"][file_id])
        else:
            candidatelist += "<td>NA</td>"

        # write unit ids
        for j in range(0,numPositions):
            candidatelist += "<td>"+candidate["unitids"][j]+"</td>"

        # write nucleotide sequence with helpful separator symbols
        sequence = ""
        for j in range(0,numPositions):
            fields = candidate["unitids"][j].split("|")
            sequence += fields[3]

            if len(fields) == 9:     # symmetry operator present
                symm = fields[8]
            else:
                symm = ""

            if j+1 < numPositions:
                # if same chain and same symmetry operator, if present
                cfields = candidate["unitids"][j+1].split("|")
                if len(cfields) == 9:     # symmetry operator present
                    csymm = cfields[8]
                else:
                    csymm = ""

                if fields[2] == cfields[2] and symm == csymm:
                    if candidate['chainindices'][j] + 1 == candidate['chainindices'][j+1]: # successive
                        sequence += "-"
                    elif candidate['chainindices'][j] + 1 < candidate['chainindices'][j+1]: # later
                        sequence += "&#8594;"
                    else:
                        sequence += "&#8592;"
                else:
                    sequence += "."
        candidatelist += "<td>"+sequence+"</td>"

        # write interactions by group
        interactions = candidate["interactions"]
        for pairType in pairTypes:
            for c in sorted(list(set(pairsToPrint[pairType]))):
                if c[0] < c[1] or (c[0] != c[1] and pairType in ["BPh","BR","sO"]) or (c[0]==c[1] and pairType in ["glycosidicBondOrientation","chiDegree"]):
                    if (c[0],c[1],pairType) in interactions and len(interactions[(c[0],c[1],pairType)]) > 0 and interactions[(c[0],c[1],pairType)][0] != "None":
                        # multiple interactions of the same type separated by commas
                        candidatelist += "<td>"+",".join(interactions[(c[0],c[1],pairType)])+"</td>"
                    else:
                        candidatelist += "<td></td>"

        candidatelist += '</tr>\n'
    candidatelist += '</table>\n'

    discrepancydata = "var data =  []\n"

    # print("  write_output: ", allvsallmatrix)

    if np.size(allvsallmatrix) > 0:
        # write discrepancy data in new 2022 list format
        # first element is a reference to the div in which the heatmap should appear
        discrepancydata = 'var data = ["#heatmap",[\n'            # start a list, start a matrix

        # second element is a matrix with the numerical values of the discrepancy
        # writing both upper and lower triangles of the matrix
        s = allvsallmatrix.shape[0]
        for c in range(0,s):
            discrepancydata += '['     # start a row of the discrepancy matrix
            ife1 = candidates[c]["unitids"][0]
            for d in range(0,s):
                ife2 = candidates[d]["unitids"][0]
                discrepancydata += "%.4f" % allvsallmatrix[c][d]  # one entry
                if d < s-1:
                    discrepancydata += ','  # commas between entries in a row
                else:
                    discrepancydata += '],\n'  # end a row, newline

        discrepancydata += '],\n'           # end the matrix, continue the list

        # third element is a list of labels of instances
        discrepancydata += '['              # start list of instances
        for c in range(0,s):
            ife1 = candidates[c]["unitids"][0]
            discrepancydata += '"' + ife1 + '"'    # write one instance name in quotes
            if c < s-1:
                discrepancydata += ","  # commas between instances
            else:
                discrepancydata += "]\n]" # end list of instances, end list of data

    # find the location of the template for the HTML file
    if "TEMPLATEFILENAME" in Q:
        template_filename = Q["TEMPLATEFILENAME"]
    else:
        # get the path of the current program
        current_path,current_program = os.path.split(os.path.abspath(__file__))
        template_filename = os.path.join(current_path,'template.html')

    # read template.html into one string
    with open(template_filename, 'r') as myfile:
        template = myfile.read()

    # replace ###PAGETITLE### with pagetitle
    template = template.replace("###PAGETITLE###",pagetitle)
    template = template.replace("###sequencecolumn###",str(sequence_column))

    if Q.get('name',None):
        name_text = "Query name: " + Q['name']
    elif Q.get('description',None):
        name_text = "Query description: " + Q['description']
    else:
        name_text = "FR3D search results"

    if len(candidates) == 1:
        s_text = ""
    else:
        s_text = "s"

    queryNote = "%s.<br>Found %d candidate%s from %d of %d files in %0.0f seconds." % (name_text,len(candidates),s_text,Q["numFilesSearched"],len(Q["searchFiles"]),Q["elapsedClockTime"])

    if "moreCandidatesThanHeatMap" in Q:
        queryNote += " " + Q["moreCandidatesThanHeatMap"]

    template = template.replace("###QUERYNAME###",str(queryNote))

    seeModifyQuery = Q.get("seeModifyQuery","")

    template = template.replace("###SEEMODIFYQUERY###",seeModifyQuery)

    if "CSVFILENAME" in Q and not Q.get("reloadOutputPage",False):
        csvlink = os.path.split(Q["CSVFILENAME"])[1]
        seeCSVOutput = '<a href="%s">See CSV output</a>' % csvlink
    else:
        seeCSVOutput = ""

    template = template.replace("###seeCSVOutput###",seeCSVOutput)

    description = """
    <br>Columns of the table show candidate number in geometric similarity order,
    checkbox to display coordinates or not, structure resolution,
    discrepancy from query in geometric or mixed searches,
    unit ids matching each position in the query,
    sequence of the units and backbone connectivity,
    glycosidic bond conformation if requested,
    basepair and stack interactions present, base-phosphate interactions,
    base-ribose interactions, oxygen stacking interactions,
    and the number of nested AU, GC, GU Watson-Crick pairs crossed by each annotated interaction.
    Click on a column header to sort by that column.
    Click on the diagonal of the heatmap to select an instance; right click to de-select.
    Click above the diagonal to select a range of instances, below the diagonal to compare two instances.
    """

    template = template.replace("###DESCRIPTION###",description)

    # replace ###CANDIDATELIST### with candidatelist
    template = template.replace("###CANDIDATELIST###",candidatelist)

    # identify the location of js folder for interactivity, coordinates
    JSLOCATION = '.'           # js folder is in HTML directory
    if Q.get('JSLOCATION',None):
        JSLOCATION = Q['JSLOCATION']
    else:
        try:
            from fr3d.search.fr3d_configuration import JSLOCATION
        except:
            pass


    JS2 = '<script src="%sjs/jquery.jmolTools.WebFR3D.js"></script>' % JSLOCATION
    if JSLOCATION == "./":
        JS1 = '<script src="%sjs/JSmol.min.nojq.js"></script>' % JSLOCATION
        JS3 = '<script src="%sjs/imagehandlinglocal.js"></script>' % JSLOCATION
    else:
        JS1 = '<script src="%sjs/jsmol/JSmol.min.nojq.js"></script>' % JSLOCATION
        JS3 = '<script src="%sjs/imagehandling.js"></script>' % JSLOCATION
    JS4 = '<script src="%sjs/jmolplugin.js" type="text/javascript"></script>' % JSLOCATION
    JS5 = '<script src="%sjs/heatmap_2024.js" type="text/javascript"></script>' % JSLOCATION
    JS6 = '<script src="%sjs/sort.table.js" type="text/javascript"></script>' % JSLOCATION
    JS7 = '<script src="%sjs/shift.click.checkbox.js" type="text/javascript"></script>' % JSLOCATION

    template = template.replace("###JS1###",JS1)
    template = template.replace("###JS2###",JS2)
    template = template.replace("###JS3###",JS3)
    template = template.replace("###JS4###",JS4)
    template = template.replace("###JS6###",JS6)
    template = template.replace("###JS7###",JS7)

    if "REFRESHTIME" in Q:
        REFRESHTIME = Q["REFRESHTIME"]
    else:
        try:
            from fr3d.search.fr3d_configuration import REFRESHTIME
        except:
            REFRESHTIME = 2

    refresh = ""
    if "reloadOutputPage" in Q and Q["reloadOutputPage"]:
        refresh = '<meta http-equiv="refresh" content="%d">' % REFRESHTIME
    template = template.replace("###REFRESH###",refresh)

    if np.size(allvsallmatrix) > 0:
        template = template.replace("###JS5###",JS5)    # include heatmap.js code
        discrepancydata = '<script type="text/javascript">\n' + discrepancydata + '\n</script>'
        template = template.replace("###DISCREPANCYDATA###",discrepancydata)
    else:
        template = template.replace("###DISCREPANCYDATA###","")
        template = template.replace("###JS5###","")    # do not display a heat map

    messages = ""
    messages += "<br>"
    if len(Q["errorMessage"]) > 0:
        messages += "Error messages:<br>"
        for line in Q["errorMessage"]:
            messages += line + "<br>"
    if len(Q["userMessage"]) > 0:
        messages += "User messages:<br>"
        for line in Q["userMessage"]:
            messages += line + "<br>"
    else:
        messages += "No error or warning messages.<br>\n"

    template = template.replace("###MESSAGES###",messages)

    if Q.get('printFileOperations',False):
        print("Writing to %s" % Q['HTMLFILENAME'])

    # could write directly to .gz when desired
    with open(Q['HTMLFILENAME'], 'w') as myfile:
        myfile.write(template)

    if "gzip" in Q and Q["gzip"]:
        os.system("rm -f %s.gz" % Q['HTMLFILENAME'])
        os.system("gzip %s" % Q['HTMLFILENAME'])


def writeCSVOutput(Q,candidates):
    """
    Write the list of candidates in comma separated value format
    """

    pairTypes = ['glycosidicBondOrientation','chiDegree','pairsStacks','BPh','BR','sO','crossingNumber']
    pairsToPrint = defaultdict(list)

    # record which of the many possible columns contain data
    for pairType in pairTypes:
        for candidate in candidates:
            interactions = candidate["interactions"]
            for a in range(0,len(candidate['indices'])):        # first position
                for b in range(0,len(candidate['indices'])):    # second position
                    if (a,b,pairType) in interactions:
                        pairsToPrint[pairType].append((a,b))

    candidatelist = ''

    numPositions = Q["numPositions"]

    # write header line
    candidatelist += "Similarity order,"
    if(Q["type"] == "geometric" or Q["type"] == "mixed"):
        candidatelist += "Discrepancy,"
    candidatelist += "Resolution,"

    for j in range(0,numPositions):
        candidatelist += "Position %d," % (j+1)
    candidatelist += "Sequence,"  # list sequence of candidate
    for pairType in pairTypes:
        for c in sorted(list(set(pairsToPrint[pairType]))):
            if c[0] < c[1] or (c[0] != c[1] and pairType in ["BPh","BR","sO"]) or (c[0]==c[1] and pairType in ["glycosidicBondOrientation","chiDegree"]):
                if c[0] == c[1] and pairType == 'glycosidicBondOrientation':
                    candidatelist += "Orient "+str(c[0]+1)+","
                elif c[0] == c[1] and pairType == 'chiDegree':
                    candidatelist += "Chi "+str(c[0]+1)+","
                elif pairType == 'crossingNumber':
                    candidatelist += "Cross "+str(c[0]+1)+"--"+str(c[1]+1)+","
                else:
                    candidatelist += str(c[0]+1)+"--"+str(c[1]+1)+","

    candidatelist += "View,Coordinates,Sequence variability"  # link to view, link for coordinates

    candidatelist += "\n"

    # write one row for each candidate
    for i in range(0,len(candidates)):
        candidate = candidates[i]
        candidatelist += str(i+1)+','
        if(Q["type"] == "geometric" or Q["type"] == "mixed"):
            candidatelist += "%0.4f," % candidate["discrepancy"]

        file_id = candidate["unitids"][0].split("|")[0]
        if "PDB_data_file" in Q and file_id in Q["PDB_data_file"]:
            candidatelist += format_resolution(Q["PDB_data_file"][file_id]) + ","
        else:
            candidatelist += "NA,"

        # write unit ids
        for j in range(0,numPositions):
            candidatelist += candidate["unitids"][j]+","

        # write nucleotide sequence with helpful separator symbols
        sequence = ""
        for j in range(0,numPositions):
            fields = candidate["unitids"][j].split("|")
            sequence += fields[3]

            if len(fields) == 9:     # symmetry operator present
                symm = fields[8]
            else:
                symm = ""

            if j+1 < numPositions:
                # if same chain and same symmetry operator, if present
                cfields = candidate["unitids"][j+1].split("|")
                if len(cfields) == 9:     # symmetry operator present
                    csymm = cfields[8]
                else:
                    csymm = ""

                if fields[2] == cfields[2] and symm == csymm:
                    if candidate['chainindices'][j] + 1 == candidate['chainindices'][j+1]: # successive
                        sequence += "--"
                    elif candidate['chainindices'][j] + 1 < candidate['chainindices'][j+1]: # later
                        sequence += "->"
                    else:
                        sequence += "<-"
                else:
                    sequence += ".."

        candidatelist += sequence + ","

        # write interactions by group
        interactions = candidate["interactions"]
        for pairType in pairTypes:
            for c in sorted(list(set(pairsToPrint[pairType]))):
                if c[0] < c[1] or (c[0] != c[1] and pairType in ["BPh","BR","sO"]) or (c[0]==c[1] and pairType in ["glycosidicBondOrientation","chiDegree"]):
                    if (c[0],c[1],pairType) in interactions and len(interactions[(c[0],c[1],pairType)]) > 0 and interactions[(c[0],c[1],pairType)][0] != "None":
                        # multiple interactions of the same type separated by commas
                        candidatelist += '"'+",".join(interactions[(c[0],c[1],pairType)])+'",'
                    else:
                        candidatelist += ","

        # make list of unit ids
        unit_id_list = ''
        for j in range(0,numPositions):
            unit_id_list += candidate["unitids"][j]
            if j < numPositions-1:
                unit_id_list += ','

        # make link to view
        if "PDB_data_file" in Q and file_id in Q["PDB_data_file"]:
            candidatelist += '"https://rna.bgsu.edu/rna3dhub/display3D/unitid/' + unit_id_list + '",'
        else:
            candidatelist += ','

        # make link to coordinates
        if "PDB_data_file" in Q and file_id in Q["PDB_data_file"]:
            candidatelist += '"https://rna.bgsu.edu/rna3dhub/rest/getCoordinates?coord=' + unit_id_list + '",'
        else:
            candidatelist += ','

        # make link to sequence variability server
        if "PDB_data_file" in Q and file_id in Q["PDB_data_file"]:
            candidatelist += '"https://rna.bgsu.edu/correspondence/variability?id=' + unit_id_list + '&format=unique",'
        else:
            candidatelist += ','

        candidatelist += '\n'

    if Q.get('printFileOperations',False):
        print("Writing to %s" % Q["CSVFILENAME"])

    with open(Q["CSVFILENAME"], 'w') as myfile:
        myfile.write(candidatelist)


# def writeCandidateOutput(candidates, Q, ifedata):
#     queryName = Q['name']
#     fileName =   '../output/' + queryName.encode('ascii','ignore') + '_python_output.txt'
#     file = open(fileName, 'w')
#     file.write("Found " + str(len(candidates)) + " candidates\n")

#     for candidate in candidates:
#         dataLine = ""

#         indices = ""
#         for index in candidate['indices']:
#             indices += "%6s" % str(index)

#         dataLine += indices

#         file.write(dataLine + '\n')
#     file.close()


