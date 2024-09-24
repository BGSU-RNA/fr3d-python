"""
set paths and parameters for FR3D to use
These can also be set in the query .json file
"""

SERVER = False
CIFPATH = "C:/Users/zirbel/Documents/FR3D/PDBFiles"
DATAPATHUNITS = "C:/Users/zirbel/Documents/FR3D/Python FR3D/data/units"
DATAPATHPAIRS = "C:/Users/zirbel/Documents/FR3D/Python FR3D/data/pairs"
OUTPUTPATH = "C:/Users/zirbel/Documents/FR3D/Python FR3D/output/"
JSONPATH = "C:/Users/zirbel/Documents/FR3D/Python FR3D/JSONqueries/"
MAXTIME = 20
MAXTIME = float('inf')
MAXCANDIDATESHEATMAP = 300
MAXCANDIDATES = 10000
REFRESHTIME = 20

# control where the output .html file looks for javascript files; optional.  Also awkward.
JS1 = '  <script src="./js/JSmol.min.nojq.js"></script>'
JS2 = '  <script src="./js/jquery.jmolTools.js"></script>'
JS3 = '  <script src="./js/imagehandlinglocal.js"></script>'
JS4 = '<script src="./js/jmolplugin.js" type="text/javascript"></script>'
JS5 = '<script type="text/javascript" src="./js/heatmap.js"></script>'
