"""
set paths and parameters for FR3D to use
These can also be set in the query .json file
"""

SERVER = True
DATAPATHUNITS = "/var/www/html/units"
DATAPATHPAIRS = "/var/www/html/pairs"
OUTPUTPATH    = "/var/www/fr3d/app/results"
# CIFPATH       = "/usr/local/pipeline/hub-core/cif-files"  # we should not need this
JSONPATH      = "/var/www/fr3d/app/results"  # where to look for JSON files by default

MAXTIME = float('inf')
MAXTIME = 20
MAXCANDIDATESHEATMAP = 300
MAXCANDIDATES = 1000
REFRESHTIME = 2

JSLOCATION = './'  # js folder is in HTML folder
