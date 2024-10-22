"""
set paths and parameters for FR3D to use
These can also be set in the query .json file
"""

SERVER = True
DATAPATHUNITS = "/var/www/html/units"
DATAPATHPAIRS = "/var/www/html/pairs"
OUTPUTPATH    = "/var/www/fr3d/app/results"
JSONPATH      = "/var/www/fr3d/app/results"  # where to look for JSON files by default
# CIFPATH       = "/usr/local/pipeline/hub-core/cif-files"  # not needed on the server

MAXTIME = float('inf')
MAXTIME = 20
MAXCANDIDATESHEATMAP = 300
MAXCANDIDATES = 1000
REFRESHTIME = 2

JSLOCATION = './'                              # standard location when running locally
JSLOCATION = 'https://rna.bgsu.edu/rna3dhub/'  # get js files from the rna3dhub server

# default values for query fields
Q = {}
Q['downloadDataFiles'] = False      # on the server, if they don't exist, can't download them
Q["PDBDATAFILEPATH"] = "/var/www/html/"