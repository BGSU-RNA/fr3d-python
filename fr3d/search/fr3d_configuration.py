import os
# try to determine where FR3D is running, so this file does not need to be modified
if "vtnguye" in os.getcwd():
	SERVER = False
	DATAPATH = "D:/vtnguye/Bowling Green State University/data/"  # where to store downloaded data files
	OUTPUTPATH = "../output/"
	JSONPATH = "../JSONqueries/"
	TEMPLATEPATH = "./"
	MAXTIME = 20
	MAXTIME = float('inf')
	MAXCANDIDATESHEATMAP = 300
	MAXCANDIDATES = 10000
	REFRESHTIME = 20

	JS1 = '  <script src="./js/JSmol.min.nojq.js"></script>'
	JS2 = '  <script src="./js/jquery.jmolTools.js"></script>'
	JS3 = '  <script src="./js/imagehandlinglocal.js"></script>'
	JS4 = '<script src="./js/jmolplugin.js" type="text/javascript"></script>'
	JS5 = '<script type="text/javascript" src="./js/heatmap.js"></script>'

elif "monke" in os.getcwd():
	SERVER = False
	DATAPATH = "D:/Users/monke/OneDrive/Desktop/oldPC/RNA/data"  # where to store downloaded data files
	OUTPUTPATH = "../output/"
	JSONPATH = "../JSONqueries/"
	TEMPLATEPATH = "./"
	MAXTIME = 20
	MAXTIME = float('inf')
	MAXCANDIDATESHEATMAP = 300
	MAXCANDIDATES = 10000
	REFRESHTIME = 20

	JS1 = '  <script src="./js/JSmol.min.nojq.js"></script>'
	JS2 = '  <script src="./js/jquery.jmolTools.js"></script>'
	JS3 = '  <script src="./js/imagehandlinglocal.js"></script>'
	JS4 = '<script src="./js/jmolplugin.js" type="text/javascript"></script>'
	JS5 = '<script type="text/javascript" src="./js/heatmap.js"></script>'

elif "zirbel" in os.getcwd():
	SERVER = False
	CIFPATH = "C:/Users/zirbel/Documents/FR3D/PDBFiles"
	DATAPATH = "C:/Users/zirbel/Documents/FR3D/Python FR3D/data/"
	OUTPUTPATH = "C:/Users/zirbel/Documents/FR3D/Python FR3D/output/"
	JSONPATH = "C:/Users/zirbel/Documents/FR3D/Python FR3D/JSONqueries/"
	TEMPLATEPATH = "./"
	MAXTIME = 20
	MAXTIME = float('inf')
	MAXCANDIDATESHEATMAP = 300
	MAXCANDIDATES = 100000
	REFRESHTIME = 20

	JS1 = '  <script src="./js/JSmol.min.nojq.js"></script>'
	JS2 = '  <script src="./js/jquery.jmolTools.js"></script>'
	JS3 = '  <script src="./js/imagehandlinglocal.js"></script>'
	JS4 = '<script src="./js/jmolplugin.js" type="text/javascript"></script>'
	JS5 = '<script type="text/javascript" src="./js/heatmap.js"></script>'
else:
	SERVER = True
	CIFPATH = "/"
	DATAPATH = "/var/www/html/"
	OUTPUTPATH = "/var/www/web-fr3d/Results/"
	JSONPATH = "/var/www/web-fr3d/Results/"
	TEMPLATEPATH = "/var/www/web-fr3d/python/"
	MAXTIME = 1200             # seconds
	MAXCANDIDATESHEATMAP = 300
	MAXCANDIDATES = 1000
	REFRESHTIME = 2

	JS1 = '  <script src="http://rna.bgsu.edu/rna3dhub/js/jsmol/JSmol.min.nojq.js"></script>'
	JS2 = '  <script src="http://rna.bgsu.edu/rna3dhub/js/jquery.jmolTools.WebFR3D.js"></script>'
	JS3 = '  <script src="http://rna.bgsu.edu/webfr3d/js/imagehandling.js"></script>'
	JS4 = '<script src="http://rna.bgsu.edu/webfr3d/js/jmolplugin.js" type="text/javascript"></script>'
	JS5 = '<script type="text/javascript" src="http://rna.bgsu.edu/webfr3d/js/heatmap.js"></script>'

