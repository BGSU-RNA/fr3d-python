# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 3:00:00 2021

@author:adamc
"""

import sys, os
import FR3D
import subprocess
# import io
# from contextlib import redirect_stdout
# from unittest.mock import patch
# from FR3d import main as __name__

# this next line makes a mock function of the builtin `print()`
# @patch('builtins.print') # bailed on idea

# def test(search = ''):
# 	f = io.StringIO() # this seems to only steal my print statement in this file
# 	with redirect_stdout(f):
# 		print(f"########### {search} ###########")
# 		os.system(f'python FR3D.py "{search}"')
# 		# FR3D.main(search) # no good
# 	return(stringy)
# 	# return(f.getvalue())

def test(search = ''):
	if sys.version_info[0] < 3:
		os.system('python FR3D.py "%s"' % search) # this line works in Python 2.7
	else:
		return(subprocess.run('python FR3D.py "%s"' % search)) # this line works in py3.5+
	# next line should be legacy
	# return(subprocess.check_output(f'python FR3D.py "{search}"'))
	


names = ["AU cWW", 
"SR triple",
"SR core",
"sarcin3geometric",
"cWW in NMR",
"cWWpairWithAminoAcids", 
"RNA-protein-cWW-minor-groove", 
"LR cWW", 
"GU 2 into junction", 
"bSS and cWW",
"stacked cWW",
"AU cWW",
"SR triple",
"SR core",
"stacked cWW",
"sarcin3geometric",
"cWWpairWithAminoAcids",
"GU 2 into junction",
"bSS and cWW",
"sarcin5mixed2",
"GU Case 1 UG cWW",
"continuity test",
"unary test",
"UGGU tandem",
"GUUG tandem",
"GGUC_GGUC tandem",
"GU tandem at end of helix",
"GU tandem at end of helix adjacent above",
"GU tandem at end of helix adjacent below",
"GU tandem at end of helix one unpaired",
"GU tandem at end of helix two unpaired",
"UG tandem at end of helix",
"UG tandem at end of helix adjacent above",
"UG tandem at end of helix adjacent below",
"UG tandem at end of helix one unpaired",
"UG tandem at end of helix two unpaired",
"AG tHS",
"tSS LR",
"bSS distance 5,6",
"bPh",
"stacked bases with aa",
"base surrounded by bases",
"base surrounded by amino acids",
"base almost surrounded by amino acids",
"junction 10",
"internal loop",
"Hairpin flanking pair",
"Hairpin interacts with single strand",
"borderSS",
"SR triple",
"cWW with amino acid in minor groove",
"UU Case 0",
"UU Case 2 UU interior",
"UU Case 2 UU cWW",
"GC Case 1b CG cWW",
"GC Case 1c CG cWW",
"GC Case 1 CG interior",
"GC Case 2 GC cWW",
"GC Case 2b GC cWW",
"GC Case 2c GC cWW",
"GC Case 2 GC interior",
"GC Case 0",
"GC Case 1 CG cWW",
"AU Case 2 AU cWW",
"AU Case 2b AU cWW",
"AU Case 2c AU cWW",
"AU Case 2 AU interior",
"AU Case 0",
"AU Case 1 UA cWW",
"AU Case 1b UA cWW",
"AU Case 1c UA cWW",
"AU Case 1 UA interior",
"GU Case 0",
"GU Case 1 UG cWW",
"GU Case 1b UG cWW",
"GU Case 1c UG cWW",
"GU Case 1 UG interior",
"GU Case 2 GU cWW",
"GU Case 2b GU cWW",
"GU Case 2c GU cWW",
"GU Case 2 GU interior",
"AU hairpin Case 1 UA cWW",
"AU hairpin Case 1 UA after cWW",
"GU Case 2b GU cWW",
"unary no pair",
"CC cWW almost",
'tHH BPh',
"kink turn 65553",
"sarcin5geometric",
"AG tHS",
"Decoding loop",
"Z step n+",
"cWW in NMR",
"cWW with amino acid in minor groove",
"RNA-protein-cWW-minor-groove",
"GNRA hairpin symbolic",
"sO4'5",
"sO4'3",
"Z step",
"IL_4Y4O_235",
"syn pair",
"syn stack",
"NUNNGN tSW",
"next",
"chi_angle",
"RNA-protein4",
'modified',
'cWW with modified',
"sO",
"Python modified bases",
"LR cWW",
"chain length",
"count BP",
"count NT",
'Compare Matlab Python',
"AA cWW",
"DNA",
"sarcin5geometric",
"sarcin5mixed",
"sarcin13mixed",
"GoU 6S0Z",
"GoG 5NJT",
"cWW with modified",
"Stacked cWW GC",
'cWW in DNA',
'cWW stack in DNA'
]


def main():
	for i in range(len(names)):
		print("\n########### %s ###########" % names[i])
		test(search = names[i])

if __name__ == "__main__":
	main()
	#breakpoint() # used for debugging

# eventually:
# for name in names:
# 	output = test(name)
# 	#find a way to CAPTURE output and compare it to known outputs
# 	#when comparisons disagree, print out which search(es) now get different output
# 	if knowns[name] != output:
# 		print("Issue on search {name}\nShould return: {knowns[name]}\nReturned:{output}\n\n")


# "sarcin5mixed2"
# "GU Case 1 UG cWW"
# "continuity test"
# "unary test"

# "UGGU tandem"
# "GUUG tandem"
# "GGUC_GGUC tandem"
# "GU tandem at end of helix"
# "GU tandem at end of helix adjacent above"
# "GU tandem at end of helix adjacent below"
# "GU tandem at end of helix one unpaired"
# "GU tandem at end of helix two unpaired"

# "UG tandem at end of helix"
# "UG tandem at end of helix adjacent above"
# "UG tandem at end of helix adjacent below"
# "UG tandem at end of helix one unpaired"
# "UG tandem at end of helix two unpaired"

# ""
# "AG tHS"
# "tSS LR"
# "bSS distance 5,6"
# "bPh"

# "RNA-protein4"
# "stacked bases with aa"
# "base surrounded by bases"
# "base surrounded by amino acids"
# "base almost surrounded by amino acids"
# "junction 10"
# "internal loop"
# "Hairpin flanking pair"
# "Hairpin interacts with single strand"
# "borderSS"
# "SR triple"
# "cWW with amino acid in minor groove"

# "UU Case 0"
# "UU Case 2 UU interior"
# "UU Case 2 UU cWW"
# ["UU Case 0","UU Case 2 UU interior","UU Case 2 UU cWW"]

# "GC Case 1b CG cWW"
# "GC Case 1c CG cWW"
# "GC Case 1 CG interior"
# "GC Case 2 GC cWW"
# "GC Case 2b GC cWW"
# "GC Case 2c GC cWW"
# "GC Case 2 GC interior"
# "GC Case 0"
# "GC Case 1 CG cWW"
# ["GC Case 0","GC Case 1 CG cWW","GC Case 1b CG cWW","GC Case 1c CG cWW","GC Case 1 CG interior","GC Case 2 GC cWW","GC Case 2b GC cWW","GC Case 2c GC cWW","GC Case 2 GC interior"]

# "AU Case 2 AU cWW"
# "AU Case 2b AU cWW"
# "AU Case 2c AU cWW"
# "AU Case 2 AU interior"
# "AU Case 0"
# "AU Case 1 UA cWW"
# "AU Case 1b UA cWW"
# "AU Case 1c UA cWW"
# "AU Case 1 UA interior"
# ["AU Case 0","AU Case 1 UA cWW","AU Case 1b UA cWW","AU Case 1c UA cWW","AU Case 1 UA interior","AU Case 2 AU cWW","AU Case 2b AU cWW","AU Case 2c AU cWW","AU Case 2 AU interior"]

# "GU Case 0"
# "GU Case 1 UG cWW"
# "GU Case 1b UG cWW"
# "GU Case 1c UG cWW"
# "GU Case 1 UG interior"
# "GU Case 2 GU cWW"
# "GU Case 2b GU cWW"
# "GU Case 2c GU cWW"
# "GU Case 2 GU interior"
# ["GU Case 0","GU Case 1 UG cWW","GU Case 1b UG cWW","GU Case 1c UG cWW","GU Case 1 UG interior","GU Case 2 GU cWW","GU Case 2b GU cWW","GU Case 2c GU cWW","GU Case 2 GU interior"]

# "AU hairpin Case 1 UA cWW"
# "AU hairpin Case 1 UA after cWW"
# "cWW with amino acid in minor groove"
# "GU Case 2b GU cWW"

# "unary no pair"

# "CC cWW almost"
# "syn pair"
# 'tHH BPh'
# "kink turn 65553"
# "syn stack"

# "Compare Matlab Python"
# "sarcin5mixed"
# "sarcin5geometric"
# "GNRA hairpin symbolic"
# "Z step"
# "AG tHS"