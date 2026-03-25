# -*- coding: utf-8 -*-
"""
Created on Sun Nov 10 3:00:00 2021

@author:adamc
"""

import sys
import os
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
		return(subprocess.run('python FR3D.py "%s"' % search))


def main():
	# read files in queries directory
	import glob
	names = glob.glob("queries/*.json")

	# put names in random order
	import random
	random.shuffle(names)

	for i,name in enumerate(names):
		if not "entire" in name:
			print("\n########### %s %d of %d ###########" % (name,i,len(names)))
			test(name)

if __name__ == "__main__":
	main()
