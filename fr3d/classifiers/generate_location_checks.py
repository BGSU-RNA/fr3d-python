# generate_location_checks.py generates Python code to check the location of a point
# relative to standard locations of atoms.
# Parameter r tells how many Angstroms it can be to the right of the line.

from fr3d.definitions import NAbasecoordinates
import math
import numpy as np

def generate_code(p1,p2,a1,a2,i,r=0):
	"""
	Generate a line of python code to tell when a generic point (x,y)
	lies in the plane to the left of the line segment from p1 to p2.
	Points p1 and p2 are two-dimensional.
	a1 and a2 are atom labels.
	i tells how many tabs to put at the start of the line
	r allows for the point (x,y) to be within r units to the right
	of the line, for near interactions.
	"""

	cx = p1[1] - p2[1]
	cy = p2[0] - p1[0]
	c  = p1[0]*p2[1] - p1[1]*p2[0] + r*math.sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

	t = '\t'*i

	if r == 0:
		t += "if %9.6f*x + %9.6f*y + %9.6f > 0:  # Left of %s-%s" % (cx,cy,c,a1,a2)
	else:
		t += "if %9.6f*x + %9.6f*y + %9.6f > 0:  # Within %9.6f Angstroms of being left of %s-%s" % (cx,cy,c,r,a1,a2)

	return t

def iterate_over_atom_lists(atom_lists,r=0):
	for sequence in sorted(atom_lists.keys()):
		print(sequence)
		atom_list = atom_lists[sequence]
		for i in range(1,len(atom_list)):
			a1 = atom_list[i-1]
			a2 = atom_list[i]
			p1 = NAbasecoordinates[sequence][a1][0:2]
			p2 = NAbasecoordinates[sequence][a2][0:2]
			print(generate_code(p1,p2,a1,a2,i,r))

def iterate_over_rings(atom_lists,r=0):
	for sequence in sorted(atom_lists.keys()):
		print(sequence)
		atom_list = atom_lists[sequence]
		seq = sequence.replace("5","").replace("6","")  # actual sequence
		s = np.array([0.0,0.0])
		for atom in atom_list:
			s += NAbasecoordinates[seq][atom][0:2]
		s = s / len(atom_list)         # center of ring
		d2 = {}
		for atom in atom_list:
			p = NAbasecoordinates[seq][atom][0]
			q = NAbasecoordinates[seq][atom][1]
			d2[atom] = math.sqrt((p-s[0])**2 + (q-s[1])**2)  # distance squared
			#print(atom, d2[atom])
		m = max([d2[atom] for atom in d2])

		# find coefficients of ellipse a(x-s[0])^2 + b(x-s[0])(y-s[0]) + (y-s[1])^2 - k = 0
		# that best approximates the points
		amin = 1.0
		bmin = 0.0
		kmin = m**2

		minscore = 100
		scoregap = 100
		for i in range(0,100000):
			a = amin + np.random.normal(0,0.01,1)
			b = bmin + np.random.normal(0,0.01,1)
			k = kmin + np.random.normal(0,0.01,1)
			score = 0
			for atom in atom_list:
				p = NAbasecoordinates[seq][atom][0]
				q = NAbasecoordinates[seq][atom][1]
				score += abs((a*(p-s[0])**2 + b*(p-s[0])*(q-s[1]) + (q-s[1])**2 - k))
			if score < minscore:
				amin = a
				bmin = b
				kmin = k
				scoregap = minscore - score
				minscore = score
				if i > 5000:
					print('i=%5d a=%0.6f b=%0.6f k=%0.6f score=%0.6f' % (i,a,b,k,score))

		a = amin
		b = bmin
		k = (math.sqrt(kmin)+r)**2  # push out by distance r
		print('if %0.6f*(x-(%0.6f))**2 + %0.6f*(x-(%0.6f))*(y-(%0.6f)) + (y-(%0.6f))**2 < %0.6f:  # %s r=%0.1f' % (a,s[0],b,s[0],s[1],s[1],k,sequence,r))

		print('if nt1_seq == "%s":' % seq)
		print('    for i in range(0,1000):       # %s' % sequence)
		print('        t = 6.28318530718*i/1000')
		print('        r = math.sqrt(%0.6f/(%0.6f*math.cos(t)**2+%0.6f*math.cos(t)*math.sin(t)+math.sin(t)**2))' % (k,a,b))
		print('        x.append(%0.6f + r*math.cos(t))' % s[0])
		print('        y.append(%0.6f + r*math.sin(t))' % s[1])

ring_atom_lists = {}
ring_atom_lists['A'] = ['C4','C5','N7','C8','N9','C4','N3','C2','N1','C6','C5']
ring_atom_lists['C'] = ['N1','C2','N3','C4','C5','C6','N1']
ring_atom_lists['G'] = ['C4','C5','N7','C8','N9','C4','N3','C2','N1','C6','C5']
ring_atom_lists['DT'] = ['N1','C2','N3','C4','C5','C6','N1']
ring_atom_lists['U'] = ['N1','C2','N3','C4','C5','C6','N1']

ring_lists = {}
ring_lists['A6'] = ['C4','N3','C2','N1','C6','C5']
ring_lists['A5'] = ['C4','C5','N7','C8','N9']
ring_lists['C'] = ['N1','C2','N3','C4','C5','C6']
ring_lists['G6'] = ['C4','N3','C2','N1','C6','C5']
ring_lists['G5'] = ['C4','C5','N7','C8','N9']
ring_lists['DT'] = ['N1','C2','N3','C4','C5','C6']
ring_lists['U'] = ['N1','C2','N3','C4','C5','C6']

convexHullAtoms = {}
convexHullAtoms['A'] = ["C1'",'N3','H2','H61','H62','H8',"C1'"] #Based on Matlab Code
convexHullAtoms['DA'] = ["C1'",'H2','N6','H8',"C1'"]
convexHullAtoms['C'] = ["C1'",'O2','N3','N4','H5','H6', "C1'"]
convexHullAtoms['DC'] = ["C1'",'O2','N3','N4','H5','H6', "C1'"]
convexHullAtoms['G'] = ["C1'",'H21','H22','H1','O6','N7','H8',"C1'"]
convexHullAtoms['DG'] = ["C1'",'H21','H22','H1','O6','N7','H8',"C1'"]
convexHullAtoms['U'] = ["C1'",'O2','H3','O4','H5','H6',"C1'"]
convexHullAtoms['DT'] = ["C1'",'O2','H3','O4','C7', 'C6', "C1'"] 

print('Code to check that an (x,y) point is inside a base ring')
iterate_over_atom_lists(ring_atom_lists,0)

print('Code to check that an (x,y) point is inside the convex hull of a base')
iterate_over_atom_lists(convexHullAtoms,0)

print('Code to check that an (x,y) point is close to being inside a base ring')
iterate_over_atom_lists(ring_atom_lists,0.4)

print('Code to check that an (x,y) point is close to the center of a ring')
iterate_over_rings(ring_lists,0.3)