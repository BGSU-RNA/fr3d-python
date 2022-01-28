# generate_location_checks.py generates Python code to check the location of a point
# relative to standard locations of atoms

from fr3d.definitions import NAbasecoordinates


def generate_code(p1,p2,a1,a2):
	"""
	Generate a line of python code to tell when a generic point (x,y)
	lies in the plane to the left of the line segment from p1 to p2.
	Points p1 and p2 are two-dimensional.
	"""

	cx = p1[1] - p2[1]
	cy = p2[0] - p1[0]
	c  = p1[0]*p2[1] - p1[1]*p2[0]

	t = "if %9.6f*x + %9.6f*y + %9.6f > 0:  # Left of %s-%s" % (cx,cy,c,a1,a2)

	return t


atom_lists = {}
atom_lists['A'] = ['C4','C5','N7','C8','N9','C4','N3','C2','N1','C6','C5']
atom_lists['C'] = ['N1','C2','N3','C4','C5','C6','N1']
atom_lists['G'] = ['C4','C5','N7','C8','N9','C4','N3','C2','N1','C6','C5']
atom_lists['U'] = ['N1','C2','N3','C4','C5','C6','N1']
atom_lists['DT'] = ['N1','C2','N3','C4','C5','C6','N1']

for sequence in atom_lists.keys():
	print(sequence)
	atom_list = atom_lists[sequence]
	for i in range(1,len(atom_list)):
		a1 = atom_list[i-1]
		a2 = atom_list[i]
		p1 = NAbasecoordinates[sequence][a1][0:2]
		p2 = NAbasecoordinates[sequence][a2][0:2]
		print(generate_code(p1,p2,a1,a2))

