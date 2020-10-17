# this program writes the list of corresponding atoms for modified and parent nucleotides as Matlab code
# that can be used with the Matlab version of FR3D

from modified_parent_mapping import modified_nucleotides

t = ""
t += "modnt = containers.Map;\n"


for key in modified_nucleotides.keys():
	t += "ParentAtomFromModAtom = containers.Map;\n"
	t += "ModAtomFromParentAtom = containers.Map;\n"
	for atom1, atom2 in modified_nucleotides[key]['atoms'].items():
		t += "ParentAtomFromModAtom('" + atom1 + "') = '" + atom2 + "';\n"
		t += "ModAtomFromParentAtom('" + atom2 + "') = '" + atom1 + "';\n"
	t += "newdict = containers.Map;\n"
	t += "newdict('standard') = '" + modified_nucleotides[key]["standard"] + "';\n"
	t += "newdict('parentfrommod') = ParentAtomFromModAtom;\n"
	t += "newdict('modfromparent') = ModAtomFromParentAtom;\n"
	t += "modnt('" + key + "') = newdict;\n"

print(t)