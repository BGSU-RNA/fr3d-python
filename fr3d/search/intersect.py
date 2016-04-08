"""
intersect.py
by Dan Schellhas
finalized February 2014

This collection of functions takes a list of lists of contraints and returns a list of
candidates that meet all of the given constraints. Contraints are given in the
form constraints[position1][position2] = [[object1, object2], [object3, object4], ...]
For example, constraints[0][1] = [[1,2],[7,11]] would mean that the constraint 
between positions 0 and 1 is only satisfied by objects 1 and 2 along with 7 and 11.

The code used to load the constraints in testing was as follows and had the position 
numbers as the first two numbers of each line in a CSV:
	def loadConstraints(filename):
		fp = open(filename)
		constraints = defaultdict(dict)
		for line in fp:
			c = line.strip().split(",")
			i = int(c[0])
			j = int(c[1])
			pairs = []
			first = None
			for x in c[2:]:
				if first:
					pairs.append((first, int(x)))
					first = None
				else:
					first = int(x)
			constraints[i][j] = pairs
		return constraints

The algorithm reorders the contraints prior to processing into a good (but not proven
to be optimal) order. It then loads the constraints into a tree for easy searching. 
Using this tree of lists, the algorithm then generates candidates through set 
intersections and returns them to the user. The runtime of candidate selection is also 
returned for help in determining reasonable query sizes.

Main function:
candidates[], runtime = intersect(constraints[][][])

Support functions:
newConstraints[][][] = getPermutation(constraints[][][], perm[])
cmap[][][], atoms = buildMap(constraints[][][])
candidates[], runtime = getCandidates(cmap[][][], atoms)
candidates[] = buildCandidate(soFar[], options[], cmap[][][], atoms)

"""
from time import time
from copy import copy
from collections import defaultdict

def intersect(constraints):
    # This section does the reordering for efficiency
    atoms = len(constraints) + 1
    constraintSizes = []
    for i in constraints:
        for j in constraints[i]:
            L = len(constraints[i][j])
            constraintSizes.append(([i,j], L))
    constraintSizes = sorted(constraintSizes, key=lambda x: x[1])
    perm = constraintSizes[0][0]
    del constraintSizes
    objects = set(range(1, atoms+1))
    objects.remove(perm[0])
    objects.remove(perm[1])
    while len(perm) < atoms:
        bestLen = float("inf")
        bestCon = -1
        for obj in objects:
            myLen = 1
            for p in perm:
                if p < obj:
                    myLen *= len(constraints[p][obj])
                else:
                    myLen *= len(constraints[obj][p])
            if myLen < bestLen:
                bestLen = myLen
                bestCon = obj
        perm.append(bestCon)
        objects.remove(bestCon)
    # This puts the constraints into that better order
    newConstraints = getPermutation(constraints, perm)
    del constraints
    # This builds the easily searchable tree
    cmap, atoms = buildMap(newConstraints)
    del newConstraints
    # This uses that tree to generate the candidates
    candidates, runtime = getCandidates(cmap, atoms)
    return candidates, runtime


def getPermutation(constraints, perm):
    newConstraints = {}
    indexOf = {}
    for i in range(0, len(perm)):
        indexOf[perm[i]] = i + 1
    for i in constraints:
        for j in constraints[i]:
            if indexOf[i] > indexOf[j]:
                pairs = []
                for pair in constraints[i][j]:
                    pair = (pair[1], pair[0])
                    pairs.append(pair)
                if indexOf[j] in newConstraints:
                    newConstraints[indexOf[j]][indexOf[i]] = pairs
                else:
                    newConstraints[indexOf[j]] = {indexOf[i]:pairs}
            else:
                if indexOf[i] in newConstraints:
                    newConstraints[indexOf[i]][indexOf[j]] = copy(constraints[i][j])
                else:
                    newConstraints[indexOf[i]] = {indexOf[j]:copy(constraints[i][j])}
    return newConstraints
					
def buildMap(constraints):
    atoms = max(k for k, v in constraints[1].iteritems() if v != 0)
    cmap = [0] * (atoms + 1)
    for i in range(1,atoms):
        cmap[i] = [0] * (atoms + 1)
        for j in range(i + 1, atoms + 1):
            cmap[i][j] = defaultdict(list)
			
    for i in range(1, atoms):
        for j in range(i + 1, atoms + 1):
            for c in constraints[i][j]:
                cmap[i][j][c[0]].append(c[1])
    return cmap, atoms
	
def getCandidates(cmap, atoms):
    startTime = time()
    candidates = []
    for start in cmap[1][2]:
        useIt = True
        for i in range(2, atoms + 1):
            if start not in cmap[1][i]:
                useIt = False
                break
        if useIt:
            candidates += buildCandidate([start], cmap[1][2][start], cmap, atoms)
    return candidates, time() - startTime

def buildCandidate(soFar, options, cmap, atoms):
    candidates = []
    if len(soFar) + 1 == atoms:
        for opt in options:
            candidates.append(soFar+[opt])
        return candidates
	
    if soFar[0] in cmap[1][len(soFar) + 2]:
        baseOptions = set(cmap[1][len(soFar) + 2][soFar[0]])
        for i in range(2, len(soFar) + 1):
            if len(baseOptions):
                if soFar[i - 1] in cmap[i][len(soFar) + 2]:
                    baseOptions = baseOptions.intersection(set(cmap[i][len(soFar) + 2][soFar[i - 1]]))
                else:
                    return []
            else:
                return []
        if len(baseOptions):
            for opt in options:
                newOptions = baseOptions
                if opt in cmap[len(soFar) + 1][len(soFar) + 2]:
                    newOptions = newOptions.intersection(set(cmap[len(soFar) + 1][len(soFar) + 2][opt]))
                else:
                    continue
                if len(newOptions): candidates += buildCandidate(soFar + [opt], newOptions, cmap, atoms)
			
    return candidates	

