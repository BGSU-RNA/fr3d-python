#!/usr/local/bin/python -O

from random import shuffle
from fr3d.ordering.metric import GivenDistances
from numpy import median
from numpy import isnan

'''
	This method inserts each point at its optimal location in a given order.
	The first two points are free since there is no context.

	m - The data of the problem in a Metric object
	w - The weight vector
	o - The order to insert points
	verbose - Print the results

	greedyInsertion returns the path that it finds as a list of points from the Metric

    Originally developed by Dan Schellhas
    Modified on 8-30-2016 by Craig Zirbel, adding pathLength
    Modified on 8-30-2016 by Craig Zirbel, adding orderWithPathLengthFromDistanceMatrix
'''

def pathLength(m, path):
  score = 0
  for i in range(0, len(path)-1):
    score += m.d(path[i], path[i+1])
  return score

def testScore(m, path, depth = False):
	if depth == False:
		depth = len(path)
	testScore = 0
	for i in range(0, len(path)):
		for j in range(i+1, min(i+depth+1, len(path))):
			testScore += m.d(path[i], path[j])/abs(j-i)
	return testScore

def orderWithPathLengthFromDistanceMatrix(distances,numReps = 1,scanForNan = False):

  if scanForNan:
    values = []
    for i in range(0,len(distances)):
      for j in range(0,len(distances[0])):
        if not isnan(distances[i,j]):
          values.append(distances[i,j])

    m = 0
    if len(values) > 0:
      m = median(values)
      for i in range(0,len(distances)):
        if isnan(distances[i,i]):
          distances[i,i] = 0
        for j in range(0,len(distances[0])):
          if isnan(distances[i,j]):
            distances[i,j] = m

  dataset = GivenDistances(distances)

  bestOrder = range(0,len(distances))

  if len(distances) < 3:
    bestScore = pathLength(dataset,bestOrder)
    return bestOrder, bestScore, distances

  bestPathLength = float("inf")

  for i in range(0,numReps):
    results = greedyInsertionPathLength(dataset)
    newPathLength = pathLength(dataset,results[0])
    if newPathLength < bestPathLength:
      bestPathLength = newPathLength
      bestOrder = results[0]

  return bestOrder, bestPathLength, distances

def greedyInsertionPathLength(m, o=[], verbose=False):

  # if no starting ordering
  if len(o) == 0:
  	o = m.points
  	shuffle(o)          # random starting ordering
  path = o[:2]          # first two points of the current ordering
  score = m.d(path[0], path[1])

  for p in range(2, len(o)):
    # score inserting point o[p] at beginning of path
    bestScore = m.d(o[p],path[0])
    bestPosition = 0

    # score inserting point o[p] at end of path
    currentScore = m.d(path[-1],o[p])
    if currentScore < bestScore:
      bestScore = m.d(path[-1],o[p])
      bestPosition = len(path)

    # score inserting point o[p] at various points within the path
    for position in range(1, len(path)):
      currentScore = m.d(path[position-1],o[p]) + m.d(o[p],path[position]) - m.d(path[position-1],path[position])

      if currentScore < bestScore:
        bestScore = currentScore
        bestPosition = position

    path.insert(bestPosition, o[p])
    score += bestScore

  return path, score

def greedyInsertionHarmonic(m, w=False, o=[], depth=False, verbose=False):
  if w == False:
    w = []
    for i in range(1,len(m.points)+1):
      if i <= depth:
        w.append(1/(i*1.0))
      else:
        w.append(0)
  if len(o) == 0:
    o = m.points
    shuffle(o)
  if depth == False:
    depth = len(o)
  path = o[:2]
  score = w[0] * m.d(path[0], path[1])

  for p in range(2, len(o)):
    distances = []
    for point in path:
      distances.append(m.d(o[p], point))
    for position in range(0, len(path)+1):
      myScore = 0
      l = min(len(path), depth)
      s = max(0, position - l)
      e = min(len(path), position + l)
      for i in range(s, e):
        d = int(abs(i+.5-position))
        myScore += distances[i] * w[d]
      if position == 0:
        tempScore = score + myScore
        bestPosition = 0
        bestScore = tempScore
        hmm = [o[p]] + path
      else:
        tempScore += myScore - myOldScore
        for i in range(max(0,s-1), position-1):
          tempScore += m.d(path[i], path[position - 1]) * (w[abs(position - 1 - i) - 1] - w[abs(position - 1 - i)])
        for i in range(position, e):
          tempScore += m.d(path[i], path[position - 1]) * (w[abs(position - 1 - i)] - w[abs(position - 1 - i) - 1])

      if tempScore < bestScore:
        bestScore = tempScore
        bestPosition = position

      myOldScore = myScore
    path.insert(bestPosition, o[p])
    score = bestScore
  return path, score
