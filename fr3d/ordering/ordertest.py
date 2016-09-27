#!/usr/local/bin/python -O

import matplotlib.pyplot as plt
import numpy as np
import random
from random import shuffle
from metric import Metric
from metric import MetricUniform2D
from metric import GivenDistances
from simHeat import simHeat
from greedyInsertion import greedyInsertion
from greedyInsertion import pathLength
from greedyInsertion import orderWithPathLengthFromDistanceMatrix

n = 20                        # number of data points
print("Generating "+str(n)+" data points")

# two-dimensional data points
m = MetricUniform2D(n,30)     # 1 by 30 rectangle in 2d
#print(m.points)               # indices 0 to n-1
#print(m.data)                 # list of 2d data points
#print(m.memo)

bestPathLength = float("inf")
numReps = 100
print("Running greedy insertion "+str(numReps)+" times")
for i in range(0,numReps):
  order = greedyInsertion(m,depth=1)
  newPathLength = pathLength(m,order[0])
  if newPathLength < bestPathLength:
    bestPathLength = newPathLength
    bestOrder = order[0]
    print("Improved path length "+str(bestPathLength)+" on replication " + str(i))
#    print(order[0])

distances = np.zeros((len(bestOrder),len(bestOrder)))
displaydistances = np.zeros((len(bestOrder),len(bestOrder)))

for i in range(0,n):
  for j in range(0,n):
    distances[i,j] = m.d(i,j)
    displaydistances[i,j] = m.d(i,j)

# replace some real distances with None to simulate a missing data problem

if 10 < 1:
  for i in range(0,n):
    for j in range(0,n):
      if random.random() < 0.05:
        distances[i,j] = None
        distances[j,i] = None
        displaydistances[i,j] = -1
        displaydistances[j,i] = -1

# replace some rows and columns with None to simulate a different missing data problem
if 0 < 1:
  for i in range(0,n):
    if random.random() < 0.1:
      for j in range(0,n):
        distances[i,j] = None
        distances[j,i] = None
        displaydistances[i,j] = -1
        displaydistances[j,i] = -1

# replace first row and column with None to simulate a different missing data problem
if 10 < 1:
  i = 1
  for j in range(0,n):
    distances[i,j] = None
    distances[j,i] = None
    displaydistances[i,j] = -1
    displaydistances[j,i] = -1

print("Repackaging data into a new data structure, ordering with 1 replication")
newData = GivenDistances(distances)
neworder = greedyInsertion(newData,depth=1)
print("Path length " + str(pathLength(newData,neworder[0])))

print("Using orderWithPathLengthFromDistanceMatrix, 1 replication")
order,score,d = orderWithPathLengthFromDistanceMatrix(distances,True)
print("Path length " + str(score))
print("Using orderWithPathLengthFromDistanceMatrix, 10 replications")
order,score,d = orderWithPathLengthFromDistanceMatrix(distances,100,True)
print("Path length " + str(score))

ordereddistances = np.zeros((len(bestOrder),len(bestOrder)))

for i in range(0,n):
  for j in range(0,n):
#    ordereddistances[i,j] = displaydistances[order[i],order[j]]
    ordereddistances[i,j] = d[order[i],order[j]]

print("Displaying heat map of distance matrix from first ordering made")
column_labels = list('ABCD')
row_labels = list('WXYZ')
data = np.random.rand(4,4)

fig, ax = plt.subplots()
heatmap = ax.pcolor(ordereddistances, cmap=plt.cm.bwr)

# put the major ticks at the middle of each cell
#ax.set_xticks(np.arange(data.shape[0])+0.5, minor=False)
#ax.set_yticks(np.arange(data.shape[1])+0.5, minor=False)

# want a more natural, table-like display
ax.invert_yaxis()
ax.xaxis.tick_top()

#ax.set_xticklabels(row_labels, minor=False)
#ax.set_yticklabels(column_labels, minor=False)
plt.show()
