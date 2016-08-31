import math
from random import random

# Metric classes for use with similarity ordering or traveling salesman problems
#
# Any metric class should define a set "points" and a method "d" that
# takes two points and returns the distance between them. Regardless of the
# file supplied, the user will get a set of points to manipulate as they
# see fit an a function that gives them the distances they need for their
# algorithm in a black box fashion.
#
# Originally designed by Dan Schellhas
# Last updated Feb 25, 2015

def Metric(filename, mode="unif2D"):
    if filename.isdigit():
      if mode == "unif2D":
        return MetricUniform2D(int(filename))
    else:
      with open(filename) as currentFile:
        curMode = "HEADER"
        curPoint = 0
        for line in currentFile:
          lineData = line.strip('\n').split(' ')
          if curMode == "HEADER":
            if lineData[0][:9] == "DIMENSION":
              size = int(lineData[len(lineData)-1].strip())
            if lineData[0][:16] == "EDGE_WEIGHT_TYPE":
              if lineData[len(lineData)-1].strip() == "EUC_2D":
                # Return the TSP-LIB 2D Euclidean metric
                if size > 1000:
                  return MetricTSPeuc2DnoMemo(filename)
                else:
                  return MetricTSPeuc2D(filename)
              else:
                # No other metrics currently defined
                print(lineData[len(lineData)-1])
                return False


# TSP-LIB 2D Euclidean metric with memoization
# A second version could be made for datasets too large for memoization
class MetricTSPeuc2D:

  points = []
  data = []
  memo = {}

  def __init__ (self, filename):
    self.points = []
    self.data = []
    self.memo = {}
    with open(filename) as currentFile:
      curMode = "HEADER"
      curPoint = 0
      for line in currentFile:
        lineData = line.strip('\n').strip('\r').split(' ')
        if curMode == "HEADER":
          if lineData[0][:16] == "EDGE_WEIGHT_TYPE":
            if lineData[len(lineData)-1].strip() != "EUC_2D":
              # Error!
              print "Invalid input file! Attempting anyway."
          if lineData[0][:18] == "NODE_COORD_SECTION": curMode = "DATA"
        elif curMode == "DATA" and lineData[0].isdigit():
          self.data.append(lineData[1:])
          self.points.append(curPoint)
          curPoint += 1

  def d(self, point1, point2):
    if not self.memo.has_key(point1):
      self.memo[point1] = {point2: math.sqrt(math.pow(float(self.data[point1][0]) - float(self.data[point2][0]),2) + math.pow(float(self.data[point1][1]) - float(self.data[point2][1]),2))}
    elif not self.memo[point1].has_key(point2):
      self.memo[point1][point2] = math.sqrt(math.pow(float(self.data[point1][0]) - float(self.data[point2][0]),2) + math.pow(float(self.data[point1][1]) - float(self.data[point2][1]),2))
    return self.memo[point1][point2]

# TSP-LIB 2D Euclidean metric without memoization
# A second version could be made for datasets too large for memoization
class MetricTSPeuc2DnoMemo:

  points = []
  data = []

  def __init__ (self, filename):
    self.points = []
    self.data = []
    with open(filename) as currentFile:
      curMode = "HEADER"
      curPoint = 0
      for line in currentFile:
        lineData = line.strip('\n').strip('\r').split(' ')
        if curMode == "HEADER":
          if lineData[0][:16] == "EDGE_WEIGHT_TYPE":
            if lineData[len(lineData)-1].strip() != "EUC_2D":
              # Error!
              print "Invalid input file! Attempting anyway."
          if lineData[0][:18] == "NODE_COORD_SECTION": curMode = "DATA"
        elif curMode == "DATA" and lineData[0].isdigit():
          self.data.append(lineData[1:])
          self.points.append(curPoint)
          curPoint += 1

  def d(self, point1, point2):
    return math.sqrt(math.pow(float(self.data[point1][0]) - float(self.data[point2][0]),2) + math.pow(float(self.data[point1][1]) - float(self.data[point2][1]),2))

# Two dimensional uniform distribution of points
#
class MetricUniform2D:

  points = []
  data = []
  memo = {}

  def __init__ (self, n, ratio=1):
    self.points = []
    self.data = []
    self.memo = {}
    curPoint = 0
    for i in range(0, n):
      self.data.append([random(),ratio*random()])
      self.points.append(curPoint)
      curPoint += 1

  def d(self, point1, point2):
    if not self.memo.has_key(point1):
      self.memo[point1] = {point2: math.sqrt(math.pow(float(self.data[point1][0]) - float(self.data[point2][0]),2) + math.pow(float(self.data[point1][1]) - float(self.data[point2][1]),2))}
    elif not self.memo[point1].has_key(point2):
      self.memo[point1][point2] = math.sqrt(math.pow(float(self.data[point1][0]) - float(self.data[point2][0]),2) + math.pow(float(self.data[point1][1]) - float(self.data[point2][1]),2))
    return self.memo[point1][point2]



# Set up a Metric object given a matrix of distances
#
class GivenDistances:

  points = []
  data = []
  memo = {}

  def __init__ (self, distances):
    self.points = []
    self.data = []
    self.memo = {}
    numPoints = len(distances[0])
    self.distances = distances
    curPoint = 0
    for i in range(0, numPoints):
      self.points.append(curPoint)
      curPoint += 1

  def d(self, point1, point2):
    return self.distances[point1,point2]
