"""
Implement tree-penalized path length ordering algorithm described in the paper
https://www.sciencedirect.com/science/article/pii/S037722172200501X by Aliyev and Zirbel
"""

import math
import numpy as np
import random

def treePenalty(distance):
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import squareform

    # set up groups to track hierarchical clustering
    group = []
    for i in range(0,distance.shape[0]):
        group.append([i])

    # hierarchical clustering with average linkage
    Z = linkage(squareform(distance), "average")

    # penalty matrix
    penalty = np.zeros(distance.shape)

    for merger in Z:
        a = int(merger[0])
        b = int(merger[1])
        group.append(group[a]+group[b])
        for i in group[a]:
            for j in group[b]:
                penalty[i][j] = merger[2]
                penalty[j][i] = merger[2]

    """
    for i in range(0,distance.shape[0]):
        for j in range(0,distance.shape[0]):
            print("%6.3f" % penalty[i][j]),
        print("")
    """

    return penalty

def normalizePenaltyMatrix(distance,penalty):
    # normalize the penalty matrix as described in the tpPL article
    dsum = 0.0
    psum = 0.0

    newpenalty = np.zeros(distance.shape)

    for i in range(0,distance.shape[0]):
        for j in range(0,distance.shape[0]):
            dsum = dsum + distance[i][j]
            psum = psum + penalty[i][j]

    if psum > 0:
        for i in range(0,distance.shape[0]):
            for j in range(0,distance.shape[0]):
                newpenalty[i][j] = penalty[i][j] * dsum / psum
    else:
        # in case of a distance matrix of all zeros
        newpenalty = penalty

    return newpenalty

def twoOptSwap(distance,order):
    # simple 2-opt swap where you check every possible swap and make the swap if it improves the score
    # do large-scale swaps first
    n = len(order)
    distance_reduction = 0
    for i in range(0,n):
        for p in range(n-i-1,1,-1):
            j = i + p
            # compute new path length minus old if you change the path
            distance_change = distance[order[i]][order[j-1]] + distance[order[i+1]][order[j]] - distance[order[i]][order[i+1]] - distance[order[j-1]][order[j]]
            if distance_change < -0.0000000000001:  # avoid changes due to roundoff error
                distance_reduction = distance_reduction + distance_change
                order = order[:i+1] + order[j-1:i:-1] + order[j:]  
    return order, distance_reduction

def greedyInsertionPathLength(distance, order=[]):
    # if no starting ordering, put points in random order
    if len(order) == 0:
        order = list(range(0,distance.shape[0]))
        random.shuffle(order)          # random starting ordering

    # make a path from the first two points of the ordering
    path = order[:2]
    score = distance[path[0], path[1]]

    # insert the remaining points into the path one by one
    for p in range(2, len(order)):
        # score inserting point order[p] at beginning of path
        bestScore = distance[order[p]][path[0]]
        bestPosition = 0

        # score inserting point order[p] at end of path
        currentScore = distance[path[-1]][order[p]]
        if currentScore < bestScore:
            bestScore = distance[path[-1]][order[p]]
            bestPosition = len(path)

        # score inserting point order[p] at various points within the path
        for position in range(1, len(path)):
            currentScore = distance[path[position-1]][order[p]] + distance[order[p]][path[position]] - distance[path[position-1]][path[position]]

            if currentScore < bestScore:
                bestScore = currentScore
                bestPosition = position

        # insert where the score is the lowest
        path.insert(bestPosition, order[p])
        score += bestScore

    return path, score

def multipleGreedyInsertionPathLength(distance, repetitions=100, seed=None):
    # repeat greedy insertion multiple times and keep the best ordering
    if seed:
        random.seed(seed)

    bestScore = float("inf")
    for rep in range(0,repetitions):
        path, score = greedyInsertionPathLength(distance)
        if score < bestScore:
            bestScore = score
            bestPath = path

    return bestPath

def multipleGreedyInsertionPathLengthTwoOpt(distance, repetitions=10, seed=None, bestScore=float("inf"), bestOrder=[]):
    # repeat greedy insertion followed by two-opt swaps multiple times and keep the best ordering

    n = distance.shape[0]

    if n == 0:
        bestOrder = []
        bestScore = 0
    elif n == 1:
        bestOrder = [0]
        bestScore = 0
    else:
        if seed:
            random.seed(seed)

        bestRep = None

        for rep in range(0,repetitions):
            order, score = greedyInsertionPathLength(distance)
            order, distance_reduction = twoOptSwap(distance,order)
            score = score + distance_reduction
            if score < bestScore:
                bestScore = score
                bestOrder = order

                # optional diagnostics
                if False:
                    print('Greedy insertion repetion %2d score %10.8f' % (rep+1,bestScore))


        bestOrder, distance_reduction = twoOptSwap(distance,bestOrder)
        bestScore = bestScore + distance_reduction

    return bestOrder, bestScore

def treePenalizedPathLength(distance,repetitions=10,seed=None,penaltyStrength=0.5):

    n = distance.shape[0]

    if n == 0:
        return []
    elif n == 1:
        return [0]
    elif n == 2:
        return [0,1]

    if penaltyStrength > 0:
        # tpPL penalized distance matrix
        penaltyMatrix = treePenalty(distance)
        penaltyMatrix = normalizePenaltyMatrix(distance,penaltyMatrix)
        penalizedDistance = distance + penaltyStrength * penaltyMatrix
        order, score = multipleGreedyInsertionPathLengthTwoOpt(penalizedDistance,repetitions,seed)
    else:
        # standard path length
        order, score = multipleGreedyInsertionPathLengthTwoOpt(distance,repetitions,seed)

    return order

def optimalLeafOrder(distance, method = "average"):
    # always produces the same ordering, but not always the best ordering
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import squareform

    Z = linkage(squareform(distance), method, optimal_ordering = True)
    dn = dendrogram(Z,no_plot=True)
    return dn['leaves']

def standardOrder(distance,order):
    # orient the path to have lower discrepancies in upper left of matrix
    n = distance.shape[0]
    if n <= 2:
        return order
    else:
        # get approximately the first third of the data
        if n == 3:
            m = 2
        elif n == 4:
            m = 3
        elif n == 5:
            m = 3
        elif n == 6:
            m = 4
        else:
            m = int(math.ceil(n/3))+1

        # sum the upper left and lower right triangles
        upperLeftSum = 0.0
        lowerRightSum = 0.0
        for i in range(0,m):
            for j in range(i+1,m):
                upperLeftSum  += distance[order[i]][order[j]]
                lowerRightSum += distance[order[n-i-1]][order[n-j-1]]

        if upperLeftSum > lowerRightSum:
            reversedorder = [order[n-i-1] for i in range(0,n)]
            return reversedorder
        else:
            return order

def setDiagonalToZero(distance):
    # linkage functions require a matrix whose diagonal is zero
    new_distance = np.copy(distance)
    for i in range(0,len(new_distance)):
        new_distance[i][i] = 0.0
    return new_distance

def reorderSymmetricMatrix(distance, newOrder):
    # apply the new ordering to a distance matrix, returning a copy
    newDistance = np.zeros(distance.shape)
    for i in range(0,newDistance.shape[0]):
        for j in range(i,newDistance.shape[1]):
            d = distance[newOrder[i]][newOrder[j]]
            newDistance[i][j] = d
            newDistance[j][i] = d
    return newDistance

def reorderList(oldList,newOrder):
    # apply the new ordering to a list, returning a copy
    newList = []
    for i in range(0,len(oldList)):
        newList.append(oldList[newOrder[i]])

    return newList

def imputeNANValues(distance):
    # replace any NAN or negative values with the maximum value in the matrix
    # This pushes the NAN values to the beginning or end of the ordering
    maxVal = 0
    for i in range(0,distance.shape[0]):
        for j in range(i+1,distance.shape[1]):
            if not math.isnan(distance[i][j]):
                maxVal = max(maxVal,distance[i][j])

    newDistance = np.zeros(distance.shape)
    for i in range(0,distance.shape[0]):
        for j in range(i+1,distance.shape[1]):
            if math.isnan(distance[i][j]) or distance[i][j] < 0:
                newDistance[i][j] = maxVal
                newDistance[j][i] = maxVal
            else:
                newDistance[i][j] = distance[i][j]
                newDistance[j][i] = distance[i][j]

        # OK to have positive values on the diagonal, but not NAN or negative
        if math.isnan(distance[i][i]) or distance[i][i] < 0:
            newDistance[i][i] = 0
        else:
            newDistance[i][i] = distance[i][i]

    return newDistance

def calculateDistanceMatrix(points):
    # calculate the distance between every pair of points
    n = len(points)
    distance = np.zeros((n,n))
    for i in range(0,n):
        for j in range(i+1,n):
            distance[i][j] = np.linalg.norm(points[i]-points[j])
            distance[j][i] = distance[i][j]
    return distance

def plotPoints(points,show=False):
    # plot the first two dimensions of points and connect with lines
    import matplotlib.pyplot as plt
    plt.plot(points[:,0],points[:,1],linewidth=1,c='cyan',zorder=1)
    plt.scatter(points[0,0],points[0,1],s=20,c='lightgreen',zorder=2)  # first point green dot
    plt.scatter(points[1:,0],points[1:,1],s=8,c='black',zorder=2)
    if show:
        plt.show()

def plotHeatMap(distance, labels=None, show=False):
    # plot the distance matrix as a heatmap
    import matplotlib.pyplot as plt
    plt.imshow(distance, cmap='viridis', interpolation='nearest')
    plt.colorbar()
    # put any labels next to the plot on x and y axes
    if labels:
        plt.xticks(range(0,len(labels)),labels,rotation=90)
        plt.yticks(range(0,len(labels)),labels)
    else:
        plt.xticks([])
        plt.yticks([])
    # if fewer than 30 points, put the actual values in the cells
    if distance.shape[0] < 30:
        for i in range(0,distance.shape[0]):
            for j in range(0,distance.shape[0]):
                plt.text(j,i,"%.1f" % distance[i][j],ha="center",va="center",color="w",fontsize=6)
    if show:
        plt.show()

def generateUniformDataset(n,dimensions=2):
    # n uniformly distributed points in d dimensions
    points = np.random.uniform(0,1,(n,dimensions))
    return points

def generateClusterDataset(num_clusters,points_per_cluster,dimensions=2):
    # generate cluster centers then points around those clusters
    points = np.zeros((num_clusters*points_per_cluster,dimensions))
    for c in range(0,num_clusters):
        center = np.random.uniform(0,1,(1,dimensions))
        for p in range(0,points_per_cluster):
            points[c*points_per_cluster+p] = center + np.random.normal(0,0.05,(1,dimensions))
    return points

def example():
    from time import time
    # generate clusters of points
    points = generateClusterDataset(6,50,2)

    # calculate the distance between every pair of points
    distance = calculateDistanceMatrix(points)

    # optionally illustrate how to deal with NAN values
    if False:
        for i in range(0,distance.shape[0]):
            for j in range(0,distance.shape[0]):
                if random.random() < 0.01:
                    distance[i][j] = float('nan')
                    distance[j][i] = float('nan')
        distance = imputeNANValues(distance)

    # get TSP ordering; set random number seed so that we get the same order every time we use the same dataset
    # TSP splits up clusters more often than tpPL, leaving a checkerboard pattern on the heat map
    startTime = time()
    TSPorder = treePenalizedPathLength(distance,20,seed=1,penaltyStrength=0.0)
    print('TSP            took %6.3f seconds to order %d points' % (time()-startTime,distance.shape[0]))
    TSPorder = standardOrder(distance,TSPorder)
    TSPdistance = reorderSymmetricMatrix(distance,TSPorder)

    # get tpPL ordering; set random number seed so that we get the same order every time we use the same dataset
    startTime = time()
    tpPLorder = treePenalizedPathLength(distance,20,seed=1,penaltyStrength=1.0)
    print('tpPL           took %6.3f seconds to order %d points' % (time()-startTime,distance.shape[0]))
    tpPLorder = standardOrder(distance,tpPLorder)
    tpPLdistance = reorderSymmetricMatrix(distance,tpPLorder)

    # optimal leaf order returns the same ordering each time
    # sometimes makes long jumps between clusters, which artificially makes the
    # clusters look more distinct than they really are
    startTime = time()
    OLOorder = optimalLeafOrder(distance)
    print('OLO            took %6.3f seconds to order %d points' % (time()-startTime,distance.shape[0]))
    # the following line fixes the long jumps in OLO, also deterministic, but sometimes slpits clusters
    # OLOorder, score_change = twoOptSwap(distance,OLOorder)
    # print('OLO plus 2-opt took %6.3f seconds to order %d points, changed score by %0.3f' % (time()-startTime,len(points),score_change))
    OLOorder = standardOrder(distance,OLOorder)
    OLOdistance = reorderSymmetricMatrix(distance,OLOorder)

    # plot data points in 2d, plot re-ordered heat map, side by side
    import matplotlib.pyplot as plt
    plt.figure(figsize=(13,8))
    plt.subplot(2,3,1)
    plotPoints(points[TSPorder],show=False)
    plt.ylabel("TSP",fontsize=10)
    plt.subplot(2,3,4)
    plotHeatMap(TSPdistance,show=False)
    plt.subplot(2,3,2)
    plotPoints(points[tpPLorder],show=False)
    plt.ylabel("tpPL",fontsize=10)
    plt.subplot(2,3,5)
    plotHeatMap(tpPLdistance,show=False)
    plt.subplot(2,3,3)
    plotPoints(points[OLOorder],show=False)
    plt.ylabel("OLO",fontsize=10)
    #plt.ylabel("OLO plus 2-opt",fontsize=10)
    plt.subplot(2,3,6)
    plotHeatMap(OLOdistance,show=False)
    plt.savefig('example.png')
    print("First point in the new orderings is shown as a green dot")
    plt.show()

if __name__ == "__main__":
    example()
