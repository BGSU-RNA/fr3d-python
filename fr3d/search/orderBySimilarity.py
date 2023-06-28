import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import random
from myTimer import myTimer
# see https://matplotlib.org/users/pyplot_tutorial.html
# from tspy import TSP

def treePenalty(distance):
    Z = linkage(squareform(distance), "average")
#    print("regular",Z)

    penalty = np.zeros(distance.shape)

    group = []
    for i in range(0,distance.shape[0]):
        group.append([i])

    for merger in Z:
        a = int(merger[0])
        b = int(merger[1])
        group.append(group[a]+group[b])
        for i in group[a]:
            for j in group[b]:
                penalty[i][j] = merger[2]
                penalty[j][i] = merger[2]

#    print(group)

    if 0 > 1:
        for i in range(0,len(distance)):
            for j in range(0,len(distance)):
                print("%6.3f" % penalty[i][j]),
            print("")

    return penalty

def treePenalizedPathLength(distance,repetitions=100,seed=None):
    penalizedMatrix = distance + 5*treePenalty(distance)
    order = multipleGreedyInsertionPathLength(penalizedMatrix,repetitions,seed)
    return order

def optimalLeafOrder(distance):
    Z = linkage(squareform(distance), "average", optimal_ordering = True)
#    print("OLO",Z)
    dn = dendrogram(Z,no_plot=True)

    return dn['leaves']

def greedyInsertionPathLength(distance, order=[], verbose=False):

    # if no starting ordering
    if len(order) == 0:
        order = list(range(0,len(distance)))
        random.shuffle(order)          # random starting ordering
    path = order[:2]            # first two points of the current ordering
    score = distance[path[0], path[1]]

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

        path.insert(bestPosition, order[p])
        score += bestScore

    return path, score

def multipleGreedyInsertionPathLength(distance, repetitions=100, seed=None):

    if seed:
        random.seed(seed)

    bestScore = float("inf")
    for rep in range(0,repetitions):
        path, score = greedyInsertionPathLength(distance)
        if score < bestScore:
            bestScore = score
            bestPath = path

    return bestPath

def reorderSymmetricMatrix(distance, newOrder):

    newDistance = np.zeros(distance.shape)
    for i in range(0,newDistance.shape[0]):
        for j in range(i+1,newDistance.shape[1]):
           newDistance[i][j] = distance[newOrder[i]][newOrder[j]]
           newDistance[j][i] = newDistance[i][j]

    return newDistance

def generateUniformDataset(n,d):
    points = np.random.uniform(0,1,(n,d))

    distance = np.zeros((n,n))
    for i in range(0,n):
        for j in range(i+1,n):
            distance[i][j] = np.linalg.norm(points[i]-points[j])
            distance[j][i] = distance[i][j]

    return points, distance


def example():
    """
    """

    size = 20

def testRuntime():
    sizes = [10, 20, 40, 80, 160, 320, 640]
    timerData = myTimer("start")
    numReplications = 10

    print("orderBySimilarity: testRuntime")

    for rep in range(0,numReplications):
        for size in sizes:
            timerData = myTimer("generate_"+str(size))
            points, distance = generateUniformDataset(size,3)

            timerData = myTimer("treepenalty_"+str(size))
            penalizedMatrix = distance + 5*treePenalty(distance)

            timerData = myTimer("OLO_"+str(size))
            order = optimalLeafOrder(distance)

            timerData = myTimer("tpTSP_"+str(size))
            penalizedMatrix = distance + 5*penalizedMatrix
            order = multipleGreedyInsertionPathLength(distance,2)

        print(myTimer("summary"))

    methods = ["generate","treepenalty","OLO","tpTSP"]
    times = {}
    for method in methods:
        times[method] = []

    minmin = 0
    maxmax = 0
    for size in sizes:
        for method in methods:
            newTime = timerData[method+"_"+str(size)]/numReplications
            times[method].append(newTime)
            maxmax = max(maxmax,newTime)
            minmin = min(minmin,newTime)

    import matplotlib.pyplot as plt

    colors = {}
    colors["generate"] = "r"
    colors["OLO"] = "b"
    colors["tpTSP"] = "g"
    colors["treepenalty"] = "k"

    for method in methods:
        plt.plot(sizes, times[method], colors[method])
        print("%s slope %8.4f" % (method,(np.log(times[method][-1])-np.log(times[method][-2]))/(np.log(sizes[-1])-np.log(sizes[-2]))))

    minmin = max(0.001,minmin)  # small numbers don't work

    plt.axis([min(sizes), max(sizes), minmin, maxmax])
    plt.xscale('log')                 # logarithmic scale to squeeze in large maximum values
    plt.yscale('log')                 # logarithmic scale to squeeze in large maximum values
    plt.xlabel('Number of points')
    plt.ylabel('Run time')
    plt.title('Log log plot of runtime versus number of points')
    plt.show()           # use plt.show() when not using repl.it

#print("testing runtime")
#testRuntime()
