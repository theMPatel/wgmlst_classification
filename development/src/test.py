
def calcname(named, tree, unnamed, distances, thresholds):

    thresholds.sort( key=lambda x: -x)



    level = 1
    while len(distances):

        closestClusters = set()
        while thresholds[level-1] >= distances[0] > thresholds[level]:

            