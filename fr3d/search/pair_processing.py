"""
Process one point on each unit to make lists of pairs of units that could possibly
be in the same motif instance.  Units that are too far apart can be excluded.
"""

import numpy as np
from collections import defaultdict


def get_pairlist(Q, models, centers, alt_index = False):
    """
    Use a fast technique to find all pairs of units whose centers have
    distances within given distance ranges.
    For each pair of positions, return a list of indices whose distances
    are within the min to max range for that pair of positions.
    """

    pairlist = defaultdict(dict)

    if Q["type"] == "mixed" or Q["type"] == "geometric":

        if centers.size > 0:
            min_coord = np.nanmin(centers)
            max_coord = np.nanmax(centers)

            # square_size is used to map cubes to integers; it's a nice integer larger than the diameter of the coordinates
            square_size = np.power(10, np.floor(np.log10(max_coord - min_coord))) * np.ceil((max_coord-min_coord)/np.power(10,np.floor(np.log10(max_coord- min_coord))))

            # get all triples of (index1,index2,distance) where distance is below Q["largestMaxRange"]
            max_pair_list = fixed_radius_search(len(centers), models, centers, Q["largestMaxRange"], square_size)

            # sort these triples by distance
            max_pair_list = sorted(max_pair_list, key = lambda x : x[2])

            # we want to return indices in both orders, so make such a list once, then pull slices from it later
            # that's much, much quicker than building lists of pairs over and over again
            both_pairs_list = []
            for a,b,c in max_pair_list:
                both_pairs_list.append((a,b))
                both_pairs_list.append((b,a))

            if alt_index == False:
                for i in range(0,Q['numpositions']):
                    for j in range(i+1,Q['numpositions']):
                        min_squared = Q["requireddistanceminimum"][i][j]**2
                        max_squared = Q["requireddistancemaximum"][i][j]**2

                        # find pair with lowest distance in the current min-max range
                        if Q["requireddistanceminimum"][i][j] <= 0:
                            # this happens sometimes, then it's clear where to start
                            pair_start = 0
                        else:
                            a = 0                    # lower limit
                            b = len(max_pair_list)   # upper limit
                            while a < b:             # far apart, haven't converged yet
                                c = (a+b)//2         # midpoint of current indices, rounded down

                                if max_pair_list[c][2] > min_squared:
                                    b = c
                                else:
                                    a = c+1

                            pair_start = b

                        #print(i,j,a,c,b,max_pair_list[b][2]>min_squared)

                        # find pair with largest distance in the current min-max range
                        a = 0
                        b = len(max_pair_list)-1
                        while a < b:
                            c = (a+b+1)//2         # midpoint of current indices, rounded
                            if max_pair_list[c][2] < max_squared:
                                a = c
                            else:
                                b = c-1

                        pair_end = b

                        #print(i,j,a,c,b,max_pair_list[b][2]<max_squared)

                        pairlist[i][j] = both_pairs_list[(2*pair_start):(2*pair_end+2)]

            else:
                for i in alt_index:
                    for j in alt_index[alt_index.index(i)+1:]:
                        min_squared = Q["requireddistanceminimum"][alt_index.index(i)][alt_index.index(j)]**2
                        max_squared = Q["requireddistancemaximum"][alt_index.index(i)][alt_index.index(j)]**2
                        pair_list = []
                        for k in max_pair_list:
                            if k[2] > min_squared and k[2] < max_squared:
                                pair_list.append((alt_index[k[0]], alt_index[k[1]]))
                                pair_list.append((alt_index[k[1]], alt_index[k[0]]))
                        pairlist[i][j] = pair_list

    else:
        if alt_index == False:
            for i in range(0,Q['numpositions']):
                for j in range(i+1,Q['numpositions']):
                    pairlist[i][j] = "full"
        else:
            for i in alt_index:
                for j in alt_index[alt_index.index(i)+1:]:
                    pairlist[i][j] = "full"

    return pairlist


def fixed_radius_search(num_points, models, positions, max_cutoff, square_size):
    """
    Find all pairs of points with the same model and whose distance
    between positions is less than max_cutoff.
    Put points into buckets, then calculate squared distances between points in
    adjacent buckets.
    Points from different models cannot make a pair.
    m = number of points
    models = list of strings
    positions = mx3 numpy array
    max_cutoff = number
    square_size = length of each side of the cubes to split into
    Returns a list of points and the squared distances between them
    """

    max_size = square_size
    x = positions[:,0]
    y = positions[:,1]
    z = positions[:,2]
    d = max_cutoff      # subsquare side length & distance you want to find from
    dSquared = d*d

    bucket = 0
    bucket_list = {}

    grid_len = int(np.ceil(max_size/(d))+2)

    # index buckets by an integer for x,y,z and by model number as a string
    for i in range(num_points):
        # skip units like GLY that don't have a center defined for the sidechain
        if not np.isnan(x[i]) and not np.isnan(y[i]) and not np.isnan(z[i]):
            bucket = (int( (x[i]//d) + (y[i]//d)*grid_len + (z[i]//d) * grid_len**2) + 1 + grid_len + grid_len**2, models[i])
            if bucket in bucket_list:
                bucket_list[bucket].append(i)
            else:
                bucket_list[bucket] = [i]

    # adjacent buckets are at the same z level or above
    # don't need to check all 27 buckets, that would check pairs twice
    dx = [1,0,-1,1,0,-1,1,0,-1,1,0,-1,1,0]
    dy = [1,1,1,0,0,0,-1,-1,-1,1,1,1,0,0]
    dz = [1,1,1,1,1,1,1,1,1,0,0,0,0,0]

    pairlist = []

    # loop over all buckets
    for first_bucket in bucket_list:
        # loop over adjacent buckets
        for j in range(13):
            # second bucket must be the same or adjacent, and must have the same model number
            second_bucket = (first_bucket[0] + dx[j] + grid_len*dy[j] + (grid_len**2)*dz[j],first_bucket[1])
            # check if there are any points in the adjacent bucket
            if second_bucket in bucket_list:
                # loop over points in the first bucket
                for pnt1 in bucket_list[first_bucket]:
                    # loop over points in the second bucket
                    for pnt2 in bucket_list[second_bucket]:
                        distanceSquared = (x[pnt1]-x[pnt2])**2 + (y[pnt1]-y[pnt2])**2+ (z[pnt1] - z[pnt2])**2 #dist(i,pnt2, x,y,z)
                        if distanceSquared < dSquared:
                            pairlist.append((pnt1,pnt2,distanceSquared))

        # loop over pairs of points in the same bucket
        for pnt1 in bucket_list[first_bucket]:
            # loop over points in the same bucket
            for pnt2 in bucket_list[first_bucket]:
                # only check in one direction
                if pnt1 < pnt2:
                    distanceSquared = (x[pnt1]-x[pnt2])**2 + (y[pnt1]-y[pnt2])**2+ (z[pnt1] - z[pnt2])**2 #dist(i,pnt2, x,y,z)
                    if distanceSquared < dSquared:
                        pairlist.append((pnt1,pnt2,distanceSquared))

    return pairlist


def get_pairlist_old(Q, models, centers, alt_index = False):
    """
    Deprecated and replaced with a newer version.

    Use a fast technique to find all pairs of units whose centers have
    distances within given distance ranges.
    This method uses append to build many short lists, but that is slow,
    so this method is decrecated.
    """

    pairlist = defaultdict(dict)

    if Q["type"] == "mixed" or Q["type"] == "geometric":

        if centers.size > 0:
            min_coord = np.nanmin(centers)
            max_coord = np.nanmax(centers)

            square_size = np.power(10, np.floor(np.log10(max_coord - min_coord))) * np.ceil((max_coord-min_coord)/np.power(10,np.floor(np.log10(max_coord- min_coord))))
            # get all triples of (index1,index2,distance) where distance is below Q["largestMaxRange"]
            max_pair_list = fixed_radius_search(len(centers), models, centers, Q["largestMaxRange"], square_size)

            if alt_index == False:
                for i in range(0,Q['numpositions']):
                    for j in range(i+1,Q['numpositions']):
                        min_squared = Q["requireddistanceminimum"][i][j]**2
                        max_squared = Q["requireddistancemaximum"][i][j]**2
                        pair_list = []
                        for k in max_pair_list:
                            if k[2] > min_squared and k[2] < max_squared:
                                pair_list.append((k[0], k[1]))
                                pair_list.append((k[1], k[0]))
                        pairlist[i][j] = pair_list
            else:
                for i in alt_index:
                    for j in alt_index[alt_index.index(i)+1:]:
                        min_squared = Q["requireddistanceminimum"][alt_index.index(i)][alt_index.index(j)]**2
                        max_squared = Q["requireddistancemaximum"][alt_index.index(i)][alt_index.index(j)]**2
                        pair_list = []
                        for k in max_pair_list:
                            if k[2] > min_squared and k[2] < max_squared:
                                pair_list.append((alt_index[k[0]], alt_index[k[1]]))
                                pair_list.append((alt_index[k[1]], alt_index[k[0]]))
                        pairlist[i][j] = pair_list

    else:
        if alt_index == False:
            for i in range(0,Q['numpositions']):
                for j in range(i+1,Q['numpositions']):
                    pairlist[i][j] = "full"
        else:
            for i in alt_index:
                for j in alt_index[alt_index.index(i)+1:]:
                    pairlist[i][j] = "full"

    return pairlist
