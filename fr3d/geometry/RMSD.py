import numpy


def RMSD(set1, set2):
    """This function calculates the root mean square difference of the
    atomic positions of set1 with set2.  In other words, this function
    calculates the root mean squared differences between two sets of a
    n 3-dimensional coordinates.  The calculation is explained here:
    wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
    """
    assert len(set1) == len(set2)
    L = len(set1)
    assert L > 0
    distance = numpy.sqrt(numpy.sum(numpy.power(set1 - set2, 2))/L)
    return distance


def sumsquarederror(set1, set2):
    """This function calculates the summed squared difference of the
    atomic positions of set1 with set2.  In other words, this function
    calculates the root mean squared differences between two sets of a
    n 3-dimensional coordinates.  The calculation is explained here:
    wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
    """
    assert len(set1) == len(set2)
    assert len(set1) > 0
    distance = numpy.sum(numpy.power(set1 - set2, 2))
    # TODO: Need to fix this for weighting.
    return distance
