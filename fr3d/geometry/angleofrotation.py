import numpy


def angle_of_rotation(rotation_matrix):
    return numpy.arccos((numpy.trace(rotation_matrix) - 1) / 2.0)


def axis_of_rotation(rotation_matrix):
    eigenvalues, eigenvector = numpy.linalg.eig(rotation_matrix)
    location = numpy.where(eigenvalues == 1)[0]
    return eigenvector[location]
