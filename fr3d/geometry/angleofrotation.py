import numpy as np


def angle_of_rotation(rotation_matrix):
    value = (np.trace(rotation_matrix) - 1.0) / 2.0
    value = np.clip(value, -1, 1)
    return np.arccos(value)


def axis_of_rotation(rotation_matrix):
    eigenvalues, eigenvector = np.linalg.eig(rotation_matrix)
    location = np.where(eigenvalues == 1)[0]
    return eigenvector[location]

def angle_between_planes(vec1, vec2):
    cosang = np.dot(vec1, vec2)
    sinang = np.linalg.norm(np.cross(vec1, vec2))
    angle = np.arctan2(sinang, cosang)
    return angle
