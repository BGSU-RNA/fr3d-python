import numpy

def angle_of_rotation(rotation_matrix):  
                   
    angle = numpy.arccos((numpy.trace(rotation_matrix)-1)/2.0)*180/numpy.pi
    
    return angle
    
def axis_of_rotation(rotation_matrix):
    eigenvalues, eigenvector = numpy.linalg.eig(rotation_matrix)

    location = numpy.where(eigenvalues==1)[0]
    axis = eigenvector[location]

    return axis
    
    