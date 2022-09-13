"""Program that generates random x,y,z coordinates to visualize the definition of the convex hull that is being used in NA_Pairwise_interactions.py"""
from draw_residues import draw_base
import matplotlib.pyplot as plt
import numpy as np
#from numpy.random import default_rng #not available in python 2.7 
from fr3d.classifiers.NA_pairwise_interactions import check_convex_hull_atoms 

def create_random_x_y_coordinates(size):
    rints_x = np.random.uniform(-5,5,size=size)
    rints_y = np.random.uniform(-5,5, size=size)
    rints_z = np.random.uniform(-5,5,size=size)
    return rints_x,rints_y,rints_z

def check_generated_coordinates(seq):
    xvalues = []
    yvalues = []
    zvalues = []

    samples = 100000
    x,y,z = create_random_x_y_coordinates(samples)
    print(x,y,z)
    for num in range(samples):
        check = check_convex_hull_atoms(x[num],y[num],z[num],seq)
        if check: 
            xvalues.append(x[num])
            yvalues.append(y[num])
            zvalues.append(z[num])
    return xvalues, yvalues, zvalues

def create_convex_hull_plot(seqlist):
    fig = plt.figure(figsize=(15.0, 6.0))
    #plt.title('Outline of Convex Hull of NA Bases',pad=20)
    sub=1
    for seq in seqlist:
        xvalues, yvalues, zvalues = check_generated_coordinates(seq)

        ax = fig.add_subplot(1, 4, sub)
        ax.axis("equal")
        ax.set_title('Random Plotted Points on Base %s' % seq) # %d %s %s' % (c,base_combination,interaction_list[0]), rotation=10)
        draw_base(seq,'default',2,ax)
        ax.scatter(xvalues,yvalues,color='c',marker=".",s=1)
        sub += 1

    #plt.savefig('/Users/jmitc/Documents/')
    plt.show()
    plt.close()

#Create Plots
seq=['A','C','G','U']
DNASeq = ['DA','DC','DG','DT']
create_convex_hull_plot(seq)
#create_convex_hull_plot(DNASeq)



