import matplotlib.pyplot as plt
import numpy as np
from draw_residues import draw_base

#Taken from zStackingOverlap.m 
A = [(-2.100463, 0.447145), (0.000000,-1.009320), (4.271447,1.317924), (1.920945,5.150733),(0.230436,4.699718),(-2.100463,0.447145)] 
C = [(-2.082733,0.123632),(0.000000,-1.010259),(2.269450,-0.120783),(1.203833,4.411996),(-0.527970,4.602202),(-2.036772, 2.647095),(-2.082733,0.123632)]
G = [(-2.101572,0.463584),(0.000000,-1.009529),(4.872516,0.097781 ),(5.295175, 1.782283),(3.613335,3.374547),(1.396986,4.394213),(-0.751391,2.120132),(-2.101572,0.463584)]
U = [(-2.082780,0.111836),(0.000000,-1.008947),(2.292490, -0.048394 ),(2.092152,2.445179),(0.177156,4.020060),(-2.124577,2.616537),(-2.082780,0.111836)]
    

def check_convex_hull_atoms(x,y,z, parent):
    near_z_cutoff = 4.5
    inside = False
    if abs(z) < near_z_cutoff:
        if parent == 'MA':
                if  1.456465*x +  2.100463*y +  8.555926 > 0:  # Left of 0-1
                        if -2.327244*x +  4.271447*y + 11.525513 > 0:  # Left of 1-2
                                if -3.832809*x + -2.350502*y +  9.820982 > 0:  # Left of 2-3
                                        if  0.451015*x + -1.690509*y +  4.463818 > 0:  # Left of 3-4
                                                if  4.252573*x + -2.330899*y +  9.350088 > 0:  # Left of 4-5
                                                        inside = True
        if parent == 'MC':
                if  1.133891*x +  2.082733*y +  6.605970 > 0:  # Left of 0-1
                        if -0.889476*x +  2.269450*y +  6.389457 > 0:  # Left of 1-2
                                if -4.532773*x + -1.065617*y +  6.350554 > 0:  # Left of 2-3
                                        if -0.190212*x + -1.731803*y +  4.412786 > 0:  # Left of 3-4
                                                if  1.955107*x + -1.508802*y +  5.771434 > 0:  # Left of 4-5
                                                        if  2.523463*x + -0.045961*y +  6.131936 > 0:  # Left of 5-6
                                                                inside = True
        if parent == 'MG':
                if  1.473113*x +  2.101572*y +  8.853837 > 0:  # Left of 0-1
                        if -1.107310*x +  4.872516*y + 13.939528 > 0:  # Left of 1-2
                                if -1.684502*x +  0.422659*y +  6.635047 > 0:  # Left of 2-3
                                        if -1.592264*x + -1.681840*y +  5.439034 > 0:  # Left of 3-4
                                                if -1.019666*x + -2.216349*y +  4.841372 > 0:  # Left of 4-5
                                                        if  2.274081*x + -2.148377*y +  4.887645 > 0:  # Left of 5-6
                                                                if  1.656548*x + -1.350181*y +  3.573758 > 0:  # Left of 6-7         
                                                                        inside = True
        
        if parent == 'MU': 
                if  1.120783*x +  2.082780*y +  6.602568 > 0:  # Left of 0-1
                        if -0.960553*x +  2.292490*y +  6.551139 > 0:  # Left of 1-2
                                if -2.493573*x + -0.200338*y +  4.495078 > 0:  # Left of 2-3
                                        if -1.574881*x + -1.914996*y +  3.661149 > 0:  # Left of 3-4
                                                if  1.403523*x + -2.301733*y +  4.892567 > 0:  # Left of 4-5
                                                        if  2.504701*x +  0.041797*y +  6.112638 > 0:  # Left of 5-6
                                                                inside = True
    return inside
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
        seq = seq.replace("M", "")
        draw_base(seq,'default',2,ax)
        ax.scatter(xvalues,yvalues,color='c',marker=".",s=1)
        sub += 1

    plt.savefig('/Users/jmitc/Documents/')
    plt.show()
    plt.close()

#Create Plots
seq=['MA','MC','MG','MU']

DNASeq = ['DA','DC','DG','DT']
create_convex_hull_plot(seq)