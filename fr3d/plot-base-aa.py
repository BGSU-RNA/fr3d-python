
def draw_base(base_seq, ax):
    """Connects atoms to draw neighboring bases and amino acids for 3D plots"""
     #creates lists of rotated base coordinates
    for basecoord_list in list_base_coord:
        new_base_x = []
        new_base_y = []
        new_base_z = []

        back_base_x = []
        back_base_y = []
        back_base_z = []


        try:
            for atomname in RNAconnections[base_seq]:
                coord_base = []
                coord_base= basecoord_list[atomname]
                new_base_x.append(coord_base[0])
                new_base_y.append(coord_base[1])
                new_base_z.append(coord_base[2])
            base_lines= ax.plot(new_base_x, new_base_y, new_base_z, label= 'Base')
            #ax.scatter(basecenter[0], basecenter[1], basecenter[2], zdir='y', color='b', marker='o')
            #ax.scatter(x = 0, y= 0, z= 0, color='b', marker='o')
            plt.setp(base_lines, 'color', 'b', 'linewidth', 1.0)

            for atomname in Ribophos_connect[base_seq]:
                back_base=[]
                back_base= basecoord_list[atomname]
                back_base_x.append(back_base[0])
                back_base_y.append(back_base[1])
                back_base_z.append(back_base[2])
            base_lines= ax.plot(back_base_x, back_base_y, back_base_z, label= 'Base')
            plt.setp(base_lines, 'color', 'g', 'linewidth', 1.0)
            #ax.text(9, 1, 1, base_residue)
        except:
            print "Missing residues"
            continue

def draw_aa(aa, ax):
    #Connects atoms to draw neighboring bases and amino acids for 3D plots
    for aacoord_list in list_aa_coord:
        new_aa_x=[]
        new_aa_y=[]
        new_aa_z=[]

        back_aa_x=[]
        back_aa_y=[]
        back_aa_z=[]

        try:
            for atomname in aa_connections[aa]:
                coord_aa=[]
                coord_aa= aacoord_list[atomname]
                new_aa_x.append(coord_aa[0])
                new_aa_y.append(coord_aa[1])
                new_aa_z.append(coord_aa[2])
            aa_lines= ax.plot(new_aa_x, new_aa_y, new_aa_z, label= 'Amino acid')
            plt.setp(aa_lines, 'color', 'r', 'linewidth', 1.0)

            for atomname in aa_backconnect[aa]:
                back_aa=[]
                back_aa= aacoord_list[atomname]
                back_aa_x.append(back_aa[0])
                back_aa_y.append(back_aa[1])
                back_aa_z.append(back_aa[2])
            aa_lines= ax.plot(back_aa_x, back_aa_y, back_aa_z, label= 'Amino acid')
            plt.setp(aa_lines, 'color', 'y', 'linewidth', 1.0)
        except:
            print "Missing residues"
            continue

def draw_aa_cent(aa, aa_part, ax):
    #Connects atoms to draw neighboring bases and amino acids for 3D plots
    for aacoord_list in list_aa_coord:
        new_aa_x=[]
        new_aa_y=[]
        new_aa_z=[]

        aa_center_x = 0
        aa_center_y = 0
        aa_center_z = 0
        n = 0

        if aa_part == 'aa_fg':
            connections = aa_connections
        elif aa_part == 'aa_backbone':
            connections = aa_backconnect
        try:
            for atomname in connections[aa]:
                coord_aa=[]
                coord_aa= aacoord_list[atomname]
                new_aa_x.append(coord_aa[0])
                new_aa_y.append(coord_aa[1])
                new_aa_z.append(coord_aa[2])

                aa_center_x = aa_center_x + coord_aa[0]
                aa_center_y = aa_center_y + coord_aa[1]
                aa_center_z = aa_center_z + coord_aa[2]
                n = n + 1
            ax.scatter(aa_center_x/n, aa_center_y/n, aa_center_z/n, c= 'r', marker = 'o')
        except:
            print "Missing residues"
            continue

def draw_base_aa_plots(list_base_aa):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for base, aa, interaction, edge in list_base_aa:
            base_seq = base.sequence
            aa= aa.sequence

            draw_base(base_seq, ax)
            draw_aa(aa, ax)
            #draw_aa_cent(aa, aa_part, ax)

            ax.set_xlabel('X Axis')
            ax.set_ylabel('Y Axis')
            ax.set_zlabel('Z Axis')
            ax.set_xlim3d(10, -15)
            ax.set_ylim3d(10, -15)
            ax.set_zlim3d(10, -15)
            plt.title('%s with ' % base_seq +'%s' % aa + ' %s' % aa_part)
            plt.show()

