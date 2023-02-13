
from fr3d.definitions import NAconnections
from fr3d.definitions import NAbasecoordinates
from fr3d.definitions import NAbasecolor

def draw_base(base_seq,colorscheme,dimensions,ax,zorder=1):
    """
    Connects atoms to draw one base in the specified number of dimensions
    base_seq is A, C, G, U, DA, DC, DG, DT
    colorscheme is 'default' or 'CPK'
    dimensions is 2 or 3
    ax is the current axis
    zorder is for 2d graphs to tell the plotting order
    """

    new_base_x = []
    new_base_y = []
    new_base_z = []

    if colorscheme == 'default':
        for atomname in NAconnections[base_seq]:
            coord_base = NAbasecoordinates[base_seq][atomname]
            new_base_x.append(coord_base[0])
            new_base_y.append(coord_base[1])
            new_base_z.append(coord_base[2])

        if dimensions == 2:
            ax.plot(new_base_x, new_base_y, color=NAbasecolor[base_seq], linewidth=2.0, zorder=zorder)
        else:
            ax.plot(new_base_x, new_base_y, new_base_z, color=NAbasecolor[base_seq], linewidth=2.0)

    else:
        for j in range(0,len(NAconnections[base_seq]),2):
            atom1 = NAconnections[base_seq][j]
            atom2 = NAconnections[base_seq][j+1]
            p = NAbasecoordinates[base_seq][atom1]
            q = NAbasecoordinates[base_seq][atom2]

            m = [(p[0]+q[0])/2,(p[1]+q[1])/2,(p[2]+q[2])/2] # midpoint

            atom_colors = {}
            atom_colors['C'] = [0.6,0.6,0.6]
            atom_colors['N'] = [0,0,1]
            atom_colors['O'] = [1,0,0]

            c1 = atom_colors[atom1[0]]
            c2 = atom_colors[atom2[0]]

            if dimensions == 2:
                ax.plot([p[0],m[0]],[p[1],m[1]], color=c1, linewidth=2.0, zorder=zorder)
                ax.plot([q[0],m[0]],[q[1],m[1]], color=c2, linewidth=2.0, zorder=zorder)
            else:
                ax.plot([q[0],m[0]],[q[1],m[1]][q[2],m[2]], color=NAbasecolor[base_seq], linewidth=2.0)



    return



    """
    from fr3d.definitions import Ribophos_connect

    for atomname in Ribophos_connect[base_seq]:
        back_base=[]
        back_base= basecoord_list[atomname]
        back_base_x.append(back_base[0])
        back_base_y.append(back_base[1])
        back_base_z.append(back_base[2])
    base_lines= ax.plot(back_base_x, back_base_y, back_base_z, label= 'Base')
    plt.setp(base_lines, 'color', 'g', 'linewidth', 1.0)
    #ax.text(9, 1, 1, base_residue)
    """