
from fr3d.definitions import RNAconnections
from fr3d.definitions import NAbasecoordinates
from fr3d.definitions import NAbasecolor


def draw_base(base_seq,dimensions,ax,zorder=1):
    """
    Connects atoms to draw one base in the specified number of dimensions
    """

    new_base_x = []
    new_base_y = []
    new_base_z = []

    for atomname in RNAconnections[base_seq]:
        coord_base= NAbasecoordinates[base_seq][atomname]
        new_base_x.append(coord_base[0])
        new_base_y.append(coord_base[1])
        new_base_z.append(coord_base[2])

    #print(new_base_x)

    if dimensions == 2:
        ax.plot(new_base_x, new_base_y, color=NAbasecolor[base_seq], linewidth=2.0, zorder=zorder)
    else:
        ax.plot(new_base_x, new_base_y, new_base_z, color=NAbasecolor[base_seq], linewidth=2.0)

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