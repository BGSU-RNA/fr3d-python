"""
refine_atom_mappings.py reads atom_mappings_provisional.txt,
superimposes the modified base on the parent base,
matches each parent base atom to the nearest modified base atom,
writes out those mappings as atom_mappings.txt,
and if there is no .png image file, it writes out a visual representation of the mappings

Review the new image files and if there is something wrong,
add lines to atom_mappings_manual.txt to tell it what to do.
Delete the image file, then run make_atom_mappings.py again and then refine_atom_mappings again
Iterate until you've got it right.
Use https://www.rcsb.org/ligand/A but substitute the non-standard nucleotide where A is

Ideas for the next version:
    Color the base letters red if the atom changes
    Choose colors closer to CPK for the base atoms; need enough reds, blues, grays.  Get close-ish
    Color connection to C1' to acknowledge that it is a connection, but don't color what C1' connects to, too complicated to view
"""

draw_figures = True      # draw new figures if they don't already exist
draw_figures = False     # don't draw new figures at all

overwrite_figures = False
overwrite_figures = True  # draw figures, overwriting existing ones

plot_standard = True     # include the standard base in the plots
plot_standard = False    # just plot the modified nucleotide

save_as_gif = False
save_as_gif = True

crop_out_white_space = False
crop_out_white_space = True

color_scheme = 'diagnostic'  # use many colors, to show the atom mappings
color_scheme = 'CPK'         # use CPK coloring

atom_label_color = 'black'
new_atom_point_color = '#80D1E3'   # blue of Argon since that probably won't be added

show_figure = True       # pause to show each figure, enable rotation of coordinates
show_figure = False

focus_list = ['YYG','MHG','1W5','8NI','TJU']
focus_list = ['1SC']
focus_list = ['MA6']
focus_list = ['OMC']
focus_list = ['G7M']
focus_list = ['XGA']
focus_list = ['1SC']
focus_list = ['0A','0C','0G','0U']
focus_list = ['F86']
focus_list = ['1CC']
focus_list = ['1FC']
focus_list = ['2GF']
focus_list = ['2DT']
focus_list = ['6MZ']
focus_list = ['5HC']
focus_list = ['6HC']
focus_list = ['C2S']
focus_list = ['0A','0C','0G','0U']  # chirality changes galore
focus_list = ['XB9']
focus_list = ['3DR']
focus_list = ['C66','E3C']
focus_list = ['16B','45A','7AT','A2P','A3P','ADS','ANC','APC','ATP','PPS']  # A updates
focus_list = ['4AC','5IC','6OO','73W','A5M','AI5','LCC','N7X','O2C']  # C updates
focus_list = ['02I','2BU','2DA','3DA','4EN','6HA','6HB','7DA','A3A','AD2','AF2','FA2','FAX','L3X','MDV','RMP','SMP','TFO','URT','XAD']  # DA updates
focus_list = []    # process all modified nucleotides

max_count = 99999      # process modified nucleotides with count at or below this number

if color_scheme == 'CPK':
    output_directory = "C:/users/zirbel/Documents/modified"
else:
    output_directory = "C:/users/zirbel/Documents/modified/diagnostic"

# from definitions import NAconnections
from definitions import NAbasecoordinates
from definitions import NAbaseheavyatoms
from definitions import NAbasehydrogens
from superpositions import besttransformation
from make_atom_mappings import read_monomer_cif

import imageio
import json
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import sys
import time
from PIL import Image

if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
    read_mode = 'rb'
    write_mode = 'w'
else:
    from urllib.request import urlretrieve as urlretrieve
    read_mode = 'rt'
    write_mode = 'wt'   # write as text


def element_to_cpk_color(element):
    # input is a text string like "C" or "O" or "BR"
    # that is an element symbol
    # output is a list of three numbers between 0 and 1 showing the red, green, blue colors

    # what other elements are needed?
    # there are 110 elements, but only 18 are common in nucleic acids
    # common elements in nucleic acids
    # C, N, O, P, S, H, F, Cl, Br, I, Mg, Na, K, Ca, Mn, Fe, Co, Zn

    # best to find the CPK colors used by Protein Data Bank (PDB) or by Mol* visualization program
    # Source of the hexadecimal string for each color: https://sciencenotes.org/molecule-atom-colors-cpk-colors/
    # Easier to read: https://jmol.sourceforge.net/jscolors/

    # rgb = [0.0,0.0,0.0]  # black, default
    rgb = None

    el = element.upper()

    if el == "H":
        rgb = [0.9,0.9,0.9] # light gray because we have a white background
        rgb = [0.8,0.8,0.8] # light gray because we have a white background
    elif el == "B":
        rgb = "FFB5B5"
    elif el == "C":
        rgb = '#909090' # gray
        rgb = '#808080' # darker gray
    elif el == "N":
        rgb = '#3050F8' # blue
    elif el == "O":
        rgb = '#FF0D0D' # red
    elif el == "F":
        rgb = '#90E050'
    elif el == "NA":
        rgb = '#AB5CF2'
    elif el == "MG":
        rgb = '#8AFF00'
    elif el == "P":
        rgb = '#FF8000'
    elif el == "S":
        rgb = '#FFFF30'
    elif el == "CL":
        rgb = '#1FF01F'
    elif el == "CA":
        rgb = '#3DFF00'
    elif el == "V":
        rgb = 'A6A6AB'
    elif el == "MN":
        rgb = '#9C7AC7'
    elif el == "FE":
        rgb = '#E06633'
    elif el == "CO":
        rgb = '#F090A0'
    elif el == "ZN":
        rgb = '#7D80B0'
    elif el == "K":
        rgb = '#8F40D4'
    elif el == "SE":
        rgb = "FFA100"
    elif el == "BR":
        rgb = '#A62929'
    elif el == "TE":
        rgb = "D47A00"
    elif el == "I":
        rgb = '#940094'
    elif el == "PT":
        rgb = "D0D0E0"

    if not rgb:
        print("Unknown element %s" % el)
        print(crashnow)

    return rgb


def get_cif_data(base):
    """
    Read the .cif file and organize its data into data structures for this program
    """

    coordinates = {}
    coordinates_ideal = {}
    connections = []
    atom_to_element = {}
    atom_to_chirality = {}
    par_comp_id = 'No par_comp_id line'
    one_letter_code = 'No one_letter_code line'

    cif_data = read_monomer_cif(base)

    par_comp_id = cif_data['chem_comp'][0].get('mon_nstd_par_comp_id',None)
    one_letter_code = cif_data['chem_comp'][0].get('one_letter_code',None)

    # apparently this does not happen
    # if par_comp_id and "," in par_comp_id:
    #     print("Multiple parent compound IDs %s for %s" % (par_comp_id,base))
    #     print(crashnow)

    for row in cif_data['chem_comp_atom']:
        atom = row['atom_id']
        x = row['model_Cartn_x']
        y = row['model_Cartn_y']
        z = row['model_Cartn_z']

        if x and y and z:
            coordinates[atom] = [float(x),float(y),float(z)]
        atom_to_element[atom] = row['type_symbol']

        x = row['pdbx_model_Cartn_x_ideal']
        y = row['pdbx_model_Cartn_y_ideal']
        z = row['pdbx_model_Cartn_z_ideal']

        if x and y and z:
            coordinates_ideal[atom] = [float(x),float(y),float(z)]

        chirality = row.get('pdbx_stereo_config',None)
        atom_to_chirality[atom] = chirality


    # not sure why but sometimes the ideal coordinates are more complete
    if len(coordinates_ideal) > len(coordinates):
        coordinates = coordinates_ideal

    # see if this fixes some H5' and H5'' labeling discrepancies
    # if base in ['5HC']:
    #     coordinates = coordinates_ideal

    for row in cif_data['chem_comp_bond']:
        atom1 = row['atom_id_1']
        atom2 = row['atom_id_2']
        connections.append((atom1,atom2))
        connections.append((atom2,atom1))

    connections = list(set(connections))

    return coordinates, connections, atom_to_element, atom_to_chirality, par_comp_id, one_letter_code


def my_norm(x,y):

    return math.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2)


def read_atom_mappings(filename):
    """
    Read a text file of four columns
    (parent,parent atom,modified,modified atom)
    """
    parent_to_modified_atom = {}
    modified_to_parent_atom = {}
    modified_base_to_parent = {}
    not_mappable = []

    with open(filename,read_mode) as f:
        lines = f.readlines()

    # store atom mappings in dictionaries to make it easy to map them around
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if len(fields) >= 3:    # allow for no trailing tabs on some lines
            parent = fields[0]
            modified = fields[2]

            if len(parent) == 0:
                not_mappable.append(modified)
                continue

            if not modified in modified_to_parent_atom:
                modified_to_parent_atom[modified] = {}
                parent_to_modified_atom[modified] = {}
                modified_base_to_parent[modified] = parent

            parent_atom = fields[1]
            if parent_atom:
                modified_atom = fields[3]
                modified_to_parent_atom[modified][modified_atom] = parent_atom
                parent_to_modified_atom[modified][parent_atom] = modified_atom

    return parent_to_modified_atom, modified_to_parent_atom, modified_base_to_parent, not_mappable


def unit_vector(v):
    return v / np.linalg.norm(v)


def pyramidal_hydrogens(P1,C,P2,bondLength=1.1):
    # return positions of hydrogens making a tetrahedron with center C and vertices P1 and P2

    # infer positions one way
    V1 = P1
    V2 = P2
    # vector from V2 to C
    u = unit_vector(C-V2)
    # construct Rodrigues rotation matrix
    # matrix to rotate 120 degrees around vector u
    W = np.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])
    R = np.identity(3) + (np.sqrt(3)/2)*W + 1.5 * np.dot(W,W)
    # remaining vertices are vector from C to V1 rotated 120 degrees in either direction
    V3 = C + bondLength * unit_vector(np.dot(R,V1-C))
    V4 = C + bondLength * unit_vector(np.dot(np.transpose(R),V1-C))

    # infer positions the other way
    V1 = P2
    V2 = P1
    # vector from V2 to C
    u = unit_vector(C-V2)
    # construct Rodrigues rotation matrix
    # matrix to rotate 120 degrees around vector u
    W = np.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])
    R = np.identity(3) + (np.sqrt(3)/2)*W + 1.5 * np.dot(W,W)
    # remaining vertices are vector from C to V1 rotated 120 degrees in either direction
    VV4 = C + bondLength * unit_vector(np.dot(R,V1-C))
    VV3 = C + bondLength * unit_vector(np.dot(np.transpose(R),V1-C))

    # average the two inferred positions
    P3 = (V3+VV3)/2
    P4 = (V4+VV4)/2

    return P3, P4


def get_mod_atom_closest_to(atom_list, parent_to_modified_atom, mod_coordinates, bond_length=1.1):
    # use tetrahedral geometry to infer the location of the hydrogen atoms
    # zzz

    c = []
    for a in atom_list:
        if a in parent_to_modified_atom:
            b = parent_to_modified_atom[a]
            c.append(np.array(mod_coordinates[b]))
        else:
            return None, None, None, None

    # compute location of H5' from C4', C5', O5', for example
    p, q = pyramidal_hydrogens(c[0],c[1],c[2],bond_length)

    mapped_to = set(parent_to_modified_atom.values())

    min_dist = 1
    min_atom = None
    # min_d = None
    for b in mod_coordinates.keys():
        if not b in mapped_to:
            d = np.array(mod_coordinates[b])
            dist = np.linalg.norm(p-d)
            if dist < min_dist:
                min_dist = dist
                min_atom = b
                # min_d = d

    # if min_atom:
    #     print('Hydrogen bond length is %8.4f' % np.linalg.norm(c[1]-min_d))

    return min_atom, min_dist, p, q


def draw_base_coordinates(base_seq,coordinates,connections,atom_to_display,ax,limits=None):
    """
    Connects atoms to draw one base
    ax is the current axis
    """

    if limits:
        xmin,xmax,ymin,ymax = limits
        # print("Limits:")
        # print(xmin,xmax,ymin,ymax)
    else:
        xmin = 100
        xmax = -100
        ymin = 100
        ymax = -100

    connections = sorted(connections)

    drawn_connections = set([])
    drawn_atoms = set([])

    if backbone:
        fs = 8     # font size for labels
        ds = 6     # display size for points
    else:
        fs = 16
        ds = 16

    for atom1,atom2 in connections:

        if atom1 in coordinates:
            p = coordinates[atom1]
        else:
            continue

        if atom2 in coordinates:
            q = coordinates[atom2]
        else:
            continue

        drawn_connections.add((atom2,atom1))
        drawn_connections.add((atom1,atom2))

        # get (bond color, point color, display mode) for each atom
        # if no color is specified, use green for the bond
        c1,pc1,d1 = atom_to_display.get(atom1,("green",new_atom_point_color,"not"))
        c2,pc2,d2 = atom_to_display.get(atom2,("green",new_atom_point_color,"not"))

        # if d1 == "thin" and d2 == "thin":
        #     # this will plot just one atom away from the mapped atoms
        #     continue

        # if one of the atoms is not to be plotted, skip this connection
        if d1 == "not" or d2 == "not":
            continue

        sh = 0.05  # shift the atom label slightly to the right and up

        if atom1 in atom_to_display and not atom1 in drawn_atoms:
            ax.text(p[0]+sh,p[1]+sh,p[2],atom1,fontsize=fs,color=atom_label_color)
            ax.scatter([p[0]],[p[1]],[p[2]], s=ds, color=pc1)
            drawn_atoms.add(atom1)

        if atom1 in atom_to_display and not atom2 in drawn_atoms:
            # plot atoms that are connected to a mapped atom
            ax.text(q[0]+sh,q[1]+sh,q[2],atom2,fontsize=fs,color=atom_label_color)
            ax.scatter([q[0]],[q[1]],[q[2]], s=ds, color=pc2)
            drawn_atoms.add(atom2)

        if atom2 in atom_to_display and not atom2 in drawn_atoms:
            ax.text(q[0]+sh,q[1]+sh,q[2],atom2,fontsize=fs,color=atom_label_color)
            ax.scatter([q[0]],[q[1]],[q[2]], s=ds, color=pc2)
            drawn_atoms.add(atom2)

        if atom2 in atom_to_display and not atom1 in drawn_atoms:
            # plot atoms that are connected to a mapped atom
            ax.text(p[0]+sh,p[1]+sh,p[2],atom1,fontsize=fs,color=atom_label_color)
            ax.scatter([p[0]],[p[1]],[p[2]], s=ds, color=pc1)

        if base_seq in ['A','C','G','U','DA','DC','DG','DT'] or (atom1 in atom_to_display and atom2 in atom_to_display):
            # count all standard base atoms toward the size of the frame, also color atoms on modified
            xmin = min(xmin,p[0],q[0])
            xmax = max(xmax,p[0],q[0])
            ymin = min(ymin,p[1],q[1])
            ymax = max(ymax,p[1],q[1])

        if d1 == "full" and d2 == "full":
            if backbone:
                lw = 6.0
                lw = 4.0
            else:
                lw = 12.0
                lw = 6.0

        else:
            # draw unmapped atoms that are connected to mapped atoms
            # draw thin lines

            if backbone:
                lw = 2.0
                lw = 1.0
            else:
                lw = 4.0
                lw = 2.0

        # color half of the bond by atom1, half of the bond by atom2
        s = [(3*p[0]+q[0])/4.0,(3*p[1]+q[1])/4.0,(3*p[2]+q[2])/4.0] # quarter point, near atom1
        m = [(p[0]+q[0])/2.0,(p[1]+q[1])/2.0,(p[2]+q[2])/2.0] # midpoint
        t = [(p[0]+3*q[0])/4.0,(p[1]+3*q[1])/4.0,(p[2]+3*q[2])/4.0] # three-quarters point, near atom2

        # points in order from atom1 to atom2 are p, s, m, t, q
        try:
            if c1 == c2:
                # draw from atom1 to atom2 with projection to make lines cover the atoms
                ax.plot3D([p[0],q[0]],[p[1],q[1]],[p[2],q[2]], color=c1, linewidth=lw, solid_capstyle='round')

            else:
                # draw from atom1 to quarter way across with projection to make lines cover the atom
                ax.plot3D([p[0],s[0]],[p[1],s[1]],[p[2],s[2]], color=c1, linewidth=lw, solid_capstyle='round')
                # draw from quarter way to midpoint with ends that don't overlap at the midpoint
                ax.plot3D([s[0],m[0]],[s[1],m[1]],[s[2],m[2]], color=c1, linewidth=lw, solid_capstyle='butt')
                # draw from midpoint to three-quarters with ends that don't overlap at the midpoint
                ax.plot3D([m[0],t[0]],[m[1],t[1]],[m[2],t[2]], color=c2, linewidth=lw, solid_capstyle='butt')
                # draw from three-quarters to atom2 with projection to make lines cover the atoms
                ax.plot3D([t[0],q[0]],[t[1],q[1]],[t[2],q[2]], color=c2, linewidth=lw, solid_capstyle='round')

                # # draw from atom1 to quarter way across with projection to make lines cover the atom
                # ax.plot3D([p[0],s[0]],[p[1],s[1]],[p[2],s[2]], color=c1, linewidth=lw)
                # # draw from quarter way to midpoint with ends that don't overlap at the midpoint
                # ax.plot3D([s[0],m[0]],[s[1],m[1]],[s[2],m[2]], color=c1, linewidth=lw)
                # # draw from midpoint to three-quarters with ends that don't overlap at the midpoint
                # ax.plot3D([m[0],t[0]],[m[1],t[1]],[m[2],t[2]], color=c2, linewidth=lw)
                # # draw from three-quarters to atom2 with projection to make lines cover the atoms
                # ax.plot3D([t[0],q[0]],[t[1],q[1]],[t[2],q[2]], color=c2, linewidth=lw)
                # pass

            xmin = min(xmin,p[0],q[0])
            xmax = max(xmax,p[0],q[0])
            ymin = min(ymin,p[1],q[1])
            ymax = max(ymax,p[1],q[1])

        except:
            print("Trouble drawing %s to %s" % (atom1,atom2))
            print(p)
            print(q)
            print(lw)
            print(c1)
            print(c2)
            continue

    return xmin, xmax, ymin, ymax


def keep_unique_changes(modified_changes):
    """
    Remove duplicate changes from the list of changes
    """

    view_order = {}
    view_order['replacement'] = 0
    view_order['addition'] = 1
    view_order['removal'] = 2
    view_order['chirality'] = 3
    view_order['chirality_reversal'] = 4
    view_order['added_bond'] = 5
    view_order['removed_bond'] = 6


    modified_changes['changes'] = sorted(modified_changes['changes'], key = lambda x : (x['change_location'],view_order[x['change_type']],x.get('parent_atom',""),x.get('modified_atom',''),x.get('new_modified_atom','')))

    seen = set()
    keep = []
    for change in modified_changes['changes']:
        s = str(change)
        if not s in seen:
            seen.add(s)
            keep.append(change)

    modified_changes['changes'] = keep

    return modified_changes


def my_bbox(image):
    # Find the box that bounds the non-white pixels of the image

    # Get the width and height of the image
    width, height = image.size

    # Initialize the bounding box coordinates
    min_x = width
    min_y = height
    max_x = 0
    max_y = 0

    # Iterate over each pixel in the image
    # left to right
    x = 0
    while x < width and min_x == width:
        for y in range(height):
            # Get the pixel value at the current position
            pixel = image.getpixel((x, y))

            # Check if the pixel is not white
            if min(pixel) < 255:
                # Update the bounding box coordinates
                min_x = x
                break
        x += 1

    # right to left
    x = width-1
    while x > 0 and max_x == 0:
        for y in range(height):
            # Get the pixel value at the current position
            pixel = image.getpixel((x, y))

            # Check if the pixel is not white
            if min(pixel) < 255:
                # Update the bounding box coordinates
                max_x = x
                break
        x -= 1

    # top to bottom
    y = 0
    while y < height and min_y == height:
        for x in range(min_x,max_x):
            # Get the pixel value at the current position
            pixel = image.getpixel((x, y))

            # Check if the pixel is not white
            if min(pixel) < 255:
                # Update the bounding box coordinates
                min_y = y
                break
        y += 1

    # bottom to top
    y = height-1
    while y > 0 and max_y == 0:
        for x in range(min_x,max_x):
            # Get the pixel value at the current position
            pixel = image.getpixel((x, y))

            # Check if the pixel is not white
            if min(pixel) < 255:
                # Update the bounding box coordinates
                max_y = y
                break
        y -= 1

    # Expand the bounding box by one pixel in each direction
    min_x = max(0,min_x-1)
    max_x = min(width,max_x+1)
    min_y = max(0,min_y-1)
    max_y = min(height,max_y+1)

    # Return the bounding box coordinates
    return (min_x, min_y, max_x, max_y)


#######################################################
# main block starts here

standard_nts = ['A','C','G','U','DA','DC','DG','DT']

# read provisional atom to atom mappings
parent_to_modified_atom, modified_to_parent_atom, modified_base_to_parent, not_mappable = read_atom_mappings("atom_mappings_provisional.txt")

# read manual mappings
parent_to_modified_atom_manual, modified_to_parent_atom_manual, modified_base_to_parent_manual, not_mappable_manual = read_atom_mappings("atom_mappings_manual.txt")

# these colors are used for diagnostics
# colors for corresponding atoms and half of their bonds
color_list = ['red','cyan','orange','blue','pink','wheat','gold','green','brown','purple','lightgrey','lime','lightblue','magenta','teal']
color_list = color_list + color_list + color_list + color_list + color_list + color_list + color_list  # never run out of colors

ribose = ["C2'","C3'","O3'","C4'","O4'","C5'"]  # for DNA
phosphate = ["O5'","P","OP1","OP2"]

# note:  C1' is not listed in ribose_full so it needs to be added when needed
ribose_full = ["C2'","C3'","O3'","C4'","O2'","O4'","C5'","H1'","H2'","H2''","H3'","H4'","H5'","H5''","HO3'","HO2'"]
phosphate_full = ["O5'","P","OP1","OP2","OP3","HOP1","HOP2","HOP3"]

# define a color scheme for parent nucleotides, which will transfer over to mapped atoms in modified nucleotides
parent_atom_to_color = {}
if color_scheme == 'diagnostic':
    for parent in NAbaseheavyatoms.keys():
        parent_atom_to_color[parent] = {}
        parent_atom_to_color[parent]["C1'"] = "tan"
        c = 0
        #for a in (NAbaseheavyatoms[parent] + NAbasehydrogens[parent] + ["C1'"]):
        for a in (NAbaseheavyatoms[parent] + NAbasehydrogens[parent]):
            parent_atom_to_color[parent][a] = color_list[c]
            c += 1

        # color backbone non-hydrogen atoms
        c = 3
        for a in ribose_full + phosphate_full:
            parent_atom_to_color[parent][a] = color_list[c]
            c += 1
        if parent in ['A','C','G','U']:
            parent_atom_to_color[parent]["O2'"] = "red"
        parent_atom_to_color[parent]["O5'"] = "magenta"
        parent_atom_to_color[parent]["P"] = "orange"

        print('Coloring for %s' % parent)
        print(parent_atom_to_color[parent])
elif color_scheme == 'CPK':
    # color atoms by CPK coloring
    for parent in NAbaseheavyatoms.keys():    # A, C, G, U, DA, DC, DG, DT
        parent_atom_to_color[parent] = {}
        for a in NAbaseheavyatoms[parent] + NAbasehydrogens[parent] + ["C1'"] + ribose_full + phosphate_full:
            # a is C7, N1, O2, etc.
            element = a[0]  # always one character for standard nucleic acids
            parent_atom_to_color[parent][a] = element_to_cpk_color(element)
else:
    print('Unknown color scheme %s' % color_scheme)

# load parent nucleotide cif data
par_coordinates = {}
par_connections = {}
par_atom_to_element = {}
par_comp_id = {}
par_base_coordinates = {}
par_atom_to_chirality = {}

for parent in standard_nts:
    # read the .cif files to get the coordinates and the atom to atom connections
    par_coord, par_conn, par_atom_to_elem, par_atom_to_chiral, par_comp_id, one_letter_code = get_cif_data(parent)

    par_coordinates[parent] = par_coord
    par_connections[parent] = par_conn
    par_atom_to_element[parent] = par_atom_to_elem
    par_atom_to_chirality[parent] = par_atom_to_chiral

    # remove hydrogens that were recorded in the QM calculations
    par_base_coordinates[parent] = NAbasecoordinates[parent]  # base coordinates in standard orientation
    if parent in ['A','G','DA','DG']:
        if "H9" in par_base_coordinates[parent]:
            del par_base_coordinates[parent]["H9"]
        if "H9'" in par_base_coordinates[parent]:
            del par_base_coordinates[parent]["H9'"]
    elif parent in ['C','U','DC']:
        if "H1" in par_base_coordinates[parent]:
            del par_base_coordinates[parent]["H1"]
        if "H1'" in par_base_coordinates[parent]:
            del par_base_coordinates[parent]["H1'"]


# get standard coordinates of parent nucleotides
par_coordinates_standard = {}

for parent in standard_nts:
    par_coordinates_standard[parent] = {}

    par_base_atoms = []
    par_atoms = []

    # assemble coordinates of mapped heavy base atoms; skip hydrogens
    for a in par_base_coordinates[parent].keys():
        if not "H" in a:
            par_base_atoms.append(par_base_coordinates[parent][a])
            par_atoms.append(par_coordinates[parent][a])

    # superimpose parent and modified base atoms
    U, new1, mean1, rmsd, sse, mean2 = besttransformation(par_base_atoms,par_atoms)

    print("%2s to %2s RMSD is %8.4f with heavy atoms" % (parent,parent,rmsd))

    # map each atom in the modified nucleotide into standard position
    for atom in par_coordinates[parent].keys():
        c = mean1 + np.dot(U,np.array(par_coordinates[parent][atom]) - mean2)
        par_coordinates_standard[parent][atom] = [c[0,0],c[0,1],c[0,2]]

# make a list of modified nucleotides that don't have mappings
not_mapped = []

# gather data about O2' mappings for RNA
if False:
    RNA_modified_url = []
    RNA_temp_mapping = []

    atom_list = ["O2'","P"]    # atoms of particular interest

    print('Modified RNA nucleotides missing one of these atoms: %s' % atom_list)
    for modified, parent in modified_base_to_parent.items():
        if parent in ['A','C','G','U']:
            RNA_modified_url.append('https://www.rcsb.org/ligand/%s' % modified)
            for atom in atom_list:
                if not modified_to_parent_atom[modified].get(atom,""):
                    print('https://www.rcsb.org/ligand/%s' % modified)
                    #RNA_temp_mapping.append('%s\t%s\t%s\t%s' % (parent,atom,modified,modified_to_parent_atom[modified].get(atom,"")))
                    map = '%s %s %s' % (parent,atom,modified)
                    RNA_temp_mapping.append(map)
                    print(map)


# keep track of changes
modified_to_changes = {}
modified_list = []

# read the list of modified nucleotides and their counts
with open('modified_nt_list.csv',read_mode) as f:
    lines = f.readlines()

# loop over modified nucleotides from most common to least
# plot the modified base
for line in lines:

    local_show_figure = show_figure

    fields = line.replace('"','').rstrip("\n").split(",")
    if len(fields) == 3:
        modified = fields[1]
        count = int(fields[2])

        if len(focus_list) > 0 and not modified in focus_list:
            continue

        if count > max_count:
            continue

        modified_list.append(modified)

        modified_to_changes[modified] = {}
        modified_to_changes[modified]['standard_base'] = []  # empty list when no standard base
        modified_to_changes[modified]['changes'] = []  # empty list when no atom changes
        modified_to_changes[modified]['count'] = count

        mod_coordinates, mod_connections, mod_atom_to_element, mod_atom_to_chirality, par_comp_id, one_letter_code = get_cif_data(modified)

        modified_to_changes[modified]['atom_count'] = len(mod_atom_to_element)

        if modified in not_mappable_manual:
            print('Skipping modified nucleotide %5s with count %4d because it is not mappable' % (modified,count))
            not_mapped.append('%5s with count %4d is not mappable' % (modified,count))
            modified_to_changes[modified]['error'] = 'not mappable'
            continue

        print("")
        print('Processing modified nucleotide %5s with count %4d' % (modified,count))

        if len(mod_coordinates) < 3:
            print('Modified nucleotide %s has only %s atoms' % (modified,len(mod_coordinates)))
            not_mapped.append('%5s with count %4d has only %d atoms' % (modified,count,len(mod_coordinates)))
            modified_to_changes[modified]['error'] = 'not enough matching atoms to map'
            continue

        # print('Found %d connections for %s' % (len(mod_connections)/2.0,modified))

        if modified_to_parent_atom_manual.get(modified,{}):
            print('Manual mappings for %s' % modified)
            print(modified_to_parent_atom_manual[modified])

        if modified in modified_base_to_parent:
            par_atoms = []
            mod_atoms = []

            par_backbone_atoms = []
            mod_backbone_atoms = []
            mod_backbone_atoms_reflected = []

            parent = modified_base_to_parent[modified]

            print('Parent nucleotide for %s is %s' % (modified,parent))

            modified_to_changes[modified]['standard_base'] = [parent]
            modified_to_changes[modified]['par_comp_id'] = par_comp_id
            modified_to_changes[modified]['one_letter_code'] = one_letter_code

            for a,b in parent_to_modified_atom[modified].items():
                # collect coordinates of mapped heavy base atoms
                if a in NAbaseheavyatoms[parent]:
                    if a in par_base_coordinates[parent] and len(par_coordinates[parent][a]) == 3:
                        if b in mod_coordinates and len(mod_coordinates[b]) == 3:
                            par_atoms.append(par_base_coordinates[parent][a])
                            mod_atoms.append(mod_coordinates[b])
                # collect coordinates of certain backbone atoms
                if a in ['P',"O5'","C5'","C2'","C3'","O3'","C4'","O4'"]:
                    if a in par_coordinates[parent] and len(par_coordinates[parent][a]) == 3:
                        if b in mod_coordinates and len(mod_coordinates[b]) == 3:
                            par_backbone_atoms.append(par_coordinates[parent][a])
                            mod_backbone_atoms.append(mod_coordinates[b])
                            x,y,z = mod_coordinates[b]
                            mod_backbone_atoms_reflected.append([-x,y,z])

            # if len(par_atoms) < 3:
            #     not_mapped.append('%5s with count %4d does not have enough atom mappings to superimpose' % (modified,count))
            #     modified_to_changes[modified]['error'] = 'not enough matching atoms to superimpose'
            #     continue

        else:
            not_mapped.append('%5s with count %4d does not have an atom mapping' % (modified,count))
            parent = 'Unknown'
            modified_to_changes[modified]['error'] = 'no known parent'
            continue

        # if len(par_atoms) < 3:
        #     # map modified base onto xy plane as well as possible, for plotting
        #     par_atoms = [[0,0,0],[0,1,0],[1,0,0]]   # points in xy plane
        #     mod_atoms = []

        #     for mod_atom in sorted(mod_coordinates.keys(), key = lambda x : len(x)):
        #         mod_atoms.append(mod_coordinates[mod_atom])
        #         if len(mod_atoms) == 3:
        #             break

        # map each atom in the modified nucleotide into standard position
        mod_coordinates_standard = {}
        if len(par_atoms) >= 3:
            # superimpose parent and modified base atoms
            U, new1, mean1, rmsd, sse, mean2 = besttransformation(par_atoms,mod_atoms)
            print("Base RMSD is %8.4f with heavy atoms" % rmsd)
            for atom in mod_coordinates.keys():
                c = mean1 + np.dot(U,np.array(mod_coordinates[atom]) - mean2)
                mod_coordinates_standard[atom] = [c[0,0],c[0,1],c[0,2]]

        else:
            # create enough information to make a graph anyway
            for atom in mod_coordinates.keys():
                mod_coordinates_standard[atom] = mod_coordinates[atom]

        if len(par_atoms) >= 3:
            # check mappings, try to improve hydrogen mappings on bases, record mappings
            for par_atom, par_coord in par_base_coordinates[parent].items():
                if par_atom in NAbaseheavyatoms[parent] or par_atom in NAbasehydrogens[parent]:
                    min_dist = 1
                    for mod_atom, mod_coord_standard in mod_coordinates_standard.items():
                        dist = my_norm(par_coord,mod_coord_standard)
                        if dist < min_dist:
                            min_dist = dist
                            nearest_mod_atom = mod_atom

                    # print("%s atom %4s is the closest to standard %s %4s distance %6.3f" % (modified,nearest_mod_atom,parent,par_atom,min_dist))

                    # find additional mappings or fix mappings for hydrogen atoms
                    if par_atom.startswith('H') and not nearest_mod_atom in modified_to_parent_atom[modified]:
                        if modified in parent_to_modified_atom_manual:
                            # this modified nucleotide has been mapped manually
                            if not par_atom in parent_to_modified_atom_manual[modified]:
                                # but this parent atom is not mapped manually
                                if not par_atom in parent_to_modified_atom[modified]:
                                    print("Mapping standard %-4s to modified %-4s distance is %8.2f" % (par_atom,nearest_mod_atom,min_dist))
                                elif not parent_to_modified_atom[modified][par_atom] == nearest_mod_atom:
                                    print("Re-mapping standard %-4s to modified %-4s distance is %8.2f" % (par_atom,nearest_mod_atom,min_dist))
                                parent_to_modified_atom[modified][par_atom] = nearest_mod_atom
                                modified_to_parent_atom[modified][nearest_mod_atom] = par_atom
                        else:
                            if not par_atom in parent_to_modified_atom[modified]:
                                print("Mapping standard %-4s to modified %-4s distance is %8.2f" % (par_atom,nearest_mod_atom,min_dist))
                            elif not parent_to_modified_atom[modified][par_atom] == nearest_mod_atom:
                                print("Re-mapping standard %-4s to modified %-4s distance is %8.2f" % (par_atom,nearest_mod_atom,min_dist))
                            parent_to_modified_atom[modified][par_atom] = nearest_mod_atom
                            modified_to_parent_atom[modified][nearest_mod_atom] = par_atom

        # attempt to map some backbone hydrogens if not already done
        hydrogen_to_heavy = {}
        hydrogen_to_heavy["H5'"]  = ["C4'","C5'","O5'"]
        hydrogen_to_heavy["H5''"] = ["O5'","C5'","C4'"]
        hydrogen_to_heavy["H2'"]  = ["C1'","C2'","C3'"]
        hydrogen_to_heavy["H2''"] = ["C3'","C2'","C1'"]

        for a, atom_list in hydrogen_to_heavy.items():
            # avoid trying to map H2'' from RNA to an atom on the modified base
            if not a in par_atom_to_element[parent]:
                continue
            if not a in parent_to_modified_atom[modified]:
                b, d, p, q = get_mod_atom_closest_to(atom_list,parent_to_modified_atom[modified],mod_coordinates)
                if b:
                    parent_to_modified_atom[modified][a] = b
                    modified_to_parent_atom[modified][b] = a
                    print("Mapping standard %-4s to modified %-4s distance %8.2f" % (a,b,d))

        # note chirality changes
        reversal = True
        num_RS_changes = 0
        for a, b in parent_to_modified_atom[modified].items():
            parent_chirality = par_atom_to_chirality[parent].get(a,"")
            modified_chirality = mod_atom_to_chirality.get(b,"")
            if parent_chirality in ['R','S'] and modified_chirality in ['N',parent_chirality]:
                reversal = False
            if modified_chirality in ['R','S'] and parent_chirality in ['N',modified_chirality]:
                reversal = False
            if parent_chirality in ['R','S'] and modified_chirality in ['R','S'] and not parent_chirality == modified_chirality:
                num_RS_changes += 1
            if parent_chirality and modified_chirality and not parent_chirality == modified_chirality:
                change_dict = {}
                change_dict['parent_atom'] = a
                change_dict['parent_chirality'] = parent_chirality
                change_dict['modified_atom'] = b
                change_dict['modified_chirality'] = modified_chirality
                change_dict['change_type'] = 'chirality'
                if a in phosphate_full:
                    change_dict['change_location'] = 'phosphate'
                elif a in ribose_full or a == "C1'":
                    change_dict['change_location'] = 'ribose'
                else:
                    change_dict['change_location'] = 'base'
                modified_to_changes[modified]['changes'].append(change_dict)
                print('Chirality change for parent %-4s %s mapped to %-4s %s' % (a,parent_chirality,b,modified_chirality))

        if num_RS_changes == 0:
            reversal = False

        if reversal:
            print('All %d chiral centers are reversed' % num_RS_changes)
            change_dict = {}
            change_dict['change_type'] = 'chirality_reversal'
            change_dict['change_location'] = 'ribose'
            modified_to_changes[modified]['changes'].append(change_dict)

        # find standard atoms not mapped
        all_parent_atoms = set(par_atom_to_element[parent].keys())
        mapped_parent_atoms = set(parent_to_modified_atom[modified].keys())
        parent_atoms_not_mapped = all_parent_atoms - mapped_parent_atoms
        print('Standard atoms not mapped: %s' % sorted(parent_atoms_not_mapped))
        for a in sorted(parent_atoms_not_mapped):
            change_dict = {}
            change_dict['parent_atom'] = a
            change_dict['parent_element'] = par_atom_to_element[parent][a]
            change_dict['change_type'] = 'removal'
            if a in phosphate_full:
                change_dict['change_location'] = 'phosphate'
            elif a in ribose_full or a == "C1'":
                change_dict['change_location'] = 'ribose'
            else:
                change_dict['change_location'] = 'base'
            modified_to_changes[modified]['changes'].append(change_dict)

        # find changes in covalent bonds between mapped atoms
        par_connections_mapped_atoms = set()
        for a1,a2 in par_connections[parent]:
            if a1 in parent_to_modified_atom[modified] and a2 in parent_to_modified_atom[modified]:
                par_connections_mapped_atoms.add((a1,a2))
                par_connections_mapped_atoms.add((a2,a1))

        mod_connections_as_parent = set()
        for (b1,b2) in mod_connections:
            if b1 in modified_to_parent_atom[modified] and b2 in modified_to_parent_atom[modified]:
                a1 = modified_to_parent_atom[modified][b1]
                a2 = modified_to_parent_atom[modified][b2]
                mod_connections_as_parent.add((a1,a2))
                mod_connections_as_parent.add((a2,a1))

        connections_in_parent_not_in_modified = par_connections_mapped_atoms - mod_connections_as_parent
        connections_in_modified_not_in_parent = mod_connections_as_parent - par_connections_mapped_atoms

        # print(connections_in_parent_not_in_modified)
        # print(connections_in_modified_not_in_parent)

        for s, t in [(connections_in_parent_not_in_modified,'removed_bond'),(connections_in_modified_not_in_parent,'added_bond')]:
            for a1,a2 in s:
                if a1 < a2:   # only list in one direction
                    change_dict = {}
                    change_dict['change_type'] = t
                    change_dict['parent_atom_1'] = a1
                    change_dict['parent_atom_2'] = a2
                    change_dict['modified_atom_1'] = parent_to_modified_atom[modified].get(a1,"")
                    change_dict['modified_atom_2'] = parent_to_modified_atom[modified].get(a2,"")
                    if a1 in phosphate_full:
                        change_dict['change_location'] = 'phosphate'
                    elif a1 in ribose_full or a1 == "C1'":
                        change_dict['change_location'] = 'ribose'
                    else:
                        change_dict['change_location'] = 'base'
                    if a2 in phosphate_full:
                        change_dict['change_location_2'] = 'phosphate'
                    elif a2 in ribose_full or a2 == "C1'":
                        change_dict['change_location_2'] = 'ribose'
                    else:
                        change_dict['change_location_2'] = 'base'
                    modified_to_changes[modified]['changes'].append(change_dict)

                    # print(change_dict)
                    # local_show_figure = True

        # record what changed and get ready to plot
        for backbone in [True, False]:
            par_atoms = []
            mod_atoms = []
            par_atom_colors = {}  # tells the bond color, point color, how to display
            mod_atom_colors = {}  # tells the bond color, point color, how to display

            base_atoms = NAbaseheavyatoms[parent] + NAbasehydrogens[parent] + ["C1'"]

            # give default colors to parent atoms
            for a in parent_atom_to_color[parent].keys():
                if backbone or a in base_atoms:
                    par_atom_colors[a] = (parent_atom_to_color[parent][a],"black","full")
                else:
                    par_atom_colors[a] = (parent_atom_to_color[parent][a],"black","not")

            # color mapped atoms of the parent and the modified nucleotide
            for a,b in sorted(parent_to_modified_atom[modified].items()):
                display = "full"
                if not backbone and a in ribose_full + phosphate_full:
                    # set display to "not" for mapped atoms that are not part of the base
                    mod_atom_colors[b] = ("","black","not")
                    continue

                if b in mod_coordinates and len(mod_coordinates[b]) == 3:
                    if color_scheme == 'CPK':
                        element = mod_atom_to_element[b]        # look up the element
                        color = element_to_cpk_color(element)   # look up the color for modified atom
                        par_color = element_to_cpk_color(par_atom_to_element[parent][a])  # look up the color for parent atom
                    elif a in parent_atom_to_color[parent]:
                        color = parent_atom_to_color[parent][a]
                        par_color = parent_atom_to_color[parent][a]
                    else:
                        print('Do not know how to color %s with parent atom %s' % (b,a))
                        color = "black"
                        par_color = "black"

                    # check if atom on modified residue is the same element as standard
                    if par_atom_to_element[parent][a] == mod_atom_to_element[b]:
                        point_color = "black"
                    else:
                        point_color = "white"

                        # parent atom gets the correct dot color
                        c,d,e = par_atom_colors[a]
                        par_atom_colors[a] = (c,point_color,e)

                        if backbone:
                            # record that the element changes, but only to that once
                            change_dict = {}
                            change_dict['parent_atom'] = a
                            change_dict['parent_element'] = par_atom_to_element[parent][a]
                            change_dict['modified_atom'] = b
                            change_dict['modified_element'] = mod_atom_to_element[b]
                            change_dict['change_type'] = 'replacement'
                            if a in phosphate_full:
                                change_dict['change_location'] = 'phosphate'
                            elif a in ribose_full or a == "C1'":
                                change_dict['change_location'] = 'ribose'
                            else:
                                change_dict['change_location'] = 'base'
                            modified_to_changes[modified]['changes'].append(change_dict)

                    if backbone or a in base_atoms:
                        mod_atom_colors[b] = (color,point_color,"full")

                    if a in base_atoms:
                        # collect atoms for superposition of bases
                        par_atoms.append(par_base_coordinates[parent][a])
                        mod_atoms.append(mod_coordinates[b])

            if not backbone:
                # compute base RMSD if possible
                if modified in parent_to_modified_atom and len(par_atoms) >= 3:
                    if len(par_atoms) >= 3:
                        U, new1, mean1, rmsd, sse, mean2 = besttransformation(par_atoms,mod_atoms)
                        #print(U,new1,mean1,rmsd,sse,mean2)

                        print("RMSD is %8.4f with hydrogens" % rmsd)
                        modified_to_changes[modified]['base_rmsd'] = rmsd

            # identify atoms connected to the modified residue that are not mapped and so are added
            for atom1,atom2 in mod_connections:
                if atom1 in mod_atom_colors and mod_atom_colors[atom1][2] == 'full' and not atom2 in mod_atom_colors:
                    if atom2 in mod_coordinates and len(mod_coordinates[atom2]) == 3:

                        if color_scheme == 'CPK':
                            element = mod_atom_to_element[atom2]
                            bond_color = element_to_cpk_color(element)
                        else:
                            bond_color = "black"

                        if not backbone and modified_to_parent_atom[modified].get(atom1,"") == "C1'":
                            mod_atom_colors[atom2] = (bond_color,new_atom_point_color,"not")
                        else:
                            mod_atom_colors[atom2] = (bond_color,new_atom_point_color,"thin")

                        # print atom2 is new, note changes
                        change_dict = {}
                        change_dict['modified_atom'] = atom1
                        change_dict['modified_element'] = mod_atom_to_element[atom1]
                        change_dict['new_modified_atom'] = atom2
                        change_dict['new_modified_element'] = mod_atom_to_element[atom2]
                        change_dict['change_type'] = 'addition'

                        parent_atom = modified_to_parent_atom[modified].get(atom1,"")

                        if parent_atom in phosphate_full:
                            change_dict['change_location'] = 'phosphate'
                        elif parent_atom in ribose_full or parent_atom == "C1'":
                            change_dict['change_location'] = 'ribose'
                        else:
                            change_dict['change_location'] = 'base'
                        modified_to_changes[modified]['changes'].append(change_dict)

            # avoid duplicate changes
            modified_to_changes[modified] = keep_unique_changes(modified_to_changes[modified])

            if backbone:
                figure_save_file = os.path.join(output_directory,"backbone_plots",'backbone_%s_%s.png' % (parent,modified))
            else:
                figure_save_file = os.path.join(output_directory,"base_plots",'base_%s_%s.png' % (parent,modified))

            figure_save_file_gif = figure_save_file.replace(".png",".gif")

            if (draw_figures or local_show_figure) and (overwrite_figures \
                or (not save_as_gif and not os.path.exists(figure_save_file    )) \
                or (    save_as_gif and not os.path.exists(figure_save_file_gif))):

                shift = 0.3
                (xmin, xmax, ymin, ymax) = (0, 1, 0, 1)

                if plot_standard:
                    # fig = plt.figure(figsize=(15.0, 9.0))
                    fig = plt.figure(figsize=(10.0, 10.0))
                    # plot parent base atoms
                    ax = fig.add_subplot(2, 2, 1, projection='3d')

                    plt.subplots_adjust(wspace=-0.20,hspace=-0.20)

                    xmin, xmax, ymin, ymax = draw_base_coordinates(parent,par_coordinates_standard[parent],par_connections[parent],par_atom_colors,ax)
                    # print('parent limits   ',xmin,xmax,ymin,ymax)

                    ax.set_aspect('equal')

                    if local_show_figure:
                        ax.view_init(elev=90, azim=-90, roll=45)
                    else:
                        ax.view_init(elev=90, azim=-90)
                    ax.axis("off")

                    # expand a bit to include full atom dots, which are cropped when outside of the axis limits
                    # ax.set_xlim(xmin-shift,xmax+shift)
                    # ax.set_ylim(ymin-shift,ymax+shift)

                    ax = fig.add_subplot(2, 2, 3, projection='3d')

                    # print the changes
                    text_list = []

                    for change in modified_to_changes[modified]['changes']:
                        if not backbone and not change['change_location'] == 'base':
                            continue
                        if change['change_type'] == 'replacement':
                            text_list.append('%s replaced with %s on %s' % (change['parent_atom'],change['modified_atom'],change['change_location']))
                        if change['change_type'] == 'addition':
                            text_list.append('%s added to %s' % (change['new_modified_atom'],change['change_location']))
                        if change['change_type'] == 'removal':
                            text_list.append('%s removed from %s' % (change['parent_atom'],change['change_location']))
                        if change['change_type'] == 'chirality':
                            text_list.append('%s chirality %s changed to %s on %s' % (change['parent_atom'],change['parent_chirality'],change['modified_chirality'],change['change_location']))
                        if change['change_type'] == 'chirality_reversal':
                            text_list.append('All chiral centers reversed')
                        if change['change_type'] == 'added_bond':
                            text_list.append('Added bond between %s and %s, locations %s and %s' % (change['modified_atom_1'],change['modified_atom_2'],change['change_location'],change['change_location_2']))
                        if change['change_type'] == 'removed_bond':
                            text_list.append('Removed bond between %s and %s, locations %s and %s' % (change['modified_atom_1'],change['modified_atom_2'],change['change_location'],change['change_location_2']))
                    y = 0.9
                    delta_y = 0.1
                    tfs = 8
                    for t in text_list:
                        ax.text(0,y,0,str(t),fontsize=tfs)
                        y -= delta_y
                    # if backbone and "chirality" in modified_to_changes[modified]:
                    #     ax.text(0,y,0,"Chirality change",fontsize=tfs)
                    #     y -= delta_y
                    # if backbone and "chirality_reversal" in modified_to_changes[modified]:
                    #     ax.text(0,y,0,"Chirality reversal",fontsize=tfs)
                    #     y -= delta_y
                    # for a in ["H5'","H5''","H2'","H2''"]:
                    #     if a in parent_to_modified_atom[modified]:
                    #         b = parent_to_modified_atom[modified][a]
                    #         if not a == b:
                    #             ax.text(0,y,0,"Backbone %s mapped to %s" % (a,parent_to_modified_atom[modified][a]),fontsize=tfs)
                    #             y -= delta_y

                    ax.set_aspect('equal')
                    ax.axis("off")
                    ax.view_init(elev=90, azim=-90)

                    ax = fig.add_subplot(2, 2, 2, projection='3d')

                else:
                    fig = plt.figure(figsize=(5.0, 6.0))
                    ax = fig.add_subplot(1, 1, 1, projection='3d')

                # plot modified nucleotide atoms
                #ax.set_title("%s count: %d" % (modified,count))

                xmin2, xmax2, ymin2, ymax2 = draw_base_coordinates(modified,mod_coordinates_standard,mod_connections,mod_atom_colors,ax)

                ax.set_aspect('equal')
                ax.axis("off")


                plt.tight_layout()
                # print('modified limits ',xmin2,xmax2,ymin2,ymax2)

                # use the axis limits from the parent, or wider if necessary
                # ax.set_xlim(min(xmin,xmin2)-shift,max(xmax,xmax2)+shift)
                # ax.set_ylim(min(ymin,ymin2)-shift,max(ymax,ymax2)+shift)

                if local_show_figure:
                    ax.view_init(elev=90, azim=-90, roll=45)
                else:
                    ax.view_init(elev=90, azim=-90)

                # plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, hspace=-0.20)
                plt.subplots_adjust(wspace=-0.20,hspace=-0.20)

                if local_show_figure:
                    plt.show()
                else:
                    plt.savefig(figure_save_file)
                    plt.close()

                if (save_as_gif or crop_out_white_space) and not local_show_figure:
                    image = imageio.imread(figure_save_file)

                    # Convert to PIL image
                    pil_image = Image.fromarray(image)

                    # Crop the image to remove white space
                    bbox = pil_image.getbbox()
                    bbox = my_bbox(pil_image)
                    # print("==============================================================")
                    # print(bbox)
                    cropped_image = pil_image.crop(bbox)

                    # Convert back to numpy array
                    cropped_array = np.array(cropped_image)

                    if save_as_gif:
                        # Save the cropped image as GIF
                        imageio.imwrite(figure_save_file.replace(".png",".gif"), cropped_array)
                    else:
                        # Save the cropped image as PNG
                        imageio.imwrite(figure_save_file, cropped_array)

if len(focus_list) == 0:
    # write out the mappings again, from most common modified nucleotide to least
    modified_written = set()
    with open('atom_mappings.txt',write_mode) as f:
        for line in lines:
            fields = line.replace('"','').rstrip("\n").split(",")
            if len(fields) == 3:
                count = int(fields[2])
                if count > 0:
                    modified = fields[1]

                    if not modified in modified_written and modified in modified_base_to_parent:
                        parent = modified_base_to_parent[modified]
                        for par_atom, mod_atom in parent_to_modified_atom[modified].items():
                            f.write('%s\t%s\t%s\t%s\n' % (parent,par_atom,modified,mod_atom))

                        if "chirality" in modified_to_changes[modified]:
                            print('Standard %-4s modified %-4s chirality reversed' % (parent,modified))

                    modified_written.add(modified)



    print('')
    print('%d messages about the mappings:' % len(not_mapped))
    print("\n".join(not_mapped))

    # list the modified nucleotides that are not mappable
    # print('%d modified nucleotides are not mappable:' % len(not_mappable))
    # print(",".join(not_mappable))

    # remove standard nucleotides from the list of changes
    for parent in standard_nts:
        if parent in modified_to_changes:
            del modified_to_changes[parent]

    changes_file = os.path.join(output_directory,'modified_to_changes.json')
    with open(changes_file, write_mode) as f:
        # write modified_to_changes to a file in json format
        f.write(json.dumps(modified_to_changes))

    print('Wrote changes to %s' % changes_file)

if draw_figures and save_as_gif:
    import glob

    file_pattern = os.path.join(output_directory,'base_plots', '*.png')
    png_files = glob.glob(file_pattern)
    for file in png_files:
        os.remove(file)

    file_pattern = os.path.join(output_directory,'backbone_plots', '*.png')
    png_files = glob.glob(file_pattern)
    for file in png_files:
        os.remove(file)
