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

# user settings below

color_scheme = 'diagnostic'  # use many colors, to check the atom mappings
color_scheme = 'CPK'         # use CPK coloring

if color_scheme == 'diagnostic':
    plot_standard = True     # include the standard base in the plots
    save_as_gif = False      # save as .png, which may be more robust
    crop_out_white_space = False

    show_figure = True       # pause to show each figure, enable rotation of coordinates
    show_figure = False

    output_directory = "diagnostic"

else:
    plot_standard = False    # just plot the modified nucleotide by itself
    save_as_gif = True       # .gif works a little better online
    crop_out_white_space = True   # read the image, crop, save again
    show_figure = False      # don't stop to show each modified nucleotide
    output_directory = ""

draw_figures = False     # don't draw new figures at all
draw_figures = True      # draw new figures if they don't already exist

overwrite_figures = True  # draw figures, overwriting existing ones.  Slow.
overwrite_figures = False # makes it easier to identify what is new

focus_list = ['MA6']  # list the ones you want to look at specifically, maybe with show_figure = True
focus_list = []       # process all modified nucleotides

# user settings above, program settings below

atom_label_color = 'black'
new_atom_point_color = '#80D1E3'   # blue of Argon since that probably won't be added

# from definitions import NAconnections
from fr3d.definitions import NAbasecoordinates
from fr3d.definitions import NAbaseheavyatoms
from fr3d.definitions import NAbasehydrogens
from fr3d.geometry.superpositions import besttransformation
from fr3d.modified.make_atom_mappings import read_monomer_cif
from fr3d.modified.make_atom_mappings import download_nakb_modified_nt_list

from collections import defaultdict
import imageio
import json
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os
import sys
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
        rgb = "#FFB5B5"
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
        rgb = '#A6A6AB'
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
        rgb = "#FFA100"
    elif el == "BR":
        rgb = '#A62929'
    elif el == "TE":
        rgb = "#D47A00"
    elif el == "I":
        rgb = '#940094'
    elif el == "PT":
        rgb = "#D0D0E0"

    if not rgb:
        print("Unknown element %s, need to know how to color it" % el)
        print(crashnow)

    return rgb


def get_cif_data_old(base):
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


def get_cif_data(base):
    """
    Read the .cif file and organize its data into data structures for this program
    """

    cif_data = read_monomer_cif(base)

    data = {}
    data['standard_base'] = []  # empty list when no standard base
    data['changes'] = []  # empty list when no atom changes

    data['pdb'] = {}
    data['pdb']['name'] = cif_data['chem_comp'][0].get('name',None)
    data['pdb']['chem_comp_type'] = cif_data['chem_comp'][0].get('type',None)
    data['pdb']['par_comp_id'] = cif_data['chem_comp'][0].get('mon_nstd_parent_comp_id','No mon_nstd_parent_comp_id line')
    data['pdb']['one_letter_code'] = cif_data['chem_comp'][0].get('one_letter_code','No one_letter_code line')
    data['pdb']['pdbx_initial_date'] = cif_data['chem_comp'][0].get('pdbx_initial_date','1900-01-01')
    data['pdb']['formula'] = cif_data['chem_comp'][0].get('formula','No formula line')

    for descriptor in cif_data['pdbx_chem_comp_descriptor']:
        t = descriptor['type']
        if t == "SMILES_CANONICAL":
            if not t in data:
                data['pdb'][t] = {}
            program = descriptor['program']
            if "openeye" in program.lower():
                data['pdb'][t][program] = {}
                data['pdb'][t][program]['program_version'] = descriptor['program_version']
                data['pdb'][t][program]['descriptor'] = descriptor['descriptor']

    # for k,v in data.items():
    #     print(k,v)
    # input("Press enter to continue")

    coordinates = {}
    coordinates_ideal = {}
    connections = []
    atom_to_element = {}
    atom_to_chirality = {}

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

    data['atom_count'] = len(atom_to_element)

    return data, coordinates, connections, atom_to_element, atom_to_chirality


def my_norm(x,y):

    return math.sqrt((x[0]-y[0])**2+(x[1]-y[1])**2+(x[2]-y[2])**2)


def get_pdb_to_modomics_mapping():
    """
    Read pdb_to_modomics.txt and modomics_id_to_png.txt
    pdb_to_modomics.txt has lines like:
    PSU	1.0000	185	pY	confirmed

    modomics_id_to_png.txt has lines like:
    2	2_dijNdIi
    https://genesilico.pl/modomics/media/mod_images/2_dijNdIi.png

    Return a dictionary mapping PDB identifier to modomics information
    """

    modomics_id_to_png = {}
    with open("modomics_id_to_png.txt",read_mode) as f:
        lines = f.readlines()
        for line in lines:
            fields = line.rstrip("\n").split("\t")
            modomics_id = fields[0]
            filename = fields[1]
            modomics_id_to_png[modomics_id] = 'https://genesilico.pl/modomics/media/mod_images/%s.png' % filename

    pdb_to_modomics_data = {}
    with open("pdb_to_modomics.txt",read_mode) as f:
        lines = f.readlines()
        for line in lines:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5 and fields[4] == 'confirmed':
                modified = fields[0]
                modomics_id = fields[2]
                modomics_short_name = fields[3]
                d = {}
                d['id'] = modomics_id
                d['short_name'] = modomics_short_name
                d['png'] = modomics_id_to_png.get(modomics_id,'')
                pdb_to_modomics_data[modified] = d

    return pdb_to_modomics_data


def read_atom_mappings(filename):
    """
    Read a text file of three or four columns
    (parent,parent atom,modified,modified atom)
    When the fourth column is missing, it means not to map the parent atom.
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
                if len(fields) == 4:
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

    c = []
    for a in atom_list:
        if a in parent_to_modified_atom:
            b = parent_to_modified_atom[a]
            c.append(np.array(mod_coordinates[b]))
        else:
            return None, None, None, None, None, None

    # compute location of H5' from C4', C5', O5', for example
    p, q = pyramidal_hydrogens(c[0],c[1],c[2],bond_length)

    left_min_dist = 1
    left_min_atom = ""
    right_min_dist = 1
    right_min_atom = ""
    for b in mod_coordinates.keys():
        d = np.array(mod_coordinates[b])
        right_dist = np.linalg.norm(p-d)
        if right_dist < right_min_dist:
            right_min_dist = right_dist
            right_min_atom = b
        left_dist = np.linalg.norm(q-d)
        if left_dist < left_min_dist:
            left_min_dist = left_dist
            left_min_atom = b

    # if min_atom:
    #     print('Hydrogen bond length is %8.4f' % np.linalg.norm(c[1]-min_d))

    return left_min_atom, left_min_dist, p, right_min_atom, right_min_dist, q

    # mapped_to = set(parent_to_modified_atom.values())

    # min_dist = 1
    # min_atom = None
    # # min_d = None
    # for b in mod_coordinates.keys():
    #     if not b in mapped_to:
    #         d = np.array(mod_coordinates[b])
    #         dist = np.linalg.norm(p-d)
    #         if dist < min_dist:
    #             min_dist = dist
    #             min_atom = b
    #             # min_d = d

    # # if min_atom:
    # #     print('Hydrogen bond length is %8.4f' % np.linalg.norm(c[1]-min_d))

    # return min_atom, min_dist, p, q


def draw_base_coordinates(base_seq,coordinates,connections,atom_to_display,backbone,ax,limits=None,shift=(0,0,0)):
    """
    Connects atoms to draw one base
    ax is the current axis
    """

    # shift coordinates to be able to display two nucleotides side by side
    xs = shift[0]
    ys = shift[1]
    zs = shift[2]

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
        fs = 10
        ds = 8

    for atom1,atom2 in connections:

        if atom1 in coordinates:
            w = coordinates[atom1]
            p = np.array([w[0] + xs, w[1] + ys, w[2] + zs])
        else:
            continue

        if atom2 in coordinates:
            w = coordinates[atom2]
            q = np.array([w[0] + xs, w[1] + ys, w[2] + zs])
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
    view_order['over_under'] = 7
    view_order['over_under_reversal'] = 8


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


def main(mod_nt=""):
    #######################################################
    # main block starts here

    standard_nts = ['A','C','G','U','DA','DC','DG','DT']

    # read provisional atom to atom mappings
    parent_to_modified_atom, modified_to_parent_atom, modified_base_to_parent, not_mappable = read_atom_mappings("atom_mappings_provisional.txt")

    # read manual mappings
    parent_to_modified_atom_manual, modified_to_parent_atom_manual, modified_base_to_parent_manual, not_mappable_manual = read_atom_mappings("atom_mappings_manual.txt")

    # download modified nucleotide counts
    mod_to_count = download_nakb_modified_nt_list()

    # get confirmed mappings from PDB identifiers to Modomics
    pdb_to_modomics_data = get_pdb_to_modomics_mapping()

    # these colors are used when color_scheme = diagnostic
    # colors for corresponding atoms and half of their bonds
    color_list = ['red','cyan','orange','blue','pink','wheat','gold','green','brown','purple','lightgrey','lime','lightblue','magenta','teal']
    color_list = ['peru','violet','fuchsia','wheat','gold','purple','brown','lightgrey','magenta','darkgoldenrod','darkkhaki','darkorchid','sienna']
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

            # make sure every atom is colored; this handles the hydrogens
            c = 0
            for a in NAbaseheavyatoms[parent] + NAbasehydrogens[parent]:
                parent_atom_to_color[parent][a] = color_list[c]
                c += 1

            # restart so every backbone hydrogen is colored the same
            c = 0
            for a in ribose_full + phosphate_full:
                parent_atom_to_color[parent][a] = color_list[c]
                c += 1

            # override the default colors for heavy atoms
            # reddish colors for oxygens, blueish colors for nitrogens, green for carbon, various for hydrogens
            parent_atom_to_color[parent]["C1'"] = "tan"
            parent_atom_to_color[parent]["C2'"] = "green"
            parent_atom_to_color[parent]["C3'"] = "springgreen"
            parent_atom_to_color[parent]["C4'"] = "teal"
            parent_atom_to_color[parent]["C5'"] = "chartreuse"
            parent_atom_to_color[parent]["O2'"] = "red"   # doesn't hurt when it's not there
            parent_atom_to_color[parent]["O3'"] = "tomato"
            parent_atom_to_color[parent]["O4'"] = "crimson"
            parent_atom_to_color[parent]["O5'"] = "deeppink"
            parent_atom_to_color[parent]["OP1"] = "lightcoral"
            parent_atom_to_color[parent]["OP2"] = "violet"
            parent_atom_to_color[parent]["OP3"] = "pink"
            parent_atom_to_color[parent]["P"] = "tab:orange"

            if parent in ['A','G','DA','DG']:
                parent_atom_to_color[parent]["N1"] = "mediumblue"
                parent_atom_to_color[parent]["N2"] = "darkslateblue"
                parent_atom_to_color[parent]["N3"] = "royalblue"
                parent_atom_to_color[parent]["N6"] = "cyan"
                parent_atom_to_color[parent]["N7"] = "deepskyblue"
                parent_atom_to_color[parent]["N9"] = "cornflowerblue"
                parent_atom_to_color[parent]["C2"] = "mediumspringgreen"
                parent_atom_to_color[parent]["C4"] = "forestgreen"
                parent_atom_to_color[parent]["C5"] = "lawngreen"
                parent_atom_to_color[parent]["C6"] = "olivedrab"
                parent_atom_to_color[parent]["C8"] = "lightgreen"
                parent_atom_to_color[parent]["O6"] = "red"
            else:
                parent_atom_to_color[parent]["N1"] = "cornflowerblue"
                parent_atom_to_color[parent]["N3"] = "royalblue"
                parent_atom_to_color[parent]["N4"] = "deepskyblue"
                parent_atom_to_color[parent]["C2"] = "mediumspringgreen"
                parent_atom_to_color[parent]["C4"] = "forestgreen"
                parent_atom_to_color[parent]["C5"] = "palegreen"
                parent_atom_to_color[parent]["C6"] = "olivedrab"
                parent_atom_to_color[parent]["C7"] = "darkgreen"
                parent_atom_to_color[parent]["O2"] = "red"
                parent_atom_to_color[parent]["O4"] = "tab:red"

            # parent_atom_to_color[parent][""] = ""

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

    # collect global plotting min and max values

    parent_to_min_max = {}
    for parent in NAbaseheavyatoms.keys():
        parent_to_min_max[parent] = {}
        parent_to_min_max[parent]['xmin'] = 100
        parent_to_min_max[parent]['xmax'] = -100
        parent_to_min_max[parent]['ymin'] = 100
        parent_to_min_max[parent]['ymax'] = -100
        parent_to_min_max[parent]['xmin2'] = 100
        parent_to_min_max[parent]['xmax2'] = -100
        parent_to_min_max[parent]['ymin2'] = 100
        parent_to_min_max[parent]['ymin3'] = 100

    # load parent nucleotide cif data
    par_coordinates = {}
    par_connections = {}
    par_atom_to_element = {}
    par_comp_id = {}
    par_base_coordinates = {}
    par_atom_to_chirality = {}

    for parent in standard_nts:
        # read the .cif files to get the coordinates and the atom to atom connections
        par_coord, par_conn, par_atom_to_elem, par_atom_to_chiral, par_comp_id, one_letter_code = get_cif_data_old(parent)
        par_data, par_coord, par_conn, par_atom_to_elem, par_atom_to_chiral = get_cif_data(parent)

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

        # make mappings of standard nucleotide atoms to standard, for making plots
        modified_base_to_parent[parent] = parent
        modified_to_parent_atom[parent] = {}
        parent_to_modified_atom[parent] = {}
        for atom in par_atom_to_elem.keys():
            modified_to_parent_atom[parent][atom] = atom
            parent_to_modified_atom[parent][atom] = atom

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
    modified_to_atoms = {}
    all_mod_atom_to_element = {}

    # DNA that have O2'
    DNA_with_O2_prime_counter = 0

    # loop over modified nucleotides
    # include standard nucleotides to also make images for
    # include nucleotides that are not mappable, and nucleotides that are not from NAKB list
    full_list = ["A","C","G","U","DA","DC","DG","DT"] + sorted(set(modified_base_to_parent.keys()) | mod_to_count.keys())
    for modified in full_list:

        if len(mod_nt) > 0 and not modified == mod_nt:
            continue

        local_show_figure = show_figure

        # when you want to redraw specific modified nucleotides, list them here
        revised_set = set(["ORP","AAB","3DR","NRI","NR1","DV3","92F","48Z","61H","HOL","HOB","NSU","T0T","I","X4A","NSU","LHO","CFV","BMN","A1LXS","A1BBA","63T","RF5","PYY","NP3","MM7","FFD","DRP","DPY","DDX","D3","ASU","D33","6U0","61H","2DF","YRR","YA4","PYP","DXD","D3N","48Z","WC7","S8U","A2M","ATP","CPN","OWR","92F","SAY","K1F","DV3","G35","I","IMP","CPN","80S","TOQ","2MA"])
        revised_set = []
        if modified in revised_set:
            redraw_figure = True
        else:
            redraw_figure = False

        # another way to focus on specifid nucleotides, put them in focus_list
        if len(focus_list) > 0 and not modified in focus_list:
            continue

        modified_list.append(modified)

        # mod_coordinates, mod_connections, mod_atom_to_element, mod_atom_to_chirality, par_comp_id, one_letter_code = get_cif_data_old(modified)
        mod_to_changes, mod_coordinates, mod_connections, mod_atom_to_element, mod_atom_to_chirality = get_cif_data(modified)

        modified_to_changes[modified] = mod_to_changes
        modified_to_changes[modified]['standard_base'] = []  # empty list when no standard base
        modified_to_changes[modified]['changes'] = []  # empty list when no atom changes
        modified_to_changes[modified]['count'] = mod_to_count.get(modified,0)
        modified_to_changes[modified]['atom_count'] = len(mod_atom_to_element)
        modified_to_atoms[modified] = set(mod_atom_to_element.keys())

        all_mod_atom_to_element[modified] = mod_atom_to_element

        if modified in pdb_to_modomics_data:
            modified_to_changes[modified]['modomics'] = pdb_to_modomics_data[modified]

        if modified in not_mappable_manual:
            print('Skipping modified nucleotide %5s because it is not mapped' % (modified))
            not_mapped.append('%5s is not mapped' % (modified))
            modified_to_changes[modified]['error'] = 'not mappable'
            # input("Press Enter to continue")
            continue

        print("")
        print('Processing nucleotide %5s' % (modified))

        if len(mod_coordinates) < 3:
            print('Modified nucleotide %s has only %s atoms' % (modified,len(mod_coordinates)))
            not_mapped.append('%5s has only %d atoms' % (modified,len(mod_coordinates)))
            modified_to_changes[modified]['error'] = 'not mappable'
            continue

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
            # modified_to_changes[modified]['par_comp_id'] = par_comp_id
            # modified_to_changes[modified]['one_letter_code'] = one_letter_code

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

        else:
            not_mapped.append('%5s does not have an atom mapping' % (modified))
            parent = 'Unknown'
            modified_to_changes[modified]['error'] = 'not mappable'
            continue

        if len(par_atoms) >= 3:
            enough_atoms_to_trust_rotation = True
        else:
            enough_atoms_to_trust_rotation = False
            # add just enough atoms to get the rotation of the glycosidic bond correct
            # like C2' and O4' if they are mapped
            for a in ["C2'","O4'","C3'","C4'"]:
                if len(par_atoms) < 3:
                    if a in par_coordinates_standard[parent] and len(par_coordinates_standard[parent][a]) == 3:
                        if a in parent_to_modified_atom[modified]:
                            b = parent_to_modified_atom[modified][a]
                            if b in mod_coordinates and len(mod_coordinates[b]) == 3:
                                par_atoms.append(par_coordinates_standard[parent][a])
                                mod_atoms.append(mod_coordinates[b])
                                print('Added %s from parent and %s from modified' % (a,b))


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
            # center the few modified atoms at the origin
            total = np.zeros(3)
            for mod_atom in mod_atoms:
                total += np.array(mod_atom)
            print('mean',total / len(mod_coordinates.keys()))
            for atom in mod_coordinates.keys():
                mod_coordinates_standard[atom] = mod_coordinates[atom] - total / len(mod_atoms)
                print(atom,mod_coordinates_standard[atom])
            # input("Press Enter to continue...  Those are the average coordinates of the modified atoms.")

        if enough_atoms_to_trust_rotation:
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

        # determine whether bonds on certain ribose atoms are "over" or "under" the ring
        over_under_change_count = 0

        target_list = []
        target_list.append((["C2'","C3'","C4'"],"O3'","O","H3'","H"))
        target_list.append((["O4'","C4'","C3'"],"C5'","C","H4'","H"))
        if parent in ['A','G','DA','DG']:
            target_list.append((["C2'","C1'","O4'"],"N9","N","H1'","H"))
        else:
            target_list.append((["C2'","C1'","O4'"],"N1","N","H1'","H"))
        if parent in ['A','C','G','U'] or "O2'" in mod_atom_to_element.keys():
            target_list.append((["C1'","C2'","C3'"],"O2'","O","H2'","H"))

        # if "O2'" in mod_atom_to_element.keys() and not parent in ['A','C','G','U']:
        #     print("Parent %s modified %s has O2'" % (parent,modified))
        #     input("Press Enter to continue")

        for parent_atom_list,left_target,left_a,right_target,right_a in target_list:
            left_atom, left_distance, left_coordinate, right_atom, right_distance, right_coordinate \
            = get_mod_atom_closest_to(parent_atom_list,parent_to_modified_atom[modified],mod_coordinates,1.25)
            if left_atom:
                print("Modified atom closest to normal %s location is %-4s distance %8.2f" % (left_target,left_atom,left_distance))
            if right_atom:
                print("Modified atom closest to normal %s location is %-4s distance %8.2f" % (right_target,right_atom,right_distance))
            # if the left atom is no longer the original element, but the right one is, an over/under change happened
            if not mod_atom_to_element.get(left_atom,"") == left_a and mod_atom_to_element.get(right_atom,"") == left_a:
                # local_show_backbone_figure = True

                target_maps_to = parent_to_modified_atom[modified].get(left_target,"")

                print("Over/under change for parent %-4s mapped to modified %-4s" % (left_target,target_maps_to))

                over_under_change_count += 1
                change_dict = {}
                change_dict['change_type'] = 'over_under'
                change_dict['change_location'] = 'ribose'
                change_dict['parent_atom'] = left_target
                change_dict['modified_atom'] = target_maps_to
                modified_to_changes[modified]['changes'].append(change_dict)

            if left_target == "O2'" and not parent in ['A','C','G','U']:
                # DNA parent but modified has O2'
                DNA_with_O2_prime_counter += 1
                t = "Number %s\n" % DNA_with_O2_prime_counter
                t += "%s has parent %s but has O2'\n" % (modified,parent)
                t += "chem_comp.type is %s\n" % modified_to_changes[modified]['pdb']['chem_comp_type']
                t += "mon_nstd_parent_comp_id is %s\n" % modified_to_changes[modified]['pdb']['par_comp_id']
                t += "one_letter_code is %s\n" % modified_to_changes[modified]['pdb']['one_letter_code']
                if not mod_atom_to_element.get(left_atom,"") == left_a and mod_atom_to_element.get(right_atom,"") == left_a:
                    t += "Over/under change for parent %-4s mapped to modified %-4s\n" % (left_target,target_maps_to)
                t += "\n"
                with open("DNA_with_O2_prime.txt","a") as f:
                    f.write(t)

        if over_under_change_count == len(target_list):

            print('All %d target ribose atoms are over/under changes' % over_under_change_count)
            change_dict = {}
            change_dict['change_type'] = 'over_under_reversal'
            change_dict['change_location'] = 'ribose'
            modified_to_changes[modified]['changes'].append(change_dict)

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
                    left_atom, left_distance, left_coordinate, right_atom, right_distance, right_coordinate \
                    = get_mod_atom_closest_to(atom_list,parent_to_modified_atom[modified],mod_coordinates)
                    if left_atom and not left_atom in parent_to_modified_atom[modified].values():
                        parent_to_modified_atom[modified][a] = left_atom
                        modified_to_parent_atom[modified][left_atom] = a
                        print("Mapping standard %-4s to modified %-4s distance %8.2f" % (a,left_atom,left_distance))

        # note chirality changes from .cif file using IUPAC definition
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
                figure_save_file = os.path.join(output_directory,"img",'backbone_%s_%s.png' % (parent,modified))
            else:
                figure_save_file = os.path.join(output_directory,"base_plots",'base_%s_%s.png' % (parent,modified))
                figure_save_file = os.path.join(output_directory,"img",'base_%s_%s.png' % (parent,modified))

            figure_save_file_gif = figure_save_file.replace(".png",".gif")

            if (draw_figures or local_show_figure) and (overwrite_figures or redraw_figure \
                or (not save_as_gif and not os.path.exists(figure_save_file    )) \
                or (    save_as_gif and not os.path.exists(figure_save_file_gif))):

                shift = 0.3
                (xmin, xmax, ymin, ymax) = (0, 1, 0, 1)
                ymin3 = 0

                if plot_standard:
                    # (width,height)
                    # fig = plt.figure(figsize=(15.0, 9.0))
                    fig = plt.figure(figsize=(12.0, 8.0))
                    ax = fig.add_subplot(1, 1, 1, projection='3d')

                    # plot parent base atoms
                    # ax = fig.add_subplot(2, 2, 1, projection='3d')

                    # plt.subplots_adjust(wspace=-0.20,hspace=-0.20)
                    if parent in ['A','G','DA','DG']:
                        if backbone:
                            parent_shift = (-12,0,0)
                        else:
                            parent_shift = (-10,0,0)
                    else:
                        if backbone:
                            parent_shift = (-10,0,0)
                        else:
                            parent_shift = (-8,0,0)
                    xmin, xmax, ymin, ymax = draw_base_coordinates(parent,par_coordinates_standard[parent],par_connections[parent],par_atom_colors,backbone,ax,shift=parent_shift)

                    # expand a bit to include full atom dots, which are cropped when outside of the axis limits
                    # ax.set_xlim(xmin-shift,xmax+shift)
                    # ax.set_ylim(ymin-shift,ymax+shift)

                    # ax = fig.add_subplot(2, 2, 3, projection='3d')

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
                            text_list.append('Added bond between %s and %s,' % (change['modified_atom_1'],change['modified_atom_2']))
                            text_list.append('locations %s and %s' % (change['change_location'],change['change_location_2']))
                        if change['change_type'] == 'removed_bond':
                            text_list.append('Removed bond between %s and %s,' % (change['modified_atom_1'],change['modified_atom_2']))
                            text_list.append('locations %s and %s' % (change['change_location'],change['change_location_2']))
                        if change['change_type'] == 'over_under':
                            text_list.append('Over/under change for %s mapped to %s' % (change['parent_atom'],change['modified_atom']))
                        if change['change_type'] == 'over_under_reversal':
                            text_list.append('All over/under situations reversed')

                    if parent in ['A','G','DA','DG']:
                        x = 1.0 + parent_shift[0] / 2.0
                    else:
                        x = 0.5 + parent_shift[0] / 2.0

                    if backbone:
                        y = -8
                    else:
                        y = -4
                    if len(text_list) > 15:
                        if backbone:
                            delta_y = 0.4
                            tfs = 8
                        else:
                            delta_y = 0.3
                            tfs = 8
                    else:
                        if backbone:
                            delta_y = 0.6
                            tfs = 12
                        else:
                            delta_y = 0.5
                            tfs = 12
                    for t in text_list:
                        ax.text(x,y,0,str(t),fontsize=tfs)
                        y -= delta_y

                    ymin3 = y

                else:
                    fig = plt.figure(figsize=(5.0, 6.0))
                    ax = fig.add_subplot(1, 1, 1, projection='3d')

                # plot modified nucleotide atoms

                xmin2, xmax2, ymin2, ymax2 = draw_base_coordinates(modified,mod_coordinates_standard,mod_connections,mod_atom_colors,backbone,ax)

                parent_to_min_max[parent]['xmin']  = min(xmin,parent_to_min_max[parent].get('xmin',xmin))
                parent_to_min_max[parent]['xmax']  = max(xmax,parent_to_min_max[parent].get('xmax',xmax))
                parent_to_min_max[parent]['xmin2'] = min(xmin2,parent_to_min_max[parent].get('xmin2',xmin2))
                parent_to_min_max[parent]['ymin']  = min(xmin,parent_to_min_max[parent].get('ymin',ymin))
                parent_to_min_max[parent]['ymin2'] = min(ymin2,parent_to_min_max[parent].get('ymin2',ymin2))
                parent_to_min_max[parent]['ymax']  = max(ymax,parent_to_min_max[parent].get('ymax',ymax))
                parent_to_min_max[parent]['xmax2'] = max(xmax2,parent_to_min_max[parent].get('xmax2',xmax2))
                parent_to_min_max[parent]['ymin3'] = min(ymin3,parent_to_min_max[parent].get('ymin3',ymin3))

                # print('x min max',min(xmin,xmin2)-shift,max(3,xmax,xmax2)+shift)
                # print('y min max',min(ymin,ymin2,ymin3)-shift,max(ymax,ymax2)+shift)
                # input("Press Enter to continue...")

                if plot_standard:
                    # use standard limits to keep the bases in the same position
                    if parent in ['A','G','DA','DG']:
                        if backbone:
                            ax.set_xlim(-19.5,6.0)
                            ax.set_ylim(-15,2.5)
                        else:
                            ax.set_xlim(-14,4.5)
                            ax.set_ylim(-9,2.5)
                    else:
                        if backbone:
                            ax.set_xlim(-14,5)
                            ax.set_ylim(-15,3.5)
                        else:
                            ax.set_xlim(-11,5)
                            ax.set_ylim(-9,3.5)

                ax.axis("off")

                ax.set_aspect('equal')

                plt.tight_layout()
                # plt.subplots_adjust(left=0.01, right=0.99, top=0.99, bottom=0.01, hspace=-0.20)
                # plt.subplots_adjust(wspace=-0.20,hspace=-0.20)

                if local_show_figure:
                    ax.view_init(elev=90, azim=-90, roll=45)
                    plt.show()

                ax.view_init(elev=90, azim=-90)
                plt.savefig(figure_save_file, bbox_inches='tight', pad_inches=0, dpi=150)
                plt.close()

                if save_as_gif or crop_out_white_space:
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

    if len(focus_list) == 0 and len(mod_nt) == 0:
        count_NAKB_C1p_name_changed = 0
        count_NAKB_same_element_name_changed = 0
        count_NAKB_atom_mappings = 0
        # write out the mappings again, from most common modified nucleotide to least
        modified_written = set()
        with open('atom_mappings.txt',write_mode) as f:
            for modified in modified_base_to_parent.keys():
                if not modified in modified_written and modified in modified_base_to_parent:
                    parent = modified_base_to_parent[modified]
                    modified_atoms_written = set()
                    # write out all parent atoms whether or not they are mapped
                    for par_atom in par_atom_to_element[parent].keys():
                        if par_atom in parent_to_modified_atom[modified]:
                            mod_atom = parent_to_modified_atom[modified][par_atom]
                            modified_atoms_written.add(mod_atom)
                            if modified in mod_to_count:
                                count_NAKB_atom_mappings += 1
                                if par_atom == "C1'" and not mod_atom == "C1'":
                                    count_NAKB_C1p_name_changed += 1
                                if not par_atom_to_element[parent][par_atom] == all_mod_atom_to_element[modified][mod_atom]:
                                    count_NAKB_same_element_name_changed += 1
                        else:
                            mod_atom = ""
                        f.write('%s\t%s\t%s\t%s\n' % (parent,par_atom,modified,mod_atom))

                    # write out any remaining modified atoms that are not mapped
                    par_atom = ""
                    for mod_atom in sorted(modified_to_atoms[modified] - modified_atoms_written):
                        f.write('%s\t%s\t%s\t%s\n' % (parent,par_atom,modified,mod_atom))

                    if "chirality" in modified_to_changes[modified]:
                        print('Standard %-4s modified %-4s chirality reversed' % (parent,modified))

                modified_written.add(modified)

        # write a small file of modified to parent mappings
        modified_to_parent = sorted(modified_base_to_parent.items(), key=lambda x : (x[1],x[0]))
        with open('nt_mappings.txt',write_mode) as f:
            for modified, parent in modified_to_parent:
                f.write("%s\t%s\n" % (modified,parent))

        print('')
        print('%d messages about the mappings:' % len(not_mapped))
        print("\n".join(not_mapped))

        print(sorted(modified_to_changes.keys()))

        # remove standard nucleotides from the list of changes
        for parent in standard_nts:
            if parent in modified_to_changes:
                del modified_to_changes[parent]

        # keep only modified nucleotides identified by NAKB non-standard residue list
        # because the .json file is for the NAKB modified nucleotide site
        # add DI and DU because they are being added to the NAKB site soon
        for modified in list(modified_to_changes.keys()):
            if not modified in mod_to_count and not modified in ['DI','DU']:
                del modified_to_changes[modified]
                print('Deleted %s from modified_to_changes because NAKB count is zero' % modified)
        print('Now there are %d modified residues in the dataset' % len(modified_to_changes.keys()))

        changes_file = 'modified_to_change_data.json'
        with open(changes_file, write_mode) as f:
            # write modified_to_changes to a file in json format
            f.write(json.dumps(modified_to_changes))

        print('Wrote changes to %s' % changes_file)

    if draw_figures and save_as_gif:
        import glob

        file_pattern = os.path.join(output_directory,'base_plots', '*.png')
        file_pattern = os.path.join(output_directory,'img', '*.png')
        png_files = glob.glob(file_pattern)
        for file in png_files:
            os.remove(file)

        # file_pattern = os.path.join(output_directory,'backbone_plots', '*.png')
        # file_pattern = os.path.join(output_directory,'img', '*.png')
        # png_files = glob.glob(file_pattern)
        # for file in png_files:
        #     os.remove(file)

    # for parent in parent_to_min_max.keys():
    #     print('Parent %s has min max %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f' % (parent,\
    #         parent_to_min_max[parent]['xmin'],\
    #         parent_to_min_max[parent]['xmax2'],\
    #         parent_to_min_max[parent]['ymin'],\
    #         parent_to_min_max[parent]['ymin2'],\
    #         parent_to_min_max[parent]['ymin3'],\
    #         parent_to_min_max[parent]['ymax']))

    year_to_mods = defaultdict(list)
    for mod, changes in modified_to_changes.items():
        year = changes['pdb']['pdbx_initial_date'].split('-')[0]
        year_to_mods[year].append(mod+"_"+str(changes['count']))

    for year, mods in sorted(year_to_mods.items()):
        print("%s\t%s\t%s" % (year, len(mods), ",".join(mods)))

if __name__ == "__main__":
    # Get mod_nt from command line argument
    if len(sys.argv) == 2:
        mod_nt = sys.argv[1]
        result = main(mod_nt)
    else:
        result = main()