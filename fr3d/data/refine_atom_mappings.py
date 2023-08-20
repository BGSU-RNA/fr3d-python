"""
refine_atom_mappings.py reads atom_mappings_provisional.txt,
superimposes the modified base on the parent base,
matches each parent base atom to the nearest modified base atom,
writes out those mappings as atom_mappings.txt,
and if there is no .png image file, it writes out a visual representation of the mappings

Review the image files and if there is something wrong, make a special case in make_atom_mappings.py
Delete the image file, then run make_atom_mappings.py again and then refine_atom_mappings again
Iterate until you've got it right.
Hopefully you only need to look at the most recent image files.
Use https://www.rcsb.org/ligand/A but substitute the non-standard nucleotide where A is

Ideas for the next version:
    Color the base letters red if the atom changes
    Choose colors closer to CPK for the base atoms; need enough reds, blues, grays.  Get close-ish
    Color connection to C1' to acknowledge that it is a connection, but don't color what C1' connects to, too complicated to view
"""

draw_figures = True
draw_figures = False

backbone = False
backbone = True

from fr3d.definitions import NAconnections
from fr3d.definitions import NAbasecoordinates
from fr3d.definitions import NAbaseheavyatoms
from fr3d.definitions import NAbasehydrogens
from fr3d.geometry.superpositions import besttransformation

import math
import matplotlib.pyplot as plt
import numpy as np
import os
import sys

if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
    read_mode = 'rb'
    write_mode = 'w'
else:
    from urllib.request import urlretrieve as urlretrieve
    read_mode = 'rt'
    write_mode = 'wt'   # write as text


def draw_base_coordinates(base_seq,coordinates,connections,atom_colors,ax,zorder=1,limits=None):
    """
    Connects atoms to draw one base
    ax is the current axis
    zorder is for 2d graphs to tell the plotting order
    """

    if limits:
        xmin,xmax,ymin,ymax = limits
        ymax_shift = 0.0  # how far to go above ymax for modified nucleotides, because there is space there
        print("Limits:")
        print(xmin,xmax,ymin,ymax)
    else:
        xmin = 100
        xmax = -100
        ymin = 100
        ymax = -100
        ymax_shift = 0.0

    print('Found %d connections for %s' % (len(connections),base_seq))

    connections = sorted(connections)

    drawn_connections = set([])
    drawn_atoms = set([])

    if backbone:
        lw = 6.0
        fs = 8
        ds = 6
    else:
        lw = 12.0
        fs = 16
        ds = 16       # need to modify to get this right

    for atom1,atom2 in connections:

        if (atom1,atom2) in drawn_connections:
            continue

        drawn_connections.add((atom2,atom1))
        drawn_connections.add((atom1,atom2))

        # print(atom_colors.get(atom1,("","black")))
        # print(atom_colors.get(atom2,("","black")))

        c1,pc1 = atom_colors.get(atom1,("","black"))
        c2,pc2 = atom_colors.get(atom2,("","black"))

        try:
            p = coordinates[atom1]
            q = coordinates[atom2]

            if atom1 in atom_colors and not atom1 in drawn_atoms:
                ax.text(p[0],p[1],atom1,fontsize=fs,color='darkgray')
                ax.scatter([p[0]],[p[1]], s=ds, color=pc1, zorder=zorder+2)
                drawn_atoms.add(atom1)

            if atom1 in atom_colors and not c1 == "tan" and not atom2 in drawn_atoms:
                # plot atoms that are connected to a mapped atom
                ax.text(q[0],q[1],atom2,fontsize=fs,color='darkgray')
                ax.scatter([q[0]],[q[1]], s=ds, color=pc2, zorder=zorder+2)
                drawn_atoms.add(atom2)

            if atom2 in atom_colors and not atom2 in drawn_atoms:
                ax.text(q[0],q[1],atom2,fontsize=fs,color='darkgray')
                ax.scatter([q[0]],[q[1]], s=ds, color=pc2, zorder=zorder+2)
                drawn_atoms.add(atom2)

            if atom2 in atom_colors and not c2 == "tan" and not atom1 in drawn_atoms:
                # plot atoms that are connected to a mapped atom
                ax.text(p[0],p[1],atom1,fontsize=fs,color='darkgray')
                ax.scatter([p[0]],[p[1]], s=ds, color=pc1, zorder=zorder+2)

            if base_seq in ['A','C','G','U','DA','DC','DG','DT'] or (atom1 in atom_colors and atom2 in atom_colors):
                # count all standard base atoms toward the size of the frame, also color atoms on modified
                xmin = min(xmin,p[0],q[0])
                xmax = max(xmax,p[0],q[0])
                ymin = min(ymin,p[1],q[1])
                ymax = max(ymax,p[1],q[1])

        except:
            print('Could not plot %s,%s' % (atom1,atom2))
            continue

        if atom1 in atom_colors and atom2 in atom_colors:
            s = [(3*p[0]+q[0])/4.0,(3*p[1]+q[1])/4.0,(3*p[2]+q[2])/4.0] # quarter point, near atom1
            m = [(p[0]+q[0])/2.0,(p[1]+q[1])/2.0,(p[2]+q[2])/2.0] # midpoint
            t = [(p[0]+3*q[0])/4.0,(p[1]+3*q[1])/4.0,(p[2]+3*q[2])/4.0] # three-quarters point, near atom2

            # points in order from atom1 to atom2 are p, s, m, t, q

            """
            # CYK coloring
            atom_colors = {}
            atom_colors['C'] = [0.6,0.6,0.6]
            atom_colors['N'] = [0,0,1]
            atom_colors['O'] = [1,0,0]
            """

            try:
                if c1 == c2:
                    # draw from atom1 to atom2 with project to make lines cover the atoms
                    ax.plot([p[0],q[0]],[p[1],q[1]], color=c1, linewidth=lw, zorder=zorder)
                else:
                    # draw from atom1 to quarter way across with projection to make lines cover the atom
                    ax.plot([p[0],s[0]],[p[1],s[1]], color=c1, linewidth=lw, zorder=zorder, solid_capstyle='round')
                    # draw from quarter way to midpoint with ends that don't overlap at the midpoint
                    ax.plot([s[0],m[0]],[s[1],m[1]], color=c1, linewidth=lw, zorder=zorder, solid_capstyle='butt')
                    # draw from midpoint to three-quarters with ends that don't overlap at the midpoint
                    ax.plot([m[0],t[0]],[m[1],t[1]], color=c2, linewidth=lw, zorder=zorder, solid_capstyle='butt')
                    # draw from three-quarters to atom2 with projection to make lines cover the atoms
                    ax.plot([t[0],q[0]],[t[1],q[1]], color=c2, linewidth=lw, zorder=zorder, solid_capstyle='round')

                xmin = min(xmin,p[0],q[0])
                xmax = max(xmax,p[0],q[0])
                ymin = min(ymin,p[1],q[1])
                ymax = max(ymax,p[1],q[1])

            except:
                continue

        elif (atom1 in atom_colors or atom2 in atom_colors) and not c1 == "tan" and not c2 == "tan":
            # C1' is tan, don't extend into the backbone atoms; fragile
            # connect to unmapped atoms that are connected to mapped atoms
            ax.plot([p[0],q[0]],[p[1],q[1]], color='gray', linewidth=2.0, zorder=zorder+1)

            xmin = min(xmin,p[0],q[0])
            xmax = max(xmax,p[0],q[0])
            ymin = min(ymin,p[1],q[1])
            ymax = max(ymax,p[1],q[1])

    return xmin, xmax, ymin, ymax


def read_monomer_coordinates(base):
    """
    Open and read the .cif file as a text file
    """

    coordinates = {}
    connections = []

    if base == 'PRN':
        filename = 'data_PRN.cif'        # crazy windows thing, can't name a file PRN so use data_PRN
    else:
        filename = base+'.cif'

    with open(os.path.join('cif',filename),read_mode) as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith(base):
            fields = line.split()
            if base in ['LHH','AFG','RFJ'] and len(fields) == 18:
                try:
                    atom = fields[1].replace('"','')
                    coordinates[atom] = [float(fields[12]),float(fields[13]),float(fields[14])]
                except:
                    continue
            elif len(fields) == 18:
                try:
                    atom = fields[1].replace('"','')
                    coordinates[atom] = [float(fields[9]),float(fields[10]),float(fields[11])]
                except:
                    continue

            if len(fields) == 7:
                atom1 = fields[1].replace('"','')
                atom2 = fields[2].replace('"','')
                connections.append((atom1,atom2))

            if base in ['LHH']:
                print(fields)

    return coordinates, connections


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
        if len(fields) >= 4:
            parent = fields[0]
            modified = fields[2]

            if not parent:
                not_mappable.append(modified)
                continue

            if not modified in modified_to_parent_atom:
                modified_to_parent_atom[modified] = {}
                parent_to_modified_atom[modified] = {}
                modified_base_to_parent[modified] = parent

            parent_atom = fields[1]
            modified_atom = fields[3]

            modified_to_parent_atom[modified][modified_atom] = parent_atom
            parent_to_modified_atom[modified][parent_atom] = modified_atom

    return parent_to_modified_atom, modified_to_parent_atom, modified_base_to_parent, not_mappable

#######################################################
# main block starts here

# read provisional atom to atom mappings
parent_to_modified_atom, modified_to_parent_atom, modified_base_to_parent, not_mappable = read_atom_mappings("atom_mappings_provisional.txt")

# read manual mappings
parent_to_modified_atom_manual, modified_to_parent_atom_manual, modified_base_to_parent_manual, not_mappable_manual = read_atom_mappings("atom_mappings_provisional.txt")

# colors for corresponding atoms and half of their bonds
color_list = ['red','cyan','orange','blue','pink','wheat','gold','green','brown','purple','lightgrey','lime','lightblue','magenta','teal']
color_list = color_list + color_list + color_list + color_list + color_list + color_list + color_list  # never run out of colors

ribose = ["C2'","C3'","O3'","C4'","O4'","C5'"]  # for DNA
phosphate = ["O5'","P","OP1","OP2"]

parent_atom_to_color = {}
for parent in NAbaseheavyatoms.keys():
    parent_atom_to_color[parent] = {}
    parent_atom_to_color[parent]["C1'"] = "tan"
    c = 0
    #for a in (NAbaseheavyatoms[parent] + NAbasehydrogens[parent] + ["C1'"]):
    for a in (NAbaseheavyatoms[parent] + NAbasehydrogens[parent]):
        parent_atom_to_color[parent][a] = color_list[c]
        c += 1

    if backbone:
        c = 3
        for a in ribose + phosphate:
            parent_atom_to_color[parent][a] = color_list[c]
            c += 1
        if parent in ['A','C','G','U']:
            parent_atom_to_color[parent]["O2'"] = "red"
        parent_atom_to_color[parent]["O5'"] = "magenta"
        parent_atom_to_color[parent]["P"] = "orange"

    print(parent_atom_to_color[parent])

not_mapped = []

# gather data about O2' mappings for RNA
if True:
    RNA_modified_url = []
    RNA_temp_mapping = []

    atom_list = ["O2'","P"]    # atoms of particular interest

    for modified, parent in modified_base_to_parent.items():
        if parent in ['A','C','G','U']:
            RNA_modified_url.append('https://www.rcsb.org/ligand/%s' % modified)
            for atom in atom_list:
                if not modified_to_parent_atom[modified].get(atom,""):
                    print('https://www.rcsb.org/ligand/%s' % modified)
                    #RNA_temp_mapping.append('%s\t%s\t%s\t%s' % (parent,atom,modified,modified_to_parent_atom[modified].get(atom,"")))
                    RNA_temp_mapping.append('%s %s %s' % (parent,atom,modified))

    # for url in RNA_modified_url:
    #     print(url)

    for map in RNA_temp_mapping:
        print(map)


# read the list of modified nucleotides and their counts
with open('modified_nt_list.csv',read_mode) as f:
    lines = f.readlines()

# loop over modified nucleotides from most common to least
# plot the parent base and the modified base
for line in lines:
    fields = line.rstrip("\n").split(",")
    if len(fields) == 3:
        count = int(fields[2])
        if count > 0:
            modified = fields[1]

            if modified in not_mappable:
                print('Skipping modified nucleotide %3s with count %4d because it is not mappable' % (modified,count))
                continue

            print("")
            print('Processing modified nucleotide %3s with count %4d' % (modified,count))

            mod_coordinates, mod_connections = read_monomer_coordinates(modified)
            mod_coordinates_standard = {}

            if len(mod_coordinates) < 3:
                print('Modified nucleotide %s has only %s atoms' % (modified,len(mod_coordinates)))
                not_mapped.append('%4s with count %4d has only %d atoms' % (modified,count,len(mod_coordinates)))
                continue

            par_atoms = []
            mod_atoms = []

            par_coordinates = {}

            if modified in modified_base_to_parent:
                parent = modified_base_to_parent[modified]

                # read the .cif files to get the coordinates and the atom to atom connections
                par_coordinates, par_connections = read_monomer_coordinates(parent)
                par_coordinates = NAbasecoordinates[parent]  # base coordinates in standard orientation

                # assemble coordinates of mapped heavy base atoms; skip hydrogens
                for a,b in parent_to_modified_atom[modified].items():
                    if a in NAbaseheavyatoms[parent]:
                        if a in par_coordinates and len(par_coordinates[a]) == 3:
                            if b in mod_coordinates and len(mod_coordinates[b]) == 3:
                                par_atoms.append(par_coordinates[a])
                                mod_atoms.append(mod_coordinates[b])

                if len(par_atoms) < 3:
                    not_mapped.append('%4s with count %4d does not have enough atom mappings to superimpose' % (modified,count))

            else:
                not_mapped.append('%4s with count %4d does not have an atom mapping' % (modified,count))
                parent = 'Unknown'

            if len(par_atoms) >= 3:
                enough_atoms_to_map = True
            else:
                enough_atoms_to_map = False

                # map modified base onto xy plane as well as possible
                par_atoms = [[0,0,0],[0,1,0],[1,0,0]]   # points in xy plane

                for mod_atom in sorted(mod_coordinates.keys(), key = lambda x : len(x)):
                    mod_atoms.append(mod_coordinates[mod_atom])
                    if len(mod_atoms) == 3:
                        break

            if len(par_atoms) != len(mod_atoms):
                print("Warning: Different number of atoms)")
                print(par_atoms)
                print(mod_atoms)
                print(mod_coordinates)

            # superimpose base atoms, fix hydrogen mappings if necessary
            U, new1, mean1, rmsd, sse, mean2 = besttransformation(par_atoms,mod_atoms)

            print("RMSD is %8.4f with heavy atoms" % rmsd)

            # map each atom in the modified nucleotide into standard position
            for atom in mod_coordinates.keys():
                c = mean1 + np.dot(U,np.array(mod_coordinates[atom]) - mean2)
                mod_coordinates_standard[atom] = [c[0,0],c[0,1],c[0,2]]

            # check mappings, try to improve hydrogen mappings on bases, record mappings
            if enough_atoms_to_map:
                for par_atom, par_coord in par_coordinates.items():
                    if par_atom in NAbaseheavyatoms[parent] or par_atom in NAbasehydrogens[parent]:
                        min_dist = 9999
                        for mod_atom, mod_coord_standard in mod_coordinates_standard.items():
                            dist = my_norm(par_coord,mod_coord_standard)
                            if dist < min_dist:
                                min_dist = dist
                                nearest_mod_atom = mod_atom

                        print("%s atom %4s is closest to parent %4s distance %6.3f" % (modified,nearest_mod_atom,par_atom,min_dist))

                        # find additional mappings or fix mappings for hydrogen atoms
                        if par_atom.startswith('H') and not nearest_mod_atom in modified_to_parent_atom[modified]:
                            if not par_atom in parent_to_modified_atom_manual[modified]:

                                print(modified_to_parent_atom_manual[modified])

                                if not par_atom in parent_to_modified_atom[modified]:
                                    print("  Mapping standard %4s to modified %4s" % (par_atom,nearest_mod_atom))
                                elif not parent_to_modified_atom[modified][par_atom] == nearest_mod_atom:
                                    print("  Re-mapping standard %4s to modified %4s" % (par_atom,nearest_mod_atom))
                                parent_to_modified_atom[modified][par_atom] = nearest_mod_atom
                                modified_to_parent_atom[modified][nearest_mod_atom] = par_atom



            if backbone:
                figure_save_file = 'backbone_%s_%s.png' % (parent,modified)
                figure_save_file = os.path.join("backbone_plots",figure_save_file)
            elif enough_atoms_to_map:
                figure_save_file = '%s_%s.png' % (parent,modified)
                figure_save_file = os.path.join("plots",figure_save_file)
            else:
                continue

            if draw_figures and not os.path.exists(figure_save_file):
                par_atoms = []
                mod_atoms = []
                par_atom_colors = {}
                mod_atom_colors = {}
                par_atoms_colored = set([])
                par_atoms_not_colored = set([])

                # color corresponding atoms
                for a,b in sorted(parent_to_modified_atom[modified].items()):
                    if a in NAbaseheavyatoms[parent] or a in NAbasehydrogens[parent] or a == "C1'":
                    # if a in NAbaseheavyatoms[parent] or a in NAbasehydrogens[parent]:
                        if a in par_coordinates and len(par_coordinates[a]) == 3:
                            if b in mod_coordinates and len(mod_coordinates[b]) == 3:
                                mod_atoms.append(mod_coordinates[b])
                                color = parent_atom_to_color[parent][a]
                                if a[0] == b[0]:
                                    # same atoms
                                    point_color = "black"
                                else:
                                    # different atoms are mapped
                                    point_color = "white"
                                mod_atom_colors[b] = (color,point_color)

                                if not a in ["H1'","H9'"]:
                                    par_atoms.append(par_coordinates[a])
                                    par_atom_colors[a] = (color,point_color)
                                    par_atoms_colored.add(a)

                par_atoms_not_colored = set(par_coordinates.keys()) - par_atoms_colored - set(["C1'","H9","H9'","H1"])

                # computed RMSD if possible
                if modified in parent_to_modified_atom and enough_atoms_to_map:
                    if len(par_atoms) >= 3:
                        U, new1, mean1, rmsd, sse, mean2 = besttransformation(par_atoms,mod_atoms)
                        #print(U,new1,mean1,rmsd,sse,mean2)

                        print("RMSD is %8.4f with hydrogens" % rmsd)

                if backbone:
                    for a,b in sorted(parent_to_modified_atom[modified].items()):
                        if not a in NAbaseheavyatoms[parent] + NAbasehydrogens[parent] + ["C1'"]:
                            if b in mod_coordinates and len(mod_coordinates[b]) == 3:
                                mod_atoms.append(mod_coordinates[b])
                                if a in parent_atom_to_color[parent]:
                                    color = parent_atom_to_color[parent][a]
                                else:
                                    color = "gainsboro"
                                if a[0] == b[0]:
                                    # same atoms
                                    point_color = "black"
                                else:
                                    # different atoms are mapped
                                    point_color = "white"
                                mod_atom_colors[b] = (color,point_color)
                                par_atoms_colored.add(a)

                # else:
                #     # color atoms of parent and modified separately
                #     print("Coloring parent and modified atoms separately")
                #     if modified in modified_base_to_parent:
                #         c = 0
                #         for a in NAbaseheavyatoms[parent] + NAbasehydrogens[parent]:
                #             par_atom_colors[a] = (color_list[c],"black")
                #             c += 1

                #     c = 0
                #     for b in mod_coordinates:
                #         mod_atom_colors[b] = (color_list[c],"black")
                #         c = min(c+1,len(color_list)-1)

                # plot the two bases, one in each subplot
                fig = plt.figure(figsize=(10.0, 6.0))

                # plot parent base atoms
                ax = fig.add_subplot(1, 2, 1)
                ax.set_aspect('equal',adjustable='box')
                #ax.axis("off")
                ax.axis("off")
                ax.set_title("%s RMSD %0.2f" % (parent,rmsd))

                xmin, xmax, ymin, ymax = draw_base_coordinates(parent,par_coordinates,par_connections,par_atom_colors,ax,1)

                # expand a bit to include full atom dots, which are cropped when outside of the axis limits
                shift = 0.2
                ax.set_xlim(xmin-shift,xmax+shift)
                ax.set_ylim(ymin-shift,ymax+shift)

                # plot modified nucleotide base atoms
                ax = fig.add_subplot(1, 2, 2)
                ax.set_aspect('equal',adjustable='box')
                ax.axis("off")
                ax.set_title("%s count: %d" % (modified,count))

                xmin2, xmax2, ymin2, ymax2 = draw_base_coordinates(modified,mod_coordinates_standard,mod_connections,mod_atom_colors,ax,1,(xmin,xmax,ymin,ymax))

                print('parent limits   ',xmin,xmax,ymin,ymax)
                print('modified limits ',xmin2,xmax2,ymin2,ymax2)

                if not parent == 'Unknown':
                    # use the axis limits from the parent, or wider if necessary
                    ax.set_xlim(min(xmin,xmin2)-shift,max(xmax,xmax2)+shift)
                    ax.set_ylim(min(ymin,ymin2)-shift,max(ymax,ymax2)+shift)

                plt.savefig(figure_save_file)
                plt.close()

# write out the mappings again, from most common modified nucleotide to least
modified_written = set([])

with open('atom_mappings.txt',write_mode) as f:

    for line in lines:
        fields = line.rstrip("\n").split(",")
        if len(fields) == 3:
            count = int(fields[2])
            if count > 0:
                modified = fields[1]

                if not modified in modified_written and modified in modified_base_to_parent:
                    parent = modified_base_to_parent[modified]
                    for par_atom, mod_atom in parent_to_modified_atom[modified].items():
                        f.write('%s\t%s\t%s\t%s\n' % (parent,par_atom,modified,mod_atom))

                modified_written.add(modified)

print('%d messages about the mappings:' % len(not_mapped))
print("\n".join(not_mapped))

print('%d modified nucleotides are not mappable:' % len(not_mappable))
print(",".join(sorted(not_mappable)))
