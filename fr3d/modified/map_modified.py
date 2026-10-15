"""
map_modified.py - Interactive 3D visualization tool for mapping modified nucleotides
"""

import sys
import os
import tkinter as tk
from tkinter import ttk
import matplotlib
matplotlib.use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import math
import re
import pdbx
from urllib.request import urlretrieve

from refine_atom_mappings import read_atom_mappings
from fr3d.definitions import NAbaseheavyatoms, NAbasehydrogens, NAbasecoordinates
from fr3d.geometry.superpositions import besttransformation


def element_to_cpk_color(element):
    # Copied from refine_atom_mappings.py and kept local so that file is independent
    rgb = None

    el = element.upper()

    if el == "H":
        rgb = [0.8,0.8,0.8]  # light gray because we have a white background
    elif el == "B":
        rgb = "#FFB5B5"
    elif el == "C":
        rgb = '#808080'  # darker gray
    elif el == "N":
        rgb = '#3050F8'  # blue
    elif el == "O":
        rgb = '#FF0D0D'  # red
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
        rgb = '#808080'

    return rgb


# Build parent_atom_to_color mapping (copied/adapted from refine_atom_mappings.py)
color_scheme = 'diagnostic'

# colors for corresponding atoms and half of their bonds
color_list = ['peru','violet','fuchsia','wheat','gold','purple','brown','lightgrey','magenta','darkgoldenrod','darkkhaki','darkorchid','sienna']
color_list = color_list + color_list + color_list + color_list + color_list + color_list + color_list

ribose = ["C2'","C3'","O3'","C4'","O4'","C5'"]
phosphate = ["O5'","P","OP1","OP2"]
ribose_full = ["C2'","C3'","O3'","C4'","O2'","O4'","C5'","H1'","H2'","H2''","H3'","H4'","H5'","H5''","HO3'","HO2'"]
phosphate_full = ["O5'","P","OP1","OP2","OP3","HOP1","HOP2","HOP3"]

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

        # override defaults for heavy/backbone atoms
        parent_atom_to_color[parent]["C1'"] = "tan"
        parent_atom_to_color[parent]["C2'"] = "green"
        parent_atom_to_color[parent]["C3'"] = "springgreen"
        parent_atom_to_color[parent]["C4'"] = "teal"
        parent_atom_to_color[parent]["C5'"] = "chartreuse"
        parent_atom_to_color[parent]["O2'"] = "red"
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
else:
    # fallback: CPK coloring for any parent
    for parent in NAbaseheavyatoms.keys():
        parent_atom_to_color[parent] = {}
        for a in NAbaseheavyatoms[parent] + NAbasehydrogens[parent] + ["C1'"] + ribose_full + phosphate_full:
            element = a[0]
            parent_atom_to_color[parent][a] = element_to_cpk_color(element)


def read_component_from_file_or_download(mod_nt, cif_dir="cif_files"):
    """
    Read one CCD component from cif/<ID>.cif, downloading it from RCSB if missing.
    """
    comp_id = mod_nt.upper()
    if comp_id == "PRN":
        # Windows reserves PRN as a device name
        filename = os.path.join(cif_dir, f"data_{comp_id}.cif")
    else:
        filename = os.path.join(cif_dir, f"{comp_id}.cif")

    if not os.path.exists(filename):
        os.makedirs(cif_dir, exist_ok=True)
        url = f"https://files.rcsb.org/ligands/download/{comp_id}.cif"
        try:
            print(f"Downloading {comp_id}.cif from RCSB...")
            urlretrieve(url, filename)
        except Exception as ex:
            print(f"Could not download {comp_id}.cif: {ex}")
            return None

    if not os.path.exists(filename):
        return None

    data = []
    with open(filename, "r", encoding="utf-8") as fh:
        reader = pdbx.reader.PdbxReader(fh)
        reader.read(data)

    if not data:
        return None

    block = data[0]
    cif_data = {}

    for category_name in block.get_object_name_list():
        category = block.get_object(category_name)
        rows = []
        for i in range(category.row_count):
            row = {}
            for item_name in category.item_name_list:
                key = item_name.split('.', 1)[1]
                row[key] = category.get_value(key, i)
            rows.append(row)
        cif_data[category_name] = rows

    return cif_data


def get_component_coordinates_and_bonds(cif_data):
    """
    Coordinate and connection extraction adapted from refine_atom_mappings.py.
    """
    coordinates = {}
    coordinates_ideal = {}
    connections = []
    atom_to_element = {}

    for row in cif_data.get('chem_comp_atom', []):
        atom = row.get('atom_id')
        if not atom:
            continue

        x = row.get('model_Cartn_x')
        y = row.get('model_Cartn_y')
        z = row.get('model_Cartn_z')

        if x and y and z and x not in ['?', '.'] and y not in ['?', '.'] and z not in ['?', '.']:
            coordinates[atom] = [float(x), float(y), float(z)]

        atom_to_element[atom] = row.get('type_symbol', 'C')

        x = row.get('pdbx_model_Cartn_x_ideal')
        y = row.get('pdbx_model_Cartn_y_ideal')
        z = row.get('pdbx_model_Cartn_z_ideal')

        if x and y and z and x not in ['?', '.'] and y not in ['?', '.'] and z not in ['?', '.']:
            coordinates_ideal[atom] = [float(x), float(y), float(z)]

    # same rule as refine_atom_mappings.py
    if len(coordinates_ideal) > len(coordinates):
        coordinates = coordinates_ideal

    for row in cif_data.get('chem_comp_bond', []):
        atom1 = row.get('atom_id_1')
        atom2 = row.get('atom_id_2')
        if atom1 and atom2:
            connections.append((atom1, atom2))
            connections.append((atom2, atom1))

    connections = list(set(connections))

    return coordinates, connections, atom_to_element


def map_atom_to_neighbors(cif_data):
    """
    Create a dictionary whose keys are atoms and whose values are a set of neighboring atoms
    (ported/adapted from make_atom_mappings.py)
    """
    atom_to_neighbor_atoms = {}

    for row in cif_data.get('chem_comp_bond', []):
        atom1 = row.get('atom_id_1')
        atom2 = row.get('atom_id_2')
        if not atom1 or not atom2:
            continue

        if atom1 not in atom_to_neighbor_atoms:
            atom_to_neighbor_atoms[atom1] = set()
        if atom2 not in atom_to_neighbor_atoms:
            atom_to_neighbor_atoms[atom2] = set()

        atom_to_neighbor_atoms[atom1].add(atom2)
        atom_to_neighbor_atoms[atom2].add(atom1)

    return atom_to_neighbor_atoms


def extend_mappings_by_neighbors(mod_cif, par_cif, initial_parent_to_mod):
    """
    Given CIF data for modified and parent components and an initial mapping
    (parent_atom -> modified_atom), attempt to extend the mapping by
    propagating along bonded neighbors where names match.
    Returns an extended parent->mod mapping (dict).
    """
    # Build neighbor maps
    mod_atom_to_neighbors = map_atom_to_neighbors(mod_cif)
    par_atom_to_neighbors = map_atom_to_neighbors(par_cif)

    # initialize mappings; ignore blank entries as they are not real mappings
    standard_to_mod = {}
    if initial_parent_to_mod:
        for k, v in initial_parent_to_mod.items():
            if k and v:
                standard_to_mod[k] = v
    mod_to_par = {v: k for k, v in standard_to_mod.items()}

    # If both have C1' and it isn't set, seed mapping
    if ("C1'" not in standard_to_mod or not standard_to_mod.get("C1'")) and "C1'" in mod_atom_to_neighbors and "C1'" in par_atom_to_neighbors:
        standard_to_mod["C1'"] = "C1'"
        mod_to_par["C1'"] = "C1'"

    # Sensible seed for terminal/bridging phosphate oxygen naming differences.
    # Some modified residues use O3A where the parent uses OP3.
    if "OP3" not in standard_to_mod and "OP3" in par_atom_to_neighbors and "O3A" in mod_atom_to_neighbors:
        standard_to_mod["OP3"] = "O3A"
        mod_to_par["O3A"] = "OP3"

    # atoms to skip mapping automatically (hydrogens and mobile hydrogens)
    skip_atoms = set(['H21','H22','H41','H42','H61','H62','H71','H72','H73', "H5'","H5''","H2''"])

    new_par_atoms = list(standard_to_mod.keys())

    if len(new_par_atoms) == 0:
        return standard_to_mod

    unresolved_par_atoms = set()

    while len(new_par_atoms) > 0:
        par_atoms = new_par_atoms
        new_par_atoms = []
        for par_atom in par_atoms:
            mod_atom = standard_to_mod.get(par_atom)
            if not mod_atom:
                continue

            mod_neighbors = set(mod_atom_to_neighbors.get(mod_atom, set())) - set(mod_to_par.keys())
            par_neighbors = set(par_atom_to_neighbors.get(par_atom, set())) - set(standard_to_mod.keys())

            # match neighbors with identical names
            new_atoms = mod_neighbors & par_neighbors
            for new_atom in new_atoms:
                # avoid auto-mapping hydrogens here; defer H mapping to
                # the specialized hydrogen-mapping routines added later
                if new_atom in skip_atoms or (isinstance(new_atom, str) and new_atom.startswith('H')):
                    continue
                standard_to_mod[new_atom] = new_atom
                mod_to_par[new_atom] = new_atom
                new_par_atoms.append(new_atom)

            # handle common naming conventions (OP1/O1P etc.)
            if "OP1" in par_neighbors and "O1P" in mod_neighbors:
                standard_to_mod["OP1"] = "O1P"
                mod_to_par["O1P"] = "OP1"
                new_par_atoms.append("OP1")

            if "OP1" in par_neighbors and "O1A" in mod_neighbors:
                standard_to_mod["OP1"] = "O1A"
                mod_to_par["O1A"] = "OP1"
                new_par_atoms.append("OP1")

            if "OP2" in par_neighbors and "O2P" in mod_neighbors:
                standard_to_mod["OP2"] = "O2P"
                mod_to_par["O2P"] = "OP2"
                new_par_atoms.append("OP2")

            if "OP2" in par_neighbors and "O2A" in mod_neighbors:
                standard_to_mod["OP2"] = "O2A"
                mod_to_par["O2A"] = "OP2"
                new_par_atoms.append("OP2")

            if "OP3" in par_neighbors and "O3P" in mod_neighbors:
                standard_to_mod["OP3"] = "O3P"
                mod_to_par["O3P"] = "OP3"
                new_par_atoms.append("OP3")

            if "OP3" in par_neighbors and "O3A" in mod_neighbors:
                standard_to_mod["OP3"] = "O3A"
                mod_to_par["O3A"] = "OP3"
                new_par_atoms.append("OP3")

            # focus on remaining differing atoms
            mod_diff_atoms = mod_neighbors - new_atoms
            par_diff_atoms = par_neighbors - new_atoms

            if len(mod_diff_atoms) == 1 and len(par_diff_atoms) == 1:
                par_diff_atom = list(par_diff_atoms)[0]
                # Avoid mapping hydrogens automatically here; let the hydrogen
                # specific mapper handle cases where two hydrogens need pairing
                if isinstance(par_diff_atom, str) and par_diff_atom.startswith('H'):
                    pass
                else:
                    mod_diff_atom = list(mod_diff_atoms)[0]
                    standard_to_mod[par_diff_atom] = mod_diff_atom
                    mod_to_par[mod_diff_atom] = par_diff_atom
                    new_par_atoms.append(par_diff_atom)
            else:
                # accumulate unresolved non-hydrogen parent atoms for potential later inspection
                if len(mod_diff_atoms) > 0 or len(par_diff_atoms) > 0:
                    pda = set([x for x in par_diff_atoms if not x.startswith('H')])
                    mda = set([x for x in mod_diff_atoms if not x.startswith('H')])
                    if len(pda) == 2 and len(mda) == 2:
                        unresolved_par_atoms |= pda

    return standard_to_mod


def orient_component_to_standard(parent, coordinates):
    """
    Orient a component's coordinates into the standard base plane using
    the `NAbasecoordinates` template and `besttransformation`.
    Returns a new coordinates dict (atom -> [x,y,z]) in standard orientation.
    """
    try:
        template = NAbasecoordinates[parent]
    except Exception:
        return coordinates

    par_base_atoms = []
    par_atoms = []

    # collect matching heavy base atoms present in both template and coordinates
    for atom_name, tpl_coord in template.items():
        if atom_name.startswith('H'):
            continue
        if atom_name in coordinates:
            par_base_atoms.append(tpl_coord)
            par_atoms.append(coordinates[atom_name])

    if len(par_base_atoms) < 3:
        return coordinates

    try:
        U, new1, mean1, rmsd, sse, mean2 = besttransformation(par_base_atoms, par_atoms)
    except Exception:
        return coordinates

    oriented = {}
    mean1_flat = np.array(mean1).reshape(3,)
    mean2_flat = np.array(mean2).reshape(3,)
    for atom, coord in coordinates.items():
        arr = np.array(coord)
        try:
            c = mean1 + np.dot(U, arr - mean2)
            c = np.array(c).reshape(3,)
            oriented[atom] = [float(c[0]), float(c[1]), float(c[2])]
        except Exception:
            oriented[atom] = coord

    # Rotate so the base plane normal aligns with screen z-axis
    try:
        # collect three non-hydrogen template base atoms to compute normal
        tpl_atoms = []
        for atom_name, tpl_coord in NAbasecoordinates[parent].items():
            if not atom_name.startswith('H'):
                tpl_atoms.append(np.array(tpl_coord))
            if len(tpl_atoms) >= 3:
                break

        if len(tpl_atoms) >= 3:
            v1 = tpl_atoms[1] - tpl_atoms[0]
            v2 = tpl_atoms[2] - tpl_atoms[0]
            n_tpl = np.cross(v1, v2)
            n_norm = np.linalg.norm(n_tpl)
            if n_norm > 1e-8:
                n_tpl = n_tpl / n_norm
                target = np.array([0.0, 0.0, 1.0])
                # If already aligned, skip
                if abs(np.dot(n_tpl, target) - 1.0) > 1e-6:
                    # rotation axis to bring n_tpl to target
                    axis = np.cross(n_tpl, target)
                    axis_norm = np.linalg.norm(axis)
                    if axis_norm > 1e-8:
                        axis = axis / axis_norm
                        angle = np.arccos(np.clip(np.dot(n_tpl, target), -1.0, 1.0))
                        # Rodrigues' rotation matrix
                        K = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
                        R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K.dot(K))

                        # apply rotation around mean1_flat
                        rotated = {}
                        for atom, coord in oriented.items():
                            vec = np.array(coord) - mean1_flat
                            rv = R.dot(vec)
                            newc = mean1_flat + rv
                            rotated[atom] = [float(newc[0]), float(newc[1]), float(newc[2])]
                        return rotated
    except Exception:
        pass

    return oriented


def center_and_orient_base_xy(parent, coordinates):
    """
    Translate and rotate coordinates so that the base heavy-atom (N1 or N9)
    is at (0,0) in the XY plane and the C1' atom lies on the negative Y axis.
    Returns a new coordinates dict.
    """
    coords = dict(coordinates)

    # choose base atom: purines use N9, pyrimidines use N1
    base_atom = None
    if parent and parent.upper()[:1] in ['A','G','D']:  # crude test; refined below
        # prefer N9 for purines
        if "N9" in coords:
            base_atom = "N9"
        elif "N1" in coords:
            base_atom = "N1"
    else:
        if "N1" in coords:
            base_atom = "N1"
        elif "N9" in coords:
            base_atom = "N9"

    if base_atom is None or base_atom not in coords:
        return coords

    base_pos = np.array(coords[base_atom], dtype=float)

    # translate so base_pos xy -> (0,0)
    translated = {}
    for a, c in coords.items():
        v = np.array(c, dtype=float) - base_pos
        translated[a] = [float(v[0]), float(v[1]), float(v[2])]

    # ensure C1' is present and compute rotation to place it on negative Y axis
    if "C1'" in translated:
        c1 = np.array(translated["C1'"])
        v_xy = np.array([c1[0], c1[1]])
        norm = np.linalg.norm(v_xy)
        if norm > 1e-8:
            # current angle of vector in XY plane
            v_angle = math.atan2(v_xy[1], v_xy[0])
            # target angle for negative Y axis is -pi/2
            target_angle = -math.pi / 2.0
            rot = target_angle - v_angle
            cos_r = math.cos(rot)
            sin_r = math.sin(rot)
            Rz = np.array([[cos_r, -sin_r, 0.0], [sin_r, cos_r, 0.0], [0.0, 0.0, 1.0]])

            rotated = {}
            for a, c in translated.items():
                vec = np.dot(Rz, np.array(c, dtype=float))
                rotated[a] = [float(vec[0]), float(vec[1]), float(vec[2])]
            return rotated

    return translated


def draw_component_bonds_and_atoms(ax, coordinates, connections, atom_to_element, bond_linewidth=4.0, atom_labelsize=8, draw_atom_markers=False, atom_markersize=50, colorfn=None, color_map=None):
    """
    Reusable bond and atom drawing routine (no scatter balls, per-atom bond coloring).
    Can be used for multiple panels with different color schemes or data.
    """
    # draw bonds with per-atom coloring
    drawn_connections = set()
    for atom1, atom2 in connections:
        if atom1 not in coordinates or atom2 not in coordinates:
            continue

        key = tuple(sorted((atom1, atom2)))
        if key in drawn_connections:
            continue
        drawn_connections.add(key)

        p = np.array(coordinates[atom1])
        q = np.array(coordinates[atom2])

        if colorfn is None:
            colorfn = element_to_cpk_color

        # prefer atom-name-based color mapping if provided
        if color_map and atom1 in color_map:
            c1 = color_map.get(atom1)
        else:
            c1 = colorfn(atom_to_element.get(atom1, 'C'))

        if color_map and atom2 in color_map:
            c2 = color_map.get(atom2)
        else:
            c2 = colorfn(atom_to_element.get(atom2, 'C'))

        # Determine linewidth per-bond: if a color_map is provided that indicates
        # mapped atoms, draw bonds between two mapped atoms with the full
        # linewidth; if either atom is unmapped, draw a thinner bond so
        # unmapped atoms are visually de-emphasized.
        if color_map:
            if atom1 in color_map and atom2 in color_map:
                lw = bond_linewidth
            else:
                # use roughly one third thickness for unmapped bonds
                lw = max(0.5, bond_linewidth / 3.0)
        else:
            lw = bond_linewidth

        # points in order from atom1 to atom2
        s = [(3*p[0]+q[0])/4.0, (3*p[1]+q[1])/4.0, (3*p[2]+q[2])/4.0]  # quarter point, near atom1
        m = [(p[0]+q[0])/2.0, (p[1]+q[1])/2.0, (p[2]+q[2])/2.0]  # midpoint
        t = [(p[0]+3*q[0])/4.0, (p[1]+3*q[1])/4.0, (p[2]+3*q[2])/4.0]  # three-quarters point, near atom2

        try:
            if c1 == c2:
                # draw from atom1 to atom2 with same color
                ax.plot3D([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], color=c1, linewidth=lw, solid_capstyle='round')
            else:
                # draw from atom1 to quarter way with atom1 color
                ax.plot3D([p[0], s[0]], [p[1], s[1]], [p[2], s[2]], color=c1, linewidth=lw, solid_capstyle='round')
                # draw from quarter way to midpoint with atom1 color (no overlap at midpoint)
                ax.plot3D([s[0], m[0]], [s[1], m[1]], [s[2], m[2]], color=c1, linewidth=lw, solid_capstyle='butt')
                # draw from midpoint to three-quarters with atom2 color (no overlap at midpoint)
                ax.plot3D([m[0], t[0]], [m[1], t[1]], [m[2], t[2]], color=c2, linewidth=lw, solid_capstyle='butt')
                # draw from three-quarters to atom2 with atom2 color
                ax.plot3D([t[0], q[0]], [t[1], q[1]], [t[2], q[2]], color=c2, linewidth=lw, solid_capstyle='round')
        except Exception as ex:
            print(f"Error drawing bond {atom1}-{atom2}: {ex}")

    # draw atom markers (optional) and labels
    if draw_atom_markers:
        xs = []
        ys = []
        zs = []
        cols = []
        for atom, coord in coordinates.items():
            xs.append(coord[0])
            ys.append(coord[1])
            zs.append(coord[2])
            if color_map and atom in color_map:
                cols.append(color_map.get(atom))
            else:
                cols.append(colorfn(atom_to_element.get(atom, 'C')))
        try:
            ax.scatter(xs, ys, zs, c=cols, s=atom_markersize, depthshade=False)
        except Exception:
            # fallback to per-point drawing if 3D scatter not supported
            for x, y, z, c in zip(xs, ys, zs, cols):
                ax.plot3D([x], [y], [z], marker='o', color=c, markersize=max(1, atom_markersize/10))

    for atom, coord in coordinates.items():
        ax.text(coord[0] + 0.08, coord[1] + 0.08, coord[2], atom, fontsize=atom_labelsize, color='black')


def configure_axes_clean(ax):
    """Configure a 3D Axes to have no labels, no ticks, no panes/planes, and hidden axis lines."""
    try:
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_zlabel('')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
    except Exception:
        pass

    # Disable panes to remove coordinate planes
    try:
        ax.xaxis.pane.fill = False
        ax.yaxis.pane.fill = False
        ax.zaxis.pane.fill = False
        ax.xaxis.pane.set_edgecolor('none')
        ax.yaxis.pane.set_edgecolor('none')
        ax.zaxis.pane.set_edgecolor('none')
    except Exception:
        pass
    # As a last resort, turn the entire axis off to ensure no coordinate axes are shown
    try:
        # This hides axis lines, ticks, and labels for the axis object
        ax.set_axis_off()
    except Exception:
        pass


def read_new_atom_mappings(filename):
    """Read the new atom mappings format:
    std_base\tstd_atom\tmod_base\tmod_atom
    (also tolerates legacy 5-column lines with a trailing date)
    Returns same 3 dicts as read_atom_mappings plus not_mappable list.
    """
    parent_to_modified_atom = {}
    modified_to_parent_atom = {}
    modified_base_to_parent = {}
    not_mappable = []

    try:
        with open(filename, 'r') as f:
            lines = f.readlines()
    except Exception:
        return parent_to_modified_atom, modified_to_parent_atom, modified_base_to_parent, not_mappable

    for line in lines:
        fields = line.rstrip('\n').split('\t')
        # pad fields to length 4 (legacy 5th date column is ignored)
        while len(fields) < 4:
            fields.append('')
        std_base = fields[0].strip()
        std_atom = fields[1].strip()
        mod_base = fields[2].strip()
        mod_atom = fields[3].strip()

        if not mod_base:
            continue

        if std_base == '':
            # modified nucleotide marked as not mappable
            not_mappable.append(mod_base)
            continue

        if mod_base not in modified_to_parent_atom:
            modified_to_parent_atom[mod_base] = {}
            parent_to_modified_atom[mod_base] = {}
            modified_base_to_parent[mod_base] = std_base

        # store mapping; allow empty std_atom or mod_atom
        parent_to_modified_atom[mod_base][std_atom] = mod_atom
        if mod_atom:
            modified_to_parent_atom[mod_base][mod_atom] = std_atom

    return parent_to_modified_atom, modified_to_parent_atom, modified_base_to_parent, not_mappable


def write_new_atom_mappings_for_mod(mod_nt, parent_to_modified_atom, modified_to_parent_atom, modified_base_to_parent):
    """Update new_atom_mappings.txt replacing entries for mod_nt with current mappings."""
    fname = 'new_atom_mappings.txt'
    try:
        try:
            with open(fname, 'r') as f:
                old_lines = f.readlines()
        except Exception:
            old_lines = []

        # Keep lines not for this modified nucleotide
        kept = []
        for line in old_lines:
            fields = line.rstrip('\n').split('\t')
            if len(fields) >= 3 and fields[2].strip() == mod_nt:
                continue
            kept.append(line)

        # Build new lines for this modified nucleotide
        new_lines = []
        parent = modified_base_to_parent.get(mod_nt, None)

        if not parent:
            # write single-line marking unmapped modified nucleotide
            new_lines.append(f"\t\t{mod_nt}\t\n")
        else:
            # Prefer parent CIF atom names as the authoritative parent atom list.
            parent_atoms = []
            try:
                cif_par = read_component_from_file_or_download(parent)
                if cif_par:
                    par_coords, _, _ = get_component_coordinates_and_bonds(cif_par)
                    parent_atoms = list(par_coords.keys())
                elif parent in NAbasecoordinates:
                    # fallback only if CIF is unavailable
                    parent_atoms = list(NAbasecoordinates[parent].keys())

                # Remove O2' for DNA parent codes (DA/DC/DG/DT)
                try:
                    if parent in ['DA','DC','DG','DT'] and "O2'" in parent_atoms:
                        parent_atoms = [a for a in parent_atoms if a != "O2'"]
                except Exception:
                    pass
            except Exception:
                parent_atoms = []

            # If still empty, use any keys from existing mappings
            if not parent_atoms:
                parent_atoms = list(parent_to_modified_atom.get(mod_nt, {}).keys())

            # Preserve any existing *non-empty* mapped parent atoms that may not
            # be listed due to naming differences, so valid correspondences are
            # not dropped while stale blank placeholders are ignored.
            try:
                for pa, ma in parent_to_modified_atom.get(mod_nt, {}).items():
                    if pa and ma and pa not in parent_atoms:
                        parent_atoms.append(pa)
            except Exception:
                pass

            # Lines for each parent atom (may have empty mod atom)
            mapped_mods = set()
            for pa in sorted(parent_atoms):
                ma = parent_to_modified_atom.get(mod_nt, {}).get(pa, '')
                if ma:
                    mapped_mods.add(ma)
                new_lines.append(f"{parent}\t{pa}\t{mod_nt}\t{ma}\n")

            # Add any modified atoms that are not mapped to a parent atom
            try:
                cif_mod = read_component_from_file_or_download(mod_nt)
                if cif_mod:
                    mod_coords, _, _ = get_component_coordinates_and_bonds(cif_mod)
                    for ma in sorted(mod_coords.keys()):
                        if ma in mapped_mods:
                            continue
                        # standard atom field empty when modified atom unmapped
                        new_lines.append(f"{parent}\t\t{mod_nt}\t{ma}\n")
            except Exception:
                pass

        # Write back file: kept old lines + new lines
        try:
            with open(fname, 'w') as f:
                for L in kept:
                    f.write(L)
                for L in new_lines:
                    f.write(L)
        except Exception:
            pass
    except Exception:
        pass

    # Hide axis lines where possible
    # axis cleanup is handled by configure_axes_clean where axes exist


def apply_roll_rotation(coordinates, roll_angle_rad, center):
    """
    Apply roll rotation around z-axis through center.
    roll_angle_rad: rotation angle in radians
    center: [x, y, z] center point for rotation
    """
    if roll_angle_rad == 0:
        return coordinates

    cos_r = np.cos(roll_angle_rad)
    sin_r = np.sin(roll_angle_rad)

    rotated = {}
    for atom, coord in coordinates.items():
        # Translate to origin
        p = np.array(coord) - np.array(center)
        # Apply rotation matrix around z-axis
        x_rot = p[0] * cos_r - p[1] * sin_r
        y_rot = p[0] * sin_r + p[1] * cos_r
        z_rot = p[2]
        # Translate back
        rotated[atom] = [x_rot + center[0], y_rot + center[1], z_rot + center[2]]

    return rotated


def center_coords_to_centroid(coordinates):
    """Return a new coordinates dict with centroid subtracted (centered at origin)."""
    if not coordinates:
        return coordinates
    try:
        arr = np.array(list(coordinates.values()), dtype=float)
        centroid = arr.mean(axis=0)
        centered = {}
        for atom, coord in coordinates.items():
            c = np.array(coord, dtype=float) - centroid
            centered[atom] = [float(c[0]), float(c[1]), float(c[2])]
        return centered
    except Exception:
        return coordinates


def center_coords_on_atoms(coordinates, atom_names):
    """Return a new coordinates dict centered on the centroid of `atom_names`.
    Falls back to the input coordinates if none of the requested atoms are present.
    """
    if not coordinates or not atom_names:
        return coordinates
    try:
        selected = [np.array(coordinates[a], dtype=float) for a in atom_names if a in coordinates]
        if len(selected) == 0:
            return coordinates
        centroid = np.mean(np.array(selected, dtype=float), axis=0)
        centered = {}
        for atom, coord in coordinates.items():
            c = np.array(coord, dtype=float) - centroid
            centered[atom] = [float(c[0]), float(c[1]), float(c[2])]
        return centered
    except Exception:
        return coordinates


def get_atom_centroid(coordinates, atom_names):
    """Return centroid of selected atoms, or None if none are present."""
    if not coordinates or not atom_names:
        return None
    try:
        selected = [np.array(coordinates[a], dtype=float) for a in atom_names if a in coordinates]
        if len(selected) == 0:
            return None
        centroid = np.mean(np.array(selected, dtype=float), axis=0)
        return np.array(centroid, dtype=float)
    except Exception:
        return None


def unit_vector(v):
    return v / np.linalg.norm(v)


def pyramidal_hydrogens(P1, C, P2, bondLength=1.1):
    """Return two inferred hydrogen positions for a tetrahedral arrangement
    using vertices P1 and P2 around center C. Adapted from refine_atom_mappings.py.
    """
    V1 = P1
    V2 = P2
    u = unit_vector(C - V2)
    W = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
    R = np.identity(3) + (np.sqrt(3) / 2) * W + 1.5 * np.dot(W, W)
    V3 = C + bondLength * unit_vector(np.dot(R, V1 - C))
    V4 = C + bondLength * unit_vector(np.dot(np.transpose(R), V1 - C))

    V1 = P2
    V2 = P1
    u = unit_vector(C - V2)
    W = np.array([[0, -u[2], u[1]], [u[2], 0, -u[0]], [-u[1], u[0], 0]])
    R = np.identity(3) + (np.sqrt(3) / 2) * W + 1.5 * np.dot(W, W)
    VV4 = C + bondLength * unit_vector(np.dot(R, V1 - C))
    VV3 = C + bondLength * unit_vector(np.dot(np.transpose(R), V1 - C))

    P3 = (V3 + VV3) / 2
    P4 = (V4 + VV4) / 2

    return P3, P4


def get_mod_atom_closest_to(atom_list, parent_to_modified_atom_map, mod_coordinates, bond_length=1.1):
    """Use tetrahedral geometry to infer hydrogen locations based on three heavy atoms.
    Returns (left_atom, left_dist, left_coord, right_atom, right_dist, right_coord)
    or (None, None, None, None, None, None) if unable to infer.
    """
    c = []
    for a in atom_list:
        if a in parent_to_modified_atom_map:
            b = parent_to_modified_atom_map[a]
            if b in mod_coordinates:
                c.append(np.array(mod_coordinates[b]))
            else:
                return None, None, None, None, None, None
        else:
            return None, None, None, None, None, None

    # compute two candidate H positions
    p, q = pyramidal_hydrogens(c[0], c[1], c[2], bond_length)

    left_min_dist = 1e6
    left_min_atom = None
    right_min_dist = 1e6
    right_min_atom = None
    for b, dcoord in mod_coordinates.items():
        d = np.array(dcoord)
        rd = np.linalg.norm(p - d)
        if rd < right_min_dist:
            right_min_dist = rd
            right_min_atom = b
        ld = np.linalg.norm(q - d)
        if ld < left_min_dist:
            left_min_dist = ld
            left_min_atom = b

    return left_min_atom, left_min_dist, p, right_min_atom, right_min_dist, q


def fill_hydrogen_mappings(self, mod_nt, parent_name, parent_to_modified_atom, modified_to_parent_atom, cif_mod=None, cif_par=None):
    """Try to fill in hydrogen mappings for a single modified nucleotide using
    the base-superposition and geometry heuristics from `refine_atom_mappings.py`.
    Updates the mapping dicts in-place and returns number of hydrogen mappings added.
    """
    added = 0
    try:
        try:
            self.changes_text.insert('end', f"Running hydrogen filler for {mod_nt} (parent={parent_name})\n")
        except Exception:
            pass
        if cif_mod is None:
            cif_mod = read_component_from_file_or_download(mod_nt)
        if cif_par is None and parent_name:
            cif_par = read_component_from_file_or_download(parent_name)
        if not cif_mod or not cif_par or not parent_name:
            return added

        mod_coords, _, mod_atom_to_element = get_component_coordinates_and_bonds(cif_mod)
        par_coords, _, par_atom_to_element = get_component_coordinates_and_bonds(cif_par)

        # Build superposition using heavy base atoms that are currently mapped
        mapped_pairs = []
        for ma, pa in modified_to_parent_atom.get(mod_nt, {}).items():
            try:
                if pa in NAbasecoordinates.get(parent_name, {}) and not pa.startswith('H') and ma in mod_coords:
                    mapped_pairs.append((pa, ma))
            except Exception:
                continue

        try:
            self.changes_text.insert('end', f"Mapped heavy pairs before H-fill: {mapped_pairs}\n")
        except Exception:
            pass

        if len(mapped_pairs) >= 3:
            par_list = [np.array(NAbasecoordinates[parent_name][pa], dtype=float) for pa, ma in mapped_pairs]
            mod_list = [np.array(mod_coords[ma], dtype=float) for pa, ma in mapped_pairs]
            try:
                U, new1, mean1, rmsd, sse, mean2 = besttransformation(par_list, mod_list)
                # transform modified coordinates into parent-based standard frame
                mod_coords_std = {}
                for atom, coord in mod_coords.items():
                    arr = np.array(coord)
                    try:
                        c = mean1 + np.dot(U, arr - mean2)
                        c = np.array(c).reshape(3,)
                        mod_coords_std[atom] = [float(c[0]), float(c[1]), float(c[2])]
                    except Exception:
                        mod_coords_std[atom] = coord
            except Exception:
                mod_coords_std = dict(mod_coords)
        else:
            # nothing to base a superposition on; use raw coords
            mod_coords_std = dict(mod_coords)

        # Also build superposition directly in parent-CIF frame for reliable H mapping.
        # This avoids relying on template hydrogen coordinates for residues where
        # template hydrogen naming/availability differs from CIF naming.
        mod_coords_parent_frame = dict(mod_coords)
        try:
            mapped_pairs_cif = []
            for pa, ma in parent_to_modified_atom.get(mod_nt, {}).items():
                if not pa or not ma:
                    continue
                if pa.startswith('H'):
                    continue
                if pa in par_coords and ma in mod_coords:
                    mapped_pairs_cif.append((pa, ma))

            if len(mapped_pairs_cif) >= 3:
                par_list_cif = [np.array(par_coords[pa], dtype=float) for pa, ma in mapped_pairs_cif]
                mod_list_cif = [np.array(mod_coords[ma], dtype=float) for pa, ma in mapped_pairs_cif]
                Uc, new1c, mean1c, rmsdc, ssec, mean2c = besttransformation(par_list_cif, mod_list_cif)
                tmp = {}
                for atom, coord in mod_coords.items():
                    arr = np.array(coord)
                    c = mean1c + np.dot(Uc, arr - mean2c)
                    c = np.array(c).reshape(3,)
                    tmp[atom] = [float(c[0]), float(c[1]), float(c[2])]
                mod_coords_parent_frame = tmp
        except Exception:
            mod_coords_parent_frame = dict(mod_coords)

        # Build neighbor maps to prefer direct neighbor-based hydrogen assignments
        mod_atom_to_neighbors = map_atom_to_neighbors(cif_mod)
        par_atom_to_neighbors = map_atom_to_neighbors(cif_par)

        pmap = parent_to_modified_atom.get(mod_nt, {})
        rmap = modified_to_parent_atom.get(mod_nt, {})

        # Map BASE hydrogens using parent-CIF neighbor topology + CIF-frame proximity.
        # This handles cases like DG:H1/H21/H22/H8 vs 1GC:HN1/HN2/HN2A/H8.
        base_heavy = set(NAbaseheavyatoms.get(parent_name, []))
        if len(base_heavy) == 0 and isinstance(parent_name, str) and parent_name.startswith('D'):
            base_heavy = set(NAbaseheavyatoms.get(parent_name[1:], []))

        # Build a reusable pool of currently-unmapped modified hydrogens
        def _unmapped_mod_hydrogens():
            return [m for m, el in mod_atom_to_element.items() if str(el).upper() == 'H' and m not in rmap]
        for par_atom, el in par_atom_to_element.items():
            if not par_atom.startswith('H'):
                continue
            # Only skip if this parent H already has a non-empty mapping.
            if pmap.get(par_atom, ''):
                continue

            par_neighbors = par_atom_to_neighbors.get(par_atom, set())
            heavy_neighbors = [n for n in par_neighbors if not n.startswith('H')]
            if len(heavy_neighbors) == 0:
                continue

            # only base-hydrogen stage: heavy neighbor must be in base heavy set
            if not any(hn in base_heavy for hn in heavy_neighbors):
                continue

            best = None
            bestd = 1e6
            for heavy in heavy_neighbors:
                if heavy not in base_heavy:
                    continue
                mod_heavy = pmap.get(heavy)
                if not mod_heavy:
                    continue
                candHs = [m for m in mod_atom_to_neighbors.get(mod_heavy, set()) if (m not in rmap) and (mod_atom_to_element.get(m, '').upper() == 'H')]
                if len(candHs) == 0:
                    continue

                if par_atom in par_coords:
                    pcoord = np.array(par_coords[par_atom], dtype=float)
                else:
                    pcoord = None

                if pcoord is None:
                    # no parent coordinate, choose first available candidate
                    cand = candHs[0]
                    d = 0.0
                    if d < bestd:
                        bestd = d
                        best = cand
                else:
                    for cand in candHs:
                        if cand not in mod_coords_parent_frame:
                            continue
                        d = np.linalg.norm(pcoord - np.array(mod_coords_parent_frame[cand], dtype=float))
                        if d < bestd:
                            bestd = d
                            best = cand

            # fallback for base H: if no bonded-candidate was found, use nearest
            # currently-unmapped modified hydrogen in parent-CIF frame.
            if best is None and par_atom in par_coords:
                pcoord = np.array(par_coords[par_atom], dtype=float)
                for cand in _unmapped_mod_hydrogens():
                    if cand not in mod_coords_parent_frame:
                        continue
                    d = np.linalg.norm(pcoord - np.array(mod_coords_parent_frame[cand], dtype=float))
                    if d < bestd:
                        bestd = d
                        best = cand

            if best is not None and bestd < 2.2:
                parent_to_modified_atom.setdefault(mod_nt, {})[par_atom] = best
                modified_to_parent_atom.setdefault(mod_nt, {})[best] = par_atom
                rmap[best] = par_atom
                added += 1
                try:
                    self.changes_text.insert('end', f"Mapped base hydrogen {par_atom}->{best} (dist {bestd:.2f})\n")
                except Exception:
                    pass

        # Map straightforward non-base/non-methylene backbone hydrogens by
        # mapped-heavy-neighbor topology, e.g. H3', H4', H1', HO3'.
        special_methylene = {"H5'", "H5''", "H2'", "H2''"}
        for par_atom, el in par_atom_to_element.items():
            if not par_atom.startswith('H'):
                continue
            if par_atom in special_methylene:
                continue
            # skip if already mapped to a non-empty modified atom
            if pmap.get(par_atom, ''):
                continue

            par_neighbors = par_atom_to_neighbors.get(par_atom, set())
            heavy_neighbors = [n for n in par_neighbors if not n.startswith('H')]
            if len(heavy_neighbors) == 0:
                continue

            # Use first heavy neighbor that is already mapped
            mod_heavy = None
            for hn in heavy_neighbors:
                if pmap.get(hn, ''):
                    mod_heavy = pmap.get(hn)
                    break
            if not mod_heavy:
                continue

            candHs = [m for m in mod_atom_to_neighbors.get(mod_heavy, set()) if (m not in rmap) and (mod_atom_to_element.get(m, '').upper() == 'H')]
            if len(candHs) == 0:
                continue

            chosen = None
            chosen_dist = None
            if len(candHs) == 1:
                chosen = candHs[0]
                chosen_dist = 0.0
            elif par_atom in par_coords:
                pcoord = np.array(par_coords[par_atom], dtype=float)
                bestd = 1e6
                best = None
                for cand in candHs:
                    if cand not in mod_coords_parent_frame:
                        continue
                    d = np.linalg.norm(pcoord - np.array(mod_coords_parent_frame[cand], dtype=float))
                    if d < bestd:
                        bestd = d
                        best = cand
                chosen = best
                chosen_dist = bestd

            if chosen:
                parent_to_modified_atom.setdefault(mod_nt, {})[par_atom] = chosen
                modified_to_parent_atom.setdefault(mod_nt, {})[chosen] = par_atom
                rmap[chosen] = par_atom
                added += 1
                try:
                    if chosen_dist is None:
                        self.changes_text.insert('end', f"Mapped backbone hydrogen {par_atom}->{chosen}\n")
                    else:
                        self.changes_text.insert('end', f"Mapped backbone hydrogen {par_atom}->{chosen} (dist {chosen_dist:.2f})\n")
                except Exception:
                    pass

        # Backbone hydrogens heuristics (H5', H5'', H2', H2'')
        hydrogen_to_heavy = {
            "H5'": ["C4'", "C5'", "O5'"],
            "H5''": ["O5'", "C5'", "C4'"],
            "H2'": ["C1'", "C2'", "C3'"],
            "H2''": ["C3'", "C2'", "C1'"],
        }

        for hyd, heavy_list in hydrogen_to_heavy.items():
            # skip if parent doesn't have that hydrogen
            if hyd not in par_atom_to_element:
                continue
            # skip only if already mapped to a non-empty modified atom
            if parent_to_modified_atom.get(mod_nt, {}).get(hyd, ''):
                continue

            # Try tetrahedral inference first (left/right candidates)
            left_atom, left_dist, left_coord, right_atom, right_dist, right_coord = get_mod_atom_closest_to(heavy_list, parent_to_modified_atom.get(mod_nt, {}), mod_coords)
            # Prefer candidate closer to parent CIF hydrogen coordinate.
            pcoord = None
            try:
                if hyd in par_coords:
                    pcoord = np.array(par_coords[hyd], dtype=float)
            except Exception:
                pcoord = None

            candidate = None
            cand_dist = None

            # For methylene hydrogens, use a deterministic geometry-side convention
            # from the pyramidal construction (right-side first). This avoids
            # ambiguous name-based swaps like H5'/H5'' and H2'/H2''.
            if hyd in ["H5'", "H5''", "H2'", "H2''"]:
                if right_atom and right_atom not in rmap:
                    candidate = right_atom
                    cand_dist = right_dist
                elif left_atom and left_atom not in rmap:
                    candidate = left_atom
                    cand_dist = left_dist

            if pcoord is not None:
                # compute distances from template hydrogen to candidate mod atom positions
                d_left = 1e6
                d_right = 1e6
                if left_atom and left_atom in mod_coords_parent_frame:
                    d_left = np.linalg.norm(pcoord - np.array(mod_coords_parent_frame[left_atom], dtype=float))
                if right_atom and right_atom in mod_coords_parent_frame:
                    d_right = np.linalg.norm(pcoord - np.array(mod_coords_parent_frame[right_atom], dtype=float))
                # choose closer candidate that isn't already used when no deterministic choice yet
                if candidate is None:
                    if left_atom and left_atom not in rmap and d_left <= d_right:
                        candidate = left_atom
                        cand_dist = d_left
                    elif right_atom and right_atom not in rmap and d_right < d_left:
                        candidate = right_atom
                        cand_dist = d_right
            else:
                # fallback: prefer left then right if available (legacy behavior)
                if candidate is None:
                    if left_atom and left_atom not in rmap:
                        candidate = left_atom
                        cand_dist = left_dist
                    elif right_atom and right_atom not in rmap:
                        candidate = right_atom
                        cand_dist = right_dist

            if candidate and cand_dist < 1.8:
                parent_to_modified_atom.setdefault(mod_nt, {})[hyd] = candidate
                modified_to_parent_atom.setdefault(mod_nt, {})[candidate] = hyd
                rmap[candidate] = hyd
                added += 1
                try:
                    self.changes_text.insert('end', f"Mapped backbone hydrogen {hyd}->{candidate} (dist {cand_dist:.2f})\n")
                except Exception:
                    pass

        try:
            self.changes_text.insert('end', f"Hydrogen filler added {added} mappings for {mod_nt}\n")
        except Exception:
            pass

    except Exception:
        pass

    return added


class MapModifiedGUI:
    def __init__(self, master, mod_nt, batch_mode=False):
        self.master = master
        self.batch_mode = batch_mode
        self.mod_nt = mod_nt
        self.std_nt = "??"  # Default standard nucleotide placeholder
        self.edit_text = self.std_nt + "\t" + self.mod_nt
        self.return_data = None
        self.return_action = None
        self.roll_angle = 0.0  # for third rotation axis via right-click
        self.right_click_start = None  # track right-click drag
        self.left_mouse_rotate_active = False
        self.pan_state = None

        # Set up the main window
        self.master.title(f"Map Modified Nucleotide: {mod_nt}")
        self.master.state('zoomed')  # Maximize window on Windows
        try:
            self.master.protocol("WM_DELETE_WINDOW", self.on_window_close)
        except Exception:
            pass

        # Create main layout
        self.create_layout()

        # Draw initial cubes (upper-right component will be drawn later when processed)
        self.draw_upper_left_cube()
        self.draw_lower_left_cube()
        self.draw_lower_right_cube()

        # Track whether upper-right component has been drawn yet
        self.upper_right_drawn = False
        # Track whether right panels have been auto-aligned to the parent orientation
        self.right_panels_aligned = False

        self.load_mod_nt(mod_nt)

    def load_mod_nt(self, mod_nt):
        """Load a new CCD into the existing GUI without recreating the window."""
        self.mod_nt = mod_nt
        self.edit_text = self.std_nt + "\t" + self.mod_nt
        self.return_data = None
        self.return_action = None
        self.right_click_start = None
        self.pan_state = None
        self.left_mouse_rotate_active = False
        self.right_panels_aligned = False
        try:
            if hasattr(self, 'right_aligned_coords'):
                del self.right_aligned_coords
        except Exception:
            pass
        try:
            if hasattr(self, 'right_alignment_signature'):
                del self.right_alignment_signature
        except Exception:
            pass
        try:
            self.master.title(f"Map Modified Nucleotide: {mod_nt}")
        except Exception:
            pass

        try:
            self.changes_text.delete('1.0', 'end')
        except Exception:
            pass

        # Show CCD self-identification hints from the CIF file in the Changes box.
        # These are informational only and are not used to auto-select the parent.
        try:
            cif_data_identity = read_component_from_file_or_download(mod_nt)
            if cif_data_identity:
                self.log_cif_identity_fields(mod_nt, cif_data_identity)
        except Exception:
            pass

        # If there is a mapping for this modified nucleotide in atom_mappings.txt,
        # populate the central mappings text area and draw the standard nucleotide
        try:
            parent_to_modified_atom, modified_to_parent_atom, modified_base_to_parent, not_mappable = read_atom_mappings("atom_mappings.txt")
        except Exception:
            parent_to_modified_atom = {}
            modified_to_parent_atom = {}
            modified_base_to_parent = {}
            not_mappable = []

        # Merge in any new mappings provided by new_atom_mappings.txt
        try:
            new_pmap, new_m2p, new_mb2p, new_not_mappable = read_new_atom_mappings("new_atom_mappings.txt")
            # overwrite or add entries from new file
            for mod_base, pmap in new_pmap.items():
                parent_to_modified_atom[mod_base] = pmap
            for mod_base, rev in new_m2p.items():
                modified_to_parent_atom[mod_base] = rev
            for mod_base, parent in new_mb2p.items():
                modified_base_to_parent[mod_base] = parent
            # merge not_mappable lists
            not_mappable = list(set(not_mappable) | set(new_not_mappable))
        except Exception:
            pass

        # store mapping tables on the GUI instance for interactive updates
        self.parent_to_modified_atom = parent_to_modified_atom
        self.modified_to_parent_atom = modified_to_parent_atom
        self.modified_base_to_parent = modified_base_to_parent

        if mod_nt in modified_base_to_parent:
            std_base = modified_base_to_parent.get(mod_nt)
        else:
            std_base = None
        # store standard base for use by other methods
        self.std_base = std_base

        # If we have a standard base, ensure the lower-left view is drawn using it
        if self.std_base:
            try:
                self.draw_lower_left_component()
            except Exception:
                pass

        # If mappings exist for this modified residue, populate central column
        if mod_nt in parent_to_modified_atom:
            mapping_lines = []
            # First line: reference nucleotide tab modified residue
            if std_base:
                mapping_lines.append(f"{std_base}\t{mod_nt}")
            else:
                mapping_lines.append(f"??\t{mod_nt}")

            # Then per-atom mappings: parent atom \t modified atom
            pmap = parent_to_modified_atom.get(mod_nt, {})
            # Sort keys for deterministic output
            for parent_atom in sorted(pmap.keys()):
                modified_atom = pmap.get(parent_atom, '')
                mapping_lines.append(f"{parent_atom}\t{modified_atom}")

            # Replace content of mappings_text
            try:
                self.mappings_text.delete('1.0', 'end')
                self.mappings_text.insert('1.0', '\n'.join(mapping_lines))
            except Exception:
                pass

            # If there is a standard nucleotide code, draw it in the upper-left window
            # Use the shared draw_upper_left_component so orientation and centering
            # match the lower-left implementation.
            if std_base:
                try:
                    self.draw_upper_left_component()
                except Exception:
                    pass

            # If mappings exist for this modified residue, also draw the modified
            # residue in the lower-right window using the mapping-aware colors
            # so the GUI shows mapped coloring on startup.
            try:
                self.draw_lower_right_component()
            except Exception:
                pass
        else:
            try:
                self.mappings_text.delete('1.0', 'end')
                if std_base:
                    self.mappings_text.insert('1.0', f"{std_base}\t{mod_nt}")
                else:
                    self.mappings_text.insert('1.0', f"??\t{mod_nt}")
            except Exception:
                pass
            try:
                if std_base:
                    self.draw_upper_left_component()
                    self.draw_lower_left_component()
                else:
                    self.draw_upper_left_cube()
                    self.draw_lower_left_cube()
            except Exception:
                pass
            try:
                self.draw_lower_right_component()
            except Exception:
                pass

        # Draw the upper-right component at startup (once)
        try:
            self.draw_upper_right_component()
            self.upper_right_drawn = True
        except Exception:
            pass

    def log_cif_identity_fields(self, mod_nt, cif_data):
        """Log parent-related identity fields from the component CIF to Changes box."""
        try:
            chem_comp_rows = cif_data.get('chem_comp', [])
            if not chem_comp_rows:
                self.changes_text.insert('end', f"{mod_nt}: no chem_comp category found in CIF\n")
                return

            row = chem_comp_rows[0] if len(chem_comp_rows) > 0 else {}

            def clean(v):
                if v is None:
                    return ''
                s = str(v).strip()
                if s in ['.', '?']:
                    return ''
                return s

            f_type = clean(row.get('type', ''))
            f_one = clean(row.get('one_letter_code', ''))
            f_parent = clean(row.get('mon_nstd_parent_comp_id', ''))
            f_pdbx_type = clean(row.get('pdbx_type', ''))

            self.changes_text.insert('end', f"CIF identity hints for {mod_nt}:\n")
            self.changes_text.insert('end', f"  _chem_comp.type = {f_type}\n")
            self.changes_text.insert('end', f"  _chem_comp.one_letter_code = {f_one}\n")
            self.changes_text.insert('end', f"  _chem_comp.mon_nstd_parent_comp_id = {f_parent}\n")
            self.changes_text.insert('end', f"  _chem_comp.pdbx_type = {f_pdbx_type}\n")
        except Exception:
            pass

    def create_layout(self):
        """Create the GUI layout with 4 graphics windows, 2 text boxes, and 3 buttons"""

        # Create main container
        main_frame = ttk.Frame(self.master)
        main_frame.pack(fill=tk.BOTH, expand=True, padx=5, pady=5)

        # Configure grid weights for responsive layout
        main_frame.columnconfigure(0, weight=2)  # Left graphics
        main_frame.columnconfigure(1, weight=1)  # Center controls
        main_frame.columnconfigure(2, weight=2)  # Right graphics
        main_frame.rowconfigure(0, weight=1)  # Upper graphics
        main_frame.rowconfigure(1, weight=1)  # Lower graphics

        # Create frames for each section
        self.upper_left_frame = ttk.Frame(main_frame, relief=tk.SUNKEN, borderwidth=2)
        self.upper_left_frame.grid(row=0, column=0, sticky='nsew', padx=2, pady=2)

        self.lower_left_frame = ttk.Frame(main_frame, relief=tk.SUNKEN, borderwidth=2)
        self.lower_left_frame.grid(row=1, column=0, sticky='nsew', padx=2, pady=2)

        self.center_frame = ttk.Frame(main_frame)
        self.center_frame.grid(row=0, column=1, rowspan=2, sticky='nsew', padx=5, pady=2)

        self.upper_right_frame = ttk.Frame(main_frame, relief=tk.SUNKEN, borderwidth=2)
        self.upper_right_frame.grid(row=0, column=2, sticky='nsew', padx=2, pady=2)

        self.lower_right_frame = ttk.Frame(main_frame, relief=tk.SUNKEN, borderwidth=2)
        self.lower_right_frame.grid(row=1, column=2, sticky='nsew', padx=2, pady=2)

        # Set up center controls
        self.setup_center_controls()

        # Create 3D graphics windows
        self.upper_left_canvas = self.create_3d_plot(self.upper_left_frame)
        self.lower_left_canvas = self.create_3d_plot(self.lower_left_frame)
        self.upper_right_canvas = self.create_3d_plot(self.upper_right_frame)
        self.lower_right_canvas = self.create_3d_plot(self.lower_right_frame)

        # Bind mouse wheel zoom and right-click roll rotation to upper right canvas
        self.upper_right_canvas.get_tk_widget().bind("<MouseWheel>", self.on_upper_right_scroll)
        self.upper_right_canvas.get_tk_widget().bind("<Button-4>", self.on_upper_right_scroll)  # Linux scroll up
        self.upper_right_canvas.get_tk_widget().bind("<Button-5>", self.on_upper_right_scroll)  # Linux scroll down
        self.upper_right_canvas.get_tk_widget().bind("<Button-3>", self.on_upper_right_button_press)  # Right-click press
        self.upper_right_canvas.get_tk_widget().bind("<B3-Motion>", self.on_upper_right_motion)  # Right-click drag

        # Bind same controls to upper-left canvas so user can zoom/roll there as well
        self.upper_left_canvas.get_tk_widget().bind("<MouseWheel>", self.on_upper_left_scroll)
        self.upper_left_canvas.get_tk_widget().bind("<Button-4>", self.on_upper_left_scroll)
        self.upper_left_canvas.get_tk_widget().bind("<Button-5>", self.on_upper_left_scroll)
        self.upper_left_canvas.get_tk_widget().bind("<Button-3>", self.on_upper_left_button_press)
        self.upper_left_canvas.get_tk_widget().bind("<B3-Motion>", self.on_upper_left_motion)

        # Bind same controls to lower-left and lower-right canvases so they have the same zoom/roll behavior
        try:
            self.lower_left_canvas.get_tk_widget().bind("<MouseWheel>", self.on_upper_left_scroll)
            self.lower_left_canvas.get_tk_widget().bind("<Button-4>", self.on_upper_left_scroll)
            self.lower_left_canvas.get_tk_widget().bind("<Button-5>", self.on_upper_left_scroll)
            self.lower_left_canvas.get_tk_widget().bind("<Button-3>", self.on_upper_left_button_press)
            self.lower_left_canvas.get_tk_widget().bind("<B3-Motion>", self.on_upper_left_motion)

            self.lower_right_canvas.get_tk_widget().bind("<MouseWheel>", self.on_upper_left_scroll)
            self.lower_right_canvas.get_tk_widget().bind("<Button-4>", self.on_upper_left_scroll)
            self.lower_right_canvas.get_tk_widget().bind("<Button-5>", self.on_upper_left_scroll)
            self.lower_right_canvas.get_tk_widget().bind("<Button-3>", self.on_upper_left_button_press)
            self.lower_right_canvas.get_tk_widget().bind("<B3-Motion>", self.on_upper_left_motion)
        except Exception:
            pass

    def on_upper_right_button_press(self, event):
        """Track right-click for panning."""
        button = getattr(event, 'num', getattr(event, 'button', None))
        if button == 3:  # Right-click (button 3)
            self.right_click_start = (event.x, event.y)
            self.pan_state = self.get_canvas_from_widget(event.widget)

    def on_upper_left_button_press(self, event):
        """Track right-click for panning on any synced panel."""
        button = getattr(event, 'num', getattr(event, 'button', None))
        if button == 3:
            self.right_click_start = (event.x, event.y)
            self.pan_state = self.get_canvas_from_widget(event.widget)

    def on_upper_right_button_release(self, event):
        """End right-click pan drag."""
        button = getattr(event, 'num', getattr(event, 'button', None))
        if button == 3:
            self.right_click_start = None
            self.pan_state = None

    def on_upper_right_motion(self, event):
        """Handle right-click drag for panning."""
        if self.right_click_start is None or not self.pan_state:
            return

        dx = event.x - self.right_click_start[0]
        dy = event.y - self.right_click_start[1]
        self.right_click_start = (event.x, event.y)
        self.pan_all_axes(self.pan_state.axes, dx, dy)

    def on_upper_left_motion(self, event):
        """Handle right-click drag for panning from any synced panel."""
        if self.right_click_start is None or not self.pan_state:
            return

        dx = event.x - self.right_click_start[0]
        dy = event.y - self.right_click_start[1]
        self.right_click_start = (event.x, event.y)
        self.pan_all_axes(self.pan_state.axes, dx, dy)

    def on_upper_right_scroll(self, event):
        """Handle scroll zoom for upper right panel"""

        ax = self.upper_right_canvas.axes
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        cur_zlim = ax.get_zlim()

        # Determine zoom direction
        if event.num == 5 or event.delta < 0:
            scale_factor = 0.8  # Zoom out
        else:
            scale_factor = 1.2  # Zoom in

        relx = (cur_xlim[1] - cur_xlim[0]) / 2
        rely = (cur_ylim[1] - cur_ylim[0]) / 2
        relz = (cur_zlim[1] - cur_zlim[0]) / 2

        ax.set_xlim([cur_xlim[0] - relx * (scale_factor - 1), cur_xlim[1] + relx * (scale_factor - 1)])
        ax.set_ylim([cur_ylim[0] - rely * (scale_factor - 1), cur_ylim[1] + rely * (scale_factor - 1)])
        ax.set_zlim([cur_zlim[0] - relz * (scale_factor - 1), cur_zlim[1] + relz * (scale_factor - 1)])

        # Disable autoscale and margin adjustment to prevent clipping
        ax.autoscale(False)
        ax.margins(0)

        # Force draw without automatic axis adjustment to prevent clipping
        # Propagate same limits to other canvases for synchronized zoom
        try:
            ul_ax = self.upper_left_canvas.axes
            ll_ax = self.lower_left_canvas.axes
            lr_ax = self.lower_right_canvas.axes
            ul_ax.set_xlim(ax.get_xlim())
            ul_ax.set_ylim(ax.get_ylim())
            ul_ax.set_zlim(ax.get_zlim())
            ul_ax.autoscale(False)
            ul_ax.margins(0)
            ll_ax.set_xlim(ax.get_xlim())
            ll_ax.set_ylim(ax.get_ylim())
            ll_ax.set_zlim(ax.get_zlim())
            ll_ax.autoscale(False)
            ll_ax.margins(0)
            lr_ax.set_xlim(ax.get_xlim())
            lr_ax.set_ylim(ax.get_ylim())
            lr_ax.set_zlim(ax.get_zlim())
            lr_ax.autoscale(False)
            lr_ax.margins(0)
            try:
                self.upper_left_canvas.draw_idle()
            except Exception:
                pass
            try:
                self.lower_left_canvas.draw_idle()
            except Exception:
                pass
            try:
                self.lower_right_canvas.draw_idle()
            except Exception:
                pass
        except Exception:
            pass

        self.upper_right_canvas.draw_idle()

    def on_upper_left_scroll(self, event):
        """Handle scroll zoom for upper left panel and sync to other panels"""

        ax = self.upper_left_canvas.axes
        cur_xlim = ax.get_xlim()
        cur_ylim = ax.get_ylim()
        cur_zlim = ax.get_zlim()

        if event.num == 5 or event.delta < 0:
            scale_factor = 0.8
        else:
            scale_factor = 1.2

        relx = (cur_xlim[1] - cur_xlim[0]) / 2
        rely = (cur_ylim[1] - cur_ylim[0]) / 2
        relz = (cur_zlim[1] - cur_zlim[0]) / 2

        ax.set_xlim([cur_xlim[0] - relx * (scale_factor - 1), cur_xlim[1] + relx * (scale_factor - 1)])
        ax.set_ylim([cur_ylim[0] - rely * (scale_factor - 1), cur_ylim[1] + rely * (scale_factor - 1)])
        ax.set_zlim([cur_zlim[0] - relz * (scale_factor - 1), cur_zlim[1] + relz * (scale_factor - 1)])

        ax.autoscale(False)
        ax.margins(0)

        # Sync these limits to all other windows so they show the same zoom
        try:
            ur_ax = self.upper_right_canvas.axes
            lr_ax = self.lower_right_canvas.axes
            ll_ax = self.lower_left_canvas.axes
            ur_ax.set_xlim(ax.get_xlim())
            ur_ax.set_ylim(ax.get_ylim())
            ur_ax.set_zlim(ax.get_zlim())
            ur_ax.autoscale(False)
            ur_ax.margins(0)
            lr_ax.set_xlim(ax.get_xlim())
            lr_ax.set_ylim(ax.get_ylim())
            lr_ax.set_zlim(ax.get_zlim())
            lr_ax.autoscale(False)
            lr_ax.margins(0)
            ll_ax.set_xlim(ax.get_xlim())
            ll_ax.set_ylim(ax.get_ylim())
            ll_ax.set_zlim(ax.get_zlim())
            ll_ax.autoscale(False)
            ll_ax.margins(0)
            try:
                self.upper_right_canvas.draw_idle()
            except Exception:
                pass
            try:
                self.lower_left_canvas.draw_idle()
            except Exception:
                pass
            try:
                self.lower_right_canvas.draw_idle()
            except Exception:
                pass
        except Exception:
            pass

        self.upper_left_canvas.draw_idle()

    def setup_center_controls(self):
        """Set up the center panel with text boxes and buttons"""

        # Configure center frame layout
        self.center_frame.rowconfigure(0, weight=3)  # Mappings text box
        self.center_frame.rowconfigure(1, weight=0)  # Buttons
        self.center_frame.rowconfigure(2, weight=1)  # Changes text box
        self.center_frame.columnconfigure(0, weight=1)

        # Mappings text box (top, editable)
        mappings_label = ttk.Label(self.center_frame, text="Mappings:")
        mappings_label.grid(row=0, column=0, sticky='nw', pady=(0, 2))

        self.mappings_text = tk.Text(self.center_frame, width=30, height=10, wrap=tk.WORD)
        self.mappings_text.configure(height=16)
        self.mappings_text.grid(row=0, column=0, sticky='nsew', pady=(20, 5))
        self.mappings_text.insert('1.0', self.edit_text)

        # Add scrollbar to mappings text
        mappings_scroll = ttk.Scrollbar(self.center_frame, orient='vertical', command=self.mappings_text.yview)
        mappings_scroll.grid(row=0, column=1, sticky='ns', pady=(20, 5))
        self.mappings_text.configure(yscrollcommand=mappings_scroll.set)

        # Buttons frame
        buttons_frame = ttk.Frame(self.center_frame)
        buttons_frame.grid(row=1, column=0, sticky='ew', pady=5)
        buttons_frame.columnconfigure(0, weight=1)
        buttons_frame.columnconfigure(1, weight=1)
        buttons_frame.rowconfigure(0, weight=1)
        buttons_frame.rowconfigure(1, weight=1)

        # Create buttons (two rows)
        process_btn = ttk.Button(buttons_frame, text="Process", command=self.on_process)
        process_btn.grid(row=0, column=0, padx=2, sticky='ew')

        save_btn = ttk.Button(buttons_frame, text="Save", command=self.on_save)
        save_btn.grid(row=0, column=1, padx=2, sticky='ew')

        skip_btn = ttk.Button(buttons_frame, text="Skip", command=self.on_skip)
        skip_btn.grid(row=1, column=0, padx=2, sticky='ew')

        unmappable_btn = ttk.Button(buttons_frame, text="Unmappable", command=self.on_unmappable)
        unmappable_btn.grid(row=1, column=1, padx=2, sticky='ew')

        # Changes text box (bottom, smaller)
        changes_label = ttk.Label(self.center_frame, text="Changes:")
        changes_label.grid(row=2, column=0, sticky='nw', pady=(5, 2))

        self.changes_text = tk.Text(self.center_frame, width=30, height=5, wrap=tk.WORD)
        self.changes_text.configure(height=3)
        self.changes_text.grid(row=2, column=0, sticky='nsew', pady=(25, 0))

        # Add scrollbar to changes text
        changes_scroll = ttk.Scrollbar(self.center_frame, orient='vertical', command=self.changes_text.yview)
        changes_scroll.grid(row=2, column=1, sticky='ns', pady=(25, 0))
        self.changes_text.configure(yscrollcommand=changes_scroll.set)

    def create_3d_plot(self, parent_frame):
        """Create a 3D matplotlib plot embedded in a tkinter frame"""

        # Create figure with square aspect ratio
        fig = Figure(figsize=(5, 5), dpi=100)
        ax = fig.add_subplot(111, projection='3d')

        # Configure axes cleanly (no labels, ticks, planes)
        configure_axes_clean(ax)

        # Create canvas
        canvas = FigureCanvasTkAgg(fig, master=parent_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

        # Store figure and axes for later access
        canvas.figure = fig
        canvas.axes = ax

        try:
            canvas.mpl_connect('button_press_event', self.on_mpl_button_press)
            canvas.mpl_connect('button_release_event', self.on_mpl_button_release)
            canvas.mpl_connect('motion_notify_event', self.on_mpl_motion_notify)
        except Exception:
            pass

        return canvas

    def get_all_canvases(self):
        return [self.upper_left_canvas, self.lower_left_canvas, self.upper_right_canvas, self.lower_right_canvas]

    def get_canvas_from_widget(self, widget):
        try:
            for canvas in self.get_all_canvases():
                if canvas.get_tk_widget() == widget:
                    return canvas
        except Exception:
            pass
        return None

    def sync_all_views_from_axes(self, source_ax):
        """Copy view orientation from one 3D axes to all other axes."""
        try:
            elev = getattr(source_ax, 'elev', None)
            azim = getattr(source_ax, 'azim', None)
            roll = getattr(source_ax, 'roll', None)
            if elev is None or azim is None:
                return
            for canvas in self.get_all_canvases():
                ax = canvas.axes
                if ax is source_ax:
                    continue
                try:
                    if roll is None:
                        ax.view_init(elev=elev, azim=azim)
                    else:
                        ax.view_init(elev=elev, azim=azim, roll=roll)
                except Exception:
                    try:
                        ax.view_init(elev=elev, azim=azim)
                    except Exception:
                        pass
                try:
                    canvas.draw_idle()
                except Exception:
                    pass
        except Exception:
            pass

    def pan_all_axes(self, source_ax, dx_pixels, dy_pixels):
        """Pan all axes together based on pixel delta in the source widget."""
        try:
            xlim = source_ax.get_xlim()
            ylim = source_ax.get_ylim()
            xspan = xlim[1] - xlim[0]
            yspan = ylim[1] - ylim[0]

            width_px = 800.0
            height_px = 800.0
            try:
                fig = source_ax.figure
                if fig and fig.bbox:
                    width_px = max(1.0, float(fig.bbox.width))
                    height_px = max(1.0, float(fig.bbox.height))
            except Exception:
                pass

            shift_x = -dx_pixels * (xspan / width_px)
            shift_y = dy_pixels * (yspan / height_px)

            for canvas in self.get_all_canvases():
                ax = canvas.axes
                try:
                    cur_xlim = ax.get_xlim()
                    cur_ylim = ax.get_ylim()
                    ax.set_xlim([cur_xlim[0] + shift_x, cur_xlim[1] + shift_x])
                    ax.set_ylim([cur_ylim[0] + shift_y, cur_ylim[1] + shift_y])
                    ax.autoscale(False)
                    ax.margins(0)
                    canvas.draw_idle()
                except Exception:
                    pass
        except Exception:
            pass

    def on_mpl_button_press(self, event):
        """Track matplotlib drag-based view rotation starts."""
        try:
            if event.button == 1 and event.inaxes in [c.axes for c in self.get_all_canvases()]:
                self.left_mouse_rotate_active = True
        except Exception:
            pass

    def on_mpl_button_release(self, event):
        """Sync views when drag-based rotation ends."""
        try:
            if event.inaxes in [c.axes for c in self.get_all_canvases()]:
                self.sync_all_views_from_axes(event.inaxes)
        except Exception:
            pass
        self.left_mouse_rotate_active = False

    def on_mpl_motion_notify(self, event):
        """Continuously sync matplotlib 3D view rotation while dragging."""
        try:
            if self.left_mouse_rotate_active and event.inaxes in [c.axes for c in self.get_all_canvases()]:
                self.sync_all_views_from_axes(event.inaxes)
        except Exception:
            pass

    def draw_cube(self, ax, color='cyan', clear=True):
        """Draw a unit cube in 3D"""
        if clear:
            ax.clear()
            ax.set_xlabel('')
            ax.set_ylabel('')
            ax.set_zlabel('')
            ax.set_xticklabels([])
            ax.set_yticklabels([])
            ax.set_zticklabels([])
            ax.grid(False)

        # Define the vertices of a unit cube
        vertices = np.array([
            [0, 0, 0], [1, 0, 0], [1, 1, 0], [0, 1, 0],  # Bottom face
            [0, 0, 1], [1, 0, 1], [1, 1, 1], [0, 1, 1]   # Top face
        ])

        # Define the edges
        edges = [
            [0, 1], [1, 2], [2, 3], [3, 0],  # Bottom face
            [4, 5], [5, 6], [6, 7], [7, 4],  # Top face
            [0, 4], [1, 5], [2, 6], [3, 7]   # Vertical edges
        ]

        # Draw edges
        for edge in edges:
            points = vertices[edge]
            ax.plot3D(*points.T, color=color, linewidth=2)

        # Set equal aspect ratio
        ax.set_xlim([-0.5, 1.5])
        ax.set_ylim([-0.5, 1.5])
        ax.set_zlim([-0.5, 1.5])

    def draw_upper_left_cube(self):
        """Draw cube in upper left window"""
        # Replace cube placeholder with a clean, empty axes (no cube, no axis text)
        ax = self.upper_left_canvas.axes
        ax.clear()
        configure_axes_clean(ax)
        self.upper_left_canvas.draw_idle()

    def draw_lower_left_cube(self):
        """Draw cube in lower left window"""
        # Draw the standard residue in lower-left using parent-atom coloring if available
        try:
            self.draw_lower_left_component()
        except Exception:
            # fallback to a clean empty axes
            ax = self.lower_left_canvas.axes
            ax.clear()
            configure_axes_clean(ax)
            self.lower_left_canvas.draw_idle()

    def draw_upper_right_cube(self):
        """Draw cube in upper right window"""
        self.draw_cube(self.upper_right_canvas.axes, color='blue')
        self.upper_right_canvas.draw()

    def draw_upper_right_component(self):
        """Draw component atom coordinates in upper right window using CPK colors and atom labels"""
        ax = self.upper_right_canvas.axes
        ax.clear()
        ax.set_xlabel('')
        ax.set_ylabel('')
        ax.set_zlabel('')
        ax.set_xticklabels([])
        ax.set_yticklabels([])
        ax.set_zticklabels([])
        ax.grid(False)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_zticks([])
        # Ensure axes are configured consistently (no labels, panes, axis lines)
        configure_axes_clean(ax)

        cif_data = read_component_from_file_or_download(self.mod_nt)
        if not cif_data:
            ax.text(0, 0, 0, f"{self.mod_nt} .cif not found/downloadable", fontsize=10, color='red')
            self.upper_right_canvas.draw()
            return

        coordinates, connections, atom_to_element = get_component_coordinates_and_bonds(cif_data)
        if len(coordinates) == 0:
            ax.text(0, 0, 0, f"No coordinates for {self.mod_nt}", fontsize=10, color='red')
            self.upper_right_canvas.draw()
            return

        # Calculate center before applying rotation
        xyz = np.array(list(coordinates.values()))
        mins = xyz.min(axis=0)
        maxs = xyz.max(axis=0)
        center = (mins + maxs) / 2.0


        # If we have good correspondences to a parent base, align this modified
        # residue to the parent's canonical orientation once (>=3 base atoms).
        try:
            modified_to_parent = getattr(self, 'modified_to_parent_atom', {}).get(self.mod_nt, {})
            parent_name = getattr(self, 'modified_base_to_parent', {}).get(self.mod_nt, None)
            if parent_name and modified_to_parent:
                # read parent CIF to obtain element types so we can exclude hydrogens
                parent_atom_to_element = {}
                try:
                    cif_par = read_component_from_file_or_download(parent_name)
                    if cif_par:
                        _, _, parent_atom_to_element = get_component_coordinates_and_bonds(cif_par)
                except Exception:
                    parent_atom_to_element = {}

                # count mapped parent atoms that appear in the base template, excluding hydrogens
                mapped_pairs = []
                for ma, pa in modified_to_parent.items():
                    try:
                        if not (parent_name in NAbasecoordinates and pa in NAbasecoordinates[parent_name]):
                            continue
                        # skip hydrogen atoms on either side
                        p_el = parent_atom_to_element.get(pa, '')
                        m_el = atom_to_element.get(ma, '')
                        if p_el and str(p_el).strip().upper().startswith('H'):
                            continue
                        if m_el and str(m_el).strip().upper().startswith('H'):
                            continue
                        mapped_pairs.append((pa, ma))
                    except Exception:
                        continue
                if len(mapped_pairs) >= 3:
                    # build lists for besttransformation: template coords then mod coords
                    par_list = [np.array(NAbasecoordinates[parent_name][pa], dtype=float) for pa, ma in mapped_pairs]
                    mod_list = [np.array(coordinates[ma], dtype=float) for pa, ma in mapped_pairs]
                    try:
                        U, new1, mean1, rmsd, sse, mean2 = besttransformation(par_list, mod_list)
                        # apply transform to all modified coords
                        transformed = {}
                        for atom, coord in coordinates.items():
                            arr = np.array(coord)
                            try:
                                c = mean1 + np.dot(U, arr - mean2)
                                c = np.array(c).reshape(3,)
                                transformed[atom] = [float(c[0]), float(c[1]), float(c[2])]
                            except Exception:
                                transformed[atom] = coord
                        # ensure base sits at origin and C1' on -Y
                        transformed = center_and_orient_base_xy(parent_name, transformed)
                        # build signature of mapping pairs used for this alignment
                        try:
                            sig = tuple(sorted([(pa, ma) for pa, ma in mapped_pairs]))
                        except Exception:
                            sig = None

                        # If we've already aligned with the same signature, skip re-applying
                        if getattr(self, 'right_panels_aligned', False) and getattr(self, 'right_alignment_signature', None) == sig:
                            # reuse existing aligned coords
                            coordinates = getattr(self, 'right_aligned_coords', transformed)
                            try:
                                self.changes_text.insert('end', f"Reused existing right-panel alignment for {self.mod_nt} (signature unchanged)\n")
                            except Exception:
                                pass
                        else:
                            # store transformed coords for rendering and mark aligned
                            coordinates = transformed
                            try:
                                self.right_aligned_coords = transformed
                            except Exception:
                                pass
                            self.right_panels_aligned = True
                            try:
                                self.right_alignment_signature = sig
                            except Exception:
                                pass
                            # Log which parent atoms were used for this alignment
                            try:
                                pairs_str = ', '.join([f"{pa}->{ma}" for pa, ma in mapped_pairs])
                                self.changes_text.insert('end', f"Aligned {self.mod_nt} to parent {parent_name} using {len(mapped_pairs)} pairs: {pairs_str}\n")
                            except Exception:
                                pass
                    except Exception:
                        pass
        except Exception:
            pass

        # If a shared aligned coords object exists, prefer it for exact matching
        if getattr(self, 'right_panels_aligned', False) and getattr(self, 'right_aligned_coords', None):
            try:
                coordinates = self.right_aligned_coords
                # recompute center for drawing
                xyz = np.array(list(coordinates.values()))
                mins = xyz.min(axis=0)
                maxs = xyz.max(axis=0)
                center = (mins + maxs) / 2.0
            except Exception:
                pass

        # Apply roll rotation if needed
        if self.roll_angle != 0:
            coordinates = apply_roll_rotation(coordinates, self.roll_angle, center)

        # Center coordinates on mapped atoms when possible so the base/backbone
        # stays central even for modified residues with long side chains.
        try:
            mapped_atoms = []
            try:
                mapped_atoms = list(getattr(self, 'modified_to_parent_atom', {}).get(self.mod_nt, {}).keys())
            except Exception:
                mapped_atoms = []
            if mapped_atoms:
                coordinates = center_coords_on_atoms(coordinates, mapped_atoms)
                center = np.array([0.0, 0.0, 0.0])
            else:
                coordinates = center_coords_to_centroid(coordinates)
                center = np.array([0.0, 0.0, 0.0])
            # recompute center for drawing after centering
            try:
                xyz = np.array(list(coordinates.values()))
                mins = xyz.min(axis=0)
                maxs = xyz.max(axis=0)
            except Exception:
                pass
        except Exception:
            pass

        # Recalculate bounds after rotation/centering and compute span for sizing
        xyz = np.array(list(coordinates.values()))
        mins = xyz.min(axis=0)
        maxs = xyz.max(axis=0)
        span = float(np.max(maxs - mins))
        if span < 1.0:
            span = 1.0

        # Draw bonds using original coordinates (no atom balls), thicker linewidth for visibility
        draw_component_bonds_and_atoms(ax, coordinates, connections, atom_to_element,
                          bond_linewidth=6.0, atom_labelsize=8, draw_atom_markers=False)

        # Compute viewing box with reasonable padding so molecule is visible and not clipped
        half = span / 2.0 + 1.0

        ax.set_xlim(center[0] - half, center[0] + half)
        ax.set_ylim(center[1] - half, center[1] + half)
        ax.set_zlim(center[2] - half, center[2] + half)

        # Set equal aspect ratio for x, y, z (1:1:1 scale)
        ax.set_box_aspect([1, 1, 1])

        # Disable matplotlib's automatic margin adjustment which can cause clipping
        ax.margins(0)
        ax.autoscale(False)

        # If we auto-aligned the right panels, make the view match the left
        try:
            if getattr(self, 'right_panels_aligned', False):
                try:
                    configure_axes_clean(ax)
                except Exception:
                    pass
                try:
                    ax.view_init(elev=90, azim=270)
                except Exception:
                    pass
        except Exception:
            pass

        self.upper_right_canvas.draw_idle()

    def draw_upper_left_component(self):
        """Draw the standard residue in upper-left using default CPK colors and shared roll."""
        if not getattr(self, 'std_base', None):
            return

        ax = self.upper_left_canvas.axes
        ax.clear()
        configure_axes_clean(ax)

        cif_data_std = read_component_from_file_or_download(self.std_base)
        if not cif_data_std:
            self.upper_left_canvas.draw_idle()
            return

        coords, conns, atom_to_element = get_component_coordinates_and_bonds(cif_data_std)

        # Orient parent coordinates into the standard/base plane using known template
        try:
            coords = orient_component_to_standard(self.std_base, coords)
        except Exception:
            pass

        # Translate/rotate so base atom sits at origin and C1' points -Y
        try:
            coords = center_and_orient_base_xy(self.std_base, coords)
        except Exception:
            pass

        # compute center and apply shared roll
        xyz = np.array(list(coords.values()))
        mins = xyz.min(axis=0)
        maxs = xyz.max(axis=0)
        center = (mins + maxs) / 2.0
        if self.roll_angle != 0:
            coords = apply_roll_rotation(coords, self.roll_angle, center)

        # Center coordinates at centroid so the molecule is centered in the view
        try:
            coords = center_coords_to_centroid(coords)
            # recompute bounds after centering
            xyz = np.array(list(coords.values()))
            mins = xyz.min(axis=0)
            maxs = xyz.max(axis=0)
            center = (mins + maxs) / 2.0
        except Exception:
            pass

        draw_component_bonds_and_atoms(ax, coords, conns, atom_to_element, bond_linewidth=6.0, atom_labelsize=8, draw_atom_markers=False)

        span = float(np.max(maxs - mins))
        if span < 1.0:
            span = 1.0
        padding = max(0.2, span * 0.12)
        half = span / 2.0 + padding

        ax.set_xlim(center[0] - half, center[0] + half)
        ax.set_ylim(center[1] - half, center[1] + half)
        ax.set_zlim(center[2] - half, center[2] + half)
        ax.set_box_aspect([1,1,1])

        # Hide axes and ticks for left panels (clean view)
        try:
            configure_axes_clean(ax)
        except Exception:
            pass

        ax.margins(0)
        ax.autoscale(False)
        # Set view to look straight down +Z and keep consistent for left panes
        try:
            ax.view_init(elev=90, azim=270)
            # ensure axes directions: positive X to right, positive Y up
            try:
                if getattr(ax, 'xaxis_inverted', None) and ax.xaxis_inverted():
                    ax.invert_xaxis()
            except Exception:
                pass
            try:
                if getattr(ax, 'yaxis_inverted', None) and ax.yaxis_inverted():
                    ax.invert_yaxis()
            except Exception:
                pass
        except Exception:
            pass
        self.upper_left_canvas.draw_idle()

    def draw_lower_right_component(self):
        """Draw the standard residue in lower-right using the alternate color scheme."""
        # Draw the modified residue (self.mod_nt) colored by any mapping to the
        # standard parent atoms. Atoms that map to a parent atom receive the
        # parent's atom color (from parent_atom_to_color); atoms that are not
        # mapped will result in thinner bonds.
        ax = self.lower_right_canvas.axes
        ax.clear()
        configure_axes_clean(ax)

        cif_data_mod = read_component_from_file_or_download(self.mod_nt)
        if not cif_data_mod:
            self.lower_right_canvas.draw_idle()
            return

        coords, conns, atom_to_element = get_component_coordinates_and_bonds(cif_data_mod)

        # compute center and potentially auto-align to parent orientation once
        xyz = np.array(list(coords.values()))
        mins = xyz.min(axis=0)
        maxs = xyz.max(axis=0)
        center = (mins + maxs) / 2.0

        try:
            modified_to_parent = getattr(self, 'modified_to_parent_atom', {}).get(self.mod_nt, {})
            parent_name = getattr(self, 'modified_base_to_parent', {}).get(self.mod_nt, None)
            if parent_name and modified_to_parent:
                # read parent CIF to obtain element types so we can exclude hydrogens
                parent_atom_to_element = {}
                try:
                    cif_par = read_component_from_file_or_download(parent_name)
                    if cif_par:
                        _, _, parent_atom_to_element = get_component_coordinates_and_bonds(cif_par)
                except Exception:
                    parent_atom_to_element = {}

                mapped_pairs = []
                for ma, pa in modified_to_parent.items():
                    try:
                        if not (parent_name in NAbasecoordinates and pa in NAbasecoordinates[parent_name]):
                            continue
                        p_el = parent_atom_to_element.get(pa, '')
                        m_el = atom_to_element.get(ma, '')
                        if p_el and str(p_el).strip().upper().startswith('H'):
                            continue
                        if m_el and str(m_el).strip().upper().startswith('H'):
                            continue
                        mapped_pairs.append((pa, ma))
                    except Exception:
                        continue
                if len(mapped_pairs) >= 3:
                    par_list = [np.array(NAbasecoordinates[parent_name][pa], dtype=float) for pa, ma in mapped_pairs]
                    mod_list = [np.array(coords[ma], dtype=float) for pa, ma in mapped_pairs]
                    try:
                        U, new1, mean1, rmsd, sse, mean2 = besttransformation(par_list, mod_list)
                        transformed = {}
                        for atom, coord in coords.items():
                            arr = np.array(coord)
                            try:
                                c = mean1 + np.dot(U, arr - mean2)
                                c = np.array(c).reshape(3,)
                                transformed[atom] = [float(c[0]), float(c[1]), float(c[2])]
                            except Exception:
                                transformed[atom] = coord
                        transformed = center_and_orient_base_xy(parent_name, transformed)
                        # compute signature for mapped pairs
                        try:
                            sig = tuple(sorted([(pa, ma) for pa, ma in mapped_pairs]))
                        except Exception:
                            sig = None

                        if getattr(self, 'right_panels_aligned', False) and getattr(self, 'right_alignment_signature', None) == sig:
                            # already aligned with same mapping; reuse existing
                            coords = getattr(self, 'right_aligned_coords', transformed)
                            try:
                                self.changes_text.insert('end', f"Reused existing right-panel alignment for {self.mod_nt} (signature unchanged)\n")
                            except Exception:
                                pass
                        else:
                            coords = transformed
                            try:
                                self.right_aligned_coords = transformed
                            except Exception:
                                pass
                            self.right_panels_aligned = True
                            try:
                                self.right_alignment_signature = sig
                            except Exception:
                                pass
                            try:
                                pairs_str = ', '.join([f"{pa}->{ma}" for pa, ma in mapped_pairs])
                                self.changes_text.insert('end', f"Aligned {self.mod_nt} to parent {parent_name} using {len(mapped_pairs)} pairs: {pairs_str}\n")
                            except Exception:
                                pass
                        # recompute center after transform
                        xyz = np.array(list(coords.values()))
                        mins = xyz.min(axis=0)
                        maxs = xyz.max(axis=0)
                        center = (mins + maxs) / 2.0
                    except Exception:
                        pass
        except Exception:
            pass

        if self.roll_angle != 0:
            coords = apply_roll_rotation(coords, self.roll_angle, center)

        # Use in-memory mappings stored on the GUI (updated by Process)
        parent_to_modified_atom = getattr(self, 'parent_to_modified_atom', {})
        modified_to_parent_atom = getattr(self, 'modified_to_parent_atom', {})
        modified_base_to_parent = getattr(self, 'modified_base_to_parent', {})

        # Determine parent for this modified nucleotide
        parent = modified_base_to_parent.get(self.mod_nt, None)
        color_map_mod = {}
        if parent:
            # handle deoxy naming variants (DA vs A etc.)
            pkey = parent
            if pkey not in parent_atom_to_color and len(pkey) > 1 and pkey.startswith('D'):
                maybe = pkey[1:]
                if maybe in parent_atom_to_color:
                    pkey = maybe

            # For each modified atom that maps to a parent atom, assign parent's color
            mod_to_parent = modified_to_parent_atom.get(self.mod_nt, {})
            for mod_atom, parent_atom in mod_to_parent.items():
                if parent in parent_atom_to_color and parent_atom in parent_atom_to_color[parent]:
                    color_map_mod[mod_atom] = parent_atom_to_color[parent][parent_atom]
                else:
                    # Fallback to element-based coloring
                    color_map_mod[mod_atom] = element_to_cpk_color(atom_to_element.get(mod_atom, 'C'))

        # Center coordinates on mapped atoms when possible so the base/backbone
        # stays central even for modified residues with long side chains.
        try:
            mapped_atoms = []
            try:
                mapped_atoms = list(self.modified_to_parent_atom.get(self.mod_nt, {}).keys())
            except Exception:
                mapped_atoms = []
            if mapped_atoms:
                coords = center_coords_on_atoms(coords, mapped_atoms)
                center = np.array([0.0, 0.0, 0.0])
            else:
                coords = center_coords_to_centroid(coords)
                center = np.array([0.0, 0.0, 0.0])
            xyz = np.array(list(coords.values()))
            mins = xyz.min(axis=0)
            maxs = xyz.max(axis=0)
        except Exception:
            pass

        # draw using mapping-aware colors. Use the same linewidth as other plots
        # so mapped bonds appear with the same prominence.
        draw_component_bonds_and_atoms(ax, coords, conns, atom_to_element, bond_linewidth=6.0, atom_labelsize=8, draw_atom_markers=False, colorfn=lambda e: 'black', color_map=color_map_mod)

        span = float(np.max(maxs - mins))
        if span < 1.0:
            span = 1.0
        padding = max(0.2, span * 0.12)
        half = span / 2.0 + padding

        ax.set_xlim(center[0] - half, center[0] + half)
        ax.set_ylim(center[1] - half, center[1] + half)
        ax.set_zlim(center[2] - half, center[2] + half)
        ax.set_box_aspect([1,1,1])
        ax.margins(0)
        ax.autoscale(False)

        # If we auto-aligned the right panels, make the view match the left
        try:
            if getattr(self, 'right_panels_aligned', False):
                try:
                    configure_axes_clean(ax)
                except Exception:
                    pass
                try:
                    ax.view_init(elev=90, azim=270)
                except Exception:
                    pass
        except Exception:
            pass

        self.lower_right_canvas.draw_idle()

    def draw_lower_left_component(self):
        """Draw the standard residue in lower-left using parent_atom_to_color mapping."""
        if not getattr(self, 'std_base', None):
            ax = self.lower_left_canvas.axes
            ax.clear()
            configure_axes_clean(ax)
            self.lower_left_canvas.draw_idle()
            return

        parent = self.std_base
        # prefer the single-letter parent if mapping uses that
        if parent not in parent_atom_to_color and len(parent) > 1 and parent.startswith('D'):
            # try removing leading D for deoxy forms
            maybe = parent[1:]
            if maybe in parent_atom_to_color:
                parent = maybe

        color_map = parent_atom_to_color.get(parent, {})

        ax = self.lower_left_canvas.axes
        ax.clear()
        configure_axes_clean(ax)

        cif_data_std = read_component_from_file_or_download(self.std_base)
        if not cif_data_std:
            self.lower_left_canvas.draw_idle()
            return


        coords, conns, atom_to_element = get_component_coordinates_and_bonds(cif_data_std)

        # Orient parent coordinates into the standard/base plane using known template
        try:
            coords = orient_component_to_standard(self.std_base, coords)
        except Exception:
            pass

        # Translate/rotate so base atom sits at origin and C1' points -Y
        try:
            coords = center_and_orient_base_xy(self.std_base, coords)
        except Exception:
            pass

        # compute center and apply shared roll
        xyz = np.array(list(coords.values()))
        mins = xyz.min(axis=0)
        maxs = xyz.max(axis=0)
        center = (mins + maxs) / 2.0
        if self.roll_angle != 0:
            coords = apply_roll_rotation(coords, self.roll_angle, center)

        # Center coordinates at centroid so the molecule is centered in the view
        try:
            coords = center_coords_to_centroid(coords)
            xyz = np.array(list(coords.values()))
            mins = xyz.min(axis=0)
            maxs = xyz.max(axis=0)
            center = (mins + maxs) / 2.0
        except Exception:
            pass

        # draw using atom-name-based colors from parent_atom_to_color
        draw_component_bonds_and_atoms(ax, coords, conns, atom_to_element, bond_linewidth=6.0, atom_labelsize=8, draw_atom_markers=False, colorfn=None, atom_markersize=50, color_map=color_map)

        span = float(np.max(maxs - mins))
        if span < 1.0:
            span = 1.0
        padding = max(0.2, span * 0.12)
        half = span / 2.0 + padding

        ax.set_xlim(center[0] - half, center[0] + half)
        ax.set_ylim(center[1] - half, center[1] + half)
        ax.set_zlim(center[2] - half, center[2] + half)
        ax.set_box_aspect([1,1,1])

        # Hide axes and ticks for left panels
        try:
            configure_axes_clean(ax)
        except Exception:
            pass

        ax.margins(0)
        ax.autoscale(False)
        self.lower_left_canvas.draw_idle()
        # Set view to look straight down +Z and keep consistent for left panes
        try:
            # rotate lower-left view 90 degrees counterclockwise relative to default
            ax.view_init(elev=90, azim=270)
            try:
                if getattr(ax, 'xaxis_inverted', None) and ax.xaxis_inverted():
                    ax.invert_xaxis()
            except Exception:
                pass
            try:
                if getattr(ax, 'yaxis_inverted', None) and ax.yaxis_inverted():
                    ax.invert_yaxis()
            except Exception:
                pass
        except Exception:
            pass

    def draw_lower_right_cube(self):
        """Draw cube in lower right window"""
        # Replace placeholder cube with either a clean empty axes or the standard residue
        try:
            self.draw_lower_right_component()
        except Exception:
            ax = self.lower_right_canvas.axes
            ax.clear()
            configure_axes_clean(ax)
            self.lower_right_canvas.draw_idle()

    def redraw_right_windows(self):
        """Redraw the two graphics windows on the right"""
        # Only draw the placeholder cube in upper-right if a component has not yet been drawn
        if not getattr(self, 'upper_right_drawn', False):
            self.draw_upper_right_cube()
        # ensure lower-right shows standard residue (or clean axes)
        try:
            self.draw_lower_right_component()
        except Exception:
            self.draw_lower_right_cube()

    def on_process(self):
        """Handle Process button click"""
        # Get text from mappings text box
        text = self.mappings_text.get('1.0', 'end-1c')

        # Parse lines robustly: support tab or space separators, allow blank mappings,
        # and capture unmatched modified-atom-only lines. If an atom appears with
        # a mapping on any line, that mapping takes precedence over any blank line.
        raw_lines = [l.rstrip('\n') for l in text.splitlines()]
        lines = [ln for ln in raw_lines if ln.strip() != '']

        new_parent = None
        parsed_mappings = {}
        unmatched_mod_user = set()
        explicitly_mapped_mod_user = set()

        if len(lines) > 0:
            # first non-empty line contains parent (and optionally modified)
            first_line = lines[0]
            first_fields = first_line.split('\t')
            if len(first_fields) >= 1 and first_fields[0].strip() != '':
                candidate_parent = first_fields[0].strip().upper()
                if candidate_parent not in ['??', 'STD']:
                    new_parent = candidate_parent
            else:
                # fallback to whitespace split
                toks = first_line.split()
                if len(toks) >= 1:
                    candidate_parent = toks[0].strip().upper()
                    if candidate_parent not in ['??', 'STD']:
                        new_parent = candidate_parent

        # helper: parent atom name list if parent known — prefer CIF atoms when available
        parent_atom_list = []
        if new_parent:
            try:
                cif_par_tmp = read_component_from_file_or_download(new_parent)
                if cif_par_tmp:
                    par_coords_tmp, _, _ = get_component_coordinates_and_bonds(cif_par_tmp)
                    parent_atom_list = list(par_coords_tmp.keys())
                elif new_parent in NAbasecoordinates:
                    parent_atom_list = list(NAbasecoordinates[new_parent].keys())
            except Exception:
                if new_parent in NAbasecoordinates:
                    parent_atom_list = list(NAbasecoordinates[new_parent].keys())

        # parse remaining lines
        for ln in lines[1:]:
            # try tab-separated first to respect generated format
            fields = ln.split('\t')
            if len(fields) >= 2:
                p = fields[0].strip()
                m = fields[1].strip()
                if p == '' and m != '':
                    # user supplied an unmatched modified atom
                    if m not in explicitly_mapped_mod_user:
                        unmatched_mod_user.add(m)
                    continue
                if p != '':
                    # Only accept parent atoms known for this parent component.
                    # This prevents stale template-only atoms (e.g., H9/H9')
                    # from persisting in the central panel for DG.
                    if parent_atom_list and p not in parent_atom_list:
                        continue
                    if m != '':
                        # explicit mapping; always override
                        parsed_mappings[p] = m
                        explicitly_mapped_mod_user.add(m)
                        if m in unmatched_mod_user:
                            unmatched_mod_user.discard(m)
                    else:
                        # explicit blank mapping; only record if no mapping present
                        if p not in parsed_mappings:
                            parsed_mappings[p] = ''
                    continue

            # fallback: whitespace split
            toks = ln.split()
            if len(toks) >= 2:
                p = toks[0].strip()
                m = toks[1].strip()
                if parent_atom_list and p not in parent_atom_list:
                    continue
                parsed_mappings[p] = m
                if m:
                    explicitly_mapped_mod_user.add(m)
                    if m in unmatched_mod_user:
                        unmatched_mod_user.discard(m)
            elif len(toks) == 1:
                token = toks[0].strip()
                # treat single token as parent atom (blank mapping) if it appears
                # in expected parent atoms, otherwise treat as unmatched modified atom
                if token in parent_atom_list:
                    if token not in parsed_mappings:
                        parsed_mappings[token] = ''
                else:
                    if token not in explicitly_mapped_mod_user:
                        unmatched_mod_user.add(token)

        # Update in-memory mappings for this modified residue (no files written)
        if not hasattr(self, 'parent_to_modified_atom'):
            self.parent_to_modified_atom = {}
        if not hasattr(self, 'modified_to_parent_atom'):
            self.modified_to_parent_atom = {}
        if not hasattr(self, 'modified_base_to_parent'):
            self.modified_base_to_parent = {}

        # Reset mappings for this modified residue and populate with parsed values.
        # Keep blanks for display, but only non-empty values are real extension seeds.
        extended_mappings = dict(parsed_mappings)
        parent_name = None
        if new_parent:
            parent_name = new_parent
        elif getattr(self, 'std_base', None):
            parent_name = self.std_base
        elif self.modified_base_to_parent.get(self.mod_nt):
            parent_name = self.modified_base_to_parent.get(self.mod_nt)

        if parent_name:
            try:
                cif_mod = read_component_from_file_or_download(self.mod_nt)
                cif_par = read_component_from_file_or_download(parent_name)
                if cif_mod and cif_par:
                    # If neither side has C1' mapped, but both CIFs contain C1', seed it
                    try:
                        coords_mod, _, _ = get_component_coordinates_and_bonds(cif_mod)
                        coords_par, _, _ = get_component_coordinates_and_bonds(cif_par)
                        if not extended_mappings.get("C1'", ""):
                            if "C1'" in coords_mod and "C1'" in coords_par:
                                extended_mappings["C1'"] = "C1'"
                                self.changes_text.insert('end', "Seeded mapping C1' -> C1'\n")
                    except Exception:
                        pass

                    # Attempt to extend mappings if we now have at least one non-empty seed
                    nonempty_seeds = {k: v for k, v in extended_mappings.items() if k and v}
                    if nonempty_seeds:
                        try:
                            ext = extend_mappings_by_neighbors(cif_mod, cif_par, nonempty_seeds)
                            # merge extended mapping but preserve any explicit user mappings
                            if ext:
                                merged = dict(extended_mappings)
                                for k, v in ext.items():
                                    # fill only currently blank/unset entries
                                    if not merged.get(k, ''):
                                        merged[k] = v
                                # explicit non-empty user mappings always win
                                for k, v in extended_mappings.items():
                                    if v:
                                        merged[k] = v
                                extended_mappings = merged
                                self.changes_text.insert('end', f"Auto-extended mappings -> {len(extended_mappings)} pairs (merged)\n")
                        except Exception:
                            pass
            except Exception:
                pass

        if extended_mappings:
            # If mappings changed for this modified residue, clear any previous alignment
            old_map = self.parent_to_modified_atom.get(self.mod_nt, {})
            if old_map != extended_mappings:
                try:
                    self.right_panels_aligned = False
                    if hasattr(self, 'right_aligned_coords'):
                        delattr = False
                        try:
                            del self.right_aligned_coords
                        except Exception:
                            pass
                        try:
                            del self.right_alignment_signature
                        except Exception:
                            pass
                except Exception:
                    pass
            # Before committing mappings, attempt to fill hydrogens using
            # specialized geometry heuristics. First build a temporary reverse map
            # and call the filler which updates the dicts in-place.
            try:
                temp_rev = {v: k for k, v in extended_mappings.items() if v}
                if parent_name:
                    try:
                        cif_mod = read_component_from_file_or_download(self.mod_nt)
                        cif_par = read_component_from_file_or_download(parent_name)
                    except Exception:
                        cif_mod = None
                        cif_par = None
                else:
                    cif_mod = None
                    cif_par = None

                try:
                    parent_to_mod_outer = {self.mod_nt: dict(extended_mappings)}
                    mod_to_parent_outer = {self.mod_nt: dict(temp_rev)}
                    added = fill_hydrogen_mappings(self, self.mod_nt, parent_name, parent_to_mod_outer, mod_to_parent_outer, cif_mod=cif_mod, cif_par=cif_par)
                    # IMPORTANT: pull back modifications from filler into extended mappings
                    extended_mappings = dict(parent_to_mod_outer.get(self.mod_nt, extended_mappings))
                    temp_rev = dict(mod_to_parent_outer.get(self.mod_nt, temp_rev))
                    if added:
                        try:
                            self.changes_text.insert('end', f"Added {added} hydrogen mappings for {self.mod_nt}\n")
                        except Exception:
                            pass
                except Exception:
                    pass

                # commit mappings (including any hydrogen additions)
                self.parent_to_modified_atom[self.mod_nt] = dict(extended_mappings)
                # rebuild reverse mapping after hydrogen filling
                rev = {}
                for pa, ma in extended_mappings.items():
                    if ma:
                        rev[ma] = pa
                self.modified_to_parent_atom[self.mod_nt] = rev
            except Exception:
                # fallback: commit simple mapping
                self.parent_to_modified_atom[self.mod_nt] = dict(extended_mappings)
                rev = {v: k for k, v in extended_mappings.items()}
                self.modified_to_parent_atom[self.mod_nt] = rev

        # If a new parent was specified, update it
        if new_parent:
            self.std_base = new_parent
            self.modified_base_to_parent[self.mod_nt] = new_parent

        # Normalize and write back the mappings_text using tabs. Show
        # all parent atoms (base first, then ribose/phosphate) and then
        # list unmatched modified atoms.
        out_lines = []
        if getattr(self, 'std_base', None):
            out_lines.append(f"{self.std_base}\t{self.mod_nt}")
            parent = self.std_base
        else:
            out_lines.append(f"??\t{self.mod_nt}")
            parent = None

        # determine canonical parent atom ordering: prefer atoms from the parent CIF
        # when available (so we only list atoms actually present), otherwise fall
        # back to the template `NAbasecoordinates`.
        parent_atoms_all = []
        par_coords = None
        try:
            if parent:
                cif_par = read_component_from_file_or_download(parent)
                if cif_par:
                    par_coords, _, _ = get_component_coordinates_and_bonds(cif_par)
                    parent_atoms_all = list(par_coords.keys())
                elif parent in NAbasecoordinates:
                    parent_atoms_all = list(NAbasecoordinates[parent].keys())
        except Exception:
            parent_atoms_all = []

        rp = list(ribose_full) + list(phosphate_full)

        # Drop parent keys that are not present in the selected parent atom list
        # (or optional rp add-ons above). This removes stale entries like H9/H9'.
        try:
            allowed_parent_atoms = set(parent_atoms_all)
            cur_map = dict(self.parent_to_modified_atom.get(self.mod_nt, {}))
            cleaned = {pa: ma for pa, ma in cur_map.items() if pa in allowed_parent_atoms}
            self.parent_to_modified_atom[self.mod_nt] = cleaned
            rev_clean = {}
            for pa, ma in cleaned.items():
                if ma:
                    rev_clean[ma] = pa
            self.modified_to_parent_atom[self.mod_nt] = rev_clean
        except Exception:
            pass

        # Remove O2' for DNA parents (DA/DC/DG/DT)
        try:
            if parent and parent in ['DA','DC','DG','DT'] and "O2'" in parent_atoms_all:
                parent_atoms_all = [a for a in parent_atoms_all if a != "O2'"]
        except Exception:
            pass

        # also ensure parent atoms explicitly mentioned by the user are preserved
        # (but avoid pulling in template-only legacy keys like H9/H9' from old maps).
        try:
            for a in parsed_mappings.keys():
                if a and a not in parent_atoms_all:
                    parent_atoms_all.append(a)
        except Exception:
            pass

        # Split into base vs ribose/phosphate and sort alphabetically within groups
        base_atoms = sorted([a for a in parent_atoms_all if a not in rp])
        rp_atoms = sorted([a for a in parent_atoms_all if a in rp])

        pmap = self.parent_to_modified_atom.get(self.mod_nt, {})
        for parent_atom in base_atoms + rp_atoms:
            out_lines.append(f"{parent_atom}\t{pmap.get(parent_atom, '')}")

        # now list unmatched modified atoms (those not mapped to any parent)
        try:
            cif_mod = read_component_from_file_or_download(self.mod_nt)
            if cif_mod:
                mod_coords, _, _ = get_component_coordinates_and_bonds(cif_mod)
                mapped_mods = set(self.modified_to_parent_atom.get(self.mod_nt, {}).keys())
                unmatched = sorted([m for m in mod_coords.keys() if m not in mapped_mods])
            else:
                unmatched = []
        except Exception:
            unmatched = []

        # include any unmatched mod atoms the user typed explicitly earlier
        try:
            unmatched_user = sorted(list(unmatched_mod_user))
            for m in unmatched_user:
                if m not in unmatched and m not in mapped_mods:
                    unmatched.append(m)
        except Exception:
            pass

        for ma in unmatched:
            out_lines.append(f"\t{ma}")

        try:
            self.mappings_text.delete('1.0', 'end')
            self.mappings_text.insert('1.0', '\n'.join(out_lines))
        except Exception:
            pass

        # Log to changes box
        self.changes_text.insert('end', f"Processed mappings for {self.mod_nt}: parent={self.std_base} mappings={len(pmap)}\n")

        # Redraw views: update left windows to reflect any parent change and
        # update lower-right to reflect new correspondences
        try:
            if getattr(self, 'std_base', None):
                self.draw_upper_left_component()
                self.draw_lower_left_component()
        except Exception:
            pass

        try:
            # Always redraw upper-right so alignment may be recomputed when Process is clicked
            self.draw_upper_right_component()
            self.upper_right_drawn = True
        except Exception:
            pass

        try:
            self.draw_lower_right_component()
        except Exception:
            self.draw_lower_right_cube()

        print(f"Process clicked. parent={self.std_base} mappings={len(pmap)}")

    def on_save(self):
        """Handle Save button click"""
        # Get text from mappings text box
        self.edit_text = self.mappings_text.get('1.0', 'end-1c')

        # Re-parse mappings_text robustly to ensure what's shown in the central box
        # is committed to the in-memory mapping tables before writing the file.
        try:
            raw_lines = [l.rstrip('\n') for l in self.edit_text.splitlines()]
            lines = [ln for ln in raw_lines if ln.strip() != '']
            new_parent = None
            parsed = {}
            unmatched_user = set()

            if len(lines) > 0:
                first_line = lines[0]
                first_fields = first_line.split('\t')
                if len(first_fields) >= 1 and first_fields[0].strip() != '':
                    candidate_parent = first_fields[0].strip().upper()
                    if candidate_parent not in ['??', 'STD']:
                        new_parent = candidate_parent
                else:
                    toks = first_line.split()
                    if len(toks) >= 1:
                        candidate_parent = toks[0].strip().upper()
                        if candidate_parent not in ['??', 'STD']:
                            new_parent = candidate_parent

            parent_atom_list = []
            if new_parent:
                try:
                    cif_par_tmp = read_component_from_file_or_download(new_parent)
                    if cif_par_tmp:
                        par_coords_tmp, _, _ = get_component_coordinates_and_bonds(cif_par_tmp)
                        parent_atom_list = list(par_coords_tmp.keys())
                    elif new_parent in NAbasecoordinates:
                        parent_atom_list = list(NAbasecoordinates[new_parent].keys())
                except Exception:
                    if new_parent in NAbasecoordinates:
                        parent_atom_list = list(NAbasecoordinates[new_parent].keys())

            for ln in lines[1:]:
                fields = ln.split('\t')
                if len(fields) >= 2:
                    p = fields[0].strip()
                    m = fields[1].strip()
                    if p == '' and m != '':
                        unmatched_user.add(m)
                        continue
                    if p != '':
                        if parent_atom_list and p not in parent_atom_list:
                            continue
                        if m != '':
                            parsed[p] = m
                        else:
                            if p not in parsed:
                                parsed[p] = ''
                        continue
                toks = ln.split()
                if len(toks) >= 2:
                    p = toks[0].strip()
                    m = toks[1].strip()
                    if parent_atom_list and p not in parent_atom_list:
                        continue
                    parsed[p] = m
                elif len(toks) == 1:
                    token = toks[0].strip()
                    if token in parent_atom_list:
                        if token not in parsed:
                            parsed[token] = ''
                    else:
                        unmatched_user.add(token)

            # Commit parsed mappings to GUI state
            try:
                if new_parent:
                    self.modified_base_to_parent[self.mod_nt] = new_parent
                    self.std_base = new_parent
                self.parent_to_modified_atom[self.mod_nt] = dict(parsed)
                rev = {}
                for pa, ma in parsed.items():
                    if ma:
                        rev[ma] = pa
                # include existing reverse entries for any mapped hydrogens added earlier
                for ma, pa in self.modified_to_parent_atom.get(self.mod_nt, {}).items():
                    if ma not in rev and pa:
                        rev[ma] = pa
                self.modified_to_parent_atom[self.mod_nt] = rev
            except Exception:
                pass

        except Exception:
            pass

        # Update new_atom_mappings.txt for this modified nucleotide
        try:
            write_new_atom_mappings_for_mod(self.mod_nt, self.parent_to_modified_atom, self.modified_to_parent_atom, self.modified_base_to_parent)
            self.changes_text.insert('end', f"Wrote new_atom_mappings.txt for {self.mod_nt}\n")
        except Exception:
            pass

        # Store return data
        self.return_action = 'save'
        self.return_data = self.edit_text

        print(f"Save clicked. Returning data: {self.return_data}")

        # Advance batch or close single-window session
        try:
            self.master.quit()
        except Exception:
            pass
        if not getattr(self, 'batch_mode', False):
            try:
                self.master.destroy()
            except Exception:
                pass

    def on_skip(self):
        """Skip current CCD without writing new_atom_mappings.txt"""
        self.return_action = 'skip'
        self.return_data = None

        print("Skip clicked. Moving to next CCD.")

        try:
            self.master.quit()
        except Exception:
            pass
        if not getattr(self, 'batch_mode', False):
            try:
                self.master.destroy()
            except Exception:
                pass

    def on_unmappable(self):
        """Mark current CCD as unmappable and move to next CCD."""
        try:
            # Remove any parent mapping so write helper emits unmappable marker line.
            try:
                if hasattr(self, 'modified_base_to_parent') and self.mod_nt in self.modified_base_to_parent:
                    del self.modified_base_to_parent[self.mod_nt]
            except Exception:
                pass

            write_new_atom_mappings_for_mod(self.mod_nt, self.parent_to_modified_atom, self.modified_to_parent_atom, self.modified_base_to_parent)
            try:
                self.changes_text.insert('end', f"Marked {self.mod_nt} as unmappable in new_atom_mappings.txt\n")
            except Exception:
                pass
        except Exception:
            pass

        self.return_action = 'unmappable'
        self.return_data = None

        print(f"Unmappable clicked for {self.mod_nt}. Moving to next CCD.")

        try:
            self.master.quit()
        except Exception:
            pass
        if not getattr(self, 'batch_mode', False):
            try:
                self.master.destroy()
            except Exception:
                pass

    def on_window_close(self):
        """Quit the entire run when the window X button is clicked."""
        self.return_action = 'quit-all'
        self.return_data = None
        print("Window close requested. Quitting entirely.")
        try:
            self.master.quit()
        except Exception:
            pass
        try:
            self.master.destroy()
        except Exception:
            pass

    def on_discard(self):
        """Backward-compatible alias for Skip."""
        self.on_skip()


def main(mod_nt):
    """Run the GUI for one CCD and return (action, data)."""
    root = tk.Tk()
    app = MapModifiedGUI(root, mod_nt)
    root.mainloop()

    action = getattr(app, 'return_action', None)
    data = getattr(app, 'return_data', None)
    return action, data


def parse_ccd_inputs(input_arg):
    """Parse CCD IDs from a literal string or from a file path.
    Accepts comma/tab/newline/space separated values (mixed separators allowed).
    """
    content = input_arg
    if os.path.isfile(input_arg):
        try:
            with open(input_arg, 'r', encoding='utf-8') as f:
                content = f.read()
        except Exception:
            content = input_arg

    parts = re.split(r'[\s,\t\r\n]+', content)
    out = []
    seen = set()
    for p in parts:
        token = p.strip().upper()
        if not token:
            continue
        if token in seen:
            continue
        seen.add(token)
        out.append(token)
    return out


def get_already_handled_ccds():
    """Return set of CCD IDs already present in atom_mappings/new_atom_mappings."""
    handled = set()

    try:
        p2m, m2p, m2parent, not_map = read_atom_mappings("atom_mappings.txt")
        handled.update([k.upper() for k in m2parent.keys()])
        handled.update([k.upper() for k in p2m.keys()])
        handled.update([k.upper() for k in not_map])
    except Exception:
        pass

    try:
        p2m_n, m2p_n, m2parent_n, not_map_n = read_new_atom_mappings("new_atom_mappings.txt")
        handled.update([k.upper() for k in m2parent_n.keys()])
        handled.update([k.upper() for k in p2m_n.keys()])
        handled.update([k.upper() for k in not_map_n])
    except Exception:
        pass

    return handled


if __name__ == "__main__":
    # Get input argument (single CCD ID or path to file containing CCD IDs)
    if len(sys.argv) < 2:
        print("Usage: python map_modified.py <CCD_or_input_file>")
        sys.exit(1)

    input_arg = sys.argv[1]
    ccd_list = parse_ccd_inputs(input_arg)
    if not ccd_list:
        print("No CCD identifiers found in input.")
        sys.exit(0)

    handled = get_already_handled_ccds()
    remaining = [ccd for ccd in ccd_list if ccd.upper() not in handled]

    if not remaining:
        print("All CCD identifiers in input are already present in atom_mappings.txt or new_atom_mappings.txt.")
        sys.exit(0)

    if len(remaining) == 1:
        mod_nt = remaining[0]
        print(f"Processing {mod_nt} ...")
        action, result = main(mod_nt)
        if result is not None:
            print(f"Returned data for {mod_nt}:\n{result}")
    else:
        root = tk.Tk()
        app = MapModifiedGUI(root, remaining[0], batch_mode=True)
        index = 0

        while index < len(remaining):
            mod_nt = remaining[index]
            print(f"Processing {mod_nt} ...")
            root.mainloop()
            action = getattr(app, 'return_action', None)
            result = getattr(app, 'return_data', None)

            if action == 'quit-all':
                print("Quit requested. Stopping batch processing.")
                break
            elif action == 'save':
                if result is not None:
                    print(f"Returned data for {mod_nt}:\n{result}")
            elif action == 'unmappable':
                print(f"Marked {mod_nt} as unmappable.")
            elif action == 'skip':
                print(f"Skipped {mod_nt}.")
            else:
                print(f"Closed {mod_nt} without explicit action.")

            handled.add(mod_nt.upper())
            index += 1
            if index >= len(remaining):
                break
            app.load_mod_nt(remaining[index])

        try:
            root.destroy()
        except Exception:
            pass
