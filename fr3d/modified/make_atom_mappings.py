"""
Read a set of manual mappings from a file atom_mappings_manual.txt
Loop over a list of modified nucleotides from CSV file modified_nt_list.csv
For each modified nucleotide, download the corresponding .cif file
Determine the parent nucleotide
Read the parent nucleotide .cif file
Map the atoms of the modified nucleotide to the parent nucleotide
"""

import os
import pdbx
import sys
from sys import argv

# import the version of urlretrieve appropriate to the Python version
if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
else:
    from urllib.request import urlretrieve as urlretrieve


def read_monomer_cif(mod_nt):
    """
    One function to read all necessary data from the .cif file for this project
    """

    # set the filename for the .cif file
    if mod_nt in ["PRN"]:
        # Windows restriction on using PRN as a filename
        filename = os.path.join("cif","data_" + mod_nt + ".cif")
    else:
        filename = os.path.join("cif",mod_nt + ".cif")

    if not os.path.exists(filename):
        # download from RCSB website and store in the cif folder
        url = "https://files.rcsb.org/ligands/download/%s.cif" % mod_nt
        urlretrieve(url,filename)

    if not os.path.exists(filename):
        print("Could not find %s" % filename)
        return None

    # Open the CIF file for the modified nucleotide
    cif = open(filename)
    data = []
    pRd = pdbx.reader.PdbxReader(cif)
    pRd.read(data)
    cif.close()
    data = data[0]

    category_names = data.get_object_name_list()
    cif_data = {}

    for category_name in category_names:
        # print('category_name',category_name)

        cif_data[category_name] = []
        category = data.get_object(category_name)

        row_count = category.row_count
        # print('  row_count',row_count)

        for i in range(row_count):
            row = {}
            for item_name in category.item_name_list:
                # print('    ',item_name)
                name = item_name.split('.')[1]
                # print('    ',name)
                # print('    ',mod_nt,category_name,i,item_name,category.get_value(name,i))
                row[name] = category.get_value(name,i)
            cif_data[category_name].append(row)
            # if category_name == "chem_comp":
            #     print(row)

    return cif_data


def map_atom_to_neighbors(cif_data):
    """
    Create a dictionary whose keys are atoms and whose values are a set of neighboring atoms
    """

    atom_to_neighbor_atoms = {}

    for row in cif_data['chem_comp_bond']:
        # get the two interacting atoms
        atom1 = row['atom_id_1']
        atom2 = row['atom_id_2']

        # add entries to the dictionary
        if not atom1 in atom_to_neighbor_atoms:    # not already a key in the dictionary
            atom_to_neighbor_atoms[atom1] = set([])  # start with an empty set
        if not atom2 in atom_to_neighbor_atoms:    # not already a key in the dictionary
            atom_to_neighbor_atoms[atom2] = set([])  # start with an empty set

        # store the connection in both directions
        atom_to_neighbor_atoms[atom1].add(atom2)
        atom_to_neighbor_atoms[atom2].add(atom1)

    return atom_to_neighbor_atoms


def reverse_map(A_to_B):
    """
    Reverse an atom mapping
    """

    B_to_A = {}
    for A,B in A_to_B.items():
        B_to_A[B] = A

    return B_to_A


def process_one_modified_nt(mod_nt,mappings=[]):
    """
    Input mod_nt is a string telling the name of a modified nucleotide like 'PSU' or '2MG'
    The function reads the corresponding .cif file, determines the parent nucleotide,
    reads the parent .cif file, and finds as many corresponding atoms as possible.

    It also reads a file of manual mappings.
    each manual mapping is a 4-tuple like ("A","N9","ZAD","N9")
    if a modified nucleotide is not mappable, the entry will be like ("","","YRR","")
    to simply indicate the mapping of modified nucleotide to standard, use ("A","","ZAD","") with no atoms
    """

    standard_nt = None

    # check if the manual mappings say that no mapping is possible
    if len(mappings) > 0:
        standard_nt = mappings[0][0]
        if len(standard_nt) == 0:
            # manually annotated that mod_nt does not map to a nucleotide
            return None

    # read the .cif file for the modified nucleotide
    mod_data = read_monomer_cif(mod_nt)

    # print(mod_data['chem_comp'][0]['mon_nstd_parent_comp_id'])
    # print(mod_data['chem_comp'][0]['one_letter_code'])

    NA_type = mod_data['chem_comp'][0].get('type',None)
    cif_standard_nt = mod_data['chem_comp'][0].get('mon_nstd_parent_comp_id',None)
    one_letter_code = mod_data['chem_comp'][0].get('one_letter_code',None)

    print("NA_type from cif file:           %s" % (NA_type))
    print("standard_nt from .cif file:      %s" % (cif_standard_nt))
    print("one_letter_code from .cif file:  %s" % (one_letter_code))
    print("standard_nt from manual mapping: %s" % (standard_nt))

    # get a list of atoms in the modified nucleotide
    mod_atoms = set([])
    for row in mod_data['chem_comp_atom']:
        mod_atoms.add(row['atom_id'])

    # use standard_nt from the .cif file
    if cif_standard_nt and not standard_nt:
        standard_nt = cif_standard_nt

    # standard nts that we are not going to model
    if standard_nt == "DI":
        standard_nt = "DG"

    if standard_nt == "DU" and "RNA" in NA_type:
        standard_nt = "U"
    elif standard_nt == "DU" and "DNA" in NA_type:
        # keep it in the DNA family, even though deoxy U is the standard nt
        standard_nt = "DT"

    # map these like the others, so their plot gets made
    if mod_nt in ["A","C","G","U","DA","DC","DG","DT"]:
        standard_nt = mod_nt

    output = ""    # start accumulating a string of output

    # if there is no parent in the .cif file, use atoms in the .cif file to guess the parent
    if not standard_nt:
        if 'N9' in mod_atoms and 'O6' in mod_atoms:
            standard_nt = 'G'
        elif 'N9' in mod_atoms and 'N6' in mod_atoms:
            standard_nt = 'A'
        elif 'N1' in mod_atoms and 'N4' in mod_atoms:
            standard_nt = 'C'
        elif 'N1' in mod_atoms and 'O4' in mod_atoms:
            standard_nt = 'U'

        old_standard_nt = None
        if standard_nt in ["A","C","G"] and not "O2'" in mod_atoms:
            old_standard_nt = "D" + standard_nt
        elif standard_nt == "U" and not "O2'" in mod_atoms:
            old_standard_nt = "DT"

        if standard_nt in ["A","C","G"] and "DNA" in NA_type:
            standard_nt = "D" + standard_nt
        elif standard_nt in ["A","C","G"] and "RNA" in NA_type:
            pass
        elif standard_nt in ["A","C","G"] and not "O2'" in mod_atoms:
            standard_nt = "D" + standard_nt
        elif standard_nt == "U" and not "O2'" in mod_atoms:
            standard_nt = "DT"

        if not standard_nt == old_standard_nt:
            print("%s changed from %s to %s" % (mod_nt,old_standard_nt,standard_nt))

        if not standard_nt:
            return None

        print('Making a guess that the parent of %s is %s' % (mod_nt,standard_nt))
        output += 'Making a guess that the parent of %s is %s\n' % (mod_nt,standard_nt)

    print("Standard nucleotide we use is:     %s" % standard_nt)

    # error messages
    # resolve them by adding the appropriate entry to atom_mappings_manual.txt
    if "," in standard_nt:
        output = "%s is two parents combined (%s) and we have no plan for which one to map\n" % (mod_nt, standard_nt)
        return output

    if not standard_nt in ['A','C','G','U','DA','DC','DG','DT']:
        output = "mon_nstd_parent_comp_id %s is not a nucleotide %s\n" % (standard_nt,mod_nt)
        return output

    if standard_nt == "None" or standard_nt == None:
        output = "No standard nucleotide identified for %s\n" % (mod_nt)
        return output

    parent_data = read_monomer_cif(standard_nt)

    parent_atoms = set([])
    for row in parent_data['chem_comp_atom']:
        parent_atoms.add(row['atom_id'])

    mod_atom_to_neighbors = map_atom_to_neighbors(mod_data) # Gives a set of atoms and its corresponding neighbors for the mod
    # print(mod_atom_to_neighbors)
    par_atom_to_neighbors = map_atom_to_neighbors(parent_data)
    # print(par_atom_to_neighbors)


    """
    for key in sorted(mod_atom_to_neighbors.keys()):
        print(key, mod_atom_to_neighbors[key])

        for par_key in par_atom_to_neighbors.keys():
            if key == par_key:
                if mod_atom_to_neighbors[key] == par_atom_to_neighbors[par_key]:
                    print("Neighbors are the same!")
                if mod_atom_to_neighbors[key] != par_atom_to_neighbors[par_key]:
                    print("Neighbors are not the same!")
    """

    # use manual mappings to get started
    # In some nucleotides, C1' attaches to an atom different from N1/N9
    # In those cases, add the next match, then the algorithm may handle the rest.
    # In other cases, C1' has a different name like C1G or C1A
    # In other cases, the atom names are all different, so the algorithm needs more help

    # keep track of parent atom to modified atom mapping
    standard_to_mod = {}
    mod_to_par = {}

    for mapping in mappings:
        if mapping[1] and mapping[3]:
            standard_to_mod[mapping[1]] = mapping[3]
            mod_to_par[mapping[3]] = mapping[1]

    if not "C1'" in standard_to_mod and "C1'" in mod_atom_to_neighbors and "C1'" in par_atom_to_neighbors:
        # The usual starting point is the C1' atom
        standard_to_mod["C1'"] = "C1'"   # map C1' to C1'
        mod_to_par["C1'"] = "C1'"        # map C1' to C1'

    # # reverse the mappings listed so far; then they don't need to be typed twice
    # mod_to_par = reverse_map(standard_to_mod)

    # new atom(s) that the algorithm should spread from
    new_par_atoms = list(standard_to_mod.keys())

    # print('new_par_atoms')
    # print(new_par_atoms)

    if len(new_par_atoms) > 0:
        while len(new_par_atoms) > 0:
            par_atoms = new_par_atoms
            new_par_atoms = []
            for par_atom in par_atoms:
                mod_atom = standard_to_mod[par_atom] # Gives back the modified atom related to the parent atom
                if mod_atom in mod_atom_to_neighbors:
                    mod_neighbors = (mod_atom_to_neighbors[mod_atom]) - set(mod_to_par.keys()) # subtracts out previously mapped atoms
                else:
                    for ma in sorted(mod_atom_to_neighbors.keys()):
                        print(ma, mod_atom_to_neighbors[ma])
                    print('Warning: Atom %s has no neighbors' % (mod_atom))
                par_neighbors = (par_atom_to_neighbors[par_atom]) - set(standard_to_mod.keys()) # subtracts out previously mapped atoms
                """
                if mod_neighbors == par_neighbors:
                    new_atoms = mod_neighbors
                    for new_atom in new_atoms: # Adds the neighbors of "atom" to standard_to_mod and mod_to_par
                        standard_to_mod[new_atom] = new_atom
                        mod_to_par[new_atom] = new_atom
                        new_par_atoms.append(new_atom)
                   # print(standard_to_mod)
                else:
                """
                # find neighboring atoms that have the same names and map them
                new_atoms = mod_neighbors & par_neighbors
                for new_atom in new_atoms:
                    if new_atom in ['H21','H22','H41','H42','H61','H62','H71','H72','H73',"H5'","H5''","H2''"]:
                        # don't map these, since they are often switched and can be inferred later
                        pass
                    elif new_atom == "H2'" and standard_nt in ["DA","DC","DG","DT"]:
                        # DNA has H2' and H2''
                        pass
                    else:
                        standard_to_mod[new_atom] = new_atom
                        mod_to_par[new_atom] = new_atom
                        new_par_atoms.append(new_atom)

                # special cases of atom naming conventions
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

                # take these out of consideration, so we can focus on unmapped atoms
                mod_diff_atoms = mod_neighbors - new_atoms
                par_diff_atoms = par_neighbors - new_atoms

                # if there is exactly one neighboring atom that does not have the same name
                if len(mod_diff_atoms) == 1 and len(par_diff_atoms) == 1:
                    par_diff_atom = list(par_diff_atoms)[0]

                    if par_diff_atom[0] == 'H' and 'H' in [x[0] for x in new_atoms]:
                        # don't map a hydrogen atom when there are other hydrogens to map here
                        pass
                    else:
                        mod_diff_atom = list(mod_diff_atoms)[0]
                        print("%-5s maps to %-5s" % (par_diff_atom, mod_diff_atom))
                        standard_to_mod[par_diff_atom] = mod_diff_atom
                        mod_to_par[mod_diff_atom] = par_diff_atom
                        new_par_atoms.append(par_diff_atom)
                elif len(mod_diff_atoms) == 2 and len(par_diff_atoms) == 2:
                    # print(par_diff_atoms)
                    # print(mod_diff_atoms)
                    p1, p2 = sorted(par_diff_atoms)
                    m1, m2 = sorted(mod_diff_atoms)
                    print('Wondering how to map %s and %s to %s and %s' % (p1,p2,m1,m2))

                    if (p1[0] == m1[0] and p1[0] != m2[0]) or (p2[0] == m2[0] and p2[0] != m1[0]):
                        print('Guessing a map of %s to %s and %s to %s but not using it' % (p1,m1,p2,m2))
                        # standard_to_mod[p1] = m1
                        # mod_to_par[m1] = p1
                        # new_par_atoms.append(p1)

                        # standard_to_mod[p2] = m2
                        # mod_to_par[m2] = p2
                        # new_par_atoms.append(p2)
                    else:
                        pass
                        #print("Neighbors %s and %s, still not sure what to do" % (par_diff_atoms,mod_diff_atoms))

                else:
                    if len(mod_diff_atoms) > 0 or len(par_diff_atoms) > 0:
                        print(str(par_diff_atoms) + " and " + str(mod_diff_atoms) + " has no atom names that match")

        #print(standard_to_mod)

        """
        # loop over key,value pairs in the dictionary standard_to_mod
        for par_atom,mod_atom in standard_to_mod.items():
            # output four columns for a spreadsheet
            output += "%s\t%s\t%s\t%s\n" % (chem_part,par_atom,mod_nt,mod_atom)
        """

        # loop over keys, get the value, make a string for m
        for par_atom in sorted(standard_to_mod.keys()):
            mod_atom = standard_to_mod[par_atom]
            # output four columns for a spreadsheet

            if mod_atom in mod_atoms:
                output += "%s\t%s\t%s\t%s\n" % (standard_nt,par_atom,mod_nt,mod_atom)
            else:
                # Over time, the atom list in the .cif file can change, making manual mappings obsolete
                output += "Warning: %s not in %s\n" % (mod_atom,mod_nt)

    else:
        print("Not sure what atom to start at with %s" % mod_nt)
        output = "Not sure what atom to start at with %s \n" % mod_nt

    return output


def map_all_modified_nucleotides():
    """
    Procedure for mapping atoms of non-standard nucleotides
    First, get a list of non-standard nucleotides
    Go to http://www.nakb.org/searchterms.html
    Check Nonstandard NA residues, copy and paste into a spreadsheet, save as CSV format
    Copy the list into atom_mappings.txt, preserving the rows for A, C, G, U, DA, DC, DG, DT
    Second, run make_atom_mappings.py which will download necessary .cif files and make a first mapping
    The output file is called atom_mappings.txt, which contains mappings and messages
    Third, run refine_atom_mappings.py to fill in missing atoms

    """

    # read manual mappings
    modified_to_mappings = {}
    with open("atom_mappings_manual.txt","r") as f:
        lines = f.readlines()
    for line in lines:
        # make it easy to comment out lines with #
        if line.startswith("#"):
            continue
        fields = line.replace("\n","").split("\t")

        if len(fields) >= 3:
            # OK to have more than 4 fields, but only the first 4 are used
            if len(fields) == 4:
                parent,par_atom,modified,mod_atom = fields[0:4]
            else:
                parent,par_atom,modified = fields[0:3]
                mod_atom = ""
            if not modified in modified_to_mappings:
                modified_to_mappings[modified] = []
            modified_to_mappings[modified].append((parent,par_atom,modified,mod_atom))

    print("Found %d modified nucleotides in manual mapping file" % len(modified_to_mappings.keys()))


    with open('modified_nt_list.csv','r') as f:
        lines = f.read()                # read the entire file
        lines = lines.replace("\r","")  # remove \r return characters
        lines = lines.replace('"','')   # remove double quotes
        lines = lines.split("\n")       # split on newline character, return a list

    # get information about each modified nucleotide
    output = ""
    skipped_mod_nt = []
    for line in lines:
        fields = line.split(",")  # split line of data into list, usually length 3
        if len(fields) == 3:
            mod_nt = fields[1]
            mod_nt_url = "https://www.rcsb.org/ligand/" + mod_nt  # for viewing
            mod_nt_count = int(fields[2])
            mod_nt_number = int(fields[0])
            print("")
            print("Processing number %3d %4s which has count %4d and url %s" % (mod_nt_number,mod_nt,mod_nt_count,mod_nt_url))

            # map the standard nucleotides so we produce images colored just like the modified ones
            # if mod_nt in ['A','C','G','U','DA','DC','DG','DT']:
            #     # no need to map these
            #     continue

            new_output = ""
            if mod_nt in modified_to_mappings:
                manual_mappings = modified_to_mappings[mod_nt]
            else:
                manual_mappings = []
            new_output = process_one_modified_nt(mod_nt,manual_mappings)
            if new_output:
                output += new_output

            if not new_output:
                skipped_mod_nt.append(mod_nt)
                svg_filename = os.path.join('skipped',mod_nt + ".svg")
                if not os.path.exists(svg_filename):
                    svg_url = "https://cdn.rcsb.org/images/ccd/unlabeled/%s/%s.svg" % (mod_nt[0],mod_nt)
                    urlretrieve(svg_url, svg_filename)


    with open("atom_mappings_provisional.txt","w") as f:
        f.write(output)

    print("Skipped the following %d modified nucleotides:" % len(skipped_mod_nt))
    print(sorted(skipped_mod_nt))

    with open("skipped/skipped.html","w") as f:
        c = 0
        for mod_nt in sorted(skipped_mod_nt):
            c += 1
            f.write("<h2>%s number %d of %d</h2>\n" % (mod_nt,c,len(skipped_mod_nt)))
            f.write('<a href="https://www.rcsb.org/ligand/%s" target = "_blank">%s in ligand explorer</a><br>\n' % (mod_nt,mod_nt))
            f.write('<a href="https://www.rcsb.org/ligand/%s" target = "_blank"><img src="%s.svg" height="300"></a>\n' % (mod_nt,mod_nt))


if __name__=="__main__":

    # Process command line arguments, if any

    print("arguments")
    print(argv)

    if len(argv) == 1:
        map_all_modified_nucleotides()
    elif len(argv) == 2:
        mod_nt = argv[1]
        print('Got modified nucleotide %s' % mod_nt)
        print(process_one_modified_nt(mod_nt))
