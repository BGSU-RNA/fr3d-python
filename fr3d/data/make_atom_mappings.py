
import os
import pdbx
import sys
from sys import argv

# import the version of urlretrieve appropriate to the Python version
if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
else:
    from urllib.request import urlretrieve as urlretrieve


def map_atom_to_neighbors(chem_comp_bond):
    """
    Create a dictionary whose keys are atoms and whose values are a set of neighboring atoms
    """
    mod_atom_to_atom = {}
    for line in chem_comp_bond:
        if not line[1] in mod_atom_to_atom:    # not already a key in the dictionary
            mod_atom_to_atom[line[1]] = set([])  # start with an empty set
        if not line[2] in mod_atom_to_atom:    # not already a key in the dictionary
            mod_atom_to_atom[line[2]] = set([])  # start with an empty set
        mod_atom_to_atom[line[1]].add(line[2])
        mod_atom_to_atom[line[2]].add(line[1])
    return mod_atom_to_atom


def reverse_map(A_to_B):
    """
    Reverse an atom mapping
    """

    B_to_A = {}
    for A,B in A_to_B.items():
        B_to_A[B] = A

    return B_to_A


def process_one_modified_nt(mod_nt):
    """
    Input mod_nt is a string telling the name of a modified nucleotide like 'PSU' or '2MG'
    The function reads the corresponding .cif file, determines the parent nucleotide,
    reads the parent .cif file, and finds as many corresponding atoms as possible.

    It also reads a file of manual mappings.
    each manual mapping is a 4-tuple like ("A","N9","ZAD","N9")
    if a modified nucleotide is not mappable, the entry will be like ("","","YRR","")
    to simply indicate the mapping of modified nucleotide to standard, use ("A","","ZAD","") with no atoms
    """

    modified_to_mappings = {}
    with open("atom_mappings_manual.txt","r") as f:
        lines = f.readlines()
    for line in lines:
        # make it easy to comment out lines with #
        if line.startswith("#"):
            continue
        fields = line.replace("\n","").split("\t")

        if len(fields) >= 4:
            # OK to have more than 4 fields, but only the first 4 are used
            parent,par_atom,modified,mod_atom = fields[0:4]
            if not modified in modified_to_mappings:
                modified_to_mappings[modified] = []
            modified_to_mappings[modified].append((parent,par_atom,modified,mod_atom))

    print("Found %d modified nucleotides in manual mapping file" % len(modified_to_mappings.keys()))

    # map the given modified nucleotide

    # try to determine the best matching standard nucleotide
    standard_nt = None
    # keep track of parent atom to modified atom mapping
    standard_to_mod = {}

    if mod_nt in modified_to_mappings:
        standard_nt = modified_to_mappings[mod_nt][0][0]
        if len(standard_nt) == 0:
            # manually annotated that mod_nt does not map to a nucleotide
            return None

    # Open the CIF file for the modified nucleotide
    if mod_nt in ["PRN"]:
        # Windows restriction on using PRN as a filename
        cif = open(os.path.join("cif","data_" + mod_nt + ".cif"))
    else:
        cif = open(os.path.join("cif",mod_nt + ".cif"))

    # A list to be propagated with data blocks
    data = []

    # Create a PdbxReader object
    pRd = pdbx.reader.PdbxReader(cif)

    # Read the CIF file, propagating the data list
    pRd.read(data)

    # Close the CIF file, as it is no longer needed
    cif.close()

    # Retrieve the first data block
    data = data[0]

    # Retrieve the struct_conn category table, which delineates connections

    # We don't actually need the atom list here
    # chem_comp_atom = data.get_object("chem_comp_atom")
    # chem_comp_atom.print_it()

    # get a list of atoms in the modified nucleotide
    chem_comp_bond = data.get_object("chem_comp_bond")
    atoms = set([])
    for line in chem_comp_bond:
        atoms.add(line[1])         # assumes that atoms will also be in these columns
        atoms.add(line[2])

    chem_comp = data.get_object("chem_comp")
    # print("chem_comp print_it")
    # chem_comp.print_it()
    # print(help(chem_comp))
    # print("item_name_list")
    # print(chem_comp.item_name_list)

    # use the .cif file chem_comp block to determine the parent, if possible
    NA_type = chem_comp.get_value('type')
    cif_standard_nt = None
    cif_standard_nt = chem_comp.get_value('mon_nstd_parent_comp_id')

    print("NA_type from cif file:           %s"     % NA_type)
    print("standard_nt from .cif file:      %s" % (cif_standard_nt))
    print("standard_nt from manual mapping: %s" % (standard_nt))

    apparently_wrong_cif_comp_id = ["T0N","T0Q","2EG","6MA"]

    if standard_nt and cif_standard_nt and not standard_nt == cif_standard_nt:
        if not "," in cif_standard_nt and not mod_nt in apparently_wrong_cif_comp_id:
            if cif_standard_nt == "DU" and standard_nt == "U":
                print("CIF says DU")
            else:
                print(crashnow)
    elif cif_standard_nt and not standard_nt:
        standard_nt = cif_standard_nt

    # standard nts that we are not going to model
    if standard_nt == "DI":
        standard_nt = "DG"

    if standard_nt == "DU" and "RNA" in NA_type:
        standard_nt = "U"
    elif standard_nt == "DU" and "DNA" in NA_type:
        # keep it in the DNA family, even though deoxy U is the parent
        standard_nt = "DT"


    output = ""    # start accumulating a string of output

    # if there is no parent in the .cif file, use atoms in the .cif file to guess the parent
    if not standard_nt:
        if 'N9' in atoms and 'O6' in atoms:
            standard_nt = 'G'
        elif 'N9' in atoms and 'N6' in atoms:
            standard_nt = 'A'
        elif 'N1' in atoms and 'N4' in atoms:
            standard_nt = 'C'
        elif 'N1' in atoms and 'O4' in atoms:
            standard_nt = 'U'

        old_standard_nt = None
        if standard_nt in ["A","C","G"] and not "O2'" in atoms:
            old_standard_nt = "D" + standard_nt
        elif standard_nt == "U" and not "O2'" in atoms:
            old_standard_nt = "DT"

        if standard_nt in ["A","C","G"] and "DNA" in NA_type:
            standard_nt = "D" + standard_nt
        elif standard_nt in ["A","C","G"] and "RNA" in NA_type:
            pass
        elif standard_nt in ["A","C","G"] and not "O2'" in atoms:
            standard_nt = "D" + standard_nt
        elif standard_nt == "U" and not "O2'" in atoms:
            standard_nt = "DT"

        if not standard_nt == old_standard_nt:
            print("%s changed from %s to %s" % (mod_nt,old_standard_nt,standard_nt))

        print('Making a guess that the parent of %s is %s' % (mod_nt,standard_nt))
        output += 'Making a guess that the parent of %s is %s\n' % (mod_nt,standard_nt)

    print("Standard nucleotide is %s" % standard_nt)

    if "," in standard_nt:
        output = "%s is two parents combined (%s) and we have no plan for which one to map\n" % (mod_nt, standard_nt)
        return output

    if standard_nt == "None" or standard_nt == None:
        output = "No standard nucleotide identified for %s\n" % (mod_nt)
        return output

    cif_parent = open(os.path.join("cif",standard_nt + ".cif"))

    # A list to be propagated with data blocks
    data_parent = []

    # Create a PdbxReader object
    pRd = pdbx.reader.PdbxReader(cif_parent)

    # Read the CIF file, propagating the data list
    pRd.read(data_parent)

    # Close the CIF file, as it is no longer needed
    cif_parent.close()

    # Retrieve the first data block
    data_parent = data_parent[0]

    #parent_chem_comp_atom = data_parent.get_object("chem_comp_atom")
    #parent_chem_comp_atom.print_it()

    parent_chem_comp_bond = data_parent.get_object("chem_comp_bond")
    """
    for line in parent_chem_comp_bond:
        print(line[1], line[2])
    """
    parent_atoms = set([])
    for line in parent_chem_comp_bond:
        atom = line[1]
        parent_atoms.add(atom)
        atom = line[2]
        parent_atoms.add(atom)

    #print(sorted(atoms))
    #print(sorted(parent_atoms))
    #print("difference")
    #print("Atoms in mod, not in par: ")
    #Modified has / parent doesn't
    #print(atoms.difference(parent_atoms))
    #print("Atoms in par, not in mod: ")
    #Modified doesn't have / parent does
    #print(parent_atoms.difference(atoms))

    """
    # difference example
    seta = set([1, 2, 4, 5])
    setb = set([1, 2 ,3 ,4])
    print(seta.difference(setb))
    print(seta-setb)
    print(setb.difference(seta))
    print(setb-seta)
    print(seta & setb) #set intersection
    print(seta | setb) #set union
    """

    #print('1a')

    mod_atom_to_neighbors = map_atom_to_neighbors(chem_comp_bond) # Gives a set of atoms and its corresponding neighbors for the mod
    #print(mod_atom_to_neighbors)
    par_atom_to_neighbors = map_atom_to_neighbors(parent_chem_comp_bond)
    #print(par_atom_to_neighbors)

    #print('2a')

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
    if mod_nt in modified_to_mappings:
        for mapping in modified_to_mappings[mod_nt]:
            if mapping[1] and mapping[3]:
                standard_to_mod[mapping[1]] = mapping[3]

    if not "C1'" in standard_to_mod and "C1'" in mod_atom_to_neighbors and "C1'" in par_atom_to_neighbors:
        # The usual starting point is the C1' atom
        standard_to_mod["C1'"] = "C1'"   # map C1' to C1'

    # reverse the mappings listed so far; then they don't need to be typed twice
    mod_to_par = reverse_map(standard_to_mod)

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
                    if not new_atom in ['H21','H22','H41','H42','H61','H62','H71','H72','H73']:
                        # don't map these, since they are often switched and can be inferred later
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
            output += "%s\t%s\t%s\t%s\n" % (standard_nt,par_atom,mod_nt,mod_atom)

    else:
        print("Not sure what atom to start at")
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

    with open('modified_nt_list.csv','r') as f:
        lines = f.read()                # read the entire file
        lines = lines.replace("\r","")  # remove \r return characters
        lines = lines.split("\n")       # split on newline character, return a list

    # get information about each modified nucleotide
    output = ""
    for line in lines:
        fields = line.split(",")  # split line of data into list, usually length 3
        if len(fields) == 3:
            mod_nt = fields[1]
            mod_nt_url = "https://www.rcsb.org/ligand/" + mod_nt  # for viewing
            mod_nt_count = int(fields[2])
            mod_nt_number = int(fields[0])
            print("")
            print("Processing number %3d %4s which has count %4d and url %s" % (mod_nt_number,mod_nt,mod_nt_count,mod_nt_url))

            mod_filename = mod_nt + ".cif"

            mod_path_and_filename = os.path.join('cif',mod_filename)  # cif folder, PSU.cif filename

            if mod_nt == 'PRN':
                # Windows restriction on files named PRN
                mod_path_and_filename = mod_path_and_filename.replace("PRN","data_PRN")

            if not os.path.exists(mod_path_and_filename):  # if file does not already exist
                # download from RCSB website and store in the cif folder
                urlretrieve("https://files.rcsb.org/ligands/download/" + mod_filename, mod_path_and_filename)

            if not mod_nt in ['A','C','G','U','DA','DC','DG','DT']:
                # run outside of try/except to catch errors
                new_output = process_one_modified_nt(mod_nt)

                try:
                    #new_output = process_one_modified_nt(mod_nt)
                    if new_output:
                        output += new_output
                    else:
                        output += "\t\t%s\t\n" % mod_nt
                except:
                    output += "\t\t%s\t\n" % mod_nt
                    # output += "Not able to process %s\n" % mod_nt

    with open("atom_mappings_provisional.txt","w") as f:
        f.write(output)


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
