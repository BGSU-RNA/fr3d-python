Methodology for adding modified nucleotides.

Optional step:
    Go to https://www.nakb.org/searchterms.html
    Click Nonstandard NA Residues, to get all modified residues that are in polymers
    Click CSV
    Save the CSV in this folder as modified_nt_list.csv

If you have never run this code, run these first to establish a baseline:
    python311 make_atom_mappings.py
    python311 refine_atom_mappings.py

The program make_atom_mappings.py will also try to download the list of modified nucleotides from NAKB.

Make provisional mappings of all modified nucleotides to standard:
    python311 make_atom_mappings.py
The program make_atom_mappings.py writes the file atom_mappings_provisional.txt.
Newly appeared modified residues might be processed correctly by this program; we'll check that later.

Have a look at skipped/skipped.html to see non-standard residues that are not being mapped even provisionally.
If you find one you can map, then add it to atom_mappings_manual.txt, following the format there.
Then run make_atom_mappings.py again.

Check the file refine_atom_mappings.py to set:
    color_scheme = 'diagnostic'  # use many colors, to show the atom mappings
    draw_figures = True      # draw new figures if they don't already exist
    overwrite_figures = False # makes it easier to identify what is new

Now process the provisional mappings into the actual mappings by running:
    python311 refine_atom_mappings.py
The program refine_atom_mappings.py writes the files atom_mappings.txt, modified_to_changes.json, and image files.
In the folder diagnostic/base_plots, sort the images by modified date to see the most recent ones.
Each image should show the standard base on the left and the modified base on the right.
Atoms that are mapped to each other are colored in the same color, even if the element has changed.
Make sure that all atoms that should be mapped are mapped, and mapped correctly.
Check the text.  If it says something like "H61 removed, HN61 added" maybe you should map H61 to HN61.
If not, type correct mappings into atom_mappings_manual.txt following the format and being sure
to separate columns with tab characters (not spaces).
You should not need to type in all atom to atom mappings, just enough to guide the program.
If the C1' atom has a new name, make sure to map the C1' atom to the corresponding atom.
Don't be fooled by a phosphate group attached to C3' or O3'; that does not correspond to the
phosphate group on a standard nucleotide.

After you add or fix mappings, delete the newest image file(s) in base_plots so they
get made again, then run:
    python311 make_atom_mappings.py
    python311 refine_atom_mappings.py
and repeat the examination until you are happy with the base mappings.
Next, check the backbone mappings by going to diagnostic/backbone_plots, then delete those plots,
and run the programs again.
It can be hard to see all the atoms in the backbone mappings; you might need to use the ligand explorer and
inspect the file atom_mappings.txt to see what mappings are actually made.

To map modified RNA nucleotides to Modomics, run map_pdb_to_modomics.py, then open pdb_to_modomics.html
to view and open pdb_to_modomics.txt to edit, perhaps with Excel since it needs to be tab delimited.
Look for rows of the table with TODO in them.
Carefully compare the 3D coordinates from PDB and from Modomics to see if there is an exact match.
Pay close attention to chirality.
Mark exact matches in pdb_to_modomics.txt as "confirmed" and incorrect matches as "no".  Save.
Run map_pdb_to_modomics.py again to update the .txt and .html files.
