# This script is for developing and testing NA_pairwise_interactions.py


from NA_pairwise_interactions import *


from fr3d.localpath import outputText
from fr3d.localpath import outputNAPairwiseInteractions
from fr3d.localpath import outputNAPickleInteractions
from fr3d.localpath import contact_list_file
from fr3d.localpath import inputPath
from fr3d.localpath import outputHTML
from fr3d.data.base import EntitySelector


PDB_list = ['5AJ3']
PDB_list = ['6hiv']
PDB_list = ['3QRQ','5J7L']
PDB_list = ['4V9F','4YBB','4Y4O','6AZ3','4P95']
PDB_list = ['3BT7']
PDB_list = ['5I4A']
PDB_list = ['6A2H']
PDB_list = ['3JB9']
PDB_list = ['1OCT']
PDB_list = ['4v9fFH.pdb']
PDB_list = ['5KCR', '4WOI', '6C4I', '5JC9', '5L3P', '5KPW', '3J9Y', '3J9Z', '6BU8', '5WF0', '4V55', '4V54', '4V57', '4V56', '4V50', '4V53', '4V52', '4WF1', '5H5U', '4V5B', '5WFS', '5O2R', '5WFK', '5LZD', '5LZA', '6O9J', '6O9K', '6ORL', '6ORE', '3R8O', '3R8N', '4V85', '5MDV', '5MDW', '4V80', '4U27', '4U26', '4U25', '4U24', '4U20', '5KPS', '6GXM', '5KPX', '4U1U', '3JBU', '4V9P', '3JBV', '6Q9A', '6DNC', '4U1V', '6GXO', '5IQR', '5NWY', '4V9C', '6OSK', '4V9D', '4V9O', '5MGP', '6Q97', '3JCJ', '5J91', '3JCD', '3JCE', '6I7V', '6GXN', '4V64', '5J7L', '5AFI', '6BY1', '6ENU', '4V7V', '4V7U', '4V7T', '4V7S', '3JA1', '6ENF', '6OUO', '6ENJ', '5JU8', '5J8A', '6GWT', '4YBB', '5NP6', '5J88', '5U9G', '5U9F', '4V6D', '4V6E', '4V6C', '5JTE', '6OT3', '5J5B', '4WWW', '6OSQ', '5U4J', '5MDZ', '5U4I', '6NQB', '5UYQ', '5UYP', '5MDY', '5WDT', '6H4N', '5UYK', '4V89', '5UYM', '5UYL', '5UYN', '5WE6', '5WE4', '5KCS', '4V4Q', '4V4H', '5IT8']
PDB_list = ['4V51','4V9K']
PDB_list = ['6WJR']
PDB_list = ['6TPQ']
PDB_list = ['4KTG']
PDB_list = ['5KCR']
PDB_list = ['7ECF']  # DNA quadruplex
PDB_list = ['5J7L']
PDB_list = ['4V9F','5J7L','4ARC']
PDB_list = ['4ARC']
PDB_list = ['4ARC']
PDB_list = ['4V9F','6ZMI','7K00']
PDB_list = ['2N1Q']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.217/3.0A/csv']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.220/1.5A/csv']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.220/2.0A/csv']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.220/2.5A/csv']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.220/3.0A/csv']
PDB_list = ['4V9F']
PDB_list = ['203D']
PDB_list = ['7k00']
PDB_list = ['4V9F','6ZMI','7K00']
PDB_list = ['6CFJ']
PDB_list = ['7K00']
PDB_list = ['4TNA']
PDB_list = ['4V9F']
from DNA_2A_list import PDB_list   # define PDB_list as a list of DNA structures
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.224/2.5A/csv']
# list needed for a WebFR3D query
PDB_list = ['7MKY', '4JF2', '5VGW', '4ENC', '5XTM', '2R8S', '4PQV', '3RW6', '4BW0', '6CB3', '4K27', '5U3G', '7OF0', '4LVW', '5D5L', '2NZ4', '3NKB', '6TFG', '2Z75', '4YAZ', '5X2G', '4V9F', '7OZQ', '4Y4O', '4WFL', '1M5K', '7K16', '5FJC', '7O7Y', '6JQ5', '6S0Z', '3P22', '7OX9', '1Q96', '6KWQ', '3LQX', '6U8D', '6SVS', '3E5C', '7RQB', '2EZ6', '6DMC', '2V3C', '5M0I', '3MXH', '4YBB', '5B2P', '4P95', '7KKV', '3NPQ', '5DDP', '4NLF', '7P7Q', '6AZ3', '7D7W', '6S0X', '7RYG', '3AM1', '4PCJ', '5UZ6', '5B2T']
PDB_list += ['5AH5', '4ENC', '7EOG', '1QU2', '2ZUE', '2QUW', '2QUS', '5KPY', '7OF0', '7EQJ', '1U0B', '5AOX', '3FOZ', '2DRA', '4YCO', '7C79', '4V9F', '4Y4O', '4WFL', '3RG5', '5UD5', '7K16', '7O7Y', '6S0Z', '4JXZ', '4J50', '3B31', '3ADD', '7RQB', '3OVB', '6UGG', '4PRF', '4YBB', '3VJR', '1QTQ', '7K98', '4P95', '2GDI', '7DCO', '7P7Q', '6AZ3', '4YYE', '6S0X', '5HR7', '7RYG', '3AM1', '2OEU', '3D2V', '1J1U']
PDB_list = list(set(PDB_list))

PDB_list = ['3AM1','4J50']  # these have symmetry operators, but no annotations there
PDB_list = ['3AM1']  # these have symmetry operators, but no annotations there
PDB_list = ['4J50']  # these have symmetry operators, but no annotations there
PDB_list = ['4RKV']

base_seq_list = ['A','U','C','G']      # for RNA
base_seq_list = ['DA','DT','DC','DG']  # for DNA
base_seq_list = []                     # for all nucleic acids, modified or not

# tell which types of interactions to annotate
categories = {}
categories['sO'] = []        # annotate all sO interactions
categories['basepair'] = []
categories['stacking'] = []

OverwriteDataFiles = False   # to save time, if a data file exists, skip annotation
OverwriteDataFiles = True    # even if a data file already exists, annotate and overwrite

ShowStructureReadingErrors = True
ShowStructureReadingErrors = False

# this path should be specified in localpath.py
# intended for writing out a .pickle file to be used by the FR3D motif search tool
unit_data_path = "C:/Users/jimitch/Documents/GitHub/fr3d-python/data/units"

# annotate all nucleotides in all chains, even when a representative set is used
annotate_entire_PDB_files = True

timerData = myTimer("start")
lastwritetime = time()

allInteractionDictionary = defaultdict(list)

timerData = myTimer("Making PDB list",timerData)

PDB_IFE_Dict = map_PDB_list_to_PDB_IFE_dict(PDB_list)

#print("PDB_IFE_Dict is %s" % PDB_IFE_Dict)

counter = 0
count_pair = 0

# loop through 3D structures and annotate interactions
PDBs = PDB_IFE_Dict.keys()
#PDBs = PDBs[::-1]  # reverse the order of the list, for debugging

if len(PDBs) > 1 and not OverwriteDataFiles:
    print("Annotating interactions if no file is found in %s" % outputNAPairwiseInteractions)
else:
    print("Annotating interactions and saving in %s" % outputNAPairwiseInteractions)


for PDB in PDBs:

    counter += 1

    outputDataFileCSV =    outputNAPairwiseInteractions + PDB + ".csv"
    outputDataFilePickle = outputNAPickleInteractions + PDB + "_RNA_pairs_exp.pickle"

    if annotate_entire_PDB_files:

        pair_file = "%s_pairs_%s.pickle" % (PDB,fr3d_classification_version)
        pair_to_data_output_file = outputNAPairwiseInteractions + pair_file

        if not os.path.exists(pair_to_data_output_file) or len(PDBs) == 1 or OverwriteDataFiles:

            print("Reading file " + PDB + ", which is number "+str(counter)+" out of "+str(len(PDB_IFE_Dict)))
            timerData = myTimer("Reading CIF files",timerData)

            if ShowStructureReadingErrors:
                # do this to make sure to see any error messages
                structure = load_structure(os.path.join(inputPath,PDB+'.cif'))
            else:
                # do it this way to suppress error messages
                try:
                    structure = load_structure(os.path.join(inputPath,PDB+'.cif'))
                except:
                    print("  Could not load structure %s" % PDB)
                    print(inputPath)

                    continue

            # write out data file of nucleotide centers and rotations that can be used by FR3D for searches
            # need to be able to identify each chain that is available
            # write_unit_data_file(PDB,unit_data_path,structure)

            #interaction_to_list_of_tuples, pair_to_interaction, pair_to_data, timerData = annotate_nt_nt_in_structure(structure,timerData)
            # annotate interactions and return pair_to_data
            interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_in_structure(structure,categories,timerData,True)

            # turn this off during development and testing
            if True:
                print("  Annotated these interactions: %s" % interaction_to_list_of_tuples.keys())
                pickle.dump(interaction_to_list_of_tuples,open(outputDataFilePickle,"wb"),2)

            timerData = myTimer("Recording interactions",timerData)
            print('Writing data file %s' % pair_to_data_output_file)
            pickle.dump(pair_to_data,open(pair_to_data_output_file,"wb"),2)

            myTimer("summary",timerData)


    else:
        # only process individual IFEs
        # this has not been tested recently and may need to be modified
        # This would be most relevant for statistical tallies, finding exemplars, etc.
        # But that could be done by just downloading the annotations, not creating them anew
        print("Reading file " + PDB + ", which is number "+str(counter)+" out of "+str(len(PDB_IFE_Dict)))
        timerData = myTimer("Reading CIF files",timerData)

        if ShowStructureReadingErrors:
            # do this to make sure to see any error messages
            structure = load_structure(os.path.join(inputPath,PDB+'.cif'))
        else:
            # do it this way to suppress error messages
            try:
                structure = load_structure(os.path.join(inputPath,PDB+'.cif'))
            except:
                print("Could not load structure %s" % PDB)
                continue

        # extract nucleotides to analyze
        IFE = PDB_IFE_Dict[PDB]          #
        if len(IFE) == 0:                # use the whole PDB file
            if base_seq_list:
                bases = structure.residues(sequence = base_seq_list)  # load just the types of bases in base_seq_list
            else:
                bases = structure.residues(type = ["RNA linking","DNA linking"])  # load all RNA/DNA nucleotides
        else:                            # use specific chains only
            chain_ids = []
            print("  Keeping only bases in chains %s" % IFE)
            chains = IFE.split("+")
            for chain in chains[1:]:            #skip element zero, leading +
                fields = chain.split("|")
                chain_ids.append(fields[2])
            if base_seq_list:
                bases = structure.residues(chain = chain_ids, sequence = base_seq_list)  # load just the types of bases in base_seq_list
            else:
                bases = structure.residues(chain = chain_ids)  # load all bases

        # ??? record which RNA/DNA chains are actually present
        # count nucleotides
        numBases = 0
        for base in bases:
            numBases += 1
        print("  Found " + str(numBases) + " bases in " + PDB)

        # build cubes to be able to find potential pairs quickly
        timerData = myTimer("Building cubes",timerData)
        print("  Building nucleotide cubes in " + PDB)
        baseCubeList, baseCubeNeighbors = make_nt_cubes(bases, nt_nt_screen_distance, nt_reference_point)

        # annotate nt-nt interactions
        timerData = myTimer("Annotating interactions",timerData)

        # annotate interactions and return pair_to_data
        interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_in_structure(structure,categories,timerData,True)

        # used to return Python_pairs, pair_to_data, timerData

        timerData = myTimer("Recording interactions",timerData)

        # write out pairs in the format that WebFR3D reads
        # accumulate list of interacting units by base, interaction type, and edges
        for nt1, nt2, interaction, edge, standard_aa, param in list_nt_nt:
            base = base_residue.unit_id()
            # skip symmetry operated instances; generally these are just duplicates anyway
            if not "||||" in str(base):
                aa = aa_residue.unit_id()
                base_component = str(base).split("|")
                aa_component = str(aa).split("|")
                key = base_component[3]+"_"+aa_component[3]+"_"+interaction+"_"+edge
                count_pair += 1
                allInteractionDictionary[key].append((base,aa,interaction,edge,standard_aa,param))  # store tuples

        myTimer("summary",timerData)

# when appropriate, write out HTML files
"""
if len(PDB_IFE_Dict) > 100:
    print("Writing " + outputDataFile)
    timerData = myTimer("Writing HTML files",timerData)
    pickle.dump((allInteractionDictionary,allAATwoBaseDictionary,PDB_list),open(outputDataFile,"wb"))
    writeInteractionsHTML(allInteractionDictionary,outputHTML,version)
"""

myTimer("summary",timerData)

