# This script is for developing and testing NA_pairwise_interactions.py

"""
When changes are made to other code in fr3d-python, as administrator:
cd c:/Users/zirbel/Documents/GitHub/fr3d-python
python27 -m pip install .
python38 -m pip install .
python311 -m pip install .

To run the code:
cd c:/Users/zirbel/Documents/GitHub/fr3d-python/fr3d/classifiers
python27 NA_pairwise_interactions.py -c basepair,stacking,sugar_ribose 4TNA
python38 NA_pairwise_interactions.py -c basepair,sugar_ribose 4TNA
python311 NA_pairwise_interactions.py -c basepair,sugar_ribose 4TNA

python38 develop_NA_pairwise_interactions.py
python311 develop_NA_pairwise_interactions.py

"""


# python NA_pairwise_interactions.py -i "C:\Users\zirbel\Documents\FR3D\PDBFiles" -o "C:\Users\zirbel\Documents\FR3D\NAPairwiseInteractions" 4TNA
# python NA_pairwise_interactions.py -i "C:\Users\zirbel\Documents\FR3D\PDBFiles" -o "C:\Users\zirbel\Documents\FR3D\NAPairwiseInteractions" 4TNA.cif.gz
# python NA_pairwise_interactions.py -i "C:\Users\zirbel\Documents\FR3D\PDBFiles" -o "C:\Users\zirbel\Documents\FR3D\NAPairwiseInteractions" -c stacking 4TNA.cif.gz
# python NA_pairwise_interactions.py -i "C:\Users\zirbel\Documents\FR3D\PDBFiles" -o "C:\Users\zirbel\Documents\FR3D\NAPairwiseInteractions" user_test.pdb
# python NA_pairwise_interactions.py -i "C:\Users\zirbel\Documents\FR3D\PDBFiles" -o "C:\Users\zirbel\Documents\FR3D\NAPairwiseInteractions" user_test.pdb.gz

from NA_pairwise_interactions import *
from NA_unit_annotation import generateUnitAnnotation

from fr3d.localpath import outputText
from fr3d.localpath import outputNAPairwiseInteractions
from fr3d.localpath import contact_list_file
from fr3d.localpath import inputPath
from fr3d.localpath import outputHTML
from fr3d.localpath import fr3d_pickle_path

from fr3d.data.base import EntitySelector

#from fr3d.pdb import pdb_reader

parser = argparse.ArgumentParser()
parser.add_argument('-c', "--category", help='Interaction category or categories (basepair,stacking,sO,basepair_detail, bphosphate)')
args = parser.parse_args()
categories = {}

Leontis_Westhof_basepairs = ['cWW', 'cSS', 'cHH', 'cHS', 'cHW', 'cSH', 'cSW', 'cWH', 'cWS', 'tSS', 'tHH', 'tHS', 'tHW', 'tSH', 'tSW', 'tWH', 'tWS', 'tWW']

if args.category:
    for category in args.category.split(","):
        categories[category] = []
else:
    # default is to annotate and write just "true" basepairs
    categories['basepair'] = Leontis_Westhof_basepairs

from hydrogen_bonds import load_ideal_basepair_hydrogen_bonds
from hydrogen_bonds import check_hydrogen_bond

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
PDB_list = ['4V9F']
PDB_list = ['203D']
PDB_list = ['7k00']
PDB_list = ['6CFJ']
PDB_list = ['7K00']
PDB_list = ['7MKY', '4JF2', '5VGW', '4ENC', '5XTM', '2R8S', '4PQV', '3RW6', '4BW0', '6CB3', '4K27', '5U3G', '7OF0', '4LVW', '5D5L', '2NZ4', '3NKB', '6TFG', '2Z75', '4YAZ', '5X2G', '4V9F', '7OZQ', '4Y4O', '4WFL', '1M5K', '7K16', '5FJC', '7O7Y', '6JQ5', '6S0Z', '3P22', '7OX9', '1Q96', '6KWQ', '3LQX', '6U8D', '6SVS', '3E5C', '7RQB', '2EZ6', '6DMC', '2V3C', '5M0I', '3MXH', '4YBB', '5B2P', '4P95', '7KKV', '3NPQ', '5DDP', '4NLF', '7P7Q', '6AZ3', '7D7W', '6S0X', '7RYG', '3AM1', '4PCJ', '5UZ6', '5B2T']
PDB_list += ['5AH5', '4ENC', '7EOG', '1QU2', '2ZUE', '2QUW', '2QUS', '5KPY', '7OF0', '7EQJ', '1U0B', '5AOX', '3FOZ', '2DRA', '4YCO', '7C79', '4V9F', '4Y4O', '4WFL', '3RG5', '5UD5', '7K16', '7O7Y', '6S0Z', '4JXZ', '4J50', '3B31', '3ADD', '7RQB', '3OVB', '6UGG', '4PRF', '4YBB', '3VJR', '1QTQ', '7K98', '4P95', '2GDI', '7DCO', '7P7Q', '6AZ3', '4YYE', '6S0X', '5HR7', '7RYG', '3AM1', '2OEU', '3D2V', '1J1U']
PDB_list = list(set(PDB_list))

PDB_list = ['3AM1','4J50']  # these have symmetry operators, but no annotations there
PDB_list = ['3AM1']  # these have symmetry operators, but no annotations there
PDB_list = ['4J50']  # these have symmetry operators, but no annotations there
PDB_list = ['4RKV','4J50','3AM1']
PDB_list = ['4TNA.cif']
PDB_list = ['5T2A']  # has a conflicting cBW annotation
PDB_list = ['5UED']
PDB_list = ['4V9F','6ZMI','7K00','4TNA']
PDB_list = ['4TNA','7QI4']
PDB_list = ['4K27']
PDB_list = ['5B2R']
PDB_list = ['1AJF', '1JTW', '1N66', '1R7W', '1S9S', '1U6P', '1ZIF', '1ZIG', '1ZIH', '2JXV', '2KPV', '2LK3', '2MXJ', '2MXL', '2N1Q', '4BY9', '6MCI', '6VU1', '6VVJ']
PDB_list = ['2GDI']
PDB_list = ['4J50']
PDB_list = ['3RG5']  # had a problem annotating U|67 basepairs in Python 3.8
PDB_list = ['4M6D']  # to see 4M6D|1|H|G|28  cSH  4M6D|1|H|U|29
PDB_list = ['124D']  # RNA-DNA duplex
PDB_list = ['1A1L']  # DNA-DNA duplex
PDB_list = ['7JQQ']  # Five RNA chains and a DNA-DNA duplex
PDB_list = ['7K00','8B0X']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.280/2.5A/csv']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.280/2.0A/csv']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.280/1.5A/csv']
PDB_list = ['1NBS','6PMO']
PDB_list = ['6DVK']
PDB_list = ['7QI4']
PDB_list = ['7O7Y']
PDB_list = ['1Q96']
PDB_list = ['7QI4','6S0Z','4C8Z','5FQ5','5D5L','5AOX','7K16','5CCB','4V8Y','8D28','4AL5','7OZQ','6YL5','4NLF','7LO9','7QIW','7QIZ','4C8Y','5F4Q','8D28']  # have hydrogens that are named wrong
PDB_list = ['7O7Y','7O7Z','8A3D']
PDB_list = ['1L2X']
PDB_list = ['7ZW0','6S0X','7UVZ']
PDB_list = ['6XU8']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.285/1.5A/csv']

from DNA_2A_list import PDB_list   # define PDB_list as a list of DNA structures

PDB_list = ['4TNA','5KFX','8B0X']
PDB_list = ['4TNA']
PDB_list = ['4V9F']
PDB_list = ['2IZN']   # gets an interaction called p_1??
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.285/3.0A/csv','8B0X','http://rna.bgsu.edu/rna3dhub/nrlist/download/3.285/2.5A/csv','http://rna.bgsu.edu/rna3dhub/nrlist/download/3.285/2.0A/csv','http://rna.bgsu.edu/rna3dhub/nrlist/download/3.285/1.5A/csv']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.339/3.0A/csv','8B0X','http://rna.bgsu.edu/rna3dhub/nrlist/download/3.339/2.5A/csv','http://rna.bgsu.edu/rna3dhub/nrlist/download/3.339/2.0A/csv','http://rna.bgsu.edu/rna3dhub/nrlist/download/3.339/1.5A/csv']

# zzz

OverwriteDataFiles = True    # even if a data file already exists, annotate and overwrite
OverwriteDataFiles = False   # to save time, if a data file exists, skip annotation

base_seq_list = ['A','U','C','G']      # for RNA
base_seq_list = ['DA','DT','DC','DG']  # for DNA
base_seq_list = []                     # for all nucleic acids, modified or not

# tell which types of interactions to annotate
categories = {}
categories['coplanar'] = []   # necessary to get all data for datapoint
categories['basepair'] = []
categories['stacking'] = []
categories['backbone'] = []
categories['sO'] = []        # annotate all sO interactions
categories['sugar_ribose']   = []

ShowStructureReadingErrors = True
ShowStructureReadingErrors = False

experimental = True          # save interactions in pairs_exp folder so they can be compared to ones from the server

# this path should be specified in localpath.py
# intended for writing out a .pickle file to be used by the FR3D motif search tool

# annotate all nucleotides in all chains, even when a representative set is used
annotate_entire_PDB_files = False
annotate_entire_PDB_files = True

timerData = myTimer("start")
lastwritetime = time()

allInteractionDictionary = defaultdict(list)

timerData = myTimer("Making PDB list",timerData)

PDB_IFE_Dict = map_PDB_list_to_PDB_IFE_dict(PDB_list)

print("PDB_IFE_Dict is %s" % PDB_IFE_Dict)

counter = 0
count_pair = 0

# loop through 3D structures and annotate interactions
PDBs = PDB_IFE_Dict.keys()
#PDBs = PDBs[::-1]  # reverse the order of the list, for debugging

print('Annotating these %d PDB files:' % len(PDBs))
print(",".join(sorted(PDBs)))

# If just a few files are requested, overwrite data files
if len(PDBs) > 10 and not OverwriteDataFiles:
    print("Annotating interactions if no file is found in %s" % outputNAPairwiseInteractions)
else:
    print("Annotating interactions and saving in %s" % outputNAPairwiseInteractions)

# restrict dictionary of cutoffs to just the basepairs needed here
Leontis_Westhof_basepairs = ['cWW', 'cSS', 'cHH', 'cHS', 'cHW', 'cSH', 'cSW', 'cWH', 'cWS', 'tSS', 'tHH', 'tHS', 'tHW', 'tSH', 'tSW', 'tWH', 'tWS', 'tWW', 'cWB', 'cBW']
focused_basepair_cutoffs = focus_basepair_cutoffs(nt_nt_cutoffs,Leontis_Westhof_basepairs)
ideal_hydrogen_bonds = load_ideal_basepair_hydrogen_bonds()

PDBs = sorted(PDBs)


# python311 develop_NA_pairwise_interactions.py

worker = [1,len(PDBs),2]     # start at 1 and do every other
worker = [0,len(PDBs),2]     # start at 0 and do every other

worker = [0,len(PDBs),1]     # process each file from 0 to end
worker = [len(PDBs)-1,0,-1]  # start at the end and work backward, one at a time

for i in range(worker[0],worker[1],worker[2]):

    PDB = PDBs[i]

    PDB_id = PDB[0:4]

    counter += 1

    outputDataFileCSV = os.path.join(outputNAPairwiseInteractions, PDB_id + ".csv")

    if not os.path.exists(outputNAPairwiseInteractions):
        os.mkdir(outputNAPairwiseInteractions)

    if experimental:
        outputDataFilePicklePath = os.path.join(fr3d_pickle_path, "pairs_exp")
    else:
        outputDataFilePicklePath = os.path.join(fr3d_pickle_path, "pairs")

    if not os.path.exists(outputDataFilePicklePath):
        os.mkdir(outputDataFilePicklePath)

    outputDataFilePickle = os.path.join(outputDataFilePicklePath, PDB_id + "_RNA_pairs.pickle")

    unit_annotation_file = os.path.join(outputNAPairwiseInteractions,"%s_glycosidic.txt" % PDB_id)
    if not os.path.exists(unit_annotation_file):
        print('Annotating units in %s, which is %d out of %d' % (PDB_id,i+1,len(PDB_IFE_Dict)))
        generateUnitAnnotation(PDB_id, '', inputPath, outputNAPairwiseInteractions, {'glycosidic':[]}, 'txt')


    #print(crashnow)

    if annotate_entire_PDB_files:

        pair_file = "%s_pairs_%s.pickle" % (PDB,fr3d_classification_version)
        pair_to_data_output_file = outputNAPairwiseInteractions + pair_file

        if not os.path.exists(pair_to_data_output_file) or len(PDBs) <= 10 or OverwriteDataFiles:

            print("Reading file %s, which is number %d out of %d" % (PDB,i+1,len(PDB_IFE_Dict)))
            timerData = myTimer("Reading CIF files",timerData)

            structure, messages = load_structure(os.path.join(inputPath,PDB),PDB)
            print(messages)

            if not structure:
                continue

            """
            print('Loading directly with the cif reader')
            rm = read_mode
            filename = os.path.join(inputPath,PDB+".cif.gz")
            if filename.lower().endswith('.cif.gz'):
                with gzip.open(filename, rm) as raw:
                    from fr3d.cif.reader import Cif
                    structure = Cif(raw).structure()
            """

            """
            for base in structure.residues(type = ["RNA linking","DNA linking"]):
                #print(base.unit_id())
                #print(base.centers['glycosidic'])
                #print(base.centers['base'])

                if base.unit_id() in ['1Q96|1|B|A|20','1Q96|1|A|A|9']:
                    print(base.unit_id())
                    print(base.centers['glycosidic'])
                    print(base.centers['base'])

                if base.unit_id() in ['7QI4|1|AA|G|1355']:
                    print('7QI4|1|AA|G|1355 H21 %s' % base.centers['H21'])
                    print('7QI4|1|AA|G|1355 H22 %s' % base.centers['H22'])
            """

            # write out data file of nucleotide centers and rotations that can be used by FR3D for searches
            # need to be able to identify each chain that is available
            write_unit_data_file(PDB,fr3d_pickle_path,structure)

            # annotate interactions and return pair_to_data
            interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_in_structure(structure,categories,focused_basepair_cutoffs,ideal_hydrogen_bonds,[],timerData,True)

            # for pair,data in pair_to_data.items():
            #     print(pair,data)

            # turn this off during development and testing
            if True:
                print("  Annotated these interactions: %s" % interaction_to_list_of_tuples.keys())
                pickle.dump(interaction_to_list_of_tuples,open(outputDataFilePickle,"wb"),2)
                print('  Wrote FR3D pair file %s' % outputDataFilePickle)

            timerData = myTimer("Recording interactions",timerData)
            pickle.dump(pair_to_data,open(pair_to_data_output_file,"wb"),2)
            print('  Wrote classification data file %s' % pair_to_data_output_file)

            write_txt_output_file(outputNAPairwiseInteractions,PDB,interaction_to_list_of_tuples,categories, category_to_interactions)
            print('  Wrote CSV file(s) to %s' % outputNAPairwiseInteractions)

            if len(PDBs) > 10:
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
            structure, messages = load_structure(os.path.join(inputPath,PDB+'.cif'),PDB)
        else:
            # do it this way to suppress error messages
            try:
                structure, messages = load_structure(os.path.join(inputPath,PDB+'.cif'))
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
                if structure:
                    bases = structure.residues(chain = chain_ids)  # load all bases
                else:
                    continue

        # ??? record which RNA/DNA chains are actually present
        # count nucleotides
        numBases = 0
        for base in bases:
            numBases += 1
        print("  Found " + str(numBases) + " bases in " + PDB)

        # build cubes to be able to find potential pairs quickly
        timerData = myTimer("Building cubes",timerData)
        print("  Building nucleotide cubes in " + PDB)
        baseCubeList, baseCubeNeighbors = make_nt_cubes_half(bases, nt_nt_screen_distance, nt_reference_point)

        # annotate nt-nt interactions
        timerData = myTimer("Annotating interactions",timerData)

        # annotate interactions and return pair_to_data
        interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_in_structure(structure,categories,focused_basepair_cutoffs,ideal_hydrogen_bonds,timerData,True)

        # used to return Python_pairs, pair_to_data, timerData

        timerData = myTimer("Recording interactions",timerData)

        # write out pairs in the format that WebFR3D reads
        # accumulate list of interacting units by base, interaction type, and edges
        write_unit_data_file(PDB,fr3d_pickle_path,structure)
        # for nt1, nt2, interaction, edge, standard_aa, param in list_nt_nt:
        #     base = base_residue.unit_id()
        #     # skip symmetry operated instances; generally these are just duplicates anyway
        #     if not "||||" in str(base):
        #         aa = aa_residue.unit_id()
        #         base_component = str(base).split("|")
        #         aa_component = str(aa).split("|")
        #         key = base_component[3]+"_"+aa_component[3]+"_"+interaction+"_"+edge
        #         count_pair += 1
        #         allInteractionDictionary[key].append((base,aa,interaction,edge,standard_aa,param))  # store tuples
            # turn this off during development and testing
        pair_file = "%s_pairs_%s.pickle" % (PDB,fr3d_classification_version)
        pair_to_data_output_file = outputNAPairwiseInteractions + pair_file

        if True:
                print("  Annotated these interactions: %s" % interaction_to_list_of_tuples.keys())
                pickle.dump(interaction_to_list_of_tuples,open(outputDataFilePickle,"wb"),2)
                print('  Wrote FR3D pair file %s' % outputDataFilePickle)

        timerData = myTimer("Recording interactions",timerData)
        pickle.dump(pair_to_data,open(pair_to_data_output_file,"wb"),2)
        print('  Wrote classification data file %s' % pair_to_data_output_file)

        write_txt_output_file(outputNAPairwiseInteractions,PDB,interaction_to_list_of_tuples,categories, category_to_interactions)
        print('  Wrote CSV file(s) to %s' % outputNAPairwiseInteractions)

        if len(PDBs) > 10:
            myTimer("summary",timerData)
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

