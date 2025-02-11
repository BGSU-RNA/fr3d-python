# This script is for developing and testing NA_pairwise_interactions.py

"""
When changes are made to other code in fr3d-python, as administrator:
cd c:/Users/zirbel/Documents/GitHub/fr3d-python
python311 -m pip install .
python27 -m pip install .
python38 -m pip install .

To run the code:
cd c:/Users/zirbel/Documents/GitHub/fr3d-python/fr3d/classifiers
python27 NA_pairwise_interactions.py -c basepair,stacking,sugar_ribose 4TNA
python38 NA_pairwise_interactions.py -c basepair,sugar_ribose 4TNA
python311 NA_pairwise_interactions.py -c basepair,sugar_ribose 4TNA

python38 develop_NA_pairwise_interactions.py
python311 develop_NA_pairwise_interactions.py 1
python311 develop_NA_pairwise_interactions.py 2
python311 develop_NA_pairwise_interactions.py 4
python311 develop_NA_pairwise_interactions.py 5

"""


# python NA_pairwise_interactions.py -i "C:\Users\zirbel\Documents\FR3D\PDBFiles" -o "C:\Users\zirbel\Documents\FR3D\NAPairwiseInteractions" 4TNA
# python NA_pairwise_interactions.py -i "C:\Users\zirbel\Documents\FR3D\PDBFiles" -o "C:\Users\zirbel\Documents\FR3D\NAPairwiseInteractions" 4TNA.cif.gz
# python NA_pairwise_interactions.py -i "C:\Users\zirbel\Documents\FR3D\PDBFiles" -o "C:\Users\zirbel\Documents\FR3D\NAPairwiseInteractions" -c stacking 4TNA.cif.gz
# python NA_pairwise_interactions.py -i "C:\Users\zirbel\Documents\FR3D\PDBFiles" -o "C:\Users\zirbel\Documents\FR3D\NAPairwiseInteractions" user_test.pdb
# python NA_pairwise_interactions.py -i "C:\Users\zirbel\Documents\FR3D\PDBFiles" -o "C:\Users\zirbel\Documents\FR3D\NAPairwiseInteractions" user_test.pdb.gz

import os
import pickle

from NA_pairwise_interactions import *
from NA_unit_annotation import generateUnitAnnotation

from fr3d.localpath import outputNAPairwiseInteractions
from fr3d.localpath import inputPath
from fr3d.localpath import fr3d_pickle_path

from hydrogen_bonds import load_ideal_basepair_hydrogen_bonds


parser = argparse.ArgumentParser()
parser.add_argument('worker', type=str, nargs='+', help='0 for all, 1 to process evens, 2 to process odds, 3 to start at end')
parser.add_argument('-c', "--category", help='Interaction category or categories (basepair,stacking,sO,basepair_detail, bphosphate)')
args = parser.parse_args()

Leontis_Westhof_basepairs = ['cWW', 'cSS', 'cHH', 'cHS', 'cHW', 'cSH', 'cSW', 'cWH', 'cWS', 'tSS', 'tHH', 'tHS', 'tHW', 'tSH', 'tSW', 'tWH', 'tWS', 'tWW']

if args.category:
    categories = {}
    for category in args.category.split(","):
        categories[category] = []
else:
    # default is to annotate and write just "true" basepairs
    categories = {}
    categories['basepair'] = Leontis_Westhof_basepairs
    # tell which types of interactions to annotate
    categories['coplanar'] = []   # necessary to get all data for datapoint
    categories['basepair'] = []
    categories['basepair_detail'] = []
    # categories['stacking'] = []
    # categories['backbone'] = []
    # categories['sO'] = []        # annotate all sO interactions
    # categories['sugar_ribose']   = []

if args.worker:
    worker = int(args.worker[0])
else:
    worker = 0

PDB_list = ['4V9F']

from DNA_2A_list import PDB_list   # define PDB_list as a list of DNA structures

PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.308/3.0A/csv']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/3.285/1.5A/csv']
PDB_list = ['4V9F','6AZ3','6GYV','7O7Y','7OYC','7QI4','7QIW','7V9E','8A98','8AZW','8GLP','5J7L','7RQB']
PDB_list = ['http://rna.bgsu.edu/rna3dhub/nrlist/download/NR/3.349/3.0A/csv','8B0X','8GLP','http://rna.bgsu.edu/rna3dhub/nrlist/download/NR/3.349/2.5A/csv','http://rna.bgsu.edu/rna3dhub/nrlist/download/NR/3.349/2.0A/csv','http://rna.bgsu.edu/rna3dhub/nrlist/download/NR/3.349/1.5A/csv']

if True:
    # read chains from datmos / nabir
    PDB_set = set()
    PDB_chain_set = set()
    for mt in ['RNA','DNA']:
        filename = '%s_reference_chains.csv' % mt
        path_filename = os.path.join('C:/Users/zirbel/Documents/FR3D/Python FR3D/data/pairs_datmos',filename)
        with open(path_filename,'rt') as f:
            lines = f.readlines()
        for line in lines:
            pdb,chain,desc,count = line.split(",")
            PDB_set.add(pdb.upper())
            PDB_chain_set.add(pdb.upper()+"|1|"+chain)
    PDB_list = list(PDB_set)

# zzz

OverwriteDataFiles = True    # even if a data file already exists, annotate and overwrite
OverwriteDataFiles = False   # to save time, if a data file exists, skip annotation

base_seq_list = ['A','U','C','G']      # for RNA
base_seq_list = ['DA','DT','DC','DG']  # for DNA
base_seq_list = []                     # for all nucleic acids, modified or not

ShowStructureReadingErrors = True
ShowStructureReadingErrors = False

experimental = True          # save interactions in pairs_exp folder so they can be compared to ones from the server
experimental = False

# this path should be specified in localpath.py
# intended for writing out a .pickle file to be used by the FR3D motif search tool

# annotate all nucleotides in all chains, even when a representative set is used
annotate_entire_PDB_files = False
annotate_entire_PDB_files = True

annotate_units = False

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

# simple parallelization
if worker == 0:     # start at 0 and process all files
    a = 0
    b = len(PDBs)
    c = 1
elif worker == 1:   # start at 0 and process even-numbered files
    a = 0
    b = len(PDBs)
    c = 2
elif worker == 2:   # start at 1 and process odd-numbered files
    a = 1
    b = len(PDBs)
    c = 2
elif worker == 3:   # start at the end and process all files
    a = len(PDBs)-1
    b = 0
    c = -1
elif worker == 4:   # start at the end and process every other
    a = len(PDBs)-1
    b = 0
    c = -2
elif worker == 5:   # start almost at the end and process every other
    a = len(PDBs)-2
    b = 0
    c = -2
else:
    a = worker      # start at indicated number and process all files
    b = len(PDBs)
    c = 1

for i in range(a,b,c):

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


    if annotate_units:
        unit_annotation_file = os.path.join(outputNAPairwiseInteractions,"%s_glycosidic.txt" % PDB_id)
        if not os.path.exists(unit_annotation_file):
            print('Annotating units in %s, which is %d out of %d' % (PDB_id,i+1,len(PDB_IFE_Dict)))
            generateUnitAnnotation(PDB_id, '', inputPath, outputNAPairwiseInteractions, {'glycosidic':[]}, 'txt')

    outputDataFilePickle = os.path.join(outputDataFilePicklePath, PDB_id + "_RNA_pairs.pickle")

    if annotate_entire_PDB_files:
        # name for file with pairs and datapoint variable about annotations
        pair_file = "%s_datapoint.pickle" % (PDB)
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
            # write_unit_data_file(PDB,fr3d_pickle_path,structure)

            # annotate interactions and return pair_to_data
            interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_in_structure(structure,categories,focused_basepair_cutoffs,ideal_hydrogen_bonds,[],timerData,True)

            # for pair,data in pair_to_data.items():
            #     print(pair,data)

            # turn this off during development and testing
            if False:
                print("  Annotated these interactions: %s" % interaction_to_list_of_tuples.keys())
                pickle.dump(interaction_to_list_of_tuples,open(outputDataFilePickle,"wb"),2)
                print('  Wrote FR3D pair file %s' % outputDataFilePickle)


            timerData = myTimer("Recording interactions",timerData)
            pickle.dump(pair_to_data,open(pair_to_data_output_file,"wb"),5)
            print('  Wrote classification data file %s' % pair_to_data_output_file)

            # write_txt_output_file(outputNAPairwiseInteractions,PDB,interaction_to_list_of_tuples,categories, category_to_interactions)
            # print('  Wrote CSV file(s) to %s' % outputNAPairwiseInteractions)

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
        # write_unit_data_file(PDB,fr3d_pickle_path,structure)

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
        pair_file = "%s_pairs.pickle" % (PDB)
        pair_to_data_output_file = outputNAPairwiseInteractions + pair_file

        if False:
            print("  Annotated these interactions: %s" % interaction_to_list_of_tuples.keys())
            pickle.dump(interaction_to_list_of_tuples,open(outputDataFilePickle,"wb"),2)
            print('  Wrote FR3D pair file %s' % outputDataFilePickle)

        timerData = myTimer("Recording interactions",timerData)
        pickle.dump(pair_to_data,open(pair_to_data_output_file,"wb"),5)
        print('  Wrote classification data file %s' % pair_to_data_output_file)

        # write_txt_output_file(outputNAPairwiseInteractions,PDB,interaction_to_list_of_tuples,categories, category_to_interactions)
        # print('  Wrote CSV file(s) to %s' % outputNAPairwiseInteractions)

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

