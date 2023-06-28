import numpy as np
import os.path
from collections import defaultdict
from time import time
import datetime
import pickle
import urllib
import sys
from fr3d_configuration import DATAPATH
from fr3d_configuration import SERVER

# import the version of urlretrieve appropriate to the Python version
if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
else:
    from urllib.request import urlretrieve as urlretrieve


def writeNAUnitData(structure,DATAPATH,messages=[]):
    """
    accumulate units, centers, rotations and write to files
    according to chain
    """

    chain_to_unit_data = {}
    chains = []

    nts = structure.residues(type = ["RNA linking","DNA linking"])

    for nt in nts:
        unit_id = nt.unit_id()
        fields = unit_id.split('|')
        chain = '|'.join(fields[0:3])

        if not chain in chain_to_unit_data:
            chain_to_unit_data[chain] = {}
            chain_to_unit_data[chain]['ids'] = []
            chain_to_unit_data[chain]['chainIndices'] = []
            chain_to_unit_data[chain]['centers'] = []
            chain_to_unit_data[chain]['rotations'] = []

        chain_to_unit_data[chain]['ids'].append(nt.unit_id())
        chain_to_unit_data[chain]['chainIndices'].append(nt.index)
        chain_to_unit_data[chain]['centers'].append(nt.centers['glycosidic'])
        chain_to_unit_data[chain]['rotations'].append(nt.rotation_matrix)

    # write chain data out to individual files
    # this code writes base centers, when we write glycosidic centers, name output file _NA.pickle
    for chain in chain_to_unit_data.keys():
        unit_info = [chain_to_unit_data[chain]['ids']]
        unit_info.append(chain_to_unit_data[chain]['chainIndices'])
        unit_info.append(chain_to_unit_data[chain]['centers'])
        unit_info.append(chain_to_unit_data[chain]['rotations'])

        chain_filename = os.path.join(DATAPATH,'units',chain.replace('|','-')+'_RNA.pickle')

        with open(chain_filename, 'wb') as fh:
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(unit_info, fh, 2)

        chains.append(chain)
        messages.append('Wrote %d nucleotides to %s' % (len(chain_to_unit_data[chain]['ids']),chain_filename))

    return chains, messages


def writeNAPairwiseInteractions(structure,DATAPATH,messages=[]):
    """
    Annotate pairwise interactions and write one file for the
    structure.
    """

    from fr3d.classifiers.NA_pairwise_interactions import annotate_nt_nt_in_structure

    # tell which categories of interactions to annotate
    categories = {}
    categories['basepair'] = []
    categories['stacking'] = []
    categories['sO'] = []

    interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_in_structure(structure,categories)

    print("  file_reading:  structure name: %s" % structure.pdb)

    # use RNA_pairs for now until _NA_pairs is established
    pdb_filename = os.path.join(DATAPATH,'pairs',structure.pdb+'_RNA_pairs.pickle')
    with open(pdb_filename, 'wb') as fh:
        # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
        pickle.dump(interaction_to_list_of_tuples, fh, 2)

    return messages


def writeNAUnitAnnotations(structure,DATAPATH,messages=[]):
    """
    Annotate chi angle and glycosidic bond conformation.
    """

    print("  file_reading:  annotating bond orientation for %s" % structure.pdb)

    bond_annotations, error_messages = annotate_bond_orientation(structure,False)

    messages += error_messages

    chain_to_unit_data = {}

    # organize annotations by chain
    for annotation in bond_annotations:
        unit_id = annotation['unit_id']
        fields = unit_id.split('|')
        chain = '|'.join(fields[0:3])

        if not chain in chain_to_unit_data:
            chain_to_unit_data[chain] = {}

        if not unit_id in chain_to_unit_data[chain]:
            chain_to_unit_data[chain][unit_id] = {}

        chain_to_unit_data[chain][unit_id]['chi_degree'] = annotation['chi_degree']
        chain_to_unit_data[chain][unit_id]['orientation'] = annotation['orientation']

    # write chain data out to individual files for each chain
    for chain in chain_to_unit_data.keys():
        chain_filename = os.path.join(DATAPATH,'units',chain.replace('|','-')+'_NA_unit_annotations.pickle')

        with open(chain_filename, 'wb') as fh:
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(chain_to_unit_data[chain], fh, 2)

        messages.append('Wrote %d nucleotides to %s' % (len(chain_to_unit_data[chain]),chain_filename))

    return messages


def writeNABackboneData(structure,DATAPATH,messages=[]):
    """
    accumulate units, phosphate, and sugar centers and write to files
    according to chain
    """

    chain_to_unit_data = {}
    chains = []

    nts = structure.residues(type = ["RNA linking","DNA linking"])

    for nt in nts:
        unit_id = nt.unit_id()
        fields = unit_id.split('|')
        chain = '|'.join(fields[0:3])

        if not chain in chain_to_unit_data:
            chain_to_unit_data[chain] = {}
            chain_to_unit_data[chain]['ids'] = []
            chain_to_unit_data[chain]['chainIndices'] = []
            chain_to_unit_data[chain]['phosphate'] = []
            chain_to_unit_data[chain]['sugar'] = []

        chain_to_unit_data[chain]['ids'].append(nt.unit_id())
        chain_to_unit_data[chain]['chainIndices'].append(nt.index)
        chain_to_unit_data[chain]['phosphate'].append(nt.centers['nt_phosphate'])
        chain_to_unit_data[chain]['sugar'].append(nt.centers['nt_sugar'])

    # write chain data out to individual files
    for chain in chain_to_unit_data.keys():
        chain_filename = os.path.join(DATAPATH,'units',chain.replace('|','-')+'_NA_phosphate_sugar.pickle')

        unit_info = [chain_to_unit_data[chain]['ids']]
        unit_info.append(chain_to_unit_data[chain]['chainIndices'])
        unit_info.append(chain_to_unit_data[chain]['phosphate'])
        unit_info.append(chain_to_unit_data[chain]['sugar'])

        with open(chain_filename, 'wb') as fh:
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(unit_info, fh, 2)

        chains.append(chain)
        messages.append('Wrote %d nucleotides to %s' % (len(chain_to_unit_data[chain]['ids']),chain_filename))

    return messages


def writeProteinUnitData(structure,DATAPATH,messages=[]):
    """
    accumulate units, indices, centers and write to one file
    """

    aas = structure.residues(type = ["L-peptide linking"])

    aa_unit_ids = []
    aa_indices  = []
    aa_centers  = []
    aa_backbones = []

    for aa in aas:
        aa_unit_ids.append(aa.unit_id())
        aa_indices.append(aa.index)
        aa_centers.append(aa.centers['aa_fg'])  # amino acid functional group center, None in GLY
        aa_backbones.append(aa.centers['aa_backbone'])  # amino acid backbone, average of ['N','CA','C','O']



        if 'GLY' in aa.unit_id():
            print('file_reading',aa.unit_id,aa.index,aa.centers['aa_fg'],aa.centers['aa_backbone'])



    # write chain data out to individual files
    protein_filename = os.path.join(DATAPATH,'units',structure.pdb+'_protein.pickle')

    with open(protein_filename, 'wb') as fh:
        # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
        pickle.dump([aa_unit_ids,aa_indices,aa_centers,aa_backbones], fh, 2)

    messages.append('Wrote %d amino acids to %s' % (len(aa_unit_ids),protein_filename))

    return messages


def processPDBFile(PDBname):
    """
    If we have to start with a .pdb or .cif file,
    load that file and write out all necessary .pickle files
    with required annotations.
    """

    from fr3d_configuration import CIFPATH
    from fr3d.classifiers.NA_pairwise_interactions import load_structure
    from fr3d.classifiers.NA_pairwise_interactions import annotate_nt_nt_in_structure
    from fr3d.classifiers.NA_unit_annotation import annotate_bond_orientation


    messages = []

    if '.cif' in PDBname.lower() or '.pdb' in PDBname.lower():
        cifPathAndFileName = PDBname
    else:
        fields = PDBname.split('|')
        cifPathAndFileName = os.path.join(CIFPATH,fields[0]+'.cif')

    # read the coordinate file
    structure, messages = load_structure(cifPathAndFileName)

    # write out centers and rotation matrices of nucleotides by chain
    chains, messages = writeNAUnitData(structure,DATAPATH,messages)

    # annotate pairwise interactions
    messages = writeNAPairwiseInteractions(structure,DATAPATH,messages)

    # annotate chi angle and glycosidic bond conformation
    messages = writeNAUnitAnnotations(structure,DATAPATH,messages)

    # write out centers and indices for amino acids for the whole structure
    messages = writeProteinUnitData(structure,DATAPATH,messages)

    # write out centers and indices for amino acids for the whole structure
    messages = writeNABackboneData(structure,DATAPATH,messages)

    for message in messages:
        print("  %s" % message)

    return chains, messages

def readPDBDatafile():
    """
    Read .pickle file containing data about each nucleic-
    acid-containing PDB file.
    """

    datafile = {}

    filename = "NA_datafile.pickle"
    pathAndFileName = os.path.join(DATAPATH,'units',filename)

    #if not os.path.exists(pathAndFileName) and not SERVER:
    # try to update the file for a local installation, but not too often
    if not SERVER:
        download_needed = True

        if os.path.exists(pathAndFileName):
            m_time = os.path.getmtime(pathAndFileName)
            dt_m = datetime.datetime.fromtimestamp(m_time)
            elapsed_days = (datetime.datetime.now() - dt_m).total_seconds()/(60*60*24.0)
            #print('%s was modified on %s' % (pathAndFileName,dt_m))
            #print('elapsed days is %0.4f' % elapsed_days)
            if elapsed_days < 7:
                download_needed = False

        if download_needed:
            try:
                print("Downloading %s" % filename)
                urlretrieve("http://rna.bgsu.edu/units/" + filename, pathAndFileName)
            except:
                print("Unable to download %s" % filename)

    if os.path.exists(pathAndFileName):
        try:
            if sys.version_info[0] < 3:
                datafile = pickle.load(open(pathAndFileName,"rb"))
            else:
                datafile = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
        except:
            print("Could not read "+filename)

    return datafile


def readNAPositionsFile(Q, chainString, starting_index):
    """
    Read .pickle file of RNA/DNA base center and rotation matrix;
    download if necessary; create if necessary
    """

    ids = []
    chainIndices = []
    centers = []
    rotations = []
    line_num = starting_index

    centers = np.zeros((0,3))
    rotations = np.zeros((0,3,3))
    id_to_index = defaultdict()
    index_to_id = defaultdict()

    filename = chainString + '_RNA.pickle'  # old, uses base center
    filename = chainString + '_NA.pickle'   # new on 2023-02-20, covers RNA and DNA, uses glycosidic atom
    pathAndFileName = os.path.join(DATAPATH,'units',filename)

    if not os.path.exists(pathAndFileName) and not SERVER:
        try:
            urlretrieve("http://rna.bgsu.edu/units/" + filename, pathAndFileName) # testing
            print("Downloaded "+filename)
        except:
            print("file_reading: Could not download %s, reading .cif" % filename)

            chains, messages = processPDBFile(chainString)
            Q['userMessage'] += messages

    if os.path.exists(pathAndFileName):
        try:
            if sys.version_info[0] < 3:
                ids, chainIndices, centers, rotations = pickle.load(open(pathAndFileName,"rb"))
            else:
                ids, chainIndices, centers, rotations = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
        except:
            print("Could not read "+filename)
            Q["errorMessage"].append("Could not retrieve RNA unit file "+filename)

        for i in range(0,len(ids)):
            id_to_index[ids[i]] = line_num
            index_to_id[line_num] = ids[i]
            line_num += 1

    else:
        print("Could not find "+filename)
        Q["errorMessage"].append("Could not retrieve RNA unit file "+filename)

    centers = np.asarray(centers)
    rotations = np.asarray(rotations)

    return Q, centers, rotations, ids, id_to_index, index_to_id, chainIndices

def readNABackboneFile(Q, chainString, starting_index):
    """
    Read .pickle file of RNA phosphate and sugar centers; download if necessary
    """

    ids = []
    chainIndices = []
    centers = []
    rotations = []
    line_num = starting_index

    phosphate = np.zeros((0,3))
    sugar = np.zeros((0,3))
    id_to_index = defaultdict()
    index_to_id = defaultdict()

    filename = chainString + "_NA_phosphate_sugar.pickle"
    pathAndFileName = os.path.join(DATAPATH,'units',filename)

    if not os.path.exists(pathAndFileName) and not SERVER:
        urlretrieve("http://rna.bgsu.edu/units/"+filename, pathAndFileName)
        print("Downloaded "+filename)

    if os.path.exists(pathAndFileName):
        if sys.version_info[0] < 3:
            try:
                ids, chainIndices, phosphate, sugar = pickle.load(open(pathAndFileName,"rb"))
            except:
                print("Could not read "+filename)
                Q["userMessage"].append("Could not retrieve RNA unit file "+filename)
        else:
            try:
                ids, chainIndices, phosphate, sugar = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
            except:
                print("Could not read "+filename)
                Q["userMessage"].append("Could not retrieve RNA unit file "+filename)

        for i in range(0,len(ids)):
            id_to_index[ids[i]] = line_num
            index_to_id[line_num] = ids[i]
            line_num += 1

    else:
        print("Could not find "+filename)
        Q["userMessage"].append("Could not retrieve RNA unit file "+filename)

    phosphate = np.asarray(phosphate)
    sugar     = np.asarray(sugar)

    return Q, phosphate, sugar, ids, id_to_index, index_to_id, chainIndices

def readProteinPositionsFile(Q, PDBID, starting_index):

    ids = []
    chainIndices = []
    centers = []
    line_num = starting_index

    centers = np.zeros((0,3))  # of lines/nucleotides in file
    rotations = np.zeros((0,3,3))
    id_to_index = defaultdict()
    index_to_id = defaultdict()

    filename = PDBID + "_protein.pickle"
    pathAndFileName = os.path.join(DATAPATH,'units',filename)

    if not os.path.exists(pathAndFileName) and not SERVER:
        urlretrieve("http://rna.bgsu.edu/units/"+filename, pathAndFileName)
        print("Downloaded "+filename)

    if os.path.exists(pathAndFileName):
        item_list = []
        if sys.version_info[0] < 3:
            try:
                item_list = pickle.load(open(pathAndFileName,"rb"))
            except:
                print("Could not read "+filename)
                Q["userMessage"].append("Could not retrieve protein unit file "+filename)
                return Q, centers, ids, id_to_index, index_to_id, chainIndices
        else:
            try:
                item_list = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
            except:
                print("Could not read "+filename)
                Q["userMessage"].append("Could not retrieve protein unit file "+filename)
                return Q, centers, ids, id_to_index, index_to_id, chainIndices

        if len(item_list) == 3:
            ids, chainIndices, centers = item_list
        elif len(item_list) == 4:
            ids, chainIndices, centers, backbones = item_list
        else:
            Q["userMessage"].append("Could not retrieve protein unit file "+filename)
            return Q, centers, ids, id_to_index, index_to_id, chainIndices

        for i in range(0,len(ids)):
            id_to_index[ids[i]] = line_num
            index_to_id[line_num] = ids[i]
            line_num += 1

    else:
        print("Could not find "+filename)
        Q["userMessage"].append("Could not retrieve protein unit file "+filename)
        return Q, centers, ids, id_to_index, index_to_id, chainIndices

    centers = np.asarray(centers)

    return Q, centers, ids, id_to_index, index_to_id, chainIndices


def readNAPairsFileRaw(Q, PDBID, alternate = ""):
    """
    Read the .pickle file that stores all pairwise interactions for a PDB id.
    Alternate could be "_exp" for experimental annotations, for example.
    """

    interactionToTriples = defaultdict(list)

    pairsFileName = PDBID + '_RNA_pairs.pickle'
    pathAndFileName = os.path.join(DATAPATH,'pairs'+alternate,pairsFileName)

    justDownloaded = False

    if not os.path.exists(pathAndFileName) and not SERVER:
        url = "http://rna.bgsu.edu/pairs" + alternate + "/"+pairsFileName
        urlretrieve(url, pathAndFileName) # testing
        print("Downloaded "+pairsFileName)

        justDownloaded = True

    if os.path.exists(pathAndFileName):
        if sys.version_info[0] < 3:
            try:
                interactionToTriples = pickle.load(open(pathAndFileName,"rb"))
            except:
                print("Could not read "+pairsFileName+", it may not be available")
                Q["userMessage"].append("Could not retrieve RNA pairs file "+pairsFileName)
                if justDownloaded:
                    os.remove(pathAndFileName)
        else:
            try:
                interactionToTriples = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
            except:
                print("Could not read "+pairsFileName+", it may not be available D*********************************")
                Q["userMessage"].append("Could not retrieve RNA pairs file "+pairsFileName)
                if justDownloaded:
                    os.remove(pathAndFileName)

    return Q, interactionToTriples


def readNAPairsFile(Q, PDBID, id_to_index, alternate = ""):
    """
    We need to read the pairwise interactions whether or not there are
    constraints, to be able to display the interactions in the output.
    alternate is for experimental or alternate lists of interacting pairs.
    The alternate .pickle files must be stored in a parallel directory.
    """

    Q, interactionToTriples = readNAPairsFileRaw(Q, PDBID, alternate)

    # interactionToTriples has this structure:
    # interactionToTriples['cWW'] is a list of triples, each triple being (unit_id_1,unit_id_2,range)
    # However, what is passed back is interactionToIndexPairs, which has a different structure
    # interactionToIndexPairs['cWW'][0] is the list of pairs
    # interactionToIndexPairs['cWW'][1] is the list of crossing numbers

    # convert pairs of unit ids to pairs of indices for use in FR3D.py
    # add alternate tag to interactions from alternate files

    interactionToIndexPairs = {}
    pairToInteractions = defaultdict(list)
    pairToCrossingNumber = {}

    for interaction in interactionToTriples:
        if interaction+alternate in Q["activeInteractions"]:
            newListOfPairs = []
            crossingNumber = []
            if "BPh" in interaction or "BR" in interaction:
                newListOfPairsSelf = []                      # 0BPh and maybe others can make self interactions
                crossingNumberSelf = []
                for pair in interactionToTriples[interaction]:
                    if pair[0] in id_to_index and pair[1] in id_to_index:
                        indexPair = (id_to_index[pair[0]],id_to_index[pair[1]])
                        pairToInteractions[indexPair].append(interaction+alternate)
                        pairToCrossingNumber[indexPair] = pair[2]
                        if pair[0] == pair[1]:
                            newListOfPairsSelf.append(indexPair)
                            crossingNumberSelf.append(pair[2])
                        else:
                            newListOfPairs.append(indexPair)
                            crossingNumber.append(pair[2])
                interactionToIndexPairs[interaction+alternate] = (newListOfPairs,crossingNumber)
                interactionToIndexPairs[interaction+alternate+"_self"] = (newListOfPairsSelf,crossingNumberSelf)
            else:
                for pair in interactionToTriples[interaction]:
                    if pair[0] in id_to_index and pair[1] in id_to_index:
                        indexPair = (id_to_index[pair[0]],id_to_index[pair[1]])
                        pairToInteractions[indexPair].append(interaction+alternate)
                        pairToCrossingNumber[indexPair] = pair[2]
                        newListOfPairs.append(indexPair)
                        crossingNumber.append(pair[2])
                interactionToIndexPairs[interaction+alternate] = (newListOfPairs,crossingNumber)
        else:
            for pair in interactionToTriples[interaction]:
                if pair[0] in id_to_index and pair[1] in id_to_index:
                    indexPair = (id_to_index[pair[0]],id_to_index[pair[1]])
                    pairToInteractions[indexPair].append(interaction+alternate)
                    pairToCrossingNumber[indexPair] = pair[2]

    return Q, interactionToIndexPairs, pairToInteractions, pairToCrossingNumber


def readUnitAnnotations(Q, ifename):
    # split up ifename and read unit annotations for each individual chain
    # unit annotations are a dictionary mapping unit id to a tuple t
    # t[0] is the chi angle in degrees
    # t[1] is anti/syn/intermediate syn
    # When more annotations are added, they can be added to the tuple

    chains =  ifename.split('+')

    unit_id_to_annotation = {}

    for chain in chains:
        filename        = "%s_NA_unit_annotations.pickle" % chain.replace("|","-")
        pathAndFileName = os.path.join(DATAPATH,'units',filename)

        if not os.path.exists(pathAndFileName) and not SERVER:
            urlretrieve("http://rna.bgsu.edu/units/" + filename, pathAndFileName)
            print("Downloaded "+filename)

        # note:  if the file is not present on the server, a text file with a 404 error will be downloaded
        # so it will look like a file was downloaded, but it's not the file you need!

        if os.path.exists(pathAndFileName):
            if sys.version_info[0] < 3:
                try:
                    new_mapping = pickle.load(open(pathAndFileName,"rb"))
                    unit_id_to_annotation.update(new_mapping)
                except:
                    print("Could not read "+filename+", it may not be available D=================================")
                    Q["userMessage"].append("Could not retrieve RNA pairs file "+filename)
                    try:
                        os.remove(pathAndFileName)
                    except:
                        Q["userMessage"].append("Could not remove RNA pairs file "+filename)
            else:
                try:
                    new_mapping = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
                    unit_id_to_annotation.update(new_mapping)
                except:
                    print("Could not read "+filename+", it may not be available D*********************************")
                    Q["userMessage"].append("Could not retrieve RNA pairs file "+filename)
                    try:
                        os.remove(pathAndFileName)
                    except:
                        Q["userMessage"].append("Could not remove RNA pairs file "+filename)


    return Q, unit_id_to_annotation
