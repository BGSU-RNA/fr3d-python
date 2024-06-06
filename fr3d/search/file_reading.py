
from collections import defaultdict
import datetime
import numpy as np
import os.path
import pickle
import sys
from time import time

from fr3d.search.fr3d_configuration import SERVER

# import the version of urlretrieve appropriate to the Python version
if sys.version_info[0] < 3:
    from urllib import urlretrieve as urlretrieve
else:
    from urllib.request import urlretrieve as urlretrieve


def checkDirectories(Q):
    """
    check for required directories for data and output
    """

    if SERVER or Q.get('server', False):
        # on the server, the directories are already set up; save time
        return Q

    if "DATAPATHUNITS" in Q:
        DATAPATHUNITS = Q["DATAPATHUNITS"]
    else:
        try:
            from fr3d_configuration import DATAPATHUNITS
            Q["DATAPATHUNITS"] = DATAPATHUNITS
        except:
            print("Error: Could not find DATAPATHUNITS in query or in fr3d_configuration.py")
            Q["errorMessage"].append("Error: Could not find DATAPATHUNITS in query or in fr3d_configuration.py")
            Q["errorStatus"] = "write and exit"
            return Q

    if "DATAPATHPAIRS" in Q:
        DATAPATHPAIRS = Q["DATAPATHPAIRS"]
    else:
        try:
            from fr3d_configuration import DATAPATHPAIRS
            Q["DATAPATHPAIRS"] = DATAPATHPAIRS
        except:
            print("Error: Could not find DATAPATHPAIRS in query or in fr3d_configuration.py")
            Q["errorMessage"].append("Error: Could not find DATAPATHPAIRS in query or in fr3d_configuration.py")
            Q["errorStatus"] = "write and exit"
            return Q

    if "OUTPUTPATH" in Q:
        OUTPUTPATH = Q["OUTPUTPATH"]
    else:
        try:
            from fr3d_configuration import OUTPUTPATH
            Q["OUTPUTPATH"] = OUTPUTPATH
        except:
            print("Error: Could not find OUTPUTPATH in query or in fr3d_configuration.py")
            Q["errorMessage"].append("Error: Could not find OUTPUTPATH in query or in fr3d_configuration.py")
            Q["errorStatus"] = "write and exit"
            return Q

    for directory in [DATAPATHUNITS, DATAPATHPAIRS, OUTPUTPATH]:
        try:
            os.stat(directory)
        except:
            try:
                os.mkdir(directory)
                print("Made " + directory + " directory")
            except:
                print("Error: Could not make " + directory + " directory")
                Q["errorMessage"].append("Error: Could not make " + directory + " directory")
                Q["errorStatus"] = "write and exit"
                return Q
    return Q


def get_DATAPATHUNITS(Q):
    """
    Figure out where to look for files containing information on individual units,
    the centers and rotation matrices.
    """

    if "DATAPATHUNITS" in Q:
        DATAPATHUNITS = Q["DATAPATHUNITS"]
    else:
        try:
            from fr3d_configuration import DATAPATHUNITS
        except:
            print("Error: Could not find DATAPATHUNITS in query or in fr3d_configuration.py")
            Q["errorMessage"].append("Error: Could not find DATAPATHUNITS in query or in fr3d_configuration.py")
            Q["errorStatus"] = "write and exit"

    directory = DATAPATHUNITS
    try:
        os.stat(directory)
    except:
        try:
            os.mkdir(directory)
            print("Made " + directory + " directory")
        except:
            print("Error: Could not make " + directory + " directory")
            Q["errorMessage"].append("Error: Could not make " + directory + " directory")
            Q["errorStatus"] = "write and exit"
            return Q

    return Q, DATAPATHUNITS


def get_DATAPATHPAIRS(Q):
    """
    Figure out where to look for files containing information on individual units,
    the centers and rotation matrices.
    """

    if "DATAPATHPAIRS" in Q:
        DATAPATHPAIRS = Q["DATAPATHPAIRS"]
    else:
        try:
            from fr3d_configuration import DATAPATHPAIRS
        except:
            print("Error: Could not find DATAPATHPAIRS in query or in fr3d_configuration.py")
            Q["errorMessage"].append("Error: Could not find DATAPATHPAIRS in query or in fr3d_configuration.py")
            Q["errorStatus"] = "write and exit"

    directory = DATAPATHPAIRS
    try:
        os.stat(directory)
    except:
        try:
            os.mkdir(directory)
            print("Made " + directory + " directory")
        except:
            print("Error: Could not make " + directory + " directory")
            Q["errorMessage"].append("Error: Could not make " + directory + " directory")
            Q["errorStatus"] = "write and exit"
            return Q

    return Q, DATAPATHPAIRS


def get_CIFPATH(Q):
    """
    When working with user-defined queries or when one wants to not download pre-computed
    annotations, we need a place to store the .cif or .pdb files.
    Check that this directory is defined and that the directory exists.
    """

    if "CIFPATH" in Q:
        CIFPATH = Q["CIFPATH"]
    else:
        try:
            from fr3d_configuration import CIFPATH
            Q["CIFPATH"] = CIFPATH
        except:
            print("Error: Could not find CIFPATH in query or in fr3d_configuration.py")
            Q["errorMessage"].append("Error: Could not find CIFPATH in query or in fr3d_configuration.py")
            Q["errorStatus"] = "write and exit"
            return Q

    directory = CIFPATH
    try:
        os.stat(directory)
    except:
        try:
            os.mkdir(directory)
            print("Made " + directory + " directory")
        except:
            print("Error: Could not make " + directory + " directory")
            Q["errorMessage"].append("Error: Could not make " + directory + " directory")
            Q["errorStatus"] = "write and exit"
            return Q

    return Q, CIFPATH


def writeNAUnitData(structure,DATAPATHUNITS,messages=[]):
    """
    accumulate units, centers, rotations and write to files according to chain
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
        chain_to_unit_data[chain]['centers'].append(nt.centers['glycosidic'])  # use glycosidic atom
        chain_to_unit_data[chain]['rotations'].append(nt.rotation_matrix)

    # write chain data out to individual files
    # this code writes base centers, when we write glycosidic centers, name output file _NA.pickle
    for chain in chain_to_unit_data.keys():
        unit_info = [chain_to_unit_data[chain]['ids']]
        unit_info.append(chain_to_unit_data[chain]['chainIndices'])
        unit_info.append(chain_to_unit_data[chain]['centers'])
        unit_info.append(chain_to_unit_data[chain]['rotations'])

        # follow format used on BGSU RNA server like 4V9F-1-0_NA.pickle
        chain_filename = os.path.join(DATAPATHUNITS,chain.replace('|','-')+'_NA.pickle')

        with open(chain_filename, 'wb') as fh:
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(unit_info, fh, 2)

        chains.append(chain)
        messages.append('file_reading: Wrote %4d nucleotides to %s' % (len(chain_to_unit_data[chain]['ids']),chain_filename))

    return chains, messages


def writeNAPairwiseInteractions(structure,DATAPATHPAIRS,messages=[]):
    """
    Annotate pairwise interactions and write one file for the
    structure.
    """

    from fr3d.classifiers.NA_pairwise_interactions import annotate_nt_nt_in_structure

    # tell which categories of interactions to annotate
    categories = {}
    categories['basepair'] = []
    categories['basepair_detail'] = []
    categories['stacking'] = []
    categories['sO'] = []
    categories['sugar_ribose'] = []
    categories['coplanar'] = []

    interaction_to_list_of_tuples, category_to_interactions, timerData, pair_to_data = annotate_nt_nt_in_structure(structure,categories)

    # use RNA_pairs for now until _NA_pairs is established

    pdb_filename = os.path.join(DATAPATHPAIRS,structure.pdb+'_RNA_pairs.pickle')
    with open(pdb_filename, 'wb') as fh:
        # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
        pickle.dump(interaction_to_list_of_tuples, fh, 2)

    return messages


def writeNAUnitAnnotations(structure,DATAPATHUNITS,messages=[]):
    """
    Annotate chi angle and glycosidic bond conformation.
    """

    from fr3d.classifiers.NA_unit_annotation import annotate_bond_orientation

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
        chain_filename = os.path.join(DATAPATHUNITS,chain.replace('|','-')+'_NA_unit_annotations.pickle')

        with open(chain_filename, 'wb') as fh:
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(chain_to_unit_data[chain], fh, 2)

        messages.append('file_reading: Wrote %4d nucleotides to %s' % (len(chain_to_unit_data[chain]),chain_filename))

    return messages


def writeNABackboneData(structure,DATAPATHUNITS,messages=[]):
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
        chain_filename = os.path.join(DATAPATHUNITS,chain.replace('|','-')+'_NA_phosphate_sugar.pickle')

        unit_info = [chain_to_unit_data[chain]['ids']]
        unit_info.append(chain_to_unit_data[chain]['chainIndices'])
        unit_info.append(chain_to_unit_data[chain]['phosphate'])
        unit_info.append(chain_to_unit_data[chain]['sugar'])

        with open(chain_filename, 'wb') as fh:
            # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
            pickle.dump(unit_info, fh, 2)

        chains.append(chain)
        messages.append('file_reading: Wrote %4d nucleotides to %s' % (len(chain_to_unit_data[chain]['ids']),chain_filename))

    return messages


def writeProteinUnitData(structure,DATAPATHUNITS,messages=[]):
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
    protein_filename = os.path.join(DATAPATHUNITS,structure.pdb+'_protein.pickle')

    with open(protein_filename, 'wb') as fh:
        # Use 2 for "HIGHEST_PROTOCOL" for Python 2.3+ compatibility.
        pickle.dump([aa_unit_ids,aa_indices,aa_centers,aa_backbones], fh, 2)

    messages.append('file_reading: Wrote %4d amino acids to %s' % (len(aa_unit_ids),protein_filename))

    return messages


def processPDBFile(Q,structure_filename,file_id=None,pairs_only=False):
    """
    Read the coordinate file (in .cif or .pdb format) and save .pickle files.
    If we have to start with a .pdb or .cif file,
    load that file and write out all necessary .pickle files with required annotations.
    structure_filename can be:
        a full path to a file with an extension .cif or .pdb
        a file_id that is presumed to be located in CIFPATH and we will add .cif or .pdb
        a chain identifier, from which we extract the file_id and look in CIFPATH
    """

    from fr3d.classifiers.NA_pairwise_interactions import load_structure

    messages = []

    if Q.get("printFileOperations",False):
        print("  file_reading: processing %s" % structure_filename)

    if os.path.exists(structure_filename):
        # full path to a file that exists, use that
        cifPathAndFileName = structure_filename
    else:
        Q, CIFPATH = get_CIFPATH(Q)

        if '|' in structure_filename:
            # chain identifier given, extract file_id and look in CIFPATH
            fields = structure_filename.split('|')
            cifPathAndFileName = os.path.join(CIFPATH,fields[0]+'.cif')
        elif '.cif' in structure_filename.lower() or '.pdb' in structure_filename.lower():
            # name ends with .cif or .pdb
            cifPathAndFileName = os.path.join(CIFPATH,structure_filename)
        else:
            # no extension provided, hopefully we can hunt it down
            cifPathAndFileName = os.path.join(CIFPATH,structure_filename)

    # file_id is the analogue of PDB id, but also allows for other coordinate file names
    if not file_id:
        file_id = os.path.basename(cifPathAndFileName).replace('.gz','').replace('.cif','').replace('.pdb','')
        print('  file_reading: file_id is %s, cifPathAndFileName is %s' % (file_id,cifPathAndFileName))

    # read the coordinate file, specifying the file_id as the preferred name
    structure, messages = load_structure(cifPathAndFileName,file_id,preferred_id = file_id)

    if not structure and Q.get("printFileOperations",False):
        print('  file_reading: Could not read %s' % cifPathAndFileName)
        print('  file_reading: messages: %s' % messages)
        return [], file_id, messages

    # use file_id to determine the name, not the name in the file itself
    if not structure.pdb == file_id and Q.get("printFileOperations",False):
        print("  file_reading: structure name is %s but using file_id %s" % (structure.pdb,file_id))
        structure.pdb = file_id

    # make sure directories exist
    Q = checkDirectories(Q)

    # annotate pairwise interactions
    if Q.get("printFileOperations",False):
        print("  file_reading: annotating pairwise interactions for %s" % structure.pdb)

    messages = writeNAPairwiseInteractions(structure,Q['DATAPATHPAIRS'],messages)

    if Q.get("printFileOperations",False):
        print("  file_reading: saving in %s" % Q['DATAPATHPAIRS'])

    if not pairs_only:
        if Q.get("printFileOperations",False):
            print("  file_reading: calculating and writing unit data for %s" % structure.pdb)

        # write out centers and rotation matrices of nucleotides by chain
        chains, messages = writeNAUnitData(structure,Q['DATAPATHUNITS'],messages)

        # annotate chi angle and glycosidic bond conformation
        messages = writeNAUnitAnnotations(structure,Q['DATAPATHUNITS'],messages)

        # write out centers and indices for amino acids for the whole structure
        messages = writeProteinUnitData(structure,Q['DATAPATHUNITS'],messages)

        # write out centers and indices for amino acids for the whole structure
        messages = writeNABackboneData(structure,Q['DATAPATHUNITS'],messages)

    else:
        chains = None

    if Q.get("printFileOperations",False):
        for message in messages:
            print("  %s" % message)

    return chains, file_id, messages


def readPDBDatafile(DATAPATHUNITS):
    """
    Read .pickle file containing data about each nucleic-
    acid-containing PDB file.
    This file is produced by pipeline stage NA_datafile
    """

    datafile = {}

    filename = "NA_datafile.pickle"
    pathAndFileName = os.path.join(DATAPATHUNITS,filename)

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
    download from BGSU RNA site if necessary; create from .cif or .pdb if necessary.
    chainString is like '4V9F|1|0'
    """

    ids = []
    chainIndices = []
    centers = []
    rotations = []
    line_num = starting_index

    id_to_index = defaultdict()
    index_to_id = defaultdict()

    fields = chainString.split('|')
    file_id = fields[0]

    # filename = chainString.replace('|','-') + '_RNA.pickle'  # old, uses base center
    filename = chainString.replace('|','-') + '_NA.pickle'   # new on 2023-02-20, covers RNA and DNA, uses glycosidic atom

    Q, DATAPATHUNITS = get_DATAPATHUNITS(Q)

    pathAndFileName = os.path.join(DATAPATHUNITS,filename)

    if not os.path.exists(pathAndFileName) and not SERVER:
        # try to download .pickle file of RNA/DNA base center and rotation matrix
        if "PDB_data_file" in Q and file_id in Q["PDB_data_file"]:
            try:
                print("Attempting to download "+filename+" from BGSU RNA site")
                urlretrieve("http://rna.bgsu.edu/units/" + filename, pathAndFileName)
                print("Downloaded "+filename)
            except:
                print("Could not download %s from BGSU RNA site" % filename)
                pass

    if not os.path.exists(pathAndFileName) and not SERVER:
        # try to download .cif file and all related .pickle files
        if Q.get("printFileOperations",False):
            print("No positions file found for %s, attempting to create it from .cif or .pdb" % chainString)
        chains, file_id, messages = processPDBFile(Q,file_id,file_id)
        Q['userMessage'] += messages

    if os.path.exists(pathAndFileName):
        try:
            if sys.version_info[0] < 3:
                ids, chainIndices, centers, rotations = pickle.load(open(pathAndFileName,"rb"))
            else:
                ids, chainIndices, centers, rotations = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
        except:
            print("Could not read "+filename)
            Q["errorMessage"].append("Could not retrieve NA unit file "+filename)

        ok_ids = []
        ok_centers = []
        ok_rotations = []

        for i in range(0,len(ids)):
            if len(centers[i]) == 3 and rotations[i].shape == (3,3):
                ok_ids.append(ids[i])
                ok_centers.append(centers[i])
                ok_rotations.append(rotations[i])
                id_to_index[ids[i]] = line_num
                index_to_id[line_num] = ids[i]
                line_num += 1

        ids = ok_ids
        centers = np.asarray(ok_centers)
        rotations = np.asarray(ok_rotations)

    else:
        print("Could not find "+filename)
        Q["errorMessage"].append("Could not retrieve NA unit file "+filename)

        centers = np.asarray(np.zeros((0,3)))
        rotations = np.asarray(np.zeros((0,3,3)))

    for center in centers:
        if not len(center) == 3:
            print("  file_reading: center is %s" % center)

    # for rotation in rotations:
    #     print(rotation)

    # centers = np.asarray(centers)
    # rotations = np.asarray(rotations)

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

    Q, DATAPATHUNITS = get_DATAPATHUNITS(Q)
    filename = chainString + "_NA_phosphate_sugar.pickle"
    pathAndFileName = os.path.join(DATAPATHUNITS,filename)

    if not os.path.exists(pathAndFileName) and not SERVER:
        urlretrieve("http://rna.bgsu.edu/units/"+filename, pathAndFileName)
        print("Downloaded "+filename)

    if os.path.exists(pathAndFileName):
        if sys.version_info[0] < 3:
            try:
                ids, chainIndices, phosphate, sugar = pickle.load(open(pathAndFileName,"rb"))
            except:
                print("Could not read "+filename)
                Q["userMessage"].append("Could not retrieve NA backbone file "+filename)
        else:
            try:
                ids, chainIndices, phosphate, sugar = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
            except:
                print("Could not read "+filename)
                Q["userMessage"].append("Could not retrieve NA backbone file "+filename)

        for i in range(0,len(ids)):
            id_to_index[ids[i]] = line_num
            index_to_id[line_num] = ids[i]
            line_num += 1

    else:
        print("Could not find "+filename)
        Q["userMessage"].append("Could not retrieve NA backbone file "+filename)

    phosphate = np.asarray(phosphate)
    sugar     = np.asarray(sugar)

    return Q, phosphate, sugar, ids, id_to_index, index_to_id, chainIndices


def readProteinPositionsFile(Q, file_id, starting_index):

    ids = []
    chainIndices = []
    centers = []
    line_num = starting_index

    centers = np.zeros((0,3))  # of lines/nucleotides in file
    rotations = np.zeros((0,3,3))
    id_to_index = defaultdict()
    index_to_id = defaultdict()

    Q, DATAPATHUNITS = get_DATAPATHUNITS(Q)
    filename = file_id + "_protein.pickle"
    pathAndFileName = os.path.join(DATAPATHUNITS,filename)

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


def readNAPairsFileRaw(Q, file_id, alternate = ""):
    """
    Read the .pickle file that stores all pairwise interactions for a PDB id.
    Alternate could be "_exp" for experimental annotations, for example.
    """

    justDownloaded = False

    Q, DATAPATHPAIRS = get_DATAPATHPAIRS(Q)

    # new standard filename does not say RNA or DNA or NA
    pairsFileName = file_id + '_pairs.pickle'
    pathAndFileName = os.path.join(DATAPATHPAIRS+alternate,pairsFileName)

    if not os.path.exists(pathAndFileName):
        # old standard was _RNA_ but that is being phased out
        pairsFileName = file_id + '_RNA_pairs.pickle'
        pathAndFileName = os.path.join(DATAPATHPAIRS+alternate,pairsFileName)

        if not os.path.exists(pathAndFileName) and not SERVER:
            if Q.get("printFileOperations",False):
                print("  file_reading: Could not find "+pairsFileName+" in "+DATAPATHPAIRS+alternate)
            if not Q.get("computePairsLocally",False) or not file_id in Q.get("PDB_data_file",[]) or alternate:
                # try to download annotations of this structure from BGSU RNA site
                url = "http://rna.bgsu.edu/pairs" + alternate + "/" + pairsFileName
                urlretrieve(url, pathAndFileName) # testing
                if Q.get("printFileOperations",False):
                    print("  file_reading: downloaded "+pairsFileName)

                justDownloaded = True
            else:
                # compute pairwise interactions locally and store locally

                # continue with old standard until we switch over completely
                pairsFileName = file_id + '_RNA_pairs.pickle'
                pathAndFileName = os.path.join(DATAPATHPAIRS+alternate,pairsFileName)

                not_chains, file_id, messages = processPDBFile(Q,file_id,file_id,pairs_only=True)
                Q['userMessage'] += messages

    if Q.get("printFileOperations",False):
        print("  file_reading: reading "+pairsFileName+" from "+DATAPATHPAIRS+alternate)

    interactionToTriples = defaultdict(list)

    if os.path.exists(pathAndFileName):
        if sys.version_info[0] < 3:
            try:
                interactionToTriples = pickle.load(open(pathAndFileName,"rb"))
            except:
                print("Could not read "+pairsFileName+", it may not be available")
                Q["userMessage"].append("Could not read RNA pairs file "+pairsFileName)
                if justDownloaded:
                    os.remove(pathAndFileName)
        else:
            try:
                interactionToTriples = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
            except:
                if Q.get("printFileOperations",False):
                    print("Could not read "+pairsFileName+", it may not be available D*********************************")
                Q["userMessage"].append("Could not read RNA pairs file "+pairsFileName)
                if justDownloaded:
                    os.remove(pathAndFileName)

    return Q, interactionToTriples


def readNAPairsFile(Q, file_id, id_to_index, alternate = ""):
    """
    We need to read the pairwise interactions whether or not there are
    constraints, to be able to display the interactions in the output.
    alternate is for experimental or alternate lists of interacting pairs.
    The alternate .pickle files must be stored in a parallel directory.
    """

    Q, interactionToTriples = readNAPairsFileRaw(Q, file_id, alternate)

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
        pathAndFileName = os.path.join(DATAPATHUNITS,filename)

        if not os.path.exists(pathAndFileName) and not SERVER:
            urlretrieve("http://rna.bgsu.edu/units/" + filename, pathAndFileName)
            if Q.get("printFileOperations",False):
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
                    Q["userMessage"].append("Could not retrieve NA unit annotation file "+filename)
                    try:
                        os.remove(pathAndFileName)
                    except:
                        Q["userMessage"].append("Could not remove NA unit annotation file "+filename)
            else:
                try:
                    new_mapping = pickle.load(open(pathAndFileName,"rb"), encoding = 'latin1')
                    unit_id_to_annotation.update(new_mapping)
                except:
                    print("Could not read "+filename+", it may not be available D*********************************")
                    Q["userMessage"].append("Could not retrieve NA unit annotation file "+filename)
                    try:
                        os.remove(pathAndFileName)
                    except:
                        Q["userMessage"].append("Could not remove NA unit annotation file "+filename)


    return Q, unit_id_to_annotation
