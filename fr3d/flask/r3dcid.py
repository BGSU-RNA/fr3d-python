# Draw circular interaction diagrams for RNA and DNA structures

# Postscript references
# https://personal.math.ubc.ca/~cass/courses/ps.html
# http://paulbourke.net/dataformats/postscript/

import math
import os
from pathlib import Path
import pickle
import requests
import shutil
import subprocess
import sys
import urllib

arc_group_to_interactions = {}
arc_group_to_interactions["nested-wc"] = ['cWW']
arc_group_to_interactions["lr-wc"] = ['cWW']
arc_group_to_interactions["lr-non-wc"]     = ['cSS', 'cHH', 'cHS', 'cWH', 'cWS', 'tSS', 'tHH', 'tHS', 'tWH', 'tWS', 'tWW', 'cWW']
arc_group_to_interactions["nested-non-wc"] = ['cSS', 'cHH', 'cHS', 'cWH', 'cWS', 'tSS', 'tHH', 'tHS', 'tWH', 'tWS', 'tWW', 'cWW']
arc_group_to_interactions["bph"] = ['0BPh', '1BPh', '2BPh', '3BPh', '4BPh', '5BPh', '6BPh', '7BPh', '8BPh', '9BPh']
arc_group_to_interactions["br"] = ['0BR', '1BR', '2BR',  '3BR', '4BR', '5BR', '6BR', '7BR', '8BR', '9BR']
arc_group_to_interactions["sr"] = ['cSR','tSR']
arc_group_to_interactions["so"] = [a+b for a in ["s3","s5"] for b in ["O2'","O3'","O4'","O5'","OP1","OP2"]]
arc_group_to_interactions["stacking"] = ['s35','s33','s55','s53']
arc_group_to_interactions["near"] = ['all']
arc_group_names = ["nested-wc","lr-wc","lr-non-wc","nested-non-wc","bph","br","sr","so","stacking","near"]

arc_group_names_by_order = ["stacking","sr","so","br","bph","nested-non-wc","lr-non-wc","nested-wc","lr-wc","near"]

arc_group_name_to_text = {}
arc_group_name_to_text["nested-wc"]     = "%d nested Watson-Crick basepairs (AU, GC, GU cWW)"
arc_group_name_to_text["lr-wc"]         = "%d long-range Watson-Crick basepairs"
arc_group_name_to_text["nested-non-wc"] = "%d nested non-Watson-Crick basepairs"
arc_group_name_to_text["bonus"]         = "(Includes AA cWW, AG cWW, UU cWW, AG tHS, etc.)"
arc_group_name_to_text["lr-non-wc"]     = "%d long-range non-Watson-Crick basepairs"
arc_group_name_to_text["stacking"]      = "%d base stacking interactions"
arc_group_name_to_text["bph"]           = "%d base-phosphate interactions"
arc_group_name_to_text["br"]            = "%d base-ribose interactions"
arc_group_name_to_text["sr"]            = "%d sugar-ribose interactions"
arc_group_name_to_text["so"]            = "%d oxygen stacking interactions"
arc_group_name_to_text["near"]          = "%d near basepairs or other interactions"

# control points to map nucleotide number to the value of each of these parameters
control_base_number_font_size=  [(35,7),(82,7),(132,7),(294,4),(585,2),(1522,1),(3000,0.5),(4000,0.4),(5000,0.4),(7000,0.25),(11662,0.15),(18000,0.1)]
control_tick_mark_length=       [(35,4),(82,3),(500, 2),(1000,0.8),(2000, 0.6),(5000,0.5),(7228,0.4),(10000,0.3),(18263,0.2)] #length of radial tick mark lines

control_linewidth=              [(35,2),(82,2),(122,1.4),(294,0.5),(585,0.3),(1522,0.1),(3000,0.1),(5000,0.1),(7000,0.1),(18000,0.045)]

control_major_number_distance=  [(35,40),(82,45),(294,40),(585,30),(1000,20),(1522,18),(2318,15),(5000,10),(18263,10)]
control_major_number_font_size= [(35,12),(82,12),(100,11),(294,10),(585,10),(1522,10),(3031,9.5),(4748,9),(7228,9),(18263,8),(40000,3)]

control_chain_gap_size=         [(10,0),(35,1),(82,2),(294,5),(585,5),(1522,20),(18263,20)]
control_chain_gap_location=     [(10,0.35),(35,0.35),(75,0.35),(102,0.35),(170,0.35),(450,0.1),(700,0.05),(1500,0),(19000,0)]

control_circle_radius=          [(35,150),(82,170),(294,180),(585,200),(1522,220),(3031,230),(4217,240)]
control_helix_size=             [(75,7),(100,6),(1500,5),(6400,2.5),(8000,2.3),(16000,2)]

rfam_order = {
    "RF00001": 1,  # 5S rRNA goes first around the circle
    "RF00002": 2,  # 5.8S rRNA goes before the LSU
    "RF02543": 3,  # LSU families
    "RF02541": 3,
    "RF02540": 3,
    "RF02546": 3,
    "RF01960": 4,  # SSU families
    "RF00177": 4,
    "RF01959": 4,
    "RF02542": 4,
    "RF02545": 4,
    None: 5,
    "RF00005": 6    # tRNA goes last
}

simplify_entity_type = {}
simplify_entity_type['Polydeoxyribonucleotide (DNA)'] = 'dna'
simplify_entity_type['Polyribonucleotide (RNA)'] = 'rna'
simplify_entity_type['polyribonucleotide'] = 'rna'
simplify_entity_type['polydeoxyribonucleotide'] = 'dna'
simplify_entity_type['polydeoxyribonucleotide/polyribonucleotide hybrid'] = 'hybrid'
simplify_entity_type['DNA/RNA Hybrid'] = 'hybrid'


def get_fr3d_interaction_to_triple_list(pdb_id,data_directory=""):
    """
    This method downloads FR3D annotations if necessary,
    then reads the .pickle data file and returns a dictionary
    which has interaction types like cWW and 5BPh as keys
    and whose values are lists of *triples* of unit ids which
    make the given interaction and the number of nested cWW
    basepairs which the interaction crosses.
    """

    pdb_id = pdb_id.upper()

    fr3d_interaction_file = "%s_NA_pairs.pickle" % pdb_id
    fr3d_interaction_to_pair_list = None
    message = ""

    if data_directory == "":
        # use relative path and try to download
        data_directory = 'r3dcid_data'
    else:
        # use absolute path and count on the file being there
        pass

    if not os.path.exists(data_directory):
        Path(data_directory).mkdir(parents=True, exist_ok=True)

    fr3d_interaction_path_file = os.path.join(data_directory,"%s_NA_pairs.pickle" % pdb_id)

    if os.path.isfile(fr3d_interaction_path_file):
        try:
            with open(fr3d_interaction_path_file, 'rb') as opener1:
                fr3d_interaction_to_pair_list = pickle.load(opener1)
            message = None
            return fr3d_interaction_to_pair_list, message
        except:
            message = "Unable to open %s, removing it" % fr3d_interaction_path_file
            os.remove(fr3d_interaction_path_file)

    print("Downloading %s to %s" % (fr3d_interaction_file,fr3d_interaction_path_file))
    # https://rna.bgsu.edu/pairs/4TNA_NA_pairs.pickle
    url = "https://rna.bgsu.edu/pairs/%s" % fr3d_interaction_file
    try:
        if sys.version_info[0] < 3:
            urllib.urlretrieve(url, fr3d_interaction_path_file)  # python 2
        else:
            urllib.request.urlretrieve(url, fr3d_interaction_path_file)  # python 3
    except:
        return None, "Unable to download %s, exiting" % fr3d_interaction_file

    try:
        with open(fr3d_interaction_path_file, 'rb') as opener1:
            fr3d_interaction_to_pair_list = pickle.load(opener1)
        message = ""
        return fr3d_interaction_to_pair_list, message
    except:
        # if os.path.exists(fr3d_interaction_path_file):
        #     os.remove(fr3d_interaction_path_file)
        return None, "Unable to load %s, exiting" % fr3d_interaction_file


def process_input_chains(input_text,params={}):
    """
    Given just a PDB id like 8GLP, look up all the RNA and DNA chains in it
    Given a PDB id and one or more models like 8QO5|3,8QO5|2, focus on those models
    Given a PDB id, wildcard model, and chain like 1ABC|*|A, use all models, use that chain
    Given a PDB id and one or more assemblies like 4V9O_2,4V9O_4, show those
    Given one or more chains, restrict to those
    """

    pdb_id = ""
    requested_models = []
    requested_model_chain = []
    requested_assemblies = []
    requested_chains = []
    requested_symmetries = []

    input_text = input_text.replace(";",",")  # allow ; but prefer ,
    pieces = input_text.split("+")            # use + to separate model/chain combinations
    for piece in pieces:
        if "_" in piece:
            # old way of specifying an assembly, like 4V9O_2
            # new way is to use the assemblies field
            fields = piece.split("_")
            if not pdb_id:
                pdb_id = fields[0].split("|")[0]  # ignore model and chain
            assemblies = fields[1].split(",")     # allow multiple assemblies
            for assembly in assemblies:
                if not assembly in requested_assemblies:
                    requested_assemblies.append(assembly)
        else:
            fields = piece.split("|")
            # use first pdb_id encountered, in case of more than one
            if not pdb_id:
                pdb_id = fields[0]
            if len(fields) > 1:
                models = fields[1].split(",")
                # use a consistent wildcard character for models
                if len(models) == 1 and (models[0] == "" or models[0] == "m"):
                    models = ["*"]
                if "*" in models:
                    models = ["*"]

            if len(fields) == 3:
                chains = fields[2].split(",")
                for m in models:
                    if not m == '*':
                        # specific model requested
                        if not m in requested_models:
                            requested_models.append(m)
                        for chain in chains:
                            requested_model_chain.append((m,chain))
                    if (m == '*' or m == '1'):
                        for chain in chains:
                            if not chain in requested_chains:
                                # specific chain requested but no model restriction
                                requested_chains.append(chain)
            elif len(fields) == 2:
                for m in models:
                    if not m in requested_models:
                        requested_models.append(m)

    if not pdb_id:
        print('No PDB id found')
        return [], "", {}

    if "assembly" in params:
        pa = params['assembly'].replace(";",",")
        requested_assemblies = []
        for x in pa.split(","):
            if not x in requested_assemblies:
                # preserve order but de-duplicate
                requested_assemblies.append(x)

    if "symmetry" in params:
        pa = params['symmetry'].replace(";",",")
        requested_symmetries = pa.split(",")

    # build the filename, reflecting user requests
    filename = 'R3DCID_%s' % pdb_id.upper()

    if len(requested_assemblies) > 0:
        filename += '_assembly_%s' % ','.join(requested_assemblies)

    if len(requested_symmetries) > 0:
        filename += '_symmetry_%s' % ','.join(requested_symmetries)

    # convert integer symmetries to ASM_ to match database; database P_ will also match
    for i,s in enumerate(requested_symmetries):
        if "_" in s:
            continue
        try:
            # convert mere integers to proper symmetry names
            k = int(s)
            requested_symmetries[i] = "ASM_" + s
        except:
            pass

    requested_models_unique = sorted(set(requested_models))

    if len(requested_models_unique) > 1 or (len(requested_models_unique) == 1 and not requested_models_unique[0] == "1"):
        # more than one model specifically requested, or model other than 1 requested
        if len(requested_model_chain) > 0:
            filename += '_model-chain_' + '_'.join([("%s-%s" % (m,c)) for (m,c) in requested_model_chain])
        else:
            filename += '_model_' + '_'.join(requested_models)
    elif len(requested_chains) > 0:
        filename += '_' + '_'.join(requested_chains)

    if len(requested_models) == 0:
        # use a placeholder for now, later we'll loop over all available models
        requested_models = ['*']

    print('pdb_id                %s' % pdb_id)
    print('requested_models      %s' % requested_models)
    print('requested_model_chain %s' % requested_model_chain)
    print('requested_chains      %s' % requested_chains)
    print('requested_assemblies  %s' % requested_assemblies)
    print('requested_symmetries  %s' % requested_symmetries)

    return pdb_id, filename, requested_assemblies, requested_models, requested_model_chain, requested_chains, requested_symmetries


def set_parameters_from_input(params,filename,pdb_id):
    # set parameters and add fields to the filename

    # process display options to make sure they make sense, before basing the filename on them
    if 'coloring' in params:
        if params['coloring'].lower() in ['default','grayscale','wong']:
            coloring = params['coloring'].lower()
        else:
            print('Coloring scheme %s not recognized; using default' % params['coloring'])
            coloring = 'default'
    else:
        coloring = 'default'
    params['coloring'] = coloring

    # process user requests
    show = params.get('show','').lower()
    hide = params.get('hide','').lower()
    dim  = params.get('dim','').lower()

    if 'wc' in show.split(","):
        show = ",".join(show.split(",") + ["nested-wc","lr-wc"])

    if 'wc' in hide.split(","):
        hide = ",".join(hide.split(",") + ["nested-wc","lr-wc"])

    if 'wc' in dim.split(","):
        dim = ",".join(dim.split(",") + ["nested-wc","lr-wc"])

    if 'all' in hide:
        show = ''
        hide = ','.join(arc_group_names)
        dim = ''

    if 'all' in dim:
        show = ''
        hide = ''
        dim = ','.join(arc_group_names)

    if 'all' in show:
        show = ''
        hide = ''
        dim = ''

    # clean up the comma-separated lists
    params['show'] = clean_comma_list(show, arc_group_names)
    params['hide'] = clean_comma_list(hide, arc_group_names, show)
    params['dim']  = clean_comma_list(dim, arc_group_names, hide+","+show)

    print(params)

    text = params.get('text','').lower()
    params['text'] = clean_comma_list(text, ["basepair","stacking","bph","br","sr","so","near","all","helix","none"])

    if "none" in params['text']:
        params['text'] = 'none'
    if "all" in params['text']:
        params['text'] = "all"

    # add user-specified display options to the filename unless they are the default
    if not coloring == "default":
        filename += '_' + coloring

    if params['show']:
        filename += '_show_' + params['show']

    if params['hide']:
        filename += '_hide_' + params['hide']

    if params['dim']:
        filename += '_dim_' + params['dim']

    if params['text']:
        filename += '_text_' + params['text']

    if 'n3d' in params:
        if type(params['n3d']) == 'string' and params['n3d'].lower() == 'false':
            params['n3d'] = False
        if params['n3d'] == False:
            filename += "_n3d-false"

    if 'helix_size' in params:
        filename += '_hs_' + str(params['helix_size'])

    if not 'header' in params:
        params['header'] = 'all'

    params['header'] = clean_comma_list(params['header'],['title','method','source','release_date','resolution','all','none'])
    if len(params['header']) == 0:
        params['header'] = 'title,method,source,release_date,resolution'
    elif 'none' in params['header']:
        del params['header']
        filename += '_header_none'
    elif 'all' in params['header']:
        params['header'] = 'title,method,source,release_date,resolution'
    elif params['header'] == 'title,method,source,release_date,resolution':
        pass
    else:
        filename += '_header_' + params['header']

    if 'description' in params and len(params['description']) > 0:
        filename += '_description'

    return params, filename


def get_header_information(params,pdb_id):
    # download and display PDB title, method, resolution, release_date, source, as requested
    # https://rna.bgsu.edu/rna3dhub/rest/getPdbInfo?pdb=4V9F
    url = "https://rna.bgsu.edu/rna3dhub/rest/getPdbInfo?pdb=%s&format=json" % pdb_id
    response = requests.get(url)
    data = response.json()

    if not type(data) is dict:
        message = 'Could not get PDB title from %s' % url
    else:
        message = None
        if 'title' in params['header']:
            params['title'] = data['title']
        if 'method' in params['header']:
            params['method'] = data['method']
        if 'release_date' in params['header']:
            params['release_date'] = data['release_date']
        if 'source' in params['header']:
            if data['source']:
                params['source'] = data['source']
            else:
                params['source'] = 'Source not available'
        if 'resolution' in params['header']:
            if data['resolution']:
                params['resolution'] = data['resolution']
            else:
                params['resolution'] = 'No resolution'

    return params, message


def get_average_column_number(chain_id):
    # retrieve alignment of this chain to the Rfam family and average the column numbers

    url = "https://rna.bgsu.edu/correspondence/align_chains?chains=column,%s" % chain_id
    response = requests.get(url)
    if not response.status_code == 200:
        print("align_chains request failed for chain %12s with status code: %s" % (chain_id,response.status_code))
        return 0
    else:
        column_unit_id = response.text.strip().split("\n")
        column_total = 0
        column_count = 0
        for line in column_unit_id:
            fields = line.split()
            try:
                column = int(fields[0])
                column_total += column
                column_count += 1
            except:
                pass

        return column_total / column_count


def order_chains_around_diagram(pdb_id, requested_assemblies, requested_models, requested_model_chain, requested_chains, requested_symmetries, interaction_to_triple_list):
    """
    Download information about the chains in pdb_id
    Determine the order of assemblies, chains, symmetries around the circle
    """

    # get list of assemblies and their chains and symmetries in this pdb_id
    # https://rna.bgsu.edu/rna3dhub/rest/getChainInfo?pdb=4V9F
    # https://rna.bgsu.edu/rna3dhub/rest/getChainInfo?pdb=8QYX
    # https://rna.bgsu.edu/rna3dhub/rest/getChainInfo?pdb=6SAE
    # https://rna.bgsu.edu/rna3dhub/rest/getChainInfo?pdb=1Z1C
    # https://rna.bgsu.edu/rna3dhub/rest/getChainInfo?pdb=4FTB
    url = "https://rna.bgsu.edu/rna3dhub/rest/getChainInfo?pdb=%s" % pdb_id
    response = requests.get(url)
    data = response.json()
    numbered_symmetries = False   # ASM_2 or P_6 and the like

    message = None

    if not type(data) is dict:
        message = 'Did not get required data from %s, exiting' % url
        return [], {}, message

    if len(data['chains']) == 0:
        message = 'Did not get any chain data from %s, exiting' % url
        return [], {}, message

    # detect symmetry lists like (1-60) and (1,2,6,10,23,24)
    expand_symmetry_lists = False
    for data_line in data['assemblies']:
        if data_line['symmetry'] and "(" in data_line['symmetry']:
            expand_symmetry_lists = True

    # ASM_ is used in all new structures, but P_ was used in older structures like 1A34
    symmetry_prefix = "ASM_"
    if requested_symmetries:
        if "P_" in requested_symmetries[0]:
            symmetry_prefix = "P_"

    # fill in null symmetry ids with something reasonable, expand (1-40) symmetries
    assembly_symmetry_count = []
    new_data_lines = []
    unique_data_strings = set()
    for data_line in data['assemblies']:
        # for whatever reason, some values of data_line are duplicated
        data_string = str(data_line)
        if data_string in unique_data_strings:
            continue
        unique_data_strings.add(data_string)

        if not 'symmetry_id' in data_line:
            data_line['symmetry_id'] = "999"

        if expand_symmetry_lists:
            numbered_symmetries = True
            if data_line['symmetry'] and "X" in data_line['symmetry']:
                # not ready to process such symmetries
                continue
            elif data_line['symmetry'] and "(" in data_line['symmetry']:
                # example (1-40) from 6SAE or (1,2,6,10,23,24) from 1Z1C
                # change this data_line and then add more copies of it

                # produce a list of all symmetry operators required
                sym_text = data_line['symmetry'].replace("(","").replace(")","")
                sym_list = []
                for sym_item in sym_text.split(","):
                    sym_limits = sym_item.split("-")
                    if len(sym_limits) == 1:
                        sym_list.append(int(sym_item))
                    elif len(sym_limits) == 2:
                        sym_list += list(range(int(sym_limits[0]),int(sym_limits[1])+1))

                assembly_symmetry_count.append((data_line['assembly_id'],len(sym_list)))

                for i in sym_list:
                    dl = {}
                    dl['assembly_id'] = data_line['assembly_id']
                    dl['chain_name']  = data_line['chain_name']
                    dl['symmetry'] = '%s%d' % (symmetry_prefix,i)
                    dl['symmetry_id'] = '%d' % i
                    new_data_lines.append(dl)
            else:
                new_data_lines.append(data_line)

        elif data_line['symmetry'] is None:
            numbered_symmetries = True
            data_line['symmetry'] = "ASM_%s" % data_line['symmetry_id']
            new_data_lines.append(data_line)
        else:
            new_data_lines.append(data_line)

    # replace previous list of data, with symmetry lists expanded
    data['assemblies'] = new_data_lines

    # sort assembly data by assembly_id
    data_lines = sorted(data['assemblies'], key = lambda x : (x['assembly_id']))

    all_assemblies = set([x['assembly_id'] for x in data_lines])

    if len(requested_assemblies) > 0:
        # focus on just the assemblies requested, in the order requested
        ok_assemblies = []
        for ra in requested_assemblies:
            if ra in all_assemblies:
                ok_assemblies.append(ra)
    else:
        # do not miss any assemblies
        ok_assemblies = sorted(all_assemblies)

    # special case, two ribosomes in one assembly
    if pdb_id.upper() == '5J7L' and not requested_chains and not requested_model_chain:
        requested_chains = 'DB,DA,AA,CB,CA,BA'.split(",")
        requested_model_chain = [('1',x) for x in requested_chains]
        requested_models = ['1']

    # accumulate data about each chunk to be shown around the circle:
    # model, assembly, rfam priority, chain name, symmetry operator priority,
    sortable_chain_data = []
    solitary_unit_ids = []
    rfam_found = False
    chain_symmetry_to_assembly = {}
    same_chain_symmetry_different_assemblies = set()
    rfam_max_length = 0

    # assemble a list of all model, assembly, chain, symmetries with nucleic acids
    for model in requested_models:
        for assembly in ok_assemblies:
            for assembly_line in data['assemblies']:
                # every chain and symmetry in each assembly is listed on a separate line
                if assembly_line['assembly_id'] == assembly:
                    chain_name = assembly_line['chain_name']

                    if (model,chain_name) in requested_model_chain or (len(requested_model_chain)== 0 and chain_name in requested_chains) \
                        or (len(requested_model_chain) == 0 and len(requested_chains) == 0) :

                        symmetry = assembly_line['symmetry']
                        symmetry_id = assembly_line['symmetry_id']

                        chain_symmetry = (chain_name,symmetry,symmetry_id)
                        if not chain_symmetry in chain_symmetry_to_assembly:
                            chain_symmetry_to_assembly[chain_symmetry] = set(assembly)
                        elif not assembly in chain_symmetry_to_assembly[chain_symmetry]:
                            for already in chain_symmetry_to_assembly[chain_symmetry]:
                                same_chain_symmetry_different_assemblies.add((assembly,already))
                                same_chain_symmetry_different_assemblies.add((already,assembly))
                            chain_symmetry_to_assembly[chain_symmetry].add(assembly)

                        # print("model::assembly::chain::symmetry_id::symmetry",model,assembly,chain_name,symmetry_id,symmetry)
                        # print(chain_symmetry)
                        # print(chain_symmetry_to_assembly[chain_symmetry])
                        # print(same_chain_symmetry_different_assemblies)

                        if requested_symmetries and not symmetry in requested_symmetries:
                            continue

                        chain_info = data['chains'].get(chain_name,{})

                        na_type = simplify_entity_type.get(chain_info.get("entity_macromolecule_type", ""),"")
                        if na_type in ["rna", "dna", "hybrid"] or "solitary" in chain_info:
                            display_name = chain_info.get("standardized_name", None)
                            if display_name == None:
                                display_name = chain_info.get("pdbx_description", None)

                            if display_name != None:
                                display_name = display_name.split(";")[0] # use long version of standardized name
                            else:
                                display_name = "Not listed"

                            if "rfam_family" in chain_info:
                                rfam_family = chain_info["rfam_family"]
                                rfam_found = True
                                try:
                                    chain_length = int(chain_info["chain_length"])
                                except:
                                    chain_length = 0
                                rfam_max_length = max(rfam_max_length,chain_length)
                            else:
                                rfam_family = None

                            new_data = {}
                            new_data['pdb_id'] = pdb_id
                            new_data['model'] = model
                            new_data['assembly'] = assembly
                            if assembly in requested_assemblies:
                                new_data['assembly_priority'] = requested_assemblies.index(assembly)
                            else:
                                new_data['assembly_priority'] = assembly
                            new_data['chain_name'] = chain_name
                            if len(requested_model_chain) > 0:
                                new_data['model_chain_priority'] = requested_model_chain.index((model,chain_name))
                            else:
                                new_data['model_chain_priority'] = 1
                            if len(requested_models) > 0:
                                new_data['chain_priority'] = requested_models.index(model)
                            else:
                                new_data['chain_priority'] = 1
                            new_data['display_name'] = display_name
                            new_data['pdbx_description'] = chain_info.get("pdbx_description", "")
                            new_data['rfam_family'] = rfam_family
                            new_data['rfam_priority'] = rfam_order.get(rfam_family,7)
                            new_data['symmetry'] = symmetry
                            new_data['symmetry_id'] = symmetry_id

                            if requested_symmetries:
                                new_data['symmetry_priority'] = requested_symmetries.index(symmetry)
                            else:
                                if symmetry == '1_555' and symmetry_id == "1":
                                    # in 8PZP, symmetry 6 is 1_555 and symmetry 1 is not
                                    new_data['symmetry_priority'] = 0
                                elif symmetry == '1_555':
                                    # cases like symmetry operator 6 in 8PZP
                                    try:
                                        new_data['symmetry_priority'] = int(symmetry_id)
                                    except:
                                        # just in case
                                        new_data['symmetry_priority'] = 998
                                elif symmetry and symmetry[0].isalpha() and len(symmetry.split("_")) >= 2:
                                    # cases like P_1, H_3, ASM_4...
                                    try:
                                        new_data['symmetry_priority'] = int(symmetry.split("_")[1])
                                    except:
                                        # things like P_P
                                        new_data['symmetry_priority'] = 998
                                elif symmetry and len(symmetry.split("_")) >= 2:
                                    # cases like 12_665
                                    try:
                                        new_data['symmetry_priority'] = int(symmetry.split("_")[0])
                                    except:
                                        new_data['symmetry_priority'] = 998

                                else:
                                    new_data['symmetry_priority'] = 999

                            new_data['chain_id'] = "%s|%s|%s" % (pdb_id, model, chain_name)

                            if "solitary" in chain_info:
                                new_data['solitary'] = chain_info['solitary']
                                solitary_unit_ids += chain_info['solitary']

                            sortable_chain_data.append(new_data)

    if len(sortable_chain_data) == 0:
        message = 'Could not find chains in %s, exiting' % pdb_id
        return [], {}, message

    if len(requested_assemblies) == 0:
        # focus on assemblies that are actually needed
        ok_assemblies = sorted(set([x['assembly'] for x in sortable_chain_data]))

    # some chains are clearly labeled but do not map to Rfam families
    # we can place some of them in the right order with ribosomal chains
    text_to_rfam_priority = {}
    text_to_rfam_priority['5s ribosomal'] = 1    # sometimes not caught by rfam mapping
    text_to_rfam_priority['5s rrna'] = 1
    text_to_rfam_priority['5s rna'] = 1
    text_to_rfam_priority['3s ribosomal'] = 1.5
    text_to_rfam_priority['3s rrna'] = 1.5
    text_to_rfam_priority['3s rna'] = 1.5
    text_to_rfam_priority['2s ribosomal'] = 2.5  # rare, goes at the end of 5.8S in drosophila
    text_to_rfam_priority['2s rrna'] = 2.5
    text_to_rfam_priority['2s rna'] = 2.5

    # try to identify ribosomes like 8APN, 6AZ3 where the LSU is broken into many chains
    # need more than two chains to match LSU or SSU families
    # if that doesn't work, will need to parse descriptions to figure it out
    # 5J7L has two ribosomes in one assembly
    # 4V3P has 23 ribosomes in 23 symmetry operators, all in one assembly
    assembly_rfam_priority_to_count = {}
    assembly_to_count = {}
    rfam_priority_to_count = {}
    for data in sortable_chain_data:
        assembly = data['assembly']
        if not assembly in assembly_to_count:
            assembly_to_count[assembly] = 0
        assembly_to_count[assembly] += 1

        # print('Fixing %s' % data)

        if rfam_found and not data.get('rfam_family',None):
            # rfam chains present but this chain is not recognized
            pdbx_description = data.get("pdbx_description","").lower()
            # print(pdbx_description)
            for text, priority in text_to_rfam_priority.items():
                # print('  ',text)
                if text in pdbx_description:
                    data['rfam_priority'] = priority
                    # print('Found %s in pdbx_description' % text)

        rfam_priority = data['rfam_priority']
        if not rfam_priority in rfam_priority_to_count:
            rfam_priority_to_count[rfam_priority] = 0
        rfam_priority_to_count[rfam_priority] += 1

        if rfam_priority in [3,4]:
            # LSU or SSU
            key = assembly + "_" + str(data['rfam_priority'])
            if not key in assembly_rfam_priority_to_count:
                assembly_rfam_priority_to_count[key] = 0
            assembly_rfam_priority_to_count[key] += 1

    verbose = 0
    if False:
        verbose = 1

    # order the chains around the circle
    if rfam_priority_to_count.get(3,0) > 10 and rfam_priority_to_count.get(4,0) > 10:
        # multiple Rfam-identified LSU or SSU chains in one assembly, like 4V3P
        if verbose > 0:
            print('More than 10 Rfam identified LSU or SSU chains, using Rfam order')
        sorted_chain_data = sorted(sortable_chain_data, key=lambda x: (x['assembly_priority'],x['model_chain_priority'],x['chain_priority'],x['symmetry_priority'],x['symmetry'],x['rfam_priority'],x['chain_name']))
    elif assembly_rfam_priority_to_count and max(assembly_rfam_priority_to_count.values()) > 1 and max(assembly_to_count.values()) > 7:
        # more than 7 chains in one assembly, and Rfam matches found, as in 9i05 Toxoplasma gondii mitochondrial ribosome
        # change rfam_priority 4 (SSU), 5 (None), 7 (other Rfam) to 3 to try to get LSU, SSU in the right place by sorting alphabetically
        # but keep 5S and 5.8S at the beginning and tRNA at the end
        for x in sortable_chain_data:
            if x['rfam_priority'] in [4,5,7]:
                x['rfam_priority'] = 3
        # sort first by assembly, then alphabetically by chain name, then by symmetry
        sorted_chain_data = sorted(sortable_chain_data, key=lambda x: (x['assembly_priority'],x['model_chain_priority'],x['chain_priority'],x['symmetry_priority'],x['symmetry'],x['rfam_priority'],x['chain_name']))
        if verbose > 0:
            print('Many chains in one assembly, plus Rfam matches')
    elif assembly_rfam_priority_to_count and not requested_chains and not requested_model_chain and max(assembly_rfam_priority_to_count.values()) > 1 and rfam_max_length > 200:
        # get some tricky ribosomes in the right order
        # more than one match to the same Rfam family in the same assembly, for example 9TVU chains 4 and 1,
        # 4V8Y chains B5 and CN
        # when more than one chain maps to the same Rfam family, tweak Rfam priority by column number in alignment
        for assembly_rfam_priority, count in assembly_rfam_priority_to_count.items():
            if count > 1:
                rfam_priority = int(assembly_rfam_priority.split("_")[1])
                for k, scd in enumerate(sortable_chain_data):
                    if scd['rfam_priority'] == rfam_priority:
                        chain_id = "%s|1|%s" % (scd['pdb_id'].upper(),scd['chain_name'])
                        scd['rfam_priority'] += get_average_column_number(chain_id) / 1000000

        sorted_chain_data = sorted(sortable_chain_data, key=lambda x: (x['assembly_priority'],x['model_chain_priority'],x['chain_priority'],x['symmetry_priority'],x['symmetry'],x['rfam_priority'],x['chain_name']))
        if verbose > 0:
            print('Multiple Rfam matches to the same family, checking columns and then using Rfam priority over chain name')
            print(assembly_rfam_priority_to_count)
    elif assembly_rfam_priority_to_count and max(assembly_rfam_priority_to_count.values()) >= 1 and rfam_max_length > 200:
        # has an SSU or LSU chain, where we have a system for ordering
        # sort first by assembly, then by symmetry, then by rfam family priority order, then alphabetically by chain name
        sorted_chain_data = sorted(sortable_chain_data, key=lambda x: (x['assembly_priority'],x['model_chain_priority'],x['chain_priority'],x['symmetry_priority'],x['symmetry'],x['rfam_priority'],x['chain_name']))
        if verbose > 0:
            print('Rfam matches, sorting by Rfam priority over chain name')
            print(assembly_rfam_priority_to_count)
    elif numbered_symmetries:
        # sort by symmetry first, to keep chains together
        sorted_chain_data = sorted(sortable_chain_data, key=lambda x: (x['assembly_priority'],x['symmetry_priority'],x['symmetry'],x['model_chain_priority'],x['chain_priority'],x['chain_name']))
    else:
        # sort first by assembly, then requested chains, then by symmetry, then alphabetically by chain name
        sorted_chain_data = sorted(sortable_chain_data, key=lambda x: (x['assembly_priority'],x['model_chain_priority'],x['chain_priority'],x['symmetry_priority'],x['symmetry'],x['chain_name']))

        # for each assembly, group chains and symmetries by cWW interactions between them
        for assembly in ok_assemblies:
            # record what cWW group each chain is in, starting with all in different groups
            chain_symmetry_to_group = {}
            for i, chain_data in enumerate(sorted_chain_data):
                if chain_data['assembly'] == assembly:
                    chain_symmetry = (chain_data['chain_name'],chain_data['symmetry'])
                    chain_symmetry_to_group[chain_symmetry] = i

            solitary_unit_ids = set(solitary_unit_ids)
            solitary_chain_interactions = {}

            # identify chains which have cWW pairs between them, and form cWW groups
            for interaction in interaction_to_triple_list:
                inter_simple = interaction.lower().replace("a","")
                if inter_simple == 'cww' or len(solitary_unit_ids) > 0:
                    for u1,u2,c in interaction_to_triple_list[interaction]:
                        fields1 = u1.split("|")
                        fields2 = u2.split("|")
                        chain1 = fields1[2]
                        chain2 = fields2[2]

                        if not chain1 == chain2:
                            if len(fields1) == 9:
                                chain_symmetry1 = (chain1,fields1[8])
                            else:
                                chain_symmetry1 = (chain1,'1_555')
                            if len(fields2) == 9:
                                chain_symmetry2 = (chain2,fields2[8])
                            else:
                                chain_symmetry2 = (chain2,'1_555')

                            # focus on chain symmetry combinations in this assembly
                            # structures with many symmetries may not have the unit id matching the chain
                            if not chain_symmetry1 in chain_symmetry_to_group:
                                continue
                            if not chain_symmetry2 in chain_symmetry_to_group:
                                continue

                            g1 = chain_symmetry_to_group[chain_symmetry1]
                            g2 = chain_symmetry_to_group[chain_symmetry2]
                            if (not g1 == g2 and inter_simple == 'cww') or (u1 in solitary_unit_ids) or (u2 in solitary_unit_ids):
                                # collapse groups g1 and g2 into group min(g1,g2)
                                for cs in chain_symmetry_to_group.keys():
                                    if chain_symmetry_to_group[cs] in [g1,g2]:
                                        chain_symmetry_to_group[cs] = min(g1,g2)

                            if u1 in solitary_unit_ids:
                                if not chain_symmetry1 in solitary_chain_interactions:
                                    solitary_chain_interactions[chain_symmetry1] = {}
                                if not chain_symmetry2 in solitary_chain_interactions[chain_symmetry1]:
                                    solitary_chain_interactions[chain_symmetry1][chain_symmetry2] = set()
                                solitary_chain_interactions[chain_symmetry1][chain_symmetry2].add(inter_simple)
                            if u2 in solitary_unit_ids:
                                if not chain_symmetry2 in solitary_chain_interactions:
                                    solitary_chain_interactions[chain_symmetry2] = {}
                                if not chain_symmetry1 in solitary_chain_interactions[chain_symmetry2]:
                                    solitary_chain_interactions[chain_symmetry2][chain_symmetry1] = set()
                                solitary_chain_interactions[chain_symmetry2][chain_symmetry1].add(inter_simple)

            # update group membership for each chain-symmetry in this assembly
            chain_symmetry_to_final_group = {}
            for i, chain_data in enumerate(sorted_chain_data):
                if chain_data['assembly'] == assembly:
                    chain_symmetry = (chain_data['chain_name'],chain_data['symmetry'])
                    if solitary_unit_ids:
                        # split the chains in each group so solitary chains can go in between
                        chain_data['final_cww_group'] = chain_symmetry_to_group[chain_symmetry] + i/1000
                    else:
                        chain_data['final_cww_group'] = chain_symmetry_to_group[chain_symmetry]
                    chain_symmetry_to_final_group[chain_symmetry] = chain_data['final_cww_group']

            # place chains with solitary nucleotides as well as possible within their cWW group
            for chain_data in sorted_chain_data:
                if chain_data['assembly'] == assembly:
                    chain_symmetry = (chain_data['chain_name'],chain_data['symmetry'])
                    if 'solitary' in chain_data:
                        if chain_symmetry in solitary_chain_interactions:
                            s = 0
                            c = 0
                            for other_chain, interactions in solitary_chain_interactions[chain_symmetry].items():
                                if verbose > 0:
                                    print('Where to put solitary chain',chain_symmetry, other_chain, interactions)
                                if other_chain in chain_symmetry_to_final_group:
                                    if 'cww' in interactions:
                                        c += 2
                                        s += 2*chain_symmetry_to_final_group[other_chain]
                                    else:
                                        c += 1
                                        s += chain_symmetry_to_final_group[other_chain]
                            if c > 0:
                                chain_data['final_cww_group'] = s/c  # weighted average, to place between chains

        # group chains having cWW pairs between them, so cWW groups don't cross each other
        # sort first by assembly, then by cWW group, then alphabetically by chain name, then by symmetry
        # within a cWW group, there may be crossing arcs, as with 5B2P
        sorted_chain_data = sorted(sorted_chain_data, key=lambda x: (x['assembly_priority'],x['model_chain_priority'],x['chain_priority'],x['final_cww_group'],x['chain_name'],x['symmetry_priority'],x['symmetry']))

    assemblies = {}
    assemblies['ok_assemblies'] = ok_assemblies
    assemblies['message'] = ''

    if expand_symmetry_lists:
        # assemblies are made up of different groupings of symmetries
        # do not show interactions across assemblies
        assemblies['valid_assembly_pairs'] = [(a,a) for a in ok_assemblies]
    else:
        # assemblies may simply be different groupings of chains as in 4V9O
        assemblies['valid_assembly_pairs'] = []
        for a in ok_assemblies:
            for b in ok_assemblies:
                if not (a,b) in same_chain_symmetry_different_assemblies:
                    assemblies['valid_assembly_pairs'].append((a,b))

    return sorted_chain_data, assemblies, message


def interaction_priority(interaction,interactions):
    if interaction[0] == "n":
        return 100
    else:
        try:
            i = interactions.index(interaction)
            return i
        except:
            return 998


def linear_interpolation(points, x, tail='constant'):
    """
    Use control (x,y) points to map input number x to corresponding y
    by interpolating between the two nearest control points
    Beyond the largest control point, the tail can be 'constant' or 'decay'
    like k/x
    """

    # sort control points by x value
    points.sort(key=lambda point: point[0])

    # find interval that given x value lies in
    i = 1
    while i < len(points)-1 and x > points[i][0]:
        i += 1

    x1, y1 = points[i-1]
    x2, y2 = points[i]

    if x > x2 and i == len(points)-1 and tail=='decay':
        # beyond the largest control point and the tail should decay
        y = y2 * x2 / x
    elif x > x2:
        # use the y value of the rightmost control point
        y = y2
    elif x < x1:
        # use the y value of the leftmost control point
        y = y1
    else:
        # interpolate between the control points
        y = y1 + 1.0 * (x - x1) * (y2 - y1) / (x2 - x1)

    return y


def short_id(unit_id, n1_chain='', chain_to_longest_base_length={}):
    """
    Produce a short version of the unit id for use around the circle
    """
    fields = unit_id.split("|")
    t = ""
    if n1_chain:
        if not fields[2] == n1_chain:
            t = fields[2] + "|"   # start base with chain identifier
    if len(fields) >= 4:
        b = fields[3]      # base sequence
        c = fields[2]      # chain
        while len(b) < chain_to_longest_base_length.get(c,1):
            b = " " + b
        t += b
    if len(fields) >= 5:
        t += fields[4]      # number
    if len(fields) >= 8:
        t += fields[7]      # insertion code
    return t


def first_five_fields(unit_id):
    """
    Extract the first five fields of the unit id and make the model 1
    """

    fields = unit_id.split("|")
    if len(fields) >= 5:
        fields[1] = "1"
        return "|".join(fields[0:5])
    else:
        return unit_id


def divide_and_round(number_of_nucleotides):
    """
    Which big numbers to show around the circle
    """
    if number_of_nucleotides >= 15000:
        result = 200
    elif number_of_nucleotides >= 2000:
        result = 100
    elif number_of_nucleotides >= 1500:
        result = 50
    elif number_of_nucleotides >= 400:
        result = 20
    elif number_of_nucleotides >= 300:
        result = 10
    else:
        result = 10
    return result


def clean_comma_list(s, valid_list, exclude_string=''):
    """
    Given a comma-separated string s and a list of valid items,
    return a comma-separated string of items that are in the valid list,
    in the order that they appear in the valid list
    """
    s_list = s.strip().split(",")
    ok_list = []
    for h in valid_list:
        if h in exclude_string:
            continue
        elif h in s_list:
            ok_list.append(h)

    return ",".join(ok_list)


def draw_arcs(pairs_and_crossing, assembly_unit_id_to_angle, arcs_drawn, crossing, base_combination, unit_id_to_standard, angle_shifter, circle_radius, unit_id_to_annotation, helix_to_number_location, helix_radius_shift, unit_id_to_assemblies, valid_assembly_pairs):
    """
    Create the PostScript and SVG commands for a set of arcs all of the same color
    """

    # Watson-Crick base combinations in terms of RNA bases
    WC_bc = set(["A,U", "U,A", "G,C", "C,G", "G,U", "U,G"])

    dna_to_rna = {'DA':'A','DC':'C','DG':'G','DT':'U','A':'A','C':'C','G':'G','U':'U'}

    # constant which divides by two and changes degrees to radians
    Scale = math.pi/360.0

    PSarcs = []
    SVGarcs = []
    num_arcs_drawn = 0

    for n1,n2,crossing_number in pairs_and_crossing:
        # don't draw self interactions
        if n1 == n2:
            continue

        # don't draw the same arc twice for the same interaction
        if (n1,n2) in arcs_drawn:
            continue

        # record that we have drawn an arc between n1 and n2; don't draw it again
        arcs_drawn.add((n1,n2))
        arcs_drawn.add((n2,n1))

        # loop over assemblies where these unit ids appear
        for a1 in unit_id_to_assemblies.get(n1,[]):
            for a2 in unit_id_to_assemblies.get(n2,[]):
                # sometimes we draw arcs between assemblies, sometimes not
                if not (a1,a2) in valid_assembly_pairs:
                    continue
                an1 = (a1,n1)
                an2 = (a2,n2)
                if not an1 in assembly_unit_id_to_angle:
                    continue
                if not an2 in assembly_unit_id_to_angle:
                    continue

                # get the correct crossing number
                cn = int(crossing_number)
                if (crossing == "all") or (crossing == "lr" and cn > 0) or (crossing == "nested" and cn == 0):
                    # for Watson-Crick pairs get AU, GC, GU base combination only
                    this_bc = dna_to_rna[unit_id_to_standard[n1]]+","+dna_to_rna[unit_id_to_standard[n2]]
                    if (base_combination == "all") or (base_combination == "wc" and this_bc in WC_bc) \
                        or (base_combination == "non_wc" and not this_bc in WC_bc):

                        num_arcs_drawn += 1

                        # map unit id to angle
                        # for each interaction we will add a special number to angle1 and angle2
                        # such that they move slightly and solve the problem of overlapping
                        angle1 = assembly_unit_id_to_angle[an1]+angle_shifter
                        angle2 = assembly_unit_id_to_angle[an2]+angle_shifter

                        # put angles in increasing order, because arcs are drawn counterclockwise
                        angle1, angle2 = sorted([angle1,angle2])

                        # use the shorter distance around the circle, less than 180 degrees
                        if angle2 - angle1 > 180:
                            angle1,angle2 = angle2,(angle1+360)  # change order, shift one by 360 degrees

                        # ArcRadius is the altitude of a right triangle
                        ArcRadius = circle_radius * math.tan((angle2-angle1)*Scale)

                        if abs(ArcRadius) > 9999:
                            # very large ArcRadius, draw a straight line across the circle
                            point1X = circle_radius*math.cos(2*angle1*Scale)
                            point1Y = circle_radius*math.sin(2*angle1*Scale)
                            point2X = circle_radius*math.cos(2*angle2*Scale)
                            point2Y = circle_radius*math.sin(2*angle2*Scale)

                            PSarcs.append('%f %f newpath moveto' % (point1X, point1Y))
                            PSarcs.append("%f %f lineto" % (point2X, point2Y))
                            PSarcs.append("stroke")

                            SVGarcs.append('<line x1="%f" y1="%f" x2="%f" y2="%f" />' % (point1X, -point1Y, point2X, -point2Y))
                        else:
                            # starting and ending angles of the arc; perpendicular to nt angles
                            ArcStartAngle = angle2 + 90
                            ArcEndAngle = 270 + angle1
                            # CenterDistance is the hypotenuse of the same right angle
                            CenterDistance = math.sqrt(circle_radius*circle_radius + ArcRadius*ArcRadius)
                            # the center of the arc is along the bisector between the nucleotides
                            CenterX = CenterDistance * math.cos((angle1+angle2)*Scale)
                            CenterY = CenterDistance * math.sin((angle1+angle2)*Scale)

                            PSarcs.append("%f %f %f %f %f arc" % (CenterX,CenterY,ArcRadius,ArcStartAngle,ArcEndAngle))
                            PSarcs.append("stroke")

                            # Convert degrees to radians
                            sa_rad = math.radians(angle1)
                            ea_rad = math.radians(angle2)

                            # radius of the circular arc
                            r = ArcRadius

                            # Compute start and end points on the outer black circle
                            x1 = circle_radius * math.cos(sa_rad)
                            y1 = circle_radius * math.sin(sa_rad)
                            x2 = circle_radius * math.cos(ea_rad)
                            y2 = circle_radius * math.sin(ea_rad)

                            # Determine if arc > 180° for large-arc-flag
                            large_arc_flag = 0 if (ArcStartAngle - ArcEndAngle) % 360 > 180 else 1
                            sweep_flag = 1

                            SVGarcs.append('<path d="M %f %f A %f %f 0 %d %d %f %f"/>' % (x1, -y1, r, r, large_arc_flag, sweep_flag, x2, -y2))

                            # accumulate helix numbers of cWW pairs in the same helix
                            if base_combination == "wc":
                                # use the first five fields, works for all symmetries
                                n1_5 = first_five_fields(n1)
                                n2_5 = first_five_fields(n2)

                                helix_number1 = unit_id_to_annotation.get(n1_5,'not1')
                                helix_number2 = unit_id_to_annotation.get(n2_5,'not2')

                                # matching helix annotation
                                if helix_number1 == helix_number2 and 'helix' in helix_number1.lower():

                                    # middle_angle = (angle1 + angle2) / 2.0
                                    middle_angle =(ArcEndAngle + ArcStartAngle)/2.0

                                    # shift slightly toward the center of the circle
                                    midpoint_x = CenterX + (ArcRadius+helix_radius_shift) * math.cos(2*(middle_angle) * Scale)
                                    midpoint_y = CenterY + (ArcRadius+helix_radius_shift) * math.sin(2*(middle_angle) * Scale)

                                    # save by helix number and model, chain, symmetry to account for structures
                                    # that have multiple copies of the same chain
                                    f1 = n1.split("|")
                                    f2 = n2.split("|")

                                    if f1[1] == f2[1]:
                                        # same model, should always be the case already
                                        if len(f1) == 9:
                                            mcs1 = "|".join([f1[i] for i in [0,1,2,8]])
                                        else:
                                            mcs1 = "|".join([f1[i] for i in [0,1,2]])
                                        if len(f2) == 9:
                                            mcs2 = "|".join([f2[i] for i in [0,1,2,8]])
                                        else:
                                            mcs2 = "|".join([f2[i] for i in [0,1,2]])

                                        mcs1, mcs2 = sorted([mcs1, mcs2])
                                        helix_triple = (a1,a2,helix_number1,mcs1,mcs2)
                                        if not helix_triple in helix_to_number_location:
                                            helix_to_number_location[helix_triple] = (midpoint_x,midpoint_y,ArcRadius)
                                        else:
                                            # Compare ArcRadius of the existing tuple with the new one
                                            if ArcRadius > helix_to_number_location[helix_triple][2]:
                                                # Replace the existing tuple with the new one if ArcRadius is bigger
                                                helix_to_number_location[helix_triple] = (midpoint_x,midpoint_y,ArcRadius)

    return PSarcs, SVGarcs, num_arcs_drawn, arcs_drawn, helix_to_number_location


def make_chain_label(assemblies,requested_models,model,chain,chain_data,symmetry):
    if len(requested_models) > 1 or not requested_models[0] == '1':
        chain_label = "Chain %s|%s" % (model,chain)
    else:
        chain_label = "Chain %s" % chain

    if len(assemblies['ok_assemblies']) > 1:
        chain_label = "A%s %s" % (chain_data['assembly'],chain_label)

    if not symmetry == '1_555':
        chain_label += " " + symmetry

    return chain_label


def effective_width(s):
    """
    A number approximately proportional to the string width on the screen
    """
    w = len(s) + 0.4*sum(1 for x in s if x.isupper()) + 0.2*sum(1 for x in s if x.isdigit())

    return w


def break_line(s,target_width,maximum_width):
    lines = []
    if effective_width(s) > target_width:
        while len(s) > 0:
            if effective_width(s) > maximum_width:
                a = 0
                i = 0
                while a < target_width and i < len(s):
                    a += 1 + 0.4*s[i].isupper() + 0.2*s[i].isdigit()
                    if target_width - a < 10 and s[i] == " ":
                        # break at a space
                        a = target_width
                    i += 1
                lines.append(s[0:i].strip())
                s = s[i:]
            else:
                if s:
                    lines.append(s.strip())
                    s = ''
    else:
        lines.append(s)

    return lines


def draw_circular_diagram(chain_info, assemblies, filename, interaction_to_triple_list, params):
    """
    Construct the PostScript and SVG strings for the circular diagram
    """

    pdb_id = chain_info[0]['pdb_id']

    all_models = set([])
    requested_models = []

    chain_to_units = {}

    # show nucleotides that do not have xyz coordinates?
    n3d = params.get('n3d',True)
    chain_to_n3d_present = {}      # track if we are displaying nucleotides like (G128)

    unit_id_to_annotation = {}
    num_distinct_chains = len(set(x['chain_name'] for x in chain_info))

    unit_id_to_assemblies = {}
    valid_assembly_pairs = assemblies['valid_assembly_pairs']

    unit_id_symmetries_observed = set()

    # loop over the chains and retrieve the original experimental sequence
    # Read sequence of each chain and store sequence id and unit id
    # example: https://rna.bgsu.edu/rna3dhub/rest/SeqtoUnitMapping?ife=1S72|1|0
    # array of tuple(sequence positions, unit id) #('1S72|Sequence|0|A|2393', '1S72|1|0|A|2394')
    print('Downloading sequence position information about %d chains' % num_distinct_chains)
    success_count = 0
    for chain_data in chain_info:
        model = chain_data['model']
        if not model == '*' and not model in requested_models:
            requested_models.append(model)

        chain = chain_data['chain_name']
        chain_to_n3d_present[chain] = False

        if not chain in chain_to_units:
            chain_to_units[chain] = []

            chain_id = "%s|1|%s" % (pdb_id,chain)
            url = 'https://rna.bgsu.edu/rna3dhub/rest/SeqtoUnitMapping?ife=' + chain_id
            response = requests.get(url)
            if not response.status_code == 200:
                print("SeqtoUnitMapping request failed for chain %12s with status code: %s" % (chain,response.status_code))
            else:
                success_count += 1
                print("Sequence position to unit id request %3d successful for chain %s" % (success_count,chain_id))
                seq_pos_to_unit_id = response.text.strip().split("</br>")

                sequence_unit_list = []
                for line in seq_pos_to_unit_id:
                    pieces = line.split(" observed_as ")
                    if len(pieces) == 2:
                        sequence_unit_list.append(pieces)
                        unit_id = pieces[1]
                        if n3d and unit_id == "NULL":
                            chain_to_n3d_present[chain] = True
                        fields = unit_id.split("|")
                        if len(fields) > 2:
                            all_models.add(fields[1])
                        if len(fields) == 9:
                            unit_id_symmetries_observed.add(fields[8])

                chain_to_units[chain] = sequence_unit_list

                # some structures like 7AS5 have too many chains for there to be nucleotide annotations
                if num_distinct_chains < 100:
                    # download annotations of each nucleotide in this chain
                    url = "https://rna.bgsu.edu/correspondence/nucleotide_annotation?chain=%s" % chain_id
                    print("Trying %s" % url)
                    response = requests.get(url)
                    if response.status_code == 200:
                        annotation_text = response.text.strip()
                        for line in annotation_text.split("\n"):
                            fields = line.split("\t")

                            if len(fields) == 2:
                                if not fields[0] in unit_id_to_annotation:
                                    # only use the first annotation, that's how this service prioritizes
                                    unit_id_5 = first_five_fields(fields[0])
                                    unit_id_to_annotation[unit_id_5] = fields[1]

            if 'solitary' in chain_data:
                chain_to_units[chain] += [(x,x) for x in chain_data['solitary']]

    if len(requested_models) > 0:
        # make sure requested models exist
        requested_models = [m for m in requested_models if m in all_models]
    else:
        # sort all models numerically
        integer_models = set([])
        for m in all_models:
            try:
                integer_models.add(int(m))
            except:
                print('Non-integer model number %s found' % m)
        integer_models = sorted(integer_models)
        requested_models = [str(m) for m in integer_models]

    # track the assemblies, models, chains around the circle
    chain_name_list = []
    sequence_position_unit_id_pairs = []  # position of each unit around the circle
    max_nts_per_group = 0

    # determine how to name symmetries in some cases
    symmetry_type = ""
    if len(unit_id_symmetries_observed) > 0:
        symmetry_types_observed = set([x.split("_")[0] for x in unit_id_symmetries_observed])
        if "P" in symmetry_types_observed and not "ASM" in symmetry_types_observed:
            # this is the only case we need to recognize and account for at the moment
            symmetry_type = "P"

    if symmetry_type == "P":
        for chain_data in chain_info:
            if chain_data['symmetry'] and chain_data['symmetry'].startswith("ASM_"):
                chain_data['symmetry'] = chain_data['symmetry'].replace("ASM_","P_")

    # loop over all potential models
    for model in requested_models:
        # loop over all assemblies being shown, because chain-symmetry may be in more than one assembly
        for assembly in assemblies['ok_assemblies']:
            # repeat all chains over all models unless specific models were requested
            for chain_data in chain_info:
                if chain_data['assembly'] == assembly:
                    chain    = chain_data['chain_name']
                    symmetry = chain_data['symmetry']
                    symmetry_id = chain_data['symmetry_id']

                    # print('model:: %10s assembly %10s chain %10s symmetry_id %10s symmetry %10s' % (model,assembly,chain,symmetry_id,symmetry))

                    # match model if it is specified
                    if chain_data['model'] == '*' or chain_data['model'] == model:
                        unit_id_no_alt_set = set([])                   # keep simplest version of units

                        num_nts_found_this_chain = 0
                        for sequence_position, unit_id in chain_to_units[chain]:
                            if unit_id == "NULL" and not chain_to_n3d_present[chain]:
                                continue

                            # track which assembly(ies) each unit_id appears in
                            if unit_id in unit_id_to_assemblies:
                                unit_id_to_assemblies[unit_id].add(assembly)
                            else:
                                unit_id_to_assemblies[unit_id] = set(assembly)

                            unit_id_fields = unit_id.split("|")

                            if len(unit_id_fields) >= 5 and unit_id_fields[1] != model:
                                # different model than the requested model
                                continue

                            if len(unit_id_fields) < 9 and not symmetry == '1_555':
                                # no symmetry in the unit id is synonymous with symmetry 1_555
                                # no symmetry in the unit id does not match ASM_1 or P_1
                                continue

                            if len(unit_id_fields) == 9:
                                unit_id_symmetry = unit_id_fields[8]

                                if symmetry and symmetry.startswith('P_') and unit_id_symmetry.startswith('ASM_'):
                                    # 1A34 and other older structures used P_ instead of ASM_
                                    unit_id_symmetry = unit_id_symmetry.replace('ASM_','P_')

                                if not symmetry == unit_id_symmetry:
                                    # different symmetry operator than required
                                    continue

                            # accumulate unit ids with best alternate id, if any
                            # the url above returns rows sorted by alt id, then symmetry operator
                            if len(unit_id_fields) >= 7:     # alternate id may be present
                                unit_id_fields[6] = ""       # imagine the unit id without alt id
                                unit_id_no_alt = "|".join(unit_id_fields)
                                if unit_id_no_alt in unit_id_no_alt_set:  # full unit id
                                    continue                 # better alt id observed
                                if len(unit_id_fields) == 7: # no insertion code, no symmetry op
                                    unit_id_no_alt = "|".join(unit_id_fields[0:5])
                                    if unit_id_no_alt in unit_id_no_alt_set:  # better alt id observed
                                        continue
                            else:
                                unit_id_no_alt = unit_id

                            unit_id_no_alt_set.add(unit_id_no_alt) # record that we saw this one

                            # combined position is used as the key to angle around the circle
                            # it must have assembly because chains may be repeated
                            # it must have both sequence and unit id, because some sequence positions are
                            # repeated around the circle in different symmetries
                            # it must also have model in it for unobserved nucleotides that
                            # do not have a unit id with a model number
                            # it must have symmetry id because some structures like 8PZP use the
                            # symmetry operator 1_555 with position 6, 31, 56
                            combined_position = "%s %s %s %s %s" % (assembly,sequence_position,model,unit_id,symmetry_id)
                            sequence_position_unit_id_pairs.append((combined_position,unit_id))
                            num_nts_found_this_chain += 1

                        if num_nts_found_this_chain > 0:
                            # label for each chain around the outside of the circle
                            chain_label = make_chain_label(assemblies,requested_models,model,chain,chain_data,symmetry)

                            chain_name_list.append(chain_label)
                            sequence_position_unit_id_pairs.append(('chain_break','NULL'))

                            max_nts_per_group = max(max_nts_per_group,num_nts_found_this_chain)

    number_of_nucleotides = 0
    number_of_breaks = len(chain_name_list)
    for position, unit_id in sequence_position_unit_id_pairs:
        if not position == 'chain_break':
            number_of_nucleotides += 1

    if number_of_nucleotides == 0:
        print("No nucleotides found for %s" % (filename))
        return "", "", 612, 792

    print("Loaded %d nucleotides from %s" % (number_of_nucleotides,pdb_id))

    # set display properties set by the user
    show = params.get('show','')
    hide = params.get('hide','')
    dim  = params.get('dim','')
    text = params.get('text','basepair')

    coloring = params.get('coloring','default')
    if coloring == "grayscale":
        group_name_to_color = {
            "scale":         "grayscale",
            "lr-wc":         (0,0,0),           # black
            "nested-wc":     (0.2,0.2,0.2),     # dark gray
            "nested-non-wc": (0.4,0.4,0.4),     #
            "lr-non-wc":     (0.5,0.5,0.5),
            "stacking":      (0.8,0.8,0.8),     # light gray
            "bph":           (0.6,0.6,0.6),
            "br":            (0.7,0.7,0.7),
            "sr":            (0.75,0.75,0.75),
            "so":            (0.85,0.85,0.85),
            "near":          (0.9,0.9,0.9),     # very light gray
            "labels":        (0.5, 0.5, 0.5)    # medium gray
        }
    elif coloring == "wong":
        group_name_to_color = {
            "scale": "wong",
            "lr-wc":         (213.0/255,94.0/255,0.0/255),    # vermillion
            "nested-wc":     (0,114.0/255,178.0/255),         # blue
            "lr-non-wc":     (0.0/255,158.0/255,115.0/255),   # bluish green
            "nested-non-wc": (86.0/255,180.0/255,233.0/255),  # sky blue
            "stacking":      (240.0/255,228.0/255,66.0/255),  # yellow
            "bph":           (204.0/255,121.0/255,167.0/255), # reddish purple
            "br":            (230.0/255,159.0/255,0.0/255),   # orange
            "sr":            (0.5,0.5,0.5),       # gray
            "so":            (0.6,0.6,0.6),       # gray
            "near":          (0.7,0.7,0.7),
            "labels":        (86.0/255,180.0/255,233.0/255)   # sky blue
        }
    else:
        group_name_to_color = {
            "scale": "default",
            "lr-wc":         (1,0,0),           # red
            "nested-wc":     (0,0,1),           # blue
            "nested-non-wc": (0,1,1),           # cyan
            "lr-non-wc":     (0,0.8,0),         # green
            "stacking":      (0.95,0.95,0.1),   # yellow
            "bph":           (1,0.05,0.95),     # purple
            "br":            (1,0.6,0.1),       # orange
            "sr":            (0.5,0,0.5),       # purple
            "so":            (0,0.5,0.5),       # teal
            "near":          (0.7,0.7,0.7),     # gray
            "labels":        (86.0/255,180.0/255,233.0/255)   # sky blue
        }

    # dim colors of entire groups, which means to make the arcs more white
    for group_name in group_name_to_color:
        if group_name in dim:
            regular_color = group_name_to_color[group_name]
            factor = 0.8
            # move each component of the color toward 1, more so if factor is near 1
            colors = tuple(((1-factor)*component + factor) for component in regular_color)
            group_name_to_color[group_name] = colors

    # set numerical parameters based on total number of nucleotides
    # size of gaps between chains, as a multiple of the space for one nucleotide
    chain_gap_size = linear_interpolation(control_chain_gap_size, number_of_nucleotides)

    # make sure spaces allocated to gaps do not overwhelm nucleotides when many symmetries or models
    chain_gap_size_alternative = 0.2 * number_of_nucleotides / number_of_breaks
    if chain_gap_size > chain_gap_size_alternative:
        print("Reducing chain_gap_size from %0f to %0f" % (chain_gap_size,chain_gap_size_alternative))
        chain_gap_size = chain_gap_size_alternative

    effective_size = number_of_nucleotides+chain_gap_size*number_of_breaks

    print('number_of_nucleotides = %d' % number_of_nucleotides)
    print('number_of_breaks      = %d' % number_of_breaks)
    print('effective_size = number_of_nucleotides+chain_gap_size*number_of_breaks = %d' % (effective_size))

    # print('chain_gap_size = %0.2f; effectively the number of nucleotide slots per gap' % chain_gap_size)
    # gap_spaces = chain_gap_size * number_of_breaks
    # print('gap_spaces = %0.2f' % gap_spaces)

    # when the number of nucleotides is very large, increase the page size for PDF
    mul = 1
    if effective_size > 120000:
        # example: 4V3P
        mul = 10
    elif effective_size > 100000:
        mul = 8
    elif effective_size > 60000:
        mul = 6
    elif effective_size > 40000:
        # example: 7AS5
        mul = 5
    elif effective_size > 20000:
        # example: 6WLN
        mul = 4
    elif effective_size > 10000:
        mul = 3
    elif effective_size > 5500:
        # eukaryotic ribosomes are generally above this limit, bacterial below
        mul = 2
    if mul > 1:
        print("Using page size multiplier %d" % mul)

    # 8.5 by 11 inch page size, 72 points per inch
    page_width  = 612 * mul
    page_height = 792 * mul

    linewidth = linear_interpolation(control_linewidth, effective_size, 'decay') * mul    # interaction arc width
    circle_radius = linear_interpolation(control_circle_radius, effective_size) * mul
    tick_mark_length = linear_interpolation(control_tick_mark_length, effective_size, 'decay') * mul
    base_distance = 0.1 * mul  # added to circle_radius + tick_mark_length

    base_number_font_size = linear_interpolation(control_base_number_font_size, effective_size, 'decay') * mul

    major_number_distance_min = linear_interpolation(control_major_number_distance, effective_size, 'decay') * mul
    major_number_nt_cutoff = 120   # above this number of nucleotides, show major numbers around the outside
    major_number_interval = divide_and_round(number_of_nucleotides)
    major_number_font_size = linear_interpolation(control_major_number_font_size, effective_size, 'decay') * mul
    chain_font_size = linear_interpolation(control_major_number_font_size, effective_size, 'decay') * mul
    major_number_aligner_angle = major_number_font_size / (15 * mul)

    chain_gap_location = linear_interpolation(control_chain_gap_location, effective_size) # angular fix

    outer_arc_line_width = linewidth * 1.5

    if 'helix_size' in params:
        helix_size = params['helix_size'] * mul
    else:
        helix_size = linear_interpolation(control_helix_size, number_of_nucleotides, 'decay') * mul
    helix_radius_shift = helix_size / 2.0

    # assign angle to each sequence id and unit id
    angle_difference = 360.0 / effective_size
    sequence_position_to_angle = {}
    assembly_unit_id_to_angle = {}

    # make first chain label horizontal
    start_angle = -(-major_number_aligner_angle+chain_gap_location*(chain_gap_size+1)*angle_difference)
    angle = start_angle

    for sequence_position, unit_id in sequence_position_unit_id_pairs:
        if sequence_position == 'chain_break':
            angle = angle - (chain_gap_size * angle_difference)
        else:
            assembly = sequence_position.split(" ")[0]
            if not unit_id == 'NULL':
                assembly_unit_id_to_angle[(assembly, unit_id)] = angle
            sequence_position_to_angle[sequence_position] = angle
            angle = angle - angle_difference

    # assign standard nucleotide to each unit id to identify Watson-Crick AU, GC, GU pairs
    unit_id_to_standard = {}
    chain_to_longest_base_length = {}
    modified_base_to_parent = {}
    for (assembly,unit_id) in assembly_unit_id_to_angle.keys():
        if not unit_id in unit_id_to_standard:
            fields = unit_id.split("|")
            chain = fields[2]
            base = fields[3]
            if not chain in chain_to_longest_base_length:
                chain_to_longest_base_length[chain] = 1

            if base in ['A','C','G','U']:
                unit_id_to_standard[unit_id] = base
            elif base in ['DA','DC','DG','DT']:
                unit_id_to_standard[unit_id] = base
                chain_to_longest_base_length[chain] = max(chain_to_longest_base_length[chain],2)
            else:
                chain_to_longest_base_length[chain] = max(chain_to_longest_base_length[chain],len(base))
                if not modified_base_to_parent:
                    # load modified nucleotide information here, if available
                    lines = []

                    # download up to date mappings from modified to standard
                    url = 'https://rna.bgsu.edu/modified/nt_mappings.txt'
                    response = requests.get(url)
                    if response.status_code == 200:
                        lines = response.text.strip().split("\n")

                        for line in lines:
                            fields = line.strip().split("\t")
                            if len(fields) == 2:
                                modified_base_to_parent[fields[0]] = fields[1]

                unit_id_to_standard[unit_id] = modified_base_to_parent.get(base,base)

    # identify basepairs actually present in the structure to save time later
    basepair_interaction = ['cWW','cWw','cwW','tWW','cWH','cHW','tWH','tHW','cHH','cHh','chH','tHH','tHh','thH','cWS','cSW','tWS','tSW','cHS','cSH','tHS','tSH','cSS','cSs','csS','tSS','tSs','tsS']
    basepairs_present = []
    for bp in basepair_interaction:
        for alternative in ["","a"]:
            if bp+alternative in interaction_to_triple_list:
                basepairs_present.append(bp+alternative)
    near_basepairs_present = []
    for bp in basepair_interaction:
        for alternative in ["","a"]:
            if "n"+bp+alternative in interaction_to_triple_list:
                near_basepairs_present.append("n"+bp+alternative)

    # draw circular arcs for interactions
    arc_commands = []
    near_commands = []

    SVG_near = []
    SVG_true_arcs = []

    # record number of arcs drawn in each group
    group_name_to_count = {}
    for group_name in arc_group_names_by_order:
        group_name_to_count[group_name] = 0

    non_near_arcs_drawn = set()

    # map helix label to outermost arc of the helix
    helix_to_number_location = {}
    unit_shift = (linewidth*360)/(2*3.14*240*4)
    # draw any dim arcs first, then regular colored arcs
    for dim_level in ['dim','not_dim']:
        # loop over groups of interactions to draw colored arcs for
        for group_name in arc_group_names_by_order:
            if group_name in hide:
                continue
            if show and not group_name in show and not group_name in dim:
                continue
            interactions = arc_group_to_interactions[group_name]
            if group_name == "near":
                # draw gray arcs for all near interactions, to indicate additional close contacts
                angle_shifter = 0
                crossing = "all"
                base_combination = "all"
                # custom list of interactions
                interactions = []
                for interaction in interaction_to_triple_list.keys():
                    if interaction.startswith("n"):
                        interactions.append(interaction)
            elif group_name == "stacking":
                angle_shifter = -1.5*unit_shift
                crossing = "all"
                base_combination = "all"
            elif group_name == "br":
                angle_shifter = -3*unit_shift
                crossing = "all"
                base_combination = "all"
            elif group_name == "bph":
                angle_shifter = 3*unit_shift
                crossing = "all"
                base_combination = "all"
            elif group_name == "sr":
                angle_shifter = -2.25*unit_shift
                crossing = "all"
                base_combination = "all"
            elif group_name == "so":
                angle_shifter = 2.25*unit_shift
                crossing = "all"
                base_combination = "all"
            elif group_name == "nested-non-wc":
                angle_shifter = 1.5*unit_shift
                crossing = "nested"
                base_combination = "depends_on_interaction"
                interactions = basepairs_present
            elif group_name == "lr-non-wc":
                angle_shifter = 1.5*unit_shift
                crossing = "lr"
                base_combination = "depends_on_interaction"
                interactions = basepairs_present
            elif group_name == "nested-wc":
                angle_shifter = 0
                crossing = "nested"
                base_combination = "wc"
                interactions = list(set(basepairs_present) & set(['cWW','cwW','cWw','cWWa','cwWa','cWwa']))
            elif group_name == "lr-wc":
                angle_shifter = 0
                crossing = "lr"
                base_combination = "wc"
                interactions = list(set(basepairs_present) & set(['cWW','cwW','cWw','cWWa','cwWa','cWwa']))

            # set the color; skip if we are dimming now but group should not be dimmed, etc.
            if dim_level == "dim" and group_name in dim:
                colors = group_name_to_color[group_name]
            elif dim_level == "not_dim" and not group_name in dim:
                colors = group_name_to_color[group_name]
            else:
                continue

            # track arcs drawn in this interaction group
            if group_name == "near":
                # do not duplicate a regular interaction arc with a near interaction arc
                arcs_drawn = non_near_arcs_drawn
            else:
                arcs_drawn = set()

            # only set the color once
            gave_color_command = False

            # loop over all types of interaction in the group
            for interaction in interactions:
                if base_combination == "depends_on_interaction":
                    if "cww" in interaction.lower():
                        bc = "non_wc"
                    else:
                        bc = "all"
                else:
                    bc = base_combination
                pairs_and_crossing = interaction_to_triple_list[interaction]
                PSarcs, SVGarcs, num_arcs_drawn, arcs_drawn, helix_to_number_location = draw_arcs(pairs_and_crossing, assembly_unit_id_to_angle, arcs_drawn, crossing, bc, unit_id_to_standard, angle_shifter, circle_radius, unit_id_to_annotation, helix_to_number_location, helix_radius_shift, unit_id_to_assemblies, valid_assembly_pairs)
                if num_arcs_drawn > 0:
                    if not gave_color_command:
                        # now that we found an arc to draw, make sure the color is set exactly once
                        PSarcs = ["%f %f %f setrgbcolor" % colors] + PSarcs

                        svg_color = [int(255*c) for c in colors]
                        SVGarcs = ['<g fill="none" stroke="rgb(%d,%d,%d)" stroke-width="%f">' % (svg_color[0],svg_color[1],svg_color[2],linewidth)] + SVGarcs

                        gave_color_command = True
                    if group_name == "near":
                        near_commands += PSarcs  # draw near first
                        SVG_near += SVGarcs
                    else:
                        arc_commands += PSarcs
                        SVG_true_arcs += SVGarcs
                        non_near_arcs_drawn = non_near_arcs_drawn.union(arcs_drawn)

                    group_name_to_count[group_name] += num_arcs_drawn

            if gave_color_command:
                if group_name == "near":
                    SVG_near.append("</g>")
                else:
                    SVG_true_arcs.append("</g>")

    # accumulate a list of PostScript commands for the page
    PSlist = []
    SVGlist = []

    # tell SVG how to scale down the image for the inital view
    target_width = 612
    target_height = 792
    SVGlist.append('<?xml version="1.0" encoding="UTF-8" standalone="no"?>')
    SVGlist.append('<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="%d" height="%d" viewBox="0 0 %d %d">' % (target_width,target_height,page_width,page_height))
    SVGlist.append('<g id="zoomLayer">')

    # collect together header lines
    header = []
    title_text = filename.replace("R3DCID_","")
    header.append(title_text)
    if 'description' in params:
        if '\\n' in params['description']:
            header += params['description'].split('\\n')
        else:
            header += break_line(params['description'],130,145)
    if len(assemblies['message']) > 0:
        header.append(assemblies['message'])
    if 'title' in params:
        header += break_line(params['title'],130,145)
    if 'method' in params:
        header.append(params['method'])
    if 'source' in params:
        header.append(params['source'])
    if 'release_date' in params:
        header.append(params['release_date'])
    if 'resolution' in params:
        try:
            r = float(params['resolution'])
            header.append(params['resolution']+"A")
        except:
            pass

    header_font_size = 10 * mul
    header_vertical = 15 * mul
    for h in header:
        # approximate width; all uppercase titles are lots wider, digits are wider than lowercase
        w = effective_width(h)
        if w > 130:
            # scale down font size; could get very small
            header_font_size = min(10 * mul * 130.0 / w, header_font_size)
            header_vertical = max(10 * mul,header_vertical * 130.0 / w)
            print('New header_font_size %6.2f' % header_font_size)

    x_header = 50 * mul
    y_header = 760 * mul

    # PSlist.append("%d %d translate" % (x_header,y_header)) # move to top left of page
    PSlist.append("<< /PageSize [%d %d] >> setpagedevice" % (page_width,page_height))
    PSlist.append("newpath")                         # new drawing starting
    PSlist.append("0 setgray")                       # black color
    PSlist.append("/Times-Roman findfont")
    PSlist.append("%f scalefont" % header_font_size)
    PSlist.append("setfont")
    PSlist.append("0 0 moveto")                      # move the text cursor to the current origin

    SVGlist.append('<g font-family="Times-Roman" font-size="%f" fill="black">' % header_font_size)

    y = y_header
    for h in header:
        PSlist.append("%d %d moveto" % (x_header,y))    # move down y points from the origin, 72 points to the inch
        PSlist.append("(%s) show" % (h))    # show the header line

        SVGlist.append('<text x="%d" y="%d">%s</text>' % (x_header,page_height-y, h))

        y -= header_vertical

    # move origin back to lower left corner
    # PSlist.append("-50 -760 translate")

    # end the font choice
    SVGlist.append('</g>')

    # draw the table of arc colors and counts
    x_table = 50 * mul
    y_table = 160 * mul
    y_delta = -12 * mul
    x_offset = 15 * mul
    table_font_size = 10 * mul

    y = y_table + y_delta

    PSlist.append("newpath")
    PSlist.append("/Times-Roman findfont")
    PSlist.append("%f scalefont" % table_font_size)
    PSlist.append("0 setgray")
    PSlist.append("setfont")

    SVGlist.append('<g font-family="Times-Roman" font-size="%f" fill="black">' % table_font_size)

    # count separately, because sometimes n1 BPh n2 and also n2 BPh n1, which is undercounted
    for group_name in ["bph","br","sr","so"]:
        total = 0
        for interaction in arc_group_to_interactions[group_name]:
            total += len(interaction_to_triple_list[interaction])
        group_name_to_count[group_name] = max(total,group_name_to_count[group_name])

    for group_name in ["nested-wc","lr-wc","nested-non-wc","bonus","lr-non-wc","stacking","bph","br","sr","so","near"]:
        if not group_name in hide and not (group_name == "bonus" and "nested-non-wc" in hide):
            if not group_name == "bonus":
                PSlist.append("%d %d moveto" % (x_table+x_offset,y))
                PSlist.append("%f %f %f setrgbcolor" % group_name_to_color[group_name])
                PSlist.append("%d %f %d %d rectfill" % (x_table, y, 10 * mul, 7 * mul))

                r, g, b = group_name_to_color[group_name]
                SVGlist.append('<rect x="%d" y="%d" width="%d" height="%d" fill="rgb(%d,%d,%d)" />'
                                % (x_table, page_height - y - 7 * mul, 10 * mul, 7 * mul, int(255*r), int(255*g), int(255*b)))

            PSlist.append("%d %f moveto" % (x_table + x_offset, y))
            PSlist.append("0 0 0 setrgbcolor")
            if group_name == "bonus":
                PSlist.append("(%s) show" % (arc_group_name_to_text[group_name]))
                SVGlist.append('<text x="%d" y="%d">%s</text>' % (x_table+x_offset,page_height-y,arc_group_name_to_text[group_name]))
            else:
                PSlist.append("(%s) show" % (arc_group_name_to_text[group_name] % (group_name_to_count[group_name])))
                SVGlist.append('<text x="%d" y="%d">%s</text>' % (x_table+x_offset,page_height-y,arc_group_name_to_text[group_name] % (group_name_to_count[group_name])))

            y = y + y_delta

    # end the font choice
    SVGlist.append('</g>')

    # collect chain identifiers and display names, and break lines
    chains_printed = set([])
    chain_lines = []
    num_distinct_chains = len(set([x['chain_name'] for x in chain_info]))
    if num_distinct_chains > 11:
        # allow a longer line length because the text will be smaller
        max_length = int(60 * num_distinct_chains / 11.0)
    else:
        # reasonable line length
        max_length = 60
    for chain_data in chain_info:
        chain_name = chain_data['chain_name']
        if not chain_name in chains_printed:
            chains_printed.add(chain_name)
            display_name = chain_data['display_name']
            chain_text = "Chain %s: %s" % (chain_name, display_name)
            chain_lines += break_line(chain_text,max_length,max_length+5)

    # font size for table of chain names and their standardized names
    y = y_table + y_delta
    chain_table_font_size = table_font_size
    if len(chain_lines) > 11:
        # reduce spacing between rows
        y_delta = y_delta * 11.0 / len(chain_lines)
        chain_table_font_size = chain_table_font_size*11.0/len(chain_lines)
        PSlist.append("/Times-Roman findfont")
        PSlist.append("%f scalefont" % (chain_table_font_size))
        PSlist.append("setfont")

    SVGlist.append('<g font-family="Times-Roman" font-size="%f" fill="black">' % chain_table_font_size)

    # print list of chains and their display names below the diagram
    x_chain_list = 320 * mul
    for chain_text in chain_lines:
        ct = chain_text.replace("(","\(").replace(")","\)")
        PSlist.append("%d %f moveto" % (x_chain_list, y))
        PSlist.append("(%s) show" % (ct))

        SVGlist.append('<text x="%d" y="%d">%s</text>' % (x_chain_list,page_height-y,ct))

        y = y + y_delta

    # end the font choice
    SVGlist.append('</g>')

    # move origin to center of the circle
    PSlist.append("%d %d translate" % (306 * mul, 445 * mul))
    PSlist.append("newpath")
    PSlist.append("%s setlinewidth" % (linewidth))

    SVGlist.append('<g transform="translate(%d, %f)">' % (306 * mul, page_height - 445 * mul))

    # plan out which interactions to show around the outside of the circle and in what order
    interactions_by_importance = []
    if "all" in text:
        interactions_by_importance = basepairs_present + arc_group_to_interactions['bph']
        interactions_by_importance += arc_group_to_interactions['br'] + arc_group_to_interactions['stacking']
        interactions_by_importance += arc_group_to_interactions['sr'] + arc_group_to_interactions['so']
        for interaction in interaction_to_triple_list.keys():
            if interaction.startswith("n"):
                interactions_by_importance.append(interaction)
    else:
        if not text or "basepair" in text:
            interactions_by_importance = basepairs_present
        if "bph" in text:
            interactions_by_importance = interactions_by_importance + arc_group_to_interactions['bph']
        if "br" in text:
            interactions_by_importance = interactions_by_importance + arc_group_to_interactions['br']
        if "stack" in text:
            interactions_by_importance = interactions_by_importance + arc_group_to_interactions['stacking']
        if "sr" in text:
            interactions_by_importance = interactions_by_importance + arc_group_to_interactions['sr']
        if "so" in text:
            interactions_by_importance = interactions_by_importance + arc_group_to_interactions['so']
        if "near" in text:
            for interaction in interaction_to_triple_list.keys():
                if interaction.startswith("n"):
                    if (not text or "basepair" in text) and interaction in near_basepairs_present:
                        interactions_by_importance.append(interaction)
                    if "stack" in text and interaction[1:] in ['s33','s35','s53','s55']:
                        interactions_by_importance.append(interaction)
                    if "bph" in text and "BPh" in interaction:
                        interactions_by_importance.append(interaction)
                    if "br" in text and "BR" in interaction:
                        interactions_by_importance.append(interaction)
                    if "sr" in text and interaction in arc_group_to_interactions["sr"]:
                        interactions_by_importance.append(interaction)
                    if "so" in text and interaction[1:] in arc_group_to_interactions["so"]:
                        interactions_by_importance.append(interaction)


    # collect annotations to put outside of the labels of the base + number for each nt
    text_interactions = {}

    for interaction in interactions_by_importance:
        for n1, n2, crossing_number in interaction_to_triple_list[interaction]:
            # model is: n1 interaction n2
            # don't show self annotations, which omits many 0BPh interactions
            if n1 == n2:
                continue

            for assembly in assemblies['ok_assemblies']:
                an1 = (assembly,n1)
                an2 = (assembly,n2)

                # skip nucleotides that have no place around the circle
                if not an1 in assembly_unit_id_to_angle:
                    continue
                if not an2 in assembly_unit_id_to_angle:
                    continue

                angle1 = assembly_unit_id_to_angle[an1]
                angle2 = assembly_unit_id_to_angle[an2]

                n1_fields = n1.split("|")
                if len(n1_fields) > 2:
                    n1_chain = n1.split("|")[2]
                else:
                    n1_chain = ''

                n2_fields = n2.split("|")
                if len(n2_fields) > 2:
                    n2_chain = n2.split("|")[2]
                else:
                    n2_chain = ''

                if abs(angle1) < 90 or abs(angle1) > 270:
                    # n1 on right side of the circle, append interaction like n1 BPh n2
                    t = interaction + " " + short_id(n2,n1_chain)
                    if an1 in text_interactions:
                        if not t in text_interactions[an1]:
                            text_interactions[an1].append(t)
                    else:
                        text_interactions[an1] = [t]

                if abs(angle2) >= 90 and abs(angle2) <= 270:
                    # n2 on left side of the circle, prepend over here
                    t = short_id(n1,n2_chain) + " " + interaction
                    if an2 in text_interactions:
                        if not t in text_interactions[an2]:
                            text_interactions[an2].insert(0,t)
                    else:
                        text_interactions[an2] = [t]

                # BPh, BR, sO, cSR, tSR are listed only once in the list of interactions
                # recognize these, produce a "reversed" version, and list on both sides
                if "BPh" in interaction or "BR" in interaction or "s3O" in interaction or "s5O" in interaction or "cSR" in interaction or "tSR" in interaction:

                    # things like 3BPh and 6BR should become 3PhB and 6RB (even though those are weird)
                    if interaction[0] == "n":
                        # for example n3BPh
                        rev_interaction = interaction[0]+interaction[1]+interaction[3:]+interaction[2]
                    else:
                        # for example tHS, 3BPh, sO4'3
                        rev_interaction = interaction[0]+interaction[2:]+interaction[1]

                    if abs(angle1) >= 90 and abs(angle1) <= 270:
                        # n1 on left side of the circle, prepend reversed interaction over here
                        t = short_id(n2,n1_chain) + " " + rev_interaction
                        if an1 in text_interactions:
                            if not t in text_interactions[an1]:
                                text_interactions[an1].insert(0,t)
                        else:
                            text_interactions[an1] = [t]

                    if abs(angle2) < 90 or abs(angle2) > 270:
                        # n2 on right side of the circle, append interaction
                        t = rev_interaction + " " + short_id(n1,n2_chain)
                        if an2 in text_interactions:
                            if not t in text_interactions[an2]:
                                text_interactions[an2].append(t)
                        else:
                            text_interactions[an2] = [t]

    # track text and lengths of individual nucleotide labels
    sequence_id_to_text = {}
    outside_arcs = []    # list of starting and stopping angles for outside black arcs for each chain
    radius_base = circle_radius + outer_arc_line_width/2 + tick_mark_length*1.2 + base_distance

    if max_nts_per_group < major_number_nt_cutoff or text == 'none':
        show_major_numbers = False
    else:
        show_major_numbers = True

    # write features around the outside of the circle
    for feature in range(3):

        PSlist.append("gsave")  # save the current rotation

        if feature == 0:
            # prepare to draw tick lines in gray
            PSlist.append("0.5 setgray")
            SVGlist.append('<g stroke="gray" stroke-width="%s">' % linewidth)

        if feature == 1:
            # modify base_number_font_size if text labels are too long
            all_text_lengths = [len(x) for x in sequence_id_to_text.values()]
            p_index = min(int(0.95 * len(all_text_lengths)),len(all_text_lengths)-1)
            text_length_95 = sorted(all_text_lengths)[p_index]
            p_index = min(int(0.98 * len(all_text_lengths)),len(all_text_lengths)-1)
            text_length_98 = sorted(all_text_lengths)[p_index]

            text_radius = circle_radius + outer_arc_line_width/2 + tick_mark_length*1.2 + base_distance
            if text_radius + 0.6 * text_length_98 * base_number_font_size > page_width / 2:
                # reduce font size for base labels
                print('base_number_font_size before %8.2f' % base_number_font_size)
                base_number_font_size = (page_width/2 - text_radius)/(text_length_98 * 0.6)
                print('base_number_font_size after  %8.2f' % base_number_font_size)

            # set the shift angle for base-number labels, to get them centered correctly
            base_number_angle = base_number_font_size / (15 * mul)

        if feature == 2:
            # draw text of base, number, helix number, interactions around the outside of the circle
            PSlist.append("0 setgray")
            # PSlist.append("/Consolas findfont")  # for writing monomers or base; Consolas N/A on server
            PSlist.append("/NimbusMonoPS-Regular findfont")  # for writing monomers or base
            PSlist.append("%f scalefont" % base_number_font_size)
            PSlist.append("setfont")

            SVGlist.append('<g font-family="Courier New, Courier, monospace" font-size="%f" fill="black" stroke="none">' % base_number_font_size)

        if feature == 1:
            # update distance of major numbers and chain labels from circle radius
            max_major_number_radius = 250 * mul

            if show_major_numbers:
                major_number_distance = max(major_number_distance_min, min(5 * mul,major_number_distance_min) + 0.6 * text_length_95 * base_number_font_size)
                # draw large numbers every 50 or 100 around the outside of the circle
                gntc = group_name_to_color["labels"]
                PSlist.append("%f %f %f setrgbcolor" % gntc)
                PSlist.append("/Times-Roman findfont")  # sets font info for numbers
                PSlist.append("%f scalefont" % major_number_font_size)
                PSlist.append("setfont")

                SVGlist.append('<g font-family="Times-Roman" font-size="%f" fill="rgb(%f,%f,%f)" stroke="none">' % (major_number_font_size,int(255*gntc[0]),int(255*gntc[1]),int(255*gntc[2])))

            elif effective_size > 220:
                major_number_distance = max(major_number_distance_min, min(5 * mul,major_number_distance_min) + 0.6 * text_length_95 * base_number_font_size)
            else:
                major_number_distance = major_number_distance_min

        # variables for perimeter circle
        current_arc_start = start_angle
        current_arc_end   = start_angle + 0.001  # just in case

        large_number_track = 0
        previous_helix_length = 2
        helix = ""

        longest_base_length = 1
        for sequence_id,unit_id in sequence_position_unit_id_pairs:
            if sequence_id == 'chain_break':
                if feature == 0:
                    # record the starting and ending angles of the chain that just ended
                    outside_arcs.append([current_arc_start,current_arc_end])
                current_arc_start = None    # next nucleotide will start the next chain
                continue

            # sequence_id contains information on assembly, model, chain, symmetry, etc.
            a = sequence_position_to_angle[sequence_id]

            if current_arc_start == None:
                current_arc_start = a
            current_arc_end = a

            if feature == 0:
                # draw tick mark at each nucleotide
                if not unit_id == "NULL":
                    co = math.cos(math.radians(a))
                    si = math.sin(math.radians(a))

                    r1 = circle_radius + outer_arc_line_width / 2
                    x1 = r1 * co
                    y1 = r1 * si

                    r2 = r1 + tick_mark_length
                    x2 = r2 * co
                    y2 = r2 * si

                    PSlist.append('%f %f newpath moveto' % (r1*co,r1*si))
                    PSlist.append("%f %f lineto" % (r2*co,r2*si))
                    PSlist.append("stroke")

                    SVGlist.append('<line x1="%f" y1="%f" x2="%f" y2="%f" />' % (x1, -y1, x2, -y2))

                # accumulate base, number, helix, interaction text so it can be scaled right
                current_chain = sequence_id.split("|")[2]
                longest_base_length = chain_to_longest_base_length.get(current_chain,1)
                n3d_present = chain_to_n3d_present.get(current_chain,False)

                if text == 'none':
                    t = ""
                elif unit_id == "NULL":
                    if n3d_present:
                        sfields = sequence_id.split(" ")[1].split("|")
                        t = "(%s%s)" % (sfields[3],sfields[4])  # display base and sequence number
                        s = longest_base_length
                        # less space needed because of (
                        while s >= 2:
                            t = " " + t
                            s = s - 1
                    else:
                        t = " "
                else:
                    t = short_id(unit_id,'',chain_to_longest_base_length)   # base, number, insertion code
                    if n3d_present:
                        # add spaces to align with parentheses in (sequence number)
                        t = " " + t + " "
                    fields = unit_id.split("|")

                # print base, number, interactions and helix number as requested
                assembly = sequence_id.split(" ")[0]
                if (assembly,unit_id) in text_interactions:
                    interactions = ", ".join(text_interactions[(assembly,unit_id)])
                else:
                    interactions = ""

                if 'helix' in text or text == 'all':
                    if unit_id_to_annotation:
                        unit_id_5 = first_five_fields(unit_id)
                        helix = unit_id_to_annotation.get(unit_id_5,None)
                        if helix == None or unit_id == "NULL":
                            helix = " " * previous_helix_length
                        else:
                            helix = helix.replace("Helix ","").replace("helix ","").strip()
                            if helix.startswith("5S"):
                                helix = helix.replace("5S","")
                            previous_helix_length = len(helix)

                if abs(float(a)) >= 90 and abs(float(a)) <= 270:
                    # left side of the circle (that's why it's interaction + " " + t + " " + helix)
                    if unit_id == "NULL" and n3d_present:
                        all_text = interactions+" "+t
                    else:
                        all_text = interactions+" "+t
                    if len(helix) > 0:
                        all_text = all_text + " " + helix
                else:
                    # right side of the circle
                    if len(helix) > 0:
                        all_text = helix+" "+t+" "+interactions
                    else:
                        all_text = t+" "+interactions

                sequence_id_to_text[sequence_id] = all_text

            if feature == 2:
                # draw base, number, helix, interactions text around the outside of the circle
                all_text = sequence_id_to_text[sequence_id]

                PSlist.append("gsave")
                if abs(float(a)) >= 90 and abs(float(a)) <= 270:
                    # left side of the circle
                    # use radius because it renders faster than stringwidth
                    radius = radius_base + 0.6*len(all_text)*base_number_font_size
                    PSlist.append("%s rotate"%(a+base_number_angle+180))  # rotate text 180 degrees
                    PSlist.append("-%0.2f 0 moveto" % (radius))  # negative radius for right location
                    PSlist.append("(%s) show" % (all_text))

                    SVGlist.append('<g transform="rotate(%f)">' % (-(a + base_number_angle + 180)))
                    SVGlist.append('<text x="%0.2f" y="0" xml:space="preserve">%s</text>' % (-radius, all_text))
                    SVGlist.append('</g>')
                else:
                    # right side of the circle; left justify text which is easier
                    radius = radius_base
                    PSlist.append("%s rotate" % (a-base_number_angle))
                    PSlist.append("%0.2f 0 moveto" % (radius))
                    PSlist.append("(%s) show" % (all_text))

                    SVGlist.append('<g transform="rotate(%f)">' % (-(a - base_number_angle)))
                    SVGlist.append('<text x="%0.2f" y="0" xml:space="preserve">%s</text>' % (radius, all_text))
                    SVGlist.append('</g>')

                PSlist.append("grestore")

            if feature == 1 and show_major_numbers and not text == 'none':
                # add large numbers around outside of the circle, all with the same font and color

                if unit_id and not unit_id == "NULL":
                    n = unit_id.split("|")[4]
                else:
                    continue

                if 0 <= current_arc_start - a < 1:
                    # within 1 degree of the start of the chain, too close, number would overlap
                    continue

                if large_number_track == n:
                    # prevents showing the same number twice in a row
                    continue

                large_number_track = n
                if n != "" and int(n) % major_number_interval == 0 :
                    radius = min(circle_radius + major_number_distance,max_major_number_radius)

                    PSlist.append("gsave")
                    PSlist.append("newpath")
                    if abs(float(a)) >= 90 and abs(float(a)) <= 270:
                        # left half of the circle, right justify
                        PSlist.append("%s rotate" % (a + major_number_aligner_angle + 180))
                        PSlist.append("%f 0 moveto" % (-radius))
                        PSlist.append("(%s) dup stringwidth pop neg 0 rmoveto show" % (n))  # right justify

                        SVGlist.append('<g transform="rotate(%f)">' % (-(a + major_number_aligner_angle + 180)))
                        SVGlist.append('<text x="%f" y="0" text-anchor="end">%s</text>' % (-radius, n))
                        SVGlist.append('</g>')
                    else:
                        # right half of the circle, left justify
                        PSlist.append("%s rotate" % (a - major_number_aligner_angle))
                        PSlist.append("%f 0 moveto" % (radius))
                        PSlist.append("(%s) show" % (n))  # left justify

                        SVGlist.append('<g transform="rotate(%f)">' % (-(a - major_number_aligner_angle)))
                        SVGlist.append('<text x="%f" y="0" text-anchor="start">%s</text>' % (radius, n))
                        SVGlist.append('</g>')
                    PSlist.append("stroke")
                    PSlist.append("grestore")

        # end line color or font size
        if not feature == 1 or show_major_numbers:
            SVGlist.append('</g>')

    # add the arcs for near interactions first so they are underneath
    PSlist += near_commands
    # draw the arcs for stacking, then base-backbone, then basepairs
    PSlist += arc_commands

    SVGlist += SVG_near
    SVGlist += SVG_true_arcs

    # outside_arcs and chain_name_list need to have the same length or the next step fails
    if not len(outside_arcs) == len(chain_name_list):
        print('Problem: outside_arcs and chain_name_list should have the same length, but do not')

    # draw chain labels outside the circle
    if not text == 'none':
        max_chain_label = max([len(x) for x in chain_name_list])
        if max_chain_label > 10:
            chain_label_font_size = chain_font_size * 10.0 / max_chain_label
        else:
            chain_label_font_size = chain_font_size

        gntc = group_name_to_color["labels"]
        PSlist.append("%f %f %f setrgbcolor" % gntc)
        PSlist.append("/Times-Roman findfont")
        PSlist.append("%f scalefont" % chain_label_font_size)
        PSlist.append("setfont")

        SVGlist.append('<g font-family="Times-Roman" font-size="%f" fill="rgb(%f,%f,%f)" stroke="none">' % (chain_label_font_size,int(255*gntc[0]),int(255*gntc[1]),int(255*gntc[2])))

        for e in range(len(outside_arcs)):
            # keep the chain label from going off the page
            radius = min(circle_radius + major_number_distance,max_major_number_radius)
            chain_label = chain_name_list[e]
            PSlist.append("gsave")
            if abs(outside_arcs[e][0]) >= 90 and abs(outside_arcs[e][0]) <= 270:
                # left side of the circle, right justify
                ang = outside_arcs[e][0]+major_number_aligner_angle+chain_gap_location*(chain_gap_size+1)*angle_difference+180
                PSlist.append("%s rotate" % (ang))
                PSlist.append("%0.3f 0 moveto" % (-radius))
                PSlist.append("(%s) dup stringwidth pop neg 0 rmoveto show" % chain_label)

                SVGlist.append('<g transform="rotate(%f)">' % (-ang))
                SVGlist.append('<text x="%f" y="0" text-anchor="end">%s</text>' % (-radius, chain_label))
                SVGlist.append('</g>')
            else:
                # right side of the circle, left justify
                ang = outside_arcs[e][0]-major_number_aligner_angle+chain_gap_location*(chain_gap_size+1)*angle_difference
                PSlist.append("%s rotate" % (ang))
                PSlist.append("%0.3f 0 moveto" % radius)
                PSlist.append("(%s) show" % chain_label)

                SVGlist.append('<g transform="rotate(%f)">' % (-ang))
                SVGlist.append('<text x="%f" y="0" text-anchor="start">%s</text>' % (radius, chain_label))
                SVGlist.append('</g>')

            PSlist.append("newpath")
            PSlist.append("grestore")

        # end chain label font
        SVGlist.append('</g>')

    if helix_size > 0:
        # draw helix numbers at the tops of the cWW arcs
        PSlist.append("0 0 0 setrgbcolor")
        PSlist.append("/Helvetica findfont %f scalefont setfont" % helix_size)  # Set font size

        SVGlist.append('<g font-family="Helvetica" font-size="%f" fill="black" stroke="none">' % (helix_size))

        for helix_triple, number_location in helix_to_number_location.items():
            a1, a2, helix_number, chain1, chain2 = helix_triple
            h_fields = helix_number.split()
            if len(h_fields) == 2:
                helix_text = h_fields[1]
                if helix_text.startswith("5S"):
                    helix_text = helix_text.replace("5S","")

                # center horizontally and vertically on midpoint location
                midpoint_x = number_location[0] - 0.5 * 0.4 * helix_size * len(helix_text)
                midpoint_y = number_location[1] - 0.5 * 0.4 * helix_size

                PSlist.append("%f %f moveto" % (midpoint_x, midpoint_y))
                PSlist.append("(%s) show" % (helix_text))
                PSlist.append("stroke")

                SVGlist.append('<text x="%f" y="%f">%s</text>' % (midpoint_x, -midpoint_y, helix_text))

        SVGlist.append('</g>')


    # draw black perimeter circular arc for each chain
    PSlist.append("newpath")
    PSlist.append("0 setgray")
    PSlist.append("%f setlinewidth" % (outer_arc_line_width))

    SVGlist.append('<g fill="none" stroke="black" stroke-width="%f">' % outer_arc_line_width)

    for start_angle, end_angle in outside_arcs:
        if start_angle is not None and end_angle is not None:
            PSlist.append("0 0 %f %f %f arc" % (circle_radius+outer_arc_line_width/4.0,end_angle - 0.25*angle_difference,start_angle + 0.25*angle_difference))
            PSlist.append("stroke")
            PSlist.append("newpath")

            sa = start_angle + 0.25 * angle_difference
            ea = end_angle - 0.25 * angle_difference

            # Convert degrees to radians for math
            sa_rad = math.radians(sa)
            ea_rad = math.radians(ea)

            r = circle_radius + outer_arc_line_width / 2.0

            # Compute start and end points
            x1 = r * math.cos(sa_rad)
            y1 = r * math.sin(sa_rad)
            x2 = r * math.cos(ea_rad)
            y2 = r * math.sin(ea_rad)

            # Determine if arc > 180° for large-arc-flag
            # Determine if arc > 180° for large-arc-flag
            large_arc_flag = 1 if (sa - ea) % 360 > 180 else 0
            sweep_flag = 1

            # SVG arc command
            SVGlist.append('<path d="M %f %f A %f %f 0 %d %d %f %f"/>' % (x1, -y1, r, r, large_arc_flag, sweep_flag, x2, -y2))

    SVGlist.append('</g>')

    PSlist.append("stroke")
    PSlist.append("showpage")

    # end the translation
    SVGlist.append('</g>')

    # end the zoomlayer
    SVGlist.append('</g>')

    # end the SVG commands
    SVGlist.append("</svg>")

    PScommands = "\n".join(PSlist)
    SVGcommands = "\n".join(SVGlist)

    return PScommands, SVGcommands, page_width, page_height


def get_filename(input_chains, params = {}):
    """
    Quickly process the inputs to produce the filename
    """

    # process pdb_id, assemblies, models, symmetries
    pdb_id, filename, requested_assemblies, requested_models, requested_model_chain, requested_chains, requested_symmetries = process_input_chains(input_chains,params)

    # process additional display parameters
    params, filename = set_parameters_from_input(params,filename,pdb_id)

    print("filename: %s" % filename)

    return filename


def main(input_chains, params = {}):
    """
    input_chains is a string like 4V9F or 4V9F|1|0 or 8GLP|1|L5+8GLP|1|L8 or 8QO5|3,8QO5|4
    input_chains is a text string telling PDB id, assemblies, chains, etc.
    output_directory tells where to write .ps and .pdf files
    params is a dictionary which optionally sets various display parameters
    """

    output_directory = params.get('output_directory','r3dcid_output')

    # make sure the output path exists
    Path(output_directory).mkdir(parents=True, exist_ok=True)

    # process pdb_id, assemblies, models
    pdb_id, filename, requested_assemblies, requested_models, requested_model_chain, requested_chains, requested_symmetries = process_input_chains(input_chains,params)

    # process additional display parameters
    params, filename = set_parameters_from_input(params,filename,pdb_id)

    if len(filename) == 0:
        print("Something went wrong with processing the chain(s)")
        return "", "Problem setting parameters from input"

    # when generating examples, remove "description" from the filename
    if True and "zirbel" in os.getcwd().lower():
        filename = filename.replace("_description","")

    # load a dictionary that maps interaction type to lists of nucleotide pairs
    interaction_to_triple_list, message = get_fr3d_interaction_to_triple_list(pdb_id,params.get('data_directory',''))

    if message:
        print(message)
        return "", message

    # get information on the chains in the structure and set the order of chains around the circle
    chain_info, assemblies, message = order_chains_around_diagram(pdb_id, requested_assemblies, requested_models, requested_model_chain, requested_chains, requested_symmetries, interaction_to_triple_list)

    if message:
        print(message)
        return "", message

    if 'header' in params:
        params, message = get_header_information(params,pdb_id)
        if message:
            print(message)

    # generate PostScript and SVG commands to make the diagram
    PScommands, SVGcommands, page_width, page_height = draw_circular_diagram(chain_info, assemblies, filename, interaction_to_triple_list, params)

    if len(PScommands) == 0:
        return "", "Empty PostScript command variable"

    if 'format' in params and 'svg' in params['format'].lower():
        path_filename = os.path.join(output_directory,filename + ".svg")

        with open(path_filename, 'w') as svgfile:
            svgfile.write(SVGcommands)
        print("Wrote output file %s.svg" % filename)

    if 'format' in params and 'pdf' in params['format'].lower():
        path_filename = os.path.join(output_directory,filename + ".ps")
        path_filename_pdf = os.path.join(output_directory,filename + ".pdf")

        with open(path_filename, 'w') as newpsfile:
            newpsfile.write(PScommands)

        try:
            if shutil.which("ps2pdf"):
                subprocess.run(
                    ["ps2pdf", path_filename, path_filename_pdf],
                    check=True,
                )
            elif shutil.which("gs"):
                subprocess.run(
                    [
                        "gs",
                        "-sDEVICE=pdfwrite",
                        "-sPAPERSIZE=custom ",
                        "-dDEVICEWIDTHPOINTS=%d " % page_width,
                        "-dDEVICEHEIGHTPOINTS=%d " % page_height,
                        "-dFIXEDMEDIA ",
                        "-dNOPAUSE",
                        "-dBATCH",
                        "-dSAFER",
                        f"-sOutputFile={path_filename_pdf}",
                        path_filename,
                    ],
                    check=True,
                )
            else:
                raise RuntimeError("Neither gs nor ps2pdf found on system")

        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"PDF conversion failed: {e}")

        print("Wrote output file %s.pdf" % filename)

        # remove the .ps file
        os.remove(path_filename)

    # return filename with no extension and empty message
    return filename, ""


if __name__ == '__main__':

    if len(sys.argv) > 1:
        # process command line arguments

        # command line examples:
        # python311 r3dcid.py 4TNA
        # python311 r3dcid.py "7JQQ|1|K+7JQQ|1|L+7JQQ|1|M+7JQQ|1|N+7JQQ|1|O"
        # python311 r3dcid.py "7K00|1|b+7K00|1|a+7K00|1|A+7K00|1|X+7K00|1|Y+7K00|1|Z+7K00|1|5"

        import argparse

        # allow user to specify input and output paths
        parser = argparse.ArgumentParser()
        parser.add_argument('structure', type=str, nargs='+', help='PDB id like 7EZ2, id with chain like 4V9F|1|9, id with model like 8QO5|2, multiple chains like 4V9F|1|0,9, multiple ids like 7EZ2;8GLP;7K00')
        parser.add_argument("--assemblies", help='Assemblies to show, comma separated list like 1,2,3')
        parser.add_argument("--symmetries", help='Symmetries to show, comma separated, values like 1_555 for the default, 2_655, P_4, ASM_6, etc.')
        parser.add_argument("--data", help="Location for data files downloaded by the program")
        parser.add_argument("--output", help="Location for diagrams produced by the program")
        parser.add_argument("--coloring", help='Arc colors, choose from ([default],wong,grayscale)')
        parser.add_argument("--format", help='Output format, comma separated list from (pdf,svg)')
        parser.add_argument("--show", help='Arcs to show, comma separated list from (nested-wc,lr-wc,wc,nested-non-wc,lr-non-wc,non-wc,stacking,bph,br,sr,so,near)')
        parser.add_argument("--hide", help='Arcs to hide, comma separated list from (nested-wc,lr-wc,wc,nested-non-wc,lr-non-wc,non-wc,stacking,bph,br,sr,so,near)')
        parser.add_argument("--dim", help='Arcs to dim, comma separated list from (nested-wc,lr-wc,wc,nested-non-wc,lr-non-wc,non-wc,stacking,bph,br,sr,so,near)')
        parser.add_argument("--text", help='Text to show outside the circle, comma separated list from ([basepair],stacking,bph,br,sr,so,near,helix,all,none)')
        parser.add_argument("--n3d", help='Show nucleotides with no 3D coordinates ([true] or false)')
        parser.add_argument("--helix_size", help='Font size for helix numbers when available, 0 for no helix numbers')
        parser.add_argument("--header", help='PDB information to show from (title,method,release_date,source,resolution,[all],none)')
        parser.add_argument("--description", help='Text description to show along the top of the diagram')
        args = parser.parse_args()

        params = {}
        if args.structure:
            requested_structure = args.structure[0].replace(" ","")
        else:
            print('Specify a structure to display, like 4TNA')

        if args.assemblies:
            params['assembly'] = args.assemblies.replace(" ","")

        if args.symmetries:
            params['symmetry'] = args.symmetries.replace(" ","")

        if args.data:
            data_path = args.data
        else:
            try:
                from r3dcid_configuration import data_directory
                data_path = data_directory
            except:
                data_path = "r3dcid_data"

        if not os.path.exists(data_path):
            print('Making data path %s' % data_path)
            Path(data_path).mkdir(parents=True, exist_ok=True)
        params['data_directory'] = data_path

        if args.output:
            output_directory = args.output
        else:
            try:
                from r3dcid_configuration import output_directory
            except:
                output_directory = "r3dcid_output"
        if not os.path.exists(output_directory):
            print('Making output path %s' % output_directory)
            Path(output_directory).mkdir(parents=True, exist_ok=True)
        params['output_directory'] = output_directory

        if args.coloring:
            params['coloring'] = args.coloring.lower()

        if args.format:
            params['format'] = args.format.lower().replace(" ","")
            if not 'pdf' in params['format'] and not 'svg' in params['format']:
                print('Producing pdf formatted output')
                params['format'] = 'pdf'
        else:
            params['format'] = 'pdf'

        if args.hide:
            params['hide'] = args.hide.lower().replace(" ","")

        if args.show:
            params['show'] = args.show.lower().replace(" ","")

        if args.dim:
            params['dim'] = args.dim.lower().replace(" ","")

        if args.text:
            if not args.text.lower() == 'basepair':
                params['text'] = args.text.lower().replace(" ","")

        if args.n3d and args.n3d.lower() == 'false':
            params['n3d'] = False

        if args.helix_size:
            try:
                hs = int(args.helix_size)
            except:
                hs = -1
            if hs >= 0:
                params['helix_size'] = hs

        if args.description:
            params['description'] = args.description

        if args.header:
            params['header'] = args.header.lower().replace("-","_").replace(" ","")

        # repeat the same parameters for semicolon-separated list of pdb ids
        for pdbmodelchain in requested_structure.split(";"):
            filename = main(pdbmodelchain, params)

            if len(filename) == 0:
                print('Something went wrong with %s' % pdbmodelchain)
