# utility functions for reading cif files for nucleotides

import pdbx
import os
import sys

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


def get_cif_data(base):
    """
    Read the .cif file and organize its data into data structures for this program
    """

    cif_data = read_monomer_cif(base)

    data = {}
    data['standard_base'] = []  # empty list when no standard base
    data['changes'] = []  # empty list when no atom changes

    data['pdb'] = {}
    data['pdb']['name'] = cif_data['chem_comp'][0].get('name',None)
    data['pdb']['chem_comp_type'] = cif_data['chem_comp'][0].get('type',None)
    data['pdb']['par_comp_id'] = cif_data['chem_comp'][0].get('mon_nstd_parent_comp_id','No mon_nstd_parent_comp_id line')
    data['pdb']['one_letter_code'] = cif_data['chem_comp'][0].get('one_letter_code','No one_letter_code line')
    data['pdb']['pdbx_initial_date'] = cif_data['chem_comp'][0].get('pdbx_initial_date','1900-01-01')

    for descriptor in cif_data['pdbx_chem_comp_descriptor']:
        t = descriptor['type']
        if t == "SMILES_CANONICAL":
            if not t in data:
                data['pdb'][t] = {}
            program = descriptor['program']
            if "openeye" in program.lower():
                data['pdb'][t][program] = {}
                data['pdb'][t][program]['program_version'] = descriptor['program_version']
                data['pdb'][t][program]['descriptor'] = descriptor['descriptor']

    # for k,v in data.items():
    #     print(k,v)
    # input("Press enter to continue")

    coordinates = {}
    coordinates_ideal = {}
    connections = []
    atom_to_element = {}
    atom_to_chirality = {}

    for row in cif_data['chem_comp_atom']:
        atom = row['atom_id']
        x = row['model_Cartn_x']
        y = row['model_Cartn_y']
        z = row['model_Cartn_z']

        if x and y and z:
            coordinates[atom] = [float(x),float(y),float(z)]
        atom_to_element[atom] = row['type_symbol']

        x = row['pdbx_model_Cartn_x_ideal']
        y = row['pdbx_model_Cartn_y_ideal']
        z = row['pdbx_model_Cartn_z_ideal']

        if x and y and z:
            coordinates_ideal[atom] = [float(x),float(y),float(z)]

        chirality = row.get('pdbx_stereo_config',None)
        atom_to_chirality[atom] = chirality


    # not sure why but sometimes the ideal coordinates are more complete
    if len(coordinates_ideal) > len(coordinates):
        coordinates = coordinates_ideal

    # see if this fixes some H5' and H5'' labeling discrepancies
    # if base in ['5HC']:
    #     coordinates = coordinates_ideal

    for row in cif_data['chem_comp_bond']:
        atom1 = row['atom_id_1']
        atom2 = row['atom_id_2']
        connections.append((atom1,atom2))
        connections.append((atom2,atom1))

    connections = list(set(connections))

    data['atom_count'] = len(atom_to_element)

    return data, coordinates, connections, atom_to_element, atom_to_chirality


def read_csv(filename):
    with open(filename,"r", encoding="utf-8") as f:
        lines = f.readlines()
    header = lines[0].split('","')
    header[0] = header[0][1:]
    header[-1] = header[-1][:-2]
    all_d = []
    for line in lines[1:]:
        fields = line.split('","')
        fields[0] = fields[0][1:]
        fields[-1] = fields[-1][:-2]
        d = {}
        for i in range(len(header)):
            d[header[i]] = fields[i]
        all_d.append(d)
    return all_d

