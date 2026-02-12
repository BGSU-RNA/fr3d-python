"""
Test r3dcid.py under a variety of inputs
"""

import random
import requests

from r3dcid import main as r3dcid
from r3dcid import arc_group_to_interactions

# examples for testing different sizes of RNA molecule

chain_strings = set()

if True:
    chain_strings.add('4QQB|1|P')
    chain_strings.add('1S72|1|9+1S72|1|0')
    chain_strings.add('3GS5|1|C')
    chain_strings.add('5VSU|1|I')
    chain_strings.add('4V8P|1|D1+4V8P|1|C2+4V8P|1|A1+4V8P|1|B2+4V8P|1|F1+4V8P|1|E2+4V8P|1|H1+4V8P|1|G2')
    chain_strings.add('4QLM|1|A')
    chain_strings.add('5T5H|1|E')
    chain_strings.add('3JCS|1|3')
    chain_strings.add('6SWD|1|2')
    chain_strings.add('6ZU5|1|S60')
    chain_strings.add('6IR9|1|P')
    chain_strings.add('4Y4O|1|2A')
    chain_strings.add('6ERI|1|AA')
    chain_strings.add('7ASE|1|0')
    chain_strings.add('4V88|1|A6')
    chain_strings.add('6ID1|1|H+6ID1|1|F')
    chain_strings.add('7O7Y|1|B5+7O7Y|1|B8')
    chain_strings.add('6UZ7|1|5+6UZ7|1|8')
    chain_strings.add('3J7Q|1|5+3J7Q|1|8')
    chain_strings.add('4V9F|1|0+4V9F|1|9')
    chain_strings.add('1R3E|1|C')
    chain_strings.add('7MLW|1|F')
    chain_strings.add('XXXX|1|B')   # made up name to test that it fails gracefully
    chain_strings.add('6ZMI|1|L7+6ZMI|1|L8+6ZMI|1|L5+6ZMI|1|S2+6ZMI|1|CC')  # Homo sapiens 5S, 5.8S, LSU, SSU, tRNA
    chain_strings.add('3CZW|1|X')
    chain_strings.add('1J5E|1|A')
    chain_strings.add('7JQQ|1|K+7JQQ|1|L+7JQQ|1|M+7JQQ|1|N+7JQQ|1|O') # Symmetrical
    chain_strings.add('5J7L|1|DA,5J7L|1|AA,5J7L|1|DB')
    chain_strings.add('4RKV|1|A,4RKV|1|B')
    chain_strings.add('7LHD|1|A')
    chain_strings.add('4V9F')
    chain_strings.add('4V9F|1|9')
    chain_strings.add("4V9O")
    chain_strings.add("6ZJ3|1|LA+6ZJ3|1|LB+6ZJ3|1|LC+6ZJ3|1|LD+6ZJ3|1|LE+6ZJ3|1|LF+6ZJ3|1|LG+6ZJ3|1|LH+6ZJ3|1|LI+6ZJ3|1|LJ+6ZJ3|1|LK+6ZJ3|1|LL+6ZJ3|1|LM+6ZJ3|1|LN+6ZJ3|1|LO")
    chain_strings.add("6SKG|1|BB")
    chain_strings.add('7K00') # E. coli 5S, LSU, SSU, mRNA, A-site, P-site, E-site
    chain_strings.add('4TNA')
    chain_strings.add("1R3E")
    chain_strings.add("6YDP") # very long sequence, not so many resolved

if False:
    urls = ["https://rna.bgsu.edu/rna3dhub/nrlist/download/dna/0.5/all/csv", "https://rna.bgsu.edu/rna3dhub/nrlist/download/rna/3.392/all/csv"]
    for url in urls:
        response = requests.get(url)
        for line in response.text.split("\n"):
            fields = line.split('","')
            if len(fields) > 2:
                pdb_id = fields[1].split("|")[0]
                chain_strings.add(pdb_id)

random.shuffle(sorted(chain_strings))
print('Processing %d PDB files' % len(chain_strings))

# loop over examples
for chain_string in chain_strings:

    params = {}

    params['output_path'] = "r3dcid_output"

    # sometimes choose display options randomly
    if random.random() < 0.1:
        # choose coloring randomly
        params['coloring'] = random.choice(["wong","grayscale"])

    if random.random() < 0.1:
        # select one or more from arc_group_to_interactions.keys()
        hide_set = []
        for a in arc_group_to_interactions.keys():
            if random.random() < 0.3:
                hide_set.append(a)
        params['hide'] = ",".join(hide_set)

    if random.random() < 0.1:
        # select one or more from arc_group_to_interactions.keys()
        dim_set = []
        for a in arc_group_to_interactions.keys():
            if random.random() < 0.3:
                dim_set.append(a)
        params['dim'] = ",".join(dim_set)

    if random.random() < 0.1:
        # select one or more from arc_group_to_interactions.keys()
        text_set = []
        for a in ["basepair","stacking","bph","br","sr","so","near","all","helix"]:
            if random.random() < 0.3:
                text_set.append(a)
        params['text'] = ",".join(text_set)

    if random.random() < 0.1:
        params['n3d'] = False

    if random.random() < 0.1:
        text_set = []
        for a in ["title","method","release_date","source","resolution","none","all"]:
            if random.random() < 0.2:
                text_set.append(a)
        if len(text_set) > 0:
            params['header'] = ",".join(text_set)

    params['format'] = "pdf,svg"

    output = r3dcid(chain_string, params)

    if len(output) == 0:
        input("Empty filename, press Enter to continue")

    print("")