"""Read a mmCIF file and get all nucleotides near an amino acid.
"""

import os
import sys
import argparse
import collections as coll
import itertools as it

import numpy as np

# We import KDTree because this a way of storing points and quickly finding all
# nearby points. If we do multiple queries for distance it should be faster.
from scipy.spatial import cKDTree as KDTree

here = os.path.abspath(os.path.dirname(__file__))
fr3d = os.path.abspath(os.path.join(here, ".."))
sys.path.insert(0, fr3d)

from fr3d.cif.reader import Cif


def generate_peptide_tree(structure, limits={}):
    """Generate a KDTree for all peptides in a structure.
    """
    aa = structure.residues(type=['L-peptide linking', 'PEPTIDE LINKING'],
                            **limits)
    data = []
    mapping = []
    for residue in aa:
        for atom in residue.atoms():
            data.append(atom.coordinates())
            mapping.append(residue)
    return mapping, KDTree(np.array(data))


def generate_rna_tree(structure, limits={}):
    """Generate a KDTree for all residues in a structure.
    """
    data = []
    mapping = []
    for residue in structure.residues(type='RNA linking', **limits):
        for atom in residue.atoms():
            data.append(atom.coordinates())
            mapping.append(residue)
    return mapping, KDTree(np.array(data))


def main(filename, limits):
    with open(filename, 'rb') as raw:
        structure = Cif(raw).structure()

    aa_mapping, aa_tree = generate_peptide_tree(structure,
                                                limits.get('residue', {}))
    rna_mapping, rna_tree = generate_rna_tree(structure, limits.get('nt', {}))

    # This produces a list of (rna_atom_index, aa_atom_index)
    # The query_ball_tree method goes through all points in rna_tree and finds
    # all points in aa_tree that have distance 5. It will return something for
    # all rna atoms, even those without any nearby aa atoms.
    results = enumerate(rna_tree.query_ball_tree(aa_tree, 5))

    # Filter out all rna atoms that have nothing near by
    filtered = it.ifilter(lambda (r, a): a, results)

    # Map the results so we transform the rna atom to a residue
    mapped = it.imap(lambda (r, _): (rna_mapping[r], _), filtered)

    # Group the list by the residue
    grouped = it.groupby(mapped, lambda (r, _): r)

    # Display what we find
    for residue, nearby in grouped:

        # This creates a dictionary we use to count the number of near by
        # atoms.
        counts = coll.defaultdict(int)

        # Near by is a iterable of the form (residue, [indexes])
        for (_, indexes) in nearby:
            for index in indexes:
                peptide = aa_mapping[index]
                counts[peptide] += 1

        for peptide, count in counts.items():
            print('%s: %s %s' % (residue.unit_id(), peptide.unit_id(), count))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--residue", default='', dest='residue',
                        help="Type of residue to limit to")
    parser.add_argument("--nucleotide", default='', dest='nt',
                        help="Type of nucleotide to limit to")
    parser.add_argument("cif", help="CIF file to read")
    args = parser.parse_args()
    limits = {}
    if args.nt:
        limits['nt'] = {'sequence': args.nt}
    if args.residue:
        limits['residue'] = {'sequence': args.residue}
    main(args.cif, limits)
