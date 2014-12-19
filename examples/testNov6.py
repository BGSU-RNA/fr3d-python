"""Iterate over parts of a cif file."""



import argparse

from fr3d.cif.reader import Cif


def main(filename):
    with open(filename, 'rb') as raw:
        structure = Cif(raw).structure()

    print('Iterating over all parts')
    for residue in structure.residues(chain='A', sequence = ['A','U','G','C']):
        print "residue center ", residue.centers['base']
            

if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("cif", help="mmCIF file to read")
    args = parser.parse_args()
    main(args.cif)
