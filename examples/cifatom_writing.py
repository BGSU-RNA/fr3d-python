import os
import sys
import argparse

here = os.path.abspath(os.path.dirname(__file__))
fr3d = os.path.abspath(os.path.join(here, ".."))
sys.path.insert(0, fr3d)

from fr3d.cif.reader import Cif as Reader
from fr3d.cif.writer import CifAtom as Writer


def main(filename):
    with open(filename, 'rb') as raw:
        structure = Reader(raw).structure()

    savename = filename + 'atoms'
    with open(savename, 'wb') as out:
        writer = Writer(out)
        writer(structure)


if __name__ == '__main__':
    parser = argparse.ArgumentParser("Create a cif atom file")
    parser.add_argument("cif", help="Cif file to read")
    args = parser.parse_args()
    main(args.cif)
