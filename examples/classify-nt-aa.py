#!/usr/bin/env python

"""A script to print the classifications of RNA/AA interactions. The
classification for each pair will be printed out.
"""

import os
import sys
import argparse

here = os.path.abspath(os.path.dirname(__file__))
fr3d = os.path.abspath(os.path.join(here, ".."))
sys.path.insert(0, fr3d)

from fr3d.cif.reader import Cif as Reader
from fr3d.classifiers.rna_protein import Classifier


def main(filename):
    with open(filename, 'rb') as raw:
        structure = Reader(raw).structure()

    classifier = Classifier()
    classifications = classifier.classify(structure)
    for classification in classifications:
        print(classification)

        # 1MZP


if __name__ == '__main__':
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("cif", help="Cif file to read")
    args = parser.parse_args()
    main(args.cif)
