#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

wget "http://mmcif.wwpdb.org/docs/sw-examples/python/src/pdbx.tar.gz"
mkdir pdbx
mv pdbx.tar.gz pdbx/
pushd pdbx
tar xvf pdbx.tar.gz
CAT<<EOS
from distutils.core import setup, , find_packages

setup(
    name='pdbx',
    version='1.0.5',
    author='NDB',
    author_email='',
    packages=find_packages(include=['pdbx', 'pdbx.*']),
    url='',
    license='LICENSE.txt',
    description='Pure python tools to work with cif data from PDB',
    long_description='Pure python tools to work with cif data from PDB',
)
EOS

python setup.py install
popd
