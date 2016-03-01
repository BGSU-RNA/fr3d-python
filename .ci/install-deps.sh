#!/usr/bin/env bash

set -euo pipefail
IFS=$'\n\t'

wget "http://mmcif.wwpdb.org/docs/sw-examples/python/src/pdbx.tar.gz"
mkdir pdbx
mv pdbx.tar.gz pdbx/
pushd pdbx
tar xvf pdbx.tar.gz
cp "../.ci/setup.py" .
python setup.py install
popd
rm -r pdbx
