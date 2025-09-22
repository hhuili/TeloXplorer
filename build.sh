#!/bin/bash

set -ex

cd ${SRC_DIR}/src

make

cp BLOOM ${PREFIX}/bin/

cd ${SRC_DIR}

${PYTHON} -m pip install . --no-deps -vv --use-pep517

