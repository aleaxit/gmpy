#!/bin/bash

DEST_DIR=$1
WHEEL=$2
LD_LIBRARY_PATH="$(pwd)/.local/lib:$LD_LIBRARY_PATH"

if [[ "$OSTYPE" == "darwin"* ]]
then
  pip install -U git+https://github.com/skirpichev/delocate.git@fix-lib-sdir
  delocate-wheel --lib-sdir .libs -w ${DEST_DIR} -v ${WHEEL}
else
  auditwheel repair -w ${DEST_DIR} ${WHEEL}
fi
