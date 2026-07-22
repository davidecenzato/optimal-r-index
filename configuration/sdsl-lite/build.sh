#!/bin/bash
CUR_DIR=`pwd`
OLD_DIR="$( cd "$( dirname "$0" )" && pwd )" # gets the directory where the script is located in
cd "${OLD_DIR}"

cmake .. -DCMAKE_POLICY_VERSION_MINIMUM=3.12 -DCMAKE_INSTALL_PREFIX=${HOME} && make install

cd ${CUR_DIR}
