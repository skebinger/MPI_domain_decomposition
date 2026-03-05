#!/bin/bash

# ==============*BASH*============== #
#   _ __ ___  _ __ (_) __| | ___| |  #
#  | '_ ` _ \| '_ \| |/ _` |/ __| |  #
#  | | | | | | |_) | | (_| | (__| |  #
#  |_| |_| |_| .__/|_|\__,_|\___|_|  #
#            |_|                     #
# ================================== #
# ================================================================= #
#  Copyright (C) 2025, Simon Kebinger
# 
#  This file is part of the MPI decomposition library "mpidcl" for 
#  structured multidmensional domains.
# 
#  This library is distributed under the BSD 3-Clause License.
# ================================================================= #

set -e
set -o pipefail

# -------------------------
# PFUnit version to install
# -------------------------
PFUNIT_TAG="v4.12.0"

# -------------------------
# Directories
# -------------------------
SRC_DIR="$HOME/opt/pfunit-src"
INSTALL_OPENMPI="$HOME/opt/pfunit-openmpi"
INSTALL_INTELMPI="$HOME/opt/pfunit-intelmpi"
BUILD_DIR="$SRC_DIR/build"

# -------------------------
# Prompt for MPI flavor
# -------------------------
echo "Which MPI version of pFUnit would you like to build?"
echo "1) OpenMPI (mpifort)"
echo "2) Intel MPI (mpiifx)"
echo "3) Both"
read -rp "Enter choice [1/2/3]: " choice

# -------------------------
# Clone PFUnit if needed
# -------------------------
if [ ! -d "$SRC_DIR" ]; then
    echo "Cloning pFUnit ${PFUNIT_TAG} into $SRC_DIR..."
    git clone --branch "$PFUNIT_TAG" --depth 1 https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git "$SRC_DIR"
fi

# -------------------------
# Build function
# -------------------------
build_pfunit() {
    local FC_COMP=$1
    local CC_COMP=$2
    local PREFIX=$3
    local LABEL=$4

    echo "Building pFUnit with $LABEL..."

    rm -rf "$BUILD_DIR"
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"

    export FC=$FC_COMP
    export CC=$CC_COMP

    cmake "$SRC_DIR" \
        -DCMAKE_INSTALL_PREFIX="$PREFIX" \
        -DMPI=YES \
        -DENABLE_MPI_F08=YES \
        -DOPENMP=NO \
        -Dcoverage=NO

    make -j$(nproc)
    make install

    echo "Installed $LABEL version to $PREFIX"
}

# -------------------------
# Build OpenMPI version
# -------------------------
if [[ "$choice" == "1" || "$choice" == "3" ]]; then
    build_pfunit "mpifort" "mpicc" "$INSTALL_OPENMPI" "OpenMPI"
fi

# -------------------------
# Build Intel MPI version
# -------------------------
if [[ "$choice" == "2" || "$choice" == "3" ]]; then
    build_pfunit "mpiifx" "mpiicx" "$INSTALL_INTELMPI" "Intel MPI"
fi

echo "pFUnit build(s) complete."