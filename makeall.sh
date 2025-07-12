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

# Usage:
# ./makeall.sh [BuildType] [SharedLib ON/OFF] [MPIFLAVOR] [BuildLib ON/OFF] [BuildTest ON/OFF]
# BuildType : Test, Debug, Release (=intended build type for use)
# SharedLib : decides if a static or shared library is created on output
# So far: Testing is restricted to a static library and compilation might fail for shared
# MPIFLAVOUR : intelmpi, openmpi or a mpi implementation of your choice given you change and
# ammend the cmake extensions in ./config
# BuildLib : On current pass build the lib or not
# BuildTest : On current pass build the testing environment; requires library to be build previously
# -> envoke tests (if compiled) in ./build_test with command <ctest>
# Example:
# ./makeall.sh Debug OFF openmpi ON OFF -> to compile without tests but debugging flags with intelmpi
# ./makeall.sh Test OFF intelmpi ON ON -> to compile with testing (requires pFUnit) with intelmpi
# ./makeall.sh Release OFF openmpi ON OFF -> to compile with testing (requires pFUnit) with openmpi

BUILD_TYPE=${1:-Debug}
BUILD_SHARED=${2:-OFF}
MPIFLAVOR=${3:-intelmpi}
BUILD_LIB=${4:-OFF}
BUILD_TESTS=${5:-OFF}

BUILD_DIR=build
TEST_BUILD_DIR=build_test

if [ "$BUILD_LIB" = "ON" ]; then

    echo "=== Building Project ==="
    echo "Build Type:        $BUILD_TYPE"
    echo "Shared Library:    $BUILD_SHARED"
    echo "MPI Flavor:        $MPIFLAVOR"
    echo "Build Directory:   $BUILD_DIR"
    echo "========================"

    # Clean build dir for a fresh build
    rm -rf $BUILD_DIR
    mkdir -p $BUILD_DIR
    cd $BUILD_DIR

    cmake .. \
        -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
        -DBUILD_SHARED_LIBS=$BUILD_SHARED \
        -DMPIFLAVOR=$MPIFLAVOR \
        #--debug-output

    make -j$(nproc)

    cd ..
fi

if [ "$BUILD_TESTS" = "ON" ]; then
    echo "=== Building Tests ==="
    echo "MPI Flavor:        $MPIFLAVOR"
    echo "Build Directory:   $TEST_BUILD_DIR"
    echo "====================="

    rm -rf $TEST_BUILD_DIR
    mkdir -p $TEST_BUILD_DIR
    cd $TEST_BUILD_DIR

    echo $(pwd)

    cmake ../test \
        -DMPIFLAVOR=$MPIFLAVOR \
        #--debug-output

    make -j$(nproc)

    cd ..
fi

echo "=== Build Finished ==="
