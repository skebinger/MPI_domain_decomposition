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
# ./makeall.sh [options]
#
# Options (named arguments):
#   --type <BuildType>         : Build type (Test, Debug, Release). Default: Release
#   --shared <ON|OFF>          : Build shared library (ON) or static library (OFF). Default: OFF
#   --mpi <MPIFLAVOR>          : MPI implementation to use (intelmpi, openmpi). Default: openmpi
#   --build-lib <ON|OFF>       : Whether to build the library. Default: ON
#   --build-tests <ON|OFF>     : Whether to build the test environment (requires library built). Default: OFF
#   --install <InstallPath>    : Custom installation directory for library and module files. Default: ./install
#   --help                     : Show this help message and exit
#
# -> envoke tests (if compiled) in ./build_test with command <ctest>
#
# Examples:
#   ./makeall.sh --type Debug
#   ./makeall.sh --mpi intelmpi --shared ON
#   ./makeall.sh --install $(HOME)/lib/mpidcl
#   ./makeall.sh --mpi intelmpi --install $(HOME)/lib/mpidcl
# DEFAULT: the script executes as:
# ./makeall.sh --type Release --shared OFF --mpi openmpi --build-lib ON --build-tests OFF --install $(pwd)/install

# ===============================
# Default values
# ===============================

BUILD_TYPE="Release"
BUILD_SHARED="OFF"
MPIFLAVOR="openmpi"
BUILD_LIB="ON"
BUILD_TESTS="OFF"
INSTALL_PREFIX="$(pwd)/install"

# ===============================
# Argument parsing
# ===============================

while [[ $# -gt 0 ]]; do
    case "$1" in
        --type)
            BUILD_TYPE="$2"
            shift 2
            ;;
        --shared)
            BUILD_SHARED="$2"
            shift 2
            ;;
        --mpi)
            MPIFLAVOR="$2"
            shift 2
            ;;
        --build-lib)
            BUILD_LIB="$2"
            shift 2
            ;;
        --build-tests)
            BUILD_TESTS="$2"
            shift 2
            ;;
        --install)
            INSTALL_PREFIX="$2"
            shift 2
            ;;
        --help)
            echo "Usage:"
            echo "./makeall.sh [options]"
            echo ""
            echo "Options:"
            echo "  --type <Release|Debug|Test>"
            echo "  --shared <ON|OFF>"
            echo "  --mpi <intelmpi|openmpi>"
            echo "  --build-lib <ON|OFF>"
            echo "  --build-tests <ON|OFF>"
            echo "  --install <path>"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

BUILD_DIR=build
TEST_BUILD_DIR=build_test

if [ "$BUILD_LIB" = "ON" ]; then

    echo "=== Building Project ==="
    echo "Build Type:        $BUILD_TYPE"
    echo "Shared Library:    $BUILD_SHARED"
    echo "MPI Flavor:        $MPIFLAVOR"
    echo "Install Prefix:    $INSTALL_PREFIX"
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
        -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX

    make -j$(nproc)

    echo "=== Installing Library ==="
    cmake --install .

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
        -DMPIFLAVOR=$MPIFLAVOR

    make -j$(nproc)

    cd ..
fi

echo "=== Build Finished ==="
