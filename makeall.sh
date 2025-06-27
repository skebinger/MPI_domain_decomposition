#!/bin/bash

# Usage:
# ./makeall.sh [BuildType] [SharedLib ON/OFF] [MPIFLAVOR]
# Example:
# ./makeall.sh Debug OFF openmpi
# ./makeall.sh Release ON intelmpi

BUILD_TYPE=${1:-Release}
BUILD_SHARED=${2:-OFF}
MPIFLAVOR=${3:-intelmpi}

BUILD_DIR=build

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
    -DMPIFLAVOR=$MPIFLAVOR

make -j$(nproc)

echo "=== Build Finished ==="
