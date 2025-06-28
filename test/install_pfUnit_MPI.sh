#!/bin/bash

set -e
set -o pipefail

# 📁 Configuration
SRC_DIR="$HOME/pfunit-src"
INSTALL_OPENMPI="$HOME/pfunit-openmpi"
INSTALL_INTELMPI="$HOME/pfunit-intelmpi"
BUILD_DIR="$SRC_DIR/build"

# --- Auto-setup OpenMPI environment if mpifort not found ---
if ! command -v mpifort >/dev/null 2>&1; then
  echo "mpifort not found in PATH, attempting to add OpenMPI paths..."

  # Common Fedora OpenMPI locations:
  OPENMPI_BIN="/usr/lib64/openmpi/bin"
  OPENMPI_LIB="/usr/lib64/openmpi/lib"

  if [ -d "$OPENMPI_BIN" ] && [ -d "$OPENMPI_LIB" ]; then
    export PATH="$OPENMPI_BIN:$PATH"
    export LD_LIBRARY_PATH="$OPENMPI_LIB:$LD_LIBRARY_PATH"
    echo "Added OpenMPI bin and lib paths to environment."

    # Verify again
    if ! command -v mpifort >/dev/null 2>&1; then
      echo "Still cannot find mpifort after adding paths. Please install OpenMPI or adjust paths."
      exit 1
    fi

  else
    echo "Could not find OpenMPI directories at $OPENMPI_BIN or $OPENMPI_LIB."
    echo "Please install OpenMPI or manually set PATH and LD_LIBRARY_PATH."
    exit 1
  fi
fi

# ⛳ Prompt user for choice
echo "Which version(s) of pFUnit would you like to build?"
echo "1) OpenMPI (mpifort)"
echo "2) Intel MPI (mpiifx)"
echo "3) Both"
read -rp "Enter choice [1/2/3]: " choice

# 📥 Clone if needed
if [ ! -d "$SRC_DIR" ]; then
    echo "Cloning pFUnit into $SRC_DIR..."
    git clone https://github.com/Goddard-Fortran-Ecosystem/pFUnit.git "$SRC_DIR"
fi

# 🔁 Function to build with given compiler
build_pfunit() {
    local FC_COMP=$1
    local CC_COMP=$2
    local PREFIX=$3
    local LABEL=$4

    echo "==> Building pFUnit with $LABEL..."

    rm -rf "$BUILD_DIR"
    mkdir -p "$BUILD_DIR"
    cd "$BUILD_DIR"

    export FC=$FC_COMP
    export CC=$CC_COMP

    cmake "$SRC_DIR" \
        -DCMAKE_INSTALL_PREFIX="$PREFIX" \
        -DMPI=YES \
        -DOPENMP=NO \
        -Dcoverage=NO

    make -j$(nproc)
    make install

    echo "✅ Installed $LABEL version to $PREFIX"
}

# 🧱 Build OpenMPI version
if [[ "$choice" == "1" || "$choice" == "3" ]]; then
    build_pfunit "mpifort" "mpicc" "$INSTALL_OPENMPI" "OpenMPI"
fi

# 🧱 Build Intel MPI version
if [[ "$choice" == "2" || "$choice" == "3" ]]; then
    echo "Loading Intel MPI environment..."
    # Adjust path to your actual Intel OneAPI setvars.sh if different
    source /opt/intel/oneapi/setvars.sh || {
        echo "❌ Failed to source Intel setvars.sh — please check Intel MPI installation."
        exit 1
    }

    build_pfunit "mpiifx" "mpiicx" "$INSTALL_INTELMPI" "Intel MPI"
fi

echo "🎉 Done. pFUnit build(s) complete."
