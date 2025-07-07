# ==============*CMAKE*============= #
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

cmake_minimum_required(VERSION 3.16)

# User configurable option to choose MPI flavor
set(MPIFLAVOR "intelmpi" CACHE STRING "MPI flavor: intelmpi or openmpi")

# Pick Fortran compiler early, before project()
if(MPIFLAVOR STREQUAL "intelmpi")
    # Intel MPI uses mpiifx
    set(CMAKE_Fortran_COMPILER "mpiifx" CACHE STRING "Intel MPI Fortran compiler" FORCE)
    set(CMAKE_C_COMPILER       "mpiicx" CACHE STRING "Intel MPI C compiler" FORCE) 
elseif(MPIFLAVOR STREQUAL "openmpi")
    # Open MPI uses mpif90 or mpifort
    set(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "Open MPI Fortran compiler" FORCE)
    set(CMAKE_C_COMPILER       "mpicc" CACHE STRING "Open MPI C compiler" FORCE)
else()
    message(FATAL_ERROR "Unsupported MPIFLAVOR: ${MPIFLAVOR}. Supported: intelmpi, openmpi")
endif()

project(mpidcl LANGUAGES Fortran C)

# Include MPI flavor specific flags and defines
include(config/${MPIFLAVOR}.cmake)

# Where to put module files
set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/mod)

# Source files
file(GLOB SRC_FILES src/*.f90)

add_library(${PROJECT_NAME} ${SRC_FILES})

# Apply compiler flags and defines from flavor file
target_compile_options(${PROJECT_NAME} PRIVATE ${MPI_COMPILE_FLAGS})
target_compile_definitions(${PROJECT_NAME} PRIVATE ${MPI_DEFINES})

# Optional: module directory include
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

# Output settings
set_target_properties(${PROJECT_NAME} PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

# Optional debug/release settings (can also be handled in flavor files)
