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


# Open MPI configuration (mpifort)

set(MPI_DEFINES "-DOPEN_MPI")

if(CMAKE_BUILD_TYPE STREQUAL "Test")
    #do nothing at all
else()
    set(MPI_BASE_FLAGS
        -cpp
        -Wall
        -Wextra
        -Wno-c-binding-type
    )
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(MPI_OPT_FLAGS -O0 -g -fbounds-check -fbacktrace)
else()
    set(MPI_OPT_FLAGS -O3)
endif()

set(MPI_COMPILE_FLAGS ${MPI_BASE_FLAGS} ${MPI_OPT_FLAGS})
