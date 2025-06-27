# Intel MPI configuration (mpiifx)

set(MPI_DEFINES "-DINTEL_MPI")

set(MPI_BASE_FLAGS
    -cpp
    -fpp
    -warn all
    -traceback
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(MPI_OPT_FLAGS -O0 -g -check all -debug extended)
else()
    set(MPI_OPT_FLAGS -O3 -xHost)
endif()

set(MPI_COMPILE_FLAGS ${MPI_BASE_FLAGS} ${MPI_OPT_FLAGS})
