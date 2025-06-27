# Open MPI configuration (mpifort)

set(MPI_DEFINES "-DOPEN_MPI")

set(MPI_BASE_FLAGS
    -cpp
    -Wall
    -Wextra
)

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(MPI_OPT_FLAGS -O0 -g -fbounds-check -fbacktrace)
else()
    set(MPI_OPT_FLAGS -O3)
endif()

set(MPI_COMPILE_FLAGS ${MPI_BASE_FLAGS} ${MPI_OPT_FLAGS})
