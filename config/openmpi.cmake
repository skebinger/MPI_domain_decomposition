# Open MPI configuration (mpifort)

set(MPI_DEFINES "-DOPEN_MPI")

if(CMAKE_BUILD_TYPE STREQUAL "Test")
    #do nothing at all
else()
    set(MPI_BASE_FLAGS
        -cpp
        -Wall
        -Wextra
    )
endif()

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(MPI_OPT_FLAGS -O0 -g -fbounds-check -fbacktrace)
else()
    set(MPI_OPT_FLAGS -O3)
endif()

set(MPI_COMPILE_FLAGS ${MPI_BASE_FLAGS} ${MPI_OPT_FLAGS})
