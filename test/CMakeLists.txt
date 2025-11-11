cmake_minimum_required(VERSION 3.16)

# Use the same MPI flavor option as parent or define again
set(MPIFLAVOR "openmpi" CACHE STRING "MPI flavor: intelmpi or openmpi")

# Pick Fortran compiler early, before project()
if(MPIFLAVOR STREQUAL "intelmpi")
    set(CMAKE_Fortran_COMPILER "mpiifx" CACHE STRING "Intel MPI Fortran compiler" FORCE)
    set(CMAKE_C_COMPILER "mpiicx" CACHE STRING "Intel MPI C compiler" FORCE)
elseif(MPIFLAVOR STREQUAL "openmpi")
    set(CMAKE_Fortran_COMPILER "mpifort" CACHE STRING "Open MPI Fortran compiler" FORCE)
    set(CMAKE_C_COMPILER "mpicc" CACHE STRING "Open MPI C compiler" FORCE)
else()
    message(FATAL_ERROR "Unsupported MPIFLAVOR: ${MPIFLAVOR}. Supported: intelmpi, openmpi")
endif()

project(mpidcl_tests LANGUAGES Fortran C)

# get the funitproc executable with full path
set(FUNITPROC_EXECUTABLE "$ENV{HOME}/pfunit-${MPIFLAVOR}/PFUNIT-4.12/bin/funitproc")

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/mod_tests")

# Find the MPI and PFUNIT packages
find_package(MPI REQUIRED Fortran C)
list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/pfunit-${MPIFLAVOR}/PFUNIT-4.12/cmake")
find_package(PFUNIT REQUIRED)

set(PFUNIT_INCLUDE_DIR "$ENV{HOME}/pfunit-${MPIFLAVOR}/PFUNIT-4.12/include")
set(PFUNIT_DRIVER_F90 "${PFUNIT_INCLUDE_DIR}/driver.F90")

# Directory where your .pf files live
set(TEST_PF_DIR "${CMAKE_SOURCE_DIR}/../test")

# Gather .pf files
file(GLOB PF_TEST_FILES "${TEST_PF_DIR}/*.pf")

# Output directory for generated Fortran files
set(GENERATED_SRC_DIR "${CMAKE_BINARY_DIR}/generated")
#set(MOD_DIR ${CMAKE_Fortran_MODULE_DIRECTORY})

file(MAKE_DIRECTORY ${GENERATED_SRC_DIR})

# Convert .pf -> .F90 using funitproc
set(GENERATED_F90_FILES "")
foreach(pf_file IN LISTS PF_TEST_FILES)
    get_filename_component(pf_name ${pf_file} NAME_WE)
    set(f90_file "${GENERATED_SRC_DIR}/${pf_name}.F90")

    add_custom_command(
        OUTPUT ${f90_file}
        COMMAND ${FUNITPROC_EXECUTABLE} ${pf_file} ${f90_file}
        DEPENDS ${pf_file}
        COMMENT "Running funitproc on ${pf_file}"
        VERBATIM
    )

    list(APPEND GENERATED_F90_FILES ${f90_file})
endforeach()

# Add the driver.F90 to sources
list(APPEND GENERATED_F90_FILES ${PFUNIT_DRIVER_F90})

add_custom_target(generate_pf_files ALL DEPENDS ${GENERATED_F90_FILES})

# Create the test executable from generated sources + driver
add_executable(pfunit_tests ${GENERATED_F90_FILES})
add_dependencies(pfunit_tests generate_pf_files)

target_include_directories(pfunit_tests PRIVATE
    ${CMAKE_Fortran_MODULE_DIRECTORY}
    "$ENV{HOME}/pfunit-${MPIFLAVOR}/FARGPARSE-1.9/include"
    "$ENV{HOME}/pfunit-${MPIFLAVOR}/GFTL_SHARED-1.10/include"
    ${PFUNIT_INCLUDE_DIR}
    "${CMAKE_SOURCE_DIR}/../build/mod"
)

target_link_libraries(pfunit_tests PRIVATE
    "$ENV{HOME}/pfunit-${MPIFLAVOR}/PFUNIT-4.12/lib/libpfunit.a"
    "$ENV{HOME}/pfunit-${MPIFLAVOR}/PFUNIT-4.12/lib/libfunit.a"
    "$ENV{HOME}/pfunit-${MPIFLAVOR}/FARGPARSE-1.9/lib/libfargparse.a"
    "$ENV{HOME}/pfunit-${MPIFLAVOR}/GFTL_SHARED-1.10/lib/libgftl-shared-v2-as-default.a"
    "$ENV{HOME}/pfunit-${MPIFLAVOR}/GFTL_SHARED-1.10/lib/libgftl-shared-v2.a"
    "$ENV{HOME}/pfunit-${MPIFLAVOR}/GFTL_SHARED-1.10/lib/libgftl-shared-v1.a"
    MPI::MPI_Fortran MPI::MPI_C mpidcl
    OpenMP::OpenMP_Fortran
)

target_link_directories(pfunit_tests PRIVATE
    "${CMAKE_SOURCE_DIR}/../build/lib"
)
target_compile_definitions(pfunit_tests PRIVATE
    _TEST_SUITES=\"${CMAKE_SOURCE_DIR}/test_suites.inc\"
)

target_compile_options(pfunit_tests PRIVATE ${MPI_COMPILE_FLAGS})
target_compile_definitions(pfunit_tests PRIVATE ${MPI_DEFINES})

set_target_properties(pfunit_tests PROPERTIES
    Fortran_MODULE_DIRECTORY "${CMAKE_Fortran_MODULE_DIRECTORY}"
)

enable_testing()
add_test(NAME pfunit_tests COMMAND mpirun -np 8 $<TARGET_FILE:pfunit_tests>)
