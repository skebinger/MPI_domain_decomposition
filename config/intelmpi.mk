FC := mpiifx
MPI_FLAGS      := -DINTEL_MPI
COMPILER_FLAGS := 

ifeq ($(BUILD_TYPE),debug)
    BUILD_FLAGS := -module mod -c -cpp -O0
else
    BUILD_FLAGS := 
endif

FCFLAGS := -module $(MODDIR) $(MPI_FLAGS) $(COMPILER_FLAGS) $(BUILD_FLAGS)
