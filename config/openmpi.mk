FC := mpif90
MPI_FLAGS      := -DOMPI
COMPILER_FLAGS := -fPIC

ifeq ($(BUILD_TYPE),debug)
    BUILD_FLAGS := -g -Wall -fcheck=all -traceback -O0
else
    BUILD_FLAGS := -O2
endif

FCFLAGS := -J$(MODDIR) -I$(MODDIR) $(MPI_FLAGS) $(COMPILER_FLAGS) $(BUILD_FLAGS)
