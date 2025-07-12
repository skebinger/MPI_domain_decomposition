# Compiler and flags
FC       = mpifort
FFLAGS   = -std=f2008 -O2 -ffree-line-length-none

# Directories
SRC_DIR   = ../src
BUILD_DIR = ../build_test_output
EXEC      = $(BUILD_DIR)/test_driver_output

# Source files
MODULE_SRC  = $(wildcard $(SRC_DIR)/*.f90)
DRIVER_SRC  = test_driver_output.f90

# Object files (built in BUILD_DIR)
MODULE_OBJS = $(patsubst $(SRC_DIR)/%.f90, $(BUILD_DIR)/%.o, $(MODULE_SRC))
DRIVER_OBJ  = $(BUILD_DIR)/test_driver_output.o

# Default target
all: $(EXEC)

# Link final executable
$(EXEC): $(MODULE_OBJS) $(DRIVER_OBJ)
	$(FC) $(FFLAGS) -o $@ $^

# Compile module sources
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILD_DIR)

# Compile test driver at root
$(BUILD_DIR)/test_driver_output.o: test_driver_output.f90 | $(BUILD_DIR)
	$(FC) $(FFLAGS) -c $< -o $@ -J$(BUILD_DIR) -I$(BUILD_DIR)

# Ensure build directory exists
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Cleanup
clean:
	rm -rf $(BUILD_DIR)

# mpirun for testing
test: all
	mpirun -np 8 ../build_test_output/test_driver_output

.PHONY: all clean
