# ========================
# Configurable Parameters
# ========================

BUILD_TYPE ?= release
# 'release' or 'debug'
LIBTYPE    ?= static
# 'static' or 'shared'
TARGET     ?= mpidcl
# output library name (no extension)
MPIFLAVOR  ?= intelmpi
# 'openmpi', 'intelmpi'

# ========================
# Directories
# ========================

SRCDIR     := src
MODDIR     := mod
OBJDIR     := build
CONFIGDIR  := config

SRC        := $(wildcard $(SRCDIR)/*.f90)
OBJ        := $(patsubst $(SRCDIR)/%.f90, $(OBJDIR)/%.o, $(SRC))

# ========================
# Load Compiler Config
# ========================

CONFIG_FILE := $(strip $(CONFIGDIR)/$(MPIFLAVOR).mk)

ifeq ($(wildcard $(CONFIG_FILE)),)
$(error "Invalid MPIFLAVOR: '$(MPIFLAVOR)'. No such config file: $(CONFIG_FILE)")
endif

include $(CONFIG_FILE)

# ========================
# Library Output Settings
# ========================

ifeq ($(LIBTYPE),shared)
    LIBEXT  := so
    LINKCMD := $(FC) -shared $(LDFLAGS) -o $(TARGET).$(LIBEXT) $(OBJ)
else
    LIBEXT  := a
    LINKCMD := ar rcs $(TARGET).$(LIBEXT) $(OBJ)
endif

# ========================
# Build Rules
# ========================

all: $(MODDIR) $(OBJDIR) $(TARGET).$(LIBEXT)

$(MODDIR) $(OBJDIR):
	mkdir -p $@

$(OBJDIR)/%.o: $(SRCDIR)/%.f90
	$(FC) $(FCFLAGS) -c $< -o $@

$(TARGET).$(LIBEXT): $(OBJ)
	$(LINKCMD)

clean:
	rm -rf $(OBJDIR) $(MODDIR) $(TARGET).a $(TARGET).so

new: clean all

.PHONY: all clean
