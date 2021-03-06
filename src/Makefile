#
#     FILE: NAG workshop makefile
#   AUTHOR: Paul Bartholomew
# COMMENTS: Currently only gfortran is supported
#

TARGET_EXE = fd1d_heat_explicit
TARGET_LIB = fmw

FC = gfortran
FFLAGS = -I.

# Compiler flags to extract dependencies
ifeq ($(FC), gfortran)
	DEPFLAGS = -M -MD -MG -cpp -MF $@.tmp $<
else ifeq($(FC), nagfor)
	DEPFLAGS = =depend -otype=make $< -o $@.tmp
endif

# Setup sources and objects variables
SOURCES := types_mod.f90 rhs_mod.f90 cfl_mod.f90 io_mod.f90 solver_mod.f90
# SOURCES := $(sort $(wildcard *_mod.f90))
OBJECTS := $(SOURCES:.f90=.o)

all: $(TARGET_EXE)

$(TARGET_EXE): $(TARGET_EXE).o lib
	$(FC) -I. $< -L. -l$(TARGET_LIB) -o $@

lib: $(OBJECTS)
	ar rcs lib$(TARGET_LIB).a $^

# Create dependencies
%.P: %.f90
	-$(FC) $(DEPFLAGS)
	grep -vi -E '(netcdf|plplot)' $@.tmp > $@
	rm -rf $@.tmp

DEPS := $(SOURCES:.f90=.P) $(TARGET_EXE).P

depends: $(DEPS)
	cat $^ > $@

# Include the dependencies
-include depends

# Actually compile a file
%.o %.mod: %.f90
	$(FC) $(FFLAGS) -c $<

.PHONY: clean

clean:
	rm -rf *.o *.mod *.tmp *.P depends $(TARGET_EXE) lib$(TARGET_LIB).a
