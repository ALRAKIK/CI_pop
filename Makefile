# Basic directories
BDIR = bin
ODIR = obj
SDIR = src
MODULE_DIR = $(SDIR)/module

# Create necessary directories
$(shell mkdir -p $(BDIR) $(ODIR)/module)

# Compiler and flags
FC = gfortran
FFLAGS = -Wall -Wno-unused -Wno-unused-dummy-argument -O3 -march=native  -g -fbacktrace  -fcheck=all -fimplicit-none -fopenmp
MODDIR = -J$(ODIR) -I$(ODIR)  # Put and find modules in obj directory
LIBS   = src/lib/libquadpack.a   -lblas -llapack -L/opt/homebrew/Cellar/gsl/2.8/lib -lgsl -lgslcblas -lm -fopenmp

# Find all module files
MODULE_SRC = $(wildcard $(MODULE_DIR)/*.f90)
MODULE_OBJ = $(patsubst $(SDIR)/%.f90,$(ODIR)/%.o,$(MODULE_SRC))

# Find all other source files excluding modules
OTHER_SRC = $(shell find $(SDIR) -name "*.f90" -not -path "$(MODULE_DIR)/*")
OTHER_OBJ = $(patsubst $(SDIR)/%.f90,$(ODIR)/%.o,$(OTHER_SRC))

# Main program file (assuming it's in the root directory)
MAIN_SRC = $(wildcard *.f90)
MAIN_OBJ = $(patsubst %.f90,$(ODIR)/%.o,$(MAIN_SRC))

# All object files
ALL_OBJ = $(MODULE_OBJ) $(OTHER_OBJ) $(MAIN_OBJ)

# Default target
all: $(BDIR)/CI_pop

# Target for modules only
modules: $(MODULE_OBJ)

# Compile the executable, ensuring modules are built first
$(BDIR)/CI_pop: modules $(OTHER_OBJ) $(MAIN_OBJ)
	$(FC) -o $@ $(ALL_OBJ) $(LIBS) 
	
# Rule to compile module files
$(ODIR)/module/%.o: $(MODULE_DIR)/%.f90
	@mkdir -p $(dir $@)
	$(FC) $(MODDIR) -c $< -o $@ $(FFLAGS)

# Rule to compile other source files
$(ODIR)/%.o: $(SDIR)/%.f90
	@mkdir -p $(dir $@)
	$(FC) $(MODDIR) -c $< -o $@ $(FFLAGS)

# Rule to compile main program files
$(ODIR)/%.o: %.f90
	$(FC) $(MODDIR) -c $< -o $@ $(FFLAGS)

# Clean everything
clean:
	rm -rf $(ODIR)/* $(BDIR)/CI_pop

.PHONY: all modules clean
