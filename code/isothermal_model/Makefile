# ==============================================================================
# Makefile — Tubular reactor prototype (MOL + ODEPACK/DLSODE)
# ==============================================================================
# Minimal build following Saša's template:
#   make          - compile
#   make run      - compile + run (output to output/output.txt)
#   make plot     - run the Python plotter
#   make runplot  - run + plot
#   make clean    - remove build artifacts
# ==============================================================================

FC            = gfortran
FFLAGS        = -O2 -Wall -Wextra -fcheck=all -fbacktrace -std=f2018
FFLAGS_LEGACY = -O2 -std=legacy -w
SRCDIR        = src
LIBDIR        = lib
BUILDDIR      = build
TARGET        = prototype

# Fortran 90 sources (dependency order: params → model, solvers → output → main)
F90_SRCS = $(SRCDIR)/params.f90 \
           $(SRCDIR)/model.f90 \
           $(SRCDIR)/solvers.f90 \
           $(SRCDIR)/output.f90 \
           $(SRCDIR)/main.f90
F90_OBJS = $(patsubst $(SRCDIR)/%.f90,$(BUILDDIR)/%.o,$(F90_SRCS))

# ODEPACK (Fortran 77, legacy)
F77_SRCS = $(LIBDIR)/opkdmain.f $(LIBDIR)/opkda1.f $(LIBDIR)/opkda2.f
F77_OBJS = $(BUILDDIR)/opkdmain.o $(BUILDDIR)/opkda1.o $(BUILDDIR)/opkda2.o

all: $(TARGET)

$(TARGET): $(F90_OBJS) $(F77_OBJS)
	$(FC) $(FFLAGS) -o $@ $^

# --- Compile F90 modules/program (with module cache in build/) ---
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -J$(BUILDDIR) -I$(BUILDDIR) -c $< -o $@

# --- Compile ODEPACK legacy Fortran 77 ---
$(BUILDDIR)/opkdmain.o: $(LIBDIR)/opkdmain.f | $(BUILDDIR)
	$(FC) $(FFLAGS_LEGACY) -c $< -o $@

$(BUILDDIR)/opkda1.o: $(LIBDIR)/opkda1.f | $(BUILDDIR)
	$(FC) $(FFLAGS_LEGACY) -c $< -o $@

$(BUILDDIR)/opkda2.o: $(LIBDIR)/opkda2.f | $(BUILDDIR)
	$(FC) $(FFLAGS_LEGACY) -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

run: $(TARGET)
	@mkdir -p output
	./$(TARGET)

plot:
	@echo "Plotting..."
	python3 plots/plot.py

runplot: run plot

clean:
	rm -rf $(BUILDDIR) $(TARGET) *.mod
	@echo "  (output/<solver>_<N>/ run directories are preserved; delete manually if desired)"

.PHONY: all run plot runplot clean
