# ==============================================================================
# Makefile — Simple A -> B PFR (Wehner-Wilhelm baseline for Ch 4 Test 2)
# ==============================================================================
# Same build convention as [1_isothermal_model]/Makefile; only TARGET differs
# so the two binaries can coexist on PATH.
#
#   make          - compile
#   make run      - compile + run (requires --solver=<name>)
#   make clean    - remove build artifacts
# ==============================================================================

FC            = gfortran
FFLAGS        = -O2 -Wall -Wextra -fcheck=all -fbacktrace -std=f2018
FFLAGS_LEGACY = -O2 -std=legacy -w
SRCDIR        = src
LIBDIR        = lib
BUILDDIR      = build
TARGET        = pfr_simple

# Fortran 90 sources (dependency order: params -> model, solvers -> output -> main)
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

$(BUILDDIR)/%.o: $(SRCDIR)/%.f90 | $(BUILDDIR)
	$(FC) $(FFLAGS) -J$(BUILDDIR) -I$(BUILDDIR) -c $< -o $@

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
	./$(TARGET) --solver=odepack

clean:
	rm -rf $(BUILDDIR) $(TARGET) *.mod
	@echo "  (output/<solver>_<N>/ run directories are preserved; delete manually if desired)"

.PHONY: all run clean
