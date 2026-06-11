# Compiler and flags
FC = ifort
FFLAGS = -O3 -qopt-matmul -axAVX -fp-model fast -static -threads -qopenmp
#FFLAGS = -O3 -qopt-matmul -xCORE-AVX2 -fp-model source -static -threads -qopenmp
LDFLAGS = -qopenmp -qopt-matmul

# Directories
SRCDIR = ..
BINDIR = bin
MODDIR = $(BINDIR)/modules

# Create the binary and module directories
$(shell mkdir -p $(BINDIR) $(MODDIR))

# Main target
all: lobbie_integration

# Compile each Fortran file explicitly, respecting the dependencies
$(BINDIR)/PosVelElements_new.o: $(SRCDIR)/PosVelElements_new.F90
	$(FC) $(FFLAGS) -module $(MODDIR) -o $@ -c $<

$(BINDIR)/ElementsPosVel_new.o: $(SRCDIR)/ElementsPosVel_new.F90
	$(FC) $(FFLAGS) -module $(MODDIR) -o $@ -c $<

$(BINDIR)/Lobbie_II.o: $(SRCDIR)/Lobbie_II.f90
	$(FC) $(FFLAGS) -module $(MODDIR) -o $@ -c $<

$(BINDIR)/Lobbie_test.o: $(SRCDIR)/Lobbie_test.f90
	$(FC) $(FFLAGS) -module $(MODDIR) -o $@ -c $<



# Link the executable
lobbie_integration: $(BINDIR)/PosVelElements_new.o \
	$(BINDIR)/ElementsPosVel_new.o \
	$(BINDIR)/Lobbie_II.o \
	$(BINDIR)/Lobbie_test.o 
	$(FC) -o ../$@ $^ $(LDFLAGS)

# Clean up the build files
clean:
	rm -rf $(BINDIR)/*
