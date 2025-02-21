FC = mpifort
#FC = mpiifx
# PRODUCTION
FCFLAGS = -O3
#FCFLAGS = -O3 -fp-model precise -fp-model source
# DEBUG
#FCFLAGS = -O0 -g -Wall -pedantic -fcheck=all -fbacktrace
#FCFLAGS = -O0 -g -Warn all -check all -traceback
CPPFLAGS = `nf-config --fflags` -DGIT_REV="'none'"
LIBS = `nf-config --flibs`

all :: preproc

OBJ = mod_vector.o mod_common.o mod_ncio.o mod_param.o mod_cellaut.o \
      mod_crtstat.o

%.o: %.F90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

%.o: %.f90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

preproc: preproc.F90 $(OBJ)
	$(FC) $(CPPFLAGS) $(FCFLAGS) $< -o $@ $(OBJ) $(LIBS)

clean:
	rm -f preproc *.o *.mod

mod_vector.o : mod_vector.F90
mod_cellaut.o : mod_cellaut.F90 mod_vector.o
mod_common.o : mod_common.F90
mod_ncio.o : mod_ncio.F90 mod_common.o
mod_param.o : mod_param.F90 mod_common.o
mod_crtstat.o : mod_crtstat.F90 mod_common.o
