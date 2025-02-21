GIT_VERSION := $(shell git describe --abbrev=4 --dirty --always --tags)

FCFLAGS += -I. $(NC_INCLUDES)
CPPFLAGS += -DGIT_REV=\"$(GIT_VERSION)\"

bin_PROGRAMS = preproc

preproc_SOURCES = mod_cellaut.F90  mod_crtstat.F90  mod_param.F90 \
                  mod_common.F90   mod_ncio.F90     mod_vector.F90 \
                  preproc.F90

FCLD = $(FC)

%.o: %.F90
	$(FC) $(CPPFLAGS) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.*

clean-local:
	rm -f *.mod *__genmod.*

mod_vector.o : mod_vector.F90
mod_cellaut.o : mod_cellaut.F90 mod_vector.o
mod_common.o : mod_common.F90
mod_ncio.o : mod_ncio.F90 mod_common.o
mod_param.o : mod_param.F90 mod_common.o mod_ncio.o
mod_crtstat.o : mod_crtstat.F90 mod_common.o
