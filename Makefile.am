#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#
#    This file is part of ICTP CHyM.
#
#    ICTP CHyM is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    ICTP CHyM is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with ICTP CHyM.  If not, see <http://www.gnu.org/licenses/>.
#
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

SUBDIRS = preproc

if CHYM_PARALLEL
rcsuffix = MPI
else
rcsuffix = Serial
endif

## FCFLAGS += -I. `nf-config --fflags`
FCFLAGS += -I. $(NC_INCLUDES)

bin_PROGRAMS = chymroff

chymroff_SOURCES = mod_param.F90 \
      mod_time.f90    \
      mod_io.f90    \
      mod_model.F90 \
      mod_iface.f90 \
      mod_global.f90 \
      mod_mpimess.f90 \
      mod_varandtypes.f90 \
      main.f90


FCLD = $(MPIFC)

%.o: %.f90
	$(FCLD) $(CPPFLAGS) $(FCFLAGS) -c $<
%.o: %.F90
	$(FCLD) $(CPPFLAGS) $(FCFLAGS) -c $<

distclean-local:
	rm -f *.mod *__genmod.*

clean-local:
	rm -f *.mod *__genmod.*

mod_param.o: mod_param.F90
mod_varandtypes.o: mod_varandtypes.f90
mod_mpimess.o: mod_mpimess.f90 mod_varandtypes.o
mod_time.o: mod_time.f90 mod_param.o
mod_io.o: mod_io.f90 mod_param.o mod_varandtypes.o mod_mpimess.o mod_time.o
mod_model.o: mod_model.F90 mod_param.o mod_varandtypes.o mod_mpimess.o
mod_iface.o: mod_iface.f90 mod_param.o mod_io.o mod_model.o mod_varandtypes.o mod_mpimess.o mod_time.o
main.o: main.f90 mod_param.o mod_iface.o mod_varandtypes.o mod_mpimess.o
