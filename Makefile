FC = ifort -O3 -fp-model precise -fp-model source -DCPL
MPI_FC = mpif90 -O3 -fp-model precise -fp-model source -DCPL
###NETCDF=/home/netapp/clima-users/users/uturunco/progs/netcdf-4.3.0
NETCDF=/opt-ictp/ESMF/7.1.0r
FCFLAGS = -I$(NETCDF)/include
CPPFLAGS = -L$(NETCDF)/lib -lnetcdf -lnetcdff
APP = main.x

SRC = mod_param.f90 \
      mod_time.f90    \
      mod_io.f90    \
      mod_model.f90 \
      mod_iface.f90 \
      mod_global.f90 \
      mod_mpimess.f90 \
      mod_varandtypes.f90 \
      main.f90

OBJ = $(SRC:.f90=.o)

$(APP): $(OBJ) mod_iface.o 
	$(MPI_FC) -o $(APP) $(OBJ) $(CPPFLAGS)
##	$(FC) -o $(APP) $(OBJ) $(CPPFLAGS)

%.o: %.f90
	$(MPI_FC) $(FCFLAGS) -c $<

%.o: %.f
	$(FC) $(FCFLAGS) -c $<

install: $(APP)

clean:
	rm -f $(APP) *.o *.mod

mod_param.o: mod_param.f90
mod_varandtypes.o: mod_varandtypes.f90
mod_mpimess.o: mod_mpimess.f90 mod_varandtypes.o
mod_time.o: mod_time.f90 mod_param.o
mod_io.o: mod_io.f90 mod_param.o mod_varandtypes.o mod_mpimess.o mod_time.o
mod_model.o: mod_model.f90 mod_param.o mod_varandtypes.o mod_mpimess.o
mod_iface.o: mod_iface.f90 mod_param.o mod_io.o mod_model.o mod_varandtypes.o mod_mpimess.o mod_time.o
main.o: main.f90 mod_param.o mod_iface.o mod_varandtypes.o mod_mpimess.o
