!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP CHyM.
!
!    ICTP CHyM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published
!    by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP CHyM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP CHyM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

program main

  use mod_iface
  use mod_mpimess
  use mod_varandtypes

  implicit none

  call mpi_init(mpierr)
  call mpi_comm_dup(mpi_comm_world,mycomm,mpierr)
  call mpi_comm_size(mycomm, nproc, mpierr)
  call mpi_comm_rank(mycomm, myid, mpierr)
  call chym_init( )
  call chym_run( )
  call chym_close( )

end program main

