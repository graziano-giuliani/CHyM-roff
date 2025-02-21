module mod_global
  use mpi
  integer :: ierr, rank, nproc
  integer, dimension(MPI_STATUS_SIZE) :: istatus
  integer, dimension(MPI_STATUS_SIZE) :: istatus1
  integer, parameter :: comm = MPI_COMM_WORLD
end module mod_global
