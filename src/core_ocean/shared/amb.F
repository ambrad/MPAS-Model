module amb
  implicit none

  integer :: mpi_rank
  logical :: amroot

  logical :: first = .true.
contains

  subroutine amb_hi()
    call init()
  end subroutine amb_hi

  subroutine init()
    use mpi
    integer :: err
    if (.not. first) return
    first = .false.
    call mpi_comm_rank(MPI_COMM_WORLD, mpi_rank, err)
    amroot = mpi_rank == 0
    if (amroot) print *,'amb> init'
  end subroutine init
end module amb
