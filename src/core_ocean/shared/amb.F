module amb
  ! call stack starting top-level loop:
  ! ~/mpas/mpas/src/core_ocean/mode_forward/mpas_ocn_forward_mode.F
  ! ~/mpas/mpas/src/core_ocean/mode_forward/mpas_ocn_time_integration.F
  ! ~/mpas/mpas/src/core_ocean/mode_forward/mpas_ocn_time_integration_split.F
  ! ~/mpas/mpas/src/core_ocean/shared/mpas_ocn_tendency.F
  ! ~/mpas/mpas/src/core_ocean/shared/mpas_ocn_tracer_advection.F
  
  implicit none

  logical :: first = .true.
  integer :: mpirank
  logical :: amroot

  integer, parameter :: supercycle = 5
  integer :: timestep

contains

  subroutine amb_set_timestep(ts)
    integer, intent(in) :: ts
    call init()
    timestep = ts
  end subroutine amb_set_timestep

  function amb_get_timestep() result(ts)
    integer :: ts
    ts = timestep
  end function amb_get_timestep

  function amb_get_supercycle() result(sc)
    integer :: sc
    sc = supercycle
  end function amb_get_supercycle

  function amb_amroot() result(ar)
    logical :: ar
    ar = amroot
  end function amb_amroot

  subroutine init()
    use mpi
    integer :: err
    if (.not. first) return
    first = .false.
    call mpi_comm_rank(MPI_COMM_WORLD, mpirank, err)
    amroot = mpirank == 0
    if (amroot) print *,'amb> init'
  end subroutine init

end module amb
