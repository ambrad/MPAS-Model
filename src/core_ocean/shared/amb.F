module amb
  use mpas_kind_types
  use ocn_config
  use mpi

  implicit none
  private

  ! settings
  integer :: supercycle
  logical :: sc_use
  ! initialization
  logical :: inited = .false.
  integer :: mpirank, comm_global
  logical :: amroot
  ! debug
  integer, parameter :: ndbgq = 3, verbose = 0
  ! state
  integer :: timestep, split_explicit_step, n_ts_iter

  public :: &
       amb_use, &
       amb_set_timestep, &
       amb_set_split_explicit_step, &
       amb_get_timestep, &
       amb_get_supercycle, &
       amb_is_tracer_step, &
       amb_set_layerThickness0, &
       amb_accum, &
       amb_fill, &
       amb_update, &
       amb_check_debug_tracers, &
       ! for now
       amb_init, amb_amroot

contains

  function amb_use() result(u)
    logical :: u
    u = sc_use
  end function amb_use

  subroutine amb_set_timestep(comm)
    integer, intent(in) :: comm
    call amb_init(comm)
    timestep = timestep + 1
  end subroutine amb_set_timestep

  subroutine amb_set_split_explicit_step(split_explicit_step_in, n_ts_iter_in)
    integer, intent(in) :: split_explicit_step_in, n_ts_iter_in
    split_explicit_step = split_explicit_step_in
    n_ts_iter = n_ts_iter_in
  end subroutine amb_set_split_explicit_step

  function amb_get_timestep() result(ts)
    integer :: ts
    ts = timestep
  end function amb_get_timestep

  function amb_get_supercycle() result(sc)
    integer :: sc
    sc = supercycle
  end function amb_get_supercycle

  function amb_is_tracer_step() result(b)
    logical :: b
    b = modulo(timestep, supercycle) == 0
  end function amb_is_tracer_step

  function amb_amroot() result(ar)
    logical :: ar
    ar = amroot
  end function amb_amroot

  function amb_get_rank() result(r)
    integer :: r
    r = mpirank
  end function amb_get_rank

  subroutine amb_init(comm)
    integer, intent(in), optional :: comm

    integer :: err

    if (inited) return
    inited = .true.

    if (present(comm)) then
       comm_global = comm
    else
       comm_global = MPI_COMM_WORLD
    end if
    timestep = 0

    call mpi_comm_rank(comm_global, mpirank, err)
    amroot = mpirank == 0

    sc_use = .false.
    supercycle = config_passive_tracers_supercycle_timestep
    sc_use = supercycle > 1 .and. &
         (    trim(config_time_integrator) == 'split_explicit'   &
         .or. trim(config_time_integrator) == 'unsplit_explicit')
    
    if (amroot) then
       print '(a,l,a,i2,a,l)', 'amb> sc use', sc_use, &
            ' sc', supercycle, ' subcyc', config_ecosysTracers_subcycle_param
    end if
  end subroutine amb_init

  subroutine amb_set_layerThickness0(tendPool, layerThickness)
    type (mpas_pool_type), intent(inout) :: tendPool
    real (kind=RKIND), intent(in) :: layerThickness(:,:)

    real (kind=RKIND), pointer :: layerThickness0(:,:)

    if (split_explicit_step /= n_ts_iter) return

    if (modulo(timestep-1, supercycle) == 0) then
       call mpas_pool_get_array(tendPool, 'scLayerThickness0', layerThickness0)
       layerThickness0(:,:) = layerThickness(:,:)
    end if
  end subroutine amb_set_layerThickness0

  subroutine amb_accum(tendPool, layerThickness, layerThicknessEdge, zMid, normalThicknessFlux, &
       vertAleTransportTop)
    type (mpas_pool_type), intent(inout) :: tendPool
    real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness, layerThicknessEdge, zMid, normalThicknessFlux, vertAleTransportTop
    real (kind=RKIND), dimension(:,:), pointer :: test, t1

    real (kind=RKIND), dimension(:,:), pointer :: accumLayerThickness, accumLayerThicknessEdge, &
         accumZMid, accumNormalThicknessFlux, accumVertAleTransportTop

    if (split_explicit_step /= n_ts_iter) return

    call mpas_pool_get_array(tendPool, 'scAccumLayerThickness', accumLayerThickness)
    call mpas_pool_get_array(tendPool, 'scAccumLayerThicknessEdge', accumLayerThicknessEdge)
    call mpas_pool_get_array(tendPool, 'scAccumZMid', accumZMid)
    call mpas_pool_get_array(tendPool, 'scAccumNormalThicknessFlux', accumNormalThicknessFlux)
    call mpas_pool_get_array(tendPool, 'scAccumVertAleTransportTop', accumVertAleTransportTop)

    if (modulo(timestep-1, supercycle) == 0) then
       if (amroot .and. verbose > 0) print *,'amb> accum zero'
       accumLayerThickness = 0
       accumLayerThicknessEdge = 0
       accumZMid = 0
       accumNormalThicknessFlux = 0
       accumVertAleTransportTop = 0
    end if

    if (amroot .and. verbose > 0) print *,'amb> accum accum',timestep,split_explicit_step
    accumZMid = accumZMid + zMid
    accumLayerThickness = accumLayerThickness + layerThickness
    accumNormalThicknessFlux = accumNormalThicknessFlux + normalThicknessFlux
    accumVertAleTransportTop = accumVertAleTransportTop + vertAleTransportTop
    accumLayerThicknessEdge = accumLayerThicknessEdge + layerThicknessEdge

    if (amb_is_tracer_step()) then
       accumLayerThickness = accumLayerThickness/supercycle
       accumZMid = accumZMid/supercycle
       accumNormalThicknessFlux = accumNormalThicknessFlux/supercycle
       accumVertAleTransportTop = accumVertAleTransportTop/supercycle
       accumLayerThicknessEdge = accumLayerThicknessEdge/supercycle
    end if
  end subroutine amb_accum

  subroutine amb_fill(tendPool, layerThickness, layerThicknessEdge, zMid, &
       normalThicknessFlux, vertAleTransportTop)
    type (mpas_pool_type), intent(inout) :: tendPool
    real (kind=RKIND), intent(out) :: normalThicknessFlux(:,:)
    real (kind=RKIND), intent(out), dimension(:,:), pointer :: layerThickness, &
         layerThicknessEdge, zMid, vertAleTransportTop

    real (kind=RKIND), dimension(:,:), pointer :: accumNormalThicknessFlux

    if (split_explicit_step /= n_ts_iter) print *, 'amb> amb_fill: should not be here'
    if (amroot .and. verbose > 0) print *,'amb> fill',timestep,split_explicit_step

    call mpas_pool_get_array(tendPool, 'scAccumNormalThicknessFlux', accumNormalThicknessFlux)
    normalThicknessFlux = accumNormalThicknessFlux
    call mpas_pool_get_array(tendPool, 'scAccumLayerThickness', layerThickness)
    call mpas_pool_get_array(tendPool, 'scAccumLayerThicknessEdge', layerThicknessEdge)
    call mpas_pool_get_array(tendPool, 'scAccumZMid', zMid)
    call mpas_pool_get_array(tendPool, 'scAccumVertAleTransportTop', vertAleTransportTop)
  end subroutine amb_fill

  subroutine amb_update(tendPool, maxLevelCell, nCells, dt, tcur, ttend, lyrnew, tnew)
    type (mpas_pool_type), intent(inout) :: tendPool
    integer, intent(in) :: maxLevelCell(:), nCells
    real (kind=RKIND), intent(in) :: dt, tcur(:,:,:), ttend(:,:,:), lyrnew(:,:)
    real (kind=RKIND), intent(out) :: tnew(:,:,:)

    real (kind=RKIND), pointer :: layerThickness0(:,:)
    integer :: ic, k

    if (modulo(timestep, supercycle) /= 0) then
       tnew(:,:,:) = tcur(:,:,:)
       return
    end if
    call mpas_pool_get_array(tendPool, 'scLayerThickness0', layerThickness0)
    do ic = 1, nCells
       do k = 1, maxLevelCell(ic)
          tnew(:,k,ic) = (tcur(:,k,ic)*layerThickness0(k,ic) + &
                          supercycle*dt*ttend(:,k,ic)) / &
                         lyrnew(k,ic)
       end do
    end do
  end subroutine amb_update

  subroutine amb_check_debug_tracers(tendPool, maxLevelCell, nCells, areaCell, h, q)
    type (mpas_pool_type), intent(inout) :: tendPool
    integer, intent(in) :: maxLevelCell(:), nCells
    real (kind=RKIND), intent(in) :: areaCell(:), q(:,:,:)
    real (kind=RKIND), intent(in), target :: h(:,:)

    integer, parameter :: conout = 999

    real (kind=RKIND), save :: overall_maxerr = 0, mass0(ndbgq)
    integer, save :: consistency_output_cnt = 0, next_timestep = 1

    integer :: idx, ic, k, iq, ierr
    real (kind=RKIND) :: maxerr, mass, send(ndbgq), recv(ndbgq), err
    real (kind=RKIND), pointer :: huse(:,:)

    if (.false. .and. modulo(timestep-1, 10) == 0) then
       ! Debug tracer 3 mixing ratio is uniformly 42. Check departure from that.
       maxerr = 0
       do ic = 1, nCells
          do k = 1, maxLevelCell(ic)
             maxerr = max(maxerr, abs(q(3,k,ic) - 42.0d0))
          end do
       end do
       if (maxerr > 1.0d-14*42.0d0 .and. maxerr > overall_maxerr .and. &
           consistency_output_cnt < conout) then
          print *, 'amb> q consistency', timestep, mpirank, maxerr
          consistency_output_cnt = consistency_output_cnt + 1
          if (consistency_output_cnt == conout .and. amroot) &
               print *,'amb> q consistency silenced'
       end if
       overall_maxerr = max(overall_maxerr, maxerr)
    end if

    if (timestep == next_timestep) then
       if (timestep == 1 .and. sc_use) then
          call mpas_pool_get_array(tendPool, 'scLayerThickness0', huse)
       else
          huse => h
       end if
       do iq = 1, ndbgq
          mass = 0
          do ic = 1, nCells
             do k = 1, maxLevelCell(ic)
                mass = mass + areaCell(ic)*huse(k,ic)*q(iq,k,ic)
             end do
          end do
          send(iq) = mass
       end do
       call mpi_allreduce(send, recv, ndbgq, MPI_REAL8, MPI_SUM, comm_global, ierr)
       if (amroot) then
          do iq = 1, ndbgq
             mass = recv(iq)
             if (timestep == 1) then
                mass0(iq) = mass
                print '(a,i2,es23.15)','amb> mass0', iq, mass
             else
                err = abs(mass - mass0(iq))/mass0(iq)
                if (err > 1e-14) then
                   print '(a,i2,i4,es11.3,es23.15)', 'amb> q mass', iq, timestep, err, mass
                end if
             end if
          end do
       end if
       next_timestep = timestep + 50
       if (supercycle > 0) then
          do while (modulo(next_timestep, supercycle) /= 0)
             next_timestep = next_timestep + 1
          end do
       end if
    end if
  end subroutine amb_check_debug_tracers

end module amb
