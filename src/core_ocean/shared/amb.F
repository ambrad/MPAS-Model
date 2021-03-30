module amb
  ! advection call stack starting with top-level loop:
  ! ~/mpas/mpas/src/core_ocean/mode_forward/mpas_ocn_forward_mode.F
  ! ~/mpas/mpas/src/core_ocean/mode_forward/mpas_ocn_time_integration.F
  ! ~/mpas/mpas/src/core_ocean/mode_forward/mpas_ocn_time_integration_split.F
  ! ~/mpas/mpas/src/core_ocean/shared/mpas_ocn_tendency.F
  ! ~/mpas/mpas/src/core_ocean/shared/mpas_ocn_tracer_advection.F
  
  use mpas_kind_types
  use mpi

  implicit none

  ! 3 for baro10km with 2hr dt
  ! 2 for baro4km with 1hr dt
  integer :: supercycle
  logical :: sc_use, subc_params, diag_carbon

  logical :: first = .true.
  integer :: mpirank, comm_global
  logical :: amroot

  integer, parameter :: ndbgq = 3, verbose = 0

  integer :: timestep, split_explicit_step, n_ts_iter

  real(RKIND), dimension(:,:), allocatable, target :: &
       layerThickness0, accumLayerThickness, accumLayerThicknessEdge, &
       accumZMid, accumVertAleTransportTop
  real(RKIND), dimension(:,:), allocatable :: &
       accumNormalThicknessFlux

  real(RKIND) :: carbon_diags(5)

contains

  function amb_use() result(u)
    logical :: u
    u = sc_use
  end function amb_use

  function amb_subc_params() result (u)
    logical :: u
    u = subc_params
  end function amb_subc_params

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
    character (len=255) :: env_str
    integer :: tmp

    if (.not. first) return
    first = .false.

    if (present(comm)) then
       comm_global = comm
    else
       comm_global = MPI_COMM_WORLD
    end if
    timestep = 0

    call mpi_comm_rank(comm_global, mpirank, err)
    amroot = mpirank == 0

    sc_use = .false.
    supercycle = 1
    call get_environment_variable("amb_sc", env_str, status=err)
    if (err /= 1) then
       read(env_str, '(i)', iostat=err) supercycle
       sc_use = supercycle > 0!1
    end if
    subc_params = sc_use
    if (sc_use) then
       call get_environment_variable("amb_subc_params", env_str, status=err)
       if (err /= 1) then
          read(env_str, '(i)', iostat=err) tmp
          subc_params = tmp /= 0
       end if
    end if
    diag_carbon = .false.
    call get_environment_variable("amb_diag_carbon", env_str, status=err)
    if (err /= 1) then
       read(env_str, '(i)', iostat=err) tmp
       diag_carbon = tmp /= 0
    end if
    
    if (amroot) then
       print '(a,l,a,l,a,i2,a,l)', 'amb> sc use', sc_use, &
            ' subc params', subc_params, ' sc', supercycle, ' diagC', diag_carbon
    end if
  end subroutine amb_init

  subroutine amb_set_layerThickness0(layerThickness)
    real(RKIND), intent(in) :: layerThickness(:,:)

    if (split_explicit_step /= n_ts_iter) return

    if (.not. allocated(layerThickness0)) then
       allocate(layerThickness0(size(layerThickness,1), size(layerThickness,2)))
    end if

    if (modulo(timestep-1, supercycle) == 0) then
       layerThickness0(:,:) = layerThickness(:,:)
    end if
  end subroutine amb_set_layerThickness0

  subroutine amb_get_layerThickness0(layerThickness)
    real(RKIND), pointer, intent(out) :: layerThickness(:,:)
    layerThickness => layerThickness0
  end subroutine amb_get_layerThickness0

  subroutine amb_accum(layerThickness, layerThicknessEdge, zMid, normalThicknessFlux, &
       vertAleTransportTop)
    real(RKIND), dimension(:,:), intent(in) :: &
         layerThickness, layerThicknessEdge, zMid, normalThicknessFlux, vertAleTransportTop

    if (split_explicit_step /= n_ts_iter) return

    if (.not. allocated(accumNormalThicknessFlux)) then
       allocate(accumLayerThickness(size(layerThickness,1), &
                                    size(layerThickness,2)))
       allocate(accumLayerThicknessEdge(size(layerThicknessEdge,1), &
                                        size(layerThicknessEdge,2)))
       allocate(accumZMid(size(zMid,1), size(zMid,2)))
       allocate(accumNormalThicknessFlux(size(normalThicknessFlux,1), &
                                         size(normalThicknessFlux,2)))
       allocate(accumVertAleTransportTop(size(vertAleTransportTop,1), &
                                         size(vertAleTransportTop,2)))
       if (amroot .and. verbose > 0) print *,'amb> accum alloc'
    end if

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

  subroutine amb_fill(layerThickness, layerThicknessEdge, zMid, normalThicknessFlux, &
       vertAleTransportTop)
    real(RKIND), intent(out) :: normalThicknessFlux(:,:)
    real(RKIND), intent(out), dimension(:,:), pointer :: layerThickness, &
         layerThicknessEdge, zMid, vertAleTransportTop

    if (split_explicit_step /= n_ts_iter) print *, 'amb> amb_fill: should not be here'

    if (amroot .and. verbose > 0) print *,'amb> fill',timestep,split_explicit_step
    normalThicknessFlux = accumNormalThicknessFlux
    zMid => accumZMid
    vertAleTransportTop => accumVertAleTransportTop
    layerThickness => accumLayerThickness
    layerThicknessEdge => accumLayerThicknessEdge
  end subroutine amb_fill

  subroutine amb_update(maxLevelCell, nCells, dt, tcur, ttend, lyrnew, tnew)
    integer, intent(in) :: maxLevelCell(:), nCells
    real(RKIND), intent(in) :: dt, tcur(:,:,:), ttend(:,:,:), lyrnew(:,:)
    real(RKIND), intent(out) :: tnew(:,:,:)

    integer :: ic, k

    if (modulo(timestep, supercycle) /= 0) then
       tnew(:,:,:) = tcur(:,:,:)
       return
    end if
    do ic = 1, nCells
       do k = 1, maxLevelCell(ic)
          tnew(:,k,ic) = (tcur(:,k,ic)*layerThickness0(k,ic) + &
                          supercycle*dt*ttend(:,k,ic)) / &
                         lyrnew(k,ic)
       end do
    end do
  end subroutine amb_update

  subroutine amb_check_debug_tracers(maxLevelCell, nCells, areaCell, h, q)
    integer, intent(in) :: maxLevelCell(:), nCells
    real(RKIND), intent(in) :: areaCell(:), q(:,:,:)
    real(RKIND), intent(in), target :: h(:,:)

    integer, parameter :: conout = 999

    real(RKIND), save :: overall_maxerr = 0, mass0(ndbgq)
    integer, save :: consistency_output_cnt = 0, next_timestep = 1

    integer :: idx, ic, k, iq, ierr
    real(RKIND) :: maxerr, mass, send(ndbgq), recv(ndbgq), err
    real(RKIND), pointer :: huse(:,:)

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
          huse => layerThickness0
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

  !!! carbon conservation
  !   Inputs to reduction routines are already reduced over nCellsSolve w/in
  ! this rank.
  !   carbon_diags = (tendencies
  !                   calcToSed + pocToSed
  !                   surfaceFlux
  !                   diag_Jint_Ctot,
  !                   diff(hnew qnew - hprev qprev)/dt)

  function amb_diag_carbon() result (u)
    logical :: u
    u = diag_carbon
  end function amb_diag_carbon

  subroutine amb_c_reduce_ecosys_tends(tot)
    real(RKIND), intent(in) :: tot
    integer :: ierr

    if (.not. diag_carbon) return
    call mpi_allreduce(tot, carbon_diags(1), 1, MPI_REAL8, MPI_SUM, comm_global, ierr)
  end subroutine amb_c_reduce_ecosys_tends

  subroutine amb_c_reduce_ecosys_diffs(tot)
    real(RKIND), intent(in) :: tot
    integer :: ierr

    if (.not. diag_carbon) return
    call mpi_allreduce(tot, carbon_diags(5), 1, MPI_REAL8, MPI_SUM, comm_global, ierr)
  end subroutine amb_c_reduce_ecosys_diffs

  subroutine amb_c_reduce_ecosys_diff_vmix(tot)
    real(RKIND), intent(in) :: tot(3)
    integer :: ierr
    real(RKIND) :: totg(3)

    if (.not. diag_carbon) return
    call mpi_allreduce(tot, totg, 3, MPI_REAL8, MPI_SUM, comm_global, ierr)
    if (amroot) print '(a,3(es23.16))', 'amb> Cvmix', totg
  end subroutine amb_c_reduce_ecosys_diff_vmix

  subroutine amb_c_reduce_ecosys_sfc_fluxes(tot)
    real(RKIND), intent(in) :: tot
    integer :: ierr

    if (.not. diag_carbon) return
    call mpi_allreduce(tot, carbon_diags(3), 1, MPI_REAL8, MPI_SUM, comm_global, ierr)
  end subroutine amb_c_reduce_ecosys_sfc_fluxes

  subroutine amb_c_reduce_ecosys_diags(calc_to_sed, poc_to_sed, diag_Jint_Ctot)
    real(RKIND), intent(in) :: calc_to_sed, poc_to_sed, diag_Jint_Ctot
    integer :: ierr
    real(RKIND) :: send(2), recv(2)

    if (.not. diag_carbon) return
    send(1) = calc_to_sed + poc_to_sed
    send(2) = diag_Jint_Ctot
    call mpi_allreduce(send, recv, 2, MPI_REAL8, MPI_SUM, comm_global, ierr)
    carbon_diags(2) = recv(1)
    carbon_diags(4) = recv(2)
  end subroutine amb_c_reduce_ecosys_diags

  subroutine amb_c_print()
    if (.not. amroot) return
    if (.not. diag_carbon) return
    print '(a,5(es11.3))', 'amb> Ctend', carbon_diags(1:4), &
         carbon_diags(1) + carbon_diags(2) - carbon_diags(3) + carbon_diags(4)
    print '(a,5(es11.3))', 'amb> Cdiff', carbon_diags(5), carbon_diags(2:4), &
         carbon_diags(5) + carbon_diags(2) - carbon_diags(3) + carbon_diags(4)
  end subroutine amb_c_print

end module amb
