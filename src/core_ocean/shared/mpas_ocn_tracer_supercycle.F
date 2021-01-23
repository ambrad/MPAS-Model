! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_supercycle
!
!> \brief  Supercycle passive tracers w.r.t. baroclinic time step.
!> \author Andrew Bradley
!> \date   January 2021
!> \details
!>  This module provides book-keeping to supercycle passive tracers w.r.t. the
!>  baroclinic time step. The parameter that controls this is TODO TODO TODO
!
!-----------------------------------------------------------------------

module ocn_tracer_supercycle

  use mpi

  implicit none
  private

  integer, parameter :: ndbgq = 3, verbose = 0

  integer :: supercycle, mpirank, timestep, split_explicit_step, n_ts_iter
  logical :: use, amroot

  logical :: first = .true.

  real(8), allocatable :: layerThickness0(:,:), accumNormalThicknessFlux(:,:), &
       accumVertAleTransportTop(:,:)

  public :: &
       ocn_tracer_supercycle_use, &
       ocn_tracer_supercycle_set_timestep, &
       ocn_tracer_supercycle_set_split_explicit_step, &
       ocn_tracer_supercycle_get_supercycle, &
       ocn_tracer_supercycle_set_layerThickness0, &
       ocn_tracer_supercycle_accum, &
       ocn_tracer_supercycle_fill, &
       ocn_tracer_supercycle_update, &
       ocn_tracer_supercycle_check_tracers

contains

  function ocn_tracer_supercycle_use() result(u)
    logical :: u
    u = use
  end function ocn_tracer_supercycle_use

  subroutine ocn_tracer_supercycle_set_timestep()
    call ots_init()
    timestep = timestep + 1
  end subroutine ocn_tracer_supercycle_set_timestep

  subroutine ocn_tracer_supercycle_set_split_explicit_step(split_explicit_step_in, &
       n_ts_iter_in)
    integer, intent(in) :: split_explicit_step_in, n_ts_iter_in
    split_explicit_step = split_explicit_step_in
    n_ts_iter = n_ts_iter_in
  end subroutine ocn_tracer_supercycle_set_split_explicit_step

  function ocn_tracer_supercycle_get_timestep() result(ts)
    integer :: ts
    ts = timestep
  end function ocn_tracer_supercycle_get_timestep

  function ocn_tracer_supercycle_get_supercycle() result(sc)
    integer :: sc
    sc = supercycle
  end function ocn_tracer_supercycle_get_supercycle

  function ocn_tracer_supercycle_is_tracer_step() result(b)
    logical :: b
    b = modulo(timestep, supercycle) == 0
  end function ocn_tracer_supercycle_is_tracer_step

  subroutine ocn_tracer_supercycle_set_layerThickness0(layerThickness)
    real(8), intent(in) :: layerThickness(:,:)

    if (split_explicit_step /= n_ts_iter) return

    if (.not. allocated(layerThickness0)) &
       allocate(layerThickness0(size(layerThickness,1), size(layerThickness,2)))

    if (modulo(timestep-1, supercycle) == 0) &
       layerThickness0(:,:) = layerThickness(:,:)
  end subroutine ocn_tracer_supercycle_set_layerThickness0

  subroutine ocn_tracer_supercycle_accum(normalThicknessFlux, vertAleTransportTop)
    real(8), intent(in) :: normalThicknessFlux(:,:), vertAleTransportTop(:,:)

    if (split_explicit_step /= n_ts_iter) return

    if (.not. allocated(accumNormalThicknessFlux)) then
       allocate(accumNormalThicknessFlux(size(normalThicknessFlux,1), &
                                         size(normalThicknessFlux,2)))
       allocate(accumVertAleTransportTop(size(vertAleTransportTop,1), &
                                         size(vertAleTransportTop,2)))
       if (amroot .and. verbose > 0) print *,'amb> accum alloc'
    end if

    if (modulo(timestep-1, supercycle) == 0) then
       if (amroot .and. verbose > 0) print *,'amb> accum zero'
       accumNormalThicknessFlux(:,:) = 0
       accumVertAleTransportTop(:,:) = 0
    end if

    if (amroot .and. verbose > 0) print *,'amb> accum accum',timestep,split_explicit_step
    accumNormalThicknessFlux = accumNormalThicknessFlux + normalThicknessFlux
    accumVertAleTransportTop = accumVertAleTransportTop + vertAleTransportTop
  end subroutine ocn_tracer_supercycle_accum

  subroutine ocn_tracer_supercycle_fill(normalThicknessFlux, vertAleTransportTop)
    real(8), intent(out) :: normalThicknessFlux(:,:), vertAleTransportTop(:,:)

    if (split_explicit_step /= n_ts_iter) &
         print *, 'amb> ocn_tracer_supercycle_fill: should not be here'

    if (amroot .and. verbose > 0) print *,'amb> fill',timestep,split_explicit_step
    normalThicknessFlux = accumNormalThicknessFlux/supercycle
    vertAleTransportTop = accumVertAleTransportTop/supercycle
  end subroutine ocn_tracer_supercycle_fill

  subroutine ocn_tracer_supercycle_update(maxLevelCell, nCells, dt, tcur, ttend, lyrnew, tnew)
    integer, intent(in) :: maxLevelCell(:), nCells
    real(8), intent(in) :: dt, tcur(:,:,:), ttend(:,:,:), lyrnew(:,:)
    real(8), intent(out) :: tnew(:,:,:)

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
  end subroutine ocn_tracer_supercycle_update

  subroutine ocn_tracer_supercycle_check_tracers(maxLevelCell, nCells, areaCell, h, q)
    integer, intent(in) :: maxLevelCell(:), nCells
    real(8), intent(in) :: areaCell(:), h(:,:), q(:,:,:)

    integer, parameter :: conout = 999

    real(8), save :: overall_maxerr = 0, mass0(ndbgq)
    integer, save :: consistency_output_cnt = 0

    integer :: idx, ic, k, iq, ierr
    real(8) :: maxerr, mass, send(ndbgq), recv(ndbgq), err

    if (modulo(timestep-1, 10) == 0) then
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

    if (modulo(timestep-1, 50) == 0) then
       do iq = 1, ndbgq
          mass = 0
          do ic = 1, nCells
             do k = 1, maxLevelCell(ic)
                mass = mass + areaCell(ic)*h(k,ic)*q(iq,k,ic)
             end do
          end do
          send(iq) = mass
       end do
       call mpi_allreduce(send, recv, ndbgq, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, ierr)
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
    end if
  end subroutine ocn_tracer_supercycle_check_tracers

  subroutine ots_init()
    integer :: err
    character (len=255) :: amb_sc

    if (.not. first) return
    first = .false.

    timestep = 0

    call mpi_comm_rank(MPI_COMM_WORLD, mpirank, err)
    amroot = mpirank == 0

    use = .false.
    supercycle = 1
    call get_environment_variable("amb_sc", amb_sc, status=err)
    if (err /= 1) then
       read(amb_sc, '(i)', iostat=err) supercycle
       use = supercycle > 0
    end if

    if (amroot) print '(a,l,a,i2)', 'amb> init use', use, ' supercycle', supercycle
  end subroutine ots_init

end module ocn_tracer_supercycle
