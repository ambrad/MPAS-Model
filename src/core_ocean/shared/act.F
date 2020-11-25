module act
  ! amb_sc   if defined, integer supercycle value
  ! amb_act  if defined, any non-0 value

  use mpas_kind_types
  use mpas_io_units
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_constants
  use mpas_stream_manager
  use mpas_dmpar

  use amb, only: amb_init, amb_amroot, amb_get_rank

  implicit none
  private

  public :: act_vel_tend

  real (RKIND), parameter :: pi = 3.141592653589793d0

  type :: act_t
     integer :: state = 0
  end type act_t

  type (act_t) :: a

contains

  subroutine act_init(domain, err)
    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: err

    type (block_type), pointer :: block_ptr
    type (mpas_pool_type), pointer :: meshPool, verticalMeshPool, statePool, &
         diagnosticsPool, forcingPool, tracersPool
    integer, pointer :: nVertLevels, nVertLevelsP1, nCellsSolve,nEdgesSolve, nEdge, &
         nVerticesSolve, index_temperature, index_salinity,index_tracer1
    integer, dimension(:), pointer :: maxLevelCell
    real (kind=RKIND), pointer :: sphere_radius
    real (kind=RKIND), dimension(:), pointer :: refBottomDepth,refZMid, lonCell, &
         latEdge, latCell, latVertex, bottomDepth, angleEdge, fCell, fEdge, fVertex, &
         xCell, yCell, interfaceLocations
    real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness, &
         normalVelocity
    real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers,debugTracers
    character (len=StrKIND), pointer :: config_init_configuration

    integer :: iCell, iEdge, iVertex, k, kML
    real (kind=RKIND) :: temp, velocity, temperature, salinity
    
    err = 0
    call amb_init()
    if (amb_amroot()) print *,'amb> act_init'

    block_ptr => domain % blocklist
    do while(associated(block_ptr))
       call mpas_pool_get_subpool(block_ptr % structs, 'mesh',meshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh',verticalMeshPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'state',statePool)
       call mpas_pool_get_subpool(block_ptr % structs, 'diagnostics',diagnosticsPool)
       call mpas_pool_get_subpool(block_ptr % structs, 'forcing',forcingPool)

       call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)

       call mpas_pool_get_dimension(meshPool, 'nVertLevels',nVertLevels)
       call mpas_pool_get_dimension(meshPool, 'nCellsSolve',nCellsSolve)
       call mpas_pool_get_dimension(meshPool, 'nEdge', nEdge)
       call mpas_pool_get_dimension(meshPool, 'nEdgesSolve',nEdgesSolve)
       call mpas_pool_get_dimension(meshPool, 'nVerticesSolve',nVerticesSolve)
       call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)

       call mpas_pool_get_dimension(tracersPool, 'index_temperature',index_temperature)
       call mpas_pool_get_dimension(tracersPool, 'index_salinity',index_salinity)
       call mpas_pool_get_dimension(tracersPool, 'index_tracer1',index_tracer1)

       call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
       call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
       call mpas_pool_get_array(meshPool, 'angleEdge', angleEdge)
       call mpas_pool_get_array(meshPool, 'latCell', latCell)
       call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
       call mpas_pool_get_array(meshPool, 'latVertex', latVertex)
       call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
       call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
       call mpas_pool_get_array(meshPool, 'fCell', fCell)
       call mpas_pool_get_array(meshPool, 'fEdge', fEdge)
       call mpas_pool_get_array(meshPool, 'fVertex', fVertex)
       call mpas_pool_get_array(meshPool, 'xCell', xCell)
       call mpas_pool_get_array(meshPool, 'yCell', yCell)


       call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
       call mpas_pool_get_array(verticalMeshPool, 'restingThickness',restingThickness)

       call mpas_pool_get_array(statePool, 'normalVelocity', normalVelocity, 1)
       call mpas_pool_get_array(tracersPool, 'activeTracers',activeTracers, 1)
       call mpas_pool_get_array(tracersPool, 'debugTracers',debugTracers, 1)
       call mpas_pool_get_array(statePool, 'layerThickness',layerThickness, 1)

       ! Set refBottomDepth and refBottomDepthTopOfCell
       do k = 1, nVertLevels
          refBottomDepth(k) = 100.0_RKIND
          refZMid(k) = 50.0_RKIND
       end do

       do iCell = 1, nCellsSolve
          if (associated(activeTracers) ) then
             activeTracers(index_temperature, :, iCell) = 15.0d0
             activeTracers(index_salinity, :, iCell) = 35.0d0
          endif

          if (associated(debugTracers)) then
             debugTracers(:,:,iCell) = 0.0_RKIND
          endif

          do k = 1, nVertLevels
             layerThickness(k, iCell) = 100.0_RKIND 
             restingThickness(k, iCell) = layerThickness(k, iCell)
          end do

          fCell(iCell) = 0.0_RKIND*omega*sin(latCell(iCell)) ! 0 Coriolis
          bottomDepth(iCell) = 100.0_RKIND
          maxLevelCell(iCell) = nVertLevels
       end do

       do iEdge = 1, nEdgesSolve
          fEdge(iEdge) = 0.0_RKIND*omega*sin(latEdge(iEdge))
       end do

       do iVertex=1, nVerticesSolve
          fVertex(iVertex) = 0.0_RKIND*omega*sin(latVertex(iVertex))
       end do

       do iEdge = 1,nEdgesSolve
          velocity = 2.0*3.14159265*sphere_radius*cos(latEdge(iEdge)) / &
               (86400.0_RKIND*12.0d0)
          normalVelocity(:,iEdge) = velocity*cos(angleEdge(iEdge))
       enddo

       block_ptr => block_ptr % next
    end do
  end subroutine act_init

  subroutine act_init_validate(cfgpl, pkgpl, iocontext, err)
    type (mpas_pool_type), intent(inout) :: cfgpl, pkgpl
    type (mpas_io_context_type), intent(inout) :: iocontext
    integer, intent(out) :: err

    integer, pointer :: config_vert_levels

    err = 0

    call mpas_pool_get_config(cfgpl, 'config_vert_levels', config_vert_levels)
    config_vert_levels = 3
  end subroutine act_init_validate

  subroutine act_vel_tend(meshPool, normalVelocity, tend_normalVelocity)
    type (mpas_pool_type), intent(in) :: meshPool
    real(RKIND), intent(in) :: normalVelocity(:,:)
    real(RKIND), intent(inout) :: tend_normalVelocity(:,:)

    integer :: err, val
    character (len=255) :: key

    if (a%state == 0) then
       a%state = 1
       call get_environment_variable("amb_act", key, status=err)
       if (err /= 1) then
          read(key, '(i)', iostat=err) val
          if (val /= 0) a%state = 2
       end if
    end if

    if (a%state /= 2) return
    if (amb_amroot()) print *,'amb> here'
  end subroutine act_vel_tend

end module act
