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

  use amb, only: amb_init, amb_amroot

  implicit none
  private

  public :: act_init_setup, act_init_validate, act_vel_tend

  real(RKIND), parameter :: pi = 3.141592653589793d0

  type :: act_t
     integer :: state = 0
     real(RKIND) :: sim_time
  end type act_t

  type (act_t) :: a

contains

  subroutine act_init_setup(domain, err)
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
    call act_init()
    if (amb_amroot()) print *,'amb> act_init_setup'

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
             do k = 1, nVertLevels
                debugTracers(1, k, iCell) = cos(2*pi*refZMid(k)/refZMid(nVertLevels))* &
                     xyztrig(lonCell(iCell), latCell(iCell))
                debugTracers(2, k, iCell) = cos(2*pi*refZMid(k)/refZMid(nVertLevels))* &
                     cosine_bells(lonCell(iCell), latCell(iCell))
                debugTracers(3, k, iCell) = 42.0d0
             end do
          endif

          do k = 1, nVertLevels
             layerThickness(k, iCell) = 100.0_RKIND 
             restingThickness(k, iCell) = layerThickness(k, iCell)
          end do

          fCell(iCell) = 0.0_RKIND ! 0 Coriolis
          bottomDepth(iCell) = 100.0_RKIND
          maxLevelCell(iCell) = nVertLevels
       end do

       do iEdge = 1, nEdgesSolve
          fEdge(iEdge) = 0.0_RKIND
       end do

       do iVertex=1, nVerticesSolve
          fVertex(iVertex) = 0.0_RKIND
       end do

       do iEdge = 1,nEdgesSolve
          normalVelocity(:,iEdge) = 0
       enddo

       block_ptr => block_ptr % next
    end do
  end subroutine act_init_setup

  subroutine act_init_validate(cfgpl, pkgpl, iocontext, err)
    type (mpas_pool_type), intent(inout) :: cfgpl, pkgpl
    type (mpas_io_context_type), intent(inout) :: iocontext
    integer, intent(out) :: err

    integer, pointer :: config_vert_levels

    err = 0

    call mpas_pool_get_config(cfgpl, 'config_vert_levels', config_vert_levels)
    config_vert_levels = 3
  end subroutine act_init_validate

  subroutine act_vel_tend(meshPool, days_since_start, dt, normalVelocity, tend_normalVelocity)
    type (mpas_pool_type), intent(in) :: meshPool
    real(RKIND), intent(in) :: days_since_start, dt, normalVelocity(:,:)
    real(RKIND), intent(inout) :: tend_normalVelocity(:,:)

    real(RKIND), dimension(:), pointer :: latEdge, lonEdge
    integer, pointer :: nEdges
    integer :: i

    call act_init()
    if (a%state /= 2) return

    if (amb_amroot()) print *,'amb> days', days_since_start, 'dt', dt
    call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
    call mpas_pool_get_array(meshPool, 'latEdge', latEdge)
    call mpas_pool_get_array(meshPool, 'lonEdge', lonEdge)
    do i = 1, nEdges
       tend_normalVelocity(:,i) = 0
    end do
  end subroutine act_vel_tend

  subroutine act_init()
    character (len=255) :: key
    integer :: err, val

    if (a%state == 2) return
    call amb_init()
    if (a%state == 0) then
       a%state = 1
       call get_environment_variable("amb_act", key, status=err)
       if (err /= 1) then
          read(key, '(i)', iostat=err) val
          if (val /= 0) a%state = 2
       end if
       if (amb_amroot()) print *,'amb> act', a%state == 2
    end if
  end subroutine act_init

  subroutine lonlat2xyz(lon, lat, x, y, z)
    real(RKIND), intent(in) :: lon, lat
    real(RKIND), intent(out) :: x, y, z
    real(RKIND) :: sinl, cosl
    sinl = sin(lat)
    cosl = cos(lat)
    x = cos(lon)*cosl
    y = sin(lon)*cosl
    z = sinl
  end subroutine lonlat2xyz

  subroutine xyz2lonlat(x, y, z, lon, lat)
    real(RKIND), intent(in) :: x, y, z
    real(RKIND), intent(out) :: lon, lat
    real(RKIND) :: r
    lon = atan2(y, x);
    r = sqrt(x*x + y*y + z*z)
    lat = asin(z/r);
  end subroutine xyz2lonlat

  function great_circle_dist_xyz(xA, yA, zA, xB, yB, zB) result(dist)
    real(RKIND), intent(in) :: xA, yA, zA, xB, yB, zB
    real(RKIND) :: dist, cp1, cp2, cp3, cpnorm, dotprod
    cp1 = yA*zB - yB*zA
    cp2 = xB*zA - xA*zB
    cp3 = xA*yB - xB*yA
    cpnorm = sqrt(cp1*cp1 + cp2*cp2 + cp3*cp3)
    dotprod = xA*xB + yA*yB + zA*zB
    dist = atan2(cpnorm, dotprod)
  end function great_circle_dist_xyz

  function great_circle_dist_ll(lon1, lat1, lon2, lat2, radius) result(dist)
    real(RKIND), intent(in) :: lon1, lat1, lon2, lat2, radius
    real(RKIND) :: dist, xA, yA, zA, xB, yB, zB, cp1, cp2, cp3, cpnorm, dotprod
    call lonlat2xyz(lon1, lat1, xA, yA, zA)
    call lonlat2xyz(lon2, lat2, xB, yB, zB)
    dist = radius*great_circle_dist_xyz(xA, yA, zA, xB, yB, zB)
  end function great_circle_dist_ll

  function cosine_bell(ri, r) result(cb)
    real(RKIND), intent(in) :: ri, r
    real(RKIND) :: cb
    real(RKIND), parameter :: hmax = 1._RKIND
    cb = 0.5_RKIND*hmax*(1 + cos(pi*ri/r))
  end function cosine_bell

  function cosine_bells(lon, lat) result (cbs)
    real(RKIND), intent(in) :: lon, lat
    real(RKIND), parameter :: r = 0.5_RKIND, b = 0.1_RKIND, c = 0.9_RKIND, &
         lon1 = 5*(pi/6), lat1 = 0, lon2 = -5*(pi/6), lat2 = 0
    real(RKIND) :: cbs, x, y, z, x1, y1, z1, x2, y2, z2, r1, r2
    call lonlat2xyz(lon, lat, x, y, z)
    call lonlat2xyz(lon1, lat1, x1, y1, z1)
    r1 = great_circle_dist_xyz(x, y, z, x1, y1, z1)
    cbs = 0
    if (r1 < r) then
       cbs = cosine_bell(r1, r)
    else
       call lonlat2xyz(lon2, lat2, x2, y2, z2)
       r2 = great_circle_dist_xyz(x, y, z, x2, y2, z2)
       if (r2 < r) cbs = cosine_bell(r2, r)
    end if
  end function cosine_bells

  function xyztrig(lon, lat) result (f)
    real(RKIND), intent(in) :: lon, lat
    real(RKIND) :: f, x, y, z
    call lonlat2xyz(lon, lat, x, y, z)
    f = 0.5_RKIND*(1 + sin(pi*x)*sin(pi*y)*sin(pi*z))
  end function xyztrig

end module act
