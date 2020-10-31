module act
  use mpas_kind_types
  use mpas_io_units
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_constants
  use mpas_stream_manager
  use mpas_dmpar

  use ocn_init_cell_markers
  use ocn_init_vertical_grids

  use amb, only: amb_init, amb_amroot, amb_get_rank

  implicit none

contains

  subroutine act_init(dmn, err)
    type (domain_type), intent(inout) :: dmn
    integer, intent(out) :: err

    err = 0
    call amb_init()
    if (amb_amroot()) print *,'amb> act_init'
  end subroutine act_init

  subroutine act_init_validate(cfgpl, pkgpl, iocontext, err)
    type (mpas_pool_type), intent(inout) :: cfgpl, pkgpl
    type (mpas_io_context_type), intent(inout) :: iocontext
    integer, intent(out) :: err

    err = 0
  end subroutine act_init_validate

end module act
