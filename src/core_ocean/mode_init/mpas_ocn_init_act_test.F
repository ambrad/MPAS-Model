module mpas_ocn_init_act_test
  use mpas_kind_types
  use mpas_io_units
  use mpas_derived_types
  use mpas_pool_routines
  use mpas_constants
  use mpas_stream_manager
  use mpas_dmpar

  use act

  implicit none

contains

  subroutine ocn_init_setup_act_test(domain, iErr)
    type (domain_type), intent(inout) :: domain
    integer, intent(out) :: iErr

    character (len=StrKIND), pointer :: config_init_configuration

    iErr = 0
    call mpas_pool_get_config(domain%configs,'config_init_configuration', &
         config_init_configuration)
    if (config_init_configuration /= trim('act')) return

    call act_init(domain, iErr)
  end subroutine ocn_init_setup_act_test

  subroutine ocn_init_validate_act_test(configPool, packagePool, iocontext, iErr)
    type (mpas_pool_type), intent(inout) :: configPool
    type (mpas_pool_type), intent(inout) :: packagePool
    type (mpas_io_context_type), intent(inout) :: iocontext
    integer, intent(out) :: iErr

    character (len=StrKIND), pointer :: config_init_configuration

    iErr = 0
    call mpas_pool_get_config(configPool,'config_init_configuration', &
         config_init_configuration)
    if (config_init_configuration /= trim('act')) return

    call act_init_validate(configPool, packagePool, iocontext, iErr)
  end subroutine ocn_init_validate_act_test

end module mpas_ocn_init_act_test
