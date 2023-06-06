MODULE mckpp_initialize_fields_mod

  USE mckpp_boundary_interpolate, ONLY: mckpp_boundary_interpolate_temp, &
       mckpp_boundary_interpolate_sal
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields, &
       mckpp_allocate_3d_fields
  USE mckpp_initialize_advection_mod, ONLY: mckpp_initialize_advection
  USE mckpp_initialize_coupling_weight_mod, ONLY: & 
       mckpp_initialize_coupling_weight
  USE mckpp_initialize_fluxes, ONLY: mckpp_initialize_fluxes_variables
  USE mckpp_initialize_geography_mod, ONLY: mckpp_initialize_geography  
  USE mckpp_initialize_landsea_mod, ONLY: mckpp_initialize_landsea
  USE mckpp_initialize_ocean, ONLY: mckpp_initialize_ocean_model
  USE mckpp_initialize_ocean_profiles_mod, ONLY: mckpp_initialize_ocean_profiles
  USE mckpp_initialize_optics_mod, ONLY: mckpp_initialize_optics
  USE mckpp_initialize_relaxtion_mod, ONLY: mckpp_initialize_relaxation
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_warning, max_message_len
  USE mckpp_physics_lookup_mod, ONLY: mckpp_physics_lookup
  USE mckpp_read_fluxes_mod, ONLY: mckpp_initialize_fluxes_file
  USE mckpp_read_heat_corrections_mod, ONLY: mckpp_read_fcorr_2d, &
       mckpp_read_fcorr_3d
  USE mckpp_read_ice_mod, ONLY: mckpp_read_ice
  USE mckpp_read_salinity_mod, ONLY: mckpp_read_salinity_3d
  USE mckpp_read_salt_corrections_mod, ONLY: mckpp_read_sfcorr_2d, & 
       mckpp_read_sfcorr_3d
  USE mckpp_read_sst_mod, ONLY: mckpp_read_sst
  USE mckpp_read_temperatures_3d_mod, ONLY: mckpp_read_temperatures_3d
  USE mckpp_read_temperatures_bottom_mod, ONLY: mckpp_read_temperatures_bottom
  USE mckpp_xios_control, ONLY: mckpp_read_restart

  IMPLICIT NONE

CONTAINS

  SUBROUTINE mckpp_initialize_fields()

    CHARACTER(LEN=23) :: routine = "MCKPP_INITIALIZE_FIELDS"
    CHARACTER(LEN=max_message_len) :: message 

    CALL mckpp_allocate_3d_fields()

     ! Initialize latitude and longitude areas and the land/sea mask
    CALL mckpp_initialize_landsea()

    ! Initialize the vertical grid
    CALL mckpp_initialize_geography()

    ! Initialize coupling weights   
    CALL mckpp_initialize_coupling_weight()    
    CALL mckpp_initialize_flags()

    ! Initialize advection options
    CALL mckpp_initialize_advection()

    ! Initialize relaxation of SST, temperature and/or salinity
    IF ( kpp_const_fields%l_relax_sst .OR. kpp_const_fields%l_relax_sal & 
         .OR. kpp_const_fields%l_relax_ocnt ) &   
       CALL mckpp_initialize_relaxation()

    ! Initialize water type for optical properties of seawater
    CALL mckpp_initialize_optics()

    ! Initialize ocean profiles
    IF ( kpp_const_fields%l_restart ) THEN     
      CALL mckpp_read_restart()
    ELSE
      CALL mckpp_initialize_ocean_profiles()
    ENDIF
 
    ! Initialize boundary conditions
    IF ( kpp_const_fields%l_climsst ) CALL mckpp_read_sst()
    IF ( kpp_const_fields%l_climice ) CALL mckpp_read_ice()

    ! Heat corrections
    IF ( kpp_const_fields%l_fcorr_withz ) THEN     
      CALL mckpp_read_fcorr_3d()

    ELSEIF ( kpp_const_fields%l_fcorr ) THEN      
      CALL mckpp_read_fcorr_2d()
    ENDIF

    ! Salt corrections
    IF ( kpp_const_fields%l_sfcorr_withz ) THEN   
       CALL mckpp_read_sfcorr_3d()

    ELSEIF ( kpp_const_fields%l_sfcorr ) THEN
      CALL mckpp_read_sfcorr_2d()
    ENDIF

    ! Bottom temperature
    IF ( kpp_const_fields%l_vary_bottom_temp ) & 
      CALL mckpp_read_temperatures_bottom()
 
   ! 3D temperature relaxation
    IF ( kpp_const_fields%l_relax_ocnt ) THEN 
      IF ( .NOT. kpp_const_fields%l_interp_ocnt ) THEN 
        CALL mckpp_read_temperatures_3d()
      ELSE
        CALL mckpp_boundary_interpolate_temp()
      END IF
    END IF 

    ! 3D salinity relaxation
    IF ( kpp_const_fields%l_relax_sal ) THEN 
      IF ( .NOT. kpp_const_fields%l_interp_sal ) THEN
        CALL mckpp_read_salinity_3d()
      ELSE
        CALL mckpp_boundary_interpolate_sal()
      END IF
    END IF 

    CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_FLUXES_VARIABLES")
    CALL MCKPP_INITIALIZE_FLUXES_VARIABLES()

    ! Isothermal detection routine requires 3D ocean temperature and salinity fields
    IF (kpp_const_fields%L_NO_ISOTHERM .AND. .NOT. kpp_const_fields%L_RELAX_SAL &
         .AND. .NOT. kpp_const_fields%L_RELAX_OCNT) THEN         
       CALL mckpp_print(routine, " Calling MCKPP_READ_TEMPERATURES_3D")
       CALL MCKPP_READ_TEMPERATURES_3D()

       CALL mckpp_print(routine, "Calling MCKPP_READ_SALINITY_3D")
       CALL MCKPP_READ_SALINITY_3D()
    ENDIF

    ! L_INTERP_OCNT and L_INTERP_SAL imply L_PERIODIC_OCNT and L_PERIODIC_SAL,
    ! respectively, to deal with times before the first time in the input file.
    IF (kpp_const_fields%L_INTERP_OCNT) kpp_const_fields%L_PERIODIC_OCNT=.TRUE.
    IF (kpp_const_fields%L_INTERP_SAL) kpp_const_fields%L_PERIODIC_SAL=.TRUE.

    IF (kpp_const_fields%L_FLUXDATA .AND. .NOT. kpp_const_fields%L_COUPLE) THEN
       CALL MCKPP_INITIALIZE_FLUXES_FILE()       
    ENDIF

    kpp_const_fields%ntime=0
    CALL mckpp_print(routine, "Calling MCKPP_PHYSICS_LOOKUP")
    CALL MCKPP_PHYSICS_LOOKUP(kpp_const_fields)

    CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_OCEAN_MODEL")
    CALL MCKPP_INITIALIZE_OCEAN_MODEL()

  END SUBROUTINE MCKPP_INITIALIZE_FIELDS


  SUBROUTINE mckpp_initialize_flags()

    ! Set initial values for flags in kpp_3d_fields, which otherwise
    ! might never be set if points are not coupled.
    kpp_3d_fields%dampu_flag(:) = 0.0
    kpp_3d_fields%dampv_flag(:) = 0.0
    kpp_3d_fields%freeze_flag(:) = 0.0
    kpp_3d_fields%reset_flag(:) = 0.0

    ! Determine whether to run KPP physics on this pt
    IF (kpp_const_fields%l_couple) THEN
      kpp_3d_fields%run_physics = kpp_3d_fields%l_ocean .AND. &
                                  kpp_3d_fields%cplwght .GT. 0.0
    ELSE 
      kpp_3d_fields%run_physics = kpp_3d_fields%l_ocean 
    END IF
    
  END SUBROUTINE mckpp_initialize_flags


END MODULE mckpp_initialize_fields_mod
