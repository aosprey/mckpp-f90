MODULE mckpp_initialize_fields_mod

  USE mckpp_boundary_interpolate, ONLY: mckpp_boundary_interpolate_temp, &
      mckpp_boundary_interpolate_sal
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields, &
      mckpp_allocate_3d_fields
  USE mckpp_initialize_advection_mod, ONLY: mckpp_initialize_advection
  USE mckpp_initialize_couplingweight_mod, ONLY: mckpp_initialize_couplingweight
  USE mckpp_initialize_fluxes, ONLY: mckpp_initialize_fluxes_variables
  USE mckpp_initialize_geography_mod, ONLY: mckpp_initialize_geography  
  USE mckpp_initialize_landsea_mod, ONLY: mckpp_initialize_landsea
  USE mckpp_initialize_ocean, ONLY: mckpp_initialize_ocean_model
  USE mckpp_initialize_ocean_profiles_mod, ONLY: mckpp_initialize_ocean_profiles
  USE mckpp_initialize_optics_mod, ONLY: mckpp_initialize_optics
  USE mckpp_initialize_relaxtion_mod, ONLY: mckpp_initialize_relaxation
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_warning, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, npts
  USE mckpp_physics_lookup_mod, ONLY: mckpp_physics_lookup
  USE mckpp_read_fluxes_mod, ONLY: mckpp_initialize_fluxes_file
  USE mckpp_read_heat_corrections_mod, ONLY: mckpp_read_fcorr_2d, mckpp_read_fcorr_3d
  USE mckpp_read_ice_mod, ONLY: mckpp_read_ice
  USE mckpp_read_salinity_mod, ONLY: mckpp_read_salinity_3d
  USE mckpp_read_salt_corrections_mod, ONLY: mckpp_read_sfcorr_2d, mckpp_read_sfcorr_3d
  USE mckpp_read_sst_mod, ONLY: mckpp_read_sst
  USE mckpp_read_temperatures_3d_mod, ONLY: mckpp_read_temperatures_3d
  USE mckpp_read_temperatures_bottom_mod, ONLY: mckpp_read_temperatures_bottom
  USE mckpp_xios_control, ONLY: mckpp_read_restart

  IMPLICIT NONE

CONTAINS

SUBROUTINE MCKPP_INITIALIZE_FIELDS()
  
  INTEGER :: iy,ix,ipt
  CHARACTER(LEN=23) :: routine = "MCKPP_INITIALIZE_FIELDS"
  CHARACTER(LEN=max_message_len) :: message 

  CALL mckpp_allocate_3d_fields()

  ! Set initial values for flags in kpp_3d_fields, which otherwise
  ! might never be set if points are not coupled.
  kpp_3d_fields%dampu_flag(:)=0.
  kpp_3d_fields%dampv_flag(:)=0.
  kpp_3d_fields%freeze_flag(:)=0.
  kpp_3d_fields%reset_flag(:)=0.
  kpp_3d_fields%cplwght(:)=0.
  
  ! Initialize latitude and longitude areas and the land/sea mask
  CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_LANDSEA")
  CALL MCKPP_INITIALIZE_LANDSEA()

  ! Initialize the vertical grid
  CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_GEOGRAPHY")
  CALL MCKPP_INITIALIZE_GEOGRAPHY()

  ! Initialize coupling weights   
  IF (kpp_const_fields%L_COUPLE .OR. kpp_const_fields%L_CPLWGHT) THEN 
    CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_COUPLINGWEIGHT") 
    CALL MCKPP_INITIALIZE_COUPLINGWEIGHT()
  ENDIF
    
  ! Initialize advection options
  CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_ADVECTION")
  CALL MCKPP_INITIALIZE_ADVECTION()

  ! Initialize relaxation of SST, temperature and/or salinity
  IF (kpp_const_fields%L_RELAX_SST .OR. kpp_const_fields%L_RELAX_SAL &
       .OR. kpp_const_fields%L_RELAX_OCNT) THEN     
     CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_RELAXATION")
     CALL MCKPP_INITIALIZE_RELAXATION()
  ENDIF

  ! Initialize water type for optical properties of seawater
  CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_OPTICS")
  CALL MCKPP_INITIALIZE_OPTICS()

  ! Initialize ocean profiles
  IF (kpp_const_fields%L_RESTART) THEN     
     CALL mckpp_print(routine, "Calling MCKPP_READ_RESTART")
     ! Still needs scattering code
     CALL mckpp_read_restart()
  ELSE
     CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_OCEAN_PROFILES")
     CALL MCKPP_INITIALIZE_OCEAN_PROFILES()
  ENDIF
  
  ! Initialize boundary conditions
  IF (kpp_const_fields%L_CLIMSST) THEN     
     CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_SST")
     CALL mckpp_read_sst()
  ENDIF

  IF (kpp_const_fields%L_CLIMICE) THEN    
     CALL mckpp_print(routine, "Calling MCKPP_READ_ICE")
     CALL MCKPP_READ_ICE()
  ENDIF
  
  IF (kpp_const_fields%L_FCORR_WITHZ) THEN     
     CALL mckpp_print(routine, "Calling MCKPP_READ_FCORR_3D")
     CALL MCKPP_READ_FCORR_3D()
  ELSEIF (kpp_const_fields%L_FCORR) THEN      
     CALL mckpp_print(routine, "Calling MCKPP_READ_FCORR_2D")
     CALL MCKPP_READ_FCORR_2D()
  ENDIF

  IF (kpp_const_fields%L_SFCORR_WITHZ) THEN   
     CALL mckpp_print(routine, "Calling MCKPP_READ_SFCORR_3D")
     CALL MCKPP_READ_SFCORR_3D()
  ELSEIF (kpp_const_fields%L_SFCORR) THEN
     CALL mckpp_print(routine, "Calling MCKPP_READ_SFCORR_2D")
     CALL MCKPP_READ_SFCORR_2D()
  ENDIF

  IF (kpp_const_fields%L_VARY_BOTTOM_TEMP) THEN     
     CALL mckpp_print(routine, "Calling MCKPP_READ_TEMPERATURES_BOTTOM")
     CALL MCKPP_READ_TEMPERATURES_BOTTOM()
  ENDIF

  !!! THESE SHOULD NOT BE CALLED IF DOING L_INTERP_OCNT / L_INTERP_SAL !!!

  IF (kpp_const_fields%L_RELAX_OCNT .AND. .NOT. kpp_const_fields%L_INTERP_OCNT) THEN 
     CALL mckpp_print(routine, "Calling MCKPP_READ_TEMPERATURES_3D")
     CALL MCKPP_READ_TEMPERATURES_3D()
  ELSEIF (kpp_const_fields%L_RELAX_OCNT .AND. kpp_const_fields%L_INTERP_OCNT) THEN     
     CALL mckpp_print(routine, "Calling MCKPP_BOUNDARY_INTERPOLATE_TEMP")
     CALL MCKPP_BOUNDARY_INTERPOLATE_TEMP()
  ENDIF

  IF (kpp_const_fields%L_RELAX_SAL .AND. .NOT. kpp_const_fields%L_INTERP_SAL) THEN
     CALL mckpp_print(routine, "Calling MCKPP_READ_SALINITY_3D")
     CALL MCKPP_READ_SALINITY_3D()
  ELSEIF (kpp_const_fields%L_RELAX_SAL .AND. kpp_const_fields%L_INTERP_SAL) THEN
     CALL mckpp_print(routine, "Calling MCKPP_BOUNDARY_INTERPOLATE_SAL")
     CALL MCKPP_BOUNDARY_INTERPOLATE_SAL()
  ENDIF

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
  
  ! Determine whether to run physics code on this column
  IF (kpp_const_fields%L_COUPLE) THEN
    kpp_3d_fields%run_physics = kpp_3d_fields%l_ocean .AND. kpp_3d_fields%cplwght(ipt) .GT. 0.0
  ELSE 
    kpp_3d_fields%run_physics = kpp_3d_fields%l_ocean 
  END IF 

  CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_OCEAN_MODEL")
  CALL MCKPP_INITIALIZE_OCEAN_MODEL()

END SUBROUTINE MCKPP_INITIALIZE_FIELDS

END MODULE mckpp_initialize_fields_mod
