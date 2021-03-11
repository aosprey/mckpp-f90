PROGRAM mckpp_ocean_model_3d

USE mckpp_data_types
USE mckpp_timer
USE mckpp_xios_control

IMPLICIT NONE

  TYPE(kpp_const_type) :: kpp_const_fields
  TYPE(kpp_3d_type) :: kpp_3d_fields

  INTEGER :: ntime
  
  ! Initialise
  WRITE(6,*) "MCKPP_OCEAN_MODEL_3D: Initialisation"

  CALL mckpp_initialize_timers()
  CALL mckpp_start_timer('Initialization')

  CALL mckpp_initialize_namelist(kpp_const_fields)
  CALL mckpp_initialize_fields(kpp_3d_fields,kpp_const_fields)
  CALL mckpp_initialize_output(kpp_3d_fields,kpp_const_fields)

  CALL mckpp_stop_timer('Initialization')

  ! Main time-stepping loop
  WRITE(6,*) "MCKPP_OCEAN_MODEL_3D: Timestepping loop"

  DO ntime=1,kpp_const_fields%nend*kpp_const_fields%ndtocn
     kpp_const_fields%ntime=ntime
     kpp_const_fields%time=kpp_const_fields%startt+(kpp_const_fields%ntime-1)*&
       kpp_const_fields%dto/kpp_const_fields%spd
     WRITE(6,*) "ntime, time = ", kpp_const_fields%ntime, kpp_const_fields%time

     ! Fluxes
     IF (MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndtocn) .EQ. 0) THEN
       CALL mckpp_start_timer("Update surface fluxes")
       CALL mckpp_fluxes(kpp_3d_fields,kpp_const_fields)
       CALL mckpp_stop_timer("Update surface fluxes")
     ENDIF

     ! Update boundary conditions
     IF (kpp_const_fields%ntime .ne. 1) THEN 
       CALL mckpp_start_timer("Update ancillaries")
       CALL mckpp_boundary_update(kpp_3d_fields,kpp_const_fields)
       CALL mckpp_stop_timer("Update ancillaries")
     ENDIF 
       
     ! Physics
     CALL mckpp_physics_driver(kpp_3d_fields,kpp_const_fields)

     ! Diagnostic output 
     CALL mckpp_start_timer("Diagnostic output") 
     CALL mckpp_output_control(kpp_3d_fields,kpp_const_fields)
     CALL mckpp_stop_timer("Diagnostic output")

     ! Restart output
     CALL mckpp_start_timer("Restart output")     
     CALL mckpp_restart_control(kpp_3d_fields,kpp_const_fields)
     CALL mckpp_stop_timer("Restart output")
     
  END DO

  ! Finalise
  WRITE(6,*) "MCKPP_OCEAN_MODEL_3D: Finalisation"
  CALL mckpp_xios_finalize() 
  CALL mckpp_print_timers()

END PROGRAM mckpp_ocean_model_3d
