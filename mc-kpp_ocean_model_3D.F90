PROGRAM mckpp_ocean_model_3d

USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
USE mckpp_timer, ONLY: mckpp_initialize_timers, mckpp_start_timer, &
    mckpp_stop_timer, mckpp_print_timers
USE mckpp_xios_control, ONLY: mckpp_initialize_output, mckpp_output_control, &
    mckpp_restart_control, mckpp_finalize_output

IMPLICIT NONE

  INTEGER :: ntime
  CHARACTER(LEN=20) :: routine = "MCKPP_OCEAN_MODEL_3D"
  CHARACTER(LEN=max_message_len) :: message
  
  ! Initialise
  CALL mckpp_print(routine, "Initialisation")

  CALL mckpp_initialize_timers()
  CALL mckpp_start_timer('Initialization')

  CALL mckpp_initialize_namelist()
  CALL mckpp_initialize_fields()
  CALL mckpp_initialize_output()

  CALL mckpp_stop_timer('Initialization')

  ! Main time-stepping loop
  CALL mckpp_print(routine, "Timestepping loop")

  DO ntime = 1, kpp_const_fields%nend*kpp_const_fields%ndtocn
     kpp_const_fields%ntime = ntime
     kpp_const_fields%time = kpp_const_fields%startt+(kpp_const_fields%ntime-1)*&
       kpp_const_fields%dto/kpp_const_fields%spd
     WRITE(message,*) "ntime, time = ", kpp_const_fields%ntime, kpp_const_fields%time
     CALL mckpp_print(routine, message)

     ! Fluxes
     IF (MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndtocn) .EQ. 0) THEN
       CALL mckpp_start_timer("Update surface fluxes")
       CALL mckpp_fluxes()
       CALL mckpp_stop_timer("Update surface fluxes")
     ENDIF

     ! Update boundary conditions
     IF (kpp_const_fields%ntime .ne. 1) THEN 
       CALL mckpp_start_timer("Update ancillaries")
       CALL mckpp_boundary_update()
       CALL mckpp_stop_timer("Update ancillaries")
     ENDIF 
       
     ! Physics
     CALL mckpp_physics_driver()

     ! Diagnostic output 
     CALL mckpp_start_timer("Diagnostic output") 
     CALL mckpp_output_control()
     CALL mckpp_stop_timer("Diagnostic output")

     ! Restart output
     CALL mckpp_start_timer("Restart output")     
     CALL mckpp_restart_control()
     CALL mckpp_stop_timer("Restart output")
     
  END DO

  ! Finalise
  CALL mckpp_print(routine, "Finalisation")
  CALL mckpp_finalize_output() 
  CALL mckpp_print_timers()

END PROGRAM mckpp_ocean_model_3d
