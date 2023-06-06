PROGRAM mckpp_ocean_model_3d

  USE mckpp_boundary_update_mod, ONLY: mckpp_boundary_update
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_fluxes_mod, ONLY: mckpp_fluxes
  USE mckpp_initialize_fields_mod, ONLY: mckpp_initialize_fields
  USE mckpp_initialize_namelist_mod, ONLY: mckpp_initialize_namelist
  USE mckpp_physics_driver_mod, ONLY: mckpp_physics_driver
  USE mckpp_time_control, ONLY: mckpp_update_time
  USE mckpp_timer, ONLY: mckpp_initialize_timers, mckpp_start_timer, &
      mckpp_stop_timer, mckpp_print_timers
  USE mckpp_xios_control, ONLY: mckpp_initialize_xios, mckpp_initialize_output, &
      mckpp_output_control, mckpp_restart_control, mckpp_finalize_xios

  IMPLICIT NONE

  INTEGER :: ntime
  CHARACTER(LEN=20) :: routine = "MCKPP_OCEAN_MODEL_3D"
  CHARACTER(LEN=max_message_len) :: message

  ! Initialize
  CALL mckpp_print(routine, "Initialization")

  CALL mckpp_initialize_timers()
  CALL mckpp_start_timer('Initialization')

  kpp_const_fields%ntime = 0 

  CALL mckpp_initialize_xios() 
  CALL mckpp_initialize_namelist()
  CALL mckpp_initialize_fields()
  CALL mckpp_initialize_output()

  CALL mckpp_stop_timer('Initialization')

  ! Main time-stepping loop
  CALL mckpp_print(routine, "Timestepping loop")

  DO ntime = 1, kpp_const_fields%num_timesteps

    ! Update time variables
    CALL mckpp_update_time(ntime) 
    WRITE(message,*) "ntime, time = ", kpp_const_fields%ntime, &
                     kpp_const_fields%time
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

  ! Finalize
  CALL mckpp_print(routine, "Finalization")
  CALL mckpp_finalize_xios() 
  CALL mckpp_print_timers()

END PROGRAM mckpp_ocean_model_3d
