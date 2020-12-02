PROGRAM mckpp_ocean_model_3d

!  USE mckpp_3d_type_mod
USE mckpp_xios_control

IMPLICIT NONE

  TYPE(kpp_timer_type) :: kpp_timer
  TYPE(kpp_const_type) :: kpp_const_fields
  TYPE(kpp_3d_type), ALLOCATABLE :: kpp_3d_fields

  INTEGER :: ntime
  
  !CALL KPP_TIMER_INIT(kpp_timer)

  ! Initialise
  ! - could have a routine to call all of these
  WRITE(6,*) "MCKPP_OCEAN_MODEL_3D: Initialisation"
  ALLOCATE(kpp_3d_fields)
  CALL mckpp_initialize_namelist(kpp_const_fields)
  CALL mckpp_initialize_fields(kpp_3d_fields,kpp_const_fields)
  CALL mckpp_initialize_output(kpp_3d_fields,kpp_const_fields)

  ! Main time-stepping loop
  ! - again this could go in a timestep routine
  WRITE(6,*) "MCKPP_OCEAN_MODEL_3D: Timestepping loop"
  DO ntime=1,kpp_const_fields%nend*kpp_const_fields%ndtocn
     kpp_const_fields%ntime=ntime
     kpp_const_fields%time=kpp_const_fields%startt+(kpp_const_fields%ntime-1)*&
       kpp_const_fields%dto/kpp_const_fields%spd

     WRITE(6,*) "ntime, time = ", kpp_const_fields%ntime, kpp_const_fields%time

     ! Fluxes
     IF (MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndtocn) .EQ. 0) THEN
       CALL mckpp_fluxes(kpp_3d_fields,kpp_const_fields,kpp_timer)
     ENDIF

     ! Update boundary conditions
     IF (kpp_const_fields%ntime .ne. 1) THEN 
       CALL mckpp_boundary_update(kpp_3d_fields,kpp_const_fields,kpp_timer)
     ENDIF 
       
     ! Physics
     CALL mckpp_physics_driver(kpp_3d_fields,kpp_const_fields,kpp_timer)

     ! Output
     CALL mckpp_output_control(kpp_3d_fields,kpp_const_fields,kpp_timer)
     
  END DO

  ! Finalise
  WRITE(6,*) "MCKPP_OCEAN_MODEL_3D: Finalisation"

  CALL mckpp_xios_finalize() 

  ! Files are opened and closed as needed
  
  !CALL KPP_TIMER_PRINT(kpp_timer)

END PROGRAM mckpp_ocean_model_3d
