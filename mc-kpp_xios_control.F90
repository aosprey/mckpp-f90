MODULE mckpp_xios_control

! Control routines to be called from main KPP code. 

USE mckpp_data_types
USE mckpp_timer
USE mckpp_xios_io

USE mpi 
USE xios 

IMPLICIT NONE 


CONTAINS 

! Initialization for diagnostic output 
SUBROUTINE mckpp_initialize_output(kpp_3d_fields,kpp_const_fields)

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  CALL mckpp_xios_initialize()
  CALL mckpp_xios_diagnostic_definition(kpp_3d_fields, kpp_const_fields)

END SUBROUTINE mckpp_initialize_output


! Diagnostic output
SUBROUTINE mckpp_output_control(kpp_3d_fields, kpp_const_fields) 

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  REAL :: restart_time
 
  ! Send diags to XIOS at every ts 
  CALL mckpp_xios_diagnostic_output(kpp_3d_fields, kpp_const_fields) 

END SUBROUTINE mckpp_output_control


! Restart output
SUBROUTINE mckpp_restart_control(kpp_3d_fields, kpp_const_fields) 

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  REAL :: restart_time
 
  ! Check if restart timestep 
  ! - always write restart at end of run 
  CALL mckpp_start_timer("Write restart output") 
  IF (MOD(kpp_const_fields%ntime,kpp_const_fields%ndt_per_restart) .EQ. 0 .OR. &
    kpp_const_fields%ntime .EQ. kpp_const_fields%nend*kpp_const_fields%ndtocn) THEN

    ! Set correct time for validity of restart file (end of this timestep = start of next timestep)
    restart_time=kpp_const_fields%time+kpp_const_fields%dto/kpp_const_fields%spd

    WRITE(6,*) "Writing restart at time ", restart_time    
    CALL mckpp_xios_write_restart(kpp_3d_fields, kpp_const_fields, restart_time) 
  END IF 
  CALL mckpp_stop_timer("Write restart output") 

END SUBROUTINE mckpp_restart_control


SUBROUTINE mckpp_xios_write_restart(kpp_3d_fields, kpp_const_fields, restart_time) 

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  REAL, INTENT(IN) :: restart_time

  CHARACTER(LEN=9) :: restart_time_str
  CHARACTER(LEN=17) :: restart_filename
  CHARACTER(LEN=21) :: context_name
  TYPE(xios_context) :: ctx_hdl_restart

  WRITE(restart_time_str,'(F9.3)') restart_time
  restart_filename = "restart_"//TRIM(ADJUSTL(restart_time_str))
  context_name = "ctx_restart_"//TRIM(ADJUSTL(restart_time_str))

  ! Define a new context for this restart file 
  CALL xios_context_initialize(context_name, xios_comm)
  CALL xios_get_handle(context_name, ctx_hdl_restart)
  CALL xios_set_current_context(ctx_hdl_restart)

  ! Define file  
  CALL mckpp_xios_restart_definition(kpp_3d_fields, kpp_const_fields, restart_filename) 

  ! Write fields 
  CALL mckpp_xios_restart_output(kpp_3d_fields, kpp_const_fields) 

  ! Now close context and switch back to diagnostic one 
  CALL xios_context_finalize()
  CALL xios_set_current_context(ctx_hdl_diags)

END SUBROUTINE


SUBROUTINE mckpp_xios_initialize()

  INTEGER :: ierr
  
  CALL mpi_init(ierr) ! This should go elsewhere but leave here for now. 
  CALL xios_initialize("client", return_comm=xios_comm)

END SUBROUTINE mckpp_xios_initialize


SUBROUTINE mckpp_xios_finalize() 

  INTEGER :: ierr

  CALL xios_context_finalize()
  CALL xios_finalize()

  CALL mpi_finalize(ierr)

END SUBROUTINE mckpp_xios_finalize


END MODULE mckpp_xios_control
