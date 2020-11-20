MODULE mckpp_xios_control

! Control routines to be called from main KPP code. 

USE mpi 
USE xios 
USE mckpp_xios_io

IMPLICIT NONE 


CONTAINS 

! Initialization of XIOS and diagnostic output 
! - equivalent to mckpp_initialize_output.
SUBROUTINE mckpp_xios_initialize_output(kpp_3d_fields,kpp_const_fields)

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  CALL mckpp_xios_initialize()
  CALL mckpp_xios_diagnostic_definition(kpp_3d_fields, kpp_const_fields)

END SUBROUTINE 


! Control of XIOS diagnostic and restart output 
! - equivalent to mckpp_output_control.
SUBROUTINE mckpp_xios_output_control(kpp_3d_fields, kpp_const_fields) 

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  REAL :: restart_time
  CHARACTER(LEN=9) :: restart_time_str
  CHARACTER(LEN=20) :: restart_filename

  ! Send diags to XIOS at every ts 
  CALL mckpp_xios_diagnostic_output(kpp_3d_fields, kpp_const_fields) 

  ! Check if restart timestep 
  ! - always write restart at end of run 
  IF (MOD(kpp_const_fields%ntime,kpp_const_fields%ndt_per_restart) .EQ. 0 .OR. &
    kpp_const_fields%ntime .EQ. kpp_const_fields%nend*kpp_const_fields%ndtocn) THEN

    ! Set correct time for validity of restart file (end of this timestep = start of next timestep)
    restart_time=kpp_const_fields%time+kpp_const_fields%dto/kpp_const_fields%spd

    WRITE(6,*) "Writing restart at time ", restart_time 

    WRITE(restart_time_str,'(F9.3)') restart_time
    restart_filename = "restart_"//TRIM(ADJUSTL(restart_time_str))//".nc" 
   
    CALL mckpp_xios_write_restart(kpp_3d_fields, kpp_const_fields, restart_filename) 
  END IF 

END SUBROUTINE mckpp_xios_output_control


SUBROUTINE mckpp_xios_write_restart(kpp_3d_fields, kpp_const_fields, restart_filename) 

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  CHARACTER(*) :: restart_filename

  CALL mckpp_xios_restart_definition(kpp_3d_fields, kpp_const_fields, restart_filename) 
  CALL mckpp_xios_restart_output(kpp_3d_fields, kpp_const_fields) 

  ! Now close context and switch back to diagnostic one 
  CALL xios_context_finalize()
  CALL xios_set_current_context(ctx_hdl_kpp)

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
