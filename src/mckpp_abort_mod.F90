MODULE mckpp_abort_mod
  
  IMPLICIT NONE

  USE mckpp_log_messages, ONLY: mckpp_print_stderr, mckpp_finalize_logs
  USE mpi 
  IMPLICIT NONE
  
CONTAINS 

  ! Print error message and shutdown model.
  ! Flush log file buffer. 
  ! The error code passed to MPI_abort doesn't mean anything.
  SUBROUTINE mckpp_abort(routine, message)

    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    INTEGER :: error_code = 101
    INTEGER :: ierr
    
    CALL mckpp_print_stderr(routine, message)
    CALL mckpp_finalize_logs()
    
    CALL mpi_abort(MPI_comm_world, error_code, ierr)
 
  END SUBROUTINE mckpp_abort
    
END MODULE mckpp_abort_mod
