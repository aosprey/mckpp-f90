! Routines to print log messages.
! Messages including warnings and errors go to a log file for
! each MPI process.
! Error message from mckpp_abort goes to stderr.
MODULE mckpp_log_messages

  IMPLICIT NONE

  PUBLIC :: mckpp_initialize_logs, mckpp_finalize_logs, &
       mckpp_print, mckpp_print_error, mckpp_print_warning, &
       mckpp_print_stderr, update_context, max_message_len, &
       nupe, nuout, nuerr, l_debug
  
  PRIVATE 

  INTEGER :: my_rank 
  INTEGER, PARAMETER :: max_message_len = 250, &
      nuout = 6, nuerr = 0, nupe = 99
  LOGICAL :: l_debug = .TRUE.  ! This could go in a namelist

CONTAINS


  ! Create a log file for each MPI process
  ! Store rank in this module, to avoid circular dependencies
  ! with MPI module. 
  SUBROUTINE mckpp_initialize_logs(rank)

    INTEGER, INTENT(IN) :: rank 
    CHARACTER(LEN=4) :: rank_str
    CHARACTER(LEN=18) :: rank_filename

    my_rank = rank 
    WRITE(rank_str, '(I4.4)') my_rank
    rank_filename = "KPP_ocean_"//rank_str//".txt"
    OPEN(nupe, FILE=rank_filename, STATUS="replace")
    
  END SUBROUTINE mckpp_initialize_logs


  ! Close log file
  SUBROUTINE mckpp_finalize_logs

    CLOSE(nupe)

  END SUBROUTINE mckpp_finalize_logs
  

  ! Write general message
  SUBROUTINE mckpp_print(routine, message)

    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    CHARACTER(LEN=max_message_len) :: print_message

    IF (routine .NE. "") THEN 
      print_message = TRIM(routine) // ": " // TRIM(ADJUSTL(message))
    ELSE
      print_message = message
    END IF 
    CALL mckpp_write(nupe, print_message) 
    
  END SUBROUTINE mckpp_print


  ! Write error to log file, split over multiple lines. 
  SUBROUTINE mckpp_print_error(routine, message)

    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    CHARACTER(LEN=max_message_len) :: print_message
    
    IF (routine .NE. "") THEN 
      print_message = "Error in " // TRIM(ADJUSTL(routine)) // ":"
      CALL mckpp_write(nupe, print_message)
    ENDIF      
    CALL mckpp_write(nupe, message)
    
  END SUBROUTINE mckpp_print_error


  ! Write an error message to stderr (called by mckpp_abort) 
  SUBROUTINE mckpp_print_stderr(routine, message)
    
    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    CHARACTER(LEN=max_message_len) :: print_message, tmp

    WRITE(tmp,*) "Error from rank ", my_rank
    IF (routine .NE. "") THEN 
       print_message = tmp // " in " // TRIM(ADJUSTL(routine)) // ":"
      CALL mckpp_write(nuerr, print_message)
    ENDIF      
    CALL mckpp_write(nuerr, message)

  END SUBROUTINE mckpp_print_stderr


  ! Write warning, split over 2 lines 
  SUBROUTINE mckpp_print_warning(routine, message)
    
    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    CHARACTER(LEN=max_message_len) :: print_message

    IF (routine .NE. "") THEN 
      print_message = "Warning in " // TRIM(ADJUSTL(routine)) // ":"
      CALL mckpp_write(nupe, print_message)
    END IF 
    CALL mckpp_write(nupe, message)
    
  END SUBROUTINE mckpp_print_warning


  ! Sometimes it is useful to have the call path of routines for log messages.
  ! This function updates string with current routine
  FUNCTION update_context(calling_routine, routine) RESULT(context)

    CHARACTER(LEN=*), INTENT(IN) :: calling_routine, routine
    CHARACTER(LEN=max_message_len) :: context
    
    WRITE(context,*) TRIM(ADJUSTL(calling_routine)), " -> ", TRIM(ADJUSTL(routine))

  END FUNCTION update_context


  ! Internal : call write
  SUBROUTINE mckpp_write(unit, string)

    INTEGER, INTENT(IN) :: unit
    CHARACTER(LEN=*) :: string

    WRITE(unit,*) TRIM(ADJUSTL(string))

  END SUBROUTINE mckpp_write
  
END MODULE mckpp_log_messages
