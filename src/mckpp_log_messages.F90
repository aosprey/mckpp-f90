! Routines to print log messages to stdout, and error/warnings to 
! stderr.
! 
! Ideas:
! - split long messages over multiple lines
! - flags to control the level of prints
! - flag to flush output after each write
! - support to redirect messages to file

MODULE mckpp_log_messages

  IMPLICIT NONE 

  PUBLIC :: mckpp_print, mckpp_print_error, mckpp_print_warning, update_context
  PUBLIC :: max_message_len

  PRIVATE
  
  INTEGER :: nuout = 6,  nuerr = 0
  INTEGER, PARAMETER :: max_message_len = 200

CONTAINS

  ! Write to stdout
  SUBROUTINE mckpp_print(routine, message)

    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    CHARACTER(LEN=max_message_len) :: print_message
    
    print_message = TRIM(routine) // ": " // TRIM(ADJUSTL(message))
    CALL mckpp_write(nuout, print_message) 
    
  END SUBROUTINE mckpp_print


  ! Write error to stderr, split over 2 lines 
  SUBROUTINE mckpp_print_error(routine, message)
    
    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    CHARACTER(LEN=max_message_len) :: print_message

    IF (routine .NE. "") THEN 
      print_message = "Error in " // TRIM(ADJUSTL(routine)) // ":"
      CALL mckpp_write(nuerr, print_message)
    ENDIF 
      
    print_message = message
    CALL mckpp_write(nuerr, print_message)
    
  END SUBROUTINE mckpp_print_error


  ! Write warning to stderr, split over 2 lines 
  SUBROUTINE mckpp_print_warning(routine, message)
    
    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    CHARACTER(LEN=max_message_len) :: print_message

    print_message = "Warning in " // TRIM(ADJUSTL(routine)) // ":"
    CALL mckpp_write(nuerr, print_message)
    
    print_message = message
    CALL mckpp_write(nuerr, print_message)
    
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
