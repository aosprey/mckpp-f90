! Routines to print log messages to stdout, and error/warnings to 
! stderr. If running in parallel (CAM) only masterproc writes messages.
! 
! Ideas:
! - split long messages over multiple lines
! - flags to control the level of prints
! - flag to flush output after each write
! - support to redirect messages to file

MODULE mckpp_log_messages

#ifdef MCKPP_CAM3
  USE pmgrid, only: masterproc
#endif  
 
  IMPLICIT NONE 

  PUBLIC :: mckpp_print, mckpp_print_error, mckpp_print_warning
  PUBLIC :: max_message_len

  PRIVATE
  
  INTEGER :: nuout = 6,  nuerr = 0
  INTEGER, PARAMETER :: max_message_len = 150, max_print_len = 200

CONTAINS

  ! Write to stdout
  SUBROUTINE mckpp_print(routine, message)

    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    CHARACTER(LEN=max_print_len) :: print_message
    
    print_message = TRIM(routine) // ": " // TRIM(ADJUSTL(message))
    CALL mckpp_write(nuout, print_message) 
    
  END SUBROUTINE mckpp_print


  ! Write error to stderr, split over 2 lines 
  SUBROUTINE mckpp_print_error(routine, message)
    
    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    CHARACTER(LEN=max_print_len) :: print_message

    print_message = "Error in " // TRIM(routine) // ":"
    CALL mckpp_write(nuerr, print_message)
    
    print_message = message
    CALL mckpp_write(nuerr, print_message)
    
  END SUBROUTINE mckpp_print_error


  ! Write warning to stderr, split over 2 lines 
  SUBROUTINE mckpp_print_warning(routine, message)
    
    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    CHARACTER(LEN=max_print_len) :: print_message

    print_message = "Warning in " // TRIM(routine) // ":"
    CALL mckpp_write(nuerr, print_message)
    
    print_message = message
    CALL mckpp_write(nuerr, print_message)
    
  END SUBROUTINE mckpp_print_warning
  

  ! Internal : call write
  SUBROUTINE mckpp_write(unit, string)

    INTEGER, INTENT(IN) :: unit
    CHARACTER(LEN=*) :: string

#ifdef MCKPP_CAM3
    IF (masterproc) THEN WRITE(unit,*) TRIM(ADJUSTL(string))
#else
    WRITE(unit,*) TRIM(ADJUSTL(string))
#endif

  END SUBROUTINE mckpp_write
  
END MODULE mckpp_log_messages 
