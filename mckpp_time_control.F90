MODULE mckpp_time_control

  ! Routines to work with timestep counter, current time, and
  ! deriving validity times for I/O routines.
  ! Time variables all stored in kpp_const_fields 

  USE mckpp_data_fields, ONLY: kpp_const_fields
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  
IMPLICIT NONE

PUBLIC
  
CONTAINS

  ! Update ntime and time in kpp_const_fields
  SUBROUTINE mckpp_update_time(ntime)

    INTEGER, INTENT(IN) :: ntime ! new timestep

    kpp_const_fields%ntime = ntime
    kpp_const_fields%time = mckpp_get_time(ntime)

  END SUBROUTINE mckpp_update_time
  

  ! Return time for a given timestep
  REAL FUNCTION mckpp_get_time(ntime)

    INTEGER, INTENT(IN) :: ntime ! timestep

    mckpp_get_time = kpp_const_fields%startt + (ntime-1) * &
       kpp_const_fields%dto / kpp_const_fields%spd
  
  END FUNCTION mckpp_get_time


  ! Work out time for ancil data reads, based on current model time, udpate freq,
  ! and time dim from file. 
  ! Return time to read in, and position in time dimension in file. 
  SUBROUTINE mckpp_get_update_time(file, time, nupdate, file_times, num_times, periodic, period, &
      update_time, update_time_pos)

    CHARACTER(LEN=*) :: file
    REAL, INTENT(IN) :: time 
    INTEGER, INTENT(IN) :: nupdate, num_times, period
    REAL*4, DIMENSION(num_times), INTENT(IN) :: file_times
    LOGICAL, INTENT(IN) :: periodic
    REAL*4, INTENT(OUT) :: update_time
    INTEGER, INTENT(OUT) :: update_time_pos

    REAL :: first_file_time, last_file_time, file_update_time
    CHARACTER(LEN=21) :: routine = "MCKPP_GET_UPDATE_TIME"
    CHARACTER(LEN=max_message_len) :: message

    first_file_time = file_times(1)
    last_file_time = file_times(num_times)

    update_time = time + 0.5 * (kpp_const_fields%dto / kpp_const_fields%spd) * nupdate
    update_time_pos = NINT( (update_time-first_file_time)*kpp_const_fields%spd / &
        (kpp_const_fields%dto*nupdate) ) + 1

    ! Check if ancil is periodic 
    IF (update_time .GT. last_file_time) THEN
      IF (periodic) THEN
         DO WHILE (update_time .GT. last_file_time) 
           update_time = update_time - period
        END DO
      ELSE
        WRITE(message,*) "Time to read exceeds the last time in the netCDF file", &
            " and periodic reads have not been specified."
        CALL mckpp_print_error(routine, message)
        WRITE(message,*) update_time, TRIM(file)
        CALL mckpp_print_error(routine, message)
        CALL mckpp_abort()
      END IF
    END IF
    
    ! Check this matches time entry in file 
    file_update_time = file_times(update_time_pos)
    IF ( ABS(file_update_time - update_time) .GT. & 
        (0.01 * kpp_const_fields%dtsec / kpp_const_fields%spd) ) THEN
      WRITE(message,*) 'Cannot find time,', update_time, 'in ', TRIM(file)
      CALL mckpp_print_error(routine, message) 
      WRITE(message,*) 'The closest I came was', file_update_time
      CALL mckpp_print_error(routine, message) 
      CALL mckpp_abort()
    END IF

  END SUBROUTINE mckpp_get_update_time 


END MODULE mckpp_time_control
