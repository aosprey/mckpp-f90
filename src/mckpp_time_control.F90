MODULE mckpp_time_control

  ! Routines to work with timestep counter, current time, and
  ! deriving validity times for I/O routines.
  ! Time variables all stored in kpp_const_fields 

  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_data_fields, ONLY: kpp_const_fields
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len

  IMPLICIT NONE

  INTEGER :: ntime 
  REAL :: time 

  PUBLIC time, ntime 
  PUBLIC mckpp_initialize_time, mckpp_update_time, mckpp_get_time, & 
         mckpp_get_update_time

  PRIVATE 

CONTAINS

  ! Initialize time variables at start of run 
  SUBROUTINE mckpp_initialize_time()

    CHARACTER(LEN=max_message_len) :: message
    CHARACTER(LEN=21) :: routine = "MCKPP_INITIALIZE_TIME"

    ntime = 0 
    time = kpp_const_fields%startt

    WRITE(message,*) "ntime, time = ", ntime, time
    CALL mckpp_print(routine, message)

  END SUBROUTINE mckpp_initialize_time 


  ! Update time based on ntime 
  SUBROUTINE mckpp_update_time(nt)

    INTEGER, INTENT(IN) :: nt
    CHARACTER(LEN=max_message_len) :: message
    CHARACTER(LEN=17) :: routine = "MCKPP_UPDATE_TIME"

    ntime = nt
    time = mckpp_get_time(nt)

    WRITE(message,*) "ntime, time = ", ntime, time
    CALL mckpp_print(routine, message)

  END SUBROUTINE mckpp_update_time


  ! Return time for a given timestep
  REAL FUNCTION mckpp_get_time(ntime)

    INTEGER, INTENT(IN) :: ntime ! timestep

    mckpp_get_time = kpp_const_fields%startt + (ntime-1) * &
        kpp_const_fields%dto / kpp_const_fields%spd

  END FUNCTION mckpp_get_time


  ! Work out time for ancil data reads, based on current model time, 
  ! udpate freq, and time dim from file.
  ! Return time to read in, and position in time dimension in file.
  ! Work out time based on method 1 or 2 - see get_update_time_[12] functions 
  ! below. 
  SUBROUTINE mckpp_get_update_time(file, time, ndt_update, file_times, & 
    num_times, periodic, period, update_time, update_time_pos, method)

    CHARACTER(LEN=*) :: file
    REAL, INTENT(IN) :: time 
    INTEGER, INTENT(IN) :: ndt_update, num_times, period, method
    REAL, DIMENSION(num_times), INTENT(IN) :: file_times
    LOGICAL, INTENT(IN) :: periodic
    REAL, INTENT(OUT) :: update_time
    INTEGER, INTENT(OUT) :: update_time_pos

    REAL :: first_file_time, last_file_time, file_update_time
    CHARACTER(LEN=21) :: routine = "MCKPP_GET_UPDATE_TIME"
    CHARACTER(LEN=max_message_len) :: message

    first_file_time = file_times(1)
    last_file_time = file_times(num_times)

    IF (method .EQ. 2) THEN 
      update_time = get_update_time_2(time, ndt_update) 
    ELSE 
      update_time = get_update_time_1(time, ndt_update)
    END IF 
      
    update_time_pos = get_update_pos(time, ndt_update, first_file_time) 

    ! Check if ancil is periodic 
    IF (update_time .GT. last_file_time) THEN
      IF (periodic) THEN
        DO WHILE (update_time .GT. last_file_time) 
          update_time = update_time - period
        END DO
      ELSE
        WRITE(message,*) "Time to read exceeds the last time in the netCDF ", & 
                         "file and periodic reads have not been specified."
        CALL mckpp_print_error(routine, message)
        WRITE(message,*) update_time, TRIM(file)
        CALL mckpp_print_error(routine, message)
        WRITE(message,*) "Error reading ", TRIM(file)
        CALL mckpp_abort(routine, message)
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
      WRITE(message,*) "Error reading ", TRIM(file)
      CALL mckpp_abort(routine, message)
    END IF

  END SUBROUTINE mckpp_get_update_time


  ! Calculate time to read data:
  ! Method 1: time + half update frequency 
  REAL FUNCTION get_update_time_1(time, ndt_update)

    REAL, INTENT(IN) :: time
    INTEGER, INTENT(IN) :: ndt_update 

    get_update_time_1 = time + 0.5 * (kpp_const_fields%dto / & 
                        kpp_const_fields%spd) * ndt_update

  END FUNCTION get_update_time_1


  ! Calculate time to read data:
  ! Method 2: ?
  REAL FUNCTION get_update_time_2(time, ndt_update)

    REAL, INTENT(IN) :: time
    INTEGER, INTENT(IN) :: ndt_update 
    REAL :: ndays_upd

    ndays_upd = ndt_update * kpp_const_fields%dto / kpp_const_fields%spd
    get_update_time_2 = ndays_upd * & 
           ( FLOOR(time,8) * NINT(kpp_const_fields%spd,8) / &
             ( ndt_update*NINT(kpp_const_fields%dto,8) ) ) + &
           ( 0.5 * kpp_const_fields%dto / kpp_const_fields%spd * ndt_update )

  END FUNCTION get_update_time_2


  ! Given a time value, work out position in time array 
  REAL FUNCTION get_update_pos(update_time, ndt_update, first_file_time)

    REAL, INTENT(IN) :: update_time, first_file_time
    INTEGER, INTENT(IN) :: ndt_update

    get_update_pos = NINT( (update_time - first_file_time) * & 
                     kpp_const_fields%spd / & 
                     ( kpp_const_fields%dto * ndt_update ) ) + 1

  END FUNCTION get_update_pos

END MODULE mckpp_time_control
