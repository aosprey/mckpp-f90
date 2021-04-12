MODULE mckpp_timer

  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  
  IMPLICIT NONE 

  PUBLIC :: mckpp_initialize_timers, mckpp_start_timer, mckpp_stop_timer, &
      mckpp_print_timers, mckpp_define_new_timer

  PRIVATE

  INTEGER, PARAMETER :: & 
      max_timers=300, & 
      max_name_length=30

  TYPE timer_type 
    REAL(KIND=8) :: elapsed_time, start_time
    LOGICAL :: running
    CHARACTER(LEN=max_name_length) :: name
  END TYPE timer_type 
  
  TYPE(timer_type), DIMENSION(max_timers) :: timers

  INTEGER :: timers_allocated, index_timer, index_total

  CHARACTER(LEN=11) :: routine = "MCKPP_TIMER"
  CHARACTER(LEN=max_message_len) :: message

CONTAINS 

  SUBROUTINE mckpp_initialize_timers

    REAL(kind=8) time
    INTEGER :: i

    time = get_current_time()

    timers_allocated = 0
    DO i = 1, max_timers 
      timers(i)%start_time = 0.
      timers(i)%elapsed_time = 0.
      timers(i)%running = .FALSE. 
      timers(i)%name = 'None'
    END DO 

    ! Built in timers 
    CALL mckpp_define_new_timer('Timer', index_timer) 
    timers(index_timer)%running = .TRUE. 
    timers(index_timer)%start_time = time

    CALL mckpp_define_new_timer('Total', index_total) 
    timers(index_total)%running = .TRUE. 
    timers(index_total)%start_time = time

    time = get_current_time()
    timers(index_timer)%elapsed_time = timers(index_timer)%elapsed_time + & 
                                      (time - timers(index_timer)%start_time) 

  END SUBROUTINE mckpp_initialize_timers


  SUBROUTINE mckpp_start_timer(name) 

    CHARACTER(LEN=*), INTENT(IN) :: name

    REAL(kind=8) time
    INTEGER :: index 

    time = get_current_time()

    index = lookup_timer_index(name) 
    IF (index .EQ. -1) THEN 
      CALL mckpp_define_new_timer(name, index) 
    END IF 

    IF (index .NE. -1) THEN     
      IF (timers(index)%running .EQV. .TRUE.) THEN
        WRITE(message,*) "Trying to start timer ", TRIM(name), " when it is already running."
        CALL mckpp_print_error(routine, message) 
      ELSE 
        timers(index)%start_time = get_current_time()
        timers(index)%running = .TRUE.
      END IF 
    END IF 

    timers(index_timer)%elapsed_time = timers(index_timer)%elapsed_time + & 
                                      (get_current_time() - time) 

  END SUBROUTINE mckpp_start_timer


  SUBROUTINE mckpp_stop_timer(name) 

    CHARACTER(LEN=*), INTENT(IN) :: name
  
    REAL(kind=8) time
    INTEGER :: index 
 
    time = get_current_time()

    index = lookup_timer_index(name)
    WRITE(message,*) "Trying to stop timer ", TRIM(name), " when it is not running."
    
    IF (index .EQ. -1) THEN  
      CALL mckpp_print_error(routine, message) 
    ELSE IF (timers(index)%running .EQV. .FALSE.) THEN 
      CALL mckpp_print_error(routine, message) 
    ELSE 
      timers(index)%elapsed_time = timers(index)%elapsed_time + & 
                                   (time - timers(index)%start_time) 
      timers(index)%running = .FALSE.
    END IF 

    timers(index_timer)%elapsed_time = timers(index_timer)%elapsed_time + & 
                                      (get_current_time() - time) 
  END SUBROUTINE mckpp_stop_timer

  
  SUBROUTINE mckpp_print_timers 
    
    REAL(KIND=8) :: time 
    INTEGER :: i

    ! Total time 
    time = get_current_time()   
    timers(index_total)%elapsed_time = time - timers(index_total)%start_time

    ! Print statistics 
    CALL mckpp_print(routine, "**** KPP TIMER STATISTICS ****")
    WRITE(message,'(A10,20X,2X,A22)') 'Timer name', 'Elapsed_time(s)' 
    CALL mckpp_print("", message) 
    DO i = 1, timers_allocated
      WRITE(message,'(A30,2X,F11.3)') timers(i)%name, timers(i)%elapsed_time
      CALL mckpp_print("", message)     
    END DO 

  END SUBROUTINE mckpp_print_timers 


  SUBROUTINE mckpp_define_new_timer(name, index)

    CHARACTER(LEN=*), INTENT(IN) :: name  
    INTEGER, INTENT(OUT) :: index 

    index = -1 
    IF (timers_allocated .GE. max_timers) THEN 
      CALL mckpp_print_error(routine, "Reached maximum number of timers") 
      RETURN
    END IF 

    index = timers_allocated + 1 
    timers_allocated = index
    timers(index)%name=name

  END SUBROUTINE mckpp_define_new_timer 


  INTEGER FUNCTION lookup_timer_index(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: i

    lookup_timer_index = -1

    ! Check length of timer name
    IF (LEN(name) .GT. max_name_length) THEN
      WRITE(message,*) 'KPP TIMER : Name of timer must not exceed ', max_name_length, ' characters: ', name
      CALL mckpp_print_error(routine, message)
      RETURN
    ENDIF
   
    ! Attempt to match name to existing timer 
    DO i = 1, timers_allocated 
      IF (name .EQ. timers(i)%name) THEN 
        lookup_timer_index = i
        RETURN 
      END IF 
    END DO 

  END FUNCTION lookup_timer_index

 
  REAL(KIND=8) FUNCTION get_current_time() 

    REAL(KIND=8) :: omp_get_wtime

#ifdef OPENMP
    get_current_time = omp_get_wtime()
#else
    CALL cpu_time(get_current_time)
#endif         

  END FUNCTION get_current_time


END MODULE mckpp_timer
