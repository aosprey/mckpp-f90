MODULE mckpp_timer

  USE mckpp_parameters
  IMPLICIT NONE 

  PUBLIC :: mckpp_initialize_timers, mckpp_start_timer, mckpp_stop_timer, mckpp_print_timers 

  PRIVATE

  INTEGER, PARAMETER :: & 
      max_timers=100, & 
      max_name_length=30

  TYPE timer_type 
    REAL(KIND=8) :: elapsed_time, start_time
    LOGICAL :: running
    CHARACTER(LEN=max_name_length) :: name
  END TYPE timer_type 
  
  TYPE(timer_type), DIMENSION(max_timers) :: timers

  INTEGER :: timers_allocated, index_timer, index_total


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
    CALL define_new_timer('Timer', index_timer) 
    timers(index_timer)%running = .TRUE. 
    timers(index_timer)%start_time = time

    CALL define_new_timer('Total', index_total) 
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
      CALL define_new_timer(name, index) 
    END IF 

    IF (index .NE. -1) THEN     
      IF (timers(index)%running .EQV. .TRUE.) THEN       
        WRITE(nuerr,*) 'KPP TIMER : Trying to start timer ', name, ' when it is already running.'
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
    CHARACTER(LEN=max_error_msg_len) :: error_msg

    time = get_current_time()

    index = lookup_timer_index(name) 
    error_msg = 'KPP TIMER : Trying to stop timer '//name//' when it is not running.'
    IF (index .EQ. -1) THEN  
      WRITE(nuerr,*) error_msg
    ELSE IF (timers(index)%running .EQV. .FALSE.) THEN 
      WRITE(nuerr,*) error_msg
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
    WRITE(nuout,*) '**** KPP TIMER STATISTICS ****'
    WRITE(nuout,'(A10,20X,2X,A22)') 'Timer name', 'Elapsed_time(s)' 
    DO i = 1, timers_allocated
      WRITE(nuout,'(A30,2X,F11.3)') timers(i)%name, timers(i)%elapsed_time
    END DO 

  END SUBROUTINE mckpp_print_timers 


  SUBROUTINE define_new_timer(name, index)

    CHARACTER(LEN=*), INTENT(IN) :: name  
    INTEGER, INTENT(OUT) :: index 

    index = -1 
    IF (timers_allocated .GE. max_timers) THEN 
      WRITE(nuerr,*) 'KPP TIMER : Reached maximum number of timers'
      RETURN
    END IF 

    index = timers_allocated + 1 
    timers_allocated = index
    timers(index)%name=name

  END SUBROUTINE define_new_timer 


  INTEGER FUNCTION lookup_timer_index(name)

    CHARACTER(LEN=*), INTENT(IN) :: name
    INTEGER :: i

    lookup_timer_index = -1

    ! Check length of timer name
    IF (LEN(name) .GT. max_name_length) THEN 
      WRITE(nuerr,*) 'KPP TIMER : Name of timer must not exceed ', max_name_length, ' characters: ', name
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
