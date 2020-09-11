      SUBROUTINE kpp_timer_init(kpp_timer)
      IMPLICIT NONE

#include "kpp_3d_type.com"
#ifdef OPENMP
      REAL(kind=8) omp_get_wtime,current_time
#else
      REAL :: current_time
#endif
      TYPE(kpp_timer_type) :: kpp_timer
      INTEGER :: i
      
#ifdef OPENMP
      current_time=OMP_GET_WTIME()
#else
      CALL CPU_TIME(current_time)
#endif
      kpp_timer%timer_number_allocated=1
      kpp_timer%timer_elapsed_time(:)=0
      kpp_timer%timer_start_time(:)=0
      kpp_timer%timer_running(:)=.FALSE.
      kpp_timer%timer_start_time(timer_max_timers)=
     +     current_time
      kpp_timer%timer_running(timer_max_timers)=.TRUE.
      DO i=1,timer_max_timers-1
         kpp_timer%timer_all_names(i)='None'
      ENDDO      
      kpp_timer%timer_all_names(timer_max_timers)='Timer'

      kpp_timer%timer_start_time(kpp_timer%timer_number_allocated)=
     +     current_time
      kpp_timer%timer_running(kpp_timer%timer_number_allocated)=.TRUE.
      kpp_timer%timer_all_names(kpp_timer%timer_number_allocated)=
     +     'Total'
#ifdef OPENMP
      current_time=OMP_GET_WTIME()
#else
      CALL CPU_TIME(current_time)
#endif
      kpp_timer%timer_elapsed_time(timer_max_timers)=
     +     current_time-
     &     kpp_timer%timer_start_time(timer_max_timers)
!      WRITE(6,*) 'Total timer time: ',
!     &     timer_elapsed_time(timer_max_timers)
      
      END

      SUBROUTINE kpp_timer_time(kpp_timer,name,action)
      IMPLICIT NONE

      CHARACTER(LEN=*),intent(in) :: name
      INTEGER,intent(in) :: action

#include "kpp_3d_type.com"            
      INTEGER :: i,current_timer,flag
#ifdef OPENMP
      REAL omp_get_wtime,current_time
#else
      REAL :: current_time
#endif
      TYPE(kpp_timer_type) :: kpp_timer
#ifdef OPENMP
      current_time=OMP_GET_WTIME()
#else
      CALL CPU_TIME(current_time)
#endif     
      kpp_timer%timer_start_time(timer_max_timers)=
     +     current_time

!     Error checking on length of timer name
      IF (LEN(name) .gt. 30) THEN 
         WRITE(6,*) 'KPP TIMER : Name of timer must '//
     &        'not exceed 30 characters: '//name
         CALL MIXED_ABORT
      ENDIF
 
!     Attempt to match name to existing timer
     
      current_timer=0
      flag=0
      i=1
      DO WHILE (flag .eq. 0 .and. i .le. timer_max_timers)
!         WRITE(6,*) name,kpp_timer%timer_all_names(i)
         IF (name.eq.kpp_timer%timer_all_names(i)) THEN 
            current_timer=i
            flag=1
         ENDIF
         i=i+1
      ENDDO
      IF (i .eq. timer_max_timers) THEN
         WRITE(6,*) 'KPP TIMER: Reached maximum number of timers'
         RETURN
      ELSEIF (flag .eq. 0) THEN
         current_timer=kpp_timer%timer_number_allocated+1
         kpp_timer%timer_number_allocated=current_timer
         !WRITE(6,*) name
         kpp_timer%timer_all_names(current_timer)=name
      ENDIF
      
!     Check we have not hit the maximum number of timers.  This is
!     defined to be timer_max_timers-2 because we reserve the last timer
!     for timing the timer.  Meta.
      IF (current_timer.eq.timer_max_timers-2) THEN
         WRITE(6,*) 'KPP TIMER : Maximum number of timers '//
     &        'reached. The integration will continue, but '//
     &        'any timers after '//name//' will not be created.'
      ELSEIF (current_timer.eq.timer_max_timers-1) THEN
         RETURN
      ENDIF
      
      IF (action.eq.1) THEN
         IF (kpp_timer%timer_running(current_timer).eqv..TRUE.) THEN 
            WRITE(6,*) 'KPP TIMER : Trying to start timer '//
     &           name//' when it is already running. Stop.'
c            CALL MIXED_ABORT
         ENDIF
#ifdef OPENMP
         current_time=OMP_GET_WTIME()
#else
         CALL CPU_TIME(current_time)
#endif     
         kpp_timer%timer_start_time(current_timer)=current_time
         kpp_timer%timer_running(current_timer)=.TRUE.
c         WRITE(6,*) 'Started timer '//name,current_timer
      ENDIF

      IF (action.eq.0) THEN
         IF (kpp_timer%timer_running(current_timer).eqv..FALSE.) THEN
            WRITE(6,*) 'KPP TIMER : Trying to stop timer '//
     &           name//' number ',current_timer,
     +           ' when it is not running. Stop.'
c            CALL MIXED_ABORT
         ENDIF
#ifdef OPENMP
         current_time=OMP_GET_WTIME()
#else
         CALL CPU_TIME(current_time)
#endif
         kpp_timer%timer_elapsed_time(current_timer)=
     &        current_time-kpp_timer%timer_start_time(current_timer)+
     &        kpp_timer%timer_elapsed_time(current_timer)
         kpp_timer%timer_running(current_timer)=.FALSE.
c         WRITE(6,*) 'Stopped timer '//name,current_timer
      ENDIF      
      
#ifdef OPENMP
      current_time=OMP_GET_WTIME()
#else
      CALL CPU_TIME(current_time)
#endif
      kpp_timer%timer_elapsed_time(timer_max_timers)=
     +     current_time-
     &     kpp_timer%timer_start_time(timer_max_timers)+
     &     kpp_timer%timer_elapsed_time(timer_max_timers)
!      WRITE(6,*) 'Total timer time: ',
!     &     timer_elapsed_time(timer_max_timers)

      RETURN
      END
      
      SUBROUTINE kpp_timer_print(kpp_timer)
      IMPLICIT NONE

#include "kpp_3d_type.com"

      INTEGER :: i
#ifdef OPENMP
      REAL omp_get_wtime,current_time,total_time
#else
      REAL :: current_time,total_time
#endif
      TYPE(kpp_timer_type) :: kpp_timer
#ifdef OPENMP
      current_time=OMP_GET_WTIME()
#else
      CALL CPU_TIME(current_time)
#endif
      kpp_timer%timer_elapsed_time(1)=current_time-
     +     kpp_timer%timer_start_time(1)

      WRITE(6,*)
      WRITE(6,*) '**** KPP TIMER STATISTICS ****'
      WRITE(6,*)

      total_time=
     +     kpp_timer%timer_elapsed_time(timer_max_timers)
      DO i=2,kpp_timer%timer_number_allocated
         total_time=total_time+kpp_timer%timer_elapsed_time(i)
      ENDDO
      WRITE(6,'(A34,F11.3)') 'Total CPU time (from timer init): ',
     &     kpp_timer%timer_elapsed_time(1)
      WRITE(6,'(A31,3X,F11.3)') 'Total CPU time that was timed: ',
     &     total_time
      WRITE(6,*)
      WRITE(6,'(A10,20X,2X,A22)') 'Timer name','Elapsed time (seconds)'
      DO i=2,kpp_timer%timer_number_allocated
         WRITE(6,1000) kpp_timer%timer_all_names(i),
     +        kpp_timer%timer_elapsed_time(i)
      ENDDO
      WRITE(6,1000) 
     +     kpp_timer%timer_all_names(timer_max_timers),
     &     kpp_timer%timer_elapsed_time(timer_max_timers)

 1000 FORMAT(A30,2X,F11.3)
      
      END
      
