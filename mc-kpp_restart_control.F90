#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE mckpp_restart_control
  USE mckpp_types, only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE shr_kind_mod, only: r8=>shr_kind_r8, r4=>shr_kind_r4
#else
SUBROUTINE mckpp_restart_control(kpp_3d_fields,kpp_const_fields)
#endif
  IMPLICIT NONE

#ifdef MCKPP_CAM3
#include <parameter.inc>
  REAL(r8) :: restart_time
  ! Parameter to force writing of restart file, no matter the time.
  ! Added to enable writing MC-KPP restart file at same time as
  ! CAM restart file.
#else
#include <mc-kpp_3d_type.com>  
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  REAL :: restart_time
#endif

  CHARACTER*5 :: restart_time_char
  
  IF (kpp_const_fields%L_RESTARTW) THEN
     ! Disable regular frequency restarts for MC-KPP when coupled to CAM3, to
     ! allow us to write restart files at the same time as the atmosphere (typically the
     ! first day of every calendar month on 365-day calendar).  Write restart file
     ! whenever this routine is called.
#ifndef MCKPP_CAM3        
     ! Write restart file either every ndt_per_restart timesteps
     ! or at the end of the simulation.
     IF (MOD(kpp_const_fields%ntime,kpp_const_fields%ndt_per_restart) .EQ. 0 .OR. &
          kpp_const_fields%ntime .EQ. kpp_const_fields%nend*kpp_const_fields%ndtocn) THEN
#endif        
        ! Set correct time for validity of restart file (end of this timestep = start of next timestep)
        restart_time=kpp_const_fields%time+kpp_const_fields%dto/kpp_const_fields%spd
        !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
        !CALL KPP_TIMER_TIME(kpp_timer,'Writing Restart File',1)
        IF (restart_time .lt. 10) THEN
           WRITE(restart_time_char,'(A4,I1)') '0000',FLOOR(restart_time)
        ELSEIF (restart_time .lt. 100) THEN
           WRITE(restart_time_char,'(A3,I2)') '000',FLOOR(restart_time)
        ELSEIF (restart_time .lt. 1000) THEN
           WRITE(restart_time_char,'(A2,I3)') '00',FLOOR(restart_time)
        ELSEIF (restart_time .lt. 10000) THEN
           WRITE(restart_time_char,'(A1,I4)') '0',FLOOR(restart_time)
        ELSE
           WRITE(restart_time_char,'(I5)') FLOOR(restart_time)
        ENDIF
        WRITE(kpp_const_fields%restart_outfile,'(A12,A5)') 'KPP.restart.',restart_time_char
        !WRITE(6,*) 'MCKPP_RESTART_CONTROL: Calling MCKPP_RESTART_IO_WRITE_NETCDF'
#ifdef MCKPP_CAM3
        CALL MCKPP_RESTART_IO_WRITE_NETCDF
#else
        CALL MCKPP_RESTART_IO_WRITE_NETCDF(kpp_3d_fields,kpp_const_fields)
#endif
        !CALL MCKPP_RESTART_IO_WRITE(kpp_3d_fields,kpp_const_fields)
        !CALL KPP_TIMER_TIME(kpp_timer,'Writing Restart File',0)
        !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
#ifndef MCKPP_CAM3
     ENDIF
#endif
  ENDIF
  
  RETURN
END SUBROUTINE mckpp_restart_control
