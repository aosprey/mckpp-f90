#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE mckpp_physics_driver
  USE mckpp_types, only: kpp_3d_fields,kpp_1d_type,kpp_const_fields
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid,only: get_ncols_p
#else
SUBROUTINE mckpp_physics_driver(kpp_3d_fields,kpp_const_fields)
  IMPLICIT NONE
#endif

#ifdef MCKPP_CAM3
#include <parameter.inc>
  INTEGER :: icol,ncol,ichnk
#else
  ! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
#endif
  
  ! Local
  TYPE(kpp_1d_type) :: kpp_1d_fields
  INTEGER, parameter :: nuout=6, nuerr=0
  INTEGER :: ipt
#ifdef OPENMP
  INTEGER :: tid,OMP_GET_THREAD_NUM
#endif
  CHARACTER(LEN=21) phys_timer_name
  CHARACTER(LEN=19) trans_timer_name
  
#ifdef MCKPP_CAM3
  !WRITE(6,*) 'Before ocnstep, U = ',kpp_3d_fields(begchunk)%U(1:ncol,1,1)
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     DO icol=1,ncol
        IF (kpp_3d_fields(ichnk)%landfrac(icol) .lt. 1.0 .and. &
             kpp_3d_fields(ichnk)%cplwght(icol) .gt. 0.0) THEN
           CALL MCKPP_FIELDS_3Dto1D(kpp_3d_fields(ichnk),icol,kpp_1d_fields)
           CALL MCKPP_PHYSICS_OCNSTEP(kpp_1d_fields,kpp_const_fields)
           CALL MCKPP_PHYSICS_OVERRIDES_CHECK_PROFILE(kpp_1d_fields,kpp_const_fields)
           CALL MCKPP_FIELDS_1Dto3D(kpp_1d_fields,icol,kpp_3d_fields(ichnk))
        ENDIF
     ENDDO
  ENDDO   
  !WRITE(6,*) 'After ocnstep, U = ',kpp_3d_fields(begchunk)%U(1:ncol,1,1)
#else
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(kpp_3d_fields,kpp_const_fields) &
!!$OMP SHARED(kpp_timer) & 
!$OMP PRIVATE(trans_timer_name,phys_timer_name,tid)
  tid=OMP_GET_THREAD_NUM()
  WRITE(trans_timer_name,'(A17,I2)') 'KPP 3D/1D thread ',tid
  WRITE(phys_timer_name,'(A19,I2)') 'KPP Physics thread ',tid
!$OMP DO SCHEDULE(dynamic)
#else
  WRITE(trans_timer_name,'(A19)') 'KPP 3D/1D thread 01'
  WRITE(phys_timer_name,'(A21)') 'KPP Physics thread 01'
#endif /*OPENMP*/
  DO ipt=1,npts
     IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
        !CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,1)
        CALL mckpp_fields_3dto1d(kpp_3d_fields,ipt,kpp_1d_fields)
        !CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,0)
        !CALL KPP_TIMER_TIME(kpp_timer,phys_timer_name,1)                
        CALL mckpp_physics_ocnstep(kpp_1d_fields,kpp_const_fields)       
        CALL mckpp_physics_overrides_check_profile(kpp_1d_fields,kpp_const_fields)        
        !CALL KPP_TIMER_TIME(kpp_timer,phys_timer_name,0)
        !CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,1)
        CALL mckpp_fields_1dto3d(kpp_1d_fields,ipt,kpp_3d_fields)
        !CALL KPP_TIMER_TIME(kpp_timer,trans_timer_name,0)
     ENDIF
  ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif /*OPENMP*/
#endif /*MCKPP_CAM3*/
  
  IF (kpp_const_fields%L_VARY_BOTTOM_TEMP) THEN
     !CALL KPP_TIMER_TIME(kpp_timer,'Top level',0)
     !CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',1)
#ifdef MCKPP_CAM3
     CALL MCKPP_PHYSICS_OVERRIDES_BOTTOMTEMP
#else
     CALL MCKPP_PHYSICS_OVERRIDES_BOTTOMTEMP(kpp_3d_fields,kpp_const_fields)
#endif
     !CALL KPP_TIMER_TIME(kpp_timer,'Update ancillaries',0)
     !CALL KPP_TIMER_TIME(kpp_timer,'Top level',1)
  ENDIF

  RETURN
END SUBROUTINE mckpp_physics_driver
