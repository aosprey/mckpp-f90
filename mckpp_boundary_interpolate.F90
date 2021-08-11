#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif

MODULE mckpp_boundary_interpolate

CONTAINS

SUBROUTINE MCKPP_BOUNDARY_INTERPOLATE_TEMP()
  
#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_3d_fields,kpp_const_fields
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid,only: get_ncols_p
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
#endif
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: npts, nzp1

  IMPLICIT NONE

  INTEGER prev_time,next_time,true_time
  REAL prev_weight,next_weight,ndays_upd_ocnT
  CHARACTER(LEN=31) :: routine = "MCKPP_BOUNDARY_INTERPOLATE_TEMP"
  CHARACTER(LEN=max_message_len) :: message

#ifdef MCKPP_CAM3
  INTEGER :: icol,ncol,ichnk
  REAL, allocatable :: prev_ocnT(:,:,:),next_ocnT(:,:)
  allocate(prev_ocnT(PCOLS,begchunk:endchunk,NZP1))
  allocate(next_ocnT(PCOLS,NZP1))
#else
  REAL, allocatable :: prev_ocnT(:,:),next_ocnT(:,:)  
  allocate(prev_ocnT(NPTS,NZP1))
  allocate(next_ocnT(NPTS,NZP1))
#endif

  true_time=kpp_const_fields%time
  ndays_upd_ocnT=kpp_const_fields%ndtupdocnT*kpp_const_fields%dto/kpp_const_fields%spd
      
  ! Read ocean temperatures for previous time
  prev_time=FLOOR((true_time+ndays_upd_ocnT/2)/ndays_upd_ocnT)*ndays_upd_ocnT-ndays_upd_ocnT*0.5
  IF (prev_time .lt. 0) THEN
     prev_weight=(ndays_upd_ocnT-ABS(true_time-prev_time))/ndays_upd_ocnT
     prev_time=prev_time+kpp_const_fields%ocnT_period
  ELSE
     prev_weight=(ndays_upd_ocnT-(true_time-prev_time))/ndays_upd_ocnT
   ENDIF

  WRITE(message,*) "true_time = ", true_time
  CALL mckpp_print(routine, message)
  WRITE(message,*) "prev_time = ", prev_time
  CALL mckpp_print(routine, message)
  WRITE(message,*) "prev_weight = ", prev_weight
  CALL mckpp_print(routine, message)
  
  kpp_const_fields%time=prev_time
  CALL MCKPP_READ_TEMPERATURES_3D()
#ifdef MCKPP_CAM3 
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     prev_ocnT(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%ocnT_clim(1:ncol,1:NZP1)
  ENDDO
#else
  prev_ocnT=kpp_3d_fields%ocnT_clim
#endif

  ! Read ocean temperatures for next time
  next_time=prev_time+ndays_upd_ocnT
  next_weight=1-prev_weight

  WRITE(message,*) "next_time = ", next_time
  CALL mckpp_print(routine, message)
  WRITE(message,*) "next_weight = ", next_weight
  CALL mckpp_print(routine, message)

  kpp_const_fields%time=next_time
  CALL MCKPP_READ_TEMPERATURES_3D()
#ifdef MCKPP_CAM3 
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     next_ocnT(1:ncol,1:NZP1)=kpp_3d_fields(ichnk)%ocnT_clim(1:ncol,1:NZP1)
     kpp_3d_fields(ichnk)%ocnT_clim(1:ncol,1:NZP1)=&
          next_ocnT(1:ncol,1:NZP1)*next_weight+prev_ocnT(1:ncol,ichnk,1:NZP1)*prev_weight
  ENDDO
#else
  next_ocnT=kpp_3d_fields%ocnT_clim  
  kpp_3d_fields%ocnT_clim=next_ocnT*next_weight+prev_ocnT*prev_weight
#endif

  kpp_const_fields%time=true_time
  deallocate(prev_ocnT)
  deallocate(next_ocnT)

END SUBROUTINE MCKPP_BOUNDARY_INTERPOLATE_TEMP


SUBROUTINE MCKPP_BOUNDARY_INTERPOLATE_SAL()
  
#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_3d_fields,kpp_const_fields
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid,only: get_ncols_p
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
#endif
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: npts, nzp1

  IMPLICIT NONE

  INTEGER prev_time,next_time,true_time
  REAL prev_weight,next_weight,ndays_upd_sal
  CHARACTER(LEN=31) :: routine = "mckpp_boundary_interpolate_temp"
  CHARACTER(LEN=max_message_len) :: message

#ifdef MCKPP_CAM3
  INTEGER :: icol,ncol,ichnk
  REAL, allocatable :: prev_sal(:,:,:),next_sal(:,:)
  allocate(prev_sal(PCOLS,begchunk:endchunk,NZP1))
  allocate(next_sal(PCOLS,NZP1))
#else
  REAL, allocatable :: prev_sal(:,:),next_sal(:,:)
  allocate(prev_sal(NPTS,NZP1))
  allocate(next_sal(NPTS,NZP1))
#endif

  true_time=kpp_const_fields%time
  ndays_upd_sal=kpp_const_fields%ndtupdsal*kpp_const_fields%dto/kpp_const_fields%spd
  
  ! Read ocean salinity for previous time
  prev_time=FLOOR((true_time+ndays_upd_sal/2)/ndays_upd_sal)*ndays_upd_sal-ndays_upd_sal*0.5
  IF (prev_time .lt. 0) THEN
     prev_weight=(ndays_upd_sal-ABS(true_time-prev_time))/ndays_upd_sal
     prev_time=prev_time+kpp_const_fields%sal_period
  ELSE
     prev_weight=(ndays_upd_sal-(true_time-prev_time))/ndays_upd_sal
  ENDIF

  WRITE(message,*) "true_time = ", true_time
  CALL mckpp_print(routine, message)
  WRITE(message,*) "prev_time = ", prev_time
  CALL mckpp_print(routine, message)
  WRITE(message,*) "prev_weight = ", prev_weight
  CALL mckpp_print(routine, message)
  
  kpp_const_fields%time=prev_time
  CALL MCKPP_READ_SALINITY_3D()
#ifdef MCKPP_CAM3 
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     prev_sal(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%sal_clim(1:ncol,1:NZP1)
  ENDDO
#else
  prev_sal=kpp_3d_fields%sal_clim
#endif
  
  ! Read ocean salinity for next time
  next_time=prev_time+ndays_upd_sal
  next_weight=1-prev_weight

  WRITE(message,*) "next_time = ", next_time
  CALL mckpp_print(routine, message)
  WRITE(message,*) "next_weight = ", next_weight
  CALL mckpp_print(routine, message)

  kpp_const_fields%time=next_time
  CALL MCKPP_READ_SALINITY_3D()
#ifdef MCKPP_CAM3 
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     next_sal(1:ncol,1:NZP1)=kpp_3d_fields(ichnk)%sal_clim(1:ncol,1:NZP1)
     kpp_3d_fields(ichnk)%sal_clim(1:ncol,1:NZP1)=&
          next_sal(1:ncol,1:NZP1)*next_weight+prev_sal(1:ncol,ichnk,1:NZP1)*prev_weight
  ENDDO
#else
  next_sal=kpp_3d_fields%sal_clim  
  kpp_3d_fields%sal_clim=next_sal*next_weight+prev_sal*prev_weight
#endif
  
  kpp_const_fields%time=true_time
  deallocate(prev_sal)
  deallocate(next_sal)
  
END SUBROUTINE MCKPP_BOUNDARY_INTERPOLATE_SAL

END MODULE mckpp_boundary_interpolate
