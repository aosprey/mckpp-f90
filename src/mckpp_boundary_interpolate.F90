MODULE mckpp_boundary_interpolate

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_mpi_control, ONLY: npts_local
  USE mckpp_parameters, ONLY: nzp1
  USE mckpp_read_temperatures_3d_mod, ONLY: mckpp_read_temperatures_3d
  USE mckpp_read_salinity_mod, ONLY: mckpp_read_salinity_3d
  USE mckpp_time_control, ONLY: time

  IMPLICIT NONE
  
CONTAINS

SUBROUTINE MCKPP_BOUNDARY_INTERPOLATE_TEMP()
  
  INTEGER prev_time,next_time,true_time
  REAL prev_weight,next_weight,ndays_upd_ocnT
  CHARACTER(LEN=31) :: routine = "MCKPP_BOUNDARY_INTERPOLATE_TEMP"
  CHARACTER(LEN=max_message_len) :: message

  REAL, allocatable :: prev_ocnT(:,:),next_ocnT(:,:)  
  allocate(prev_ocnT(npts_local, nzp1))
  allocate(next_ocnT(npts_local, nzp1))

  true_time=time
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
  
  time=prev_time
  CALL MCKPP_READ_TEMPERATURES_3D()
  prev_ocnT=kpp_3d_fields%ocnT_clim

  ! Read ocean temperatures for next time
  next_time=prev_time+ndays_upd_ocnT
  next_weight=1-prev_weight

  WRITE(message,*) "next_time = ", next_time
  CALL mckpp_print(routine, message)
  WRITE(message,*) "next_weight = ", next_weight
  CALL mckpp_print(routine, message)

  time=next_time
  CALL MCKPP_READ_TEMPERATURES_3D()
  next_ocnT=kpp_3d_fields%ocnT_clim  
  kpp_3d_fields%ocnT_clim=next_ocnT*next_weight+prev_ocnT*prev_weight

  time=true_time
  deallocate(prev_ocnT)
  deallocate(next_ocnT)

END SUBROUTINE MCKPP_BOUNDARY_INTERPOLATE_TEMP


SUBROUTINE MCKPP_BOUNDARY_INTERPOLATE_SAL()
  
  INTEGER prev_time,next_time,true_time
  REAL prev_weight,next_weight,ndays_upd_sal
  CHARACTER(LEN=31) :: routine = "mckpp_boundary_interpolate_temp"
  CHARACTER(LEN=max_message_len) :: message

  REAL, allocatable :: prev_sal(:,:),next_sal(:,:)
  allocate(prev_sal(npts_local, nzp1))
  allocate(next_sal(npts_local, nzp1))

  true_time=time
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
  
  time=prev_time
  CALL MCKPP_READ_SALINITY_3D()
  prev_sal=kpp_3d_fields%sal_clim
  
  ! Read ocean salinity for next time
  next_time=prev_time+ndays_upd_sal
  next_weight=1-prev_weight

  WRITE(message,*) "next_time = ", next_time
  CALL mckpp_print(routine, message)
  WRITE(message,*) "next_weight = ", next_weight
  CALL mckpp_print(routine, message)

  time=next_time
  CALL MCKPP_READ_SALINITY_3D()
  next_sal=kpp_3d_fields%sal_clim  
  kpp_3d_fields%sal_clim=next_sal*next_weight+prev_sal*prev_weight
  
  time=true_time
  deallocate(prev_sal)
  deallocate(next_sal)
  
END SUBROUTINE MCKPP_BOUNDARY_INTERPOLATE_SAL

END MODULE mckpp_boundary_interpolate
