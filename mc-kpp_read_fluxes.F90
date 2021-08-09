SUBROUTINE MCKPP_READ_FLUXES(taux, tauy, swf, lwf, lhf, shf, rain, snow)

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, npts
  USE mckpp_time_control, ONLY: mckpp_get_update_time

  IMPLICIT NONE

#include <netcdf.inc>

  REAL, DIMENSION(npts), INTENT(OUT) :: taux, tauy, swf, lwf, lhf, shf, rain, snow
  
  REAL*4, DIMENSION(nx,ny) :: var_in
  REAL*4, DIMENSION(:), ALLOCATABLE :: file_times
  REAL*4 :: first_timein, last_timein, time, update_time
  INTEGER :: ipt, ix, iy
  INTEGER :: status, flx_ncid, time_varid, num_times
  INTEGER, DIMENSION(3) :: count, start

  CHARACTER(LEN=17) :: routine = "MCKPP_READ_FLUXES"
  CHARACTER(LEN=max_message_len) :: message

  status=NF_OPEN(kpp_const_fields%forcing_file,0,flx_ncid)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  count = (/nx,ny,1/)
  start = (/1,1,1/)

  ! Boundaries of data in file
  CALL mckpp_print(routine, "Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES")
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(flx_ncid,'fluxes','latitude','longitude','time',kpp_3d_fields%dlon(1),&
       kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid,num_times)

  ALLOCATE(file_times(num_times))
  status=NF_GET_VARA_REAL(flx_ncid,time_varid,1,num_times,file_times)
   
  ! Work out time to read and check valid
  CALL mckpp_get_update_time(kpp_const_fields%forcing_file, kpp_const_fields%time, & 
      kpp_const_fields%ndtocn, file_times, num_times, .FALSE., 0, update_time, start(3))
  WRITE(message,*) 'Reading fluxes for time ', update_time
  CALL mckpp_print(routine, message)
  WRITE(message,*) 'Reading fluxes from time point ',start(3)
  CALL mckpp_print(routine, message)
  
  CALL mckpp_print(routine, "Read taux")
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(1),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      taux(ipt)=var_in(ix,iy)
    ENDDO
  ENDDO

  CALL mckpp_print(routine, "Read tauy")
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(2),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      tauy(ipt)=var_in(ix,iy)
    ENDDO
  ENDDO
 
  CALL mckpp_print(routine, "Read swf")
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(3),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      swf(ipt)=var_in(ix,iy)
    ENDDO
  ENDDO
 
  CALL mckpp_print(routine, "Read lwf")
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(4),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      lwf(ipt)=var_in(ix,iy)
    ENDDO
  ENDDO

  CALL mckpp_print(routine, "Read lhf")
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(5),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      lhf(ipt)=var_in(ix,iy)
    ENDDO
  ENDDO
 
  CALL mckpp_print(routine, "Read shf")
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(6),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      shf(ipt)=var_in(ix,iy)
    ENDDO
  ENDDO
 
  CALL mckpp_print(routine, "Read rain")
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(7),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      rain(ipt)=var_in(ix,iy)
    ENDDO
 ENDDO
 
  CALL mckpp_print(routine, "Set snow to zero")
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      snow(ipt)=0.0
    ENDDO
  ENDDO
 
  CALL mckpp_print(routine, "Finished reading fluxes") 
  
END SUBROUTINE MCKPP_READ_FLUXES
