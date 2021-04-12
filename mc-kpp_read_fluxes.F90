SUBROUTINE MCKPP_READ_FLUXES(taux, tauy, swf, lwf, lhf, shf, rain, snow)

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_parameters, ONLY: nx, ny, npts, nuout, nuerr

  IMPLICIT NONE

#include <netcdf.inc>

  REAL, DIMENSION(npts), INTENT(OUT) :: taux, tauy, swf, lwf, lhf, shf, rain, snow
  
  REAL*4, DIMENSION(nx,ny) :: var_in
  REAL*4 :: first_timein, last_timein, time, time_in
  INTEGER :: ipt, ix, iy
  INTEGER :: status, flx_ncid, time_varid
  INTEGER, DIMENSION(3) :: count, start

  status=NF_OPEN(kpp_const_fields%forcing_file,0,flx_ncid)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  count = (/nx,ny,1/)
  start = (/1,1,1/)

  ! Boundaries of data in file
  WRITE(nuout,*) 'MCKPP_READ_FLUXES: Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES'
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(flx_ncid,'fluxes','latitude','longitude','time',kpp_3d_fields%dlon(1),&
       kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)

  ! Work out time and check valid 
  time=kpp_const_fields%time+0.5*kpp_const_fields%dtsec/kpp_const_fields%spd
  WRITE(nuout,*) 'MCKPP_READ_FLUXES: Reading fluxes for time ',&
       time,kpp_const_fields%time,kpp_const_fields%dtsec,kpp_const_fields%spd
  
  start(3)=MAX(NINT((time-first_timein)*kpp_const_fields%spd/kpp_const_fields%dtsec)+1,1)
  WRITE(nuout,*) 'MCKPP_READ_FLUXES: Reading time from time point ',start(3)

  status=NF_GET_VAR1_REAL(flx_ncid,time_varid,start(3),time_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (abs(time_in-time) .GT. 0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
    WRITE(nuerr,*) 'MCKPP_READ_FLUXES: Cannot find time,',time,'in fluxes file'
    WRITE(nuerr,*) 'MCKPP_READ_FLUXES: The closest I came was',time_in
    CALL MCKPP_ABORT()
  ENDIF

  WRITE(nuout,*) 'MCKPP_READ_FLUXES: Reading fluxes from time point ',start(3)
  WRITE(nuout,*) 'MCKPP_READ_FLUXES: Read taux'
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(1),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      taux(ipt)=var_in(ix,iy)
    ENDDO
  ENDDO

  WRITE(nuout,*) 'MCKPP_READ_FLUXES: Read tauy'
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(2),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      tauy(ipt)=var_in(ix,iy)
    ENDDO
  ENDDO
 
  WRITE(nuout,*) 'MCKPP_READ_FLUXES: Read swf'
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(3),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      swf(ipt)=var_in(ix,iy)
    ENDDO
  ENDDO
 
  WRITE(nuout,*) 'MCKPP_READ_FLUXES: Read lwf'
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(4),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      lwf(ipt)=var_in(ix,iy)
    ENDDO
  ENDDO

  WRITE(nuout,*) 'MCKPP_READ_FLUXES: Read lhf'
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(5),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      lhf(ipt)=var_in(ix,iy)
    ENDDO
  ENDDO
 
  WRITE(nuout,*) 'MCKPP_READ_FLUXES: Read shf'
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(6),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      shf(ipt)=var_in(ix,iy)
    ENDDO
  ENDDO
 
  WRITE(nuout,*) 'MCKPP_READ_FLUXES: Read rain'
  status=NF_GET_VARA_REAL(flx_ncid,kpp_const_fields%flx_varin_id(7),start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      rain(ipt)=var_in(ix,iy)
    ENDDO
 ENDDO
 
  WRITE(nuout,*) 'MCKPP_READ_FLUXES: Set snow to zero'
  DO iy=1,ny
    DO ix=1,nx
      ipt=(iy-1)*nx+ix
      snow(ipt)=0.0
    ENDDO
  ENDDO
 
  WRITE(nuout,*) 'MCKPP_READ_FLUXES: Finished reading fluxes'  
  
END SUBROUTINE MCKPP_READ_FLUXES
