#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif 

MODULE mckpp_read_sst_mod

  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_netcdf_subs

CONTAINS

SUBROUTINE MCKPP_READ_SST()

#ifdef MCKPP_CAM3  
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_types,only: kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
#endif
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, nx_globe, ny_globe

  IMPLICIT NONE
#include <netcdf.inc>

#ifdef MCKPP_CAM3
  REAL(r8) :: sst_temp(PLON,PLAT), sst_chunk(PCOLS,begchunk:endchunk)
  INTEGER :: ichnk,ncol,icol
#endif

  INTEGER :: sst_nx, sst_ny
  REAL :: offset_sst
  INTEGER :: status,ncid
  REAL*4 :: time_in,first_timein,last_timein,sstclim_time
  REAL*4, ALLOCATABLE :: var_in(:,:,:), longitudes(:),latitudes(:)
  INTEGER :: varid, time_varid,lat_varid,lon_varid,lon_dimid,lat_dimid,time_dimid
  INTEGER count(3),start(3)
  INTEGER ix,iy,nlat_file,nlon_file,ntime_file
  CHARACTER(LEN=30) tmp_name

  CHARACTER(LEN=14) :: routine = "MCKPP_READ_SST"
  CHARACTER(LEN=max_message_len) :: message
  
#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif

#ifdef MCKPP_COUPLE
  sst_nx = nx_globe 
  sst_ny = ny_globe 
#else
  sst_nx = nx
  sst_ny = ny 
#endif

  ALLOCATE( var_in(sst_nx,sst_ny,1) ) 
  ALLOCATE( longitudes(NX_GLOBE) ) 
  ALLOCATE( latitudes(NY_GLOBE) )

  status=NF_OPEN(kpp_const_fields%sst_file,0,ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)

  count=(/sst_nx,sst_ny,1/)
  start=(/1,1,1/)
  
  status=NF_INQ_VARID(ncid,'sst',varid)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
#ifndef MCKPP_COUPLE
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(ncid,'SST climatology','latitude','longitude','t',kpp_3d_fields%dlon(1),&
       kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)
#else
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(ncid,'t',time_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIMID(ncid,'t',time_dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIM(ncid,time_dimid,tmp_name,ntime_file)
  status=NF_GET_VAR1_REAL(ncid,time_varid,start(3),first_timein)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VAR1_REAL(ncid,time_varid,ntime_file,last_timein)
#endif
  sstclim_time=kpp_const_fields%time+0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdsst
  IF (sstclim_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_CLIMSST) THEN
        DO WHILE (sstclim_time .gt. last_timein)
           sstclim_time=sstclim_time-kpp_const_fields%climsst_period
        ENDDO
     ELSE
        WRITE(message,*) "Time for which to read SST exceeds the last time in the netCDF file ", & 
            "and L_PERIODIC_CLIMSST has not been specified."
        CALL mckpp_print_error(routine, message) 
        WRITE(message,*) "Attempting to read SST will lead to an error, so aborting now ..."
        CALL mckpp_print_error(routine, message) 
        CALL MCKPP_ABORT()
     ENDIF
  ENDIF
  WRITE(message,*) 'Reading climatological SST for time ',sstclim_time
  CALL mckpp_print(routine, message)
  start(3)=NINT((sstclim_time-first_timein)*kpp_const_fields%spd/&
       (kpp_const_fields%dto*kpp_const_fields%ndtupdsst))+1
  WRITE(message,*) 'Reading climatological SST from position ',start(3)
  CALL mckpp_print(routine, message)

  status=NF_GET_VAR1_REAL(ncid,time_varid,start(3),time_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  IF (abs(time_in-sstclim_time) .GT. 0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     WRITE(message,*) 'Cannot find time,',sstclim_time,'in SST climatology file'
     CALL mckpp_print_error(routine, message) 
     WRITE(message,*) 'The closest I came was',time_in
     CALL mckpp_print_error(routine, message) 
     CALL MCKPP_ABORT()
  ENDIF
 
  status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)  
  status=NF_CLOSE(ncid)

  ! KPP expects temperatures in CELSIUS.  If climatological SSTs are
  ! in Kelvin, subtract 273.15.
  offset_sst=0.
  DO ix=1,sst_nx
     DO iy=1,sst_ny
        IF (var_in(ix,iy,1) .gt. 200 .and. var_in(ix,iy,1) .lt. 400) &
             offset_sst = 273.15
     END DO
  END DO
  
#ifdef MCKPP_CAM3
  sst_temp=var_in(:,:,1)
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,1,PLON,sst_temp,sst_chunk(1,begchunk))
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%sst(1:ncol)=sst_chunk(1:ncol,ichnk)     
     IF (.NOT. kpp_const_fields%L_CLIMICE) kpp_3d_fields(ichnk)%iceconc(:)=0.0
     IF (.NOT. kpp_const_fields%L_CLIMCURR) THEN
        kpp_3d_fields(ichnk)%usf(:)=0.0
        kpp_3d_fields(ichnk)%vsf(:)=0.0
     ENDIF
  ENDDO
#else
  DO ix=1,sst_nx
     DO iy=1,sst_ny
        kpp_3d_fields%sst(ix,iy)=var_in(ix,iy,1)-offset_sst          
        IF (.NOT. kpp_const_fields%L_CLIMICE) kpp_3d_fields%iceconc(ix,iy)=0.0
        IF (.NOT. kpp_const_fields%L_CLIMCURR) THEN
           kpp_3d_fields%usf(ix,iy)=0.0
           kpp_3d_fields%vsf(ix,iy)=0.0
        ENDIF
     ENDDO
  ENDDO
#endif
  
END SUBROUTINE mckpp_read_sst

END MODULE mckpp_read_sst_mod
