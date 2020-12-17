#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE MCKPP_READ_TEMPERATURES_3D  
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_parameters
  USE mckpp_types,only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
#else
SUBROUTINE MCKPP_READ_TEMPERATURES_3D(kpp_3d_fields,kpp_const_fields)  
  USE mckpp_data_fields
#endif

  IMPLICIT NONE
  INTEGER,parameter :: nuout=6,nuerr=0
#include <netcdf.inc>

#ifdef MCKPP_CAM3
  REAL(r8) :: ocnT_temp(PLON,PLAT,NZP1), ocnT_chunk(PCOLS,begchunk:endchunk,NZP1)
  INTEGER :: ichnk,icol,ncol
  INTEGER,parameter :: my_nx=NX_GLOBE,my_ny=NY_GLOBE
#else
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER,parameter :: my_nx=NX,my_ny=NY
#endif

  INTEGER ix,iy,iz,ipoint,ocnT_varid,status,lat_varid,lon_varid,z_varid,z_dimid,time_varid,&
       ocnT_ncid,k,lat_dimid,lon_dimid,time_dimid,nlon_file,nlat_file,ntime_file,nz_file,prev_start,&
       start(4),count(4)
  REAL*4, allocatable :: ocnT_in(:,:,:,:),latitudes(:),longitudes(:),z(:)
  REAL :: ocnT_time
  REAL*4 ixx,jyy,first_timein,time_in,ndays_upd_ocnT,last_timein
  CHARACTER(LEN=30) tmp_name

#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif

  count=(/my_nx,my_ny,NZP1,1/)
  start=(/1,1,1,1/)
  allocate(ocnT_in(MY_NX,MY_NY,NZP1,1))
  allocate(longitudes(NX_GLOBE))
  allocate(latitudes(NY_GLOBE))
  allocate(z(NZP1))
  
  status=NF_OPEN(kpp_const_fields%ocnT_file,0,ocnT_ncid)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)

  status=NF_INQ_VARID(ocnT_ncid,'z',z_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)  
  status=NF_INQ_DIMID(ocnT_ncid,'z',z_dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIM(ocnT_ncid,z_dimid,tmp_name,nz_file)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (NZP1.ne.nz_file) THEN
     WRITE(nuout,*) 'MCKPP_READ_TEMPERATURES_3D: Input file for ocean temperature climatology &
          & does not have the correct number of vertical levels. &
          & It should have ',NZP1,' but instead has ',nz_file
     CALL MCKPP_ABORT
  ELSE
     status=NF_GET_VAR_REAL(ocnT_ncid,z_varid,z)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     WRITE(nuout,*) 'Read in depths from the ocean temperature climatology input file'
  ENDIF

#ifdef MCKPP_CAM3
  WRITE(6,*) 'MCKPP_READ_TEMPERATURES_3D: Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES'  
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(ocnT_ncid,'ocean temp clim','latitude','longitude','t',&
       kpp_global_fields%longitude(1),kpp_global_fields%latitude(1),start(1),start(2),first_timein,&
       last_timein,time_varid)
  WRITE(6,*) 'MCKPP_READ_TEMPERATURES_3D: Returned from MCKPP_DETERMINE_NETCDF_BOUNDARIES'
#else
  WRITE(6,*) 'MCKPP_READ_TEMPERATURES_3D: Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES'
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(ocnT_ncid,'ocean temp clim','latitude','longitude','t',kpp_3d_fields%dlon(1),&
       kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)
  WRITE(6,*) 'MCKPP_READ_TEMPERATURES_3D: Returned from MCKPP_DETERMINE_NETCDF_BOUNDARIES'
#endif
  
  status=NF_INQ_VARID(ocnT_ncid,'temperature',ocnT_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  ndays_upd_ocnT = kpp_const_fields%ndtupdocnT*kpp_const_fields%dto/kpp_const_fields%spd
  ocnT_time=(ndays_upd_ocnT)*(FLOOR(kpp_const_fields%time,8)*NINT(kpp_const_fields%spd,8)/&
       (kpp_const_fields%ndtupdocnT*NINT(kpp_const_fields%dto,8)))+&
       (0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdocnT)
  WRITE(nuout,*) ocnT_time,last_timein
     
  IF (ocnT_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_OCNT) THEN 
        DO WHILE (ocnT_time .gt. last_timein)
           ocnT_time=ocnT_time-kpp_const_fields%ocnT_period
        ENDDO
     ELSE
        WRITE(nuout,*) 'MCKPP_READ_TEMPERATURES_3D: Time for which to read the ocean temperatures exceeds &
             & the last time in the netCDF file and L_PERIODIC_OCNT has not been specified.  Attempting to &
             & read ocean temperatures will lead to an error, so aborting now ...'
        CALL MCKPP_ABORT
     ENDIF
  ENDIF

  write(nuout,*) 'MCKPP_READ_TEMPERATUERS_3D: Reading ocean temperature for time ',ocnT_time
  start(4)=NINT((ocnT_time-first_timein)*kpp_const_fields%spd/(kpp_const_fields%dto*kpp_const_fields%ndtupdocnT))+1
  status=NF_GET_VAR1_REAL(ocnT_ncid,time_varid,start(4),time_in)
  
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (abs(time_in-ocnT_time) .GT. 0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     write(nuerr,*) 'MCKPP_READ_TEMPERATURES_3D: Cannot find time ',ocnT_time,&
          ' in ocean temperature climatology input file'
     write(nuerr,*) 'MCKPP_READ_TEMPERATURES_3D: The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
  status=NF_GET_VARA_REAL(ocnT_ncid,ocnT_varid,start,count,ocnT_in)
  status=NF_CLOSE(ocnT_ncid)
  
#ifdef MCKPP_CAM3
  ocnT_temp=ocnT_in(:,:,:,1)
  deallocate(ocnT_in)
  deallocate(longitudes)
  deallocate(latitudes)
  deallocate(z)
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,NZP1,PLON,ocnT_temp,ocnT_chunk(1,begchunk,1))
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%ocnT_clim(1:ncol,1:NZP1)=ocnT_chunk(1:ncol,ichnk,1:NZP1)
  ENDDO
#else
  ! Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
  ! into one long array with dimension NPTS.       
  DO ix=1,NX
     DO iy=1,NY
        ipoint=(iy-1)*nx+ix
        DO k=1,NZP1
           kpp_3d_fields%ocnT_clim(ipoint,k)=ocnT_in(ix,iy,k,1)
        ENDDO
     ENDDO
  ENDDO
  deallocate(ocnT_in)
  deallocate(longitudes)
  deallocate(latitudes)
  deallocate(z)
#endif
  
  RETURN
END SUBROUTINE MCKPP_READ_TEMPERATURES_3D

#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE MCKPP_READ_TEMPERATURES_BOTTOM
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_parameters
  USE mckpp_types,only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
#else
SUBROUTINE MCKPP_READ_TEMPERATURES_BOTTOM(kpp_3d_fields,kpp_const_fields)
  USE mckpp_data_fields
#endif

  IMPLICIT NONE
  INTEGER,parameter :: nuout=6,nuerr=0
  
#include <netcdf.inc>

#ifdef MCKPP_CAM3
  REAL(r8) :: bottom_temp(PLON,PLAT), bottom_chunk(PCOLS,begchunk:endchunk)
  INTEGER :: ichnk,icol,ncol
  INTEGER, parameter :: my_nx=NX_GLOBE,my_ny=NY_GLOBE
#else
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER, parameter :: my_nx=NX, my_ny=NY
#endif

  INTEGER status,ncid
  REAL :: bottomclim_time
  REAL*4 var_in(my_nx,my_ny,1),time_in,first_timein,latitudes(NY_GLOBE),&
       longitudes(NX_GLOBE), last_timein, offset_temp
  INTEGER varid,time_varid,lat_varid,lon_varid,time_dimid,lat_dimid,lon_dimid,&
       nlat_file,nlon_file,ntime_file,count(3),start(3),ix,iy,ipoint
  CHARACTER(LEN=30) tmp_name
  
#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif

  count(:)=(/my_nx,my_ny,1/)
  start(:)=(/1,1,1/)

  status=NF_OPEN(kpp_const_fields%bottom_file,0,ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)  
  status=NF_INQ_VARID(ncid,'T',varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  WRITE(6,*) 'MCKPP_READ_TEMPERATURES_BOTTOM: Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES'    
#ifdef MCKPP_CAM3 
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(ncid,'bottom temp climatology','latitude','longitude',&
       't',kpp_global_fields%longitude(1),kpp_global_fields%latitude(1),start(1),start(2),&
       first_timein,last_timein,time_varid)
#else
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(ncid,'bottom temp climatology','latitude','longitude',&
       't',kpp_3d_fields%dlon(1),kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)
#endif
  WRITE(6,*) 'MCKPP_READ_TEMPERATURES_BOTTOM: Returned from MCKPP_DETERMINE_NETCDF_BOUNDARIES'
  
  bottomclim_time=kpp_const_fields%time+0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdbottom
  IF (bottomclim_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_BOTTOM_TEMP) THEN 
        DO WHILE (bottomclim_time .gt. last_timein)
           bottomclim_time=bottomclim_time-kpp_const_fields%bottom_temp_period
        ENDDO
     ELSE
        WRITE(nuout,*) 'MCKPP_READ_TEMPERATURES_BOTTOM: Time for which to read bottom temperature exceeds &
             & the last time in the netCDF file and L_PERIODIC_BOTTOM_TEMP has not been specified.  &
             & Attempting to read bottom temperature will lead to an error, so aborting now...'
        CALL MCKPP_ABORT
     ENDIF
  ENDIF
 
  write(nuout,*) 'MCKPP_READ_TEMPERATURES_BOTTOM: Reading climatological bottom temp for time ',bottomclim_time
  start(3)=NINT((bottomclim_time-first_timein)*kpp_const_fields%spd/&
       (kpp_const_fields%dto*kpp_const_fields%ndtupdbottom))+1
  
  status=NF_GET_VAR1_REAL(ncid,time_varid,start(3),time_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (abs(time_in-bottomclim_time) .GT. 0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     write(nuerr,*) 'MCKPP_READ_TEMPERATURES_BOTTOM: Cannot find time',bottomclim_time,&
          'in bottom temperature climatology file'
     write(nuerr,*) 'MCKPP_READ_TEMPERATURES_BOTTOM: The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
  status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
!     KPP expects temperatures in CELSIUS.  If climatological bottom 
!     temperatures are in Kelvin, subtract 273.15.     
  offset_temp=0.
  ix=1
  iy=1
  DO WHILE (offset_temp.EQ.0.AND.ix.LE.my_nx)
     DO iy=1,my_ny
        IF (var_in(ix,iy,1) .gt. 200 .and. var_in(ix,iy,1) .lt. 400) offset_temp = 273.15    
     END DO
     ix=ix+1
  ENDDO
  status=NF_CLOSE(ncid)

#ifdef MCKPP_CAM3
  bottom_temp=var_in(:,:,1)
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,1,PLON,bottom_temp,bottom_chunk(1,begchunk))
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%bottom_temp(1:ncol)=bottom_chunk(1:ncol,ichnk)
  ENDDO
#else
  DO ix=1,NX
     DO iy=1,NY
        ipoint=(iy-1)*NX+ix
        kpp_3d_fields%bottom_temp(ipoint) = var_in(ix,iy,1)-offset_temp        
     ENDDO
  ENDDO
#endif
  
  RETURN
END SUBROUTINE MCKPP_READ_TEMPERATURES_BOTTOM
