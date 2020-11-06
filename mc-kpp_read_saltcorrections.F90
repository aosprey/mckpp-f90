#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE MCKPP_READ_SFCORR_2D
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_types,only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
#else
SUBROUTINE MCKPP_READ_SFCORR_2D(kpp_3d_fields,kpp_const_fields)
#endif

  IMPLICIT NONE
  INTEGER nuout,nuerr,start(3),count(3)
  INTEGER ix,iy,ipoint,sfcorr_varid,status,lat_varid,lon_varid,time_varid,&
       lat_dimid,lon_dimid,time_dimid,sfcorr_ncid,k,nlat_file,nlon_file,ntime_file
  PARAMETER (nuout=6,nuerr=0)     
#include <netcdf.inc>

#ifdef MCKPP_CAM3
#include <parameter.inc>
  REAL(r8) :: sfcorr_temp(PLON,PLAT), sfcorr_chunk(PCOLS,begchunk:endchunk)
  INTEGER :: ichnk,icol,ncol
  INTEGER,parameter :: my_nx=NX_GLOBE,my_ny=NY_GLOBE
#else
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER,parameter :: my_nx=NX,my_ny=NY_GLOBE
#endif

  REAL :: sfcorr_time
  REAL*4 ixx,jyy,sfcorr_twod_in(my_nx,my_ny,1),latitudes(my_ny),longitudes(my_nx),z(NZP1),&
       first_timein,time_in,ndays_upd_sfcorr,last_timein
  CHARACTER(LEN=30) tmp_name
  
  ! Read in a NetCDF file containing a time-varying salinity correction
  ! at the surface only.  Frequency of read is controlled by ndtupdsfcorr
  ! in the namelist
  ! NPK 29/06/08

#ifdef MCKPP_CAM3
  IF (masterproc) THEN     
#endif

  status = NF_OPEN(kpp_const_fields%sfcorr_file,0,sfcorr_ncid)
  IF (status.NE.0) CALL MCKPP_HANDLE_ERR(status)

  start(:)=(/1,1,1/)
  count(:)=(/my_nx,my_ny,1/)

#ifdef MCKPP_CAM3  
  !WRITE(6,*) 'MCKPP_READ_SFCORR_2D: Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES'    
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(sfcorr_ncid,'salinity correction','latitude','longitude',&
       't',kpp_global_fields%longitude(1),kpp_global_fields%latitude(1),start(1),start(2),&
       first_timein,last_timein,time_varid)
  !WRITE(6,*) 'MCKPP_READ_SFCORR_2D: Returned from MCKPP_DETERMINE_NETCDF_BOUNDARIES'        
#else
  !WRITE(6,*) 'MCKPP_READ_SFCORR_2D: Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES'    
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(sfcorr_ncid,'salinity correction','latitude','longitude',&
       't',kpp_3d_fields%dlon(1),kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)
  !WRITE(6,*) 'MCKPP_READ_SFCORR_2D: Returned from MCKPP_DETERMINE_NETCDF_BOUNDARIES'            
#endif
 
  status=NF_INQ_VARID(sfcorr_ncid,'sfcorr',sfcorr_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
       
  ndays_upd_sfcorr = kpp_const_fields%ndtupdsfcorr*kpp_const_fields%dto/kpp_const_fields%spd
  sfcorr_time=(ndays_upd_sfcorr)*(FLOOR(kpp_const_fields%time,8)*NINT(kpp_const_fields%spd,8)/&
       (kpp_const_fields%ndtupdsfcorr*NINT(kpp_const_fields%dto,8)))+&
       (0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdsfcorr)
  
  IF (sfcorr_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_SFCORR) THEN 
        DO WHILE (sfcorr_time .gt. last_timein)
           sfcorr_time=sfcorr_time-kpp_const_fields%sfcorr_period
        ENDDO
     ELSE
        WRITE(nuout,*) 'MCKPP_READ_SFCORR_2D: Time for which to read the salinity corrections &
             & exceeds the last time in the netCDF file and L_PERIODIC_SFCORR has not been specified. &
             & Attempting to read salinity corrections will lead to an error, so aborting now ...'
        CALL MCKPP_ABORT
     ENDIF
  ENDIF

  WRITE(nuout,*) 'MCKPP_READ_SFCORR_2D: Reading salinity correction for time ',sfcorr_time
  start(3)=NINT((sfcorr_time-first_timein)*kpp_const_fields%spd/(kpp_const_fields%dto*kpp_const_fields%ndtupdsfcorr))+1
  status=NF_GET_VAR1_REAL(sfcorr_ncid,time_varid,start(3),time_in)
  
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (abs(time_in-sfcorr_time) .GT. 0.03*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     write(nuerr,*) 'MCKPP_READ_SFCORR_2D: Cannot find time',sfcorr_time,'in flux-correction input file'
     write(nuerr,*) 'MCKPP_READ_SFCORR_3D: The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
  status=NF_GET_VARA_REAL(sfcorr_ncid,sfcorr_varid,start,count,sfcorr_twod_in)
  status=NF_CLOSE(sfcorr_ncid)

#ifdef MCKPP_CAM3
  sfcorr_temp=sfcorr_twod_in(:,:,1)
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,1,PLON,sfcorr_temp,sfcorr_chunk(1,begchunk))
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%sfcorr_twod(1:ncol)=sfcorr_chunk(1:ncol,ichnk)
  ENDDO
#else
  !    Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
  !    into one long array with dimension NPTS.         
  DO ix=1,NX
     DO iy=1,NY
        ipoint=(iy-1)*nx+ix
        kpp_3d_fields%sfcorr_twod(ipoint)=sfcorr_twod_in(ix,iy,1)
     ENDDO
  ENDDO
#endif

  RETURN
END SUBROUTINE MCKPP_READ_SFCORR_2D

#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE MCKPP_READ_SFCORR_3D
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_types,only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
#else
SUBROUTINE MCKPP_READ_SFCORR_3D(kpp_3d_fields,kpp_const_fields)
#endif
        
  IMPLICIT NONE
  INTEGER, parameter :: nuout=6,nuerr=0
#include <netcdf.inc>
  
#ifdef MCKPP_CAM3
#include <parameter.inc>
  REAL(r8) :: sfcorr_temp(PLON,PLAT,NZP1), sfcorr_chunk(PCOLS,begchunk:endchunk,NZP1)
  INTEGER :: icol,ncol,ichnk
  INTEGER,parameter :: my_nx=NX_GLOBE,my_ny=NY_GLOBE
#else
#include <mc-kpp_3d_type.com>
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER,parameter :: my_nx=NX,my_ny=NY_GLOBE
#endif

  INTEGER ix,iy,iz,ipoint,sfcorr_varid,status,lat_varid,lon_varid,z_varid,z_dimid,time_varid,&
       sfcorr_ncid,k,lat_dimid,lon_dimid,time_dimid,nlon_file,nlat_file,ntime_file,nz_file,&
       start(4),count(4)
  REAL :: sfcorr_time
  REAL*4 ixx,jyy,first_timein,time_in,ndays_upd_sfcorr,last_timein
  CHARACTER(LEN=30) tmp_name
  REAL*4, allocatable :: sfcorr_in(:,:,:,:),longitudes(:),latitudes(:),z(:)
  
  ! Read in a NetCDF file containing a 
  ! time-varying salinity correction at every model vertical level.
  ! Frequency of read is controlled by ndtupdsfcorr in the namelist
  ! NPK 12/02/08  

#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif

  count=(/my_nx,my_ny,NZP1,1/)
  start=(/1,1,1,1/)

  allocate(sfcorr_in(my_nx,my_ny,NZP1,1))
  allocate(longitudes(NX_GLOBE))
  allocate(latitudes(NY_GLOBE))
  allocate(z(NZP1))

  status=NF_OPEN(kpp_const_fields%sfcorr_file,0,sfcorr_ncid)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)

  status=NF_INQ_VARID(sfcorr_ncid,'z',z_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)  
  status=NF_INQ_DIMID(sfcorr_ncid,'z',z_dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIM(sfcorr_ncid,z_dimid,tmp_name,nz_file)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (NZP1.ne.nz_file) THEN
     WRITE(nuout,*) 'MCKPP_READ_SFCORR_3D: Input file for salinity corrections does not have the &
          & correct number of vertical levels. ',&
          'It should have ',NZP1,' but instead has ',nz_file
     CALL MCKPP_ABORT
  ELSE
     status=NF_GET_VAR_REAL(sfcorr_ncid,z_varid,z)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ENDIF

#ifdef MCKPP_CAM3
  !WRITE(6,*) 'MCKPP_READ_SFCORR_3D: Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES'    
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(sfcorr_ncid,'salinity correction','latitude','longitude',&
       't',kpp_global_fields%longitude(1),kpp_global_fields%latitude(1),start(1),start(2),&
       first_timein,last_timein,time_varid)
  !WRITE(6,*) 'MCKPP_READ_SFCORR_3D: Returned from MCKPP_DETERMINE_NETCDF_BOUNDARIES'        

#else
  !WRITE(6,*) 'MCKPP_READ_SFCORR_3D: Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES'    
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(sfcorr_ncid,'salinity correction','latitude','longitude',&
       't',kpp_3d_fields%dlon(1),kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)
  !WRITE(6,*) 'MCKPP_READ_SFCORR_3D: Returned from MCKPP_DETERMINE_NETCDF_BOUNDARIES'        
#endif
  status=NF_INQ_VARID(sfcorr_ncid,'sfcorr',sfcorr_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
!  ndays_upd_sfcorr = kpp_const_fields%ndtupdsfcorr*kpp_const_fields%dto/kpp_const_fields%spd
!  WRITE(nuout,*) ndays_upd_sfcorr,FLOOR(kpp_const_fields%time,8)*NINT(kpp_const_fields%spd,8),&
!       kpp_const_fields%ndtupdsfcorr*NINT(kpp_const_fields%dto,8),&
!       0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdsfcorr
!  sfcorr_time=(ndays_upd_sfcorr)*FLOOR(kpp_const_fields%time,8)*NINT(kpp_const_fields%spd,8)/&
!       FLOAT(kpp_const_fields%ndtupdsfcorr*NINT(kpp_const_fields%dto,8))+&
!       (0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdsfcorr)  
  sfcorr_time=kpp_const_fields%time+0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdsfcorr

  IF (sfcorr_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_SFCORR) THEN 
        DO WHILE (sfcorr_time .gt. last_timein)
           sfcorr_time=sfcorr_time-kpp_const_fields%sfcorr_period
        ENDDO
     ELSE
        WRITE(nuout,*) 'MCKPP_READ_SFCORR_3D: Time for which to read the salinity corrections exceeds the &
             & last time in the netCDF file and L_PERIODIC_SFCORR has not been specified. &
             & Attempting to read salinity corrections will lead to an error, so aborting now ...'
        CALL MCKPP_ABORT
     ENDIF
  ENDIF
  
  write(nuout,*) 'MCKPP_READ_SFCORR_3D: Reading salinity correction for time ',sfcorr_time
  start(4)=NINT((sfcorr_time-first_timein)*kpp_const_fields%spd/&
       (kpp_const_fields%dto*kpp_const_fields%ndtupdsfcorr))+1  
  status=NF_GET_VAR1_REAL(sfcorr_ncid,time_varid,start(4),time_in)
      
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (abs(time_in-sfcorr_time) .GT. 0.01) THEN
     write(nuerr,*) 'MCKPP_READ_SFCORR_3D: Cannot find time',sfcorr_time,'in flux-correction input file'
     write(nuerr,*) 'MCKPP_READ_SFCORR_3D: The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
  status=NF_GET_VARA_REAL(sfcorr_ncid,sfcorr_varid,start,count,sfcorr_in)
  status=NF_CLOSE(sfcorr_ncid)

#ifdef MCKPP_CAM3
  sfcorr_temp=sfcorr_in(:,:,:,1) 
  deallocate(sfcorr_in)
  deallocate(longitudes)
  deallocate(latitudes)
  deallocate(z)
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,NZP1,PLON,sfcorr_temp,sfcorr_chunk(1,begchunk,1))
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%sfcorr_withz(1:ncol,1:NZP1) = sfcorr_chunk(1:ncol,ichnk,1:NZP1)
  ENDDO
#else  
  !     Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
  !     into one long array with dimension NPTS.         
  DO ix=1,NX
     DO iy=1,NY
        ipoint=(iy-1)*nx+ix
        DO k=1,NZP1
           kpp_3d_fields%sfcorr_withz(ipoint,k)=sfcorr_in(ix,iy,k,1)
        ENDDO
     ENDDO
  ENDDO  
  deallocate(sfcorr_in)
  deallocate(longitudes)
  deallocate(latitudes)
  deallocate(z)
#endif
  
  RETURN
END SUBROUTINE MCKPP_READ_SFCORR_3D
