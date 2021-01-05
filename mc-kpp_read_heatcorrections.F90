#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE MCKPP_READ_FCORR_2D
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_paramaters 
  USE mckpp_types,only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
#else
SUBROUTINE MCKPP_READ_FCORR_2D(kpp_3d_fields,kpp_const_fields)
  USE mckpp_data_fields
#endif

  IMPLICIT NONE
  INTEGER start(3),count(3)
  INTEGER ix,iy,ipoint,fcorr_varid,status,lat_varid,lon_varid,time_varid,&
       lat_dimid,lon_dimid,time_dimid,fcorr_ncid,k,nlat_file,nlon_file,ntime_file
#include <netcdf.inc>

#ifdef MCKPP_CAM3
  REAL (r8) :: fcorr_temp(PLON,PLAT), fcorr_chunk(PCOLS,begchunk:endchunk)
  INTEGER :: ichnk,icol,ncol
#else
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
#endif

  INTEGER :: my_nx, my_ny
  REAL :: fcorr_time
  REAL*4 ixx,jyy,first_timein,time_in,ndays_upd_fcorr,last_timein
  REAL*4, ALLOCATABLE :: fcorr_twod_in(:,:,:), latitudes(:), longitudes(:), z(:)
  CHARACTER(LEN=30) tmp_name
  
  ! Read in a NetCDF file containing a time-varying flux correction
  ! at the surface only.  Frequency of read is controlled by ndtupdfcorr
  ! in the namelist
  ! NPK 29/06/08

#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif

#ifdef MCKPP_CAM3
  my_nx = nx_globe 
  my_ny = ny_globe 
#else
  my_nx = nx
  my_ny = ny 
#endif

  ALLOCATE( fcorr_twod_in(MY_NX,MY_NY,1) ) 
  ALLOCATE( latitudes(MY_NY) ) 
  ALLOCATE( longitudes(MY_NX) ) 
  ALLOCATE( z(NZP1) )

  status = NF_OPEN(kpp_const_fields%fcorr_file,0,fcorr_ncid)
  IF (status.NE.0) CALL MCKPP_HANDLE_ERR(status)

  start(:)=(/1,1,1/)  
  count(:)=(/my_nx,my_ny,1/)

#ifdef MCKPP_CAM3  
  !WRITE(6,*) 'MCKPP_READ_FCORR_2D: Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES'    
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(fcorr_ncid,'flux correction','latitude','longitude',&
       't',kpp_global_fields%longitude(1),kpp_global_fields%latitude(1),start(1),start(2),&
       first_timein,last_timein,time_varid)
  !WRITE(6,*) 'MCKPP_READ_FCORR_2D: Returned from MCKPP_DETERMINE_NETCDF_BOUNDARIES'
#else
  !WRITE(6,*) 'MCKPP_READ_FCORR_2D: Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES'    
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(fcorr_ncid,'flux correction','latitude','longitude',&
       't',kpp_3d_fields%dlon(1),kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)  
  !WRITE(6,*) 'MCKPP_READ_FCORR_2D: Returned from MCKPP_DETERMINE_NETCDF_BOUNDARIES'
#endif

  status=NF_INQ_VARID(fcorr_ncid,'fcorr',fcorr_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
       
  ndays_upd_fcorr = kpp_const_fields%ndtupdfcorr*kpp_const_fields%dto/kpp_const_fields%spd
  fcorr_time=(ndays_upd_fcorr)*(FLOOR(kpp_const_fields%time,8)*NINT(kpp_const_fields%spd,8)/&
       (kpp_const_fields%ndtupdfcorr*NINT(kpp_const_fields%dto,8)))+&
       (0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdfcorr)
  
  IF (fcorr_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_FCORR) THEN 
        DO WHILE (fcorr_time .gt. last_timein)
           fcorr_time=fcorr_time-kpp_const_fields%fcorr_period
        ENDDO
     ELSE
        WRITE(nuout,*) 'MCKPP_READ_FCORR_2D: &
             & Time for which to read the flux corrections exceeds the last time in the netCDF file &
             & and L_PERIODIC_FCORR has not been specified. Attempting to read flux corrections will lead to &
             & an error, so aborting now ...'
        CALL MCKPP_ABORT
     ENDIF
  ENDIF

  write(nuout,*) 'MCKPP_READ_FCORR_2D: Reading flux correction for time ',fcorr_time
  start(3)=NINT((fcorr_time-first_timein)*kpp_const_fields%spd/(kpp_const_fields%dto*kpp_const_fields%ndtupdfcorr))+1  
  status=NF_GET_VAR1_REAL(fcorr_ncid,time_varid,start(3),time_in)
  
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (abs(time_in-fcorr_time) .GT. 0.03*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     write(nuerr,*) 'MCKPP_READ_FCORR_2D: Cannot find time',fcorr_time,'in flux-correction input file'
     write(nuerr,*) 'MCKPP_READ_FCORR_2D: The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
  status=NF_GET_VARA_REAL(fcorr_ncid,fcorr_varid,start,count,fcorr_twod_in)
  status=NF_CLOSE(fcorr_ncid)
  
#ifdef MCKPP_CAM3
  fcorr_temp=fcorr_twod_in(:,:,1)
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,1,PLON,fcorr_temp,fcorr_chunk(1,begchunk))
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%fcorr_twod(1:ncol)=fcorr_chunk(1:ncol,ichnk)
  ENDDO
#else    
  !    Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
  !    into one long array with dimension NPTS.         
  DO ix=1,NX
     DO iy=1,NY
        ipoint=(iy-1)*nx+ix
        kpp_3d_fields%fcorr_twod(ipoint)=fcorr_twod_in(ix,iy,1)
     ENDDO
  ENDDO
#endif
  
  RETURN
END SUBROUTINE MCKPP_READ_FCORR_2D

#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE MCKPP_READ_FCORR_3D
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_parameters
  USE mckpp_types,only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
#else
SUBROUTINE MCKPP_READ_FCORR_3D(kpp_3d_fields,kpp_const_fields)
  USE mckpp_data_fields
#endif
  
  IMPLICIT NONE
#include <netcdf.inc>

#ifdef MCKPP_CAM3
  REAL(r8) :: fcorr_temp(PLON,PLAT,NZP1), fcorr_chunk(PCOLS,begchunk:endchunk,NZP1)
  INTEGER :: icol,ncol,ichnk
#else
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
#endif

  INTEGER :: my_nx, my_ny
  INTEGER ix,iy,iz,ipoint,fcorr_varid,status,lat_varid,lon_varid,z_varid,z_dimid,time_varid,&
       fcorr_ncid,k,lat_dimid,lon_dimid,time_dimid,nlon_file,nlat_file,ntime_file,nz_file,&
       start(4),count(4)
  REAL :: fcorr_time
  REAL*4 ixx,jyy,first_timein,time_in,ndays_upd_fcorr,last_timein
  CHARACTER(LEN=30) tmp_name
  REAL*4, allocatable :: fcorr_in(:,:,:,:),longitudes(:),latitudes(:),z(:)
  
  ! Read in a NetCDF file containing a 
  ! time-varying flux correction at every model vertical level.
  ! Frequency of read is controlled by ndtupdfcorr in the namelist
  ! NPK 12/02/08

#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif

#ifdef MCKPP_CAM3
  my_nx = nx_globe 
  my_ny = ny_globe 
#else
  my_nx = nx
  my_ny = ny 
#endif

  count=(/my_nx,my_ny,NZP1,1/)
  start=(/1,1,1,1/)

  allocate(fcorr_in(MY_NX,MY_NY,NZP1,1))
  allocate(longitudes(NX_GLOBE))
  allocate(latitudes(NY_GLOBE))
  allocate(z(NZP1))
  
  status=NF_OPEN(kpp_const_fields%fcorr_file,0,fcorr_ncid)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)  
    
  status=NF_INQ_VARID(fcorr_ncid,'z',z_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)  
  status=NF_INQ_DIMID(fcorr_ncid,'z',z_dimid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_DIM(fcorr_ncid,z_dimid,tmp_name,nz_file)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (NZP1.ne.nz_file) THEN
     WRITE(nuout,*) 'MCKPP_READ_FCORR_3D: Input file for flux corrections does &
          & not have the correct number of vertical levels. ',&
          'It should have ',NZP1,' but instead has ',nz_file
     CALL MCKPP_ABORT
  ELSE
     status=NF_GET_VAR_REAL(fcorr_ncid,z_varid,z)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ENDIF
  
#ifdef MCKPP_CAM3  
  !WRITE(6,*) 'MCKPP_READ_FCORR_3D: Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES'    
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(fcorr_ncid,'flux correction','latitude','longitude',&
       't',kpp_global_fields%longitude(1),kpp_global_fields%latitude(1),start(1),start(2),&
       first_timein,last_timein,time_varid)
  !WRITE(6,*) 'MCKPP_READ_FCORR_3D: Returned from MCKPP_DETERMINE_NETCDF_BOUNDARIES'    
#else
  !WRITE(6,*) 'MCKPP_READ_FCORR_3D: Calling MCKPP_DETERMINE_NETCDF_BOUNDARIES'    
  CALL MCKPP_DETERMINE_NETCDF_BOUNDARIES(fcorr_ncid,'flux correction','latitude','longitude',&
       't',kpp_3d_fields%dlon(1),kpp_3d_fields%dlat(1),start(1),start(2),first_timein,last_timein,time_varid)
  !WRITE(6,*) 'MCKPP_READ_FCORR_3D: Returned from MCKPP_DETERMINE_NETCDF_BOUNDARIES'    
#endif
  status=NF_INQ_VARID(fcorr_ncid,'fcorr',fcorr_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  !ndays_upd_fcorr = kpp_const_fields%ndtupdfcorr*kpp_const_fields%dto/kpp_const_fields%spd
  !WRITE(nuout,*) ndays_upd_fcorr,FLOOR(kpp_const_fields%time,8)*NINT(kpp_const_fields%spd,8),&
  !     kpp_const_fields%ndtupdfcorr*NINT(kpp_const_fields%dto,8),&
  !     0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdfcorr
  !fcorr_time=(ndays_upd_fcorr)*FLOOR(kpp_const_fields%time,8)*NINT(kpp_const_fields%spd,8)/&
  !     FLOAT(kpp_const_fields%ndtupdfcorr*NINT(kpp_const_fields%dto,8))+&
  !     (0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdfcorr)
  fcorr_time=kpp_const_fields%time+0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdfcorr
   
  IF (fcorr_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_FCORR) THEN 
        DO WHILE (fcorr_time .gt. last_timein)
           fcorr_time=fcorr_time-kpp_const_fields%fcorr_period
        ENDDO
     ELSE
        WRITE(nuout,*) 'MCKPP_READ_FCORR_3D: Time for which to read the flux corrections exceeds the &
             & last time in the netCDF file and L_PERIODIC_FCORR has not been specified. &
             & Attempting to read flux corrections will lead to an error, so aborting now ...'
        CALL MCKPP_ABORT
     ENDIF
  ENDIF
  
  write(nuout,*) 'MCKPP_READ_FCORR_3D: Reading flux correction for time ',fcorr_time
  start(4)=NINT((fcorr_time-first_timein)*kpp_const_fields%spd/&
       (kpp_const_fields%dto*kpp_const_fields%ndtupdfcorr))+1
  status=NF_GET_VAR1_REAL(fcorr_ncid,time_varid,start(4),time_in)
      
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  IF (abs(time_in-fcorr_time) .GT. 0.01) THEN
     write(nuerr,*) 'MCKPP_READ_FCORR_3D: Cannot find time',fcorr_time,'in flux-correction input file'
     write(nuerr,*) 'The closest I came was',time_in
     CALL MCKPP_ABORT
  ENDIF
  status=NF_GET_VARA_REAL(fcorr_ncid,fcorr_varid,start,count,fcorr_in)
  status=NF_CLOSE(fcorr_ncid)

#ifdef MCKPP_CAM3
  fcorr_temp=fcorr_in(:,:,:,1)
  deallocate(fcorr_in)
  deallocate(longitudes)
  deallocate(latitudes)
  deallocate(z)
  ENDIF ! End of masterproc section  
  CALL scatter_field_to_chunk(1,1,NZP1,PLON,fcorr_temp,fcorr_chunk(1,begchunk,1))
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%fcorr_withz(1:ncol,1:NZP1) = fcorr_chunk(1:ncol,ichnk,1:NZP1)
  ENDDO
#else
  !     Convert from REAL*4 to REAL*(default precision). Put all (NX,NY) points
  !     into one long array with dimension NPTS.         
  DO ix=1,NX
     DO iy=1,NY
        ipoint=(iy-1)*nx+ix
        DO k=1,NZP1
           kpp_3d_fields%fcorr_withz(ipoint,k)=fcorr_in(ix,iy,k,1)
        ENDDO
     ENDDO
  ENDDO  
  deallocate(fcorr_in)
  deallocate(longitudes)
  deallocate(latitudes)
  deallocate(z)
#endif  

  RETURN
END SUBROUTINE MCKPP_READ_FCORR_3D
