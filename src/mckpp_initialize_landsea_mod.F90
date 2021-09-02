#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif

MODULE mckpp_initialize_landsea_mod

#ifdef MCKPP_CAM3 
  USE shr_kind_mod, only : r8=>shr_kind_r8
  USE mckpp_types, only: kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, get_ncols_p
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
#endif
  USE mckpp_netcdf_read, ONLY: mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_var
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_warning, max_message_len
  USE mckpp_parameters, ONLY: npts, nx, ny
 
  IMPLICIT NONE

CONTAINS

SUBROUTINE mckpp_initialize_landsea()
 
#ifdef MCKPP_CAM3
  REAL(r8) :: landsea(PLON,PLAT)
  REAL(r8) :: ocdepth(PLON,PLAT)
  REAL(r8) :: landsea_chunk(PCOLS,begchunk:endchunk),ocdepth_chunk(PCOLS,begchunk:endchunk)
  INTEGER :: ichnk, icol, ncol
#else
  REAL, DIMENSION(npts) :: landsea, ocdepth
#endif

  INTEGER :: ncid, offset_lon, offset_lat, i, j, ipt
  INTEGER, DIMENSION(2) :: start, count
  REAL, DIMENSION(nx) :: lon_in
  REAL, DIMENSION(ny) :: lat_in
  CHARACTER(LEN=24) :: routine = "MCKPP_INITIALIZE_LANDSEA"
  CHARACTER(LEN=max_message_len) :: message

#ifdef MCKPP_CAM3
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     CALL GET_RLON_ALL_P(ichnk,ncol,clon1)
     CALL GET_RLAT_ALL_P(ichnk,ncol,clat1)
     kpp_3d_fields(ichnk)%dlon(:)=clon1*360./kpp_const_fields%twopi    
     kpp_3d_fields(ichnk)%dlat(:)=clat1*360./kpp_const_fields%twopi
  ENDDO
#else
  kpp_3d_fields%dlon(1)=kpp_const_fields%alon
  kpp_3d_fields%dlat(1)=kpp_const_fields%alat
#endif

  IF (kpp_const_fields%L_LANDSEA) THEN     
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
       WRITE(message,*) "Reading", kpp_const_fields%landsea_file
       CALL mckpp_print(routine, message)
       CALL mckpp_netcdf_open(routine, kpp_const_fields%landsea_file, ncid)
       CALL mckpp_netcdf_determine_boundaries(routine, kpp_const_fields%landsea_file, ncid, &
           kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), offset_lon, offset_lat)
       CALL mckpp_netcdf_get_var(routine, kpp_const_fields%landsea_file, ncid, &
           "longitude", lon_in, offset_lon)
       CALL mckpp_netcdf_get_var(routine, kpp_const_fields%landsea_file, ncid, &
           "latitude", lat_in, offset_lat)
       start(1) = offset_lon
       start(2) = offset_lat
       count(1) = nx
       count(2) = ny
       CALL mckpp_netcdf_get_var(routine, kpp_const_fields%landsea_file, ncid, &
           "lsm", landsea, start, count, 2)
       CALL mckpp_netcdf_get_var(routine, kpp_const_fields%landsea_file, ncid, &
           "maxdepth", ocdepth, start, count, 2)       
       CALL mckpp_netcdf_close(routine, kpp_const_fields%landsea_file, ncid)
       
#ifdef MCKPP_CAM3
     ENDIF    
     CALL scatter_field_to_chunk(1,1,1,PLON,ocdepth,ocdepth_chunk(1,begchunk))
     CALL scatter_field_to_chunk(1,1,1,PLON,landsea,landsea_chunk(1,begchunk))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%ocdepth(:)=ocdepth_chunk(:,ichnk)
        kpp_3d_fields(ichnk)%landfrac(:)=landsea_chunk(:,ichnk)
        DO icol=1,ncol
           IF (landsea_chunk(icol,ichnk) .EQ. 1.0) THEN
              kpp_3d_fields(ichnk)%L_OCEAN(icol)=.FALSE.
           ELSE
              kpp_3d_fields(ichnk)%L_OCEAN(icol)=.TRUE.
           ENDIF
        ENDDO
      ENDDO
   ELSE
     DO ipt=1,npts
        kpp_3d_fields%L_OCEAN(ipt)=.FALSE.
     ENDDO
  ENDIF
     
#else
     DO j = 1, ny 
       DO i =1, nx
         ipt = (j-1)*nx + i
         kpp_3d_fields%dlon(ipt) = lon_in(j) 
         kpp_3d_fields%dlat(ipt) = lat_in(i)
         IF (landsea(ipt) .EQ. 1.0) THEN
           kpp_3d_fields%L_OCEAN(ipt)=.FALSE.
         ELSE
           kpp_3d_fields%L_OCEAN(ipt)=.TRUE.
         ENDIF  
       ENDDO
     ENDDO
     kpp_3d_fields%ocdepth = ocdepth
     
   ELSEIF (kpp_const_fields%L_REGGRID) THEN
     DO j = 1, ny
        DO i = 1, nx
           ipt = (j-1)*nx + i
           kpp_3d_fields%dlon(ipt) = kpp_const_fields%alon+(i-1)*kpp_const_fields%delta_lon
           kpp_3d_fields%dlat(ipt) = kpp_const_fields%alat+(j-1)*kpp_const_fields%delta_lat
           kpp_3d_fields%ocdepth(ipt) = -10000.
           kpp_3d_fields%L_OCEAN(ipt) = .TRUE.
        ENDDO
      ENDDO
      
   ELSEIF (.NOT. kpp_const_fields%L_REGGRID .AND. .NOT. kpp_const_fields%L_LANDSEA) THEN
     WRITE(message,*) "If you set L_REGGRID=.FALSE., you must specify a land-sea mask file from which", & 
         " to read the locations of the gridpoints in the horizontal."
     CALL mckpp_print_warning(routine, message)
  ENDIF      
#endif        
 
END SUBROUTINE mckpp_initialize_landsea

END MODULE mckpp_initialize_landsea_mod

