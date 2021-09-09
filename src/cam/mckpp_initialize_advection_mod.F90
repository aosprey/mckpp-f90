#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>

MODULE mckpp_initialize_advection_mod

  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_types,only: kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_var
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, npts, maxmodeadv
  
  IMPLICIT NONE

CONTAINS 

SUBROUTINE MCKPP_INITIALIZE_ADVECTION()
  
  REAL(r8) :: advection_chunk(PCOLS,begchunk:endchunk,2)
  INTEGER :: nmodeadv_temp(PLON,PLAT,2),modeadv_temp(PLON,PLAT,maxmodeadv,2),&
       nmodeadv_chunk(PCOLS,begchunk:endchunk,2)
  REAL(r8) :: advection_temp(PLON,PLAT,maxmodeadv,2)
  INTEGER :: ichnk,ncol,icol
  INTEGER, DIMENSION(3) :: start
  INTEGER, DIMENSION(2) :: shape
  INTEGER :: i, ipt, ivar, ncid, offset_lon, offset_lat
  CHARACTER(LEN=max_nc_filename_len) :: file
  CHARACTER(LEN=26) :: routine = "MCKPP_INITIALIZE_ADVECTION"
  CHARACTER(LEN=max_message_len) :: message
  
  IF (kpp_const_fields%L_ADVECT) THEN
    IF (masterproc) THEN
      file = kpp_const_fields%advect_file
      WRITE(message,*) "Reading advection from file ", TRIM(file)
      CALL mckpp_print(routine, message)
      CALL mckpp_netcdf_open(routine, file, ncid)
      
      CALL mckpp_netcdf_determine_boundaries(routine, file, ncid, &
          kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), offset_lon, offset_lat)
      start(1) = offset_lon
      start(2) = offset_lat
      start(3) = 1

      CALL mckpp_netcdf_get_var(routine, file, ncid, "nmode_tadv", nmodeadv_temp(:,:,1), start(1:2)) 
      CALL mckpp_netcdf_get_var(routine, file, ncid, "nmode_sadv", nmodeadv_temp(:,:,2), start(1:2))
      CALL mckpp_netcdf_get_var(routine, file, ncid, "mode_tadv", modeadv_temp(:,:,:,1), start) 
      CALL mckpp_netcdf_get_var(routine, file, ncid, "mode_sadv", modeadv_temp(:,:,:,2), start)
      CALL mckpp_netcdf_get_var(routine, file, ncid, "tadv", advection_temp(:,:,:,1), start) 
      CALL mckpp_netcdf_get_var(routine, file, ncid, "sadv", advection_temp(:,:,:,2), start)

      CALL mckpp_netcdf_close(routine, file, ncid)
     ENDIF

     CALL scatter_field_to_chunk_int(1,1,2,PLON,nmodeadv_temp,nmodeadv_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%nmodeadv(1:ncol,:)=nmodeadv_chunk(1:ncol,ichnk,:)
     ENDDO
     DO i=1,maxmodeadv 
        CALL scatter_field_to_chunk_int(1,1,2,PLON,modeadv_temp(:,:,i,:),nmodeadv_chunk(1,begchunk,1))
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           kpp_3d_fields(ichnk)%modeadv(1:ncol,i,:)=nmodeadv_chunk(1:ncol,ichnk,:)
        ENDDO
     ENDDO
     DO i=1,maxmodeadv
        CALL scatter_field_to_chunk(1,1,2,PLON,advection_temp(:,:,i,:),advection_chunk(1,begchunk,1))
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           kpp_3d_fields(ichnk)%advection(1:ncol,i,:)=advection_chunk(1:ncol,ichnk,:)
        ENDDO
     ENDDO

   ELSE
     CALL mckpp_print(routine, "No advection has been specified.")
     
     DO ichnk = begchunk, endchunk
       ncol=get_ncols_p(ichnk)
       kpp_3d_fields(ichnk)%nmodeadv(1:ncol,:)=0
     ENDDO
  ENDIF

END SUBROUTINE mckpp_initialize_advection

END MODULE mckpp_initialize_advection_mod
#endif
