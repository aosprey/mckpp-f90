#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif

MODULE mckpp_initialize_optics_mod
  
#ifdef MCKPP_CAM3
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_types,only: kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
#endif  
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_var
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, npts

  IMPLICIT NONE

CONTAINS

SUBROUTINE mckpp_initialize_optics()
  
#ifdef MCKPP_CAM3
  INTEGER :: jerlov_temp(PLON,PLAT),jerlov_chunk(PCOLS,begchunk:endchunk)
  INTEGER :: ichnk,ncol,icol
#else
  INTEGER :: jerlov_temp(nx,ny)
#endif
  INTEGER :: ipt, ncid, offset_lon, offset_lat
  INTEGER, DIMENSION(2) :: start
  INTEGER, DIMENSION(1) :: shape
  CHARACTER(LEN=max_nc_filename_len) :: file
  CHARACTER(LEN=23) :: routine = "MCKPP_INITIALIZE_OPTICS"
  CHARACTER(LEN=max_message_len) :: message
  
  IF (kpp_const_fields%L_JERLOV) THEN
#ifdef MCKPP_CAM3
    IF (masterproc) THEN
#endif
      file = kpp_const_fields%paras_file
      WRITE(message,*) "Reading optical properties of seawater from ", TRIM(file)
      CALL mckpp_print(routine, message)
      
      CALL mckpp_netcdf_open(routine, file, ncid)
      CALL mckpp_netcdf_determine_boundaries(routine, file, ncid, &
          kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), offset_lon, offset_lat)
      start(1) = offset_lon
      start(2) = offset_lat
      CALL mckpp_netcdf_get_var(routine, file, ncid, "jerlov", jerlov_temp, start)
      CALL mckpp_netcdf_close(routine, file, ncid)

#ifdef MCKPP_CAM3
     ENDIF
     CALL scatter_field_to_chunk_int(1,1,1,PLON,jerlov_temp,jerlov_chunk(1,begchunk))     
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%jerlov(1:ncol)=jerlov_chunk(1:ncol,ichnk)
     ENDDO
#else
     shape(1) = npts
     kpp_3d_fields%jerlov = RESHAPE(jerlov_temp, shape) 
#endif
     
  ELSE ! no optics file 
#ifdef MCKPP_CAM3
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%jerlov(1:ncol)=3
     ENDDO
#else
     DO ipt=1,npts
        kpp_3d_fields%jerlov(ipt)=3
     ENDDO
#endif
  ENDIF
  
END SUBROUTINE mckpp_initialize_optics

END MODULE mckpp_initialize_optics_mod
