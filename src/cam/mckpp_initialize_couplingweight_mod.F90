#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>

MODULE mckpp_initialize_couplingweight_mod

  USE shr_kind_mod,only: r8=>shr_kind_r8, r4=>shr_kind_r4
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid, only: scatter_field_to_chunk,get_ncols_p
  USE mckpp_types,only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE mckpp_netcdf_read, ONLY: mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_get_var, max_nc_filename_len
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: nx_globe, ny_globe

  IMPLICIT NONE

CONTAINS

SUBROUTINE MCKPP_INITIALIZE_COUPLINGWEIGHT()
  
  REAL(r8) :: cplwght(PLON,PLAT),cplwght_chunk(PCOLS,begchunk:endchunk)
  INTEGER :: icol,ncol,ichnk,lat_varid,lon_varid
  REAL(r4) :: latitude_in(PLAT),longitude_in(PLON)
  INTEGER ix, jy, ipoint_globe, ncid 
  REAL cplwght_in(NX_GLOBE,NY_GLOBE)
  CHARACTER(LEN=max_nc_filename_len) :: file
  CHARACTER(LEN=31) :: routine = "MCKPP_INITIALIZE_COUPLINGWEIGHT"
  CHARACTER(LEN=max_message_len) :: message

!  If L_CPLWGHT has been set, then we will use the
!  NetCDF file to set values of cplwght over the
!  entire globe.
!  Otherwise, we will set the values ourselves, based
!  on the coupling region.
!  NPK 10/9/07 - R1
  
  IF (masterproc) THEN
    file = kpp_const_fields%cplwght_file
    WRITE(message,*) "Reading coupling weight (alpha) from file ", TRIM(file)
    CALL mckpp_print(routine, message)
    CALL mckpp_netcdf_open(routine, file, ncid)
    CALL mckpp_netcdf_get_var(routine, file, ncid, "alpha", cplwght_in)
    
    ! Use coupling weight to get global latitude and longitude grid
    CALL mckpp_netcdf_get_var(routine, file, ncid, "latitude", kpp_global_fields%latitude)
    CALL mckpp_netcdf_get_var(routine, file, ncid, "longitude", kpp_global_fields%longitude)

    CALL mckpp_netcdf_close(routine, file, ncid)
  ENDIF ! End of masterproc section
  
  cplwght(:,:)=cplwght_in(:,:)
  CALL scatter_field_to_chunk(1,1,1,PLON,cplwght,cplwght_chunk(1,begchunk))
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%cplwght(1:ncol)=cplwght_chunk(1:ncol,ichnk)        
  ENDDO
  
END SUBROUTINE mckpp_initialize_couplingweight

END MODULE mckpp_initialize_couplingweight_mod
#endif
