#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif

MODULE mckpp_initialize_couplingweight_mod

  USE mckpp_netcdf_subs

CONTAINS

SUBROUTINE MCKPP_INITIALIZE_COUPLINGWEIGHT()
  
#ifdef MCKPP_CAM3
  USE shr_kind_mod,only: r8=>shr_kind_r8, r4=>shr_kind_r4
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid, only: scatter_field_to_chunk,get_ncols_p
  USE mckpp_types,only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
#endif
  USE mckpp_netcdf_read, ONLY: mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_get_var
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: nx_globe, ny_globe

  IMPLICIT NONE

#ifdef MCKPP_CAM3
  REAL(r8) :: cplwght(PLON,PLAT),cplwght_chunk(PCOLS,begchunk:endchunk)
  INTEGER :: icol,ncol,ichnk,lat_varid,lon_varid
  REAL(r4) :: latitude_in(PLAT),longitude_in(PLON)
#endif

  INTEGER ix, jy, ipoint_globe, ncid 
  REAL cplwght_in(NX_GLOBE,NY_GLOBE)
  CHARACTER(LEN=31) :: routine = "MCKPP_INITIALIZE_COUPLINGWEIGHT"
  CHARACTER(LEN=max_message_len) :: message

!  If L_CPLWGHT has been set, then we will use the
!  NetCDF file to set values of cplwght over the
!  entire globe.
!  Otherwise, we will set the values ourselves, based
!  on the coupling region.
!  NPK 10/9/07 - R1
  
#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif
    CALL mckpp_print(routine, "Reading coupling weight (alpha)")
    CALL mckpp_netcdf_open(routine, kpp_const_fields%cplwght_file, ncid)
    CALL mckpp_netcdf_get_var(routine, kpp_const_fields%cplwght_file, ncid, &
        "alpha", cplwght_in)
    
#ifdef MCKPP_CAM3
    ! Use coupling weight to get global latitude and longitude grid
    CALL mckpp_netcdf_get_var(routine, kpp_const_fields%cplwght_file, ncid, &
        "latitude", kpp_global_fields%latitude)
    CALL mckpp_netcdf_get_var(routine, kpp_const_fields%cplwght_file, ncid, &
        "longitude", kpp_global_fields%longitude)
#endif

    CALL mckpp_netcdf_close(routine, kpp_const_fields%cplwght_file, ncid)
#ifdef MCKPP_CAM3  
  ENDIF ! End of masterproc section
  
  cplwght(:,:)=cplwght_in(:,:)
  CALL scatter_field_to_chunk(1,1,1,PLON,cplwght,cplwght_chunk(1,begchunk))
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%cplwght(1:ncol)=cplwght_chunk(1:ncol,ichnk)        
  ENDDO
#else
  DO ix=1,NX_GLOBE
     DO jy=1,NY_GLOBE
        ipoint_globe=(jy-1)*NX_GLOBE+ix
        kpp_3d_fields%cplwght(ipoint_globe)=cplwght_in(ix,jy)
     ENDDO
  ENDDO  
#endif
  
END SUBROUTINE mckpp_initialize_couplingweight

END MODULE mckpp_initialize_couplingweight_mod
