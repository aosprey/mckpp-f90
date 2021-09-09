MODULE mckpp_initialize_couplingweight_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_netcdf_read, ONLY: mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_get_var, max_nc_filename_len
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: nx_globe, ny_globe

  IMPLICIT NONE

CONTAINS

  SUBROUTINE mckpp_initialize_couplingweight()

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

    file = kpp_const_fields%cplwght_file
    WRITE(message,*) "Reading coupling weight (alpha) from file ", TRIM(file)
    CALL mckpp_print(routine, message)
    CALL mckpp_netcdf_open(routine, file, ncid)
    CALL mckpp_netcdf_get_var(routine, file, ncid, "alpha", cplwght_in)
    CALL mckpp_netcdf_close(routine, file, ncid)
    DO ix=1,NX_GLOBE
       DO jy=1,NY_GLOBE
          ipoint_globe=(jy-1)*NX_GLOBE+ix
          kpp_3d_fields%cplwght(ipoint_globe)=cplwght_in(ix,jy)
       ENDDO
    ENDDO

  END SUBROUTINE mckpp_initialize_couplingweight

END MODULE mckpp_initialize_couplingweight_mod
