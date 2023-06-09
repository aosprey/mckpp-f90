MODULE mckpp_initialize_coupling_weight_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_netcdf_read, ONLY: mckpp_netcdf_open, mckpp_netcdf_close, &
        max_nc_filename_len, mckpp_netcdf_determine_boundaries, & 
        mckpp_netcdf_get_var
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_mpi_control, ONLY: mckpp_scatter_field, l_root, root
  USE mckpp_parameters, ONLY: nx, ny, npts

  IMPLICIT NONE

CONTAINS

  ! ** We currently don't have any coupling **
  ! Just read in over KPP domain for now 
  ! - as we are using cplwght as a local variable
  ! 
  ! If l_cplwght has been set, then we will use the NetCDF file to set
  ! values of cplwght over the entire globe.
  ! Otherwise, we will set the values ourselves, based on the coupling region.
  ! NPK 10/9/07 - R1
  SUBROUTINE mckpp_initialize_coupling_weight()

    INTEGER :: ncid, offset_lat, offset_lon
    INTEGER, DIMENSION(2) :: start, count 
    REAL, DIMENSION(npts) :: cplwght_global
    CHARACTER(LEN=max_nc_filename_len) :: file
    CHARACTER(LEN=31) :: routine = "MCKPP_INITIALIZE_COUPLINGWEIGHT"
    CHARACTER(LEN=max_message_len) :: message

    CALL mckpp_print(routine, "")

    IF (kpp_const_fields%l_couple .OR. kpp_const_fields%l_cplwght) THEN  

      IF (l_root) THEN 

        file = kpp_const_fields%cplwght_file
        WRITE(message,*) "Reading coupling weight (alpha) from file ", TRIM(file)
        CALL mckpp_print(routine, message)
        CALL mckpp_netcdf_open(routine, file, ncid)

        CALL mckpp_netcdf_determine_boundaries( & 
          routine, file, ncid, kpp_const_fields%alon, kpp_const_fields%alat, & 
          offset_lon, offset_lat )

        start = (/ offset_lon, offset_lat /) 
        count = (/ nx,  ny/)
        CALL mckpp_netcdf_get_var( routine, file, ncid, "alpha", & 
                                   cplwght_global, start, count, 2 )

        CALL mckpp_netcdf_close(routine, file, ncid)

      END IF
  
      CALL mckpp_scatter_field(cplwght_global, kpp_3d_fields%cplwght, root)

    ELSE
      kpp_3d_fields%cplwght(:) = 0.0
    END IF

  END SUBROUTINE mckpp_initialize_coupling_weight

END MODULE mckpp_initialize_coupling_weight_mod
