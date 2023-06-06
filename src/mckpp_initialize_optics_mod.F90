MODULE mckpp_initialize_optics_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, &
        mckpp_netcdf_close, &
        mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_var
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, npts

  IMPLICIT NONE

CONTAINS

  SUBROUTINE mckpp_initialize_optics()

    INTEGER :: ncid, offset_lon, offset_lat
    INTEGER, DIMENSION(2) :: start, count
    CHARACTER(LEN=max_nc_filename_len) :: file
    CHARACTER(LEN=23) :: routine = "MCKPP_INITIALIZE_OPTICS"
    CHARACTER(LEN=max_message_len) :: message

    CALL mckpp_print(routine, "")

    IF (kpp_const_fields%l_jerlov) THEN 

      file = kpp_const_fields%paras_file
      WRITE(message,*) "Reading optical properties of seawater from ", TRIM(file)
      CALL mckpp_print(routine, message)
      CALL mckpp_netcdf_open(routine, file, ncid)

      CALL mckpp_netcdf_determine_boundaries( & 
          routine, file, ncid, &
          kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), offset_lon, offset_lat)

      start = (/ offset_lon, offset_lat /) 
      count = (/ nx, ny /)
      CALL mckpp_netcdf_get_var( routine, file, ncid, "jerlov", & 
                                 kpp_3d_fields%jerlov, start, count, 2)
      CALL mckpp_netcdf_close(routine, file, ncid)

    ELSE ! no optics file 

      kpp_3d_fields%jerlov(:) = 3

    ENDIF

  END SUBROUTINE mckpp_initialize_optics

END MODULE mckpp_initialize_optics_mod
