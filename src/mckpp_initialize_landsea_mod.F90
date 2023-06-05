MODULE mckpp_initialize_landsea_mod

  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_netcdf_read, ONLY: mckpp_netcdf_open, mckpp_netcdf_close, &
        mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_var, & 
        max_nc_filename_len
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_parameters, ONLY: npts, nx, ny

  IMPLICIT NONE

CONTAINS

  SUBROUTINE mckpp_initialize_landsea()

    INTEGER :: ncid, offset_lon, offset_lat, i, j, ipt
    INTEGER, DIMENSION(2) :: start, count
    REAL, DIMENSION(nx) :: lon_in
    REAL, DIMENSION(ny) :: lat_in
    REAL, DIMENSION(npts) :: landsea
    CHARACTER(LEN=max_nc_filename_len) :: file
    CHARACTER(LEN=24) :: routine = "MCKPP_INITIALIZE_LANDSEA"
    CHARACTER(LEN=max_message_len) :: message

    CALL mckpp_print(routine, "")

    kpp_3d_fields%dlon(1) = kpp_const_fields%alon
    kpp_3d_fields%dlat(1) = kpp_const_fields%alat

    ! Read from file 
    IF (kpp_const_fields%l_landsea) THEN
      
      file = kpp_const_fields%landsea_file
      WRITE(message,*) "Reading land sea mask from file ", TRIM(file)
      CALL mckpp_print(routine, message)
      CALL mckpp_netcdf_open(routine, file, ncid)
      
      CALL mckpp_netcdf_determine_boundaries( & 
        routine, file, ncid, & 
        kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), offset_lon, offset_lat)

      CALL mckpp_netcdf_get_var( routine, file, ncid, "longitude", & 
                                 lon_in, offset_lon)
      CALL mckpp_netcdf_get_var( routine, file, ncid, "latitude", & 
                                 lat_in, offset_lat)

      start = (/ offset_lon, offset_lat /)
      count = (/ nx, ny /) 
      CALL mckpp_netcdf_get_var( routine, file, ncid, "lsm", landsea, & 
                                 start, count, 2)
      CALL mckpp_netcdf_get_var( routine, file, ncid, "max_depth", & 
                                 kpp_3d_fields%ocdepth, start, count, 2)       
      CALL mckpp_netcdf_close(routine, file, ncid)

      ! Expand lons and lats to full grid
      ipt = 1
      DO j = 1, ny
        DO i = 1, nx
          kpp_3d_fields%dlon(ipt) = lon_in(i)
          kpp_3d_fields%dlat(ipt) = lat_in(j)
          ipt = ipt + 1 
        END DO
      END DO
      
      ! Generate logical land sea mask
      kpp_3d_fields%l_ocean = landsea .NE. 1.0 

    ! Generate a regularly spaced grid  
    ELSEIF (kpp_const_fields%l_reggrid) THEN    

      ipt = 1
      DO j = 1, ny
        DO i = 1, nx
          kpp_3d_fields%dlon(ipt) = kpp_const_fields%alon + (i-1) * & 
                                    kpp_const_fields%delta_lon
          kpp_3d_fields%dlat(ipt) = kpp_const_fields%alat + (j-1) * & 
                                    kpp_const_fields%delta_lat
          ipt = ipt + 1 
        ENDDO
      ENDDO

      kpp_3d_fields%ocdepth = -10000.
      kpp_3d_fields%l_ocean = .TRUE.

    ! Catch all other cases 
    ELSE

      WRITE(message,*) "If you set l_reggrid = .FALSE., you must specify a ", & 
                       "land-sea mask file from which to read the locations ", &
                       "of the gridpoints in the horizontal."
      CALL mckpp_print_error(routine, message)
      CALL mckpp_abort()

    ENDIF

  END SUBROUTINE mckpp_initialize_landsea

END MODULE mckpp_initialize_landsea_mod

