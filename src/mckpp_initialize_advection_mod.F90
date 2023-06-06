MODULE mckpp_initialize_advection_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, &
        mckpp_netcdf_close, mckpp_netcdf_determine_boundaries, & 
        mckpp_netcdf_get_var
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, npts, maxmodeadv

  IMPLICIT NONE

CONTAINS 

  SUBROUTINE mckpp_initialize_advection()

    REAL, DIMENSION(nx, ny, 2) :: nmodeadv_temp
    REAL, DIMENSION(nx, ny, maxmodeadv, 2) :: modeadv_temp, advection_temp 
    INTEGER, DIMENSION(3) :: start, count
    INTEGER, DIMENSION(2) :: shape
    INTEGER :: i, ipt, ivar, ncid, offset_lon, offset_lat
    CHARACTER(LEN=max_nc_filename_len) :: file
    CHARACTER(LEN=26) :: routine = "MCKPP_INITIALIZE_ADVECTION"
    CHARACTER(LEN=max_message_len) :: message

    CALL mckpp_print(routine, "") 
    
    IF (kpp_const_fields%l_advect) THEN

      file = kpp_const_fields%advect_file
      WRITE(message,*) "Reading advection from file ", TRIM(file)
      CALL mckpp_print(routine, message)
      CALL mckpp_netcdf_open(routine, file, ncid)

      CALL mckpp_netcdf_determine_boundaries( & 
        routine, file, ncid,                  &
        kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), offset_lon, offset_lat)

      start = (/ offset_lon, offset_lat, 1 /)
      count = (/ nx, ny, maxmodeadv /)

      CALL mckpp_netcdf_get_var( routine, file, ncid, "nmode_tadv", & 
                                 kpp_3d_fields%nmodeadv(:,1),       & 
                                 start(1:2), count(1:2), 2) 
      CALL mckpp_netcdf_get_var( routine, file, ncid, "nmode_sadv", & 
                                 kpp_3d_fields%nmodeadv(:,2),       &
                                 start(1:2), count(1:2), 2)     
      CALL mckpp_netcdf_get_var( routine, file, ncid, "mode_tadv",  & 
                                 kpp_3d_fields%modeadv(:,:,1),      &
                                 start, count, 3) 
      CALL mckpp_netcdf_get_var( routine, file, ncid, "mode_sadv",  & 
                                 kpp_3d_fields%modeadv(:,:,2),      &
                                 start, count, 3) 
      CALL mckpp_netcdf_get_var( routine, file, ncid, "tadv",       & 
                                 kpp_3d_fields%advection(:,:,1),    & 
                                 start, count, 3) 
      CALL mckpp_netcdf_get_var( routine, file, ncid, "sadv",       & 
                                 kpp_3d_fields%advection(:,:,2),    & 
                                 start, count, 3) 

      CALL mckpp_netcdf_close(routine, file, ncid)
 
    ELSE

      CALL mckpp_print(routine, "No advection has been specified.")     
      kpp_3d_fields%nmodeadv(:,:) = 0.0 
    
    ENDIF

  END SUBROUTINE mckpp_initialize_advection

END MODULE mckpp_initialize_advection_mod
