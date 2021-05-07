! Routines for reading netcdf data
! With eror handling. 
MODULE mckpp_netcdf_read

  USE mckpp_log_messages, ONLY: mckpp_print_error, max_message_len
  
  USE netcdf

  IMPLICIT NONE

  PUBLIC :: max_nc_filename_len
  PUBLIC :: mckpp_netcdf_open, mckpp_netcdf_close, mckpp_netcdf_determine_boundaries, &
      mckpp_netcdf_get_coord, mckpp_netcdf_get_var

  PRIVATE 

  INTEGER, PARAMETER :: max_nc_filename_len = 200

  INTERFACE mckpp_netcdf_get_var
    MODULE PROCEDURE mckpp_netcdf_get_var_real4_1d, mckpp_netcdf_get_var_real4_3d
  END INTERFACE mckpp_netcdf_get_var

CONTAINS

  ! Open netCDF file for reading.
  ! Return ncid 
  SUBROUTINE mckpp_netcdf_open(calling_routine, file_name, ncid) 

    CHARACTER(LEN=*), INTENT(IN) ::  calling_routine, file_name
    INTEGER, INTENT(OUT) :: ncid
  
    CHARACTER(LEN=17) :: routine = "MCKPP_NETCDF_OPEN"
    CHARACTER(LEN=max_message_len) :: context, message

    context = update_context(calling_routine, routine) 
    WRITE(message, *) "Opening file ", TRIM(file_name)
    CALL check( context, message, NF90_open(file_name, NF90_nowrite, ncid) )
    
  END SUBROUTINE mckpp_netcdf_open


  ! Close netcdf file
  SUBROUTINE mckpp_netcdf_close(calling_routine, file_name, ncid)

    CHARACTER(LEN=*), INTENT(IN) ::  calling_routine, file_name
    INTEGER, INTENT(IN) :: ncid

    CHARACTER(LEN=18) :: routine = "MCKPP_NETCDF_CLOSE"
    CHARACTER(LEN=max_message_len) :: context, message
    
    context = update_context(calling_routine, routine) 
    WRITE(message, *) "Closing file ", TRIM(file_name)
    CALL check( context, message, NF90_close(ncid) )
    
  END SUBROUTINE mckpp_netcdf_close


  ! Return bounds of the dimensions, and the starting lats and lons that match model grid.
  ! Assume dimensions are named "latitude", "longitude", and "time", but these could be
  ! optional arguments. 
  SUBROUTINE mckpp_netcdf_determine_boundaries(calling_routine, file_name, ncid, &
      start_lon, start_lat, offset_lon, offset_lat, first_time, last_time) 

    CHARACTER(LEN=*), INTENT(IN) :: calling_routine, file_name
    INTEGER, INTENT(IN) :: ncid
    REAL, INTENT(IN) :: start_lon, start_lat
    INTEGER, INTENT(OUT) :: offset_lon, offset_lat
    REAL, INTENT(OUT), OPTIONAL :: first_time, last_time

    CHARACTER(LEN=max_message_len) :: context, message
    INTEGER :: nlon, nlat, ntime, lon_varid, lat_varid, time_varid
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: longitudes, latitudes

    CHARACTER(LEN=33) :: routine = "MCKPP_NETCDF_DETERMINE_BOUNDARIES"
    CHARACTER(LEN=20) :: lat_name="latitude", lon_name="longitude", time_name="time"
    REAL :: tol = 1.e-3

    context = update_context(calling_routine, routine)
    
    WRITE(message, *) "Reading longitude from file ", TRIM(file_name)
    CALL mckpp_netcdf_get_coord(context, file_name, ncid, lon_name, nlon, lon_varid)
    ALLOCATE(longitudes(nlon))
    CALL check( context, message, NF90_get_var(ncid, lon_varid, longitudes) )

    CALL get_start_position(start_lon, nlon, longitudes, tol, offset_lon) 
    IF (offset_lon .EQ. 0 ) THEN 
      CALL mckpp_print_error(context, message)
      WRITE(message,*) "Could not find starting longitude", start_lon
      CALL mckpp_print_error("", message) 
      CALL mckpp_abort()
    END IF
    
    WRITE(message, *) "Reading latitude from file ", TRIM(file_name)
    CALL mckpp_netcdf_get_coord(context, file_name, ncid, lat_name, nlat, lat_varid)
    ALLOCATE(latitudes(nlat))
    CALL check( context, message, NF90_get_var(ncid, lat_varid, latitudes) )

    CALL get_start_position(start_lat, nlat, latitudes, tol, offset_lat) 
    IF (offset_lat .EQ. 0 ) THEN 
      CALL mckpp_print_error(context, message) 
      WRITE(message,*) "Could not find starting latitude", start_lat
      CALL mckpp_print_error("", message) 
      CALL mckpp_abort()
    END IF

    IF ( PRESENT(first_time) .AND. PRESENT(last_time) ) THEN 
      WRITE(message, *) "Reading time from file ", TRIM(file_name)
      CALL mckpp_netcdf_get_coord(context, file_name, ncid, time_name, ntime, time_varid)
      CALL check( context, message, NF90_get_var(ncid, time_varid, first_time, (/1/)) )
      CALL check( context, message, NF90_get_var(ncid, time_varid, last_time, (/ntime/)) )
    END IF
    
  END SUBROUTINE mckpp_netcdf_determine_boundaries


  ! Return size and varid for a coordinate variable. 
  SUBROUTINE mckpp_netcdf_get_coord(calling_routine, file_name, ncid, dim_name, dim_len, dim_varid)

    CHARACTER(LEN=*), INTENT(IN) :: calling_routine, file_name, dim_name
    INTEGER, INTENT(IN) :: ncid
    INTEGER, INTENT(OUT) :: dim_len
    INTEGER, INTENT(OUT), OPTIONAL :: dim_varid

    INTEGER :: dimid, varid
    CHARACTER(LEN=max_message_len) :: context, message
    CHARACTER(LEN=22) :: routine = "MCKPP_NETCDF_GET_COORD"

    context = update_context(calling_routine, routine)
    WRITE(message, *) "Reading ", TRIM(dim_name), " from file ", TRIM(file_name)
    
    CALL check( context, message, NF90_inq_dimid(ncid, dim_name, dimid) )
    CALL check( context, message, NF90_inquire_dimension(ncid, dimid, len=dim_len) )
    CALL check( context, message, NF90_inq_varid(ncid, dim_name, varid) )

    IF (PRESENT(dim_varid)) dim_varid = varid
   
  END SUBROUTINE mckpp_netcdf_get_coord


  ! 4-byte real with 1 dimension
  SUBROUTINE mckpp_netcdf_get_var_real4_1d(calling_routine, file_name, ncid, var_name, array)

    CHARACTER(LEN=*), INTENT(IN) :: calling_routine, file_name, var_name
    INTEGER, INTENT(IN) :: ncid
    REAL(KIND=4), DIMENSION(:), INTENT(OUT) :: array

    INTEGER :: varid
    CHARACTER(LEN=max_message_len) :: context, message
    CHARACTER(LEN=20) :: routine = "MCKPP_NETCDF_GET_VAR"

    context = update_context(calling_routine, routine)
    WRITE(message, *) "Reading ", TRIM(var_name), " from file ", TRIM(file_name)

    CALL check( context, message, NF90_inq_varid(ncid, var_name, varid) )
    CALL check( context, message, NF90_get_var(ncid, varid, array) ) 

  END SUBROUTINE mckpp_netcdf_get_var_real4_1d

  
  ! Read variable
  ! 4-byte real with 3 dimension
  ! Start and count not optional. 
  SUBROUTINE mckpp_netcdf_get_var_real4_3d(calling_routine, file_name, &
      ncid, var_name, array, start, count)

    CHARACTER(LEN=*), INTENT(IN) :: calling_routine, file_name, var_name
    INTEGER, INTENT(IN) :: ncid
    REAL(KIND=4), DIMENSION(:,:,:), INTENT(OUT) :: array
    INTEGER, DIMENSION(3), INTENT(IN) :: start, count

    INTEGER :: varid
    CHARACTER(LEN=max_message_len) :: context, message
    CHARACTER(LEN=20) :: routine = "MCKPP_NETCDF_GET_VAR"

    context = update_context(calling_routine, routine)
    WRITE(message, *) "Reading ", TRIM(var_name), " from file ", TRIM(file_name)

    CALL check( context, message, NF90_inq_varid(ncid, var_name, varid) )
    CALL check( context, message, NF90_get_var(ncid, varid, array, start, count) ) 

  END SUBROUTINE mckpp_netcdf_get_var_real4_3d
  

  ! Find first point in array that is within tolerance of a given target value,
  ! and return array position. 
  ! If none found, return 0
  SUBROUTINE get_start_position(val, size, array, tolerance, position)
    
    REAL, INTENT(IN) :: val, tolerance
    INTEGER, INTENT(IN) :: size
    REAL(KIND=4), DIMENSION(size) :: array
    INTEGER, INTENT(OUT) :: position

    INTEGER :: n
    
    position = 0 
    DO n = 1, size
      IF ( ABS(array(n) - val) .LT. tolerance ) THEN ! match 
        position = n
        EXIT
      END IF    
    ENDDO
    
  END SUBROUTINE get_start_position


  ! Update call tree with current routine
  FUNCTION update_context(calling_routine, routine) RESULT(context)

    CHARACTER(LEN=*), INTENT(IN) :: calling_routine, routine
    CHARACTER(LEN=max_message_len) :: context
    
    WRITE(context,*) TRIM(ADJUSTL(calling_routine)), " -> ", TRIM(ADJUSTL(routine))

  END FUNCTION update_context


  ! Check netcdf routine error codes.
  ! Print error message, with some context from calling routines.
  ! Then abort
  SUBROUTINE check(routine, message, status)
    
    CHARACTER(LEN=*), INTENT(IN) :: routine, message
    INTEGER, INTENT(IN) :: status
 
    IF (status /= NF90_noerr) THEN 
      CALL mckpp_print_error(routine, message) 
      CALL mckpp_print_error("",NF90_strerror(status))
      CALL mckpp_abort()
    END IF 
    
  END SUBROUTINE check 

END MODULE mckpp_netcdf_read
