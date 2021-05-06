! Routines for reading netcdf data 
MODULE mckpp_netcdf_read

  USE mckpp_log_messages, ONLY: mckpp_print_error, max_message_len
  
  USE netcdf

  IMPLICIT NONE

  PUBLIC :: max_nc_filename_len
  PUBLIC :: mckpp_netcdf_open, mckpp_netcdf_close, mckpp_netcdf_determine_boundaries

  PRIVATE 

  INTEGER, PARAMETER :: max_nc_filename_len = 200 

CONTAINS

  ! Open netCDF file for reading.
  ! Return ncid 
  SUBROUTINE mckpp_netcdf_open(calling_routine, filename, ncid) 

    CHARACTER(LEN=*), INTENT(IN) ::  calling_routine
    CHARACTER(LEN=max_nc_filename_len), INTENT(IN) :: filename
    INTEGER, INTENT(OUT) :: ncid
  
    CHARACTER(LEN=17) :: routine = "MCKPP_NETCDF_OPEN"
    CHARACTER(LEN=max_message_len) :: context, message

    WRITE(context, *) calling_routine, " -> ", routine
    WRITE(message, *) "Opening file ", TRIM(filename)
    CALL check( context, message, NF90_open(filename, NF90_nowrite, ncid) )
    
  END SUBROUTINE mckpp_netcdf_open


  ! Close netcdf file
  SUBROUTINE mckpp_netcdf_close(calling_routine, filename, ncid)

    CHARACTER(LEN=*), INTENT(IN) ::  calling_routine
    CHARACTER(LEN=max_nc_filename_len), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: ncid

    CHARACTER(LEN=18) :: routine = "MCKPP_NETCDF_CLOSE"
    CHARACTER(LEN=max_message_len) :: context, message
    
    WRITE(context, *) calling_routine, " -> ", routine
    WRITE(message, *) "Closing file ", TRIM(filename)
    CALL check( context, message, NF90_close(ncid) )
    
  END SUBROUTINE mckpp_netcdf_close


  ! Return bounds of the dimensions, and the starting lats and lons that match model grid.
  ! Assume dimensions are named "latitude", "longitude", and "time", but these could be
  ! optional arguments. 
  SUBROUTINE mckpp_netcdf_determine_boundaries(calling_routine, filename, ncid, &
      start_lon, start_lat, offset_lon, offset_lat, first_time, last_time) 

    CHARACTER(LEN=*), INTENT(IN) ::  calling_routine
    CHARACTER(LEN=max_nc_filename_len), INTENT(IN) :: filename
    INTEGER, INTENT(IN) :: ncid
    REAL, INTENT(IN) :: start_lon, start_lat
    INTEGER, INTENT(OUT) :: offset_lon, offset_lat
    REAL, INTENT(OUT), OPTIONAL :: first_time, last_time

    CHARACTER(LEN=max_message_len) :: context, message
    INTEGER :: nlon, nlat, ntime, lon_varid, lat_varid, time_varid
    REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: longitudes, latitudes

    CHARACTER(LEN=33) :: routine = "MCKPP_DETERMINE_NETCDF_BOUNDARIES"
    CHARACTER(LEN=20) :: lat_name="latitude", lon_name="longitude", time_name="time"
    REAL :: tol = 1.e-3

    WRITE(context, *) calling_routine, " -> ", routine

    WRITE(message, *) "Reading longitude from file ", TRIM(filename)
    CALL get_coord(context, message, ncid, lon_name, nlon, lon_varid)
    ALLOCATE(longitudes(nlon))
    CALL check( context, message, NF90_get_var(ncid, lon_varid, longitudes) )

    CALL get_start_position(start_lon, nlon, longitudes, tol, offset_lon) 
    IF (offset_lon .EQ. 0 ) THEN 
      CALL mckpp_print_error(context, message) 
      CALL mckpp_print_error("", "Could not find starting longitude") 
      CALL mckpp_abort()
    END IF
    
    WRITE(message, *) "Reading latitude from file ", TRIM(filename)
    CALL get_coord(context, message, ncid, lat_name, nlat, lat_varid)
    ALLOCATE(latitudes(nlat))
    CALL check( context, message, NF90_get_var(ncid, lat_varid, latitudes) )

    CALL get_start_position(start_lat, nlat, latitudes, tol, offset_lat) 
    IF (offset_lat .EQ. 0 ) THEN 
      CALL mckpp_print_error(context, message) 
      CALL mckpp_print_error("", "Could not find starting latitude") 
      CALL mckpp_abort()
    END IF

    IF ( PRESENT(first_time) .AND. PRESENT(last_time) ) THEN 
      WRITE(message, *) "Reading time from file ", TRIM(filename)
      CALL get_coord(context, message, ncid, time_name, ntime, time_varid)
      CALL check( context, message, NF90_get_var(ncid, time_varid, first_time, (/1/)) )
      CALL check( context, message, NF90_get_var(ncid, time_varid, last_time, (/ntime/)) )
    END IF
    
  END SUBROUTINE mckpp_netcdf_determine_boundaries


  ! Return size and varid for a coordinate variable. 
  SUBROUTINE get_coord(context, message, ncid, dim_name, dim_len, dim_varid)

    CHARACTER(LEN=max_message_len), INTENT(IN) :: context, message
    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: dim_name
    INTEGER, INTENT(OUT) :: dim_len, dim_varid

    INTEGER :: dimid

    CALL check( context, message, NF90_inq_dimid(ncid, dim_name, dimid) )
    CALL check( context, message, NF90_inquire_dimension(ncid, dimid, len=dim_len) )
    CALL check( context, message, NF90_inq_varid(ncid, dim_name, dim_varid) )
   
  END SUBROUTINE get_coord


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

  END SUBROUTINE 


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
