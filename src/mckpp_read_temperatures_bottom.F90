MODULE mckpp_read_temperatures_bottom_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, & 
       mckpp_netcdf_close, mckpp_netcdf_determine_boundaries, & 
        mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_parameters, ONLY: nx, ny
  USE mckpp_time_control, ONLY: mckpp_get_update_time

  IMPLICIT NONE

  PUBLIC mckpp_read_temperatures_bottom

  PRIVATE

  LOGICAL :: l_initialized = .FALSE. 
  INTEGER :: num_times
  INTEGER, DIMENSION(3) :: start, count
  REAL, DIMENSION(:), ALLOCATABLE :: file_times

CONTAINS

  SUBROUTINE initialize_temperatures(file, ncid)

    CHARACTER(LEN=max_nc_filename_len), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: ncid
    INTEGER :: offset_lon, offset_lat
    CHARACTER(LEN=23) :: routine = "INITIALIZE_TEMPERATURES"
    CHARACTER(LEN=max_message_len) :: message

    ! Work out start and count for each time entry
    CALL mckpp_netcdf_determine_boundaries( & 
      routine, file, ncid, kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), & 
      offset_lon, offset_lat, num_times)

    start = (/ offset_lon, offset_lat, 1 /)
    count = (/ nx, ny, 1 /)

    ! Read in time field
    ALLOCATE(file_times(num_times)) 
    CALL mckpp_netcdf_get_var(routine, file, ncid, "t", file_times)

    l_initialized = .TRUE.

  END SUBROUTINE initialize_temperatures


  SUBROUTINE mckpp_read_temperatures_bottom()

    CHARACTER(LEN=max_nc_filename_len) :: file
    REAL :: update_time, offset_temp
    INTEGER :: ncid
    CHARACTER(LEN=30) :: routine = "MCKPP_READ_TEMPERATURES_BOTTOM"
    CHARACTER(LEN=max_message_len) :: message

    file = kpp_const_fields%bottom_file
    CALL mckpp_netcdf_open(routine, file, ncid)

    ! On first call, get file dimensions
    IF (.NOT. l_initialized) CALL initialize_temperatures(file, ncid)

    ! Work out time to read and check against times in file
    CALL mckpp_get_update_time( & 
      file, kpp_const_fields%time, kpp_const_fields%ndtupdbottom, &
      file_times, num_times, kpp_const_fields%l_periodic_bottom_temp, & 
      kpp_const_fields%bottom_temp_period, update_time, start(3), method=1 )

    WRITE(message,*) 'Reading climatological bottom temp for time ', update_time
    CALL mckpp_print(routine, message)
    WRITE(message,*) 'Reading climatological bottom temp from position ', & 
                     start(3)
    CALL mckpp_print(routine, message)

    ! Read data 
    CALL mckpp_netcdf_get_var( routine, file, ncid, "T", & 
                               kpp_3d_fields%bottom_temp, start, count, 3) 

    CALL mckpp_netcdf_close(routine, file, ncid)

    ! If Kelvin, convert to Celsius
    IF ( ANY( kpp_3d_fields%bottom_temp .GT. 200 .AND. & 
              kpp_3d_fields%bottom_temp .LT. 400) ) & 
      kpp_3d_fields%bottom_temp = kpp_3d_fields%bottom_temp - & 
                                  kpp_const_fields%tk0
    
  END SUBROUTINE mckpp_read_temperatures_bottom

END MODULE mckpp_read_temperatures_bottom_mod
