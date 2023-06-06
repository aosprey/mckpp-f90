MODULE mckpp_read_sst_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, & 
        mckpp_netcdf_close, mckpp_netcdf_determine_boundaries, & 
        mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_parameters, ONLY: nx, ny
  USE mckpp_time_control, ONLY: mckpp_get_update_time

  IMPLICIT NONE

  PUBLIC mckpp_read_sst

  PRIVATE

  LOGICAL :: l_initialized = .FALSE. 
  INTEGER :: num_times
  INTEGER, DIMENSION(3) :: start
  REAL, DIMENSION(:), ALLOCATABLE :: file_times

CONTAINS

  SUBROUTINE initialize_sst(file, ncid)

    CHARACTER(LEN=max_nc_filename_len), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: ncid
    INTEGER :: offset_lon, offset_lat
    CHARACTER(LEN=14) :: routine = "INITIALIZE_SST"
    CHARACTER(LEN=max_message_len) :: message

    WRITE(message,*) "Initializing ", TRIM(file)
    CALL mckpp_print(routine, message)

    ! Work out start and count for each time entry
    CALL mckpp_netcdf_determine_boundaries( & 
      routine, file, ncid, kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), & 
      offset_lon, offset_lat, num_times)

    start = (/ offset_lon, offset_lat, 1 /)

    ! Read in time field
    ALLOCATE(file_times(num_times)) 
    CALL mckpp_netcdf_get_var(routine, file, ncid, "t", file_times)

    l_initialized = .TRUE.

  END SUBROUTINE initialize_sst


  SUBROUTINE mckpp_read_sst()

    CHARACTER(LEN=max_nc_filename_len) :: file
    REAL :: update_time
    INTEGER :: ncid, i, j
    REAL :: offset_sst

    CHARACTER(LEN=14) :: routine = "MCKPP_READ_SST"
    CHARACTER(LEN=max_message_len) :: message

    file = kpp_const_fields%sst_file
    CALL mckpp_netcdf_open(routine, file, ncid)

    ! On first call, get file dimensions
    IF (.NOT. l_initialized) CALL initialize_sst(file, ncid)

    ! Work out time to read and check against times in file
    CALL mckpp_get_update_time( & 
      file, kpp_const_fields%time, kpp_const_fields%ndtupdsst, file_times, & 
      num_times, kpp_const_fields%l_periodic_climsst, &  
      kpp_const_fields%climsst_period, update_time, start(3), method=1 )

    WRITE(message,*) 'Reading climatological SST for time ', update_time
    CALL mckpp_print(routine, message)
    WRITE(message,*) 'Reading climatological SST from position ',start(3)
    CALL mckpp_print(routine, message)

    ! Read data 
    CALL mckpp_netcdf_get_var( routine, file, ncid, "sst", kpp_3d_fields%sst, & 
                               start) 
    CALL mckpp_netcdf_close(routine, file, ncid)

    ! KPP expects temperatures in CELSIUS.  If climatological SSTs are
    ! in Kelvin, subtract 273.15.
    IF ( ANY( kpp_3d_fields%sst .GT. 200 .AND. &
              kpp_3d_fields%sst .LT. 400) ) & 
      kpp_3d_fields%sst = kpp_3d_fields%sst - kpp_const_fields%tk0
    
    IF ( .NOT. kpp_const_fields%l_climice ) THEN
      kpp_3d_fields%iceconc = 0.0
    ENDIF

    IF ( .NOT. kpp_const_fields%l_climcurr ) THEN
      kpp_3d_fields%usf = 0.0
      kpp_3d_fields%vsf = 0.0
    ENDIF

  END SUBROUTINE mckpp_read_sst

END MODULE mckpp_read_sst_mod
