MODULE mckpp_read_salinity_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, & 
        mckpp_netcdf_close, mckpp_netcdf_determine_boundaries, & 
        mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_parameters, ONLY: nx, ny, nzp1
  USE mckpp_time_control, ONLY: mckpp_get_update_time, time

  IMPLICIT NONE

  PUBLIC mckpp_read_salinity_3d

  PRIVATE

  LOGICAL :: l_initialized = .FALSE. 
  INTEGER :: num_times
  INTEGER, DIMENSION(4) :: start, count
  REAL, DIMENSION(:), ALLOCATABLE :: file_times

CONTAINS

  SUBROUTINE initialize_salinity(file, ncid)

    CHARACTER(LEN=max_nc_filename_len), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: ncid
    INTEGER :: offset_lon, offset_lat, nz_in
    CHARACTER(LEN=19) :: routine = "INITIALIZE_SALINITY"
    CHARACTER(LEN=max_message_len) :: message

    WRITE(message,*) "Initializing ", TRIM(file)
    CALL mckpp_print(routine, message)

    ! Work out start and count for each time entry
    CALL mckpp_netcdf_determine_boundaries( & 
      routine, file, ncid, kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), & 
      offset_lon, offset_lat, num_times )

    start = (/ offset_lon, offset_lat, 1, 1 /)
    count = (/ nx, ny, nzp1, 1 /)
 
    ! Read in time field
    ALLOCATE( file_times(num_times) ) 
    CALL mckpp_netcdf_get_var(routine, file, ncid, "t", file_times)

    ! Check vertical levels
    CALL mckpp_netcdf_get_coord(routine, file, ncid, "z", nz_in)
    IF (nz_in .NE. nzp1) THEN
      WRITE(message,*) "Salinity climatology file does not have the correct ", &
          "number of vertical levels."
      CALL mckpp_print_error(routine, message)
      WRITE(message,*) "It should have ", nzp1, " but instead has ", nz_in
      CALL mckpp_print_error(routine, message)
      CALL mckpp_abort()
    END IF

    l_initialized = .TRUE.

  END SUBROUTINE initialize_salinity


  ! Read in a NetCDF file containing a 
  ! time-varying salinity field at every model vertical level.
  ! Frequency of read is controlled by ndtupdsal in the namelist
  ! NPK 12/02/08
  SUBROUTINE mckpp_read_salinity_3d()

    CHARACTER(LEN=max_nc_filename_len) :: file
    REAL :: update_time
    INTEGER :: ncid, k
    CHARACTER(LEN=22) :: routine = "MCKPP_READ_SALINITY_3D"
    CHARACTER(LEN=max_message_len) :: message

    file = kpp_const_fields%sal_file
    CALL mckpp_netcdf_open(routine, file, ncid)

    ! On first call, get file dimensions
    IF (.NOT. l_initialized) CALL initialize_salinity(file, ncid)
 
    ! Work out time to read and check against times in file
    CALL mckpp_get_update_time( & 
      file, time, kpp_const_fields%ndtupdsal, file_times, & 
      num_times, kpp_const_fields%L_PERIODIC_SAL, kpp_const_fields%sal_period, &
      update_time, start(4), method=2 )
    WRITE(message,*) 'Reading ocean temperature for time ', update_time
    CALL mckpp_print(routine, message)
    WRITE(message,*) 'Reading ocean temperature from position ',start(4)
    CALL mckpp_print(routine, message)

    ! Read data 
    CALL mckpp_netcdf_get_var( routine, file, ncid, "salinity", & 
                               kpp_3d_fields%sal_clim, start, count, 4 ) 
    CALL mckpp_netcdf_close(routine, file, ncid)

    ! Subtract reference salinity from climatology, for compatability with
    ! salinity values stored in main model.  
    DO k = 1, nzp1
      kpp_3d_fields%sal_clim(:,k) = kpp_3d_fields%sal_clim(:,k) - & 
                                    kpp_3d_fields%sref(:)
    ENDDO

  END SUBROUTINE mckpp_read_salinity_3d

END MODULE mckpp_read_salinity_mod
