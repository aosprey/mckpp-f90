MODULE mckpp_read_temperatures_bottom_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_parameters, ONLY: nx, ny, nx_globe, ny_globe, nzp1
  USE mckpp_time_control, ONLY: mckpp_get_update_time

  IMPLICIT NONE

  PUBLIC mckpp_read_temperatures_bottom

  PRIVATE

  LOGICAL :: l_initialized = .FALSE. 
  INTEGER :: my_nx, my_ny, num_times
  INTEGER, DIMENSION(3) :: start, count
  REAL, DIMENSION(:), ALLOCATABLE :: file_times

CONTAINS

  SUBROUTINE initialize_temperatures(file, ncid)

    CHARACTER(LEN=max_nc_filename_len), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: ncid

    REAL :: start_lat, start_lon
    CHARACTER(LEN=23) :: routine = "INITIALIZE_TEMPERATURES"
    CHARACTER(LEN=max_message_len) :: message

    my_nx = nx
    my_ny = ny 

    ! Work out start and count for each time entry
    count = (/my_nx,my_ny,1/)
    start = (/1,1,1/)
    start_lon = kpp_3d_fields%dlon(1)
    start_lat = kpp_3d_fields%dlat(1)
    CALL mckpp_netcdf_determine_boundaries(routine, file, ncid, &
        start_lon, start_lat, start(1), start(2), num_times)

    ! Read in time field
    ALLOCATE(file_times(num_times)) 
    CALL mckpp_netcdf_get_var(routine, file, ncid, "t", file_times)

    l_initialized = .TRUE.

  END SUBROUTINE initialize_temperatures

  SUBROUTINE mckpp_read_temperatures_bottom()

    CHARACTER(LEN=max_nc_filename_len) :: file
    REAL :: update_time, offset_temp
    INTEGER :: ncid, ix, iy, ipt
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: var_in

    CHARACTER(LEN=30) :: routine = "MCKPP_READ_TEMPERATURES_BOTTOM"
    CHARACTER(LEN=max_message_len) :: message

    file = kpp_const_fields%bottom_file
    CALL mckpp_netcdf_open(routine, file, ncid)

    ! On first call, get file dimensions
    IF (.NOT. l_initialized) CALL initialize_temperatures(file, ncid)
    ALLOCATE( var_in(my_nx,my_ny,1) )

    ! Work out time to read and check against times in file
    CALL mckpp_get_update_time(file, kpp_const_fields%time, kpp_const_fields%ndtupdbottom, &
        file_times, num_times, kpp_const_fields%L_PERIODIC_BOTTOM_TEMP, kpp_const_fields%bottom_temp_period, &
        update_time, start(3), method=1)
    WRITE(message,*) 'Reading climatological bottom temp for time ', update_time
    CALL mckpp_print(routine, message)
    WRITE(message,*) 'Reading climatological bottom temp from position ', start(3)
    CALL mckpp_print(routine, message)

    ! Read data 
    CALL mckpp_netcdf_get_var(routine, file, ncid, "T", var_in, start) 
    CALL mckpp_netcdf_close(routine, file, ncid)

    offset_temp = 0.
    ix = 1
    iy = 1
    DO WHILE (offset_temp .EQ. 0 .AND. ix .LE. my_nx)
      DO iy = 1,my_ny
        IF (var_in(ix,iy,1) .gt. 200 .and. var_in(ix,iy,1) .lt. 400) offset_temp = 273.15    
      END DO
      ix = ix+1
    ENDDO

    DO ix = 1,NX
      DO iy = 1,NY
        ipt = (iy-1)*NX+ix
        kpp_3d_fields%bottom_temp(ipt) = var_in(ix,iy,1)-offset_temp        
      ENDDO
    ENDDO

  END SUBROUTINE mckpp_read_temperatures_bottom

END MODULE mckpp_read_temperatures_bottom_mod
