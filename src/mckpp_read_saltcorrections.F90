MODULE mckpp_read_salt_corrections_mod 

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_parameters, ONLY: nx, ny, nx_globe, ny_globe, nzp1
  USE mckpp_time_control, ONLY: mckpp_get_update_time

  IMPLICIT NONE

  PUBLIC mckpp_read_sfcorr_2d, mckpp_read_sfcorr_3d

  PRIVATE

  LOGICAL :: l_initialized = .FALSE. 
  INTEGER :: my_nx, my_ny, num_times
  INTEGER, DIMENSION(:), ALLOCATABLE :: start, count
  REAL, DIMENSION(:), ALLOCATABLE :: file_times

CONTAINS

  SUBROUTINE initialize_sfcorr(file, ncid, ndims)

    CHARACTER(LEN=max_nc_filename_len), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: ncid, ndims

    REAL :: start_lat, start_lon
    INTEGER :: nz_in
    CHARACTER(LEN=17) :: routine = "INITIALIZE_SFCORR"
    CHARACTER(LEN=max_message_len) :: message

    ! num_dims must be 2 or 3
    IF (ndims .LT. 2 .OR. ndims .GT. 3) THEN
      WRITE(message,*) "Correction must be 2d or 3d, but num_dims = ", ndims
      CALL mckpp_print_error(routine, message)
      CALL MCKPP_ABORT
    END IF

    WRITE(message,*) "Initializing ", TRIM(file)
    CALL mckpp_print(routine, message)

    my_nx = nx
    my_ny = ny 

    ! Work out start and count for each time entry
    ALLOCATE(start(ndims+1))
    ALLOCATE(count(ndims+2))
    start = 1
    IF (ndims .EQ. 2) THEN
      count = (/my_nx,my_ny,1/)
    ELSE IF (ndims .EQ. 3) THEN 
      count = (/my_nx,my_ny,nzp1,1/)
    END IF

    start_lon = kpp_3d_fields%dlon(1)
    start_lat = kpp_3d_fields%dlat(1)
    CALL mckpp_netcdf_determine_boundaries(routine, file, ncid, &
        start_lon, start_lat, start(1), start(2), num_times)

    ! Read in time field
    ALLOCATE(file_times(num_times)) 
    CALL mckpp_netcdf_get_var(routine, file, ncid, "t", file_times)

    ! Check vertical levels
    IF (ndims .EQ. 3) THEN
      CALL mckpp_netcdf_get_coord(routine, file, ncid, "z", nz_in)
      IF (nz_in .NE. nzp1) THEN
        WRITE(message,*) "Salinity corrections file does not have the correct ", &
            "number of vertical levels."
        CALL mckpp_print_error(routine, message)
        WRITE(message,*) "It should have ", NZP1, " but instead has ", nz_in
        CALL mckpp_print_error(routine, message)
        CALL mckpp_abort()
      END IF
    END IF
    l_initialized = .TRUE.

  END SUBROUTINE initialize_sfcorr


  ! Read in a NetCDF file containing a time-varying salinity correction
  ! at the surface only.  Frequency of read is controlled by ndtupdsfcorr
  ! in the namelist
  ! NPK 29/06/08
  SUBROUTINE mckpp_read_sfcorr_2d()

    CHARACTER(LEN=max_nc_filename_len) :: file
    REAL :: update_time
    INTEGER :: ncid, ix, iy, ipt
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: var_in

    CHARACTER(LEN=20) :: routine = "MCKPP_READ_SFCORR_2D"
    CHARACTER(LEN=max_message_len) :: message

    file = kpp_const_fields%sfcorr_file
    CALL mckpp_netcdf_open(routine, file, ncid)

    ! On first call, get file dimensions
    IF (.NOT. l_initialized) CALL initialize_sfcorr(file, ncid, ndims=2)
    ALLOCATE( var_in(my_nx,my_ny,1) )

    ! Work out time to read and check against times in file
    CALL mckpp_get_update_time(file, kpp_const_fields%time, kpp_const_fields%ndtupdsfcorr, &
        file_times, num_times, kpp_const_fields%l_periodic_sfcorr, kpp_const_fields%sfcorr_period, &
        update_time, start(3), method=2)
    WRITE(message,*) 'Reading salinity correction for time ', update_time
    CALL mckpp_print(routine, message)
    WRITE(message,*) 'Reading salinity correction from position ',start(3)
    CALL mckpp_print(routine, message)

    ! Read data 
    CALL mckpp_netcdf_get_var(routine, file, ncid, "sfcorr", var_in, start) 
    CALL mckpp_netcdf_close(routine, file, ncid)

    !  Put all (NX,NY) points into one long array with dimension NPTS.         
    DO ix = 1,NX
      DO iy = 1,NY
        ipt = (iy-1)*nx+ix
        kpp_3d_fields%sfcorr_twod(ipt) = var_in(ix,iy,1)
      ENDDO
    ENDDO
    DEALLOCATE(var_in)

  END SUBROUTINE MCKPP_READ_SFCORR_2D


  ! Read in a NetCDF file containing a 
  ! time-varying salinity correction at every model vertical level.
  ! Frequency of read is controlled by ndtupdsfcorr in the namelist
  ! NPK 12/02/08  
  SUBROUTINE MCKPP_READ_SFCORR_3D()

    CHARACTER(LEN=max_nc_filename_len) :: file
    REAL :: update_time
    INTEGER :: ncid, ix, iy, iz, ipt
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: var_in

    CHARACTER(LEN=20) :: routine = "MCKPP_READ_SFCORR_3D"
    CHARACTER(LEN=max_message_len) :: message

    file = kpp_const_fields%sfcorr_file
    CALL mckpp_netcdf_open(routine, file, ncid)

    ! On first call, get file dimensions
    IF (.NOT. l_initialized) CALL initialize_sfcorr(file, ncid, ndims=3)
    ALLOCATE( var_in(my_nx,my_ny,nzp1,1) )

    ! Work out time to read and check against times in file
    CALL mckpp_get_update_time(file, kpp_const_fields%time, kpp_const_fields%ndtupdsfcorr, &
        file_times, num_times, kpp_const_fields%l_periodic_sfcorr, kpp_const_fields%sfcorr_period, &
        update_time, start(4), method=1)
    WRITE(message,*) 'Reading salinity correction for time ', update_time
    CALL mckpp_print(routine, message)
    WRITE(message,*) 'Reading salinity correction from position ',start(4)
    CALL mckpp_print(routine, message)

    ! Read data 
    CALL mckpp_netcdf_get_var(routine, file, ncid, "sfcorr", var_in, start) 
    CALL mckpp_netcdf_close(routine, file, ncid)

    ! Put all (NX,NY) points into one long array with dimension NPTS.         
    DO ix = 1,NX
      DO iy = 1,NY
        ipt = (iy-1)*nx+ix
        DO iz = 1,NZP1
          kpp_3d_fields%sfcorr_withz(ipt,iz) = var_in(ix,iy,iz,1)
        ENDDO
      ENDDO
    ENDDO

    DEALLOCATE(var_in)

  END SUBROUTINE MCKPP_READ_SFCORR_3D

END MODULE mckpp_read_salt_corrections_mod
