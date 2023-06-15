MODULE mckpp_read_ice_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_mpi_control, ONLY: l_root, root, mckpp_scatter_field
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, & 
        mckpp_netcdf_close, mckpp_netcdf_determine_boundaries, & 
        mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_parameters, ONLY: nx, ny, npts
  USE mckpp_time_control, ONLY: mckpp_get_update_time, time

  IMPLICIT NONE

  PUBLIC mckpp_read_ice

  PRIVATE

  LOGICAL :: l_initialized = .FALSE. 
  INTEGER :: num_times
  INTEGER, DIMENSION(3) :: start, count
  REAL, DIMENSION(:), ALLOCATABLE :: file_times

CONTAINS


  SUBROUTINE initialize_ice(file, ncid)

    CHARACTER(LEN=max_nc_filename_len), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: ncid
    INTEGER :: offset_lon, offset_lat
    CHARACTER(LEN=14) :: routine = "INITIALIZE_ICE"
    CHARACTER(LEN=max_message_len) :: message

    WRITE(message,*) "Initializing ", TRIM(file)
    CALL mckpp_print(routine, message)

    ! Work out start and count for each time entry
    CALL mckpp_netcdf_determine_boundaries( & 
      routine, file, ncid, kpp_const_fields%alon, kpp_const_fields%alat, & 
      offset_lon, offset_lat, num_times )

    start = (/ offset_lon, offset_lat, 1 /)
    count = (/ nx, ny, 1 /)

    ! Read in time field
    ALLOCATE(file_times(num_times)) 
    CALL mckpp_netcdf_get_var(routine, file, ncid, "t", file_times)

    l_initialized = .TRUE.

  END SUBROUTINE initialize_ice

  ! Read in ice concentrations from a user-provided netCDF file.
  ! Called _only_ if L_CLIMICE is TRUE in 3D_ocn.nml.
  ! Probably only necessary when coupling the model to an atmosphere
  ! that requires realistic ice concentrations.
  ! Written by Nick Klingaman, 11/01/08.
  SUBROUTINE mckpp_read_ice()

    IMPLICIT NONE

    CHARACTER(LEN=max_nc_filename_len) :: file
    REAL :: update_time
    INTEGER :: ncid
    REAL, DIMENSION(npts) :: var_global
!    REAL :: max_ice, min_ice
!    INTEGER :: i, j

    CHARACTER(LEN=14) :: routine = "MCKPP_READ_ICE"
    CHARACTER(LEN=max_message_len) :: message

    IF (l_root) THEN 
 
      file = kpp_const_fields%ice_file
      CALL mckpp_netcdf_open(routine, file, ncid)

      ! On first call, get file dimensions
      IF ( .NOT. l_initialized ) CALL initialize_ice(file, ncid)

      ! Work out time to read and check against times in file
      CALL mckpp_get_update_time( & 
        file, time, kpp_const_fields%ndtupdice, file_times, & 
        num_times, kpp_const_fields%L_periodic_climice, & 
        kpp_const_fields%climice_period,  update_time, start(3), method=1 )
   
      WRITE(message,*) 'Reading climatological ice concentration for time ', & 
                        update_time
      CALL mckpp_print(routine, message)

      CALL mckpp_netcdf_get_var( routine, file, ncid, "iceconc", & 
                                 var_global, start, count, 3 ) 

    END IF 
    CALL mckpp_scatter_field( var_global, kpp_3d_fields%iceconc, root )

!     max_ice = -1000.
!     min_ice = 1000.
!     DO j = 1, ny 
!       DO i = 1, nx
!        IF ( kpp_3d_fields%iceconc(i, j) .GT. max_ice ) & 
!          max_ice = kpp_3d_fields%iceconc(i, j)
!         IF ( kpp_3d_fields%iceconc(i, j) .LT. min_ice ) & 
!           min_ice = kpp_3d_fields%iceconc(i, j)
!       ENDDO
!     ENDDO

    IF ( kpp_const_fields%l_clim_ice_depth ) THEN
      WRITE(message,*) 'Reading climatological ice depth for time ', & 
                       update_time
      CALL mckpp_print(routine, message)

      IF (l_root) & 
        CALL mckpp_netcdf_get_var( routine, file, ncid, "icedepth", & 
                                   var_global, start, count, 3 )
      CALL mckpp_scatter_field( var_global, kpp_3d_fields%icedepth, root )
    END IF

    IF ( kpp_const_fields%l_clim_snow_on_ice ) THEN
      WRITE(message,*) 'Reading climatological snow depth for time ', & 
                       update_time
      CALL mckpp_print(routine, message)

      IF (l_root) &
        CALL mckpp_netcdf_get_var( routine, file, ncid, "snowdepth", & 
                                   var_global, start, count, 3 )
      CALL mckpp_scatter_field( var_global, kpp_3d_fields%snowdepth, root )

    ENDIF

    IF (l_root) CALL mckpp_netcdf_close(routine, file, ncid)

  END SUBROUTINE mckpp_read_ice

END MODULE mckpp_read_ice_mod
