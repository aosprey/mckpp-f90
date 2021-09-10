MODULE mckpp_initialize_landsea_mod

  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_netcdf_read, ONLY: mckpp_netcdf_open, mckpp_netcdf_close, &
        mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_var, & 
        max_nc_filename_len
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, & 
        max_message_len, nupe 
  USE mckpp_mpi_control, ONLY: mckpp_broadcast_field, mckpp_scatter_field, &
        l_root_proc, root_proc, npts_local, offset_global, start_global, & 
        end_global 
  USE mckpp_parameters, ONLY: npts, nx, ny

  IMPLICIT NONE

CONTAINS

  SUBROUTINE mckpp_initialize_landsea()

    INTEGER :: ncid, offset_lon, offset_lat, i, j, ipt
    INTEGER, DIMENSION(2) :: start, count
    REAL, DIMENSION(nx) :: lon_in
    REAL, DIMENSION(ny) :: lat_in
    REAL, DIMENSION(npts) :: landsea_global, ocedepth_global 
    REAL, DIMENSION(npts_local) :: landsea
    CHARACTER(LEN=max_nc_filename_len) :: file
    CHARACTER(LEN=24) :: routine = "MCKPP_INITIALIZE_LANDSEA"
    CHARACTER(LEN=max_message_len) :: message

    CALL mckpp_print(routine, "")

    ! Read from file 
    IF (kpp_const_fields%l_landsea) THEN

      IF (l_root_proc) THEN

        file = kpp_const_fields%landsea_file
        WRITE(message,*) "Reading land sea mask from file ", TRIM(file)
        CALL mckpp_print(routine, message)
        CALL mckpp_netcdf_open(routine, file, ncid)

        CALL mckpp_netcdf_determine_boundaries( & 
          routine, file, ncid, kpp_const_fields%alon, kpp_const_fields%alat, & 
          offset_lon, offset_lat)

        CALL mckpp_netcdf_get_var( routine, file, ncid, "longitude", & 
                                   lon_in, offset_lon )
        CALL mckpp_netcdf_get_var( routine, file, ncid, "latitude", & 
                                   lat_in, offset_lat )

        start = (/ offset_lon, offset_lat /)
        count = (/ nx, ny /) 
        CALL mckpp_netcdf_get_var( routine, file, ncid, "lsm", & 
                                   landsea_global, start, count, 2 )
        CALL mckpp_netcdf_get_var( routine, file, ncid, "max_depth", & 
                                   ocdepth_global, start, count, 2 )       
        CALL mckpp_netcdf_close(routine, file, ncid)

      ENDIF 

      ! Broadcast lons and lats 
      CALL mckpp_broadcast_field(lon_in, nx, root_proc)
      CALL mckpp_broadcast_field(lat_in, ny, root_proc) 

      ! Expand to full grid
      ipt = 1
      DO j = 1, ny
        DO i = 1, nx
          kpp_3d_fields%dlon_all(ipt) = lon_in(i)
          kpp_3d_fields%dlat_all(ipt) = lat_in(j)
          ipt = ipt + 1 
        END DO
      END DO

      ! Scatter ocdepth and lsm 
      CALL mckpp_scatter_field(ocdepth_global, kpp_3d_fields%ocdepth, root_proc)
      CALL mckpp_scatter_field(landsea_global, landsea, root_proc)

      ! Generate logical lsm 
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
      kpp_3d_fields%L_OCEAN = .TRUE.

    ! Not method specified
    ELSE

      WRITE(message,*) "If you set L_REGGRID=.FALSE., you must specify a ", & 
                       "land-sea mask file from which to read the locations ", &
                       "of the gridpoints in the horizontal."
      CALL mckpp_abort(routine, message)
 
    ENDIF

    ! Get local lat, lons 
    kpp_3d_fields%dlat = kpp_3d_fields%dlat_all(start_global:end_global) 
    kpp_3d_fields%dlon = kpp_3d_fields%dlon_all(start_global:end_global) 

    WRITE(nupe,*) "kpp_3d_fields%dlat_all = ", kpp_3d_fields%dlat_all
    WRITE(nupe,*) "kpp_3d_fields%dlon_all = ", kpp_3d_fields%dlon_all
    WRITE(nupe,*) "kpp_3d_fields%dlat = ", kpp_3d_fields%dlat
    WRITE(nupe,*) "kpp_3d_fields%dlon = ", kpp_3d_fields%dlon
    WRITE(nupe,*) "kpp_3d_fields%l_ocean = ", kpp_3d_fields%l_ocean
    WRITE(nupe,*) "kpp_3d_fields%ocdepth = ", kpp_3d_fields%ocdepth

  END SUBROUTINE mckpp_initialize_landsea

END MODULE mckpp_initialize_landsea_mod

