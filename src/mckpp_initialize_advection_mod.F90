MODULE mckpp_initialize_advection_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, &
        mckpp_netcdf_close, mckpp_netcdf_determine_boundaries, & 
        mckpp_netcdf_get_var
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len, nupe
  USE mckpp_mpi_control, ONLY: l_root, root, mckpp_scatter_field
  USE mckpp_parameters, ONLY: nx, ny, npts, maxmodeadv

  IMPLICIT NONE

CONTAINS 

  SUBROUTINE mckpp_initialize_advection()

    INTEGER, DIMENSION(npts, 2) :: nmodeadv_tmp
    INTEGER, DIMENSION(npts, maxmodeadv, 2) :: modeadv_tmp
    REAL, DIMENSION(npts, maxmodeadv, 2) :: advection_tmp 
    INTEGER, DIMENSION(3) :: start, count
    INTEGER, DIMENSION(2) :: shape
    INTEGER :: m, p, ncid, offset_lon, offset_lat
    CHARACTER(LEN=max_nc_filename_len) :: file
    CHARACTER(LEN=26) :: routine = "MCKPP_INITIALIZE_ADVECTION"
    CHARACTER(LEN=max_message_len) :: message

    CALL mckpp_print(routine, "") 

    IF (kpp_const_fields%l_advect) THEN

      IF (l_root) THEN 

        file = kpp_const_fields%advect_file
        WRITE(message,*) "Reading advection from file ", TRIM(file)
        CALL mckpp_print(routine, message)
        CALL mckpp_netcdf_open(routine, file, ncid)

        CALL mckpp_netcdf_determine_boundaries( & 
          routine, file, ncid, kpp_const_fields%alon, kpp_const_fields%alat, & 
          offset_lon, offset_lat )

        start = (/ offset_lon, offset_lat, 1 /)
        count = (/ nx, ny, maxmodeadv /)

        CALL mckpp_netcdf_get_var( routine, file, ncid, "nmode_tadv", & 
                                   nmodeadv_tmp(:,1), start(1:2), count(1:2), 2 )
        CALL mckpp_netcdf_get_var( routine, file, ncid, "nmode_sadv", & 
                                   nmodeadv_tmp(:,2), start(1:2), count(1:2), 2 )
        CALL mckpp_netcdf_get_var( routine, file, ncid, "mode_tadv",  & 
                                   modeadv_tmp(:,:,1), start, count, 3) 
        CALL mckpp_netcdf_get_var( routine, file, ncid, "mode_sadv",  & 
                                   modeadv_tmp(:,:,2), start, count, 3) 
        CALL mckpp_netcdf_get_var( routine, file, ncid, "tadv", & 
                                   advection_tmp(:,:,1), start, count, 3) 
        CALL mckpp_netcdf_get_var( routine, file, ncid, "sadv", & 
                                   advection_tmp(:,:,2), start, count, 3) 

        CALL mckpp_netcdf_close(routine, file, ncid)

      END IF 

      ! We need to scatter 1 layer at a time (for now) 
      DO p = 1, 2
        CALL mckpp_scatter_field( nmodeadv_tmp(:,p), &
                                  kpp_3d_fields%nmodeadv(:,p), root )
        DO m = 1, maxmodeadv
          CALL mckpp_scatter_field( modeadv_tmp(:,m,p), & 
                                    kpp_3d_fields%modeadv(:,m,p), root )
          CALL mckpp_scatter_field( advection_tmp(:,m,p), & 
                                    kpp_3d_fields%advection(:,m,p), root )
        END DO 
      END DO 

      WRITE(nupe,*) "kpp_3d_fields%nmodeadv(:,1) = ", & 
                     kpp_3d_fields%nmodeadv(:,1)   
      WRITE(nupe,*) "kpp_3d_fields%nmodeadv(:,2) = ", & 
                     kpp_3d_fields%nmodeadv(:,2) 
      WRITE(nupe,*) "kpp_3d_fields%modeadv(:,:,1) = ", & 
                     kpp_3d_fields%modeadv(:,:,1) 
      WRITE(nupe,*) "kpp_3d_fields%modeadv(:,:,2) = ", & 
                     kpp_3d_fields%modeadv(:,:,2) 
      WRITE(nupe,*) "kpp_3d_fields%advection(:,:,1) = ", & 
                     kpp_3d_fields%advection(:,:,1)
      WRITE(nupe,*) "kpp_3d_fields%advection(:,:,2) = ", & 
                     kpp_3d_fields%advection(:,:,2)   
    ELSE

      CALL mckpp_print(routine, "No advection has been specified.")     
      kpp_3d_fields%nmodeadv(:,:) = 0.0 

    ENDIF

  END SUBROUTINE mckpp_initialize_advection

END MODULE mckpp_initialize_advection_mod
