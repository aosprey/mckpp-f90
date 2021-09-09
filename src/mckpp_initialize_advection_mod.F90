MODULE mckpp_initialize_advection_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_var
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, npts, maxmodeadv

  IMPLICIT NONE

CONTAINS 

  SUBROUTINE mckpp_initialize_advection()

    REAL, DIMENSION(nx, ny, 2) :: nmodeadv_temp
    REAL, DIMENSION(nx, ny, maxmodeadv, 2) :: modeadv_temp, advection_temp 
    INTEGER, DIMENSION(3) :: start
    INTEGER, DIMENSION(2) :: shape
    INTEGER :: i, ipt, ivar, ncid, offset_lon, offset_lat
    CHARACTER(LEN=max_nc_filename_len) :: file
    CHARACTER(LEN=26) :: routine = "MCKPP_INITIALIZE_ADVECTION"
    CHARACTER(LEN=max_message_len) :: message

    IF (kpp_const_fields%L_ADVECT) THEN

      file = kpp_const_fields%advect_file
      WRITE(message,*) "Reading advection from file ", TRIM(file)
      CALL mckpp_print(routine, message)
      CALL mckpp_netcdf_open(routine, file, ncid)

      CALL mckpp_netcdf_determine_boundaries(routine, file, ncid, &
          kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), offset_lon, offset_lat)
      start(1) = offset_lon
      start(2) = offset_lat
      start(3) = 1

      CALL mckpp_netcdf_get_var(routine, file, ncid, "nmode_tadv", nmodeadv_temp(:,:,1), start(1:2)) 
      CALL mckpp_netcdf_get_var(routine, file, ncid, "nmode_sadv", nmodeadv_temp(:,:,2), start(1:2))
      CALL mckpp_netcdf_get_var(routine, file, ncid, "mode_tadv", modeadv_temp(:,:,:,1), start) 
      CALL mckpp_netcdf_get_var(routine, file, ncid, "mode_sadv", modeadv_temp(:,:,:,2), start)
      CALL mckpp_netcdf_get_var(routine, file, ncid, "tadv", advection_temp(:,:,:,1), start) 
      CALL mckpp_netcdf_get_var(routine, file, ncid, "sadv", advection_temp(:,:,:,2), start)

      CALL mckpp_netcdf_close(routine, file, ncid)
      shape(1) = npts
      shape(2) = maxmodeadv
      DO ivar = 1, 2 
        kpp_3d_fields%nmodeadv(:,ivar) = RESHAPE(nmodeadv_temp(:,:,ivar), shape(1:1)) 
        kpp_3d_fields%modeadv(:,:,ivar) = RESHAPE(modeadv_temp(:,:,:,ivar), shape) 
        kpp_3d_fields%advection(:,:,ivar) = RESHAPE(advection_temp(:,:,:,ivar), shape)
      END DO

    ELSE
      CALL mckpp_print(routine, "No advection has been specified.")     
      DO ipt = 1, npts
        DO ivar = 1, 2
          kpp_3d_fields%nmodeadv(ipt,ivar)=0
        ENDDO
      ENDDO
    ENDIF

  END SUBROUTINE mckpp_initialize_advection

END MODULE mckpp_initialize_advection_mod
