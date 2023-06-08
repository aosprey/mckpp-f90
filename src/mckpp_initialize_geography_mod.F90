MODULE mckpp_initialize_geography_mod

  USE mckpp_data_fields, ONLY: kpp_const_fields, kpp_3d_fields
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len, nupe
  USE mckpp_mpi_control, ONLY: mckpp_broadcast_field, l_root, root, &
        npts_local
  USE mckpp_netcdf_read, ONLY: mckpp_netcdf_open, mckpp_netcdf_close, & 
        mckpp_netcdf_get_var, max_nc_filename_len
  USE mckpp_parameters, ONLY: nz, nzp1, npts

  IMPLICIT NONE

CONTAINS

  SUBROUTINE mckpp_initialize_geography()

    INTEGER :: i, ipt, ncid
    REAL :: sumh, hsum, dfac, sk
    REAL, DIMENSION(nz) :: vgrid_in
    CHARACTER(LEN=max_nc_filename_len) :: file
    CHARACTER(LEN=26) :: routine = "MCKPP_INITIALIZE_GEOGRAPHY"
    CHARACTER(LEN=max_message_len) :: message

    CALL mckpp_print(routine, "") 

    ! Read vertical grid fields
    IF (kpp_const_fields%l_vgrid_file) THEN

      IF (l_root) THEN 

        file = kpp_const_fields%vgrid_file
        WRITE(message,*) "Reading vertical grid from file ", TRIM(file)
        CALL mckpp_print(routine, message)

        CALL mckpp_netcdf_open(routine, file, ncid)     
        CALL mckpp_netcdf_get_var( routine, file, ncid, "d", &
                                   kpp_const_fields%dm(1:nz) )
        CALL mckpp_netcdf_get_var( routine, file, ncid, "h", & 
                                   kpp_const_fields%hm(1:nz) )
        CALL mckpp_netcdf_get_var( routine, file, ncid, "z", & 
                                   kpp_const_fields%zm(1:nz) )
        CALL mckpp_netcdf_close(routine, file, ncid)

      END IF

      CALL mckpp_broadcast_field(kpp_const_fields%dm(1:nz), nz, root) 
      CALL mckpp_broadcast_field(kpp_const_fields%hm(1:nz), nz, root) 
      CALL mckpp_broadcast_field(kpp_const_fields%zm(1:nz), nz, root) 

      kpp_const_fields%dmax = -1. * ( kpp_const_fields%zm(nz) - & 
                              kpp_const_fields%hm(nz) )

    ELSE

      IF (kpp_const_fields%l_stretchgrid) THEN
        sumh = 0.0
        dfac = 1.0 - EXP(-kpp_const_fields%dscale)
        DO i = 1, nz
          sk = - (float(i)-0.5)/float(nz)
          kpp_const_fields%hm(i) = kpp_const_fields%dmax * dfac / &
              float(nz) / kpp_const_fields%dscale / ( 1.0 + sk*dfac )
          sumh = sumh + kpp_const_fields%hm(i)
        ENDDO
      ENDIF

      ! layer thickness h, layer grids zgrid, interface depths d
      hsum = 0.0
      DO i = 1, nz
        IF (kpp_const_fields%l_stretchgrid) THEN 
          kpp_const_fields%hm(i) = kpp_const_fields%hm(i) * & 
              kpp_const_fields%dmax / sumh 
        ELSE
          kpp_const_fields%hm(i) = kpp_const_fields%dmax / real(nz) 
        END IF
        kpp_const_fields%zm(i) =  0.0 - (hsum + 0.5 * kpp_const_fields%hm(i) )
        hsum = hsum + kpp_const_fields%hm(i)
        kpp_const_fields%dm(i) = hsum
      END DO

    END IF

    kpp_const_fields%dm(0) = 0.0
    kpp_const_fields%hm(nzp1) = 1.e-10 
    kpp_const_fields%zm(nzp1) = -kpp_const_fields%dmax

    ! Calculate Coriolis parameter
    ! Enforce minimum value of Coriolis parameter equal to 2.5 degrees latitude
    DO ipt = 1, npts_local
      IF (ABS(kpp_3d_fields%dlat(ipt)) .LT. 2.5) THEN    
        kpp_3d_fields%f(ipt) = 2. * (kpp_const_fields%twopi / 86164.) * &
            SIN(2.5 * kpp_const_fields%twopi / 360.) * & 
            SIGN(1., kpp_3d_fields%dlat(ipt))
      ELSE
        kpp_3d_fields%f(ipt) = 2. * (kpp_const_fields%twopi / 86164.) * &
            SIN(kpp_3d_fields%dlat(ipt) * & 
            kpp_const_fields%twopi / 360.)
      END IF
    END DO

    WRITE(nupe,*) "kpp_const_fields%dm = ", kpp_const_fields%dm
    WRITE(nupe,*) "kpp_const_fields%hm = ", kpp_const_fields%hm
    WRITE(nupe,*) "kpp_const_fields%zm = ", kpp_const_fields%zm

  END SUBROUTINE mckpp_initialize_geography

END MODULE mckpp_initialize_geography_mod
