MODULE mckpp_initialize_geography_mod

  USE mckpp_data_fields, ONLY: kpp_const_fields, kpp_3d_fields
  USE mckpp_netcdf_read, ONLY: mckpp_netcdf_open, mckpp_netcdf_close, mckpp_netcdf_get_var, &
      max_nc_filename_len
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: nz, nzp1, npts

  IMPLICIT NONE

CONTAINS

  SUBROUTINE mckpp_initialize_geography()

    ! Local Variables
    REAL :: sumh, hsum, dfac, sk
    REAL, DIMENSION(nz) :: vgrid_in
    INTEGER :: i, ipt, ncid
    CHARACTER(LEN=max_nc_filename_len) :: file
    CHARACTER(LEN=26) :: routine = "MCKPP_INITIALIZE_GEOGRAPHY"
    CHARACTER(LEN=max_message_len) :: message

    ! Read vertical grid fields
    IF (kpp_const_fields%L_VGRID_FILE) THEN
      file = kpp_const_fields%vgrid_file
      WRITE(message,*) "Reading vertical grid from file ", TRIM(file)
      CALL mckpp_print(routine, message)

      CALL mckpp_netcdf_open(routine, file, ncid)     
      CALL mckpp_netcdf_get_var(routine, file, ncid, "d", kpp_const_fields%dm(1:NZ))
      CALL mckpp_netcdf_get_var(routine, file, ncid, "h", kpp_const_fields%hm(1:NZ))
      CALL mckpp_netcdf_get_var(routine, file, ncid, "z", kpp_const_fields%zm(1:NZ))
      CALL mckpp_netcdf_close(routine, file, ncid)

      kpp_const_fields%DMAX=-1.*(kpp_const_fields%zm(NZ)-kpp_const_fields%hm(NZ))
    ELSE     
      IF (kpp_const_fields%L_STRETCHGRID) THEN
        sumh = 0.0
        dfac = 1.0 - exp(-kpp_const_fields%dscale)
        DO i = 1,NZ
          sk = - (float(i)-0.5)/float(NZ)
          kpp_const_fields%hm(i) = kpp_const_fields%DMAX*dfac/&
              float(NZ)/kpp_const_fields%dscale / ( 1.0 + sk*dfac )
          sumh = sumh + kpp_const_fields%hm(i)
        ENDDO
      ENDIF

      ! layer thickness h, layer grids zgrid, interface depths d
      hsum = 0.0
      DO i=1,NZ
        if(kpp_const_fields%L_STRETCHGRID) then
          kpp_const_fields%hm(i) = kpp_const_fields%hm(i) * kpp_const_fields%DMAX / sumh 
        else   
          kpp_const_fields%hm(i) = kpp_const_fields%DMAX / real(NZ) 
        endif
        kpp_const_fields%zm(i) =  0.0 - (hsum + 0.5 * kpp_const_fields%hm(i) )
        hsum = hsum + kpp_const_fields%hm(i)
        kpp_const_fields%dm(i) = hsum
      ENDDO
    ENDIF
    kpp_const_fields%dm(0) = 0.0
    kpp_const_fields%hm(nzp1) = 1.e-10 
    kpp_const_fields%zm(nzp1) = -kpp_const_fields%DMAX

    ! Calculate Coriolis parameter
    ! Enforce minimum value of Coriolis parameter equal to 2.5 degrees latitude

    DO ipt=1,npts
      if(abs(kpp_3d_fields%dlat(ipt)).lt.2.5) then    
        kpp_3d_fields%f(ipt) = 2. * (kpp_const_fields%twopi/86164.) * &
            sin(2.5*kpp_const_fields%twopi/360.)*SIGN(1.,kpp_3d_fields%dlat(ipt))
      else  
        kpp_3d_fields%f(ipt) = 2. * (kpp_const_fields%twopi/86164.) * &
            sin(kpp_3d_fields%dlat(ipt)*kpp_const_fields%twopi/360.)
      endif
    ENDDO

  END SUBROUTINE mckpp_initialize_geography

END MODULE mckpp_initialize_geography_mod
