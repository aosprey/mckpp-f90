MODULE mckpp_initialize_ocean_profiles_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, & 
        mckpp_netcdf_close, mckpp_netcdf_determine_boundaries, & 
        mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, & 
        max_message_len, nupe
  USE mckpp_mpi_control, ONLY: l_root, root, npts_local, mckpp_scatter_field
  USE mckpp_parameters, ONLY: nx, ny, npts, nzp1

  IMPLICIT NONE

  PRIVATE :: mckpp_read_and_vinterp, mckpp_profiles_vinterp

CONTAINS


  SUBROUTINE mckpp_initialize_ocean_profiles()

    INTEGER :: k, ncid, offset_lat, offset_lon
    INTEGER, DIMENSION(3) :: start
    CHARACTER(LEN=max_nc_filename_len) :: file
    CHARACTER(LEN=31) :: routine = "MCKPP_INITIALIZE_OCEAN_PROFILES"
    CHARACTER(LEN=max_message_len) :: message

    CALL mckpp_print(routine, "")
 
    IF ( .NOT. kpp_const_fields%l_initdata ) & 
      CALL mckpp_abort(routine, "No code for l_initdata=.FALSE.")
    IF ( .NOT. kpp_const_fields%l_interpinit ) &
      CALL mckpp_abort(routine, "You have to interpolate initial profiles")

    IF (l_root) THEN 

      file = kpp_const_fields%initdata_file
      WRITE(message,*) "Reading ", TRIM(file)
      CALL mckpp_print(routine, message)
      CALL mckpp_netcdf_open(routine, file, ncid)

      CALL mckpp_netcdf_determine_boundaries( & 
        routine, file, ncid, &
        kpp_const_fields%alon, kpp_const_fields%alat, offset_lon, offset_lat)
      start = (/ offset_lon, offset_lat, 1 /)

    END IF 
    
    ! Velocity 
    CALL mckpp_read_and_vinterp( ncid, file, "zvel", "u", start, & 
                                 kpp_3d_fields%u(:,:,1) )
    CALL mckpp_print(routine, "Initialized zonal velocity")

    CALL mckpp_read_and_vinterp( ncid, file, "zvel", "v", start, & 
                                 kpp_3d_fields%u(:,:,2) )

    ! Save initial currents in case they are needed to reinitalise
    ! dodgy profiles (see resetting routines in mc-kpp_physics_overrides)
    kpp_3d_fields%u_init = kpp_3d_fields%u
    CALL mckpp_print(routine, "Initialized meridional velocity") 

    ! Temperature
    CALL mckpp_read_and_vinterp( ncid, file, "ztemp", "temp", start, & 
                                 kpp_3d_fields%x(:,:,1) )

    ! KPP requires temperatures in celsius. If initial conditions
    ! are in Kelvin, subtract 273.15
    IF ( ANY( kpp_3d_fields%x(:,:,1) .GT. 200 .AND. &
              kpp_3d_fields%x(:,:,1) .LT. 400) ) & 
      kpp_3d_fields%x(:,:,1) = kpp_3d_fields%x(:,:,1) - kpp_const_fields%tk0
    CALL mckpp_print(routine, "Initialized temperature") 

    ! Salinity 
    CALL mckpp_read_and_vinterp( ncid, file, "zsal", "sal", start, & 
                                 kpp_3d_fields%x(:,:,2) )
    CALL mckpp_print(routine, "Initialized salinity") 

    ! Calculate and remove reference salinity
    kpp_3d_fields%sref(:) = ( kpp_3d_fields%x(:,1,2) + & 
                              kpp_3d_fields%x(:,nzp1,2) ) / 2.
    kpp_3d_fields%ssref(:) = kpp_3d_fields%sref(:)
    DO k = 1, nzp1
      kpp_3d_fields%x(:,k,2) = kpp_3d_fields%x(:,k,2) - kpp_3d_fields%sref(:)
    ENDDO

    IF (kpp_const_fields%l_ssref) THEN
      kpp_3d_fields%ssurf(:) = kpp_3d_fields%ssref(:)
    ELSE
      kpp_3d_fields%ssurf(:) = kpp_3d_fields%x(:,1,2) + kpp_3d_fields%sref(:)
    ENDIF

    ! Initial surface temp
    kpp_3d_fields%tref(:) = kpp_3d_fields%x(:,1,1)

    WRITE(nupe,*) "kpp_3d_fields%u(:,2,1) = ", kpp_3d_fields%u(:,2,1)
    WRITE(nupe,*) "kpp_3d_fields%u(:,2,2) = ", kpp_3d_fields%u(:,2,2)
    WRITE(nupe,*) "kpp_3d_fields%x(:,2,1) = ", kpp_3d_fields%x(:,2,1)
    WRITE(nupe,*) "kpp_3d_fields%x(:,2,2) = ", kpp_3d_fields%x(:,2,2)

  END SUBROUTINE mckpp_initialize_ocean_profiles


  ! Read in variable, vertically interpolate and scatter. 
  ! MPI prototype - not the most efficient. 
  ! We may be reading in the vertical dimensions more 
  ! than once. Also should scatter and then interpolate locally. 
  SUBROUTINE mckpp_read_and_vinterp(ncid, file, z_var_name, var_name, start, & 
                                    var_out_local) 

    INTEGER, INTENT(IN) :: ncid
    CHARACTER(LEN=*), INTENT(IN) :: file, z_var_name, var_name 
    INTEGER, DIMENSION(3), INTENT(IN) :: start 
    REAL, DIMENSION(npts_local, nzp1), INTENT(OUT) :: var_out_local
    INTEGER :: nz_in 
    REAL, ALLOCATABLE, DIMENSION(:) :: z_in
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: var_in
    REAL, DIMENSION(npts, nzp1) :: var_out_global
    CHARACTER(LEN=22) :: routine = "MCKPP_READ_AND_VINTERP"

    IF (l_root) THEN 

      ! Get size of z dimension from file 
      CALL mckpp_netcdf_get_coord(routine, file, ncid, z_var_name, nz_in)

      ALLOCATE( z_in(nz_in) )
      ALLOCATE( var_in(nx, ny, nz_in) )

      ! Read in z dimension 
      CALL mckpp_netcdf_get_var(routine, file, ncid, z_var_name, z_in)  

      ! Read variable into temp global buffer
      CALL mckpp_netcdf_get_var(routine, file, ncid, var_name, var_in, start)
  
      ! Vertical interpolation 
      CALL mckpp_profiles_vinterp( var_in, z_in, nz_in, var_out_global )

    END IF 

    ! Scatter to local domains 
    CALL mckpp_scatter_field(var_out_global, var_out_local, root)  

  END SUBROUTINE mckpp_read_and_vinterp


  SUBROUTINE mckpp_profiles_vinterp(var_in, var_z, nz_in, var_out)

    INTEGER, INTENT(in) :: nz_in
    REAL, DIMENSION(nx, ny, nz_in), INTENT(IN) :: var_in
    REAL, DIMENSION(npts, nzp1), INTENT(OUT) :: var_out
    REAL, DIMENSION(nz_in), INTENT(IN) :: var_z
    REAL :: deltaz, deltavar
    INTEGER :: i, j, k, ipt, kin

    DO j = 1, ny
      DO i = 1, nx
        ipt = (j-1)*nx+i
        kin = 1 

        DO k = 1, nzp1
          IF ( kpp_const_fields%zm(k) .GT. var_z(1) ) THEN
            var_out(ipt,k) = var_in(i,j,1)

          ELSEIF ( kpp_const_fields%zm(k) .LT. var_z(nz_in) ) THEN
            var_out(ipt,k) = var_in(i,j,nz_in)

          ELSE
            DO WHILE ( var_z(kin+1) .GT. kpp_const_fields%zm(k) )
              kin=kin+1
            ENDDO
            deltaz = var_z(kin) - var_z(kin+1)
            deltavar = var_in(i,j,kin) - var_in(i,j,kin+1)
            var_out(ipt,k) = var_in(i,j,kin) + deltavar * & 
                             (kpp_const_fields%zm(k) - var_z(kin) ) / deltaz
          ENDIF
        ENDDO

      ENDDO
    ENDDO

  END SUBROUTINE mckpp_profiles_vinterp

END MODULE mckpp_initialize_ocean_profiles_mod
