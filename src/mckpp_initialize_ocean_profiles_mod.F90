MODULE mckpp_initialize_ocean_profiles_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, mckpp_netcdf_close, &
       mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, npts, nx_globe, ny_globe, nzp1

  IMPLICIT NONE

CONTAINS

  SUBROUTINE mckpp_initialize_ocean_profiles()

    ! local
    REAL, ALLOCATABLE, DIMENSION(:) :: z_in
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: var_in
    INTEGER :: my_nx, my_ny, i, ichnk, ix, iy
    INTEGER :: ncid, offset_lat, offset_lon, nz_in
    INTEGER, DIMENSION(3) :: start, count
    REAL :: offset_sst

    CHARACTER(LEN=max_nc_filename_len) :: file
    CHARACTER(LEN=31) :: routine = "MCKPP_INITIALIZE_OCEAN_PROFILES"
    CHARACTER(LEN=max_message_len) :: message

    my_nx = nx
    my_ny = ny 

    IF (.NOT. kpp_const_fields%L_INITDATA) THEN
       CALL mckpp_print_error(routine, "No code for L_INITDATA=.FALSE.")
       CALL MCKPP_ABORT()
    ENDIF
    IF (.NOT. kpp_const_fields%L_INTERPINIT) THEN
       CALL mckpp_print_error(routine, "You have to interpolate initial profiles")
       CALL mckpp_abort()
    ENDIF

    file = kpp_const_fields%initdata_file
    WRITE(message,*) "Reading ", TRIM(file)
    CALL mckpp_print(routine, message)

    CALL mckpp_netcdf_open(routine, file, ncid)
    CALL mckpp_netcdf_determine_boundaries(routine, file, ncid, &
         kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), offset_lon, offset_lat)     
    start(1) = offset_lon
    start(2) = offset_lat
    start(3) = 1

    CALL mckpp_netcdf_get_coord(routine, file, ncid, "zvel", nz_in)
    ALLOCATE(z_in(nz_in))
    ALLOCATE(var_in(my_nx, my_ny, nz_in))
    CALL mckpp_netcdf_get_var(routine, file, ncid, "zvel", z_in)     
    CALL mckpp_netcdf_get_var(routine, file, ncid, "u", var_in, start)

    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,kpp_3d_fields%U(:,:,1))
    CALL mckpp_print(routine, "Initialized zonal velocity")

    CALL mckpp_netcdf_get_var(routine, file, ncid, "v", var_in, start)
    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,kpp_3d_fields%U(:,:,2))
    ! Save initial currents in case they are needed to reinitalise
    ! dodgy profiles (see resetting routines in mc-kpp_physics_overrides)
    kpp_3d_fields%U_init = kpp_3d_fields%U
    CALL mckpp_print(routine, "Initialized meridional velocity") 

    DEALLOCATE(z_in)
    DEALLOCATE(var_in)     
    CALL mckpp_netcdf_get_coord(routine, file, ncid, "ztemp", nz_in)
    ALLOCATE(z_in(nz_in))
    ALLOCATE(var_in(my_nx, my_ny, nz_in))
    CALL mckpp_netcdf_get_var(routine, file, ncid, "ztemp", z_in)
    CALL mckpp_netcdf_get_var(routine, file, ncid, "temp", var_in, start)

    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,kpp_3d_fields%X(:,:,1))

    ! KPP requires temperatures in CELSIUS.  If initial conditions
    ! are in Kelvin, subtract 273.15
    offset_sst = 0.
    DO ix=1,nx
       DO iy=1,ny
          IF (var_in(ix,iy,1) .gt. 200 .and. var_in(ix,iy,1) .lt. 400) offset_sst = kpp_const_fields%TK0
       END DO
    END DO
    kpp_3d_fields%X(:,:,1) = kpp_3d_fields%X(:,:,1) - offset_sst
    CALL mckpp_print(routine, "Initialized temperature") 

    DEALLOCATE(z_in)
    DEALLOCATE(var_in)     
    CALL mckpp_netcdf_get_coord(routine, file, ncid, "zsal", nz_in)
    ALLOCATE(z_in(nz_in))
    ALLOCATE(var_in(my_nx, my_ny, nz_in))
    CALL mckpp_netcdf_get_var(routine, file, ncid, "zsal", z_in)
    CALL mckpp_netcdf_get_var(routine, file, ncid, "sal", var_in, start)
    CALL mckpp_netcdf_close(routine, file, ncid)

    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,kpp_3d_fields%X(:,:,2))
    CALL mckpp_print(routine, "Initialized salinity") 

    ! Calculate and remove reference salinity

    kpp_3d_fields%Sref(:)=(kpp_3d_fields%X(:,1,2)+kpp_3d_fields%X(:,nzp1,2))/2.
    kpp_3d_fields%Ssref(:)=kpp_3d_fields%Sref(:)
    do i=1,nzp1
       kpp_3d_fields%X(:,i,2)=kpp_3d_fields%X(:,i,2)-kpp_3d_fields%Sref(:)
    enddo
    ! Initial surface temp
    kpp_3d_fields%Tref(:) = kpp_3d_fields%X(:,1,1)
    IF (kpp_const_fields%L_SSref) THEN
       kpp_3d_fields%Ssurf(:)=kpp_3d_fields%SSref(:)
    ELSE
       kpp_3d_fields%Ssurf(:)=kpp_3d_fields%X(:,1,2)+kpp_3d_fields%Sref(:)
    ENDIF

  END SUBROUTINE mckpp_initialize_ocean_profiles


  SUBROUTINE mckpp_initialize_ocean_profiles_vinterp(var_in,var_z,nz_in,model_z,var_out)

    INTEGER, intent(in) :: nz_in
    REAL,dimension(NX,NY,NZ_IN),intent(in) :: var_in
    REAL,dimension(NPTS,NZP1),intent(out) :: var_out

    INTEGER :: my_nx, my_ny
    REAL,dimension(NZ_IN), intent(in) :: var_z
    REAL,dimension(NZP1), intent(in) :: model_z
    REAL :: deltaz,deltavar
    INTEGER :: ix,iy,iz,ipt,kin

    my_nx=NX
    my_ny=NY

    DO iy=1,my_ny
       DO ix=1,my_nx
          ipt=(iy-1)*my_nx+ix
          kin=1
          DO iz=1,NZP1
             IF (model_z(iz) .GT. var_z(1)) THEN
                var_out(ipt,iz)=var_in(ix,iy,1)
             ELSEIF (model_z(iz) .LT. var_z(nz_in)) THEN
                var_out(ipt,iz)=var_in(ix,iy,nz_in)
             ELSE
                DO WHILE (var_z(kin+1) .GT. model_z(iz))
                   kin=kin+1
                ENDDO
                deltaz=var_z(kin)-var_z(kin+1)
                deltavar=var_in(ix,iy,kin)-var_in(ix,iy,kin+1)
                var_out(ipt,iz)=var_in(ix,iy,kin)+deltavar*(model_z(iz)-var_z(kin))/deltaz
             ENDIF
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE mckpp_initialize_ocean_profiles_vinterp

END MODULE mckpp_initialize_ocean_profiles_mod
