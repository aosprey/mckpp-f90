#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif

MODULE mckpp_initialize_ocean

#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_global_fields,kpp_3d_fields,kpp_const_fields,kpp_1d_type
  USE shr_kind_mod, only: r8=>shr_kind_r8, r4=>shr_kind_r4
  USE pmgrid, only: masterproc
  USE ppgrid, only: pcols,begchunk,endchunk,pcols
  USE phys_grid, only: get_ncols_p, scatter_field_to_chunk
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields, kpp_1d_type
#endif
  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_parameters
  USE mckpp_types_transfer, ONLY: mckpp_fields_3dto1d, mckpp_fields_1dto3d
  USE mckpp_physics_verticalmixing_mod, ONLY: mckpp_physics_verticalmixing
  
  IMPLICIT NONE

CONTAINS

SUBROUTINE mckpp_initialize_ocean_profiles()
  
#ifdef MCKPP_CAM3
  INTEGER :: icol,ncol
  REAL(r8) :: temp_init(PLON,PLAT,NZP1),init_chunk(PCOLS,begchunk:endchunk,NZP1)
#endif

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

#ifdef MCKPP_CAM3
  my_nx = nx_globe 
  my_ny = ny_globe 
#else
  my_nx = nx
  my_ny = ny 
#endif  

  IF (.NOT. kpp_const_fields%L_INITDATA) THEN
    CALL mckpp_print_error(routine, "No code for L_INITDATA=.FALSE.")
    CALL MCKPP_ABORT()
  ENDIF   
  IF (.NOT. kpp_const_fields%L_INTERPINIT) THEN
    CALL mckpp_print_error(routine, "You have to interpolate initial profiles")
    CALL mckpp_abort()
  ENDIF
    
#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif
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
    
#ifdef MCKPP_CAM3
    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,temp_init(:,:,:))
#else
    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,kpp_3d_fields%U(:,:,1))
#endif
#ifdef MCKPP_CAM3
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_init,init_chunk(1,begchunk,1))
  DO ichnk=begchunk,endchunk
    ncol=get_ncols_p(ichnk)
    kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,1)=init_chunk(1:ncol,ichnk,1:NZP1)
    kpp_3d_fields(ichnk)%U_init(1:ncol,1:NZP1,1)=init_chunk(1:ncol,ichnk,1:NZP1)        
  ENDDO
#endif
  CALL mckpp_print(routine, "Initialized zonal velocity")

#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif
    CALL mckpp_netcdf_get_var(routine, file, ncid, "v", var_in, start)
#ifdef MCKPP_CAM3
    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,temp_init(:,:,:))
#else
    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,kpp_3d_fields%U(:,:,2))
#endif
#ifdef MCKPP_CAM3
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_init,init_chunk(1,begchunk,1))
  DO ichnk=begchunk,endchunk
    ncol=get_ncols_p(ichnk)
    kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,2)=init_chunk(1:ncol,ichnk,1:NZP1)
    kpp_3d_fields(ichnk)%U_init(1:ncol,1:NZP1,2)=init_chunk(1:ncol,ichnk,1:NZP1)        
  ENDDO     
#else
  ! Save initial currents in case they are needed to reinitalise
  ! dodgy profiles (see resetting routines in mc-kpp_physics_overrides)
  kpp_3d_fields%U_init = kpp_3d_fields%U
#endif
  CALL mckpp_print(routine, "Initialized meridional velocity") 

#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif
    DEALLOCATE(z_in)
    DEALLOCATE(var_in)     
    CALL mckpp_netcdf_get_coord(routine, file, ncid, "ztemp", nz_in)
    ALLOCATE(z_in(nz_in))
    ALLOCATE(var_in(my_nx, my_ny, nz_in))
    CALL mckpp_netcdf_get_var(routine, file, ncid, "ztemp", z_in)
    CALL mckpp_netcdf_get_var(routine, file, ncid, "temp", var_in, start)
    
#ifdef MCKPP_CAM3
    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,temp_init(:,:,:))
#else
    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,kpp_3d_fields%X(:,:,1))
#endif
 
    ! KPP requires temperatures in CELSIUS.  If initial conditions
    ! are in Kelvin, subtract 273.15
    offset_sst = 0.
    DO ix=1,nx
      DO iy=1,ny
        IF (var_in(ix,iy,1) .gt. 200 .and. var_in(ix,iy,1) .lt. 400) offset_sst = kpp_const_fields%TK0
      END DO
    END DO
#ifdef MCKPP_CAM3
    temp_init=temp_init-offset_sst
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_init,init_chunk(1,begchunk,1))
  DO ichnk=begchunk,endchunk
    ncol=get_ncols_p(ichnk)
    kpp_3d_fields(ichnk)%X(1:ncol,1:NZP1,1)=init_chunk(1:ncol,ichnk,1:NZP1)
  ENDDO
#else
  kpp_3d_fields%X(:,:,1) = kpp_3d_fields%X(:,:,1) - offset_sst
#endif
  CALL mckpp_print(routine, "Initialized temperature") 

#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif
    DEALLOCATE(z_in)
    DEALLOCATE(var_in)     
    CALL mckpp_netcdf_get_coord(routine, file, ncid, "zsal", nz_in)
    ALLOCATE(z_in(nz_in))
    ALLOCATE(var_in(my_nx, my_ny, nz_in))
    CALL mckpp_netcdf_get_var(routine, file, ncid, "zsal", z_in)
    CALL mckpp_netcdf_get_var(routine, file, ncid, "sal", var_in, start)
    CALL mckpp_netcdf_close(routine, file, ncid)
  
#ifdef MCKPP_CAM3
    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,temp_init(:,:,:))
#else
    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,kpp_3d_fields%X(:,:,2))
#endif
#ifdef MCKPP_CAM3
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_init,init_chunk(1,begchunk,1))
  DO ichnk=begchunk,endchunk
    ncol=get_ncols_p(ichnk)
    kpp_3d_fields(ichnk)%X(1:ncol,1:NZP1,2)=init_chunk(1:ncol,ichnk,1:NZP1)
  ENDDO
#endif
  CALL mckpp_print(routine, "Initialized salinity") 
 
  ! Calculate and remove reference salinity
  
#ifdef MCKPP_CAM3
  DO ichnk=begchunk,endchunk
    kpp_3d_fields(ichnk)%Sref(:)=(kpp_3d_fields(ichnk)%X(:,1,2)+kpp_3d_fields(ichnk)%X(:,nzp1,2))/2.
    kpp_3d_fields(ichnk)%Ssref(:)=kpp_3d_fields(ichnk)%Sref(:)
    do i=1,nzp1
      kpp_3d_fields(ichnk)%X(:,i,2)=kpp_3d_fields(ichnk)%X(:,i,2)-kpp_3d_fields(ichnk)%Sref(:)
    enddo
    ! Initial surface temp
    kpp_3d_fields(ichnk)%Tref(:) = kpp_3d_fields(ichnk)%X(:,1,1)
    IF (kpp_const_fields%L_SSref) THEN
      kpp_3d_fields(ichnk)%Ssurf(:)=kpp_3d_fields(ichnk)%SSref(:)
    ELSE
      kpp_3d_fields(ichnk)%Ssurf(:)=kpp_3d_fields(ichnk)%X(:,1,2)+kpp_3d_fields(ichnk)%Sref(:)
    ENDIF
  ENDDO
#else
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
#endif
  
END SUBROUTINE mckpp_initialize_ocean_profiles


SUBROUTINE mckpp_initialize_ocean_profiles_vinterp(var_in,var_z,nz_in,model_z,var_out)
 
  INTEGER, intent(in) :: nz_in

#ifdef MCKPP_CAM3
  REAL(r8),dimension(PLON,PLAT,NZ_IN), intent(in) :: var_in
  REAL(r8),dimension(PLON,PLAT,NZP1), intent(out) :: var_out
#else
  REAL,dimension(NX,NY,NZ_IN),intent(in) :: var_in
  REAL,dimension(NPTS,NZP1),intent(out) :: var_out
#endif
 
  INTEGER :: my_nx, my_ny
  REAL,dimension(NZ_IN), intent(in) :: var_z
  REAL,dimension(NZP1), intent(in) :: model_z
  REAL :: deltaz,deltavar
  INTEGER :: ix,iy,iz,ipt,kin

#ifdef MCKPP_CAM3
  my_nx=PLON
  my_ny=PLAT
#else
  my_nx=NX
  my_ny=NY
#endif

  DO iy=1,my_ny
     DO ix=1,my_nx
        ipt=(iy-1)*my_nx+ix
        kin=1
        DO iz=1,NZP1
           IF (model_z(iz) .GT. var_z(1)) THEN
#ifdef MCKPP_CAM3
              var_out(ix,iy,iz)=var_in(ix,iy,1)
#else              
              var_out(ipt,iz)=var_in(ix,iy,1)
#endif
           ELSEIF (model_z(iz) .LT. var_z(nz_in)) THEN
#ifdef MCKPP_CAM3
              var_out(ix,iy,iz)=var_in(ix,iy,nz_in)
#else
              var_out(ipt,iz)=var_in(ix,iy,nz_in)
#endif
           ELSE
              DO WHILE (var_z(kin+1) .GT. model_z(iz))
                 kin=kin+1
              ENDDO
              deltaz=var_z(kin)-var_z(kin+1)
              deltavar=var_in(ix,iy,kin)-var_in(ix,iy,kin+1)
#ifdef MCKPP_CAM3
              var_out(ix,iy,iz)=var_in(ix,iy,kin)+deltavar*(model_z(iz)-var_z(kin))/deltaz
#else
              var_out(ipt,iz)=var_in(ix,iy,kin)+deltavar*(model_z(iz)-var_z(kin))/deltaz
#endif              
           ENDIF
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE mckpp_initialize_ocean_profiles_vinterp


! Initialize ocean model:
! Set coefficients for tridiagonal matrix solver.
! Compute hmix and diffusivity profiles for initial profile.
! Prepare for first time step.
SUBROUTINE MCKPP_INITIALIZE_OCEAN_MODEL()
  
#ifdef MCKPP_CAM3
  INTEGER :: icol,ncol,ichnk
#endif
  
  ! Local
  TYPE(kpp_1d_type) :: kpp_1d_fields
  real dzb(NZ)              ! diff. between grid-levels below z(j)
  integer k,kmix0,n,l,ipt
  real hmix0,deltaz

  ! Compute factors for coefficients of tridiagonal matrix elements.
  ! tri(0     ,1,.........) : dt/h(1) factor for rhs flux
  ! tri(k=1:NZ,0,.........) : dt/h(k)/ {dzb(k-1)=z(k-1)-z(k)=dzabove}
  ! tri(k=1:NZ,1,.........) : dt/h(k)/ {dzb(k  )=z(k)-z(k+1)=dzbelow}

  DO k=1,NZ
     dzb(k) = kpp_const_fields%zm(k) - kpp_const_fields%zm(k+1)
  ENDDO
  
  kpp_const_fields%tri(0,1,1) = kpp_const_fields%dto/kpp_const_fields%hm(1)
  kpp_const_fields%tri(1,1,1) = kpp_const_fields%dto/kpp_const_fields%hm(1)/dzb(1)
  DO k=2,NZ
     kpp_const_fields%tri(k,1,1) = kpp_const_fields%dto/kpp_const_fields%hm(k)/dzb(k)
     kpp_const_fields%tri(k,0,1) = kpp_const_fields%dto/kpp_const_fields%hm(k)/dzb(k-1)
  ENDDO
  
  IF ( .NOT. kpp_const_fields%L_RESTART) THEN    
     ! Determine hmix for initial profile:

#ifdef MCKPP_CAM3
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        DO icol=1,ncol
           IF (kpp_3d_fields(ichnk)%L_OCEAN(icol) .and. &
                kpp_3d_fields(ichnk)%cplwght(icol) .gt. 0.0) THEN
              CALL mckpp_fields_3dto1d(kpp_3d_fields(ichnk),icol,kpp_1d_fields)
#else
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(none) &
!$OMP SHARED(kpp_3d_fields, kpp_const_fields) &
!$OMP SHARED(nz, nzp1, nx, ny, npts, nvel, nsclr, nvp1, nsp1, itermax) &
!$OMP SHARED(hmixtolfrac, nztmax, nzp1tmax, nsflxs, njdt, maxmodeadv) &
!$OMP PRIVATE(ipt, k, l, deltaz, kpp_1d_fields, hmix0, kmix0)
!$OMP DO SCHEDULE(dynamic)
#endif
     DO ipt=1,npts
#ifdef MCKPP_COUPLE
        IF (kpp_3d_fields%L_OCEAN(ipt) .and. kpp_3d_fields%cplwght(ipt) .gt. 0.0) THEN
#else       
        IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
#endif
           CALL mckpp_fields_3dto1d(kpp_3d_fields,ipt,kpp_1d_fields)
#endif
           kpp_1d_fields%L_INITFLAG=.TRUE.
           CALL MCKPP_PHYSICS_VERTICALMIXING(kpp_1d_fields,kpp_const_fields,hmix0,kmix0)
           kpp_1d_fields%L_INITFLAG=.FALSE.
           kpp_1d_fields%hmix = hmix0
           kpp_1d_fields%kmix = kmix0
           kpp_1d_fields%Tref = kpp_1d_fields%X(1,1)
           ! Evaluate initial fluxes (to write to output data file)
           DO k=1,NZ
              deltaz = 0.5*(kpp_const_fields%hm(k)+kpp_const_fields%hm(k+1))
              DO n=1,NSCLR
                 kpp_1d_fields%wX(k,n)=-kpp_1d_fields%difs(k)*&
                      ((kpp_1d_fields%X(k,n)-kpp_1d_fields%X(k+1,n))/deltaz-&
                      kpp_1d_fields%ghat(k)*kpp_1d_fields%wX(0,n))
              ENDDO
              IF(kpp_const_fields%LDD) kpp_1d_fields%wX(k,1)=-kpp_1d_fields%dift(k)*&
                   ((kpp_1d_fields%X(k,1)-kpp_1d_fields%X(k+1,1))/deltaz-kpp_1d_fields%ghat(k)*kpp_1d_fields%wX(0,1))
              kpp_1d_fields%wX(k,nsp1)= kpp_const_fields%grav * (kpp_1d_fields%talpha(k)*kpp_1d_fields%wX(k,1) - &
                   kpp_1d_fields%sbeta(k) * kpp_1d_fields%wX(k,2))
              DO  n=1,NVEL
                 kpp_1d_fields%wU(k,n)= -kpp_1d_fields%difm(k)*&
                      (kpp_1d_fields%U(k,n)-kpp_1d_fields%U(k+1,n))/deltaz
              ENDDO
           ENDDO
           
           ! Prepare for first time step
           
           ! indices for extrapolation
           kpp_1d_fields%old = 0
           kpp_1d_fields%new = 1               
           ! initialize array for extrapolating hmixd,Us,Xs
           kpp_1d_fields%hmixd(0) = kpp_1d_fields%hmix
           kpp_1d_fields%hmixd(1) = kpp_1d_fields%hmix
           DO k=1,NZP1
              DO l=1,NVEL
                 kpp_1d_fields%Us(k,l,0)=kpp_1d_fields%U(k,l)
                 kpp_1d_fields%Us(k,l,1)=kpp_1d_fields%U(k,l)
              ENDDO
              DO l=1,NSCLR
                 kpp_1d_fields%Xs(k,l,0)=kpp_1d_fields%X(k,l)
                 kpp_1d_fields%Xs(k,l,1)=kpp_1d_fields%X(k,l)
              ENDDO
           ENDDO
#ifdef MCKPP_CAM3
           CALL mckpp_fields_1dto3d(kpp_1d_fields,icol,kpp_3d_fields(ichnk))
        ENDIF
     ENDDO
  ENDDO
#else
           CALL mckpp_fields_1dto3d(kpp_1d_fields,ipt,kpp_3d_fields)
        ENDIF
     ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif      
#endif

  ENDIF

END SUBROUTINE mckpp_initialize_ocean_model

END MODULE mckpp_initialize_ocean
