#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE mckpp_initialize_ocean_profiles
  USE mckpp_types, only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE shr_kind_mod, only: r8=>shr_kind_r8
  USE pmgrid, only: masterproc
  USE ppgrid, only: pcols,begchunk,endchunk
  USE phys_grid, only: get_ncols_p, scatter_field_to_chunk
#else
SUBROUTINE mckpp_initialize_ocean_profiles(kpp_3d_fields,kpp_const_fields)
#endif

  IMPLICIT NONE
  INTEGER, parameter :: nuout=6,nuerr=0
#include <netcdf.inc>  

#ifdef MCKPP_CAM3
#include <parameter.inc>
  INTEGER :: icol,ncol
  INTEGER, parameter :: my_nx=NX_GLOBE,my_ny=NY_GLOBE
  REAL(r8) :: temp_init(PLON,PLAT,NZP1),init_chunk(PCOLS,begchunk:endchunk,NZP1)
#else
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>  
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER, parameter :: my_nx=NX,my_ny=NY
#endif

! local
  INTEGER n,i,ipt,ichnk
  INTEGER status,ncid
  REAL*4, allocatable :: var_in(:,:,:),z_in(:),x_in(:),y_in(:)
  
  INTEGER varid,dimid
  INTEGER nz_in,nx_in,ny_in
  
  INTEGER ix,iy,count(3),start(3)
  INTEGER kin,k
  REAL deltaz,deltavar,offset_sst
  
  IF ( kpp_const_fields%L_INITDATA ) THEN
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
     allocate(var_in(my_nx,my_ny,200))
     allocate(z_in(200))
     allocate(x_in(NX_GLOBE))
     allocate(y_in(NY_GLOBE))
     
     status=NF_OPEN(kpp_const_fields%initdata_file,0,ncid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)     
     status=NF_INQ_DIMID(ncid,'longitude',dimid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_DIMLEN (ncid,dimid,nx_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'longitude',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,x_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     ix=1
#ifdef MCKPP_CAM3
     DO WHILE (abs(x_in(ix)-kpp_global_fields%longitude(1)) .GT. 1e-3)
#else
     DO WHILE (abs(x_in(ix)-kpp_3d_fields%dlon(1)) .GT. 1.e-3)
#endif
        ix=ix+1
        IF (ix .GE. nx_in) THEN
           WRITE(nuerr,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: Error reading initial conditions'
#ifdef MCKPP_CAM3
           WRITE(nuerr,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: Cannot find longitude ',&
                kpp_global_fields%longitude(1),' in range ',x_in(1),x_in(nx_in)
#else
           WRITE(nuerr,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: Cannot find longitude ',&
                kpp_3d_fields%dlon(1),' in range ',x_in(1),x_in(nx_in)
#endif
           CALL MCKPP_ABORT
        ENDIF
     ENDDO
     start(1)=ix
     count(1)=my_nx

     status=NF_INQ_DIMID(ncid,'latitude',dimid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_DIMLEN (ncid,dimid,ny_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'latitude',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,y_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     iy=1
#ifdef MCKPP_CAM3
     DO WHILE (abs(y_in(iy)-kpp_global_fields%latitude(1)) .GT. 1e-3)
#else
     DO WHILE (abs(y_in(iy)-kpp_3d_fields%dlat(1)) .GT. 1.e-3)
#endif
        iy=iy+1
        IF (iy .GE. ny_in) THEN
           write(nuerr,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: Error reading initial conditions'
#ifdef MCKPP_CAM3
           WRITE(nuerr,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: Cannot find latitude ',&
                kpp_global_fields%latitude(1),' in range ',y_in(1),y_in(ny_in)
#else
           write(nuerr,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: Cannot find latitude ',&
                kpp_3d_fields%dlat(1),' in range ',y_in(1),y_in(ny_in)
#endif
           CALL MCKPP_ABORT
        ENDIF
     ENDDO
     start(2)=iy
     count(2)=my_ny
     
     status=NF_INQ_DIMID(ncid,'zvel',dimid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_DIMLEN (ncid,dimid,nz_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     start(3)=1
     count(3)=nz_in
     
     status=NF_INQ_VARID(ncid,'zvel',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,z_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'u',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     
     IF (kpp_const_fields%L_INTERPINIT) THEN
#ifdef MCKPP_CAM3
        CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,temp_init(:,:,:))
#else
        CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,kpp_3d_fields%U(:,:,1))
#endif
     ELSE
        write(nuerr,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: You have to interpolate'
     ENDIF
#ifdef MCKPP_CAM3
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_init,init_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,1)=init_chunk(1:ncol,ichnk,1:NZP1)
        kpp_3d_fields(ichnk)%U_init(1:ncol,1:NZP1,1)=init_chunk(1:ncol,ichnk,1:NZP1)        
     ENDDO
     IF (masterproc) THEN
        WRITE(6,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: Initialized zonal velocity'
#endif
     status=NF_INQ_VARID(ncid,'v',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)

     IF (kpp_const_fields%L_INTERPINIT) THEN
#ifdef MCKPP_CAM3
        CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,temp_init(:,:,:))
#else
        CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,kpp_3d_fields%U(:,:,2))
#endif
     ELSE
        write(nuerr,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: You have to interpolate'
     ENDIF
#ifdef MCKPP_CAM3
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_init,init_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,2)=init_chunk(1:ncol,ichnk,1:NZP1)
        kpp_3d_fields(ichnk)%U_init(1:ncol,1:NZP1,2)=init_chunk(1:ncol,ichnk,1:NZP1)        
     ENDDO     
     IF (masterproc) THEN 
        WRITE(nuout,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: Initialized meridional velocity'
#else
     ! Save initial currents in case they are needed to reinitalise
     ! dodgy profiles (see resetting routines in mc-kpp_physics_overrides)
     kpp_3d_fields%U_init=kpp_3d_fields%U
#endif     
     status=NF_INQ_DIMID(ncid,'ztemp',dimid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_DIMLEN (ncid,dimid,nz_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     start(3)=1
     count(3)=nz_in
     
     status=NF_INQ_VARID(ncid,'ztemp',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,z_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'temp',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     
     IF (kpp_const_fields%L_INTERPINIT) THEN
#ifdef MCKPP_CAM3
        CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,temp_init(:,:,:))
#else
        CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,kpp_3d_fields%X(:,:,1))
#endif
     ELSE
        WRITE(nuerr,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: You have to interpolate'
     ENDIF

     ! KPP requires temperatures in CELSIUS.  If initial conditions
     ! are in Kelvin, subtract 273.15
     offset_sst = 0.
     DO ix=1,nx
        DO iy=1,ny
           IF (var_in(ix,iy,1) .gt. 200 .and. var_in(ix,iy,1) .lt. 400) offset_sst = kpp_const_fields%TK0
        END DO
     END DO
     temp_init=temp_init-offset_sst
#ifdef MCKPP_CAM3
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_init,init_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%X(1:ncol,1:NZP1,1)=init_chunk(1:ncol,ichnk,1:NZP1)
     ENDDO
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: Initialized temperature'
#else
     kpp_3d_fields%X(:,:,1) = kpp_3d_fields%X(:,:,1) - offset_sst     
#endif
          
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
     status=NF_INQ_DIMID(ncid,'zsal',dimid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_DIMLEN (ncid,dimid,nz_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     start(3)=1
     count(3)=nz_in
     
     status=NF_INQ_VARID(ncid,'zsal',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,z_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'sal',varid)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VARA_REAL(ncid,varid,start,count,var_in)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     
     IF (kpp_const_fields%L_INTERPINIT) THEN
#ifdef MCKPP_CAM3
        CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,temp_init(:,:,:))
#else
        CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,kpp_3d_fields%X(:,:,2))
#endif
     ELSE
        WRITE(nuerr,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: You have to interpolate'
     ENDIF
#ifdef MCKPP_CAM3
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_init,init_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%X(1:ncol,1:NZP1,2)=init_chunk(1:ncol,ichnk,1:NZP1)
     ENDDO
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_OCEAN_PROFILES: Initialized salinity' 
#endif
  ELSE
     WRITE(nuerr,*) "MCKPP_INITIALIZE_OCEAN_PROFILES: No code for L_INITDATA=.FALSE."
     CALL MCKPP_ABORT
  ENDIF

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
  
  RETURN
END SUBROUTINE mckpp_initialize_ocean_profiles

#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif
SUBROUTINE mckpp_initialize_ocean_profiles_vinterp(var_in,var_z,nz_in,model_z,var_out)
#ifdef MCKPP_CAM3
  USE shr_kind_mod, only : r8=>shr_kind_r8, r4=>shr_kind_r4
#endif  
  IMPLICIT NONE
  
  INTEGER,parameter :: nuout=6,nuerr=0
#include <parameter.inc>
  INTEGER, intent(in) :: nz_in

#ifdef MCKPP_CAM3
  REAL(r4),dimension(PLON,PLAT,NZ_IN), intent(in) :: var_in
  REAL(r8),dimension(PLON,PLAT,NZP1), intent(out) :: var_out
  INTEGER, parameter :: my_nx=PLON,my_ny=PLAT
#else
  REAL*4,dimension(NX,NY,NZ_IN),intent(in) :: var_in
  REAL,dimension(NPTS,NZP1),intent(out) :: var_out
  INTEGER, parameter :: my_nx=NX,my_ny=NY
#endif

  REAL*4,dimension(NZ_IN), intent(in) :: var_z
  REAL,dimension(NZP1), intent(in) :: model_z
  REAL :: deltaz,deltavar
  INTEGER :: ix,iy,iz,ipt,kin

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
  RETURN
END SUBROUTINE mckpp_initialize_ocean_profiles_vinterp

#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE MCKPP_INITIALIZE_OCEAN_MODEL
  USE shr_kind_mod, only: r8=>shr_kind_r8
  USE mckpp_types, only: kpp_3d_fields,kpp_const_fields,kpp_1d_type
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid, only: get_ncols_p  
#else
SUBROUTINE mckpp_initialize_ocean_model(kpp_3d_fields,kpp_const_fields)
#endif

  ! Initialize ocean model:
  ! Set coefficients for tridiagonal matrix solver.
  ! Compute hmix and diffusivity profiles for initial profile.
  ! Prepare for first time step.
  
  IMPLICIT NONE
  INTEGER, parameter :: nuout=6,nuerr=0

#ifdef MCKPP_CAM3
#include <parameter.inc>
  INTEGER :: icol,ncol,ichnk
#else
  ! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
  ! Input
  TYPE(kpp_const_type) :: kpp_const_fields  
  ! Output
  TYPE(kpp_3d_type) :: kpp_3d_fields
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
!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(kpp_3d_fields,kpp_const_fields)
!$OMP DO SCHEDULE(dynamic)
#endif
     DO ipt=1,npts
        IF (kpp_3d_fields%L_OCEAN(ipt) .and. kpp_3d_fields%cplwght(ipt) .gt. 0.0) THEN
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
  RETURN
END SUBROUTINE mckpp_initialize_ocean_model
