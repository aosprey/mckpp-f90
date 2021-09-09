#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>

MODULE mckpp_initialize_ocean_profiles_mod

  USE mckpp_types, only: kpp_3d_fields,kpp_const_fields
  USE shr_kind_mod, only: r8=>shr_kind_r8
  USE pmgrid, only: masterproc
  USE ppgrid, only: pcols,begchunk,endchunk,pcols
  USE phys_grid, only: get_ncols_p, scatter_field_to_chunk
  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, npts, nx_globe, ny_globe, nzp1
  
  IMPLICIT NONE

CONTAINS

SUBROUTINE mckpp_initialize_ocean_profiles()
  
  ! local
  INTEGER :: icol,ncol
  REAL(r8) :: temp_init(PLON,PLAT,NZP1),init_chunk(PCOLS,begchunk:endchunk,NZP1)
  REAL, ALLOCATABLE, DIMENSION(:) :: z_in
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: var_in
  INTEGER :: my_nx, my_ny, i, ichnk, ix, iy
  INTEGER :: ncid, offset_lat, offset_lon, nz_in
  INTEGER, DIMENSION(3) :: start, count
  REAL :: offset_sst
  
  CHARACTER(LEN=max_nc_filename_len) :: file
  CHARACTER(LEN=31) :: routine = "MCKPP_INITIALIZE_OCEAN_PROFILES"
  CHARACTER(LEN=max_message_len) :: message

  my_nx = nx_globe 
  my_ny = ny_globe 

  IF (.NOT. kpp_const_fields%L_INITDATA) THEN
    CALL mckpp_print_error(routine, "No code for L_INITDATA=.FALSE.")
    CALL MCKPP_ABORT()
  ENDIF   
  IF (.NOT. kpp_const_fields%L_INTERPINIT) THEN
    CALL mckpp_print_error(routine, "You have to interpolate initial profiles")
    CALL mckpp_abort()
  ENDIF
    
  IF (masterproc) THEN
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
    
    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,temp_init(:,:,:))
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_init,init_chunk(1,begchunk,1))
  DO ichnk=begchunk,endchunk
    ncol=get_ncols_p(ichnk)
    kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,1)=init_chunk(1:ncol,ichnk,1:NZP1)
    kpp_3d_fields(ichnk)%U_init(1:ncol,1:NZP1,1)=init_chunk(1:ncol,ichnk,1:NZP1)        
  ENDDO
  CALL mckpp_print(routine, "Initialized zonal velocity")

  IF (masterproc) THEN
    CALL mckpp_netcdf_get_var(routine, file, ncid, "v", var_in, start)
    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,temp_init(:,:,:))
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_init,init_chunk(1,begchunk,1))
  DO ichnk=begchunk,endchunk
    ncol=get_ncols_p(ichnk)
    kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,2)=init_chunk(1:ncol,ichnk,1:NZP1)
    kpp_3d_fields(ichnk)%U_init(1:ncol,1:NZP1,2)=init_chunk(1:ncol,ichnk,1:NZP1)        
  ENDDO     
  CALL mckpp_print(routine, "Initialized meridional velocity") 

  IF (masterproc) THEN
    DEALLOCATE(z_in)
    DEALLOCATE(var_in)     
    CALL mckpp_netcdf_get_coord(routine, file, ncid, "ztemp", nz_in)
    ALLOCATE(z_in(nz_in))
    ALLOCATE(var_in(my_nx, my_ny, nz_in))
    CALL mckpp_netcdf_get_var(routine, file, ncid, "ztemp", z_in)
    CALL mckpp_netcdf_get_var(routine, file, ncid, "temp", var_in, start)
    
    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,temp_init(:,:,:))
 
    ! KPP requires temperatures in CELSIUS.  If initial conditions
    ! are in Kelvin, subtract 273.15
    offset_sst = 0.
    DO ix=1,nx
      DO iy=1,ny
        IF (var_in(ix,iy,1) .gt. 200 .and. var_in(ix,iy,1) .lt. 400) offset_sst = kpp_const_fields%TK0
      END DO
    END DO
    temp_init=temp_init-offset_sst
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_init,init_chunk(1,begchunk,1))
  DO ichnk=begchunk,endchunk
    ncol=get_ncols_p(ichnk)
    kpp_3d_fields(ichnk)%X(1:ncol,1:NZP1,1)=init_chunk(1:ncol,ichnk,1:NZP1)
  ENDDO
  CALL mckpp_print(routine, "Initialized temperature") 

  IF (masterproc) THEN
    DEALLOCATE(z_in)
    DEALLOCATE(var_in)     
    CALL mckpp_netcdf_get_coord(routine, file, ncid, "zsal", nz_in)
    ALLOCATE(z_in(nz_in))
    ALLOCATE(var_in(my_nx, my_ny, nz_in))
    CALL mckpp_netcdf_get_var(routine, file, ncid, "zsal", z_in)
    CALL mckpp_netcdf_get_var(routine, file, ncid, "sal", var_in, start)
    CALL mckpp_netcdf_close(routine, file, ncid)
  
    CALL MCKPP_INITIALIZE_OCEAN_PROFILES_VINTERP(var_in,z_in,nz_in,kpp_const_fields%zm,temp_init(:,:,:))
  ENDIF ! End of masterproc section
  CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_init,init_chunk(1,begchunk,1))
  DO ichnk=begchunk,endchunk
    ncol=get_ncols_p(ichnk)
    kpp_3d_fields(ichnk)%X(1:ncol,1:NZP1,2)=init_chunk(1:ncol,ichnk,1:NZP1)
  ENDDO
  CALL mckpp_print(routine, "Initialized salinity") 
 
  ! Calculate and remove reference salinity
  
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
  
END SUBROUTINE mckpp_initialize_ocean_profiles


SUBROUTINE mckpp_initialize_ocean_profiles_vinterp(var_in,var_z,nz_in,model_z,var_out)
 
  INTEGER, intent(in) :: nz_in

  REAL(r8),dimension(PLON,PLAT,NZ_IN), intent(in) :: var_in
  REAL(r8),dimension(PLON,PLAT,NZP1), intent(out) :: var_out
 
  INTEGER :: my_nx, my_ny
  REAL,dimension(NZ_IN), intent(in) :: var_z
  REAL,dimension(NZP1), intent(in) :: model_z
  REAL :: deltaz,deltavar
  INTEGER :: ix,iy,iz,ipt,kin

  my_nx=PLON
  my_ny=PLAT

  DO iy=1,my_ny
     DO ix=1,my_nx
        ipt=(iy-1)*my_nx+ix
        kin=1
        DO iz=1,NZP1
           IF (model_z(iz) .GT. var_z(1)) THEN
              var_out(ix,iy,iz)=var_in(ix,iy,1)
           ELSEIF (model_z(iz) .LT. var_z(nz_in)) THEN
              var_out(ix,iy,iz)=var_in(ix,iy,nz_in)
           ELSE
              DO WHILE (var_z(kin+1) .GT. model_z(iz))
                 kin=kin+1
              ENDDO
              deltaz=var_z(kin)-var_z(kin+1)
              deltavar=var_in(ix,iy,kin)-var_in(ix,iy,kin+1)
              var_out(ix,iy,iz)=var_in(ix,iy,kin)+deltavar*(model_z(iz)-var_z(kin))/deltaz
           ENDIF
        ENDDO
     ENDDO
  ENDDO

END SUBROUTINE mckpp_initialize_ocean_profiles_vinterp

END MODULE mckpp_initialize_ocean_profiles_mod
#endif
