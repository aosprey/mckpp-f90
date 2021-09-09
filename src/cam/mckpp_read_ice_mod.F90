#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>

MODULE mckpp_read_ice_mod

  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_types,only: kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_parameters, ONLY: nx, ny, nx_globe, ny_globe
  USE mckpp_time_control, ONLY: mckpp_get_update_time

  IMPLICIT NONE

  PUBLIC mckpp_read_ice

  PRIVATE

  LOGICAL :: l_initialized = .FALSE. 
  INTEGER :: ice_nx, ice_ny, num_times
  INTEGER, DIMENSION(3) :: start, count
  REAL, DIMENSION(:), ALLOCATABLE :: file_times

CONTAINS


  SUBROUTINE initialize_ice(file, ncid)

    CHARACTER(LEN=max_nc_filename_len), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: ncid

    CHARACTER(LEN=14) :: routine = "INITIALIZE_ICE"
    CHARACTER(LEN=max_message_len) :: message

    IF (masterproc) THEN
      WRITE(message,*) "Initializing ", TRIM(file)
      CALL mckpp_print(routine, message)

      ice_nx=NX_GLOBE
      ice_ny=NY_GLOBE

      ! Work out start and count for each time entry
      count=(/ice_nx,ice_ny,1/)
      start=(/1,1,1/)

      ! Read in time field
      ALLOCATE(file_times(num_times)) 
      CALL mckpp_netcdf_get_var(routine, file, ncid, "t", file_times)

    ENDIF ! End of masterproc section
    l_initialized = .TRUE.

  END SUBROUTINE initialize_ice


  SUBROUTINE mckpp_read_ice()

    ! Read in ice concentrations from a user-provided netCDF file.
    ! Called _only_ if L_CLIMICE is TRUE in 3D_ocn.nml.
    ! Probably only necessary when coupling the model to an atmosphere
    ! that requires realistic ice concentrations.
    ! Written by Nick Klingaman, 11/01/08.

    IMPLICIT NONE

    REAL(r8) :: ice_temp(PLON,PLAT), ice_chunk(PCOLS,begchunk:endchunk)
    INTEGER :: ichnk,ncol,ico
    CHARACTER(LEN=max_nc_filename_len) :: file
    REAL :: update_time, max_ice, min_ice
    INTEGER :: ncid
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: var_in

    CHARACTER(LEN=14) :: routine = "MCKPP_READ_ICE"
    CHARACTER(LEN=max_message_len) :: message

    IF (masterproc) THEN
      file = kpp_const_fields%ice_file
      CALL mckpp_netcdf_open(routine, file, ncid)

      ! On first call, get file dimensions
      IF (.NOT. l_initialized) CALL initialize_ice(file, ncid)
      ALLOCATE( var_in(ice_nx,ice_ny,1) )

      ! Work out time to read and check against times in file
      CALL mckpp_get_update_time(file, kpp_const_fields%time, kpp_const_fields%ndtupdice, &
          file_times, num_times, kpp_const_fields%L_PERIODIC_CLIMICE, kpp_const_fields%climice_period, &
          update_time, start(3), method=1)
      WRITE(message,*) 'Reading climatological ice concentration for time ', update_time
      CALL mckpp_print(routine, message)

      ! Read and process data
      CALL mckpp_netcdf_get_var(routine, file, ncid, "iceconc", var_in, start) 

      ice_temp=var_in(:,:,1)
    ENDIF ! End of masterproc section
    CALL scatter_field_to_chunk(1,1,1,PLON,ice_temp,ice_chunk(1,begchunk))
    DO ichnk=begchunk,endchunk
      ncol=get_ncols_p(ichnk)
      kpp_3d_fields(ichnk)%iceconc(1:ncol)=ice_chunk(1:ncol,ichnk)
    ENDDO
    kpp_3d_fields%iceconc(:,:) = var_in(:,:,1)

    IF (kpp_const_fields%L_CLIM_ICE_DEPTH) THEN 
      IF (masterproc) THEN
        WRITE(message,*) 'Reading climatological ice depth for time ', update_time
        CALL mckpp_print(routine, message)
        CALL mckpp_netcdf_get_var(routine, file, ncid, "icedepth", var_in, start) 

        ice_temp=var_in(:,:,1)
      ENDIF ! End of masterproc section
      CALL scatter_field_to_chunk(1,1,1,PLON,ice_temp,ice_chunk(1,begchunk))
      DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%icedepth(1:ncol)=ice_chunk(1:ncol,ichnk)
      ENDDO
    ENDIF

    IF (kpp_const_fields%L_CLIM_SNOW_ON_ICE) THEN
      IF (masterproc) THEN
        WRITE(message,*) 'Reading climatological snow depth for time ', update_time
        CALL mckpp_print(routine, message)
        CALL mckpp_netcdf_get_var(routine, file, ncid, "snowdepth", var_in, start)

        CALL mckpp_netcdf_close(routine, file, ncid)
        ice_temp=var_in(:,:,1)
      ENDIF
      CALL scatter_field_to_chunk(1,1,1,PLON,ice_temp,ice_chunk(1,begchunk))
      DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%snowdepth(1:ncol)=ice_chunk(1:ncol,ichnk)
      ENDDO
    ENDIF
    
  END SUBROUTINE mckpp_read_ice

END MODULE mckpp_read_ice_mod
#endif
