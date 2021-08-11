#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif

MODULE mckpp_read_salinity_mod

#ifdef MCKPP_CAM3 
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_types,only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, get_ncols_p
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
#endif
  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_parameters, ONLY: nx, ny, nzp1, nx_globe, ny_globe
  USE mckpp_time_control, ONLY: mckpp_get_update_time

  IMPLICIT NONE

  PUBLIC mckpp_read_salinity_3d

  PRIVATE

  LOGICAL :: l_initialized = .FALSE. 
  INTEGER :: my_nx, my_ny, num_times
  INTEGER, DIMENSION(4) :: start, count
  REAL, DIMENSION(:), ALLOCATABLE :: file_times

CONTAINS

  SUBROUTINE initialize_salinity(file, ncid)

    CHARACTER(LEN=max_nc_filename_len), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: ncid

    REAL :: start_lat, start_lon
    INTEGER :: nz_in
    CHARACTER(LEN=19) :: routine = "INITIALIZE_SALINITY"
    CHARACTER(LEN=max_message_len) :: message

#ifdef MCKPP_CAM3
    IF (masterproc) THEN
#endif
      WRITE(message,*) "Initializing ", TRIM(file)
      CALL mckpp_print(routine, message)

#ifdef MCKPP_CAM3
      my_nx = nx_globe 
      my_ny = ny_globe 
#else
      my_nx = nx
      my_ny = ny 
#endif

      ! Work out start and count for each time entry
      count=(/my_nx,my_ny,NZP1,1/)
      start=(/1,1,1,1/)
#ifdef MCKPP_CAM3
      start_lon = kpp_global_fields%longitude(1)
      start_lat = kpp_global_fields%latitude(1)
#else
      start_lon = kpp_3d_fields%dlon(1)
      start_lat = kpp_3d_fields%dlat(1)
#endif
      CALL mckpp_netcdf_determine_boundaries(routine, file, ncid, &
          start_lon, start_lat, start(1), start(2), num_times)

      ! Read in time field
      ALLOCATE(file_times(num_times)) 
      CALL mckpp_netcdf_get_var(routine, file, ncid, "t", file_times)

      ! Check vertical levels
      CALL mckpp_netcdf_get_coord(routine, file, ncid, "z", nz_in)
      IF (nz_in .NE. nzp1) THEN
        WRITE(message,*) "Salinity climatology file does not have the correct ", &
            "number of vertical levels."
        CALL mckpp_print_error(routine, message)
        WRITE(message,*) "It should have ", NZP1, " but instead has ", nz_in
        CALL mckpp_print_error(routine, message)
        CALL mckpp_abort()
      END IF
#ifdef MCKPP_CAM3
    ENDIF ! End of masterproc section
#endif
    l_initialized = .TRUE.

  END SUBROUTINE initialize_salinity


  ! Read in a NetCDF file containing a 
  ! time-varying salinity field at every model vertical level.
  ! Frequency of read is controlled by ndtupdsal in the namelist
  ! NPK 12/02/08
  SUBROUTINE mckpp_read_salinity_3d()

#ifdef MCKPP_CAM3
    REAL(r8) :: sal_temp(PLON,PLAT,NZP1), sal_chunk(PCOLS,begchunk:endchunk,NZP1)
    INTEGER :: ichnk,icol,ncol
#endif
    CHARACTER(LEN=max_nc_filename_len) :: file
    REAL :: update_time
    INTEGER :: ncid, ix, iy, iz, ipt
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: var_in

    CHARACTER(LEN=22) :: routine = "MCKPP_READ_SALINITY_3D"
    CHARACTER(LEN=max_message_len) :: message

#ifdef MCKPP_CAM3
    IF (masterproc) THEN
#endif
      file = kpp_const_fields%sal_file
      CALL mckpp_netcdf_open(routine, file, ncid)

      ! On first call, get file dimensions
      IF (.NOT. l_initialized) CALL initialize_salinity(file, ncid)
      ALLOCATE( var_in(my_nx,my_ny,nzp1,1) )

      ! Work out time to read and check against times in file
      CALL mckpp_get_update_time(file, kpp_const_fields%time, kpp_const_fields%ndtupdsal, &
          file_times, num_times, kpp_const_fields%L_PERIODIC_SAL, kpp_const_fields%sal_period, &
          update_time, start(4), method=2)
      WRITE(message,*) 'Reading ocean temperature for time ', update_time
      CALL mckpp_print(routine, message)
      WRITE(message,*) 'Reading ocean temperature from position ',start(4)
      CALL mckpp_print(routine, message)

      ! Read data 
      CALL mckpp_netcdf_get_var(routine, file, ncid, "salinity", var_in, start) 
      CALL mckpp_netcdf_close(routine, file, ncid)
     
#ifdef MCKPP_CAM3
      sal_temp = var_in(:,:,:,1)
    ENDIF ! End of masterproc section
    CALL scatter_field_to_chunk(1,1,NZP1,PLON,sal_temp,sal_chunk(1,begchunk,1))
    DO ichnk=begchunk,endchunk
      ncol=get_ncols_p(ichnk)     
      ! Subtract reference salinity from climatology, for compatability with
      ! salinity values stored in main model.
      DO k=1,NZP1
        kpp_3d_fields(ichnk)%sal_clim(1:ncol,k)=sal_chunk(1:ncol,ichnk,k)-&
            kpp_3d_fields(ichnk)%Sref(1:ncol)
      ENDDO
    ENDDO
#else
    ! Put all (NX,NY) points into one long array with dimension NPTS.      
    DO ix = 1,NX
      DO iy = 1,NY
        ipt = (iy-1)*nx+ix
        DO iz = 1,NZP1
          ! Subtract reference salinity from climatology, for compatability with
          ! salinity values stored in main model.  
          kpp_3d_fields%sal_clim(ipt,iz) = var_in(ix,iy,iz,1) - kpp_3d_fields%Sref(ipt)
        ENDDO
      ENDDO
    ENDDO
#endif
    DEALLOCATE(var_in)

  END SUBROUTINE mckpp_read_salinity_3d

END MODULE mckpp_read_salinity_mod
