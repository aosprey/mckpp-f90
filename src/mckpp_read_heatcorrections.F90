#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif

MODULE mckpp_read_heat_corrections_mod

#ifdef MCKPP_CAM3
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_types,only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
#endif
  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, mckpp_netcdf_close, &
      mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_parameters, ONLY: nx, ny, nx_globe, ny_globe, nzp1
  USE mckpp_time_control, ONLY: mckpp_get_update_time

  IMPLICIT NONE

  PUBLIC mckpp_read_fcorr_2d, mckpp_read_fcorr_3d

  PRIVATE

  LOGICAL :: l_initialized = .FALSE. 
  INTEGER :: my_nx, my_ny, num_times
  INTEGER, DIMENSION(:), ALLOCATABLE :: start, count
  REAL, DIMENSION(:), ALLOCATABLE :: file_times

CONTAINS

  SUBROUTINE initialize_fcorr(file, ncid, ndims)

    CHARACTER(LEN=max_nc_filename_len), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: ncid, ndims

    REAL :: start_lat, start_lon
    INTEGER :: nz_in
    CHARACTER(LEN=22) :: routine = "MCKPP_INITIALIZE_FCORR"
    CHARACTER(LEN=max_message_len) :: message

    ! num_dims must be 2 or 3
    IF (ndims .LT. 2 .OR. ndims .GT. 3) THEN
      WRITE(message,*) "Correction must be 2d or 3d, but num_dims = ", ndims
      CALL mckpp_print_error(routine, message)
      CALL MCKPP_ABORT
    END IF

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
      ALLOCATE(start(ndims+1))
      ALLOCATE(count(ndims+2))
      start = 1
      IF (ndims .EQ. 2) THEN
        count = (/my_nx,my_ny,1/)
      ELSE IF (ndims .EQ. 3) THEN 
        count = (/my_nx,my_ny,nzp1,1/)
      END IF

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
      CALL mckpp_netcdf_get_var(routine, file, ncid, "time", file_times)

      ! Check vertical levels
      IF (ndims .EQ. 3) THEN
        CALL mckpp_netcdf_get_coord(routine, file, ncid, "z", nz_in)
        IF (nz_in .NE. nzp1) THEN
          WRITE(message,*) "Heat corrections file does not have the correct ", &
              "number of vertical levels."
          CALL mckpp_print_error(routine, message)
          WRITE(message,*) "It should have ", NZP1, " but instead has ", nz_in
          CALL mckpp_print_error(routine, message)
          CALL mckpp_abort()
        END IF
      END IF

#ifdef MCKPP_CAM3
    ENDIF ! End of masterproc section
#endif
    l_initialized = .TRUE.

  END SUBROUTINE initialize_fcorr


  ! Read in a NetCDF file containing a time-varying flux correction
  ! at the surface only.  Frequency of read is controlled by ndtupdfcorr
  ! in the namelist
  ! NPK 29/06/08
  SUBROUTINE mckpp_read_fcorr_2d()

#ifdef MCKPP_CAM3
    REAL (r8) :: fcorr_temp(PLON,PLAT), fcorr_chunk(PCOLS,begchunk:endchunk)
    INTEGER :: ichnk,icol,ncol
#endif
    CHARACTER(LEN=max_nc_filename_len) :: file
    REAL :: update_time
    INTEGER :: ncid, ix, iy, ipt
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: var_in

    CHARACTER(LEN=19) :: routine = "MCKPP_READ_FCORR_2D"
    CHARACTER(LEN=max_message_len) :: message

#ifdef MCKPP_CAM3
    IF (masterproc) THEN
#endif
      file = kpp_const_fields%fcorr_file
      CALL mckpp_netcdf_open(routine, file, ncid)

      ! On first call, get file dimensions
      IF (.NOT. l_initialized) CALL initialize_fcorr(file, ncid, ndims=2)
      ALLOCATE( var_in(my_nx,my_ny,1) )

      ! Work out time to read and check against times in file
      CALL mckpp_get_update_time(file, kpp_const_fields%time, kpp_const_fields%ndtupdfcorr, &
          file_times, num_times, kpp_const_fields%l_periodic_fcorr, kpp_const_fields%fcorr_period, &
          update_time, start(3), method=2)
      WRITE(message,*) 'Reading heat correction for time ', update_time
      CALL mckpp_print(routine, message)
      WRITE(message,*) 'Reading heat correction from position ',start(3)
      CALL mckpp_print(routine, message)

      ! Read data 
      CALL mckpp_netcdf_get_var(routine, file, ncid, "fcorr", var_in, start) 
      CALL mckpp_netcdf_close(routine, file, ncid)

#ifdef MCKPP_CAM3
      fcorr_temp = var_in(:,:,1)
    ENDIF ! End of masterproc section
    CALL scatter_field_to_chunk(1,1,1,PLON,fcorr_temp,fcorr_chunk(1,begchunk))
    DO ichnk=begchunk,endchunk
      ncol=get_ncols_p(ichnk)
      kpp_3d_fields(ichnk)%fcorr_twod(1:ncol)=fcorr_chunk(1:ncol,ichnk)
    ENDDO
#else    
    ! Put all (NX,NY) points into one long array with dimension NPTS.         
    DO ix = 1,NX
      DO iy = 1,NY
        ipt = (iy-1)*nx+ix
        kpp_3d_fields%fcorr_twod(ipt) = var_in(ix,iy,1)
      ENDDO
    ENDDO
#endif
    DEALLOCATE(var_in)

  END SUBROUTINE mckpp_read_fcorr_2d

  ! Read in a NetCDF file containing a 
  ! time-varying flux correction at every model vertical level.
  ! Frequency of read is controlled by ndtupdfcorr in the namelist
  ! NPK 12/02/08
  SUBROUTINE mckpp_read_fcorr_3d()

#ifdef MCKPP_CAM3
    REAL(r8) :: fcorr_temp(PLON,PLAT,NZP1), fcorr_chunk(PCOLS,begchunk:endchunk,NZP1)
    INTEGER :: icol,ncol,ichnk
#endif
    CHARACTER(LEN=max_nc_filename_len) :: file
    REAL :: update_time
    INTEGER :: ncid, ix, iy, iz, ipt
    REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: var_in

    CHARACTER(LEN=19) :: routine = "MCKPP_READ_FCORR_2D"
    CHARACTER(LEN=max_message_len) :: message

#ifdef MCKPP_CAM3
    IF (masterproc) THEN
#endif
      file = kpp_const_fields%fcorr_file
      CALL mckpp_netcdf_open(routine, file, ncid)

      ! On first call, get file dimensions
      IF (.NOT. l_initialized) CALL initialize_fcorr(file, ncid, ndims=3)
      ALLOCATE( var_in(my_nx,my_ny,nzp1,1) )

      ! Work out time to read and check against times in file
      CALL mckpp_get_update_time(file, kpp_const_fields%time, kpp_const_fields%ndtupdfcorr, &
          file_times, num_times, kpp_const_fields%l_periodic_fcorr, kpp_const_fields%fcorr_period, &
          update_time, start(4), method=1)
      WRITE(message,*) 'Reading heat correction for time ', update_time
      CALL mckpp_print(routine, message)
      WRITE(message,*) 'Reading heat correction from position ',start(3)
      CALL mckpp_print(routine, message)

      ! Read data 
      CALL mckpp_netcdf_get_var(routine, file, ncid, "fcorr", var_in, start) 
      CALL mckpp_netcdf_close(routine, file, ncid)

#ifdef MCKPP_CAM3
      fcorr_temp = var_in(:,:,:,1)
    ENDIF ! End of masterproc section  
    CALL scatter_field_to_chunk(1,1,NZP1,PLON,fcorr_temp,fcorr_chunk(1,begchunk,1))
    DO ichnk=begchunk,endchunk
      ncol=get_ncols_p(ichnk)
      kpp_3d_fields(ichnk)%fcorr_withz(1:ncol,1:NZP1) = fcorr_chunk(1:ncol,ichnk,1:NZP1)
    ENDDO
#else
    ! Put all (NX,NY) points into one long array with dimension NPTS.         
    DO ix = 1,NX
      DO iy = 1,NY
        ipt = (iy-1)*nx+ix
        DO iz = 1,NZP1
          kpp_3d_fields%fcorr_withz(ipt,iz) = var_in(ix,iy,iz,1)
        ENDDO
      ENDDO
    ENDDO
#endif
    DEALLOCATE(var_in)

  END SUBROUTINE mckpp_read_fcorr_3d

END MODULE mckpp_read_heat_corrections_mod
>>>>>>> Rewrite read_heatcorrections module.:mc-kpp_read_heatcorrections.F90
