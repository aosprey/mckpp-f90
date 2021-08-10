#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif

MODULE mckpp_read_sst_mod

#ifdef MCKPP_CAM3  
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_types,only: kpp_3d_fields,kpp_const_fields
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
  USE mckpp_parameters, ONLY: nx, ny, nx_globe, ny_globe
  USE mckpp_time_control, ONLY: mckpp_get_update_time

  IMPLICIT NONE

  PUBLIC mckpp_read_sst

  PRIVATE

  LOGICAL :: l_initialized = .FALSE. 
  INTEGER :: sst_nx, sst_ny, num_times
  INTEGER, DIMENSION(3) :: start, count
  REAL, DIMENSION(:), ALLOCATABLE :: file_times

CONTAINS

  SUBROUTINE initialize_sst(file, ncid)

    CHARACTER(LEN=max_nc_filename_len), INTENT(IN) :: file
    INTEGER, INTENT(IN) :: ncid

    CHARACTER(LEN=20) :: routine = "MCKPP_INITIALIZE_SST"
    CHARACTER(LEN=max_message_len) :: message

#ifdef MCKPP_CAM3
    IF (masterproc) THEN
#endif
      WRITE(message,*) "Initializing ", TRIM(file)
      CALL mckpp_print(routine, message)
      
#ifdef MCKPP_COUPLE
      sst_nx = nx_globe 
      sst_ny = ny_globe 
#else
      sst_nx = nx
      sst_ny = ny 
#endif
      
      ! Work out start and count for each time entry
      count=(/sst_nx,sst_ny,1/)
      start=(/1,1,1/)
#ifndef MCKPP_COUPLE
      CALL mckpp_netcdf_determine_boundaries(routine, file, ncid, &
          kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), start(1), start(2))
#endif

      ! Read in time field
      CALL mckpp_netcdf_get_coord(routine, file, ncid, "time", num_times)
      ALLOCATE(file_times(num_times)) 
      CALL mckpp_netcdf_get_var(routine, file, ncid, "time", file_times)

#ifdef MCKPP_CAM3
    ENDIF ! End of masterproc section
#endif
    l_initialized = .TRUE.

  END SUBROUTINE initialize_sst


  SUBROUTINE mckpp_read_sst()

#ifdef MCKPP_CAM3
    REAL(r8) :: sst_temp(PLON,PLAT), sst_chunk(PCOLS,begchunk:endchunk)
    INTEGER :: ichnk,ncol,icol
#endif
    CHARACTER(LEN=max_nc_filename_len) :: file
    REAL :: update_time
    INTEGER :: ncid, ix, iy, ipt
    REAL :: offset_sst
    REAL, ALLOCATABLE, DIMENSION(:,:,:) :: var_in

    CHARACTER(LEN=14) :: routine = "MCKPP_READ_SST"
    CHARACTER(LEN=max_message_len) :: message

#ifdef MCKPP_CAM3
    IF (masterproc) THEN
#endif
      file = kpp_const_fields%sst_file
      CALL mckpp_netcdf_open(routine, file, ncid)

      ! On first call, get file dimensions
      IF (.NOT. l_initialized) CALL initialize_sst(file, ncid)
      ALLOCATE( var_in(sst_nx,sst_ny,1) )

      ! Work out time to read and check against times in file
      CALL mckpp_get_update_time(file, kpp_const_fields%time, kpp_const_fields%ndtupdsst, &
          file_times, num_times, kpp_const_fields%L_PERIODIC_CLIMSST, kpp_const_fields%climsst_period, &
          update_time, start(3))
      WRITE(message,*) 'Reading climatological SST for time ', update_time
      CALL mckpp_print(routine, message)
      WRITE(message,*) 'Reading climatological SST from position ',start(3)
      CALL mckpp_print(routine, message)

      ! Read data 
      CALL mckpp_netcdf_get_var(routine, file, ncid, "sst", var_in, start) 
      CALL mckpp_netcdf_close(routine, file, ncid)

      ! KPP expects temperatures in CELSIUS.  If climatological SSTs are
      ! in Kelvin, subtract 273.15.
      offset_sst = 0.
      DO ix = 1, sst_nx
        DO iy = 1, sst_ny
          IF (var_in(ix,iy,1) .gt. 200 .and. var_in(ix,iy,1) .lt. 400) &
              offset_sst = 273.15
        END DO
      END DO

#ifdef MCKPP_CAM3
      sst_temp=var_in(:,:,1)-offset_sst
    ENDIF ! End of masterproc section
    CALL scatter_field_to_chunk(1,1,1,PLON,sst_temp,sst_chunk(1,begchunk))
    DO ichnk=begchunk,endchunk
      ncol=get_ncols_p(ichnk)
      kpp_3d_fields(ichnk)%sst(1:ncol)=sst_chunk(1:ncol,ichnk)     
      IF (.NOT. kpp_const_fields%L_CLIMICE) kpp_3d_fields(ichnk)%iceconc(:)=0.0
      IF (.NOT. kpp_const_fields%L_CLIMCURR) THEN
        kpp_3d_fields(ichnk)%usf(:)=0.0
        kpp_3d_fields(ichnk)%vsf(:)=0.0
      ENDIF
    ENDDO
#else
    kpp_3d_fields%sst = var_in(:,:,1) - offset_sst          
    IF (.NOT. kpp_const_fields%L_CLIMICE) THEN
      kpp_3d_fields%iceconc = 0.0
    ENDIF
    IF (.NOT. kpp_const_fields%L_CLIMCURR) THEN
      kpp_3d_fields%usf = 0.0
      kpp_3d_fields%vsf = 0.0
    ENDIF
#endif

  END SUBROUTINE mckpp_read_sst

END MODULE mckpp_read_sst_mod
