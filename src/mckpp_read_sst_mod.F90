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

  IMPLICIT NONE

  PUBLIC mckpp_initialize_sst, mckpp_read_sst

  PRIVATE

  INTEGER :: sst_nx, sst_ny 
  INTEGER, DIMENSION(3) :: start, count
  REAL, DIMENSION(:), ALLOCATABLE :: time_in
  REAL :: first_timein, last_timein
  LOGICAL :: file_open = .FALSE.
 
CONTAINS

SUBROUTINE mckpp_initialize_sst()

  INTEGER :: ncid, ntime_in
  CHARACTER(LEN=max_nc_filename_len) :: file

  CHARACTER(LEN=20) :: routine = "MCKPP_INITIALIZE_SST"
  CHARACTER(LEN=max_message_len) :: message
  
#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif

#ifdef MCKPP_COUPLE
  sst_nx = nx_globe 
  sst_ny = ny_globe 
#else
  sst_nx = nx
  sst_ny = ny 
#endif

  file = kpp_const_fields%sst_file
  WRITE(message,*) "Initializing ", TRIM(file)
  CALL mckpp_print(routine, message)
  CALL mckpp_netcdf_open(routine, file, ncid)
  file_open = .TRUE. 
  
  ! Work out start and count for each time entry
  count=(/sst_nx,sst_ny,1/)
  start=(/1,1,1/)
#ifndef MCKPP_COUPLE
  CALL mckpp_netcdf_determine_boundaries(routine, file, ncid, &
       kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), start(1), start(2))
#endif
  
  ! Read in time field
  CALL mckpp_netcdf_get_coord(routine, file, ncid, "time", ntime_in)
  ALLOCATE(time_in(ntime_in)) 
  CALL mckpp_netcdf_get_var(routine, file, ncid, "time", time_in)
  first_timein = time_in(1)
  last_timein = time_in(ntime_in)

#ifdef MCKPP_CAM3
  ENDIF ! End of masterproc section
#endif
  
  ! Now read first entry
  ! This should close file at end
  CALL mckpp_read_sst()
 
END SUBROUTINE mckpp_initialize_sst


SUBROUTINE mckpp_read_sst()

  IMPLICIT NONE

  CHARACTER(LEN=max_nc_filename_len) :: file
  REAL :: time, file_time, sstclim_time
  INTEGER :: ncid, ix, iy, ipt
  REAL :: offset_sst
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: var_in

#ifdef MCKPP_CAM3
  REAL(r8) :: sst_temp(PLON,PLAT), sst_chunk(PCOLS,begchunk:endchunk)
  INTEGER :: ichnk,ncol,icol
#endif
 
  CHARACTER(LEN=14) :: routine = "MCKPP_READ_SST"
  CHARACTER(LEN=max_message_len) :: message
  
#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif
     
  ALLOCATE( var_in(sst_nx,sst_ny,1) )

  file = kpp_const_fields%sst_file
  IF (.NOT. file_open) THEN
     CALL mckpp_netcdf_open(routine, file, ncid)
     file_open = .TRUE. 
  ENDIF 

  ! Work out time to read and check against times in file 
  sstclim_time = kpp_const_fields%time+0.5*kpp_const_fields%dto/kpp_const_fields%spd*kpp_const_fields%ndtupdsst
  IF (sstclim_time .gt. last_timein) THEN
     IF (kpp_const_fields%L_PERIODIC_CLIMSST) THEN
        DO WHILE (sstclim_time .gt. last_timein)
           sstclim_time=sstclim_time-kpp_const_fields%climsst_period
        ENDDO
     ELSE
        WRITE(message,*) "Time for which to read SST exceeds the last time in the netCDF file ", & 
            "and L_PERIODIC_CLIMSST has not been specified."
        CALL mckpp_print_error(routine, message) 
        WRITE(message,*) "Attempting to read SST will lead to an error, so aborting now ..."
        CALL mckpp_print_error(routine, message) 
        CALL MCKPP_ABORT()
     ENDIF
  ENDIF
  WRITE(message,*) 'Reading climatological SST for time ', sstclim_time
  CALL mckpp_print(routine, message)
  start(3) = NINT((sstclim_time-first_timein)*kpp_const_fields%spd/&
       (kpp_const_fields%dto*kpp_const_fields%ndtupdsst))+1
  WRITE(message,*) 'Reading climatological SST from position ',start(3)
  CALL mckpp_print(routine, message)

  file_time = time_in(start(3)) 
  IF (abs(file_time-sstclim_time) .GT. 0.01*kpp_const_fields%dtsec/kpp_const_fields%spd) THEN
     WRITE(message,*) 'Cannot find time,',  sstclim_time, 'in SST climatology file'
     CALL mckpp_print_error(routine, message) 
     WRITE(message,*) 'The closest I came was', file_time
     CALL mckpp_print_error(routine, message) 
     CALL MCKPP_ABORT()
  ENDIF

  CALL mckpp_netcdf_get_var(routine, file, ncid, "sst", var_in, start) 
 
  CALL mckpp_netcdf_close(routine, file, ncid)
  file_open = .FALSE.
  
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
