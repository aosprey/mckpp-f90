MODULE mckpp_read_fluxes_mod

  IMPLICIT NONE
  
  PUBLIC mckpp_initialize_fluxes_file, mckpp_read_fluxes

  PRIVATE

  INTEGER, DIMENSION(3) :: start, count
  REAL, DIMENSION(:), ALLOCATABLE :: time_in

CONTAINS
  
! Work out and store start positions of data in file.
! Read and store the time dimension.
SUBROUTINE mckpp_initialize_fluxes_file()

  USE mckpp_data_fields, ONLY: kpp_const_fields, kpp_3d_fields
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, mckpp_netcdf_close, &
       mckpp_netcdf_determine_boundaries, mckpp_netcdf_get_coord, mckpp_netcdf_get_var
  USE mckpp_parameters, ONLY: nx, ny

  IMPLICIT NONE

  INTEGER :: ncid, ntime_in
  CHARACTER(LEN=max_nc_filename_len) :: file
  
  CHARACTER(LEN=28) :: routine = "MCKPP_INITIALIZE_FLUXES_FILE"
  CHARACTER(LEN=max_message_len) :: message

  file = kpp_const_fields%forcing_file
  WRITE(message,*) "Initializing ", TRIM(file)
  CALL mckpp_print(routine, message)
  CALL mckpp_netcdf_open(routine, file, ncid)

  ! Work out start and count for each time entry
  CALL mckpp_netcdf_determine_boundaries(routine, file, ncid, &
       kpp_3d_fields%dlon(1), kpp_3d_fields%dlat(1), start(1), start(2))
  start(3) = 1
  count = (/nx, ny, 1/)

  ! Read in time field
  CALL mckpp_netcdf_get_coord(routine, file, ncid, "time", ntime_in)
  ALLOCATE(time_in(ntime_in)) 
  CALL mckpp_netcdf_get_var(routine, file, ncid, "time", time_in)
  
  CALL mckpp_netcdf_close(routine, file, ncid)
              
END SUBROUTINE mckpp_initialize_fluxes_file


! Read flux data.
! This routine is called according to a specific update frequency, which should match spacing of data in
! the file. 
SUBROUTINE mckpp_read_fluxes(taux, tauy, swf, lwf, lhf, shf, rain, snow)

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, mckpp_netcdf_open, mckpp_netcdf_close, &
       mckpp_netcdf_get_var
  USE mckpp_parameters, ONLY: nx, ny, npts
  USE mckpp_time_control, ONLY: mckpp_get_update_time

  IMPLICIT NONE

  REAL, DIMENSION(npts), INTENT(OUT) :: taux, tauy, swf, lwf, lhf, shf, rain, snow
  
  CHARACTER(LEN=max_nc_filename_len) :: file  
  REAL :: time, file_time
  INTEGER :: ncid, ix, iy, ipt
  REAL, DIMENSION(nx, ny, 1) :: var_in
  
  CHARACTER(LEN=17) :: routine = "MCKPP_READ_FLUXES"
  CHARACTER(LEN=max_message_len) :: message

  file = kpp_const_fields%forcing_file
  CALL mckpp_netcdf_open(routine, file, ncid)
 
  ! Work out time to read in from 
  CALL mckpp_get_update_time(file, kpp_const_fields%time, kpp_const_fields%ndtocn, file_times, num_times, &
      .FALSE., 0, update_time, start(3))
  WRITE(message,*) 'Reading fluxes for time ', update_time
  CALL mckpp_print(routine, message)
  WRITE(message,*) 'Reading fluxes from time point ',start(3)
  CALL mckpp_print(routine, message)

  ! Read fileds
  CALL mckpp_netcdf_get_var(routine, file, ncid, "taux", taux, start, count)
  CALL mckpp_netcdf_get_var(routine, file, ncid, "tauy", tauy, start, count)   
  CALL mckpp_netcdf_get_var(routine, file, ncid, "swf", swf, start, count)   
  CALL mckpp_netcdf_get_var(routine, file, ncid, "lwf", lwf, start, count)   
  CALL mckpp_netcdf_get_var(routine, file, ncid, "lhf", lhf, start, count)   
  CALL mckpp_netcdf_get_var(routine, file, ncid, "shf", shf, start, count)   
  CALL mckpp_netcdf_get_var(routine, file, ncid, "precip", rain, start, count)   
  snow = 0.0
 
  CALL mckpp_netcdf_close(routine, file, ncid)

END SUBROUTINE mckpp_read_fluxes

END MODULE mckpp_read_fluxes_mod
