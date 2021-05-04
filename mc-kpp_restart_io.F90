#ifndef MCKPP_CAM3
SUBROUTINE MCKPP_RESTART_IO_READ()

  USE mckpp_data_fields, ONLY: kpp_3d_fields,kpp_const_fields
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_parameters, ONLY: npts, nzp1
  
  IMPLICIT NONE
  
  CHARACTER(LEN=21) :: routine = "MCKPP_RESTART_IO_READ"
  CHARACTER(LEN=max_message_len) :: message
  
  WRITE(message,*) 'Total number of points = ',REAL(NPTS)*REAL(NZP1)
  CALL mckpp_print(routine, message)
  IF ( REAL(NPTS)*REAL(NZP1) .LT. 3000000. ) THEN   
     OPEN(30,FILE=kpp_const_fields%restart_infile,status='unknown',form='unformatted')
     READ(30) kpp_const_fields%time,kpp_3d_fields%U,kpp_3d_fields%X,kpp_3d_fields%CP,&
          kpp_3d_fields%rho,kpp_3d_fields%hmix,kpp_3d_fields%kmix,kpp_3d_fields%Sref,&
          kpp_3d_fields%SSref,kpp_3d_fields%Ssurf,kpp_3d_fields%Tref,kpp_3d_fields%old,&
          kpp_3d_fields%new,kpp_3d_fields%Us,kpp_3d_fields%Xs,kpp_3d_fields%hmixd
     CLOSE(30)
  ELSE
     OPEN(30,FILE=kpp_const_fields%restart_infile//'.1',status='unknown',form='unformatted')
     READ(30) kpp_const_fields%time,kpp_3d_fields%U,kpp_3d_fields%X,kpp_3d_fields%CP,&
          kpp_3d_fields%rho,kpp_3d_fields%hmix,kpp_3d_fields%kmix,kpp_3d_fields%Sref,&
          kpp_3d_fields%SSref,kpp_3d_fields%Ssurf,kpp_3d_fields%Tref,kpp_3d_fields%old,&
          kpp_3d_fields%new
     CLOSE(30)
     OPEN(31,FILE=kpp_const_fields%restart_infile//'.2',status='unknown',form='unformatted')
     READ(31) kpp_3d_fields%Us,kpp_3d_fields%Xs,kpp_3d_fields%hmixd
     CLOSE(31)
  ENDIF

  IF (abs(kpp_const_fields%time-kpp_const_fields%startt) .GT. 1.e-4) THEN 
     WRITE(message,*) 'Start time doesn''t match the restart record'
     CALL mckpp_print_error(routine, message) 
     WRITE(message,*) 'Start time in restart record = ',kpp_const_fields%time
     CALL mckpp_print_error(routine, message) 
     WRITE(message,*) 'Start time in namelist = ',kpp_const_fields%startt
     CALL mckpp_print_error(routine, message) 
     CALL MCKPP_ABORT()
  ENDIF
  
END SUBROUTINE MCKPP_RESTART_IO_READ


SUBROUTINE MCKPP_RESTART_IO_WRITE()
 
  USE mckpp_data_fields, ONLY: kpp_3d_fields,kpp_const_fields
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: npts, nzp1
     
  IMPLICIT NONE
  
  CHARACTER(LEN=22) :: routine = "MCKPP_RESTART_IO_WRITE"
  CHARACTER(LEN=max_message_len) :: message
    
  ! When the number of points in the model (NX*NY*NZP1) becomes
  ! quite large, we exceed the maximum size for Fortran unformatted
  ! binary files on certain machines.  The IF test below
  ! works around this by splitting the restart file in two.
  ! %Us and %Xs are the largest fields, so they get their own file.

  WRITE(message,*) 'Total number of points = ',REAL(NPTS)*REAL(NZP1)
  CALL mckpp_print(routine, message)
  IF ( REAL(NPTS)*REAL(NZP1) .LT. 3000000. ) THEN
     OPEN(31,FILE=kpp_const_fields%restart_outfile,status='unknown',form='unformatted')
     WRITE(31) kpp_const_fields%time,kpp_3d_fields%U,kpp_3d_fields%X,kpp_3d_fields%CP,&
          kpp_3d_fields%rho,kpp_3d_fields%hmix,kpp_3d_fields%kmix,kpp_3d_fields%Sref,&
          kpp_3d_fields%SSref,kpp_3d_fields%Ssurf,kpp_3d_fields%Tref,kpp_3d_fields%old,&
          kpp_3d_fields%new,kpp_3d_fields%Us,kpp_3d_fields%Xs,kpp_3d_fields%hmixd
     CLOSE(31)
  ELSE
     OPEN(31,FILE=kpp_const_fields%restart_outfile//'.1',status='unknown',form='unformatted')
     WRITE(31) kpp_const_fields%time,kpp_3d_fields%U,kpp_3d_fields%X,kpp_3d_fields%CP,&
          kpp_3d_fields%rho,kpp_3d_fields%hmix,kpp_3d_fields%kmix,kpp_3d_fields%Sref,&
          kpp_3d_fields%SSref,kpp_3d_fields%Ssurf,kpp_3d_fields%Tref,kpp_3d_fields%old,&
          kpp_3d_fields%new
     CLOSE(31)
     OPEN(32,FILE=kpp_const_fields%restart_outfile//'.2',status='unknown',form='unformatted')
     WRITE(32) kpp_3d_fields%Us,kpp_3d_fields%Xs,kpp_3d_fields%hmixd
     CLOSE(32)
  ENDIF
  
  RETURN
END SUBROUTINE MCKPP_RESTART_IO_WRITE
#endif


#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif
SUBROUTINE MCKPP_RESTART_IO_WRITE_NETCDF()
  
#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE shr_kind_mod, only: r8=>shr_kind_r8, r4=>shr_kind_r4
  USE ppgrid, only: pcols,begchunk,endchunk
  USE phys_grid,only: get_ncols_p,gather_chunk_to_field
  USE pmgrid, only: masterproc
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
#endif
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, npts, nzp1, nzp1tmax

  IMPLICIT NONE
#include <netcdf.inc>

#ifdef MCKPP_CAM3
  INTEGER :: icol,ncol,ichnk
  REAL(r8) :: temp_global_3d(PLON,PLAT,NZP1,2),temp_global_2d(PLON,PLAT),&
       temp_chunk(PCOLS,begchunk:endchunk,NZP1)
#else
  LOGICAL, ALLOCATABLE :: ALL_OCEAN(:)
#endif
  
  INTEGER :: my_nx, my_ny
  CHARACTER(LEN=20) :: netcdf_restart_outfile
  INTEGER, parameter :: nvars=20, ndims=6
  INTEGER, dimension(nvars) :: varids
  INTEGER, dimension(ndims) :: dim_dimids,dim_varids
  INTEGER :: ncid,ivar,status,i
  REAL*4, allocatable, dimension(:,:,:,:) :: temp_output
  !REAL*4,dimension(my_nx,my_ny,1:NZP1tmax+1,2) :: temp_output_large
  REAL*4, allocatable, dimension(:,:) :: temp_twod
  !REAL*4,dimension(NPTS,1:NZP1tmax+1) :: temp_twod_large
  REAL*4, allocatable :: lon_out(:),lat_out(:),z_out(:)
  REAL*4 :: delta

  CHARACTER(LEN=29) :: routine = "MCKPP_RESTART_IO_WRITE_NETCDF"
  CHARACTER(LEN=max_message_len) :: message
  

#ifdef MCKPP_CAM3
  IF (masterproc) THEN

  my_nx = PLON
  my_ny = PLAT
#else
  my_nx = nx
  my_ny = ny
  ALLOCATE( ALL_OCEAN(NPTS) ) 
  ! Make a fake land/sea mask that is all ocean, to avoid land points being masked out in the restart file
  ALL_OCEAN(1:NPTS) = .TRUE.
#endif
  ALLOCATE( temp_output(my_nx,my_ny,1:NZP1,2) ) 
  ALLOCATE( temp_twod(NPTS,1:NZP1) ) 
  ALLOCATE( lon_out(my_nx) ) 
  ALLOCATE( lat_out(my_ny) ) 
  ALLOCATE( z_out(1:NZP1tmax+1) ) 

  WRITE(netcdf_restart_outfile,'(A17,A3)') kpp_const_fields%restart_outfile,'.nc'

  status=NF_CREATE(netcdf_restart_outfile,NF_CLOBBER,ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR
  WRITE(message,*) 'Created file ', netcdf_restart_outfile
  CALL mckpp_print(routine, message)
#ifdef MCKPP_CAM3
  IF (my_nx .gt. 1) delta=kpp_global_fields%longitude(2)-kpp_global_fields%longitude(1)
#else
  IF (my_nx .gt. 1) delta=kpp_3d_fields%dlon(2)-kpp_3d_fields%dlon(1)
#endif
  CALL MCKPP_NCDF_DEF_DIM(ncid,dim_dimids(1),my_nx,dim_varids(1),'longitude','deg',delta,' ')
#ifdef MCKPP_CAM3
  IF (my_ny .gt. 1) delta=kpp_global_fields%latitude(2)-kpp_global_fields%latitude(1)
#else
  IF (my_ny .gt. 1) delta=kpp_3d_fields%dlat(NX+1)-kpp_3d_fields%dlat(1)
#endif
  CALL MCKPP_NCDF_DEF_DIM(ncid,dim_dimids(2),my_ny,dim_varids(2),'latitude','deg',delta,' ')
  CALL MCKPP_NCDF_DEF_DIM(ncid,dim_dimids(3),NZP1,dim_varids(3),'z','m',0.0,' ')
  CALL MCKPP_NCDF_DEF_DIM(ncid,dim_dimids(4),NZP1tmax+1,dim_varids(4),'z_1','number_of_gridpoint',0.0,' ')
  CALL MCKPP_NCDF_DEF_DIM(ncid,dim_dimids(5),1,dim_varids(5),'t','days',0.0,' ')
  CALL MCKPP_NCDF_DEF_DIM(ncid,dim_dimids(6),2,dim_varids(6),'intcnt','unitless',0.0,' ')

  DO ivar=1,nvars
     SELECT CASE(ivar)
     CASE (1)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),1,(/dim_dimids(5)/),'time','days','time')
     CASE (2)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),3,(/dim_dimids(1),dim_dimids(2),dim_dimids(3)/),'uvel','m/s',&
             'zonal velocity')
     CASE (3)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),3,(/dim_dimids(1),dim_dimids(2),dim_dimids(3)/),'vvel','m/s',&
             'meridional velocity')
     CASE (4)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),3,(/dim_dimids(1),dim_dimids(2),dim_dimids(3)/),'T',&
             'degrees Celsius','temperature')
     CASE (5)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),3,(/dim_dimids(1),dim_dimids(2),dim_dimids(3)/),'S',&
             'o/oo','salinity')
     CASE (6)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),3,(/dim_dimids(1),dim_dimids(2),dim_dimids(3)/),'CP',&
             'J/kg/K','specific heat capacity')
     CASE (7)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),3,(/dim_dimids(1),dim_dimids(2),dim_dimids(3)/),'rho',&
             'kg/m^3','density')
     CASE (8)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),2,(/dim_dimids(1),dim_dimids(2)/),'hmix','m',&
             'mixed-layer depth')
     CASE (9)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),2,(/dim_dimids(1),dim_dimids(2)/),'kmix','m',&
             'mixed-layer depth')
     CASE (10)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),2,(/dim_dimids(1),dim_dimids(2)/),'Sref','o/oo',&
             'reference salinity')
     CASE (11)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),2,(/dim_dimids(1),dim_dimids(2)/),'SSref','o/oo',&
             'reference surface salinity')
     CASE (12)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),2,(/dim_dimids(1),dim_dimids(2)/),'Ssurf','o/oo',&
             'surface salinity')
     CASE (13)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),2,(/dim_dimids(1),dim_dimids(2)/),'Tref','o/oo',&
             'surface salinity')
     CASE (14)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),2,(/dim_dimids(1),dim_dimids(2)/),'old','unitless',&
             'integration counter')
     CASE (15)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),2,(/dim_dimids(1),dim_dimids(2)/),'new','unitless',&
             'integration counter')
     CASE (16)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),4,(/dim_dimids(1),dim_dimids(2),dim_dimids(3),dim_dimids(6)/),'Us','m/s',&
             'zonal velocity in integration')
     CASE (17)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),4,(/dim_dimids(1),dim_dimids(2),dim_dimids(3),dim_dimids(6)/),'Vs','m/s',&
             'meridional velocity in integration')
     CASE (18)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),4,(/dim_dimids(1),dim_dimids(2),dim_dimids(3),dim_dimids(6)/),'Ts','m/s',&
             'Temperature in integration')
     CASE (19)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),4,(/dim_dimids(1),dim_dimids(2),dim_dimids(3),dim_dimids(6)/),'Ss','m/s',&
             'Salinity in integration')
     CASE (20)
        CALL MCKPP_NCDF_DEF_VAR(ncid,varids(ivar),3,(/dim_dimids(1),dim_dimids(2),dim_dimids(6)/),'hmixd','m','mixed-layer depth')
     END SELECT
  ENDDO
  
  status=NF_ENDDEF(ncid)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
#ifdef MCKPP_CAM3
  lon_out=kpp_global_fields%longitude
  lat_out=kpp_global_fields%latitude
#else
  lon_out=kpp_3d_fields%dlon(1:NX)
  lat_out=kpp_3d_fields%dlat(1::NX)
#endif
  status=NF_PUT_VAR_REAL(ncid,dim_varids(1),lon_out)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_PUT_VAR_REAL(ncid,dim_varids(2),lat_out)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)  
  z_out(1:NZP1)=kpp_const_fields%zm
  status=NF_PUT_VAR_REAL(ncid,dim_varids(3),z_out(1:NZP1))
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  DO i=1,NZP1tmax+1
     z_out(i)=i
  ENDDO
  status=NF_PUT_VAR_REAL(ncid,dim_varids(4),z_out(1:NZP1tmax+1))
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ! Set correct time for validity of restart fields 
  ! (end of this timestep = start of next timestep)
  z_out(1)=kpp_const_fields%time !+kpp_const_fields%dto/kpp_const_fields%spd
  status=NF_PUT_VAR_REAL(ncid,dim_varids(5),z_out(1))
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  z_out(1:2)=(/0,1/)
  status=NF_PUT_VAR_REAL(ncid,dim_varids(6),z_out(1:2))
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)

#ifdef MCKPP_CAM3
ENDIF ! End of masterproc section
#endif
  DO ivar=1,nvars
     SELECT CASE(ivar)
     CASE(1)
#ifdef MCKPP_CAM3
        IF (masterproc) THEN
#endif           
        z_out(1)=kpp_const_fields%time+kpp_const_fields%dto/kpp_const_fields%spd
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),z_out(1))
#ifdef MCKPP_CAM3
        ENDIF ! End of masterproc section
#endif
     CASE(2)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,1)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk,NZP1,1e20,temp_chunk)
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,:,1))
        IF (masterproc) THEN
           temp_output(:,:,1:NZP1,1)=temp_global_3d(:,:,1:NZP1,1)
#else
        temp_twod(:,1:NZP1)=kpp_3d_fields%U(:,:,1)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,1:NZP1),NZP1,ALL_OCEAN,2e20,temp_output(:,:,1:NZP1,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,1:NZP1,1))
#ifdef MCKPP_CAM3
        ENDIF ! End of masterproc section
#endif
     CASE(3)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,2)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk,NZP1,1e20,temp_chunk)
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,:,1))
        IF (masterproc) THEN
           temp_output(:,:,1:NZP1,1)=temp_global_3d(:,:,1:NZP1,1)
#else
        temp_twod(:,1:NZP1)=kpp_3d_fields%U(:,:,2)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,1:NZP1),NZP1,ALL_OCEAN,2e20,temp_output(:,:,1:NZP1,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,1:NZP1,1))
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(4)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%X(1:ncol,1:NZP1,1)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk,NZP1,1e20,temp_chunk)
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,:,1))
        IF (masterproc) THEN
           temp_output(:,:,1:NZP1,1)=temp_global_3d(:,:,1:NZP1,1)
#else
        temp_twod(:,1:NZP1)=kpp_3d_fields%X(:,:,1)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,1:NZP1),NZP1,ALL_OCEAN,2e20,temp_output(:,:,1:NZP1,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,1:NZP1,1))
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(5)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%X(1:ncol,1:NZP1,2)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk,NZP1,1e20,temp_chunk)
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d)
        IF (masterproc) THEN
           temp_output(:,:,1:NZP1,1)=temp_global_3d(:,:,1:NZP1,1)
#else
        temp_twod(:,1:NZP1)=kpp_3d_fields%X(:,:,2)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,1:NZP1),NZP1,ALL_OCEAN,2e20,temp_output(:,:,1:NZP1,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,1:NZP1,1))
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(6)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%cp(1:ncol,1:NZP1)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk,NZP1,1e20,temp_chunk)
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,1:NZP1,1))
        IF (masterproc) THEN
           temp_output(:,:,1:NZP1,1)=temp_global_3d(:,:,1:NZP1,1)
#else
        temp_twod=kpp_3d_fields%cp(:,1:NZP1)        
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,:),NZP1,ALL_OCEAN,2e20,&
             temp_output(:,:,:,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,:,1))
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(7)        
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%rho(1:ncol,1:NZP1)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk,NZP1,1e20,temp_chunk)
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,1:NZP1,1))
        IF (masterproc) THEN
           temp_output(:,:,1:NZP1,1)=temp_global_3d(:,:,1:NZP1,1)
#else
        temp_twod=kpp_3d_fields%rho(:,1:NZP1)        
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,:),NZP1,ALL_OCEAN,2e20,&
             temp_output(:,:,:,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,:,1))
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(8)
#ifdef MCKPP_CAM3        
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%hmix(1:ncol)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,1),1,1e20,temp_chunk(:,:,1))
        CALL gather_chunk_to_field(1,1,1,PLON,temp_chunk(1,begchunk,1),temp_global_2d(:,:))
        IF (masterproc) THEN
           temp_output(:,:,1,1)=temp_global_2d(:,:)
#else
        temp_twod(:,1)=kpp_3d_fields%hmix
        CALL MCKPP_REFORMAT_MASK_OUTPUT_1D(temp_twod(:,1),ALL_OCEAN,2e20,temp_output(:,:,1,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,1,1))
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(9)
#ifdef MCKPP_CAM3        
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%kmix(1:ncol)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,1),1,1e20,temp_chunk(:,:,1))
        CALL gather_chunk_to_field(1,1,1,PLON,temp_chunk(1,begchunk,1),temp_global_2d)
        IF (masterproc) THEN
           temp_output(:,:,1,1)=temp_global_2d
#else
        temp_twod(:,1)=kpp_3d_fields%kmix
        CALL MCKPP_REFORMAT_MASK_OUTPUT_1D(temp_twod(:,1),ALL_OCEAN,2e20,temp_output(:,:,1,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,1,1))
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(10)    
#ifdef MCKPP_CAM3        
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%Sref(1:ncol)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,1),1,1e20,temp_chunk(:,:,1))
        CALL gather_chunk_to_field(1,1,1,PLON,temp_chunk(1,begchunk,1),temp_global_2d)
        IF (masterproc) THEN
           temp_output(:,:,1,1)=temp_global_2d
#else
        temp_twod(:,1)=kpp_3d_fields%Sref
        CALL MCKPP_REFORMAT_MASK_OUTPUT_1D(temp_twod(:,1),ALL_OCEAN,2e20,temp_output(:,:,1,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,1,1))
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(11)        
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%SSref(1:ncol)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,1),1,1e20,temp_chunk(:,:,1))
        CALL gather_chunk_to_field(1,1,1,PLON,temp_chunk(1,begchunk,1),temp_global_2d)
        IF (masterproc) THEN
           temp_output(:,:,1,1)=temp_global_2d
#else
        temp_twod(:,1)=kpp_3d_fields%SSref
        CALL MCKPP_REFORMAT_MASK_OUTPUT_1D(temp_twod(:,1),ALL_OCEAN,2e20,temp_output(:,:,1,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,1,1))
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(12)      
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%Ssurf(1:ncol)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,1),1,1e20,temp_chunk(:,:,1))
        CALL gather_chunk_to_field(1,1,1,PLON,temp_chunk(1,begchunk,1),temp_global_2d)
        IF (masterproc) THEN
           temp_output(:,:,1,1)=temp_global_2d
#else
        temp_twod(:,1)=kpp_3d_fields%Ssurf
        CALL MCKPP_REFORMAT_MASK_OUTPUT_1D(temp_twod(:,1),ALL_OCEAN,2e20,temp_output(:,:,1,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,1,1))        
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(13)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%Tref(1:ncol)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,1),1,1e20,temp_chunk(:,:,1))
        CALL gather_chunk_to_field(1,1,1,PLON,temp_chunk(1,begchunk,1),temp_global_2d)
        IF (masterproc) THEN
           temp_output(:,:,1,1)=temp_global_2d
#else
        temp_twod(:,1)=kpp_3d_fields%Tref
        CALL MCKPP_REFORMAT_MASK_OUTPUT_1D(temp_twod(:,1),ALL_OCEAN,2e20,temp_output(:,:,1,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,1,1))  
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(14) 
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%old(1:ncol)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,1),1,1e20,temp_chunk(:,:,1))
        CALL gather_chunk_to_field(1,1,1,PLON,temp_chunk(1,begchunk,1),temp_global_2d)
        IF (masterproc) THEN
           temp_output(:,:,1,1)=temp_global_2d
#else           
        temp_twod(:,1)=kpp_3d_fields%old
        CALL MCKPP_REFORMAT_MASK_OUTPUT_1D(temp_twod(:,1),ALL_OCEAN,2e20,temp_output(:,:,1,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,1,1))  
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(15)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%old(1:ncol)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,1),1,1e20,temp_chunk(:,:,1))
        CALL gather_chunk_to_field(1,1,1,PLON,temp_chunk(1,begchunk,1),temp_global_2d)
        IF (masterproc) THEN
           temp_output(:,:,1,1)=temp_global_2d
#else
        temp_twod(:,1)=kpp_3d_fields%new
        CALL MCKPP_REFORMAT_MASK_OUTPUT_1D(temp_twod(:,1),ALL_OCEAN,2e20,temp_output(:,:,1,1))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,1,1))  
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(16)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%Us(1:ncol,1:NZP1,1,0)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,:),NZP1,1e20,temp_chunk(:,:,:))
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,1:NZP1,1))
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%Us(1:ncol,1:NZP1,1,1)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,:),NZP1,1e20,temp_chunk(:,:,:))
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,1:NZP1,2))
        IF (masterproc) THEN
           temp_output(:,:,:,:)=temp_global_3d(:,:,:,:)
#else
        temp_twod(:,:)=kpp_3d_fields%Us(:,:,1,0)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,:),NZP1,ALL_OCEAN,2e20,temp_output(:,:,:,1))
        temp_twod(:,:)=kpp_3d_fields%Us(:,:,1,1)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,:),NZP1,ALL_OCEAN,2e20,temp_output(:,:,:,2))        
#endif        
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,:,:))
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(17)        
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%Us(1:ncol,1:NZP1,2,0)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,:),NZP1,1e20,temp_chunk(:,:,:))
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,1:NZP1,1))
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%Us(1:ncol,1:NZP1,2,1)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,:),NZP1,1e20,temp_chunk(:,:,:))
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,1:NZP1,2))
        IF (masterproc) THEN
           temp_output(:,:,:,:)=temp_global_3d(:,:,:,:)
#else
        temp_twod(:,:)=kpp_3d_fields%Us(:,:,2,0)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,:),NZP1,ALL_OCEAN,2e20,temp_output(:,:,:,1))
        temp_twod(:,:)=kpp_3d_fields%Us(:,:,2,1)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,:),NZP1,ALL_OCEAN,2e20,temp_output(:,:,:,2))        
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,:,:))
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(18)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%Xs(1:ncol,1:NZP1,1,0)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,:),NZP1,1e20,temp_chunk(:,:,:))
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,1:NZP1,1))
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%Xs(1:ncol,1:NZP1,1,1)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,:),NZP1,1e20,temp_chunk(:,:,:))
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,1:NZP1,2))
        IF (masterproc) THEN
           temp_output(:,:,:,:)=temp_global_3d(:,:,:,:)
#else
        temp_twod(:,:)=kpp_3d_fields%Xs(:,:,1,0)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,:),NZP1,ALL_OCEAN,2e20,temp_output(:,:,:,1))
        temp_twod(:,:)=kpp_3d_fields%Xs(:,:,1,1)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,:),NZP1,ALL_OCEAN,2e20,temp_output(:,:,:,2))        
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,:,:))
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(19)        
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%Xs(1:ncol,1:NZP1,2,0)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,:),NZP1,1e20,temp_chunk(:,:,:))
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,1:NZP1,1))
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%Xs(1:ncol,1:NZP1,2,1)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,:),NZP1,1e20,temp_chunk(:,:,:))
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,1:NZP1,2))
        IF (masterproc) THEN
           temp_output(:,:,:,:)=temp_global_3d(:,:,:,:)
#else
        temp_twod(:,:)=kpp_3d_fields%Xs(:,:,2,0)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,:),NZP1,ALL_OCEAN,2e20,temp_output(:,:,:,1))
        temp_twod(:,:)=kpp_3d_fields%Xs(:,:,2,1)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_twod(:,:),NZP1,ALL_OCEAN,2e20,temp_output(:,:,:,2))        
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,:,:))
#ifdef MCKPP_CAM3
        ENDIF
#endif
     CASE(20)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%hmixd(1:ncol,0)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(1:ncol,:,1),1,1e20,temp_chunk(1:ncol,ichnk,1))
        CALL gather_chunk_to_field(1,1,1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,1,1))
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%hmixd(1:ncol,1)
        ENDDO
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(1:ncol,:,1),1,1e20,temp_chunk(1:ncol,ichnk,1))
        CALL gather_chunk_to_field(1,1,1,PLON,temp_chunk(1,begchunk,1),temp_global_3d(:,:,1,2))
        IF (masterproc) THEN
           temp_output(:,:,1,1:2)=temp_global_3d(:,:,1,1:2)
#else
        temp_twod(:,1)=kpp_3d_fields%hmixd(:,0)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_1D(temp_twod(:,1),ALL_OCEAN,2e20,temp_output(:,:,1,1))
        temp_twod(:,1)=kpp_3d_fields%hmixd(:,1)
        CALL MCKPP_REFORMAT_MASK_OUTPUT_1D(temp_twod(:,1),ALL_OCEAN,2e20,temp_output(:,:,1,2))
#endif
        status=NF_PUT_VAR_REAL(ncid,varids(ivar),temp_output(:,:,1,:))
#ifdef MCKPP_CAM3
        ENDIF ! End of masterproc section
#endif
     END SELECT
#ifdef MCKPP_CAM3
     IF (masterproc .and. status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
#else
     IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
#endif
  ENDDO

#ifdef MCKPP_CAM3 
  IF (masterproc) THEN
#endif
  status=NF_CLOSE(ncid)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
#ifdef MCKPP_CAM3
  ENDIF
#endif

END SUBROUTINE MCKPP_RESTART_IO_WRITE_NETCDF


SUBROUTINE MCKPP_RESTART_IO_READ_NETCDF()

#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE shr_kind_mod, only: r8=>shr_kind_r8, r4=>shr_kind_r4
  USE ppgrid, only: pcols,begchunk,endchunk
  USE phys_grid,only: get_ncols_p,scatter_field_to_chunk
#ifdef SPMD
  USE mpishorthand
#endif
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
#endif
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, npts, nzp1

  IMPLICIT NONE
#include <netcdf.inc>

#ifdef MCKPP_CAM3
  INTEGER :: icol,ncol,ichnk
  REAL(r8) :: temp_global_3d(PLON,PLAT,NZP1,2),temp_global_2d(PLON,PLAT),&
       temp_chunk(PCOLS,begchunk:endchunk,NZP1),double_time
#else
  LOGICAL, ALLOCATABLE :: ALL_OCEAN(:)  
#endif

  INTEGER :: my_nx, my_ny
  CHARACTER(LEN=20) :: netcdf_restart_infile
  INTEGER, parameter :: nvars=20, ndims=6
  INTEGER, dimension(nvars) :: varids
  INTEGER, dimension(ndims) :: dimids
  REAL*4, allocatable, dimension(:,:,:,:) :: temp_input
  REAL*4, allocatable, dimension(:,:) :: temp_twod
  INTEGER :: status,ncid,ivar,file_nlon,file_nlat,file_nz
  CHARACTER(LEN=5) :: varnames(nvars) = &
       (/'time ','uvel ','vvel ','T    ','S    ','CP   ',&
         'rho  ','hmix ','kmix ','Sref ','SSref','Ssurf',&
         'Tref ','old  ','new  ','Us   ','Vs   ','Ts   ','Ss   ','hmixd'/)

  CHARACTER(LEN=28) :: routine = "MCKPP_RESTART_IO_READ_NETCDF"
  CHARACTER(LEN=max_message_len) :: message

#ifdef MCKPP_CAM3
  IF (masterproc) THEN

  my_nx=PLON
  my_ny=PLAT
#else
  my_nx=nx
  my_ny=ny
  ALLOCATE( all_ocean(npts) ) 
  ALL_OCEAN(:)=.TRUE.
#endif  
  ALLOCATE( temp_input(my_nx,my_ny,1:NZP1,2) ) 
  ALLOCATE( temp_twod(my_nx*my_ny,1:NZP1) )
     
  netcdf_restart_infile = TRIM(ADJUSTL(kpp_const_fields%restart_infile)) // '.nc'
  status=NF_OPEN(netcdf_restart_infile,NF_NOWRITE,ncid)
  
  ! Conduct safety checks on the restart file for number of longitudes, latitudes and vertical points
  status=NF_INQ_DIMID(ncid,'longitude',dimids(1))
  status=NF_INQ_DIMLEN(ncid,dimids(1),file_nlon)
  IF (file_nlon .NE. my_nx) THEN
     WRITE(message,*) 'Number of longitudes in restart file is ', file_nlon, &
          ' but number of longitudes in the model is ' ,my_nx
     CALL mckpp_print_error(routine, message) 
     WRITE(message,*) 'Please correct either the restart file or the model in order to continue.'
     CALL mckpp_print_error(routine, message) 
     ! CALL MCKPP_ABORT()
  ENDIF
  status=NF_INQ_DIMID(ncid,'latitude',dimids(1))
  status=NF_INQ_DIMLEN(ncid,dimids(1),file_nlat)
  IF (file_nlat .NE. my_ny) THEN
     WRITE(message,*) 'Number of latitudes in restart file is ', file_nlat, &
         ' but number of latitudes in the model is ', my_ny
     CALL mckpp_print_error(routine, message) 
     WRITE(message,*) 'Please correct either the restart file or the model in order to continue.'
     CALL mckpp_print_error(routine, message) 
     CALL MCKPP_ABORT()
  ENDIF
  status=NF_INQ_DIMID(ncid,'z',dimids(1))
  status=NF_INQ_DIMLEN(ncid,dimids(1),file_nz)
  IF (file_nz .NE. NZP1) THEN
     WRITE(message,*) 'Number of depths in restart file is ', file_nz, &
         ' but number of depths in the model is ', NZP1
     CALL mckpp_print_error(routine, message) 
     WRITE(message,*) 'Please correct either the restart file or the model in order to continue.'
     CALL mckpp_print_error(routine, message) 
     CALL MCKPP_ABORT()
  ENDIF

  DO ivar=1,nvars
     status=NF_INQ_VARID(ncid,varnames(ivar),varids(ivar))
  ENDDO
  
#ifdef MCKPP_CAM3
ENDIF ! End of masterproc section
#endif

DO ivar=1,nvars
   SELECT CASE(ivar)
   CASE (1)
#ifdef MCKPP_CAM3
      IF (masterproc) THEN
#endif
      status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(1,1,1,1))
#ifdef MCKPP_CAM3
      double_time=FLOOR(temp_input(1,1,1,1))
      ENDIF ! End of masterproc section
      CALL MPIBCAST(double_time,1,mpir8,0,mpicom)
      kpp_const_fields%time=double_time
#else
      kpp_const_fields%time=temp_input(1,1,1,1)
#endif
#ifdef MCKPP_CAM3
      IF (masterproc) THEN
#endif
      IF (ABS(kpp_const_fields%time-kpp_const_fields%startt) .GT. 1e-4) THEN
         WRITE(message,*) 'Start time does not match the restart record.'
         CALL mckpp_print_error(routine, message) 
         WRITE(message,*) 'Start time in restart record = ',kpp_const_fields%time
         CALL mckpp_print_error(routine, message) 
         WRITE(message,*) 'Start time in namelist = ',kpp_const_fields%startt
         CALL mckpp_print_error(routine, message) 
         ! CALL MCKPP_ABORT
      ENDIF
#ifdef MCKPP_CAM3
      ENDIF ! End of masterproc section
#endif
   CASE (2)
#ifdef MCKPP_CAM3
      IF (masterproc) THEN
#endif
      status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,:,1))
#ifdef MCKPP_CAM3
      temp_global_3d(:,:,1:NZP1,1)=temp_input(:,:,:,1)
      ENDIF ! End of masterproc section
      CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,1),temp_chunk(1,begchunk,1))
      DO ichnk=begchunk,endchunk
         ncol=get_ncols_p(ichnk)
         kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,1)=temp_chunk(1:ncol,ichnk,1:NZP1)
      ENDDO
#else
      CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,1),NZP1,ALL_OCEAN,2e20,temp_twod)
      kpp_3d_fields%U(:,:,1)=temp_twod
#endif
   CASE (3)        
#ifdef MCKPP_CAM3
      IF (masterproc) THEN
#endif
      status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,:,1))
#ifdef MCKPP_CAM3
      temp_global_3d(:,:,1:NZP1,1)=temp_input(:,:,:,1)
      ENDIF ! End of masterproc section
      CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,1),temp_chunk(1,begchunk,1))
      DO ichnk=begchunk,endchunk
         ncol=get_ncols_p(ichnk)
         kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,2)=temp_chunk(1:ncol,ichnk,1:NZP1)
      ENDDO
#else
      CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,1),NZP1,ALL_OCEAN,2e20,temp_twod)
      kpp_3d_fields%U(:,:,2)=temp_twod
#endif
   CASE (4)
#ifdef MCKPP_CAM3
      IF (masterproc) THEN
#endif
      status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,:,1))
#ifdef MCKPP_CAM3
      temp_global_3d(:,:,1:NZP1,1)=temp_input(:,:,:,1)
      ENDIF ! End of masterproc section
      CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,1),temp_chunk(1,begchunk,1))
      DO ichnk=begchunk,endchunk
         ncol=get_ncols_p(ichnk)
         kpp_3d_fields(ichnk)%X(1:ncol,1:NZP1,1)=temp_chunk(1:ncol,ichnk,1:NZP1)
      ENDDO
#else
      CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,1),NZP1,ALL_OCEAN,2e20,temp_twod)
      kpp_3d_fields%X(:,:,1)=temp_twod
#endif
   CASE (5)
#ifdef MCKPP_CAM3
      IF (masterproc) THEN
#endif
      status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,:,1))
#ifdef MCKPP_CAM3
      temp_global_3d(:,:,1:NZP1,1)=temp_input(:,:,:,1)
      ENDIF ! End of masteproc section
      CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,1),temp_chunk(1,begchunk,1))
      DO ichnk=begchunk,endchunk
         ncol=get_ncols_p(ichnk)
         kpp_3d_fields(ichnk)%X(1:ncol,1:NZP1,2)=temp_chunk(1:ncol,ichnk,1:NZP1)
      ENDDO
#else
      CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,1),NZP1,ALL_OCEAN,2e20,temp_twod)
      kpp_3d_fields%X(:,:,2)=temp_twod
#endif
   CASE (6)
#ifdef MCKPP_CAM3
      IF (masterproc) THEN
#endif
      status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,:,1))
#ifdef MCKPP_CAM3
      temp_global_3d(:,:,1:NZP1,1)=temp_input(:,:,:,1)
      ENDIF ! End of masterproc section
      CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,1),temp_chunk(1,begchunk,1))
      DO ichnk=begchunk,endchunk
         ncol=get_ncols_p(ichnk)
         kpp_3d_fields(ichnk)%cp(1:ncol,1:NZP1)=temp_chunk(1:ncol,ichnk,1:NZP1)
      ENDDO
#else
      CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,1),NZP1,ALL_OCEAN,2e20,temp_twod)
      kpp_3d_fields%cp(:,1:NZP1)=temp_twod
#endif
   CASE (7)
#ifdef MCKPP_CAM3
      IF (masterproc) THEN
#endif
      status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,:,1))
#ifdef MCKPP_CAM3
      temp_global_3d(:,:,1:NZP1,1)=temp_input(:,:,:,1)
      ENDIF ! End of masterproc section
      CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,1),temp_chunk(1,begchunk,1))
      DO ichnk=begchunk,endchunk
         ncol=get_ncols_p(ichnk)
         kpp_3d_fields(ichnk)%rho(1:ncol,1:NZP1)=temp_chunk(1:ncol,ichnk,1:NZP1)
      ENDDO
#else
      CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,1),NZP1,ALL_OCEAN,2e20,temp_twod)
      kpp_3d_fields%rho(:,1:NZP1)=temp_twod    
#endif    
   CASE (8)
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
     status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,1,1))
#ifdef MCKPP_CAM3
     temp_global_2d(:,:)=temp_input(:,:,1,1)
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,1,PLON,temp_global_2d(:,:),temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%hmix(1:ncol)=temp_chunk(1:ncol,ichnk,1)
     ENDDO
#else
     CALL MCKPP_REFORMAT_MASK_INPUT_1D(temp_input(:,:,1,1),ALL_OCEAN,2e20,temp_twod(:,1))
     kpp_3d_fields%hmix=temp_twod(:,1)
#endif
  CASE (9)
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif        
     status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,1,1))
#ifdef MCKPP_CAM3
     temp_global_2d(:,:)=temp_input(:,:,1,1)
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,1,PLON,temp_global_2d,temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%kmix(1:ncol)=temp_chunk(1:ncol,ichnk,1)
     ENDDO
#else    
     CALL MCKPP_REFORMAT_MASK_INPUT_1D(temp_input(:,:,1,1),ALL_OCEAN,2e20,temp_twod(:,1))
     kpp_3d_fields%kmix=temp_twod(:,1)
#endif
  CASE (10)
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif        
     status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,1,1))
#ifdef MCKPP_CAM3
     temp_global_2d(:,:)=temp_input(:,:,1,1)
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,1,PLON,temp_global_2d,temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%Sref(1:ncol)=temp_chunk(1:ncol,ichnk,1)
     ENDDO
#else 
     CALL MCKPP_REFORMAT_MASK_INPUT_1D(temp_input(:,:,1,1),ALL_OCEAN,2e20,temp_twod(:,1))
     kpp_3d_fields%Sref=temp_twod(:,1)
#endif
  CASE (11)
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
     status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,1,1))
#ifdef MCKPP_CAM3
     temp_global_2d(:,:)=temp_input(:,:,1,1)
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,1,PLON,temp_global_2d,temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%SSref(1:ncol)=temp_chunk(1:ncol,ichnk,1)
     ENDDO
#else
     CALL MCKPP_REFORMAT_MASK_INPUT_1D(temp_input(:,:,1,1),ALL_OCEAN,2e20,temp_twod(:,1))
     kpp_3d_fields%SSref=temp_twod(:,1)
#endif
  CASE (12)    
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
     status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,1,1))
#ifdef MCKPP_CAM3
     temp_global_2d(:,:)=temp_input(:,:,1,1)
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,1,PLON,temp_global_2d,temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%Ssurf(1:ncol)=temp_chunk(1:ncol,ichnk,1)
     ENDDO
#else
     CALL MCKPP_REFORMAT_MASK_INPUT_1D(temp_input(:,:,1,1),ALL_OCEAN,2e20,temp_twod(:,1))
     kpp_3d_fields%Ssurf=temp_twod(:,1)
#endif
  CASE (13)
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
     status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,1,1))
#ifdef MCKPP_CAM3
     temp_global_2d(:,:)=temp_input(:,:,1,1)
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,1,PLON,temp_global_2d,temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%Tref(1:ncol)=temp_chunk(1:ncol,ichnk,1)
     ENDDO
#else
     CALL MCKPP_REFORMAT_MASK_INPUT_1D(temp_input(:,:,1,1),ALL_OCEAN,2e20,temp_twod(:,1))
     kpp_3d_fields%Tref=temp_twod(:,1)
#endif
  CASE (14)        
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
     status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,1,1))
#ifdef MCKPP_CAM3
     temp_global_2d(:,:)=temp_input(:,:,1,1)
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,1,PLON,temp_global_2d,temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%old(1:ncol)=temp_chunk(1:ncol,ichnk,1)
     ENDDO
#else
     CALL MCKPP_REFORMAT_MASK_INPUT_1D(temp_input(:,:,1,1),ALL_OCEAN,2e20,temp_twod(:,1))
     kpp_3d_fields%old=temp_twod(:,1)
#endif
  CASE (15)
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
     status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,1,1))
#ifdef MCKPP_CAM3
     temp_global_2d(:,:)=temp_input(:,:,1,1)
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,1,PLON,temp_global_2d,temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%old(1:ncol)=temp_chunk(1:ncol,ichnk,1)
     ENDDO
#else
     CALL MCKPP_REFORMAT_MASK_INPUT_1D(temp_input(:,:,1,1),ALL_OCEAN,2e20,temp_twod(:,1))
     kpp_3d_fields%new=temp_twod(:,1)
#endif
  CASE (16)
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
     status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,1:NZP1,:))
#ifdef MCKPP_CAM3
     temp_global_3d(:,:,:,:)=temp_input(:,:,:,:)
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,1),temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%Us(1:ncol,1:NZP1,1,0)=temp_chunk(1:ncol,ichnk,1:NZP1)
     ENDDO
     CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,2),temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%Us(1:ncol,1:NZP1,1,1)=temp_chunk(1:ncol,ichnk,1:NZP1)
     ENDDO
#else
     CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,1),NZP1,ALL_OCEAN,2e20,temp_twod)
     kpp_3d_fields%Us(:,:,1,0)=temp_twod
     CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,2),NZP1,ALL_OCEAN,2e20,temp_twod)
     kpp_3d_fields%Us(:,:,1,1)=temp_twod  
#endif     
  CASE (17)
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
     status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,1:NZP1,:))
#ifdef MCKPP_CAM3
     temp_global_3d(:,:,:,:)=temp_input(:,:,:,:)
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,1),temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%Us(1:ncol,1:NZP1,2,0)=temp_chunk(1:ncol,ichnk,1:NZP1)
     ENDDO
     CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,2),temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%Us(1:ncol,1:NZP1,2,1)=temp_chunk(1:ncol,ichnk,1:NZP1)
     ENDDO
#else
     CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,1),NZP1,ALL_OCEAN,2e20,temp_twod)
     kpp_3d_fields%Us(:,:,2,0)=temp_twod
     CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,2),NZP1,ALL_OCEAN,2e20,temp_twod)
     kpp_3d_fields%Us(:,:,2,1)=temp_twod
#endif
  CASE (18)
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
     status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,1:NZP1,:))
#ifdef MCKPP_CAM3
     temp_global_3d(:,:,:,:)=temp_input(:,:,:,:)
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,1),temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%Xs(1:ncol,1:NZP1,1,0)=temp_chunk(1:ncol,ichnk,1:NZP1)
     ENDDO
     CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,2),temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%Xs(1:ncol,1:NZP1,1,1)=temp_chunk(1:ncol,ichnk,1:NZP1)
     ENDDO
#else
     CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,1),NZP1,ALL_OCEAN,2e20,temp_twod)
     kpp_3d_fields%Xs(:,:,1,0)=temp_twod
     CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,2),NZP1,ALL_OCEAN,2e20,temp_twod)
     kpp_3d_fields%Xs(:,:,1,1)=temp_twod     
#endif
  CASE (19)
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
     status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,1:NZP1,:))
#ifdef MCKPP_CAM3
     temp_global_3d(:,:,:,:)=temp_input(:,:,:,:)
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,1),temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%Xs(1:ncol,1:NZP1,2,0)=temp_chunk(1:ncol,ichnk,1:NZP1)
     ENDDO
     CALL scatter_field_to_chunk(1,1,NZP1,PLON,temp_global_3d(:,:,:,2),temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%Xs(1:ncol,1:NZP1,2,1)=temp_chunk(1:ncol,ichnk,1:NZP1)
     ENDDO
#else
     CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,1),NZP1,ALL_OCEAN,2e20,temp_twod)
     kpp_3d_fields%Xs(:,:,2,0)=temp_twod
     CALL MCKPP_REFORMAT_MASK_INPUT_2D(temp_input(:,:,:,2),NZP1,ALL_OCEAN,2e20,temp_twod)
     kpp_3d_fields%Xs(:,:,2,1)=temp_twod   
#endif  
  CASE (20)
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
     status=NF_GET_VAR_REAL(ncid,varids(ivar),temp_input(:,:,1,:))
#ifdef MCKPP_CAM3
     temp_global_3d(:,:,1,:)=temp_input(:,:,1,:)
     ENDIF ! End of masterproc section
     CALL scatter_field_to_chunk(1,1,1,PLON,temp_global_3d(:,:,1,1),temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%hmixd(1:ncol,0)=temp_chunk(1:ncol,ichnk,1)
     ENDDO
     CALL scatter_field_to_chunk(1,1,1,PLON,temp_global_3d(:,:,1,2),temp_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%hmixd(1:ncol,1)=temp_chunk(1:ncol,ichnk,1)
     ENDDO
#else
     CALL MCKPP_REFORMAT_MASK_INPUT_1D(temp_input(:,:,1,1),ALL_OCEAN,2e20,temp_twod(:,1))
     kpp_3d_fields%hmixd(:,0)=temp_twod(:,1)
     CALL MCKPP_REFORMAT_MASK_INPUT_1D(temp_input(:,:,1,2),ALL_OCEAN,2e20,temp_twod(:,1))
     kpp_3d_fields%hmixd(:,1)=temp_twod(:,1)
#endif
  END SELECT
ENDDO

#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif
  status=NF_CLOSE(ncid)
  IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
#ifdef MCKPP_CAM3
  ENDIF ! End of masterproc section
#endif

END SUBROUTINE MCKPP_RESTART_IO_READ_NETCDF

