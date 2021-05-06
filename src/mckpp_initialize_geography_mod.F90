MODULE mckpp_initialize_geography_mod

  USE mckpp_netcdf_subs

CONTAINS

SUBROUTINE mckpp_initialize_geography()

#ifdef MCKPP_CAM3
  USE shr_kind_mod, only: r8=>shr_kind_r8
  USE mckpp_types, only: kpp_const_fields,kpp_3d_fields
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid,only: get_ncols_p
#else
  USE mckpp_data_fields, ONLY: kpp_const_fields, kpp_3d_fields
#endif /*MCKPP_CAM3*/
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: nz, nzp1, npts

  IMPLICIT NONE

#ifdef MCKPP_CAM3
  INTEGER :: ichnk,icol,ncol
#endif  
  
#include <netcdf.inc>
  
  ! Local Variables
  REAL sumh,hsum,dfac,sk
  REAL*4 vgrid_in(NZ)
  INTEGER i,ipt,ncid,status,dimid,varid
  CHARACTER(LEN=26) :: routine = "MCKPP_INITIALIZE_GEOGRAPHY"
  CHARACTER(LEN=max_message_len) :: message

  ! define vertical grid fields
  IF (kpp_const_fields%L_VGRID_FILE) THEN 
     ! WRITE(message,*) "Reading vertical grid from file ", kpp_const_fields%vgrid_file
     ! CALL mckpp_print(routine, message)
     status=NF_OPEN(kpp_const_fields%vgrid_file,0,ncid)
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_INQ_VARID(ncid,'d',varid)
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,vgrid_in)
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     kpp_const_fields%dm(1:NZ)=vgrid_in
     status=NF_INQ_VARID(ncid,'h',varid)
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,vgrid_in)
     kpp_const_fields%hm(1:NZ)=vgrid_in
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)     
     status=NF_INQ_VARID(ncid,'z',varid)
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_GET_VAR_REAL(ncid,varid,vgrid_in)
     kpp_const_fields%zm(1:NZ)=vgrid_in
     IF (status.ne.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     status=NF_CLOSE(ncid)
     kpp_const_fields%DMAX=-1.*(kpp_const_fields%zm(NZ)-kpp_const_fields%hm(NZ))
     ! WRITE(message,*) 'hm = ', kpp_const_fields%hm, 'zm = ', kpp_const_fields%zm,' &
     !     dm = ', kpp_const_fields%dm,' DMAX = ', kpp_const_fields%DMAX
     ! CALL mckpp_print(routine, message)
  ELSE     
     IF (kpp_const_fields%L_STRETCHGRID) THEN
        sumh = 0.0
        dfac = 1.0 - exp(-kpp_const_fields%dscale)
        DO i = 1,NZ
           sk = - (float(i)-0.5)/float(NZ)
           kpp_const_fields%hm(i) = kpp_const_fields%DMAX*dfac/&
                float(NZ)/kpp_const_fields%dscale / ( 1.0 + sk*dfac )
           sumh = sumh + kpp_const_fields%hm(i)
        ENDDO
     ENDIF
        
     ! layer thickness h, layer grids zgrid, interface depths d
     hsum = 0.0
     DO i=1,NZ
        if(kpp_const_fields%L_STRETCHGRID) then
           kpp_const_fields%hm(i) = kpp_const_fields%hm(i) * kpp_const_fields%DMAX / sumh 
        else   
           kpp_const_fields%hm(i) = kpp_const_fields%DMAX / real(NZ) 
        endif
        kpp_const_fields%zm(i) =  0.0 - (hsum + 0.5 * kpp_const_fields%hm(i) )
        hsum = hsum + kpp_const_fields%hm(i)
        kpp_const_fields%dm(i) = hsum
     ENDDO
  ENDIF
  kpp_const_fields%dm(0) = 0.0
  kpp_const_fields%hm(nzp1) = 1.e-10 
  kpp_const_fields%zm(nzp1) = -kpp_const_fields%DMAX
  
  ! Calculate Coriolis parameter
  ! Enforce minimum value of Coriolis parameter equal to 2.5 degrees latitude

#ifdef MCKPP_CAM3
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     DO icol=1,ncol
        IF (ABS(kpp_3d_fields(ichnk)%dlat(icol)) .lt. 2.5) THEN
           kpp_3d_fields(ichnk)%f(icol) = 2. * (kpp_const_fields%twopi/86164.) * &
                sin(2.5*kpp_const_fields%twopi/360.)*SIGN(1.,kpp_3d_fields(ichnk)%dlat(icol))
        ELSE
           kpp_3d_fields(ichnk)%f(icol) = 2. * (kpp_const_fields%twopi/86164.) * &
                sin(kpp_3d_fields(ichnk)%dlat(icol)*kpp_const_fields%twopi/360.)
        ENDIF
     ENDDO
  ENDDO
#else
  DO ipt=1,npts
     if(abs(kpp_3d_fields%dlat(ipt)).lt.2.5) then    
        kpp_3d_fields%f(ipt) = 2. * (kpp_const_fields%twopi/86164.) * &
             sin(2.5*kpp_const_fields%twopi/360.)*SIGN(1.,kpp_3d_fields%dlat(ipt))
     else  
        kpp_3d_fields%f(ipt) = 2. * (kpp_const_fields%twopi/86164.) * &
             sin(kpp_3d_fields%dlat(ipt)*kpp_const_fields%twopi/360.)
     endif
  ENDDO  
#endif
  
END SUBROUTINE mckpp_initialize_geography

END MODULE mckpp_initialize_geography_mod