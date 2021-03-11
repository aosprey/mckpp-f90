#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE MCKPP_INITIALIZE_ADVECTION
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_parameters
  USE mckpp_types,only: kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
#else
SUBROUTINE mckpp_initialize_advection(kpp_3d_fields,kpp_const_fields)
  USE mckpp_data_types
#endif
  
  IMPLICIT NONE
  
#include <netcdf.inc>

#ifdef MCKPP_CAM3
  REAL(r8) :: advection_chunk(PCOLS,begchunk:endchunk,2)
  INTEGER :: nmodeadv_temp(PLON,PLAT,2),modeadv_temp(PLON,PLAT,maxmodeadv,2),&
       nmodeadv_chunk(PCOLS,begchunk:endchunk,2)
  REAL(r8) :: advection_temp(PLON,PLAT,maxmodeadv,2)
  INTEGER :: ichnk,ncol,icol
#else
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
#endif

  INTEGER nmode(npts)  
  INTEGER i,ipt,ivar,imode,status,ncid_advec
  
  IF (kpp_const_fields%L_ADVECT) THEN
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
#endif
        status=NF_OPEN(kpp_const_fields%advect_file,0,ncid_advec)
        IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
#ifdef MCKPP_CAM3       
        CALL MCKPP_READ_IPAR(ncid_advec,'nmode_tadv',1,1,nmodeadv_temp(:,:,1))
        CALL MCKPP_READ_IPAR(ncid_advec,'nmode_sadv',1,1,nmodeadv_temp(:,:,2))
     ENDIF
     CALL scatter_field_to_chunk_int(1,1,2,PLON,nmodeadv_temp,nmodeadv_chunk(1,begchunk,1))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%nmodeadv(1:ncol,:)=nmodeadv_chunk(1:ncol,ichnk,:)
     ENDDO
     IF (masterproc) THEN
        CALL MCKPP_READ_IPAR(ncid_advec,'mode_tadv',maxmodeadv,1,modeadv_temp(:,:,:,1))
        CALL MCKPP_READ_IPAR(ncid_advec,'mode_sadv',maxmodeadv,1,modeadv_temp(:,:,:,2))
     ENDIF
     DO i=1,maxmodeadv 
        CALL scatter_field_to_chunk_int(1,1,2,PLON,modeadv_temp(:,:,i,:),nmodeadv_chunk(1,begchunk,1))
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           kpp_3d_fields(ichnk)%modeadv(1:ncol,i,:)=nmodeadv_chunk(1:ncol,ichnk,:)
        ENDDO
     ENDDO
     IF (masterproc) THEN        
        CALL MCKPP_READ_PAR(ncid_advec,'tadv',maxmodeadv,1,advection_temp(:,:,:,1))
        CALL MCKPP_READ_PAR(ncid_advec,'sadv',maxmodeadv,1,advection_temp(:,:,:,2))
     ENDIF
     DO i=1,maxmodeadv
        CALL scatter_field_to_chunk(1,1,2,PLON,advection_temp(:,:,i,:),advection_chunk(1,begchunk,1))
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           kpp_3d_fields(ichnk)%advection(1:ncol,i,:)=advection_chunk(1:ncol,ichnk,:)
        ENDDO
     ENDDO
#else
     call MCKPP_READ_IPAR(kpp_3d_fields,ncid_advec,'nmode_tadv',1,1,kpp_3d_fields%nmodeadv(:,1))
     call MCKPP_READ_IPAR(kpp_3d_fields,ncid_advec,'mode_tadv',maxmodeadv,1,kpp_3d_fields%modeadv(:,:,1))
     call MCKPP_READ_PAR(kpp_3d_fields,ncid_advec,'tadv',maxmodeadv,1,kpp_3d_fields%advection(:,:,1))
     call MCKPP_READ_IPAR(kpp_3d_fields,ncid_advec,'nmode_sadv',1,1,kpp_3d_fields%nmodeadv(:,2))
     call MCKPP_READ_IPAR(kpp_3d_fields,ncid_advec,'mode_sadv',maxmodeadv,1,kpp_3d_fields%modeadv(:,:,2))
     call MCKPP_READ_PAR(kpp_3d_fields,ncid_advec,'sadv',maxmodeadv,1,kpp_3d_fields%advection(:,:,2))  
#endif
     status=NF_CLOSE(ncid_advec)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ELSE
#ifdef MCKPP_CAM3
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%nmodeadv(1:ncol,:)=0
     ENDDO
#else
     DO ipt=1,npts
        DO ivar=1,2
           kpp_3d_fields%nmodeadv(ipt,ivar)=0
        ENDDO
     ENDDO
#endif
  ENDIF

  RETURN
END SUBROUTINE mckpp_initialize_advection
