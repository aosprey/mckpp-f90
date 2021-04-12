SUBROUTINE mckpp_initialize_fluxes_variables()
  
#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_3d_fields
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid,only : get_ncols_p
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields
#endif  
  USE mckpp_parameters, ONLY: npts, nsflxs

! Set up parameters for calculating fluxes and initialize fluxes.
! intermediate values computed every ndtld
  IMPLICIT NONE
  
#ifdef MCKPP_CAM3
  INTEGER :: ichnk,ncol
#endif
  integer i,ipt

! extrapolation parameters for fluxes
! Initialize flux arrays
#ifdef MCKPP_CAM3
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%wU(1:ncol,:,:)=0.0
     kpp_3d_fields(ichnk)%wX(1:ncol,:,:)=0.0
     kpp_3d_fields(ichnk)%wXNT(1:ncol,:,:)=0.0
     kpp_3d_fields(ichnk)%sflux(1:ncol,:,:,:)=0.0
     kpp_3d_fields(ichnk)%sflux(:,:,5,0)=1e-20
  ENDDO
#else
  kpp_3d_fields%wU=0.0
  kpp_3d_fields%wX=0.0
  kpp_3d_fields%wXNT=0.0
  kpp_3d_fields%sflux=0.0
  DO ipt=1,npts
     do i=1,nsflxs
        kpp_3d_fields%sflux(ipt,i,5,0)=1e-20
     enddo
  ENDDO
#endif

end SUBROUTINE mckpp_initialize_fluxes_variables

! No support for data atmosphere when coupled to CAM3
#ifndef MCKPP_CAM3 
SUBROUTINE mckpp_initialize_fluxes_file()
  
  USE mckpp_data_fields, ONLY: kpp_const_fields

  IMPLICIT NONE
  
#include <netcdf.inc>
  
  INTEGER status,index(3)
  
  index(1)=1
  index(2)=1
  index(3)=1
  
  WRITE(6,*) 'MCKPP_INITIALIZE_FLUXES: Opening file ',kpp_const_fields%forcing_file
  status=NF_OPEN(kpp_const_fields%forcing_file,0,kpp_const_fields%flx_ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'time',kpp_const_fields%flx_timein_id)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'taux',kpp_const_fields%flx_varin_id(1))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'tauy',kpp_const_fields%flx_varin_id(2))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'swf',kpp_const_fields%flx_varin_id(3))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'lwf',kpp_const_fields%flx_varin_id(4))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'lhf',kpp_const_fields%flx_varin_id(5))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'shf',kpp_const_fields%flx_varin_id(6))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_INQ_VARID(kpp_const_fields%flx_ncid,'precip',kpp_const_fields%flx_varin_id(7))
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  status=NF_GET_VAR1_REAL(kpp_const_fields%flx_ncid,kpp_const_fields%flx_timein_id,&
       index,kpp_const_fields%flx_first_timein)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
            
  RETURN
END SUBROUTINE mckpp_initialize_fluxes_file
#endif
