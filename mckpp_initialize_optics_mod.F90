#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif

MODULE mckpp_initialize_optics_mod

CONTAINS

SUBROUTINE mckpp_initialize_optics()
  
#ifdef MCKPP_CAM3
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_types,only: kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
#endif  
  USE mckpp_parameters, ONLY: npts

  IMPLICIT NONE
#include <netcdf.inc>

#ifdef MCKPP_CAM3
  INTEGER :: jerlov_temp(PLON,PLAT),jerlov_chunk(PCOLS,begchunk:endchunk)
  INTEGER :: ichnk,ncol,icol
#else
  integer jerlov(npts)
#endif
  integer ipt, status, ncid_paras
  
  IF (kpp_const_fields%L_JERLOV) THEN
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
        status=NF_OPEN(kpp_const_fields%paras_file,0,ncid_paras)
        IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
        CALL MCKPP_READ_IPAR(ncid_paras,'jerlov',1,1,jerlov_temp)
        status=NF_CLOSE(ncid_paras)
        IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     ENDIF
     CALL scatter_field_to_chunk_int(1,1,1,PLON,jerlov_temp,jerlov_chunk(1,begchunk))     
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%jerlov(1:ncol)=jerlov_chunk(1:ncol,ichnk)
     ENDDO
#else     
     status=NF_OPEN(kpp_const_fields%paras_file,0,ncid_paras)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     call MCKPP_READ_IPAR(ncid_paras,'jerlov',1,1,kpp_3d_fields%jerlov)
     status=NF_CLOSE(ncid_paras)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
#endif
  ELSE
#ifdef MCKPP_CAM3
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%jerlov(1:ncol)=3
     ENDDO
#else
     DO ipt=1,npts
        kpp_3d_fields%jerlov(ipt)=3
     ENDDO
#endif
  ENDIF
  
END SUBROUTINE mckpp_initialize_optics

END MODULE mckpp_initialize_optics_mod
