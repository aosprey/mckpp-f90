#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE mckpp_initialize_landsea
  USE shr_kind_mod, only : r8=>shr_kind_r8
  USE mckpp_types, only: kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, get_ncols_p
#else
SUBROUTINE mckpp_initialize_landsea(kpp_3d_fields,kpp_const_fields)
#endif
 
  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)
 
#include <netcdf.inc>
#ifdef MCKPP_CAM3
#include <parameter.inc>
  REAL(r8) :: landsea(PLON,PLAT)
  REAL(r8) :: ocdepth(PLON,PLAT)
  REAL(r8) :: landsea_chunk(PCOLS,begchunk:endchunk),ocdepth_chunk(PCOLS,begchunk:endchunk)
  INTEGER :: ichnk, icol, ncol
#else
! Automatically includes parameter.inc
#include <mc-kpp_3d_type.com>
  TYPE (kpp_3d_type) :: kpp_3d_fields
  TYPE (kpp_const_type) :: kpp_const_fields
  REAL :: landsea(npts)
#endif

  INTEGER :: ipt,status,ncid_landsea
  
  IF (kpp_const_fields%L_LANDSEA) THEN
     WRITE(6,*) kpp_const_fields%landsea_file
     status=NF_OPEN(kpp_const_fields%landsea_file,0,ncid_landsea)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)    
     WRITE(6,*) 'Initialized landsea file',ncid_landsea
     
#ifdef MCKPP_CAM3
     IF (masterproc) THEN 
        CALL MCKPP_READ_PAR(ncid_landsea,'lsm',1,1,landsea)        
        CALL MCKPP_READ_PAR(ncid_landsea,'max_depth',1,1,ocdepth)
     ENDIF
     CALL scatter_field_to_chunk(1,1,1,PLON,ocdepth,ocdepth_chunk(1,begchunk))
     CALL scatter_field_to_chunk(1,1,1,PLON,landsea,landsea_chunk(1,begchunk))
     WRITE(6,*) ocdepth_chunk(2,begchunk),landsea_chunk(2,begchunk)
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%ocdepth(:)=ocdepth_chunk(:,ichnk)
        kpp_3d_fields(ichnk)%landfrac(:)=landsea_chunk(:,ichnk)
        DO icol=1,ncol
           IF (landsea_chunk(icol,ichnk) .EQ. 1.0) THEN
              kpp_3d_fields(ichnk)%L_OCEAN(icol)=.FALSE.
           ELSE
              kpp_3d_fields(ichnk)%L_OCEAN(icol)=.TRUE.
           ENDIF
        ENDDO
     ENDDO
#else
     call MCKPP_READ_PAR(kpp_3d_fields,ncid_landsea,'lsm',1,1,landsea)
     DO ipt=1,npts
        IF (landsea(ipt) .EQ. 1.0) THEN
           kpp_3d_fields%L_OCEAN(ipt)=.FALSE.
        ELSE
           kpp_3d_fields%L_OCEAN(ipt)=.TRUE.
        ENDIF
     ENDDO     
     call MCKPP_READ_PAR(kpp_3d_fields,ncid_landsea,'max_depth',1,1,kpp_3d_fields%ocdepth)
#endif     
     WRITE(6,*) 'Read landsea mask'     
     status=NF_CLOSE(ncid_landsea)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ELSE
     DO ipt=1,npts
        kpp_3d_fields%L_OCEAN(ipt)=.FALSE.
     ENDDO
  ENDIF
  
  RETURN
END SUBROUTINE mckpp_initialize_landsea
