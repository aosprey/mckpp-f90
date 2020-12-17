#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE MCKPP_INITIALIZE_RELAXATION
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_parameters
  USE mckpp_types,only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk, endchunk, pcols
  USE phys_grid, only: scatter_field_to_chunk, scatter_field_to_chunk_int, get_ncols_p
#else
SUBROUTINE mckpp_initialize_relaxation(kpp_3d_fields,kpp_const_fields)
  USE mckpp_data_fields
#endif

! Re-write logic to allow for relaxing either SST or
! salinity - NPK 24/08/11

  IMPLICIT NONE
  INTEGER nuout,nuerr
  PARAMETER (nuout=6,nuerr=0)

#ifdef MCKPP_CAM3
  REAL(r8) :: relax_chunk(PCOLS,begchunk:endchunk)
  INTEGER :: ichnk,icol,ncol
#else
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
#endif

!#include <constants.com>
!#include <ocn_advec.com>
!#include <couple.com>
!#include <sstclim.com>
!#include <relax_3d.com>
  
  INTEGER ix,iy,ipoint,my_ny
  
!  REAL sst_in(NX_GLOBE,NY_GLOBE,1)
!  COMMON /save_sstin/ sst_in

#ifdef MCKPP_CAM3  
  my_ny=NY_GLOBE
#else
  my_ny=NY
#endif

  IF (kpp_const_fields%L_RELAX_SST) THEN 
     DO iy=1,my_ny
        IF (kpp_const_fields%relax_sst_in(iy) .EQ. 0.0) THEN
#ifdef MCKPP_CAM3
           IF (masterproc) &
                kpp_global_fields%relax(:,iy)=0.0
#else
           DO ix=1,nx
              ipoint=(iy-1)*nx+ix
              kpp_3d_fields%relax_sst(ipoint)=0.0
           ENDDO
#endif
        ELSE
#ifdef MCKPP_CAM3
           IF (masterproc) &
                kpp_global_fields%relax(:,iy)=1./(kpp_const_fields%relax_sst_in(iy)*kpp_const_fields%spd)
#else
           DO ix=1,nx
              ipoint=(iy-1)*nx+ix
              kpp_3d_fields%relax_sst(ipoint)=1./(kpp_const_fields%relax_sst_in(iy)*kpp_const_fields%spd)
           ENDDO
#endif        
        ENDIF
     ENDDO
#ifdef MCKPP_CAM3
     CALL scatter_field_to_chunk(1,1,1,PLON,kpp_global_fields%relax,relax_chunk(1,begchunk))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%relax_sst(1:ncol)=relax_chunk(1:ncol,ichnk)
     ENDDO
#endif
  ENDIF
  IF (kpp_const_fields%L_RELAX_SAL) THEN
     DO iy=1,my_ny
        IF (kpp_const_fields%relax_sal_in(iy) .EQ. 0.0) THEN
#ifdef MCKPP_CAM3
           IF (masterproc) &
                kpp_global_fields%relax(:,iy)=0.0
#else           
           DO ix=1,nx
              ipoint=(iy-1)*nx+ix
              kpp_3d_fields%relax_sal(ipoint)=0.0
           ENDDO
#endif
        ELSE
#ifdef MCKPP_CAM3
           IF (masterproc) THEN
              WRITE(6,*) iy,kpp_const_fields%relax_sal_in(iy)
              kpp_global_fields%relax(:,iy)=1./(kpp_const_fields%relax_sal_in(iy)*kpp_const_fields%spd)           
           ENDIF
#else
           DO ix=1,nx
              ipoint=(iy-1)*nx+ix
              kpp_3d_fields%relax_sal(ipoint)=1./(kpp_const_fields%relax_sal_in(iy)*kpp_const_fields%spd)
           ENDDO
#endif
        ENDIF
     ENDDO
#ifdef MCKPP_CAM3
     CALL scatter_field_to_chunk(1,1,1,PLON,kpp_global_fields%relax,relax_chunk(1,begchunk))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%relax_sal(1:ncol)=relax_chunk(1:ncol,ichnk)
     ENDDO
#endif
  ENDIF
  IF (kpp_const_fields%L_RELAX_OCNT) THEN
     DO iy=1,my_ny
        IF (kpp_const_fields%relax_ocnt_in(iy) .EQ. 0.0) THEN
#ifdef MCKPP_CAM3
           IF (masterproc) &
                kpp_global_fields%relax(:,iy)=0.0
#else           
           DO ix=1,nx
              ipoint=(iy-1)*nx+ix
              kpp_3d_fields%relax_ocnT(ipoint)=0.0
           ENDDO
#endif
        ELSE
#ifdef MCKPP_CAM3
           IF (masterproc) &
                kpp_global_fields%relax(:,iy)=1./(kpp_const_fields%relax_ocnT_in(iy)*kpp_const_fields%spd)
#else           
           DO ix=1,nx
              ipoint=(iy-1)*nx+ix
              kpp_3d_fields%relax_ocnT(ipoint)=1./(kpp_const_fields%relax_ocnT_in(iy)*kpp_const_fields%spd)
           ENDDO
#endif
        ENDIF
     ENDDO
#ifdef MCKPP_CAM3
     CALL scatter_field_to_chunk(1,1,1,PLON,kpp_global_fields%relax,relax_chunk(1,begchunk))
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%relax_ocnT(1:ncol)=relax_chunk(1:ncol,ichnk)
     ENDDO
#endif
  ENDIF
  CALL MCKPP_PHYSICS_OVERRIDES_SST0(kpp_3d_fields,kpp_const_fields)

! Do we need these initialization statements?
!  DO iy=1,ny
!     DO ix=1,nx
!        ipoint=(iy-1)*nx+ix
!        kpp_3d_fields%fcorr(ipoint)=0.0
!        kpp_3d_fields%scorr(ipoint,:)=0.0
!     ENDDO
!  ENDDO
  
  write(6,*) 'MCKPP_INITIALIZE_RELAXATION: Calculated SST0, fcorr and scorr'
  RETURN
END SUBROUTINE mckpp_initialize_relaxation
