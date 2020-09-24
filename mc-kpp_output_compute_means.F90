#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE mckpp_output_compute_means
  USE mckpp_types, only: kpp_3d_fields,kpp_const_fields
  USE phys_grid, only: get_ncols_p
  USE ppgrid, only: begchunk,endchunk,pcols
#else
SUBROUTINE mckpp_output_compute_means(kpp_3d_fields,kpp_const_fields)
#endif
  
  IMPLICIT NONE
  INTEGER, parameter :: nuout=6,nuerr=0

#ifdef MCKPP_CAM3
#include <parameter.inc>
#else
  ! Automatically includes parameter.inc
#include <mc-kpp_3d_type.com>
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
#endif

  REAL, allocatable :: field(:,:),vec(:)
  INTEGER :: i,j,k,ivar,ix,iy,ipt,ipt_globe

#ifdef MCKPP_CAM3
  INTEGER :: ichnk,icol,ncol
  allocate(field(PCOLS,NZP1))
  allocate(vec(PCOLS))
#else
  allocate(field(NPTS,NZP1))
  allocate(vec(NPTS))
#endif

#ifdef MCKPP_CAM3
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
#endif
  i=1
  DO ivar=1,N_VAROUTS
     IF (kpp_const_fields%ndt_varout_mean(ivar) .gt. 0) THEN
        !WRITE(6,*) 'Computing means for ivar = ',ivar,'i=',i
        !WRITE(6,*) 'ndt_varout_mean(ivar)=',kpp_const_fields%ndt_varout_mean(ivar)
        SELECT CASE (ivar)
        CASE(1)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,1)
#else
           field(:,:)=kpp_3d_fields%U(:,:,1)
#endif
        CASE(2)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,2)
#else
           field(:,:)=kpp_3d_fields%U(:,:,2)
#endif
        CASE(3)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%X(1:ncol,1:NZP1,1)
#else
           field(:,:)=kpp_3d_fields%X(:,:,1)
#endif
        CASE(4)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%X(1:ncol,1:NZP1,2)
#else
           field(:,:)=kpp_3d_fields%X(:,:,2)
#endif
        CASE(5)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%buoy(1:ncol,1:NZP1)
#else
           field(:,:)=kpp_3d_fields%buoy(:,1:NZP1)
#endif
        CASE(6)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%wu(:,0:NZ,1)
#else
           field(:,:)=kpp_3d_fields%wu(:,0:NZ,1)
#endif
        CASE(7)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%wu(:,0:NZ,2)
#else
           field(:,:)=kpp_3d_fields%wu(:,0:NZ,2)
#endif
        CASE(8)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%wx(:,0:NZ,1)
#else
           field(:,:)=kpp_3d_fields%wx(:,0:NZ,1)
#endif
        CASE(9)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%wx(:,0:NZ,2)
#else
           field(:,:)=kpp_3d_fields%wx(:,0:NZ,2)
#endif
        CASE(10)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%wx(:,0:NZ,NSP1)
#else
           field(:,:)=kpp_3d_fields%wx(:,0:NZ,NSP1)
#endif
        CASE(11)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%wXNT(:,0:NZ,1)
#else
           field(:,:)=kpp_3d_fields%wXNT(:,0:NZ,1)
#endif
        CASE(12)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%difm(:,1:NZP1)
#else
           field(:,:)=kpp_3d_fields%difm(:,1:NZP1)
#endif
        CASE(13)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%dift(:,1:NZP1)
#else
           field(:,:)=kpp_3d_fields%dift(:,1:NZP1)
#endif
        CASE(14)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%difs(:,1:NZP1)
#else
           field(:,:)=kpp_3d_fields%difs(:,1:NZP1)
#endif           
        CASE(15)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%rho(:,1:NZP1)
#else
           field(:,:)=kpp_3d_fields%rho(:,1:NZP1)
#endif
        CASE(16)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%cp(:,1:NZP1)
#else
           field(:,:)=kpp_3d_fields%cp(:,1:NZP1)
#endif
        CASE(17)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%scorr(:,:)
#else
           field(:,:)=kpp_3d_fields%scorr(:,:)
#endif
        CASE(18)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%Rig(:,:)
#else
           field(:,:)=kpp_3d_fields%Rig(:,:)
#endif
        CASE(19)
#ifdef MCKPP_CAM3
           field(:,1:NZ)=kpp_3d_fields(ichnk)%dbloc(:,1:NZ)
#else           
           field(:,1:NZ)=kpp_3d_fields%dbloc(:,1:NZ)
#endif 
           field(:,NZP1)=0.
        CASE(20)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%Shsq(:,:)
#else
           field(:,:)=kpp_3d_fields%Shsq(:,:)
#endif
        CASE(21)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%Tinc_fcorr(:,:)
#else
           field(:,:)=kpp_3d_fields%Tinc_fcorr(:,:)
#endif
        CASE(22)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%ocnTcorr(:,:)
#else
           field(:,:)=kpp_3d_fields%ocnTcorr(:,:)
#endif
        CASE(23)
#ifdef MCKPP_CAM3
           field(:,:)=kpp_3d_fields(ichnk)%sinc_fcorr(:,:)
#else
           field(:,:)=kpp_3d_fields%sinc_fcorr(:,:)
#endif
        END SELECT
#ifdef MCKPP_CAM3
        DO j=1,PCOLS
           IF (kpp_3d_fields(ichnk)%L_OCEAN(j)) THEN
              DO k=1,NZP1
                 kpp_3d_fields(ichnk)%VEC_mean(j,k,i)=&                      
                      field(j,k)/real(kpp_const_fields%ndt_varout_mean(ivar)) + kpp_3d_fields(ichnk)%VEC_mean(j,k,i)
#else
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,kpp_const_fields) &
!$OMP SHARED(i,field,ivar)
!$OMP DO SCHEDULE(static)
#endif
        DO j=1,NPTS
           IF (kpp_3d_fields%L_OCEAN(j)) THEN            
              DO k=1,NZP1
                 kpp_3d_fields%VEC_mean(j,k,i)=&
                      field(j,k)/real(kpp_const_fields%ndt_varout_mean(ivar)) + kpp_3d_fields%VEC_mean(j,k,i)
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif         
#endif
              ENDDO
           ENDIF
        ENDDO
        i=i+1
     ENDIF
  ENDDO
#ifdef MCKPP_CAM3
  ENDDO
  DO ichnk=begchunk,endchunk
#endif
  i=1
  DO ivar=1,N_SINGOUTS         
     !         WRITE(6,*) 'Means with ivar=',ivar,'i=',i
     IF (kpp_const_fields%ndt_singout_mean(ivar) .gt. 0) THEN
        SELECT CASE (ivar)
        CASE(1)
#ifdef MCKPP_CAM3
           vec(:)=kpp_3d_fields(ichnk)%hmix(1:PCOLS)
#else
           vec(:)=kpp_3d_fields%hmix(:)
#endif
        CASE(2)
#ifdef MCKPP_CAM3
           vec(:)=kpp_3d_fields(ichnk)%fcorr(1:PCOLS)
#else
           vec(:)=kpp_3d_fields%fcorr(:)
#endif
        CASE(3)
#ifdef MCKPP_CAM3
           vec(:)=kpp_3d_fields(ichnk)%sflux(1:PCOLS,1,5,0)
#else
           vec(:)=kpp_3d_fields%sflux(:,1,5,0)
#endif
        CASE(4)
#ifdef MCKPP_CAM3
           vec(:)=kpp_3d_fields(ichnk)%sflux(1:PCOLS,2,5,0)
#else
           vec(:)=kpp_3d_fields%sflux(:,2,5,0)
#endif
        CASE(5)
#ifdef MCKPP_CAM3
           vec(:)=kpp_3d_fields(ichnk)%sflux(1:PCOLS,3,5,0)
#else
           vec(:)=kpp_3d_fields%sflux(:,3,5,0)
#endif
        CASE(6)
#ifdef MCKPP_CAM3
           vec(:)=kpp_3d_fields(ichnk)%sflux(1:PCOLS,4,5,0)
#else
           vec(:)=kpp_3d_fields%sflux(:,4,5,0)
#endif
        CASE(7)
#ifdef MCKPP_CAM3
           vec(:)=kpp_3d_fields(ichnk)%sflux(1:PCOLS,6,5,0)
#else
           vec(:)=kpp_3d_fields%sflux(:,6,5,0)
#endif
        CASE(8)
#ifdef MCKPP_CAM3
           vec(:)=kpp_3d_fields(ichnk)%cplwght(1:PCOLS)
#else
           DO ix=kpp_const_fields%ifirst,kpp_const_fields%ilast
              DO iy=kpp_const_fields%jfirst,kpp_const_fields%jlast
                 ipt_globe=(iy-1)*NX_GLOBE+ix
                 ipt=(iy-kpp_const_fields%jfirst)*nx+(ix-kpp_const_fields%ifirst+1)
                 vec(ipt)=kpp_3d_fields%cplwght(ipt_globe)
              ENDDO
           ENDDO
#endif
        CASE(9)
#ifdef MCKPP_CAM3
           vec(:)=kpp_3d_fields(ichnk)%freeze_flag(1:PCOLS)
#else
           vec(:)=kpp_3d_fields%freeze_flag(:)
#endif
        CASE(10)
#ifdef MCKPP_CAM3
           vec(:)=kpp_3d_fields(ichnk)%reset_flag(1:PCOLS)
#else
           vec(:)=kpp_3d_fields%reset_flag(:)
#endif
        CASE(11)
#ifdef MCKPP_CAM3
           vec(:)=kpp_3d_Fields(ichnk)%dampu_flag(1:PCOLS)
#else
           vec(:)=kpp_3d_fields%dampu_flag(:)
#endif
        CASE(12)
#ifdef MCKPP_CAM3
           vec(:)=kpp_3d_fields(ichnk)%dampv_flag(1:PCOLS)
#else
           vec(:)=kpp_3d_fields%dampv_flag(:)
#endif
        END SELECT
#ifdef MCKPP_CAM3
        DO j=1,PCOLS
           IF (kpp_3d_fields(ichnk)%L_OCEAN(j)) &
                kpp_3d_fields(ichnk)%SCLR_mean(j,i)=vec(j)/real(kpp_const_fields%ndt_singout_mean(ivar))+&
                kpp_3d_fields(ichnk)%SCLR_mean(j,i)
#else
        DO j=1,NPTS
           IF (kpp_3d_fields%L_OCEAN(j)) &
                kpp_3d_fields%SCLR_mean(j,i)=vec(j)/real(kpp_const_fields%ndt_singout_mean(ivar))+&
                kpp_3d_fields%SCLR_mean(j,i)
#endif
        ENDDO
        i=i+1
     ENDIF
  ENDDO
#ifdef MCKPP_CAM3
  ENDDO
#endif
  deallocate(field)
  deallocate(vec)
  
  RETURN
END SUBROUTINE mckpp_output_compute_means
