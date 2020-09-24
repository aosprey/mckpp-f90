#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE mckpp_output_inst(diag_num)
  USE mckpp_types, only: kpp_3d_fields,kpp_const_fields
  USE shr_kind_mod, only: r8=>shr_kind_r8, r4=>shr_kind_r4
  USE ppgrid, only: pcols,begchunk,endchunk
  USE phys_grid,only: get_ncols_p,gather_chunk_to_field
  USE pmgrid, only: masterproc
#else
SUBROUTINE mckpp_output_inst(kpp_3d_fields,kpp_const_fields,diag_num)      
#endif

  IMPLICIT NONE
  INTEGER,parameter :: nuout=6,nuerr=0
#include <netcdf.inc>
  INTEGER,intent(in) :: diag_num

#ifdef MCKPP_CAM3
#include <parameter.inc>
  INTEGER, parameter :: my_nx=PLON, my_ny=PLAT
  INTEGER :: icol,ncol,ichnk
  REAL(r8) :: temp_global_3d(my_nx,my_ny,NZP1),temp_global_2d(my_nx,my_ny),&
       temp_zprof(PCOLS,begchunk:endchunk,NZP1)
  REAL(r8) :: temp_chunk(PCOLS,begchunk:endchunk,NZP1)
  REAL(r4) :: varout(my_nx,my_ny,NZP1,1),singout(my_nx,my_ny)
#else
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER, parameter :: my_nx=NX, my_ny=NY
  REAL*4, allocatable :: temp_2d(:,:),temp_1d(:),temp_zprof(:,:)
  REAL*4, allocatable :: varout(:,:,:,:), singout(:,:)
#endif
  REAL*4, parameter :: missval=1.e20
  INTEGER start(4),count(4),status
  INTEGER i,j,k,ivar
  INTEGER ix,iy,ipt,zprof
  
  count(1:2)=(/my_nx,my_ny/)
  start(1:2)=(/1,1/)
          
  IF (diag_num .le. N_VAROUTS) THEN 
     zprof=kpp_const_fields%zprof_varout_inst(diag_num)
     start(3)=1
     start(4)=kpp_const_fields%ntout_vec_inst(diag_num)
     count(3)=kpp_const_fields%zprofs_nvalid(zprof)
     count(4)=1          
     
#ifndef MCKPP_CAM3
     allocate(temp_2d(NPTS,NZP1))
     allocate(varout(my_nx,my_ny,kpp_const_fields%zprofs_nvalid(zprof),1))
#endif
     SELECT CASE (diag_num)
     CASE (1)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,1)           
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%U(:,:,1)
#endif
     CASE (2)        
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)        
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%U(1:ncol,1:NZP1,2)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%U(:,:,2)
#endif
     CASE (3)          
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%X(1:ncol,1:NZP1,1)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%X(:,:,1)
#endif
     CASE (4)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           DO k=1,NZP1
              temp_chunk(1:ncol,ichnk,k)=kpp_3d_fields(ichnk)%X(1:ncol,k,2)+&
                   kpp_3d_fields(ichnk)%Sref(1:ncol)
           ENDDO
        ENDDO
#else
        DO k=1,NZP1
           temp_2d(:,k)=kpp_3d_fields%X(:,k,2)+kpp_3d_fields%Sref(:)
        ENDDO
#endif
     CASE(5)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%buoy(1:ncol,1:NZP1)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%buoy(:,1:NZP1)
#endif
     CASE(6)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%wU(1:ncol,0:NZ,1)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%wU(:,0:NZ,1)
#endif
     CASE(7)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%wU(1:ncol,0:NZ,2)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%wU(:,0:NZ,2)
#endif
     CASE(8)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%wX(1:ncol,0:NZ,1)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%wX(:,0:NZ,1)
#endif
     CASE(9)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%wX(1:ncol,0:NZ,2)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%wX(:,0:NZ,2)
#endif
     CASE(10)        
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%wX(1:ncol,0:NZ,NSP1)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%wX(:,0:NZ,NSP1)
#endif
     CASE(11)      
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%wXNT(1:ncol,0:NZ,1)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%wXNT(:,0:NZ,1)
#endif
     CASE(12)      
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,2:NZP1)=kpp_3d_fields(ichnk)%difm(1:ncol,1:NZ)
           temp_chunk(1:ncol,ichnk,1)=0.0
        ENDDO
#else        
        temp_2d(:,1)=0.0
        temp_2d(:,2:NZP1)=kpp_3d_fields%difm(:,1:NZ)
#endif
     CASE(13)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,2:NZP1)=kpp_3d_fields(ichnk)%dift(1:ncol,1:NZ)
           temp_chunk(1:ncol,ichnk,1)=0.0
        ENDDO
#else
        temp_2d(:,1)=0.0
        temp_2d(:,2:NZP1)=kpp_3d_fields%dift(:,1:NZ)
#endif
     CASE(14)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,2:NZP1)=kpp_3d_fields(ichnk)%difs(1:ncol,1:NZ)
           temp_chunk(1:ncol,ichnk,1)=0.0
        ENDDO
#else
        temp_2d(:,1)=0.0
        temp_2d(:,2:NZP1)=kpp_3d_fields%difs(:,1:NZ)
#endif
     CASE(15)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%rho(1:ncol,1:NZP1)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%rho(:,1:NZP1)
#endif
     CASE(16)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%cp(1:ncol,1:NZP1)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%cp(:,1:NZP1)
#endif
     CASE(17)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%scorr(1:ncol,1:NZP1)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%scorr(:,:)
#endif
     CASE(18)        
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%Rig(1:ncol,1:NZP1)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%Rig(:,:)
#endif
     CASE(19)   
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZ)=kpp_3d_fields(ichnk)%dbloc(1:ncol,1:NZ)
           temp_chunk(1:ncol,ichnk,NZP1)=0.0
        ENDDO
#else
        temp_2d(:,1:NZ)=kpp_3d_fields%dbloc(:,1:NZ)
        temp_2d(:,NZP1)=0.0
#endif
     CASE(20)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%Shsq(1:ncol,1:NZP1)
        ENDDO
#else        
        temp_2d(:,:)=kpp_3d_fields%Shsq(:,:)
#endif
     CASE(21)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%Tinc_fcorr(1:ncol,1:NZP1)
        ENDDO
#else        
        temp_2d(:,:)=kpp_3d_fields%Tinc_fcorr(:,:)
#endif
     CASE(22)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%ocnTcorr(1:ncol,1:NZP1)
        ENDDO
#else        
        temp_2d(:,:)=kpp_3d_fields%ocnTcorr(:,:)
#endif
     CASE(23)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%Sinc_fcorr(1:ncol,1:NZP1)
        ENDDO
#else        
        temp_2d(:,:)=kpp_3d_fields%Sinc_fcorr(:,:)
#endif
     END SELECT
     
     IF (zprof .gt. 0) THEN
        j=1
#ifdef MCKPP_CAM3
        !allocate(temp_zprof(NCOL,begchunk:endchunk,kpp_const_fields%zprofs_nvalid(zprof)))
#else
        allocate(temp_zprof(NPTS,kpp_const_fields%zprofs_nvalid(zprof)))
#endif
        DO i=1,NZP1
           IF (kpp_const_fields%zprofs_mask(i,zprof)) THEN
#ifdef MCKPP_CAM3
              temp_zprof(:,:,j)=temp_chunk(:,:,i)
#else
              temp_zprof(:,j)=temp_2d(:,i)
#endif
              j=j+1
           ENDIF
        ENDDO
#ifdef MCKPP_CAM3
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_zprof,kpp_const_fields%zprofs_nvalid(zprof),&
             1e20,temp_zprof)
        !WRITE(6,*) my_nx,my_ny,kpp_const_fields%zprofs_nvalid(zprof)
!        IF (masterproc) &
!             allocate(temp_global_3d(my_nx,my_ny,kpp_const_fields%zprofs_nvalid(zprof)))
        CALL gather_chunk_to_field(1,1,kpp_const_fields%zprofs_nvalid(zprof),PLON,&
             temp_zprof(1,begchunk,1),temp_global_3d)
#else
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_zprof,kpp_const_fields%zprofs_nvalid(zprof),&
             kpp_3d_fields%L_OCEAN,missval,varout)
        deallocate(temp_zprof)
#endif
     ELSE       
#ifdef MCKPP_CAM3
        CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk,NZP1,1e20,temp_chunk)
        !IF (masterproc) &
        !     allocate(temp_global_3d(my_nx,my_ny,NZP1))
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d)
#else
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_2d,NZP1,kpp_3d_fields%L_OCEAN,missval,varout)
#endif
     ENDIF
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
     !varout(1:my_nx,1:my_ny,1:kpp_const_fields%zprofs_nvalid(zprof),1)=temp_global_3d
     varout(:,:,:,1)=temp_global_3d
     !WRITE(nuout,*) 'MCKPP_OUTPUT_INST: Writing output for diagnostic ',diag_num,&
     !     'to start=',start,'count=',count,'zprof=',zprof,'ncid=',kpp_const_fields%ncid_out,&
     !     'varid=',kpp_const_fields%varid_vec_inst(diag_num)
     !WRITE(nuout,*) VAROUT(:,:,38,1)
#endif
     status=NF_PUT_VARA_REAL(kpp_const_fields%ncid_out,kpp_const_fields%varid_vec_inst(diag_num),&
          start,count,VAROUT(:,:,1:kpp_const_fields%zprofs_nvalid(zprof),1))
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
#ifdef MCKPP_CAM3
     ENDIF ! End of masterproc section
     !deallocate(temp_global_3d)
#else
     deallocate(temp_2d)
     deallocate(varout)
#endif
     kpp_const_fields%ntout_vec_inst(diag_num)=kpp_const_fields%ntout_vec_inst(diag_num)+1
  ELSE
#ifdef MCKPP_CAM3
     !allocate(temp_global_2d(my_nx,my_ny))
     !allocate(singout(my_nx,my_ny))
#else
     allocate(singout(NX,NY))
     allocate(temp_1d(NPTS))
#endif
     start(3)=kpp_const_fields%ntout_sing_inst(diag_num-N_VAROUTS)
     count(3)=1
     SELECT CASE (diag_num-N_VAROUTS)
     CASE (1)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%hmix(1:ncol)
        ENDDO
#else
        temp_1d(:)=kpp_3d_fields%hmix(:)
#endif
     CASE (2)        
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%fcorr(1:ncol)
        ENDDO
#else
        temp_1d(:)=kpp_3d_fields%fcorr(:)
#endif
     CASE (3:6)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%sflux(1:ncol,diag_num-N_VAROUTS-2,5,0)
        ENDDO
#else
        temp_1d(:)=kpp_3d_fields%sflux(:,diag_num-2,5,0)
#endif
     CASE (7)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%sflux(1:ncol,6,5,0)
        ENDDO
#else
        temp_1d(:)=kpp_3d_fields%sflux(:,6,5,0)
#endif
     CASE (8)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%cplwght(1:ncol)
        ENDDO
#else
        DO ix=ifirst,ilast
           DO iy=jfirst,jlast
              ipt=(iy-1)*NX_GLOBE+ix
              temp_1d((iy-jfirst)*NX+ix-ifirst+1)=kpp_3d_fields%cplwght(ipt)                     
           ENDDO
        ENDDO
#endif
     CASE (9)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%freeze_flag(1:ncol)
        ENDDO
#else
        temp_1d(:)=kpp_3d_fields%freeze_flag(:)
#endif
     CASE (10)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%reset_flag(1:ncol)
        ENDDO
#else
        temp_1d(:)=kpp_3d_fields%reset_flag(:)
#endif
     CASE (11)        
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%dampu_flag(1:ncol)
        ENDDO
#else
        temp_1d(:)=kpp_3d_fields%dampu_flag(:)
#endif
     CASE (12)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=kpp_3d_fields(ichnk)%dampv_flag(1:ncol)
        ENDDO
#else
        temp_1d(:)=kpp_3d_fields%dampv_flag(:)
#endif
     CASE DEFAULT
        WRITE(6,*) 'MCKPP_OUTPUT_INST: Invalid diagnostic specified.'
     END SELECT

#ifdef MCKPP_CAM3
     CALL MCKPP_REFORMAT_MASK_OUTPUT_CHUNK(temp_chunk(:,:,1),1,1e20,temp_chunk(:,:,1))
     CALL gather_chunk_to_field(1,1,1,PLON,temp_chunk(1,begchunk,1),temp_global_2d)
     IF (masterproc) THEN
        SINGOUT=temp_global_2d
#else
     CALL MCKPP_REFORMAT_MASK_OUTPUT_1D(temp_1d,kpp_3d_fields%L_OCEAN,missval,singout)
#endif  
     
     status=NF_PUT_VARA_REAL(kpp_const_fields%ncid_out,kpp_const_fields%varid_sing_inst(diag_num-N_VAROUTS),&
          start,count,SINGOUT)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
#ifdef MCKPP_CAM3
     ENDIF ! End of masterproc section
#endif
     kpp_const_fields%ntout_sing_inst(diag_num-N_VAROUTS)=&
          kpp_const_fields%ntout_sing_inst(diag_num-N_VAROUTS)+1
  ENDIF

  RETURN
END SUBROUTINE mckpp_output_inst
