#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE mckpp_output_mean(diag_num)
  USE mckpp_types, only: kpp_3d_fields,kpp_const_fields
  USE shr_kind_mod, only: r8=>shr_kind_r8, r4=>shr_kind_r4
  USE ppgrid, only: pcols,begchunk,endchunk
  USE phys_grid,only: get_ncols_p,gather_chunk_to_field
  USE pmgrid, only: masterproc
#else
SUBROUTINE mckpp_output_mean(kpp_3d_fields,kpp_const_fields,diag_num)
#endif

  IMPLICIT NONE
      
  INTEGER,parameter :: nuout=6,nuerr=0
#include <netcdf.inc>
  INTEGER, intent(in) :: diag_num

#ifdef MCKPP_CAM3
#include <parameter.inc>
  INTEGER, parameter :: my_nx=PLON, my_ny=PLAT
  INTEGER :: icol,ncol,ichnk
  REAL(r8) :: temp_global_3d(my_nx,my_ny,NZP1),temp_global_2d(my_nx,my_ny),&
       temp_zprof(PCOLS,begchunk:endchunk,NZP1),temp_chunk(PCOLS,begchunk:endchunk,NZP1)
  REAL(r4) :: varout(my_nx,my_ny,NZP1),singout(my_nx,my_ny)
#else
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER, parameter :: my_nx=NX, my_ny=NY
  REAL*4,allocatable :: temp_2d(:,:),temp_1d(:),temp_zprof(:,:)
  REAL*4,allocatable :: VAROUT(:,:,:), SINGOUT(:,:)
#endif
  REAL*4, parameter :: missval=1.e20
  INTEGER i,j,ivar,ipt,ix,iy,start(4),count(4),k,status,mean_num,zprof
  
  count(1:2)=(/my_nx,my_ny/)
  start(1:2)=(/1,1/)
  
  IF (diag_num .le. N_VAROUTS) THEN
     zprof=kpp_const_fields%zprof_varout_mean(diag_num)
     start(3)=1
     start(4)=kpp_const_fields%ntout_vec_mean(diag_num)
     count(3)=kpp_const_fields%zprofs_nvalid(zprof)
     count(4)=1
     
#ifndef MCKPP_CAM3
     allocate(temp_2d(NPTS,NZP1))
     allocate(VAROUT(my_nx,my_ny,kpp_const_fields%zprofs_nvalid(zprof)))
#endif
     mean_num=1
     DO i=1,diag_num-1
        IF (kpp_const_fields%ndt_varout_mean(i) .gt. 0) mean_num=mean_num+1
     ENDDO
     SELECT CASE (diag_num)
     CASE (4) 
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           DO k=1,NZP1
              temp_chunk(1:PCOLS,ichnk,k)=kpp_3d_fields(ichnk)%VEC_mean(1:ncol,k,mean_num)+&
                   kpp_3d_fields(ichnk)%Sref(:)
           ENDDO
        ENDDO
#else
        DO k=1,NZP1
           temp_2d(:,k)=kpp_3d_fields%VEC_mean(:,k,mean_num)+kpp_3d_fields%Sref(:)
        ENDDO
#endif
     CASE (12,13,14)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1)=0.
           temp_chunk(1:ncol,ichnk,2:NZP1)=kpp_3d_fields(ichnk)%VEC_mean(1:ncol,2:NZP1,mean_num)
        ENDDO
#else
        temp_2d(:,1)=0.
        temp_2d(:,2:NZP1)=kpp_3d_fields%VEC_mean(:,2:NZP1,mean_num)
#endif
     CASE (18,19)
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZ)=kpp_3d_fields(ichnk)%VEC_mean(1:ncol,1:NZ,mean_num)
           temp_chunk(:,ichnk,NZP1)=0.
        ENDDO
#else
        temp_2d(:,1:NZ)=kpp_3d_fields%VEC_mean(:,1:NZ,mean_num)
        temp_2d(:,NZP1)=0.               
#endif
     CASE DEFAULT
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%VEC_mean(1:ncol,1:NZP1,mean_num)
        ENDDO
#else
        temp_2d(:,:)=kpp_3d_fields%VEC_mean(:,:,mean_num)
#endif
     END SELECT

     IF (zprof .gt. 0) THEN
        j=1
#ifdef MCKPP_CAM3
        !allocate(temp_zprof(PCOLS,begchunk:endchunk,kpp_const_fields%zprofs_nvalid(zprof)))
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
        !allocate(temp_global_3d(my_nx,my_ny,kpp_const_fields%zprofs_nvalid(zprof)))
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
        !allocate(temp_global_3d(my_nx,my_ny,NZP1))
        CALL gather_chunk_to_field(1,1,NZP1,PLON,temp_chunk(1,begchunk,1),temp_global_3d)        
#else
        CALL MCKPP_REFORMAT_MASK_OUTPUT_2D(temp_2d,NZP1,kpp_3d_fields%L_OCEAN,missval,varout)
#endif
     ENDIF
     
#ifdef MCKPP_CAM3
     IF (masterproc) THEN
        varout(:,:,:)=temp_global_3d
        !WRITE(nuout,*) 'MCKPP_OUTPUT_MEAN: Writing output for diagnostic ',diag_num,&
        !     'to start=',start,'count=',count,'zprof=',zprof
#endif     
     status=NF_PUT_VARA_REAL(kpp_const_fields%mean_ncid_out,kpp_const_fields%varid_vec_mean(diag_num),&
          start,count,VAROUT(:,:,1:kpp_const_fields%zprofs_nvalid(zprof)))
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
#ifdef MCKPP_CAM3
     ENDIF ! End of masterproc section
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%VEC_mean(1:ncol,1:NZP1,mean_num)=0.
     ENDDO
     !deallocate(temp_global_3d)
#else     
     kpp_3d_fields%VEC_mean(:,:,mean_num)=0.
     deallocate(temp_2d)
#endif
     kpp_const_fields%ntout_vec_mean(diag_num)=kpp_const_fields%ntout_vec_mean(diag_num)+1
     i=i+1
  ELSE
#ifdef MCKPP_CAM3
!     allocate(temp_global_2d(my_nx,my_ny))
!     allocate(singout(my_nx,my_ny))
     temp_global_2d(:,:)=missval
#else     
     allocate(SINGOUT(NX,NY))
     SINGOUT(:,:) = missval
#endif     
     start(3)=kpp_const_fields%ntout_sing_mean(diag_num-N_VAROUTS)
     count(3)=1 
     mean_num=1
     DO i=1,(diag_num-N_VAROUTS-1)
        IF (kpp_const_fields%ndt_singout_mean(i) .gt. 0) mean_num=mean_num+1
     ENDDO
#ifdef MCKPP_CAM3
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        temp_chunk(1:ncol,ichnk,1)=missval
        DO j=1,ncol
           IF (kpp_3d_fields(ichnk)%L_OCEAN(j)) &
                temp_chunk(j,ichnk,1)=kpp_3d_fields(ichnk)%SCLR_mean(j,mean_num)
        ENDDO
     ENDDO
     CALL gather_chunk_to_field(1,1,1,PLON,temp_chunk(1,begchunk,1),temp_global_2d)
     IF (masterproc) THEN
        SINGOUT=temp_global_2d
#else   
     DO ix=1,nx
        DO iy=1,ny
           ipt=(iy-1)*nx+ix
           IF (kpp_3d_fields%L_OCEAN(ipt)) SINGOUT(ix,iy)=kpp_3d_fields%SCLR_mean(ipt,mean_num)
        ENDDO
     ENDDO
#endif
     !WRITE(nuout,*) 'In write_means for singout, i=',mean_num,'diag_num=',diag_num
     status=NF_PUT_VARA_REAL(kpp_const_fields%mean_ncid_out,&
          kpp_const_fields%varid_sing_mean(diag_num-N_VAROUTS),start,count,SINGOUT)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
#ifdef MCKPP_CAM3
     ENDIF ! End of masterproc section
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%SCLR_mean(1:ncol,mean_num)=0.
     ENDDO
     !deallocate(temp_global_2d)
#else
     kpp_3d_fields%SCLR_mean(:,mean_num)=0.
#endif     
     kpp_const_fields%ntout_sing_mean(diag_num-N_VAROUTS)=&
          kpp_const_fields%ntout_sing_mean(diag_num-N_VAROUTS)+1
     i=i+1         
  ENDIF
  start(3)=1

  RETURN
END SUBROUTINE mckpp_output_mean
