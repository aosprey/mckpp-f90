#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE mckpp_output_range(diag_num)
  USE mckpp_types, only: kpp_3d_fields,kpp_const_fields
  USE shr_kind_mod, only: r8=>shr_kind_r8, r4=>shr_kind_r4
  USE ppgrid, only: pcols,begchunk,endchunk
  USE phys_grid, only: get_ncols_p,gather_chunk_to_field
  USE pmgrid, only: masterproc
#else
SUBROUTINE mckpp_output_range(kpp_3d_fields,kpp_const_fields,diag_num)
#endif

  IMPLICIT NONE

  INTEGER,parameter :: nuout=6,nuerr=0
  INTEGER,intent(in) :: diag_num
#include <netcdf.inc>      

#ifdef MCKPP_CAM3
#include <parameter.inc>
  INTEGER, parameter :: my_nx=PLON, my_ny=PLAT
  INTEGER :: icol,ncol,ichnk
  REAL(r8) :: temp_global_3d(my_nx,my_ny,NZP1),temp_global_2d(my_nx,my_ny),&
       temp_zprof(PCOLS,begchunk:endchunk,NZP1),temp_chunk(PCOLS,begchunk:endchunk,NZP1)
  REAL(r4) :: varout(my_nx,my_ny,NZP1),SINGOUT(my_nx,my_ny)
#else
! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER, parameter :: my_nx=NX, my_ny=NY
  REAL*4, allocatable :: temp_2d(:,:),temp_zprof(:,:),temp_1d(:)
  REAL*4,allocatable :: VAROUT(:,:,:), SINGOUT(:,:)
#endif
  REAL*4, parameter :: missval=1.e20
  INTEGER i,j,ivar,ipt,ix,iy,start(4),count(4),k,status,range_num,zprof
      
  count(1:2)=(/my_nx,my_ny/)
  start(1:2)=(/1,1/)

#ifndef MCKPP_CAM3       
  allocate(VAROUT(my_nx,my_ny,NZP1))
  allocate(SINGOUT(my_nx,my_ny))
#endif

  DO j=1,2
     IF (diag_num .le. N_VAROUTS) THEN 
        zprof=kpp_const_fields%zprof_varout_range(diag_num)
        start(3)=1
        start(4)=kpp_const_fields%ntout_vec_range(diag_num)
        count(3)=kpp_const_fields%zprofs_nvalid(zprof)
        count(4)=1
#ifndef MCKPP_CAM3
        allocate(temp_2d(NPTS,NZP1))
#endif
        range_num=1
        DO i=1,diag_num-1
           IF (kpp_const_fields%ndt_varout_range(i) .gt. 0) range_num=range_num+1
        ENDDO
        SELECT CASE (diag_num)
        CASE (4) 
#ifdef MCKPP_CAM3
           DO ichnk=begchunk,endchunk
              ncol=get_ncols_p(ichnk)
              DO k=1,NZP1
                 temp_chunk(1:PCOLS,ichnk,k)=kpp_3d_fields(ichnk)%VEC_range(1:ncol,k,range_num,j)+&
                      kpp_3d_fields(ichnk)%Sref(1:ncol)
              ENDDO
           ENDDO
#else
           DO k=1,NZP1
              temp_2d(:,k)=kpp_3d_fields%VEC_range(:,k,range_num,j)+kpp_3d_fields%Sref(:)
           ENDDO
#endif
        CASE (12,13,14)
#ifdef MCKPP_CAM3
           DO ichnk=begchunk,endchunk
              ncol=get_ncols_p(ichnk)
              temp_chunk(1:ncol,ichnk,1)=0.
              temp_chunk(1:ncol,ichnk,2:NZP1)=kpp_3d_fields(ichnk)%VEC_range(1:ncol,2:NZP1,range_num,j)
           ENDDO
#else
           temp_2d(:,1)=0.
           temp_2d(:,2:NZP1)=kpp_3d_fields%VEC_range(:,2:NZP1,range_num,j)
#endif
        CASE (18,19)
#ifdef MCKPP_CAM3
           DO ichnk=begchunk,endchunk
              ncol=get_ncols_p(ichnk)
              temp_chunk(1:ncol,ichnk,1:NZ)=kpp_3d_fields(ichnk)%VEC_range(1:ncol,1:NZ,range_num,j)
              temp_chunk(1:ncol,ichnk,NZP1)=0.
           ENDDO
#else
           temp_2d(:,1:NZ)=kpp_3d_fields%VEC_range(:,1:NZ,range_num,j)
           temp_2d(:,NZP1)=0.               
#endif
        CASE DEFAULT
#ifdef MCKPP_CAM3
           DO ichnk=begchunk,endchunk
              ncol=get_ncols_p(ichnk)
              temp_chunk(1:ncol,ichnk,1:NZP1)=kpp_3d_fields(ichnk)%VEC_range(1:ncol,1:NZP1,range_num,j)
           ENDDO
#else
           temp_2d(:,:)=kpp_3d_fields%VEC_range(:,:,range_num,j)
#endif
        END SELECT
        
        IF (zprof .gt. 0) THEN
           k=1
#ifdef MCKPP_CAM3
           !IF (j .eq. 1) &
           !     allocate(temp_zprof(PCOLS,begchunk:endchunk,kpp_const_fields%zprofs_nvalid(zprof)))
#else
           allocate(temp_zprof(NPTS,kpp_const_fields%zprofs_nvalid(zprof)))
#endif
           DO i=1,NZP1
              IF (kpp_const_fields%zprofs_mask(i,zprof)) THEN
#ifdef MCKPP_CAM3
                 temp_zprof(:,:,k)=temp_chunk(:,:,i)
#else
                 temp_zprof(:,k)=temp_2d(:,i)
#endif
                 k=k+1
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
           varout=temp_global_3d
#endif
        SELECT CASE (j)
        CASE (1)                              
           !WRITE(6,*) 'MCKPP_OUTPUT_RANGE: Writing to ncid=',kpp_const_fields%min_ncid_out,&
           !     ' varid=',kpp_const_fields%varid_vec_range(diag_num),&
           !     ' with start =',start,' and count =',count,' diag_num = ',diag_num,' range_num = ',range_num,&
           !     ' nz for this zprof = ',kpp_const_fields%zprofs_nvalid(zprof)
           status=NF_PUT_VARA_REAL(kpp_const_fields%min_ncid_out,kpp_const_fields%varid_vec_range(diag_num),&
                start,count,VAROUT(:,:,1:kpp_const_fields%zprofs_nvalid(zprof)))
           IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
           !WRITE(6,*) 'Written successfully'
        CASE (2)
           !WRITE(6,*) 'MCKPP_OUTPUT_RANGE: Writing to ncid=',kpp_const_fields%max_ncid_out,&
           !     ' varid=',kpp_const_fields%varid_vec_range(diag_num),&
           !     ' with start =',start,' and count =',count,' diag_num = ',diag_num,' range_num = ',range_num,&
           !     ' nz for this zprof = ',kpp_const_fields%zprofs_nvalid(zprof)
           status=NF_PUT_VARA_REAL(kpp_const_fields%max_ncid_out,kpp_const_fields%varid_vec_range(diag_num),&
                start,count,VAROUT(:,:,1:kpp_const_fields%zprofs_nvalid(zprof)))
           IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)       
           !WRITE(6,*) 'Written successfully'
        END SELECT
#ifdef MCKPP_CAM3
        ENDIF ! End of masterproc section
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           IF (j .eq. 1) THEN
              kpp_3d_fields(ichnk)%VEC_range(1:ncol,1:NZP1,range_num,j)=2e20
           ELSEIF (j .eq. 2) THEN
              kpp_3d_fields(ichnk)%VEC_range(1:ncol,1:NZP1,range_num,j)=-2e20
           ENDIF
        ENDDO
        IF (j .eq. 2) &
           kpp_const_fields%ntout_vec_range(diag_num)=kpp_const_fields%ntout_vec_range(diag_num)+1 
!        deallocate(temp_global_3d)
#else
        IF (j .eq. 1) THEN
           kpp_3d_fields%VEC_range(:,:,range_num,j)=2e20
        ELSEIF (j .eq. 2) THEN
           kpp_3d_fields%VEC_range(:,:,range_num,j)=-2e20
           kpp_const_fields%ntout_vec_range(diag_num)=kpp_const_fields%ntout_vec_range(diag_num)+1
        ENDIF
        deallocate(temp_2d)
#endif
     ELSE

#ifdef MCKPP_CAM3
!        allocate(temp_global_2d(my_nx,my_ny))
!	allocate(singout(my_nx,my_ny))
        temp_chunk(:,:,:)=missval
#else
        singout(:,:)=missval
#endif
        start(3)=kpp_const_fields%ntout_sing_range(diag_num-N_VAROUTS)
        count(3)=1
        range_num=1
        DO i=1,(diag_num-N_VAROUTS-1)
           IF (kpp_const_fields%ndt_singout_range(i) .gt. 0) range_num=range_num+1
        ENDDO
#ifdef MCKPP_CAM3
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           DO k=1,ncol
              IF (kpp_3d_fields(ichnk)%L_OCEAN(k)) &
                   temp_chunk(k,ichnk,1)=kpp_3d_fields(ichnk)%SCLR_range(k,range_num,j)
           ENDDO
        ENDDO
        !WRITE(6,*) temp_chunk(:,:,1)
        CALL gather_chunk_to_field(1,1,1,PLON,temp_chunk(1,begchunk,1),temp_global_2d)
        IF (masterproc) THEN
           SINGOUT=temp_global_2d
#else
        DO ix=1,nx
           DO iy=1,ny
              ipt=(iy-1)*nx+ix
              IF (kpp_3d_fields%L_OCEAN(ipt)) SINGOUT(ix,iy)=kpp_3d_fields%SCLR_range(ipt,range_num,j)
           ENDDO
        ENDDO
#endif
        SELECT CASE (j)
        CASE (1)               
           status=NF_PUT_VARA_REAL(kpp_const_fields%min_ncid_out,&
                kpp_const_fields%varid_sing_range(diag_num-N_VAROUTS),start,count,SINGOUT)
           IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
        CASE (2)
           status=NF_PUT_VARA_REAL(kpp_const_fields%max_ncid_out,&
                kpp_const_fields%varid_sing_range(diag_num-N_VAROUTS),start,count,SINGOUT)
           IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
        END SELECT
#ifdef MCKPP_CAM3
        ENDIF ! End of masterproc section
        DO ichnk=begchunk,endchunk
           ncol=get_ncols_p(ichnk)
           IF (j .eq. 1) THEN
              kpp_3d_fields(ichnk)%SCLR_range(1:ncol,range_num,j)=2e20
           ELSEIF (j .eq. 2) THEN
              kpp_3d_fields(ichnk)%SCLR_range(1:ncol,range_num,j)=-2e20
           ENDIF
        ENDDO
!        deallocate(temp_global_2d)
#else
        IF (j .eq. 1) THEN
           kpp_3d_fields%SCLR_range(:,range_num,j)=2e20
        ELSEIF (j .eq. 2) THEN
           kpp_3d_fields%SCLR_range(:,range_num,j)=-2e20
        ENDIF
#endif
        IF (j .eq. 2) THEN
           kpp_const_fields%ntout_sing_range(diag_num-N_VAROUTS)=&
                kpp_const_fields%ntout_sing_range(diag_num-N_VAROUTS)+1
        ENDIF
     ENDIF
  ENDDO
  start(3)=1
  
  RETURN
END SUBROUTINE mckpp_output_range
