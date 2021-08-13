MODULE mckpp_reformat_mask_output

CONTAINS 

#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE mckpp_reformat_mask_output_chunk(threed_in,nz_in,missval,threed_out)
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid, only: get_ncols_p
  USE shr_kind_mod, only: r8=>shr_kind_r8
  USE mckpp_types, only: kpp_3d_fields
  
  IMPLICIT NONE
  
  INTEGER, intent(in) :: nz_in
  REAL(r8), intent(in) :: threed_in(PCOLS,begchunk:endchunk,nz_in)
  REAL(r8), intent(out) :: threed_out(PCOLS,begchunk:endchunk,nz_in)
  REAL(r8), intent(in) :: missval

  INTEGER :: icol,ichnk,ncol

  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     DO icol=1,ncol
        IF (kpp_3d_fields(ichnk)%L_OCEAN(icol)) THEN
           threed_out(icol,ichnk,:)=threed_in(icol,ichnk,:)
        ELSE
           threed_out(icol,ichnk,:)=missval
        ENDIF
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE mckpp_reformat_mask_output_chunk

#else

SUBROUTINE mckpp_reformat_mask_input_2d(threed_in,nz_in,mask,missval,twod_out)

  USE mckpp_parameters, ONLY: nx, ny, npts

  IMPLICIT NONE

  INTEGER, intent(in) :: nz_in
  REAL*4,intent(in) :: threed_in(NX,NY,nz_in)
  REAL*8 :: missval
  REAL*4,intent(out) :: twod_out(NPTS,nz_in)
  LOGICAL,intent(in) :: mask(NPTS)
  INTEGER :: ix,iy,ipt

  DO ix=1,NX
     DO iy=1,NY
        ipt=(iy-1)*NX+ix
        IF (mask(ipt)) THEN
           twod_out(ipt,:)=threed_in(ix,iy,:)
        ELSE
           twod_out(ipt,:)=missval
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE mckpp_reformat_mask_input_2d


SUBROUTINE mckpp_reformat_mask_input_1d(twod_in,mask,missval,oned_out)

  USE mckpp_parameters, ONLY: nx, ny, npts

  IMPLICIT NONE

  REAL*4,intent(in) :: twod_in(NX,NY)
  REAL*8 :: missval
  REAL*4,intent(out) :: oned_out(NPTS)
  LOGICAL,intent(in) :: mask(NPTS)
  INTEGER :: ix,iy,ipt

  DO ix=1,NX
     DO iy=1,NY
        ipt=(iy-1)*NX+ix
        IF (mask(ipt)) THEN
           oned_out(ipt)=twod_in(ix,iy)
        ELSE
           oned_out(ipt)=missval
        ENDIF
     ENDDO
  ENDDO

  RETURN
END SUBROUTINE mckpp_reformat_mask_input_1d


SUBROUTINE mckpp_reformat_mask_output_1d(oned_in,mask,missval,twod_out)

  USE mckpp_parameters, ONLY: nx, ny, npts

  IMPLICIT NONE

  REAL*4,intent(in) :: oned_in(NPTS)
  REAL*8, INTENT(IN) :: missval
  REAL*4,intent(out) :: twod_out(NX,NY)
  LOGICAL,intent(in) :: mask(NPTS)
  INTEGER :: i,j,ipt
  
  DO i=1,NX
     DO j=1,NY
        ipt=(j-1)*NX+i
        IF (mask(ipt)) THEN
           twod_out(i,j)=oned_in(ipt)
        ELSE
           twod_out(i,j)=missval
        ENDIF
     ENDDO
  ENDDO

END SUBROUTINE mckpp_reformat_mask_output_1d


SUBROUTINE mckpp_reformat_mask_output_2d(twod_in,nz_in,mask,missval,threed_out)

  USE mckpp_parameters, ONLY: nx, ny, npts
  
  IMPLICIT NONE
  
  INTEGER,intent(in) :: nz_in
  REAL*4,intent(in) :: twod_in(NPTS,nz_in)
  REAL*8, INTENT(IN) :: missval
  REAL*4,intent(out) :: threed_out(NX,NY,nz_in)
  LOGICAL,intent(in) :: mask(NPTS)
  INTEGER :: i,j,ipt
  
  DO i=1,NX
     DO j=1,NY
        ipt=(j-1)*NX+i
        IF (mask(ipt)) THEN
           threed_out(i,j,:)=twod_in(ipt,:)
        ELSE
           threed_out(i,j,:)=missval
        ENDIF
     ENDDO
  ENDDO
  
END SUBROUTINE mckpp_reformat_mask_output_2d
#endif

END MODULE mckpp_reformat_mask_output
