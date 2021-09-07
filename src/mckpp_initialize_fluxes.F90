MODULE mckpp_initialize_fluxes

#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_3d_fields
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid,only : get_ncols_p
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields
#endif  
  USE mckpp_parameters, ONLY: npts, nsflxs

  IMPLICIT NONE

CONTAINS

! Set up parameters for calculating fluxes and initialize fluxes.
! intermediate values computed every ndtld
SUBROUTINE mckpp_initialize_fluxes_variables()
    
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

END SUBROUTINE mckpp_initialize_fluxes_variables

END MODULE mckpp_initialize_fluxes
