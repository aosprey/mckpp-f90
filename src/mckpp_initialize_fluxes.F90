MODULE mckpp_initialize_fluxes

  USE mckpp_data_fields, ONLY: kpp_3d_fields
  USE mckpp_parameters, ONLY: npts, nsflxs

  IMPLICIT NONE

CONTAINS

! Set up parameters for calculating fluxes and initialize fluxes.
! intermediate values computed every ndtld
SUBROUTINE mckpp_initialize_fluxes_variables()
    
  integer i,ipt

! extrapolation parameters for fluxes
! Initialize flux arrays
  kpp_3d_fields%wU=0.0
  kpp_3d_fields%wX=0.0
  kpp_3d_fields%wXNT=0.0
  kpp_3d_fields%sflux=0.0
  DO ipt=1,npts
     do i=1,nsflxs
        kpp_3d_fields%sflux(ipt,i,5,0)=1e-20
     enddo
  ENDDO

END SUBROUTINE mckpp_initialize_fluxes_variables

END MODULE mckpp_initialize_fluxes
