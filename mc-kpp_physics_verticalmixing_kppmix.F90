SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_KPPMIX(km,kmp1,dVsq,ustar,Bo,Bosol,alphaDT,betaDS,&
     Ritop, hbl , kbl, kpp_1d_fields,kpp_const_fields)
#ifdef MCKPP_CAM3
  USE mckpp_parameters
  USE mckpp_types, only: kpp_1d_type,kpp_const_type
#else 
  USE mckpp_data_types
#endif
  IMPLICIT NONE

!.......................................................................
!
!     Main driver subroutine for kpp vertical mixing scheme and 
!     interface to greater ocean model
!
!     written by : bill large,   june,  6, 1994
!     modified by: jan morzel,   june, 30, 1994
!                  bill large, august, 11, 1994
!                  bill large, november 1994,   for 1d code
!
!.......................................................................
!     

  ! input
  integer, intent(in) :: km,kmp1
  real, intent(in) :: &
       ustar, &             ! surface friction velocity       (m/s)
       Bo, &                ! surface turbulent buoy. forcing (m^2/s^3)
       Bosol, &             ! radiative buoyancy forcing      (m^2/s^3) 
       dVsq(kmp1), &        ! (velocity shear re sfc)^2      (m/s)^2
       alphaDT(kmp1), &     ! alpha * DT  across interfaces
       betaDS(kmp1), &      ! beta  * DS  across interfaces
       Ritop(km)             ! numerator of bulk Richardson Number (m/s)^2
  type(kpp_1d_type), intent(inout) :: kpp_1d_fields
  type(kpp_const_type), intent(in) :: kpp_const_fields

  ! output
  ! visc replaced by kpp_1d_fields%difm (NPK 8/2/2013)
  real, intent(out) :: &
       hbl, &               ! boundary layer depth (m)
       kbl                  ! index of first grid level below hbl     
  
  ! local
  integer mdiff,ki
  parameter (mdiff = 3)     ! number of diffusivities for local arrays
  real bfsfc                ! surface buoyancy forcing        (m^2/s^3)
  real ws                   ! momentum velocity scale
  real wm                   ! scalar   velocity scale 
  real caseA                ! = 1 in case A; =0 in case B
  real stable               ! = 1 in stable forcing; =0 in unstable
  real dkm1(mdiff)          ! boundary layer difs at kbl-1 level
  real gat1(mdiff)          ! shape function at sigma=1
  real dat1(mdiff)          ! derivative of shape function at sigma=1
  real, allocatable :: & 
       blmc(:,:)            ! boundary layer mixing coefficients
  real sigma                ! normalized depth (d / hbl)
  real Rib(2)               ! bulk Richardson number

  ALLOCATE( blmc(km,mdiff) )
      
  ! zero the mixing coefficients 
  DO ki=0,km
     kpp_1d_fields%difm(ki) = 0.0
     kpp_1d_fields%difs(ki) = 0.0
     kpp_1d_fields%dift(ki) = 0.0
  END DO
      
  ! compute RI and IW interior diffusivities everywhere
  IF(kpp_const_fields%LRI) THEN
     call MCKPP_PHYSICS_VERTICALMIXING_RIMIX(km,kmp1,kpp_1d_fields,kpp_const_fields)
  ENDIF
     
  ! add double diffusion if desired
  IF(kpp_const_fields%LDD) THEN
     call MCKPP_PHYSICS_VERTICALMIXING_DDMIX(km,kmp1,alphaDT,betaDS,kpp_1d_fields)
  ENDIF

  ! fill the bottom kmp1 coefficients for blmix      
  kpp_1d_fields%difm(kmp1) = kpp_1d_fields%difm(km)
  kpp_1d_fields%difs(kmp1) = kpp_1d_fields%difs(km)
  kpp_1d_fields%dift(kmp1) = kpp_1d_fields%dift(km)
  
  ! compute boundary layer mixing coefficients ??
  IF(kpp_const_fields%LKPP) THEN
     ! diagnose the new boundary layer depth
     call  MCKPP_PHYSICS_VERTICALMIXING_BLDEPTH (km, kmp1, dVsq, Ritop, ustar, Bo, Bosol, & 
          hbl, bfsfc, stable, caseA, kbl, Rib, sigma, wm, &
          ws, kpp_1d_fields, kpp_const_fields)
     
     ! compute boundary layer diffusivities
     call MCKPP_PHYSICS_VERTICALMIXING_BLMIX(km, mdiff, ustar, bfsfc, hbl, stable, caseA, &
          kbl, gat1 , dat1 , dkm1, blmc, sigma, wm, ws, &
          kpp_1d_fields,kpp_const_fields)
                  
     ! enhance diffusivity at interface kbl - 1
     call MCKPP_PHYSICS_VERTICALMIXING_ENHANCE (km, mdiff, dkm1, hbl, kbl, caseA, &
          blmc, kpp_1d_fields,kpp_const_fields)
         
     ! combine interior and boundary layer coefficients and nonlocal term
     do ki= 1,km
        if(ki.lt.kbl) then
           kpp_1d_fields%difm(ki)=blmc(ki,1)
           kpp_1d_fields%difs(ki)=blmc(ki,2)
           kpp_1d_fields%dift(ki)=blmc(ki,3)
        else
           kpp_1d_fields%ghat(ki)=0.
        endif
     enddo
         
!     NPK 25/9/08.  Trap for negative values of diffusivities.
!     If negative, set to a background value of 1E-05.
!     
!     DO 205 ki= 1,km
!     DO 206 i = ipt,ipt
!     IF (dift(i,ki) .LT. 0) dift(i,ki)=1E-05
!     IF (difs(i,ki) .LT. 0) difs(i,ki)=1E-05
!     IF (visc(i,ki) .LT. 0) visc(i,ki)=1E-05
!     206     continue
!     205  continue
     
  ENDIF                     ! of LKPP
      
  return
end SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_KPPMIX



