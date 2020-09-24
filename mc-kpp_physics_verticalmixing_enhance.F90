SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_ENHANCE(km, mdiff, dkm1, hbl, kbl,caseA, blmc,&
     kpp_1d_fields, kpp_const_fields)
#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_1d_type,kpp_const_type
#endif
  IMPLICIT NONE

  ! enhance the diffusivity at the kbl-.5 interface

#ifdef MCKPP_CAM3
#include <parameter.inc>
#else
  ! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
#endif

  !     Necessary for IMPLICIT NONE (NPK 11/2/13)
  real dkmp5,dstar
  integer ki
  integer km                !,kmp1           ! number of vertical levels  
  integer mdiff             ! number of viscosities + diffusivities
  integer kbl               ! grid above hbl
  real hbl                  ! boundary layer depth             (m)
  real dkm1(mdiff)          ! bl diffusivity at kbl-1 grid level
  real caseA                ! = 1 in caseA, = 0 in case B
 
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  ! output
  real blmc(km,mdiff)       ! enhanced bound. layer mixing coeff.

  ! local
  real delta                ! fraction hbl lies beteen zgrid neighbors

  do ki=1,km-1         
     if(ki .eq. (kbl - 1) ) then            
        delta = (hbl+kpp_const_fields%zm(ki)) / (kpp_const_fields%zm(ki)-kpp_const_fields%zm(ki+1))
            
        dkmp5 = caseA * kpp_1d_fields%difm(ki) + (1.-caseA) * blmc(ki,1)
        dstar = (1.-delta)**2 * dkm1(1) + delta**2 * dkmp5      
        blmc(ki,1) = (1.-delta) * kpp_1d_fields%difm(ki) + delta * dstar
        
        dkmp5 = caseA * kpp_1d_fields%difs(ki) + (1.-caseA) * blmc(ki,2)
        dstar = (1.-delta)**2 * dkm1(2) + delta**2 * dkmp5    
        blmc(ki,2) = (1.-delta) * kpp_1d_fields%difs(ki) + delta * dstar
        
        dkmp5 = caseA * kpp_1d_fields%dift(ki) + (1.-caseA) * blmc(ki,3)
        dstar = (1.-delta)**2 * dkm1(3) + delta**2 * dkmp5     
        blmc(ki,3) = (1.-delta) * kpp_1d_fields%dift(ki) + delta * dstar
        
        kpp_1d_fields%ghat(ki) = (1.-caseA) * kpp_1d_fields%ghat(ki)            
     endif
  enddo
  
  return
end SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_ENHANCE
