SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_DDMIX(km, kmp1, alphaDT,betaDS,kpp_1d_fields)

#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_1d_type
#else 
  USE mckpp_data_fields, ONLY: kpp_1d_type
#endif

  IMPLICIT NONE
  
  ! Rho dependent interior flux parameterization.
  ! Add double-diffusion diffusivities to Ri-mix values at blending
  ! interface and below.
  
  ! Necessary for IMPLICIT NONE (NPK 11/2/13)
  integer km,kmp1,ki
  real dsfmax,rrho0
  
  real alphaDT(kmp1)  ! alpha * DT  across interfaces
  real betaDS(kmp1)   ! beta  * DS  across interfaces

  TYPE(kpp_1d_type) :: kpp_1d_fields      

  !local
  real Rrho              ! dd parameter
  real diffdd            ! double diffusion diffusivity scale
  real prandtl           ! prandtl number

  data Rrho0  /  1.9   / ! Rp=(alpha*delT)/(beta*delS)
  data dsfmax / 1.0e-4 / ! .0001 m2/s

  DO ki= 1, km           
     ! salt fingering case
     if((alphaDT(ki).gt.betaDS(ki)).and.(betaDS(ki).gt.0.)) then
        Rrho  = MIN(alphaDT(ki) / betaDS(ki) , Rrho0)
        diffdd     =         1.0-((Rrho-1)/(Rrho0-1))**2
        diffdd     = dsfmax*diffdd*diffdd*diffdd
        kpp_1d_fields%dift(ki) = kpp_1d_fields%dift(ki) + diffdd * 0.8 / Rrho
        kpp_1d_fields%difs(ki) = kpp_1d_fields%difs(ki) + diffdd
        
        ! diffusive convection
     else if ((alphaDT(ki).lt.0.0).and.(betaDS(ki).lt.0.0).and. &
          (alphaDT(ki).lt.betaDS(ki)) ) then
        Rrho    = alphaDT(ki) / betaDS(ki) 
        diffdd  = 1.5e-6*9.0*0.101*exp(4.6*exp(-0.54*(1/Rrho-1)))
        prandtl = 0.15*Rrho
        if (Rrho.gt.0.5) prandtl = (1.85-0.85/Rrho)*Rrho
        kpp_1d_fields%dift(ki) = kpp_1d_fields%dift(ki) + diffdd
        kpp_1d_fields%difs(ki) = kpp_1d_fields%difs(ki) + prandtl*diffdd
        
     endif
  ENDDO
  
END SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_DDMIX
