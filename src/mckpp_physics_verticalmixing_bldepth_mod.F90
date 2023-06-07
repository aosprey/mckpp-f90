MODULE mckpp_physics_verticalmixing_bldepth_mod

  USE mckpp_data_fields, ONLY: kpp_1d_type, kpp_const_type
  USE mckpp_physics_swfrac_mod, ONLY: mckpp_physics_swfrac, mckpp_physics_swfrac_opt
  USE mckpp_physics_verticalmixing_wscale_mod, ONLY: mckpp_physics_verticalmixing_wscale
  USE mckpp_time_control, ONLY: ntime

  IMPLICIT NONE

CONTAINS

!     the oceanic planetray boundary layer depth, hbl, is determined as
!     the shallowest depth where the bulk richardson number is
!     equal to the critical value, Ricr.
!     
!     bulk richardson numbers are evaluated by computing velocity and
!     buoyancy differences between values at zgrid(kl) < 0 and surface
!     reference values.
!     in this configuration, the reference values are equal to the
!     values in the surface layer.  
!     when using a very fine vertical grid, these values should be 
!     computed as the vertical average of velocity and buoyancy from 
!     the surface down to epsilon*zgrid(kl).
!
!     when the bulk richardson number at k exceeds Ricr, hbl is
!     linearly interpolated between grid levels zgrid(k) and zgrid(k-1).
!
!     The water column and the surface forcing are diagnosed for 
!     stable/ustable forcing conditions, and where hbl is relative 
!     to grid points (caseA), so that conditional branches can be 
!     avoided in later subroutines.
SUBROUTINE mckpp_physics_verticalmixing_bldepth (km, kmp1, dVsq, Ritop, ustar, Bo, Bosol, hbl, & 
    bfsfc, stable, caseA, kbl, Rib, sigma, wm, ws, kpp_1d_fields, kpp_const_fields)
    
  real bvsq,cekman,cmonob,cs,cv,epsilon,fekman,fmonob,&
       hbf,hekman,hmin,hmin2,hmonob,hri,Ricr,vtc,vtsq,epsln
  integer ka,ksave,ku,kl

  integer km,kmp1           ! number of vertical levels
  
  ! input
  real dVsq(kmp1)           ! (velocity shear re sfc)^2      (m/s)^2
  real Ritop(km)            ! numerator of bulk Richardson Number (m/s)^2
  !     Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc(m/s^2)
  real ustar                ! surface friction velocity         (m/s)
  real Bo                   ! surface turbulent buoyancy forcing(m^2/s^3)
  real Bosol                ! radiative buoyancy forcing        (m^2/s^3)
  type(kpp_1d_type) :: kpp_1d_fields
  type(kpp_const_type) :: kpp_const_fields
     
  ! output
  real hbl                  ! boundary layer depth              (m)
  real bfsfc                ! Bo+radiation absorbed to d=hbf*hbl(m^2/s^3)
  real stable               ! =1 in stable forcing; =0 unstable
  real caseA                ! =1 in case A, =0 in case B 
  integer kbl               ! index of first grid level below hbl 
       
  ! local
  real Rib(2)               ! Bulk Richardson number
  real sigma                ! normalized depth (d/hbl)
  real wm,ws                ! turbulent velocity scales         (m/s)
  real dmo(2)               ! Monin-Obukhov Depth
  real hek                  ! Ekman depth
  logical LEK, LMO          ! flags for MO and Ekman depth checks
      
  epsln   = 1.e-16
  Ricr    = 0.30
  epsilon = 0.1 
  cekman  = 0.7 
  cmonob  = 1.0 
  cs      = 98.96
  cv      = 1.6
  hbf     = 1.0     

  ! Set MO and Ekman depth flags
  LEK = .true. 
  LMO = .true. 
  
  ! find bulk Richardson number at every grid level find implied hri
  ! find Monin-Obukvov depth at every grid level find L
  ! Compare hri, L and hek to give hbl and kbl
  !   
  ! note: the reference depth is -epsilon/2.*zgrid(k), but the reference
  ! u,v,t,s values are simply the surface layer values,
  ! and not the averaged values from 0 to 2*ref.depth,
  ! which is necessary for very fine grids(top layer < 2m thickness)
  ! note: max values when Ricr never satisfied are
  ! kbl(i)=km and hbl(i) -zgrid(km)
  ! min values are                 kbl(i)=2      hbl(i) -zgrid(1)
      
  Vtc =  cv * sqrt(0.2/cs/epsilon) / kpp_const_fields%vonk**2 / Ricr
      
  ! indices for array Rib(i,k), the bulk Richardson number.
  ka = 1
  ku = 2
  
      
  ! initialize hbl and kbl to bottomed out values
  Rib(ka) = 0.0
  dmo(ka) = -kpp_const_fields%zm(kmp1)
  kbl = km
  hbl = -kpp_const_fields%zm(km)
  hek =  cekman * ustar / (abs(kpp_1d_fields%f) + epsln)
      
  do kl = 2,km
     ! compute bfsfc = sw fraction at hbf * zgrid
     ! To optimize the code choose the "swfrac_opt" ?
     ! call swfrac(imt,hbf,zgrid(kl),jwtype,bfsfc)
         
     ! Replaces the IF test at the beginning of the swfrac_opt
     ! subroutine.  The value will be stored in kpp_1d_field%swfrac
     ! thereafter (see below).
     if (ntime .le. 1 .and. kl.eq.2) then
        call MCKPP_PHYSICS_SWFRAC_OPT(hbf,kpp_1d_fields,kpp_const_fields)
     endif
                    
     IF(kbl.ge.km) THEN
        ! use caseA as temporary array for next call to wscale
        caseA = -kpp_const_fields%zm(kl)
            
        ! compute bfsfc= Bo + radiative contribution down to hbf * hbl
        bfsfc  = Bo + Bosol * (1. - kpp_1d_fields%swfrac(kl))
        stable = 0.5 + SIGN( 0.5, bfsfc+epsln )
        sigma  = stable * 1. + (1.-stable) * epsilon
     ENDIF

     ! compute velocity scales at sigma, for hbl= caseA = -zgrid(kl)
     call MCKPP_PHYSICS_VERTICALMIXING_WSCALE(sigma, caseA, ustar, bfsfc, wm,ws,kpp_const_fields)
     
     IF(kbl.ge.km) THEN        
        ! compute the turbulent shear contribution to Rib
        bvsq =0.5*(kpp_1d_fields%dbloc(kl-1) / (kpp_const_fields%zm(kl-1)-kpp_const_fields%zm(kl))+& 
             kpp_1d_fields%dbloc(kl) / (kpp_const_fields%zm(kl)-kpp_const_fields%zm(kl+1)))
        Vtsq = -kpp_const_fields%zm(kl) * ws * sqrt(abs(bvsq)) * Vtc
        !     compute bulk Richardson number at new level, dunder
        Rib(ku) = Ritop(kl) / (dVsq(kl)+Vtsq+epsln)
        Rib(ku) = MAX( Rib(ku), Rib(ka) + epsln)
        !     linear interpolate to find hbl where Rib = Ricr
        hri   = -kpp_const_fields%zm(kl-1) + (kpp_const_fields%zm(kl-1)-kpp_const_fields%zm(kl)) * &
             (Ricr - Rib(ka)) / (Rib(ku)-Rib(ka))
        
        ! compute the Monin Obukov length scale 
        ! fmonob    = stable(i) * LMO
        fmonob    = stable * 1.0
        dmo(ku) = cmonob * ustar * ustar * ustar &
             / kpp_const_fields%vonk / (abs(bfsfc) + epsln)
        dmo(ku) = fmonob * dmo(ku) - (1.-fmonob) * kpp_const_fields%zm(kmp1) 
        if(dmo(ku).le.(-kpp_const_fields%zm(kl))) then
           hmonob =(dmo(ku)-dmo(ka))/(kpp_const_fields%zm(kl-1)-kpp_const_fields%zm(kl))
           hmonob =(dmo(ku)+hmonob*kpp_const_fields%zm(kl)) / (1.-hmonob)
        else
           hmonob = -kpp_const_fields%zm(kmp1)
        endif
        
        ! compute the Ekman depth
        ! fekman  =  stable(i) * LEK
        fekman  =  stable * 1.0
        hekman  = fekman * hek - (1.-fekman) * kpp_const_fields%zm(kmp1)
            
        ! compute boundary layer depth
        hmin  = MIN(hri, hmonob,  hekman, -kpp_1d_fields%ocdepth)
        if(hmin .lt. -kpp_const_fields%zm(kl) ) then
     
           ! Code below added by SJW 09/07/04 to solve problems where hek 
           ! less than zgrid(kl-1) giving negative diffusions
           ! if this occurs to often then we need to rethink this fix
     
           ! Modified by NPK 25/09/08 to include Monin-Obukov depth as well.
           ! Richardson depth is sometimes calculated to be huge (> 1E10) when 
           ! Ritop is negative or very small and so is not always helpful in 
           ! this scenario.
           
           if (.not. kpp_1d_fields%l_initflag) then
              if (hmin .lt. -kpp_const_fields%zm(kl-1)) then
                 hmin2=MIN(hri,hmonob,-kpp_1d_fields%ocdepth)
                 if (hmin2 .lt. -kpp_const_fields%zm(kl)) THEN
                    hmin=hmin2
                 endif
              endif
           endif
           
           hbl = hmin
           kbl = kl
        endif
        
     ENDIF                  !kbl
     
     ksave = ka
     ka    = ku
     ku    = ksave     
  ENDDO
  
  call MCKPP_PHYSICS_SWFRAC(-1.0,hbl,kpp_1d_fields%jerlov,bfsfc)
  
  bfsfc  = Bo + Bosol * (1. - bfsfc)
  stable = 0.5 + SIGN( 0.5, bfsfc)
  bfsfc  = bfsfc + stable * epsln !ensures bfsfc never=0
  
  
  ! determine caseA and caseB
  caseA  = 0.5 + SIGN( 0.5,-kpp_const_fields%zm(kbl) -0.5 * kpp_const_fields%hm(kbl) -hbl)

end SUBROUTINE mckpp_physics_verticalmixing_bldepth

END MODULE mckpp_physics_verticalmixing_bldepth_mod
