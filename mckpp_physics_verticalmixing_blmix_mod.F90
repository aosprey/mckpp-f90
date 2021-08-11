MODULE mckpp_physics_verticalmixing_blmix_mod

CONTAINS

SUBROUTINE mckpp_physics_verticalmixing_blmix(km, mdiff, ustar, bfsfc, hbl, stable, caseA, kbl, &
    gat1, dat1, dkm1, blmc, sigma, wm, ws, kpp_1d_fields,kpp_const_fields)
  
#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_1d_type,kpp_const_type
#else
  USE mckpp_data_fields, ONLY: kpp_1d_type,kpp_const_type
#endif

  ! mixing coefficients within boundary layer depend on surface
  ! forcing and the magnitude and gradient of interior mixing below
  ! the boundary layer ("matching").

  IMPLICIT NONE

  ! CAUTION if mixing bottoms out at hbl = -zgrid(km) THEN
  ! fictious layer kmp1 is needed with small but finite width (eg. 1.e-10)

  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  integer km                !,kmp1        ! number of vertical levels
  integer mdiff             ! number of viscosities + diffusivities
  
  !    input
  real ustar                ! surface friction velocity         (m/s)
  real bfsfc                ! surface buoyancy forcing        (m^2/s^3)
  real hbl                  ! boundary layer depth              (m)
  real stable               ! = 1 in stable forcing
  real caseA                ! = 1 in case A
  integer kbl               ! index of first grid level below hbl
     
  ! output
  real gat1(mdiff)
  real dat1(mdiff)
  real dkm1(mdiff)          ! boundary layer difs at kbl-1 level
  real blmc(km,mdiff)       ! boundary layer mixing coefficients(m^2/s)
  
  !  local
  real sigma                ! normalized depth (d / hbl)
  real ws, wm               ! turbulent velocity scales         (m/s)
  !  None of these were previously declared ... (NPK 6/2/13)
  real a1,a2,a3,am,as,c1,c2,c3,cg,cm,cs,cstar,delhat,difsh,&
       difsp,difth,diftp,dvdzup,epsln,f1,gm,gs,&
       dvdzdn,epsilon,gt,r,visch,viscp,zetam,sig,&
       zetas
  integer ki,kn
 
  data epsln             /   1.e-20 /
  data epsilon           /   0.1    /
  data c1                /   5.0    /
  data am,cm,c2,zetam    /   1.257  ,  8.380, 16.0, - 0.2 / !7-24-92
  data as,cs,c3,zetas    / -28.86   , 98.96 , 16.0, - 1.0 /
  data cstar             /    5.    /     

  cg = cstar * kpp_const_fields%vonk * (cs * kpp_const_fields%vonk * epsilon)**(1./3.)
      
  ! compute velocity scales at hbl
  sigma = stable * 1.0 + (1.-stable) * epsilon
      
  CALL MCKPP_PHYSICS_VERTICALMIXING_WSCALE(sigma, hbl, ustar, bfsfc,wm,ws,kpp_const_fields)
  kn = ifix(caseA+epsln) *(kbl -1) + (1-ifix(caseA+epsln)) * kbl
      
  ! find the interior viscosities and derivatives at hbl(i) 
  delhat = 0.5*kpp_const_fields%hm(kn)-kpp_const_fields%zm(kn) - hbl
  R      = 1.0 - delhat / kpp_const_fields%hm(kn)
  dvdzup = (kpp_1d_fields%difm(kn-1) - kpp_1d_fields%difm(kn)) / kpp_const_fields%hm(kn) 
  dvdzdn = (kpp_1d_fields%difm(kn)   - kpp_1d_fields%difm(kn+1)) / kpp_const_fields%hm(kn+1)
  viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup)) + R * (dvdzdn + abs(dvdzdn)) )
      
  dvdzup = (kpp_1d_fields%difs(kn-1) - kpp_1d_fields%difs(kn)) / kpp_const_fields%hm(kn) 
  dvdzdn = (kpp_1d_fields%difs(kn)   - kpp_1d_fields%difs(kn+1)) / kpp_const_fields%hm(kn+1)
  difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup)) + R * (dvdzdn + abs(dvdzdn)) )
      
  dvdzup = (kpp_1d_fields%dift(kn-1) - kpp_1d_fields%dift(kn)) / kpp_const_fields%hm(kn) 
  dvdzdn = (kpp_1d_fields%dift(kn)   - kpp_1d_fields%dift(kn+1)) / kpp_const_fields%hm(kn+1)
  diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup)) + R * (dvdzdn + abs(dvdzdn)) )
     
  visch  = kpp_1d_fields%difm(kn) + viscp * delhat
  difsh  = kpp_1d_fields%difs(kn) + difsp * delhat
  difth  = kpp_1d_fields%dift(kn) + diftp * delhat
      
  f1 = stable * c1 * bfsfc / (ustar**4+epsln) 
  gat1(1) = visch / hbl / (wm+epsln)
  dat1(1) = -viscp / (wm+epsln) + f1 * visch
  dat1(1) = min(dat1(1),0.) 
  
  gat1(2) = difsh  / hbl / (ws+epsln)
  dat1(2) = -difsp / (ws+epsln) + f1 * difsh 
  dat1(2) = min(dat1(2),0.) 
  
  gat1(3) = difth /  hbl / (ws+epsln)
  dat1(3) = -diftp / (ws+epsln) + f1 * difth 
  dat1(3) = min(dat1(3),0.) 
      
!     Turn off interior matching here
!     gat1(i,1) = 0.0001
!     gat1(i,2) = 0.00001
!     gat1(i,3) = 0.00001
!     do m=1,3
!       dat1(i,m) = 0.0
!       enddo     

  DO ki = 1,km       
     ! compute turbulent velocity scales on the interfaces
     sig = (-kpp_const_fields%zm(ki) + 0.5 * kpp_const_fields%hm(ki)) / hbl
     sigma   = stable*sig + (1.-stable)*AMIN1(sig,epsilon)
     CALL MCKPP_PHYSICS_VERTICALMIXING_WSCALE(sigma, hbl, ustar, bfsfc,wm,ws,kpp_const_fields)
     
     ! compute the dimensionless shape functions at the interfaces     
     sig = (-kpp_const_fields%zm(ki) + 0.5 * kpp_const_fields%hm(ki)) / hbl
     a1 = sig - 2.
     a2 = 3.-2.*sig
     a3 = sig - 1.
   
     Gm = a1 + a2 * gat1(1) + a3 * dat1(1) 
     Gs = a1 + a2 * gat1(2) + a3 * dat1(2)
     Gt = a1 + a2 * gat1(3) + a3 * dat1(3)
   
     !  compute boundary layer diffusivities at the interfaces 
     blmc(ki,1) = hbl * wm * sig * (1. + sig * Gm)
     blmc(ki,2) = hbl * ws * sig * (1. + sig * Gs)
     blmc(ki,3) = hbl * ws * sig * (1. + sig * Gt)
     
     ! nonlocal transport term = ghats * <ws>o
     kpp_1d_fields%ghat(ki) = (1.-stable) * cg / (ws*hbl+epsln)
  ENDDO
 
  ! find diffusivities at kbl-1 grid level 
  sig   =  -kpp_const_fields%zm(kbl-1)  / hbl
  sigma =  stable * sig + (1.-stable) * AMIN1(sig,epsilon)
  
  CALL MCKPP_PHYSICS_VERTICALMIXING_WSCALE(sigma, hbl, ustar, bfsfc, wm, ws, kpp_const_fields)
  sig = -kpp_const_fields%zm(kbl-1) / hbl
  a1= sig - 2.
  a2 = 3.-2.*sig
  a3 = sig - 1.
  Gm = a1 + a2 * gat1(1) + a3 * dat1(1)
  Gs = a1 + a2 * gat1(2) + a3 * dat1(2)
  Gt = a1 + a2 * gat1(3) + a3 * dat1(3)
  dkm1(1) = hbl * wm * sig * (1. + sig * Gm)
  dkm1(2) = hbl * ws * sig * (1. + sig * Gs)
  dkm1(3) = hbl * ws * sig * (1. + sig * Gt)
  
END SUBROUTINE mckpp_physics_verticalmixing_blmix

END MODULE mckpp_physics_verticalmixing_blmix_mod
