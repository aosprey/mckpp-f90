MODULE mckpp_physics_verticalmixing_wscale_mod

  USE mckpp_data_fields, ONLY: kpp_const_type

  IMPLICIT NONE

CONTAINS

! compute turbulent velocity scales.
! use a 2D-lookup table for wm and ws as functions of ustar and
! zetahat (=vonk*sigma*hbl*bfsfc).
SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_WSCALE(sigma, hbl, ustar, bfsfc, wm, ws, kpp_const_fields)
  
  INTEGER i,iz,izp1,j,ju,jup1
  REAL am,as,c1,c2,c3,cm,cs,epsln,fzfrac,ucube,udiff,ufrac,&
       usta,wam,was,wbm,wbs,zdiff,zetas,zfrac,zetam
  
  TYPE(kpp_const_type) :: kpp_const_fields

  ! lookup table
  INTEGER :: ni = 890, &   ! number of values for zehat
             nj = 48       ! number of values for ustar

  real deltaz               ! delta zehat in table
  real deltau               ! delta ustar in table
  real zmin, zmax           ! zehat limits for table
  real umin, umax           ! ustar limits for table
     
  ! input
  real sigma                ! normalized depth (d/hbl)
  real hbl                  ! boundary layer depth (m)
  real ustar                ! surface friction velocity         (m/s)
  real bfsfc                ! total surface buoyancy flux       (m^2/s^3)
      
  ! output
  real wm,ws                ! turbulent velocity scales at sigma
      
  ! local
  real zehat                ! = zeta *  ustar**3
  real zeta                 ! = stability parameter d/L

  zmin   = -4.e-7
  zmax   = 0.0
  umin   = 0.0
  umax   = 0.04
  epsln  = 1.0e-20
  c1     = 5.0
  am     = 1.257
  cm     = 8.380
  c2     = 16.0
  zetam  = -0.2
  as     = -28.86
  cs     = 98.96
  c3     = 16.0
  zetas  = -1.0

  deltaz = (zmax-zmin)/(ni+1) 
  deltau = (umax-umin)/(nj+1)
  
  ! use lookup table for zehat < zmax  ONLY;  otherwise use stable formulae
  zehat = kpp_const_fields%vonk * sigma * hbl * bfsfc
  
  IF (zehat .le. zmax) THEN
     zdiff  = zehat-zmin
     iz = int( zdiff/deltaz )
     iz = min( iz , ni )
     iz = max( iz , 0  )
     izp1=iz+1
     
     udiff  = ustar-umin
     ju = int( udiff/deltau)
     ju = min( ju , nj )
     ju = max( ju , 0  )
     jup1=ju+1
     
     zfrac = zdiff/deltaz - float(iz)
     ufrac = udiff/deltau - float(ju)
     
     fzfrac= 1.-zfrac
     wam   = (fzfrac)  * kpp_const_fields%wmt(iz,jup1) + &
          zfrac*kpp_const_fields%wmt(izp1,jup1)
     wbm   = (fzfrac)  * kpp_const_fields%wmt(iz,ju  ) + &
          zfrac*kpp_const_fields%wmt(izp1,ju  )
     wm    = (1.-ufrac)* wbm          + ufrac*wam
         
     was   = (fzfrac)  * kpp_const_fields%wst(iz,jup1) + &
          zfrac*kpp_const_fields%wst(izp1,jup1)
     wbs   = (fzfrac)  * kpp_const_fields%wst(iz,ju  ) + &
          zfrac*kpp_const_fields%wst(izp1,ju  )
     ws    = (1.-ufrac)* wbs          + ufrac*was       
  ELSE
     ucube = ustar**3
     wm = kpp_const_fields%vonk * ustar * ucube / (ucube + c1 * zehat)
     ws = wm     
  ENDIF
  
END SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_WSCALE

END MODULE mckpp_physics_verticalmixing_wscale_mod
