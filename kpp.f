      subroutine kppmix(
     $     dVsq,
     $     ustar , Bo    , Bosol , alphaDT , betaDS ,
     $     Ritop, 
     $     hbl , kbl, 
     $     kpp_2d_fields,kpp_const_fields)

      IMPLICIT NONE

c.......................................................................
c
c     Main driver subroutine for kpp vertical mixing scheme and 
c     interface to greater ocean model
c
c     written by : bill large,   june,  6, 1994
c     modified by: jan morzel,   june, 30, 1994
c                  bill large, august, 11, 1994
c                  bill large, november 1994,   for 1d code
c
c.......................................................................
c     
! Automatically includes parameter.inc
#include "kpp_3d_type.com"

      integer km,kmp1,mdiff,ki

      parameter (km = NZ, kmp1 = nzp1)!, imt = 1) !NX*NY)
      parameter (mdiff = 3)  ! number of diffusivities for local arrays
c   include 'local_pt.com'
c  #include "landsea.com>
c
c input
c      real zgrid(kmp1)          ! vertical grid (<= 0)            (m)
c      real hwide(kmp1)          ! layer thicknesses               (m)
c      real Shsq_oned(kmp1)      ! (local velocity shear)^2       (m/s)^2
      real dVsq(kmp1)           ! (velocity shear re sfc)^2      (m/s)^2
      real ustar                ! surface friction velocity       (m/s)
      real Bo                   ! surface turbulent buoy. forcing (m^2/s^3)
      real Bosol                ! radiative buoyancy forcing      (m^2/s^3)
      real alphaDT(kmp1)        ! alpha * DT  across interfaces
      real betaDS(kmp1)         ! beta  * DS  across interfaces
c      real dbloc_oned(km)        ! local delta buoyancy            (m/s^2)
      real Ritop(km)            ! numerator of bulk Richardson Number (m/s)^2
c     Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc(m/s^2)
c     real Coriol               ! Coriolis parameter              (s^{-1})
c     integer jwtype            ! Jerlov water type               (1 -- 5)
c      logical LRI, LDD, LKPP    ! mixing process switches
c      real ocdepth
      type(kpp_2d_type) :: kpp_2d_fields
      type(kpp_const_type) :: kpp_const_fields
c
c output
c     visc replaced by kpp_2d_fields%difm (NPK 8/2/2013)
c      real visc(0:kmp1)         ! vertical viscosity coefficient  (m^2/s)
c      real difs(0:kmp1)         ! vertical scalar diffusivity     (m^2/s)
c      real dift(0:kmp1)         ! vertical temperature diffusivity(m^2/s)
c      real ghats(km)            ! nonlocal transport              (s/m^2)
      real hbl                  ! boundary layer depth (m)
      integer kbl               ! index of first grid level below hbl
c      real Rig_oned(km)         ! local Richardson number (NPK diagnostic)
c     
c     local
      real bfsfc                ! surface buoyancy forcing        (m^2/s^3)
      real ws                   ! momentum velocity scale
      real wm                   ! scalar   velocity scale 
      real caseA                ! = 1 in case A; =0 in case B
      real stable               ! = 1 in stable forcing; =0 in unstable
      real dkm1(mdiff)          ! boundary layer difs at kbl-1 level
      real gat1(mdiff)          ! shape function at sigma=1
      real dat1(mdiff)          ! derivative of shape function at sigma=1
      real blmc(km,mdiff)       ! boundary layer mixing coefficients
      real sigma                ! normalized depth (d / hbl)
      real Rib(2)               ! bulk Richardson number
      
c      WRITE(6,*) 'In kppmix, hwide= ',hwide

c     zero the mixing coefficients 
      DO ki=0,km
         kpp_2d_fields%difm(ki) = 0.0
         kpp_2d_fields%difs(ki) = 0.0
         kpp_2d_fields%dift(ki) = 0.0
      END DO
      
c     compute RI and IW interior diffusivities everywhere
      IF(kpp_const_fields%LRI) THEN
         call ri_iwmix (km,kmp1,kpp_2d_fields,kpp_const_fields)
      ENDIF
      

c     add double diffusion if desired
      IF(kpp_const_fields%LDD) THEN
         call ddmix ( km  , kmp1,
     $        alphaDT, betaDS , kpp_2d_fields)
      ENDIF

c fill the bottom kmp1 coefficients for blmix      
      kpp_2d_fields%difm(kmp1) = kpp_2d_fields%difm(km)
      kpp_2d_fields%difs(kmp1) = kpp_2d_fields%difs(km)
      kpp_2d_fields%dift(kmp1) = kpp_2d_fields%dift(km)
 
c  compute boundary layer mixing coefficients ??
      IF(kpp_const_fields%LKPP) THEN
c diagnose the new boundary layer depth

         call  bldepth (
     +        km   , kmp1, dVsq,   
     $        Ritop, ustar , Bo   , Bosol, 
     +        hbl  , bfsfc, stable, caseA, kbl  ,
     $        Rib  , sigma, wm    , ws   , kpp_2d_fields, 
     +        kpp_const_fields)
                           
c     compute boundary layer diffusivities
         call blmix   (km   , mdiff,
     $        ustar, bfsfc, hbl , stable, caseA, 
     $        kbl, 
     $        gat1 , dat1 , dkm1, blmc, 
     $        sigma, wm   , ws, kpp_2d_fields,kpp_const_fields)
         
         
         
c     enhance diffusivity at interface kbl - 1
         call enhance (km   , mdiff , dkm1,
     $        hbl , kbl   , caseA,
     $        blmc, kpp_2d_fields,kpp_const_fields)
         

c     combine interior and boundary layer coefficients and nonlocal term
         do ki= 1,km
            if(ki.lt.kbl) then
               kpp_2d_fields%difm(ki)=blmc(ki,1)
               kpp_2d_fields%difs(ki)=blmc(ki,2)
               kpp_2d_fields%dift(ki)=blmc(ki,3)
            else
               kpp_2d_fields%ghat(ki)=0.
            endif
         enddo
        
c     For slab, set all values to small numbers
         IF (kpp_const_fields%L_SLAB) THEN
            kpp_2d_fields%difm(1:NZtmax)=1e-20
            kpp_2d_fields%dift(1:NZtmax)=1e-20
            kpp_2d_fields%difs(1:NZtmax)=1e-20
         ENDIF
         

c     
c     NPK 25/9/08.  Trap for negative values of diffusivities.
c     If negative, set to a background value of 1E-05.
c     
c     DO 205 ki= 1,km
c     DO 206 i = ipt,ipt
c     IF (dift(i,ki) .LT. 0) dift(i,ki)=1E-05
c     IF (difs(i,ki) .LT. 0) difs(i,ki)=1E-05
c     IF (visc(i,ki) .LT. 0) visc(i,ki)=1E-05
c     206     continue
c     205  continue
         
         
      ENDIF                     ! of LKPP
      
      return
      end
      
c ********************************************************************

      subroutine  bldepth (
     $     km   , kmp1 ,  dVsq , 
     $     Ritop, ustar , Bo   , Bosol, 
     $     hbl  , bfsfc, stable, caseA, kbl  ,
     $     Rib  , sigma, wm    , ws   , kpp_2d_fields, kpp_const_fields)

      IMPLICIT NONE
c     
c     the oceanic planetray boundary layer depth, hbl, is determined as
c     the shallowest depth where the bulk richardson number is
c     equal to the critical value, Ricr.
c     
c     bulk richardson numbers are evaluated by computing velocity and
c     buoyancy differences between values at zgrid(kl) < 0 and surface
c     reference values.
c     in this configuration, the reference values are equal to the
c     values in the surface layer.  
c     when using a very fine vertical grid, these values should be 
c     computed as the vertical average of velocity and buoyancy from 
c     the surface down to epsilon*zgrid(kl).
c
c     when the bulk richardson number at k exceeds Ricr, hbl is
c     linearly interpolated between grid levels zgrid(k) and zgrid(k-1).
c
c     The water column and the surface forcing are diagnosed for 
c     stable/ustable forcing conditions, and where hbl is relative 
c     to grid points (caseA), so that conditional branches can be 
c     avoided in later subroutines.
c
c
c     model  
c      include 'local_pt.com'
c      include 'times.com'
#include "kpp_3d_type.com"

c     Necessary for IMPLICIT NONE
      real bvsq,cekman,cmonob,cs,cv,epsilon,fekman,fmonob,
     +     hbf,hekman,hmin,hmin2,hmonob,hri,Ricr,vtc,vtsq,epsln
      integer ka,ksave,ku,kl

      integer km,kmp1           ! number of vertical levels
c      integer imt               ! number of horizontal grid points
c      real zgrid(kmp1)          ! vertical grid (<= 0)              (m)
c      real hwide(kmp1)          ! layer thicknesses                 (m)
c     
c     input
      real dVsq(kmp1)           ! (velocity shear re sfc)^2      (m/s)^2
c      real dbloc(km)            ! local delta buoyancy              (m/s^2)
      real Ritop(km)            ! numerator of bulk Richardson Number (m/s)^2
c     Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc(m/s^2)
      real ustar                ! surface friction velocity         (m/s)
      real Bo                   ! surface turbulent buoyancy forcing(m^2/s^3)
      real Bosol                ! radiative buoyancy forcing        (m^2/s^3)
c     real Coriol               ! Coriolis parameter                (1/s)
c      integer jwtype            ! Jerlov water type                 (1 to 5)
c      real ocdepth
      type(kpp_2d_type) :: kpp_2d_fields
      type(kpp_const_type) :: kpp_const_fields
c     
c     output
      real hbl                  ! boundary layer depth              (m)
      real bfsfc                ! Bo+radiation absorbed to d=hbf*hbl(m^2/s^3)
      real stable               ! =1 in stable forcing; =0 unstable
      real caseA                ! =1 in case A, =0 in case B 
      integer kbl               ! index of first grid level below hbl 
c     
c     local
      real Rib(2)               ! Bulk Richardson number
      real sigma                ! normalized depth (d/hbl)
      real wm,ws                ! turbulent velocity scales         (m/s)
      real dmo(2)               ! Monin-Obukhov Depth
      real hek                  ! Ekman depth
      logical LEK, LMO          ! flags for MO and Ekman depth checks
      
c      LOGICAL L_INITFLAG
c      common /initflag/L_INITFLAG
      
c      save epsln,Ricr,epsilon,cekman,cmonob,cs,cv,hbf
c     
      data epsln           /  1.e-16 /
c      data DelVmin         /  .005   /
      data Ricr            /  0.30   /
      data epsilon         /  0.1    /
      data cekman          /  0.7    /
      data cmonob          /  1.0    /
      data cs              / 98.96   /
      data cv              /  1.6    /
      data hbf             /  1.0    /     

c     Set MO and Ekman depth flags
      LEK = .true. 
      LMO = .true. 
 
c     find bulk Richardson number at every grid level find implied hri
c     find Monin-Obukvov depth at every grid level find L
c     Compare hri, L and hek to give hbl and kbl
c     
c     note: the reference depth is -epsilon/2.*zgrid(k), but the reference
c     u,v,t,s values are simply the surface layer values,
c     and not the averaged values from 0 to 2*ref.depth,
c     which is necessary for very fine grids(top layer < 2m thickness)
c     note: max values when Ricr never satisfied are
c     kbl(i)=km and hbl(i) -zgrid(km)
c     min values are                 kbl(i)=2      hbl(i) -zgrid(1)
      
      Vtc =  cv * sqrt(0.2/cs/epsilon) / kpp_const_fields%vonk**2 / Ricr
      
c     indices for array Rib(i,k), the bulk Richardson number.
      ka = 1
      ku = 2
      
      
c     initialize hbl and kbl to bottomed out values
      Rib(ka) = 0.0
      dmo(ka) = -kpp_2d_fields%zm(kmp1)
      kbl    = km
      hbl    = -kpp_2d_fields%zm(km)
c     Coriol(i) = 2. * (twopi/86164.) * sin(2.5*twopi/360.)
      hek =  cekman * ustar / 
     +     (abs(kpp_2d_fields%f) + epsln)
      
      do 200 kl = 2,km
c     compute bfsfc = sw fraction at hbf * zgrid
c     To optimize the code choose the "swfrac_opt" ?
c     call swfrac(imt,hbf,zgrid(kl),jwtype,bfsfc)
         
c     Replaces the IF test at the beginning of the swfrac_opt
c     subroutine.  The value will be stored in kpp_2d_field%swfrac
c     thereafter (see below).
         if (kpp_const_fields%ntime .le. 1 .and. kl.eq.2) then
            call swfrac_opt(hbf,
     +           kpp_2d_fields,kpp_const_fields)
         endif
                    
         IF(kbl.ge.km) THEN
c            WRITE(6,*) 'kbl = ',kbl
c     use caseA as temporary array for next call to wscale
            caseA = -kpp_2d_fields%zm(kl)
            
c     compute bfsfc= Bo + radiative contribution down to hbf * hbl
            bfsfc  = Bo
     $           + Bosol * (1. - kpp_2d_fields%swfrac(kl))
            stable = 0.5 + SIGN( 0.5, bfsfc+epsln )
            sigma  = stable * 1. + (1.-stable) * epsilon
c            WRITE(6,*) 'sigma = ',sigma
c            WRITE(6,*) 'hbl = ',hbl
c            WRITE(6,*) 'ustar = ',ustar
c            WRITE(6,*) 'bfsfc = ',bfsfc
         ENDIF

c        compute velocity scales at sigma, for hbl= caseA = -zgrid(kl)
c         WRITE(6,*) 'kbl = ',kbl
c         WRITE(6,*) 'sigma = ',sigma
c         WRITE(6,*) 'hbl = ',hbl
c         WRITE(6,*) 'ustar = ',ustar
c         WRITE(6,*) 'bfsfc = ',bfsfc
c         WRITE(6,*) 'wscale(',sigma,hbl,ustar,bfsfc

         call wscale(sigma, caseA, ustar, bfsfc, wm,ws,kpp_const_fields)
         
         IF(kbl.ge.km) THEN
            
c     compute the turbulent shear contribution to Rib
            bvsq =0.5*
     $           (kpp_2d_fields%dbloc(kl-1) / 
     $           (kpp_2d_fields%zm(kl-1)-kpp_2d_fields%zm(kl))+ 
     $           kpp_2d_fields%dbloc(kl) / 
     $           (kpp_2d_fields%zm(kl)-kpp_2d_fields%zm(kl+1)))
            Vtsq = -kpp_2d_fields%zm(kl) * ws * sqrt(abs(bvsq)) * Vtc
c     compute bulk Richardson number at new level, dunder
            Rib(ku) = Ritop(kl) / (dVsq(kl)+Vtsq+epsln)
            Rib(ku) = MAX( Rib(ku), Rib(ka) + epsln)
c     linear interpolate to find hbl where Rib = Ricr
            hri   = -kpp_2d_fields%zm(kl-1) + 
     +           (kpp_2d_fields%zm(kl-1)-kpp_2d_fields%zm(kl)) *
     $           (Ricr - Rib(ka)) / (Rib(ku)-Rib(ka))
            
c     compute the Monin Obukov length scale 
c     fmonob    = stable(i) * LMO
            fmonob    = stable * 1.0
            dmo(ku) = cmonob * ustar * ustar 
     +           * ustar
     >           / kpp_const_fields%vonk / (abs(bfsfc) + epsln)
            dmo(ku) = fmonob * dmo(ku) - (1.-fmonob) *
     +           kpp_2d_fields%zm(kmp1) 
            if(dmo(ku).le.(-kpp_2d_fields%zm(kl))) then
               hmonob =(dmo(ku)-dmo(ka))/(kpp_2d_fields%zm(kl-1)-
     +              kpp_2d_fields%zm(kl))
               hmonob =(dmo(ku)+hmonob*kpp_2d_fields%zm(kl)) / 
     +              (1.-hmonob)
            else
               hmonob = -kpp_2d_fields%zm(kmp1)
            endif
              
c     compute the Ekman depth
c     fekman  =  stable(i) * LEK
            fekman  =  stable * 1.0
            hekman  = fekman * hek - (1.-fekman) * 
     +           kpp_2d_fields%zm(kmp1)
            
c     compute boundary layer depth
            hmin  = MIN(hri, hmonob,  hekman, -kpp_2d_fields%ocdepth)
c     WRITE(6,*) 'hri=',hri,'hmonob=',hmonob,'hekman=',hekman
            if(hmin .lt. -kpp_2d_fields%zm(kl) ) then
c     
c     Code below added by SJW 09/07/04 to solve problems where hek 
c     less than zgrid(kl-1) giving negative diffusions
c     if this occurs to often then we need to rethink this fix
c     
c     Modified by NPK 25/09/08 to include Monin-Obukov depth as well.
c     Richardson depth is sometimes calculated to be huge (> 1E10) when 
c     Ritop is negative or very small and so is not always helpful in 
c     this scenario.
c     
               if (.not. kpp_2d_fields%l_initflag) then
                  if (hmin .lt. -kpp_2d_fields%zm(kl-1)) then
                     hmin2=MIN(hri,hmonob,-kpp_2d_fields%ocdepth)
                     if (hmin2 .lt. -kpp_2d_fields%zm(kl)) THEN
c     write(6,*) 'Setting hmin=',
c     &                         hmin2,'from hek? ',hekman,
c     &                         'hri=',hri,'hmonob=',hmonob
                        hmin=hmin2
c     else
c     write(6,*) 'Leaving hmin=hek ',hekman,
c     &                         'as hmonob= ',hmonob,'and hri= ',hri
                     endif
                  endif
               endif
               
               hbl = hmin
               kbl = kl
c     if(hmin.ge. hmonob) write(6,*) 'MO ',kl,hmin, hmonob,
c     &                hri
c     if(hmin.ge. hekman) write(6,*) 'EK ',kl,hmin, hekman,
c     &                hri
            endif
            
         ENDIF                  !kbl
         
!     Save hekman for use with Ekman pumping code in ocn.f
!         IF (kpp_const_fields%L_EKMAN_PUMP) 
!     c        kpp_2d_fields%hekman = hekman

         ksave = ka
         ka    = ku
         ku    = ksave
         
 200  continue
      
      call swfrac(-1.0,hbl,kpp_2d_fields,bfsfc)     
      
      bfsfc  = Bo + Bosol * (1. - bfsfc)
      stable = 0.5 + SIGN( 0.5, bfsfc)
      bfsfc  = bfsfc + stable * epsln !ensures bfsfc never=0
      
 
c determine caseA and caseB
      caseA  = 0.5 + 
     $     SIGN( 0.5,-kpp_2d_fields%zm(kbl) -0.5*
     +     kpp_2d_fields%hm(kbl) -hbl)


      return
      end

c *********************************************************************

      subroutine wscale(sigma, hbl, ustar, bfsfc,
     $     wm , ws, kpp_const_fields)
      IMPLICIT NONE
c
c     compute turbulent velocity scales.
c     use a 2D-lookup table for wm and ws as functions of ustar and
c     zetahat (=vonk*sigma*hbl*bfsfc).
c
c
! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
c
c     Necessary for IMPLICIT NONE (NPK 11/2/13)
      INTEGER ni,nj,i,iz,izp1,j,ju,jup1
      REAL am,as,c1,c2,c3,cm,cs,epsln,fzfrac,ucube,udiff,ufrac,
     +     usta,wam,was,wbm,wbs,zdiff,zetas,zfrac,zetam
c
      TYPE(kpp_const_type) :: kpp_const_fields
c
c lookup table
      parameter ( ni = 890,     ! number of values for zehat
     $     nj = 48)             ! number of values for ustar

c      real wmt(0:ni+1,0:nj+1)   ! lookup table for wm
c      real wst(0:ni+1,0:nj+1)   ! lookup table for ws
      real deltaz               ! delta zehat in table
      real deltau               ! delta ustar in table
      real zmin,zmax            ! zehat limits for table
      real umin,umax            ! ustar limits for table
c      logical firstf
c      save deltaz,deltau,zmin,zmax,umin,umax,firstf
c     
      data zmin,zmax  / -4.e-7, 0.0   / ! m3/s3
      data umin,umax  /  0.   , .04   / ! m/s
c      data firstf     / .true.        /
      
c     model
c      include 'local_pt.com'
c      integer imt               ! number of horizontal grid points
      
c     input
      real sigma                ! normalized depth (d/hbl)
      real hbl                  ! boundary layer depth (m)
      real ustar                ! surface friction velocity         (m/s)
      real bfsfc                ! total surface buoyancy flux       (m^2/s^3)
      
c     output
      real wm,ws                ! turbulent velocity scales at sigma
      
c     local
      real zehat                ! = zeta *  ustar**3
      real zeta                 ! = stability parameter d/L
      
c      save epsln,c1,am,cm,c2,zetam,as,cs,c3,zetas,vonk
      
      data epsln           /   1.0e-20/
      data c1              /   5.0   /
      data am,cm,c2,zetam  /   1.257 ,  8.380 , 16.0 , - 0.2  /
      data as,cs,c3,zetas  / -28.86  , 98.96  , 16.0 , - 1.0  /
c      data vonk            /   0.40  /
c     
c     construct the wm and ws lookup tables
c     
c     Moved to a separate subroutine for openMP compatability
c     NPK 13/2/2013.
c     
      deltaz = (zmax-zmin)/(ni+1) 
      deltau = (umax-umin)/(nj+1)
 
c use lookup table for zehat < zmax  ONLY;  otherwise use stable formulae
c      WRITE(6,*) 'vonk = ',vonk
c      WRITE(6,*) 'sigma = ',sigma
c      WRITE(6,*) 'hbl = ',hbl
c      WRITE(6,*) 'bfsfc = ',bfsfc
c      WRITE(6,*) 'zehat = ',vonk,sigma,hbl,bfsfc
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
         wam   = (fzfrac)  * kpp_const_fields%wmt(iz,jup1) + 
     +        zfrac*kpp_const_fields%wmt(izp1,jup1)
         wbm   = (fzfrac)  * kpp_const_fields%wmt(iz,ju  ) + 
     +        zfrac*kpp_const_fields%wmt(izp1,ju  )
         wm    = (1.-ufrac)* wbm          + ufrac*wam
         
         was   = (fzfrac)  * kpp_const_fields%wst(iz,jup1) + 
     +        zfrac*kpp_const_fields%wst(izp1,jup1)
         wbs   = (fzfrac)  * kpp_const_fields%wst(iz,ju  ) + 
     +        zfrac*kpp_const_fields%wst(izp1,ju  )
         ws    = (1.-ufrac)* wbs          + ufrac*was
         
      ELSE
         
         ucube = ustar**3
         wm = kpp_const_fields%vonk * ustar * ucube / 
     +        (ucube + c1 * zehat)
         ws = wm
         
      ENDIF         
      
      return
      end
      
c     **********************************************************************
      
      subroutine ri_iwmix (km,kmp1,kpp_2d_fields,kpp_const_fields)

      implicit none
c     
c     compute interior viscosity diffusivity coefficients due to
c     shear instability (dependent on a local richardson number)
c     and due to background internal wave activity.
c     
      
c     input
c      include 'local_pt.com'
! Automatically includes parameter.inc!
#include "kpp_3d_type.com"

      TYPE(kpp_2d_type) :: kpp_2d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

      integer km,kmp1           ! number of vertical levels
c      integer imt               ! number of horizontal grid points
c     real Shsq(kmp1)           ! (local velocity shear)^2          (m/s)^2
c     real dbloc(km)            ! local delta buoyancy              (m/s^2)
c     real zgrid(kmp1)          ! vertical grid (<= 0)              (m)
      
c     output
c     real visc(0:kmp1)         ! vertical viscosivity coefficient  (m^2/s)
c     real difs(0:kmp1)         ! vertical scalar diffusivity       (m^2/s)
c     real dift(0:kmp1)         ! vertical temperature diffusivity  (m^2/s)
c     real Rig(km)              ! local Richardson number
      
c     local variables
c     real Rig,Rigg             ! local richardson number
      real Rigg
      real fri,fcon             ! function of Rig
      real ratio
      real epsln,Riinfty,Ricon,difm0,difs0,difmiw,difsiw,difmcon,
     &     difscon,c1,c0
c      integer i,ki,mr,mRi
      INTEGER ki,mRi,j

c      save epsln,Riinfty,Ricon,difm0,difs0,difmcon,difscon,
c     &     difmiw,difsiw,c1
      
      data  epsln   / 1.e-16 /  ! a small number          
      data  Riinfty /  0.8     / ! LMD default was = 0.7
      data  Ricon   / -0.2    / ! note: exp was repl by multiplication
      IF (kpp_const_fields%L_SLAB) THEN
         data  difm0   / 0.000005  / ! max visc due to shear instability
         data  difs0   / 0.000005  / ! max diff ..  .. ..    ..
         data  difmiw  / 0.000001  / ! background/internal waves visc(m^2/s)
         data  difsiw  / 0.000001 / ! ..         ..       ..    diff(m^2/s)
      ELSE
         data difm0   / 0.005   /
         data difs0   / 0.005   /
         data difmiw  / 0.0001  /
         data difsiw  / 0.00001 /
      ENDIF
      data  difmcon / 0.0000   / ! max visc for convection  (m^2/s)
      data  difscon / 0.0000   / ! max diff for convection  (m^2/s)
      data  c1/ 1.0/
      data  c0/ 0.0/
      data  mRi/ 1 /            ! number of vertical smoothing passes
c     
c     compute interior gradient Ri at all interfaces, except surface
      
c-----------------------------------------------------------------------
c     compute interior gradient Ri at all interfaces ki=1,km, (not surface)
c     use visc(imt,ki=1,km) as temporary storage to be smoothed
c     use dift(imt,ki=1,km) as temporary storage of unsmoothed Ri
c     use difs(imt,ki=1,km) as dummy in smoothing call
      
      do 110 ki = 1, km
c     WRITE(6,*) ki
c     WRITE(6,*) 'Shsq(ki) = ',Shsq(ki)
         kpp_2d_fields%Rig(ki)  = kpp_2d_fields%dbloc(ki) * 
     +        (kpp_2d_fields%zm(ki)-kpp_2d_fields%zm(ki+1))/
     $        (kpp_2d_fields%Shsq(ki) + epsln)
         kpp_2d_fields%dift(ki) = kpp_2d_fields%Rig(ki)
         kpp_2d_fields%difm(ki) = kpp_2d_fields%dift(ki)                        
 110  continue
 
c-----------------------------------------------------------------------
c     vertically smooth Ri mRi times
      do j = 1,mRi
        call z121(kmp1,c0,Riinfty,kpp_2d_fields%difm,
     +        kpp_2d_fields%difs)
      enddo

c-----------------------------------------------------------------------
c                           after smoothing loop
      DO ki = 1, km
c  evaluate f of unsmooth Ri (fri) for convection        store in fcon
c  evaluate f of   smooth Ri (fri) for shear instability store in fri
 
         Rigg  = AMAX1( kpp_2d_fields%dift(ki) , Ricon )
         ratio = AMIN1( (Ricon-Rigg)/Ricon , c1 )
         fcon  = (c1 - ratio*ratio)
         fcon  = fcon * fcon * fcon
         
         Rigg  = AMAX1( kpp_2d_fields%difm(ki) , c0 )
         ratio = AMIN1( Rigg/Riinfty , c1 )
         fri   = (c1 - ratio*ratio)
         fri   = fri * fri * fri
         
c            if(i.eq.1)write(6,*)ki,dift(i,ki),visc(i,ki),Rigg,ratio,fri,
c     +           fcon
	     
c ************************   Overwrite with Gent's PP **********
c           fcon = 0.0
c           Rigg  = AMAX1( dift(i,ki) , c0 )
c           fri   = c1 / (c1 + 10. * Rigg )
c           difm0 = 0.1 * fri
c           difs0 = 0.1 * fri * fri

c  ************************   Overwrite with original PP
c           fcon = 0.0
c           Rigg  = AMAX1( dift(i,ki) , c0 )
c           fri   = c1 / (c1 +  5. * Rigg )
c           difm0 = 0.01 * fri
c           difs0 = (difmiw + fri * difm0)

c ----------------------------------------------------------------------
c            evaluate diffusivities and viscosity
c    mixing due to internal waves, and shear and static instability
 
         kpp_2d_fields%difm(ki) = 
     +        (difmiw + fcon * difmcon + fri * difm0)
         kpp_2d_fields%difs(ki) = 
     +        (difsiw + fcon * difscon + fri * difs0)
         kpp_2d_fields%dift(ki) = kpp_2d_fields%difs(ki)
      END DO

c ------------------------------------------------------------------------
c         set surface values to 0.0
 
      kpp_2d_fields%difm(0) = c0
      kpp_2d_fields%dift(0) = c0
      kpp_2d_fields%difs(0) = c0      

      return
      end

c *********************************************************************
      Subroutine z121 (kmp1,vlo,vhi,V,w)
      IMPLICIT NONE

c     Necessary for IMPLICIT NONE (NPK 11/2/13)
      INTEGER kmp1
      REAL vlo,vhi,tmp,wait
      INTEGER k,km
      
c    Apply 121 smoothing in k to 2-d array V(i,k=1,km)
c     top (0) value is used as a dummy
c     bottom (kmp1) value is set to input value from above.
 
c  input
c      include 'local_pt.com'

      real V(0:kmp1)  ! 2-D array to be smoothed in kmp1 direction
      real w(0:kmp1)  ! 2-D array of internal weights to be computed
 
      km  = kmp1 - 1
 
      w(0)    =   0.0         
      w(kmp1) =   0.0             
      V(0)    =   0.0      
      V(kmp1) =   0.0
      
      do k=1,km  
         if((V(k).lt.vlo).or.(V(k).gt.vhi)) then
            w(k) = 0.0
c     w(i,k) = 1.0
         else 
            w(k) = 1.0
         endif
         
      enddo
 
      do k=1,km
         tmp    = V(k)
         V(k) = w(k-1)*V(0)+2.*V(k)+w(k+1)*V(k+1)
         wait   = w(k-1) + 2.0 + w(k+1)
         V(k) = V(k) / wait             
         V(0) = tmp
      enddo
       
      return
      end
c     *********************************************************************

      subroutine ddmix (km, kmp1,
     +     alphaDT,betaDS,kpp_2d_fields)
      IMPLICIT NONE
c
c     Rrho dependent interior flux parameterization.
c     Add double-diffusion diffusivities to Ri-mix values at blending
c     interface and below.
c
 
c input    
! Automatically includes parameter.inc!
#include "kpp_3d_type.com"     
c      include 'local_pt.com'

c     Necessary for IMPLICIT NONE (NPK 11/2/13)
      integer km,kmp1,ki
      real dsfmax,rrho0
      
      real alphaDT(kmp1)  ! alpha * DT  across interfaces
      real betaDS(kmp1)   ! beta  * DS  across interfaces
c     real zgrid(kmp1)

      TYPE(kpp_2d_type) :: kpp_2d_fields      

c output
c      real visc(0:kmp1)  ! interior viscosity           (m^2/s)
c      real dift(0:kmp1)  ! interior thermal diffusivity (m^2/s)
c      real difs(0:kmp1)  ! interior scalar  diffusivity (m^2/s)
c
c local
      real Rrho              ! dd parameter
      real diffdd            ! double diffusion diffusivity scale
      real prandtl           ! prandtl number

c      save Rrho0,dsfmax

      data Rrho0  /  1.9   / ! Rp=(alpha*delT)/(beta*delS)
      data dsfmax / 1.0e-4 / ! .0001 m2/s


      DO ki= 1, km           

c     salt fingering case
         if((alphaDT(ki).gt.betaDS(ki)).and.(betaDS(ki).gt.0.)) then
            Rrho  = MIN(alphaDT(ki) / betaDS(ki) , Rrho0)
            diffdd     =         1.0-((Rrho-1)/(Rrho0-1))**2
            diffdd     = dsfmax*diffdd*diffdd*diffdd
            kpp_2d_fields%dift(ki) = kpp_2d_fields%dift(ki) + diffdd * 
     +           0.8 / Rrho
            kpp_2d_fields%difs(ki) = kpp_2d_fields%difs(ki) + diffdd
            
c     diffusive convection
         else if ((alphaDT(ki).lt.0.0).and.(betaDS(ki).lt.0.0).and.
     $           (alphaDT(ki).lt.betaDS(ki)) ) then
            Rrho    = alphaDT(ki) / betaDS(ki) 
            diffdd  = 1.5e-6*9.0*0.101*exp(4.6*exp(-0.54*(1/Rrho-1)))
            prandtl = 0.15*Rrho
            if (Rrho.gt.0.5) prandtl = (1.85-0.85/Rrho)*Rrho
            kpp_2d_fields%dift(ki) = kpp_2d_fields%dift(ki) + diffdd
            kpp_2d_fields%difs(ki) = kpp_2d_fields%difs(ki) + 
     +           prandtl*diffdd
            
         endif
      ENDDO

      return
      end

c *********************************************************************

      subroutine blmix 
     $     (km   , mdiff , 
     $     ustar, bfsfc, hbl , stable, caseA, 
     $     kbl   , 
     $     gat1 , dat1 , dkm1, blmc, 
     $     sigma, wm, ws, kpp_2d_fields,kpp_const_fields)
c     mixing coefficients within boundary layer depend on surface
c     forcing and the magnitude and gradient of interior mixing below
c     the boundary layer ("matching").

      IMPLICIT NONE

CAUTION if mixing bottoms out at hbl = -zgrid(km) THEN
c   fictious layer kmp1 is needed with small but finite width (eg. 1.e-10)
c model
c      include 'local_pt.com'

! Automatically includes parameter.inc!
#include "kpp_3d_type.com"

      TYPE(kpp_2d_type) :: kpp_2d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      
      integer km                !,kmp1        ! number of vertical levels
c     integer imt            ! number of horizontal grid points
      integer mdiff             ! number of viscosities + diffusivities
c     real zgrid(kmp1)   ! vertical grid (<=0)               (m)
c     real hwide(kmp1)   ! layer thicknesses                 (m)
c     
c     input
      real ustar                ! surface friction velocity         (m/s)
      real bfsfc                ! surface buoyancy forcing        (m^2/s^3)
      real hbl                  ! boundary layer depth              (m)
      real stable               ! = 1 in stable forcing
      real caseA                ! = 1 in case A
c      real visc(imt,0:kmp1)     ! vertical viscosity coefficient    (m^2/s)
c      real difs(imt,0:kmp1)     ! vertical scalar diffusivity       (m^2/s)
c      real dift(imt,0:kmp1)     ! vertical temperature diffusivity  (m^2/s)
      integer kbl               ! index of first grid level below hbl
c     
c output
      real gat1(mdiff)
      real dat1(mdiff)
      real dkm1(mdiff)          ! boundary layer difs at kbl-1 level
      real blmc(km,mdiff)       ! boundary layer mixing coefficients(m^2/s)
c      real ghats(km)            ! nonlocal scalar transport
c
c  local
      real sigma                ! normalized depth (d / hbl)
      real ws, wm               ! turbulent velocity scales         (m/s)
c  None of these were previously declared ... (NPK 6/2/13)
      real a1,a2,a3,am,as,c1,c2,c3,cg,cm,cs,cstar,delhat,difsh,
     +     difsp,difth,diftp,dvdzup,epsln,f1,gm,gs,
     +     dvdzdn,epsilon,gt,r,visch,viscp,zetam,sig,
     +     zetas
      integer ki,kn
 
c      save epsln,epsilon,c1,am,cm,c2,zetam,as,cs,c3,zetas,
c     $     cstar

      data epsln             /   1.e-20 /
      data epsilon           /   0.1    /
      data c1                /   5.0    /
      data am,cm,c2,zetam    /   1.257  ,  8.380, 16.0, - 0.2 / !7-24-92
      data as,cs,c3,zetas    / -28.86   , 98.96 , 16.0, - 1.0 /
      data cstar             /    5.    /
c     
      cg = cstar * kpp_const_fields%vonk * 
     +     (cs * kpp_const_fields%vonk * epsilon)**(1./3.)
      
c compute velocity scales at hbl
      sigma = stable * 1.0 + (1.-stable) * epsilon
      
c      WRITE(6,*) 'wscale(',sigma,hbl,ustar,bfsfc
      call wscale(sigma, hbl, ustar, bfsfc,wm,ws,kpp_const_fields)
      kn    = ifix(caseA+epsln) *(kbl -1) +
     $     (1-ifix(caseA+epsln)) * kbl
      
c     find the interior viscosities and derivatives at hbl(i) 
      delhat = 0.5*kpp_2d_fields%hm(kn)-kpp_2d_fields%zm(kn) - 
     +     hbl
      R      = 1.0 - delhat / kpp_2d_fields%hm(kn)
c      WRITE(6,*) 'kn = ',kn
c      WRITE(6,*) 'kpp_2d_fields%difm(kn-1) =',kpp_2d_fields%difm(kn-1)
c      WRITE(6,*) 'kpp_2d_fields%difm(kn) =',kpp_2d_fields%difm(kn)
c      WRITE(6,*) 'kpp_2d_fields%hm(kn) = ',kpp_2d_fields%hm(kn)
c      WRITE(6,*) 'kpp_2d_fields%hm(kn+1) =',kpp_2d_fields%hm(kn+1)

      dvdzup = (kpp_2d_fields%difm(kn-1) - kpp_2d_fields%difm(kn)) / 
     +     kpp_2d_fields%hm(kn) 
      dvdzdn = (kpp_2d_fields%difm(kn)   - kpp_2d_fields%difm(kn+1)) 
     +     / kpp_2d_fields%hm(kn+1)
      viscp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $     R  * (dvdzdn + abs(dvdzdn)) )
      
      dvdzup = (kpp_2d_fields%difs(kn-1) - kpp_2d_fields%difs(kn)) / 
     +     kpp_2d_fields%hm(kn) 
      dvdzdn = (kpp_2d_fields%difs(kn)   - kpp_2d_fields%difs(kn+1)) 
     +     / kpp_2d_fields%hm(kn+1)
      difsp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $     R  * (dvdzdn + abs(dvdzdn)) )
      
      dvdzup = (kpp_2d_fields%dift(kn-1) - kpp_2d_fields%dift(kn)) / 
     +     kpp_2d_fields%hm(kn) 
      dvdzdn = (kpp_2d_fields%dift(kn)   - kpp_2d_fields%dift(kn+1)) 
     +     / kpp_2d_fields%hm(kn+1)
      diftp  = 0.5 * ( (1.-R) * (dvdzup + abs(dvdzup))+
     $     R  * (dvdzdn + abs(dvdzdn)) )
c     
      visch  = kpp_2d_fields%difm(kn) + viscp * delhat
      difsh  = kpp_2d_fields%difs(kn) + difsp * delhat
      difth  = kpp_2d_fields%dift(kn) + diftp * delhat
c     
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
      
c     Turn off interior matching here
c     gat1(i,1) = 0.0001
c     gat1(i,2) = 0.00001
c     gat1(i,3) = 0.00001
c     do m=1,3
c       dat1(i,m) = 0.0
c       enddo

c     
      do 300 ki = 1,km       
c
c     compute turbulent velocity scales on the interfaces
c     
         sig     = (-kpp_2d_fields%zm(ki) + 0.5 * 
     +        kpp_2d_fields%hm(ki)) / hbl
         sigma   = stable*sig + (1.-stable)*AMIN1(sig,epsilon)
         call wscale(sigma, hbl, ustar, bfsfc,wm,ws,kpp_const_fields)
c
c     compute the dimensionless shape functions at the interfaces
c
         sig = (-kpp_2d_fields%zm(ki) + 0.5 * 
     +        kpp_2d_fields%hm(ki)) / hbl
         a1 = sig - 2.
         a2 = 3.-2.*sig
         a3 = sig - 1.
c     
         Gm = a1 + a2 * gat1(1) + a3 * dat1(1) 
         Gs = a1 + a2 * gat1(2) + a3 * dat1(2)
         Gt = a1 + a2 * gat1(3) + a3 * dat1(3)
c     
c     compute boundary layer diffusivities at the interfaces
c     
         blmc(ki,1) = hbl * wm * sig * (1. + sig * Gm)
         blmc(ki,2) = hbl * ws * sig * (1. + sig * Gs)
         blmc(ki,3) = hbl * ws * sig * (1. + sig * Gt)
c     
c     nonlocal transport term = ghats * <ws>o
         kpp_2d_fields%ghat(ki) = (1.-stable) * cg / (ws*hbl+epsln)
 300  continue
 
c find diffusivities at kbl-1 grid level 
      sig   =  -kpp_2d_fields%zm(kbl-1)  / hbl
      sigma =  stable * sig + (1.-stable) * AMIN1(sig,epsilon)
c
      call wscale(sigma, hbl, ustar, bfsfc,   wm, ws,kpp_const_fields)
c     
      sig = -kpp_2d_fields%zm(kbl-1) / hbl
      a1= sig - 2.
      a2 = 3.-2.*sig
      a3 = sig - 1.
      Gm = a1 + a2 * gat1(1) + a3 * dat1(1)
      Gs = a1 + a2 * gat1(2) + a3 * dat1(2)
      Gt = a1 + a2 * gat1(3) + a3 * dat1(3)
      dkm1(1) = hbl * wm * sig * (1. + sig * Gm)
      dkm1(2) = hbl * ws * sig * (1. + sig * Gs)
      dkm1(3) = hbl * ws * sig * (1. + sig * Gt)

      return
      end

c ******************************************************************

      subroutine enhance (km   ,  mdiff , dkm1,
     &                    hbl   , kbl  , caseA,
     &                    blmc, kpp_2d_fields, kpp_const_fields)
      IMPLICIT NONE
c
c enhance the diffusivity at the kbl-.5 interface
c
c
! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
c     input
c     include 'local_pt.com'

c     Necessary for IMPLICIT NONE (NPK 11/2/13)
      real dkmp5,dstar
      integer ki

      integer km                !,kmp1           ! number of vertical levels
c      integer imt               ! number of horizontal grid points
      integer mdiff             ! number of viscosities + diffusivities
      integer kbl               ! grid above hbl
      real hbl                  ! boundary layer depth             (m)
      real dkm1(mdiff)          ! bl diffusivity at kbl-1 grid level
c     real zgrid(kmp1)          ! vertical grid (<= 0)             (m)
c     real visc(imt,0:kmp1)     ! enhanced viscosity               (m^2/s) 
c     real difs(imt,0:kmp1)     ! enhanced thermal diffusivity     (m^2/s)
c     real dift(imt,0:kmp1)     ! enhanced scalar  diffusivity     (m^2/s)
      real caseA                ! = 1 in caseA, = 0 in case B
 
      TYPE(kpp_2d_type) :: kpp_2d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

c input/output
c     real ghats(km)            ! nonlocal transport               (s/m**2)
c                               ! modified ghats at kbl(i)-1 interface
c output
      real blmc(km,mdiff)       ! enhanced bound. layer mixing coeff.
c
c local
      real delta                ! fraction hbl lies beteen zgrid neighbors
c
      do ki=1,km-1         
         if(ki .eq. (kbl - 1) ) then            
            delta = (hbl+kpp_2d_fields%zm(ki)) / 
     +           (kpp_2d_fields%zm(ki)-kpp_2d_fields%zm(ki+1))
            
            dkmp5 = caseA * kpp_2d_fields%difm(ki) + (1.-caseA) * 
     +           blmc(ki,1)
            dstar = (1.-delta)**2 * dkm1(1) + delta**2 * dkmp5      
            blmc(ki,1) = (1.-delta) * kpp_2d_fields%difm(ki) + delta * 
     +           dstar

            dkmp5 = caseA * kpp_2d_fields%difs(ki) + (1.-caseA) * 
     +           blmc(ki,2)
            dstar = (1.-delta)**2 * dkm1(2) + delta**2 * dkmp5    
            blmc(ki,2) = (1.-delta) * kpp_2d_fields%difs(ki) + delta * 
     +           dstar

            dkmp5 = caseA * kpp_2d_fields%dift(ki) + (1.-caseA) * 
     +           blmc(ki,3)
            dstar = (1.-delta)**2 * dkm1(3) + delta**2 * dkmp5     
            blmc(ki,3) = (1.-delta) * kpp_2d_fields%dift(ki) + delta * 
     +           dstar
            
            kpp_2d_fields%ghat(ki) = (1.-caseA) * kpp_2d_fields%ghat(ki)            
         endif
      enddo

      return
      end

c***********************************************************************

      SUBROUTINE wm_ws_lookup(kpp_const_fields)

      IMPLICIT NONE
      
#include "kpp_3d_type.com"
      TYPE(kpp_const_type) :: kpp_const_fields

      real zmin,zmax,umin,umax,usta,zeta,zehat,epsln,
     +     am,cm,c1,c2,zetam,as,cs,c3,zetas,cstar,deltau,deltaz
      integer i,j,ni,nj
      parameter ( ni = 890,     ! number of values for zehat
     $     nj = 48)             ! number of values for ustar      
      data epsln             /   1.e-20 /
      data c1                /   5.0    /
      data zmin,zmax  / -4.e-7, 0.0   / ! m3/s3
      data umin,umax  /  0.   , .04   / ! m/s
      data am,cm,c2,zetam    /   1.257  ,  8.380, 16.0, - 0.2 / !7-24-92
      data as,cs,c3,zetas    / -28.86   , 98.96 , 16.0, - 1.0 /
      data cstar             /    5.    /

      deltaz = (zmax-zmin)/(ni+1) 
      deltau = (umax-umin)/(nj+1)
      
      do 100 i=0,ni+1
         zehat = deltaz*(i) + zmin
         do 90 j=0,nj+1
            usta = deltau*(j) + umin
            zeta = zehat/(usta**3+epsln)
            
            if(zehat.ge.0.) then
               kpp_const_fields%wmt(i,j) = kpp_const_fields%vonk*usta/
     +              (1.+c1*zeta)
               kpp_const_fields%wst(i,j) = kpp_const_fields%wmt(i,j)
            else
               if(zeta.gt.zetam) then
                  kpp_const_fields%wmt(i,j) = kpp_const_fields%vonk* 
     +                 usta * (1.-c2*zeta)**(1./4.)
               else
                  kpp_const_fields%wmt(i,j) = kpp_const_fields%vonk* 
     +                 (am*usta**3 - cm*zehat)**(1./3.)
               endif
               if(zeta.gt.zetas) then
                  kpp_const_fields%wst(i,j) = kpp_const_fields%vonk* 
     +                 usta * (1.-c3*zeta)**(1./2.)
               else
                  kpp_const_fields%wst(i,j) = kpp_const_fields%vonk* 
     +                 (as*usta**3 - cs*zehat)**(1./3.)
               endif
            endif   
 90      continue   
 100  continue
      
      END
