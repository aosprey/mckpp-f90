c Test test test
*********************************************************************
      subroutine vmix(kpp_2d_fields,kpp_const_fields,hmixn,kmixn)
c  Interface between 1-d model and vertical mixing

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include "kpp_3d_type.com"

c inputs including those from common.inc and parameter.inc
      type(kpp_2d_type) :: kpp_2d_fields
      type(kpp_const_type) :: kpp_const_fields
      real B0,B0sol,ustar

c outputs including those to common.inc
      real hmixn                ! boundary layer depth (m)
      integer kmixn       
      real rhob

c local
      real dVsq(nzp1)    ! (velocity shear re sfc)^2      (m/s)^2
      real Ritop(nz)     ! numerator of bulk Richardson Number (m/s)^2
c     Ritop = (-z - -zref)* delta buoyancy w/ respect to sfc(m/s^2)
      real alphaDT(nz)   ! alpha * DT across interfaces
      real betaDS(nz)    ! beta  * DS across interfaces
      real epsilon,epsln
      real alpha,beta,exppr
      real sigma,sigma0
      real cpsw
      real tau
      real zref,wz,bref
      integer k,n,kl
      real del
      real dlimit,vlimit
      integer jerl(12)
c          month  1   2   3   4   5   6   7   8   9   10  11  12
      data jerl / 2 , 2 , 2 , 3 , 3 , 3 , 4 , 4 , 4 , 4 , 3 , 2 /

      epsilon = 0.1
      epsln   = 1.e-20
c     LKPP = .true.
c                        ensure that bottom layer isn't zero thickness
c     Moved to init_env (init.f) for OPENMP compatability - NPK
c     kpp_const_fields%hm(nzp1) = AMAX1(kpp_const_fields%hm(nzp1),epsln)

c                        find the jerlov water type for this month 
c      if(kpp_2d_fields%jerlov.lt.1) then
c        write(nuerr,*) 'Make sure you have got time in days'
c        write(nuerr,*) 'jerlov = ',kpp_2d_fields%jerlov
c       swtime = time
c11     if(swtime.ge.dpy) then
c        swtime = swtime-dpy
c        goto 11
c       endif  
c       mon = int(dble(swtime/dpy)*12.) + 1 !(swtime/30.147) with dpy=365.
c       mon = MIN0(mon,12)
c      jwtype(ipt) = jerl(mon)
c      else
c         jwtype = kpp_2d_fields%jerlov
c      endif


c calculate density of fresh water and brine in surface layer
      alpha = 1.
      beta  = 1.
      exppr = 0.0
      sigma0=0
      sigma=0
c      WRITE(6,*) 'Before ABK80, sst = ',kpp_2d_fields%X(1,1)
      call ABK80(0.0,kpp_2d_fields%X(1,1),-kpp_2d_fields%zm(1),
     +     alpha,beta,exppr,sigma0,sigma)
c      WRITE(6,*) 'After ABK80, sst = ',kpp_2d_fields%X(1,1)
c      IF (beta .lt. 0) THEN
c         WRITE(6,*) 'sigma0=',sigma0,kpp_2d_fields%X(1,1),alpha,beta,
c     +        exppr,sigma
c      ENDIF
      kpp_2d_fields%rhoh2o = 1000. + sigma0
      call ABK80(kpp_const_fields%SICE,kpp_2d_fields%X(1,1),
     +     -kpp_2d_fields%zm(1),alpha,beta,exppr,sigma0,sigma)
      rhob      = 1000. + sigma0
 
c     calculate temperature and salt contributions of buoyancy gradients  
c     calculate buoyancy profile (m/s**2) on gridlevels

      do 10 k=1,nzp1
         call ABK80(kpp_2d_fields%X(k,2)+kpp_2d_fields%Sref,
     +        kpp_2d_fields%X(k,1),-kpp_2d_fields%zm(k),
     +        alpha,beta,exppr,sigma0,sigma)
         kpp_2d_fields%rho(k)= 1000. + sigma0
         kpp_2d_fields%CP(k) = CPSW(kpp_2d_fields%X(k,2)+
     +        kpp_2d_fields%Sref,kpp_2d_fields%X(k,1),
     +        -kpp_2d_fields%zm(k))
         kpp_2d_fields%talpha(k) = alpha
         kpp_2d_fields%sbeta(k)  = beta
         kpp_2d_fields%buoy(k) = -kpp_const_fields%grav * sigma0 / 1000.
 10   continue
 
      kpp_2d_fields%rho(0) = kpp_2d_fields%rho(1) 
      kpp_2d_fields%CP(0) = kpp_2d_fields%CP(1)   
      kpp_2d_fields%talpha(0) = kpp_2d_fields%talpha(1)   
      kpp_2d_fields%sbeta(0)  = kpp_2d_fields%sbeta(1)           

c Call to ntflx, put here to allow removal of diagnostic call to vmix
c and to ensure the most recent cp,rho used (consistent with other 
c surface fluxes?)
      call ntflx(kpp_2d_fields,kpp_const_fields)

c calculate kinematic surface momentum fluxes
      kpp_2d_fields%wU(0,1) = -kpp_2d_fields%sflux(1,5,0) / 
     +     kpp_2d_fields%rho(0)
      kpp_2d_fields%wU(0,2) = -kpp_2d_fields%sflux(2,5,0) / 
     +     kpp_2d_fields%rho(0)
      tau     = sqrt( kpp_2d_fields%sflux(1,5,0)**2 +
     +     kpp_2d_fields%sflux(2,5,0)**2 ) 
     &     +1.e-16
c  1.e-16 added to stop subsequent division by zero if tau=0.0
      ustar = sqrt( tau / kpp_2d_fields%rho(0) )
 
c total turbulent kinematic temperature flux (C m/s)
      kpp_2d_fields%wX(0,1)  = -kpp_2d_fields%sflux(4,5,0) / 
     +     kpp_2d_fields%rho(0) / kpp_2d_fields%CP(0)
 
c total turbulent kinematic salinity flux (o/oo m/s)
      kpp_2d_fields%wX(0,2) = kpp_2d_fields%Ssurf*
     +     kpp_2d_fields%sflux(6,5,0)/
     +     kpp_2d_fields%rhoh2o+(kpp_2d_fields%Ssurf-
     +     kpp_const_fields%SICE)*
     +     kpp_2d_fields%sflux(5,5,0)/rhob

c calculate total kinematic surface buoyancy flux (m**2/s**3)
      B0 = -kpp_const_fields%grav*(kpp_2d_fields%talpha(0)*
     +     kpp_2d_fields%wX(0,1) - kpp_2d_fields%sbeta(0)*
     +     kpp_2d_fields%wX(0,2) )
      kpp_2d_fields%wX(0,NSP1) =  - B0
      B0sol = kpp_const_fields%grav * kpp_2d_fields%talpha(0) * 
     +     kpp_2d_fields%sflux(3,5,0) / 
     +     (kpp_2d_fields%rho(0) * kpp_2d_fields%CP(0))
c     calculate temperature and salt contributions of buoyancy gradients
c               on interfaces for double diffusion      
      do 105 n = 1,nz
         alphaDT(n) =0.5 *(kpp_2d_fields%talpha(n)+
     +        kpp_2d_fields%talpha(n+1)) *
     +        (kpp_2d_fields%X(n,1) - kpp_2d_fields%X(n+1,1))
         betaDS(n)  =0.5 *(kpp_2d_fields%sbeta(n) + 
     +        kpp_2d_fields%sbeta(n+1)) *
     +        (kpp_2d_fields%X(n,2) - kpp_2d_fields%X(n+1,2))
 105  continue
      
c     compute buoyancy and shear profiles
      do 115  n = 1,nz
         zref =  epsilon * kpp_2d_fields%zm(n)
c     compute reference buoyancy and velocity
         wz    = AMAX1(kpp_2d_fields%zm(1),zref) 
         kpp_2d_fields%uref  = kpp_2d_fields%U(1,1) * wz / zref
         kpp_2d_fields%vref  = kpp_2d_fields%U(1,2) * wz / zref
         bref  = kpp_2d_fields%buoy(1)* wz / zref
         do 125 kl = 1,nz
            IF(zref.ge.kpp_2d_fields%zm(kl)) go to 126
            wz = AMIN1(kpp_2d_fields%zm(kl)-
     +           kpp_2d_fields%zm(kl+1),kpp_2d_fields%zm(kl)-zref) 
            del = 0.5 * wz / (kpp_2d_fields%zm(kl) - 
     +           kpp_2d_fields%zm(kl+1))
c            WRITE(6,*) wz, del
            kpp_2d_fields%uref = kpp_2d_fields%uref - 
     +           wz*( kpp_2d_fields%U(kl,1) + del *
     +           ( kpp_2d_fields%U(kl+1,1)- kpp_2d_fields%U(kl,1))) 
     +           /zref
            kpp_2d_fields%vref = kpp_2d_fields%vref - 
     +           wz*( kpp_2d_fields%U(kl,2) + del *
     +           ( kpp_2d_fields%U(kl+1,2)- kpp_2d_fields%U(kl,2))) 
     +           /zref
            bref=bref -wz*(kpp_2d_fields%buoy(kl) + 
     &           del *(kpp_2d_fields%buoy(kl+1)-kpp_2d_fields%buoy(kl))) 
     +           /zref
 125     continue
 126     continue

         Ritop(n) = (zref - kpp_2d_fields%zm(n)) * 
     +        (bref - kpp_2d_fields%buoy(n))
c     NPK Additions (25/9/2008). Prevent Ritop from going negative.
c     IF (Ritop(ipt,n) .lt. 0) Ritop(ipt,n) = epsln
         kpp_2d_fields%dbloc(n) = kpp_2d_fields%buoy(n) - 
     +        kpp_2d_fields%buoy(n+1)
         dVsq(n)  = (kpp_2d_fields%Uref - kpp_2d_fields%U(n,1))**2 + 
     +        (kpp_2d_fields%Vref - kpp_2d_fields%U(n,2))**2
         kpp_2d_fields%shsq(n)  = (kpp_2d_fields%U(n,1)-
     +        kpp_2d_fields%U(n+1,1))**2 + 
     +        (kpp_2d_fields%U(n,2)-kpp_2d_fields%U(n+1,2))**2         
 115  continue

      call kppmix(    
     $     dVsq,
     $     ustar , B0    , B0sol , alphaDT, betaDS ,
     $     Ritop,
     $     hmixn, kmixn, 
     $     kpp_2d_fields,kpp_const_fields)
      
      
c      limit the bottom diffusity and viscosity
c      zero diffusivities for no bottom flux option
c      if(LNBFLX) then 
c        dlimit = 0.0
c        vlimit = 0.0
c        do n=1,nsclr
c        wxNT(ipt,nz,n) = 0.0
c        enddo
c      else

      IF (kpp_const_fields%L_SLAB) THEN
c     Set very small background diffusivity values for slab ocean
c     Should we set wxNT(NZ,:) to zero as well for the slab? This is
c     done above for LNBFLX.
         dlimit=0.0
         vlimit=0.0
         DO n=1,nsclr
            kpp_2d_fields%wXNT(NZ,n)=0.0
         ENDDO
      ELSE
         dlimit = 0.00001
         vlimit = 0.0001
      ENDIF
      do k=nz,nzp1
         kpp_2d_fields%difm(k) = vlimit
         kpp_2d_fields%difs(k) = dlimit
         kpp_2d_fields%dift(k) = dlimit
      enddo
      kpp_2d_fields%ghat(nz) = 0.0
      
      return
      end
**********************************************************************

