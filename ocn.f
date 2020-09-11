      SUBROUTINE  ocnstep (kpp_2d_fields,kpp_const_fields)
c-----------------------------------------------------------------------
c Note in this version:
c   -  ADVECTIVE CORRECTIONS AVAILABLE
c   - read forcing, but compute SW from okta model   (fread in atmrad.f)
c   - use okta model for PAPA          (use cloudpapa in solar in bio.f)
c   - Jerlov water type II                            (SWDK in fluxes.f)
c   - albedo for ocean is set to 0.06, which is used for QSW from fcomp
c        or fread when rad/conv is not running:
c        albocn=0.06                             (init cnsts in input.f)
c   - no net fresh water flux into ocean when forcing with state
c        variables (laflx >= 1):
c        sflux(7,2,jptr) = - sflux(6,2,jptr)        (atmflx in fluxes.f)
c   - use psnow flux data to input observed SST,
c                                                    (fread in atmrad.f)
c  ALSO :
c         General ability to read large scale forcing from euc and doc LES
c ----------------------------------------------------------------------
c Originally ~/KPP/LARGE/ocn.f (from Bill Large)
c  Modified by SJW to remove those bits which a specific to Bill's
c work and just leave the bit we need
c Started 18/03/02
c-----------------------------------------------------------------------

c     Main driver for ocean module.
c     Integration is performed only on the permanent grid
c     Written   3 Mar 1991 - WGL
c     Modified  5 Jun 1992 - jan : implicit scheme
c              16 Nov      - jan : latest version
c              16 Nov 1994 - wgl : new KPP codes no temporary grid


      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
#include "ocn_energy.com"

c Input/Output
c      real U(NPTS,NZP1,NVEL), X(NPTS,NZP1,NSCLR), Rig(npts,nz),
c     +     dbloc(npts,nz), shsq(npts,nz)
c      REAL U_oned(NZP1,NVEL), X_oned(NZP1,NSCLR), Rig_oned(NZP1),
c     +     dbloc_oned(NZ),shsq_oned(NZP1)

      TYPE(kpp_2d_type) :: kpp_2d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
c
c Local Common Blocks
c
c      real hmixd_oned(0:1),     ! storage arrays for extrapolations
c     +     Us_oned(NZP1,NVEL ,0:1), ! ..      ..     ..  ..
c     +     Xs_oned(NZP1,NSCLR,0:1)  ! ..      ..     ..  ..
c      integer old,new           ! extrapolation index for Us,Xs,hmixd
c      common/ saveUXh /
c     +     old,new,Us,Xs,hmixd
c
c Local
c
c      real Un(NZP1,NVEL),       ! new profiles
c     +     Xn(NZP1,NSCLR),      ! ..  ..
c     +     hmixe,               ! estimated hmix (integration input )
c     +     hmixn,               ! new comp. hmix (    ..      output)
c     +     tol                  ! tolerance in hmix iteration
      real hmixe,hmixn,tol
      real Uo(NZP1,NVEL),
     +     Xo(NZP1,NSCLR)
      real Ux(NZP1,NVEL),       ! Additional variables to provide
     +     Xx(NZP1,NSCLR)       ! smoothing in the iteration.
      real Ui                   ! Ui used in damping (LH 8/08/2013)
      real dampU(NVEL)          ! dampU used to flag which Ui chosen in damping (LH 29/08/2013)
      real lambda               ! Factor to control smoothing
      integer
     +     iter,iconv                ! number of iterations
c      integer kmixe(npts),kmixn(npts)
      integer kmixe,kmixn
c
c More Local Variables (to make implicit none)
c
      real deltaz,rhonot,a,b, dzb(NZ), dm(0:NZ), hm(NZP1), zm(NZP1)
      integer k,l,n
c Number of iterations for computational instability
      integer comp_iter_max
      real rmsd(4),rmsd_threshold(4)
      data comp_iter_max /10/

      data lambda /0.5/

c Initialise grid
      kpp_2d_fields%zm = kpp_const_fields%zm
      kpp_2d_fields%dm = kpp_const_fields%dm
      kpp_2d_fields%hm = kpp_const_fields%hm
      kpp_2d_fields%tri = kpp_const_fields%tri

c Critical depth-integrated RMS difference between old and new profiles
c for repeating integration, for (/U,V,T,S/). Typical (stable)
c values are O(10^-2) for U and V, O(10^-3) for T and O(10^-4) for S.
c NPK 17/5/13
      IF (kpp_const_fields%L_SLAB) THEN
         rmsd_threshold = (/4,4,100,100/)
      ELSE
         rmsd_threshold = (/4,4,100,100/)
      ENDIF      
c Change slab depth for Columbia ITCZ experiments within specified bounds
      IF (kpp_const_fields%L_SLAB .and.
     +     kpp_const_fields%L_COLUMBIA_LAND) THEN
         IF (kpp_2d_fields%dlat .ge. -30 .and.
     +        kpp_2d_fields%dlat .le. 30 .and.
     +        kpp_2d_fields%dlon .ge. 0 .and.
     +        kpp_2d_fields%dlon .le. 45) THEN
            kpp_2d_fields%dm(1)=0.1
            kpp_2d_fields%zm(1)=-0.05
            kpp_2d_fields%hm(1)=0.1
            rmsd_threshold = (/4,4,100,100/)
c Force recomputation of tridiagonal matrix coefficients
            do k=1,NZ
               dzb(k) = kpp_2d_fields%zm(k) - kpp_2d_fields%zm(k+1)
            enddo
            kpp_2d_fields%tri(0,1,1) = kpp_const_fields%dto/
     +           kpp_2d_fields%hm(1)
            kpp_2d_fields%tri(1,1,1) = kpp_const_fields%dto/
     +           kpp_2d_fields%hm(1)/dzb(1)
            do k=2,NZ
               kpp_2d_fields%tri(k,1,1) = kpp_const_fields%dto/
     +              kpp_2d_fields%hm(k)/dzb(k)
               kpp_2d_fields%tri(k,0,1) = kpp_const_fields%dto/
     +              kpp_2d_fields%hm(k)/dzb(k-1)
            enddo
         ENDIF
      ENDIF
         
      Uo=kpp_2d_fields%U(:,:)
      Xo=kpp_2d_fields%X(:,:)
      kpp_2d_fields%comp_flag=.TRUE.
      kpp_2d_fields%reset_flag=0.0

c      WRITE(6,*) 'Beginning of timestep:'
c      WRITE(6,*) 'U=',kpp_2d_fields%U(:,1)
c      WRITE(6,*) 'V=',kpp_2d_fields%U(:,2)
c      WRITE(6,*) 'T=',kpp_2d_fields%X(:,1)
c      WRITE(6,*) 'S=',kpp_2d_fields%X(:,2)

      DO WHILE (kpp_2d_fields%comp_flag .and.
     +     kpp_2d_fields%reset_flag .le. comp_iter_max)
c     Estimate new profiles by  extrapolation
         do 20 k=1,NZP1
            do 22 l=1,NVEL
               kpp_2d_fields%U(k,l)=2.*
     +              kpp_2d_fields%Us(k,l,kpp_2d_fields%new)-
     +              kpp_2d_fields%Us(k,l,kpp_2d_fields%old)
               Ux(k,l)=kpp_2d_fields%U(k,l)
 22         continue
            do 24 l=1,NSCLR
               kpp_2d_fields%X(k,l)=2.*
     +              kpp_2d_fields%Xs(k,l,kpp_2d_fields%new)-
     +              kpp_2d_fields%Xs(k,l,kpp_2d_fields%old)
               Xx(k,l)=kpp_2d_fields%X(k,l)
 24         continue
 20      continue

c     Iteration loop for semi-implicit integration
c     Reset iteration counter
         iter=0
         iconv=0
c     This loop controls the number of compulsory iterations
c     The original value was 2, using (the upper value of the loop +2)
c     added by SJW (17 Jan 03) to try an alleviate some non-convergences
         DO iter=0,2
            DO k=1,NZP1
               DO l=1,NVEL
                  kpp_2d_fields%U(k,l)=lambda*Ux(k,l)+(1-lambda)*
     +                 kpp_2d_fields%U(k,l)
                  Ux(k,l)=kpp_2d_fields%U(k,l)
               ENDDO
               DO l=1,NSCLR
                  kpp_2d_fields%X(k,l)=lambda*Xx(k,l)+(1-lambda)*
     +                 kpp_2d_fields%X(k,l)
                  Xx(k,l)=kpp_2d_fields%X(k,l)
               ENDDO
            ENDDO
             call vmix(kpp_2d_fields,kpp_const_fields,hmixe,kmixe)
c     Overwrite mixing depth when using slab ocean
            IF (kpp_const_fields%L_SLAB) THEN
               IF (kpp_const_fields%L_COLUMBIA_LAND) THEN
                  IF (kpp_2d_fields%dlat .ge. -30 .and.
     +                 kpp_2d_fields%dlat .le. 30 .and.
     +                 kpp_2d_fields%dlon .ge. 0 .and.
     +                 kpp_2d_fields%dlon .le. 45) THEN
c                     WRITE(6,*) 'Change slab depth at ',
c     +                    kpp_2d_fields%dlat,kpp_2d_fields%dlon
                     hmixe=0.1
                  ELSE
                     hmixe=kpp_const_fields%slab_depth
                  ENDIF
               ELSE
                  hmixe=kpp_const_fields%slab_depth
               ENDIF
               kmixe=1
            ENDIF
            call ocnint(kpp_2d_fields,kpp_const_fields,1,kmixe,Uo,Xo)
            IF (kpp_const_fields%L_SLAB) kmixe=1
         ENDDO
c     The original code can be restored by reseting iter=1 and removing the
c     above loop
c     iter=1

         IF (kpp_const_fields%LKPP) THEN
 45         continue
            DO k=1,NZP1
               DO l=1,NVEL
                  kpp_2d_fields%U(k,l)=lambda*Ux(k,l)+(1-lambda)*
     +                 kpp_2d_fields%U(k,l)
                  Ux(k,l)=kpp_2d_fields%U(k,l)
               ENDDO
               DO l=1,NSCLR
                  kpp_2d_fields%X(k,l)=lambda*Xx(k,l)+(1-lambda)*
     +                 kpp_2d_fields%X(k,l)
                  Xx(k,l)=kpp_2d_fields%X(k,l)
               ENDDO
            ENDDO
            call vmix(kpp_2d_fields,kpp_const_fields,hmixn,kmixn)
c     Overwrite mixing depth when using slab ocean
            IF (kpp_const_fields%L_SLAB) THEN
               IF (kpp_const_fields%L_COLUMBIA_LAND) THEN
                  IF (kpp_2d_fields%dlat .ge. -30 .and.
     +                 kpp_2d_fields%dlat .le. 30 .and.
     +                 kpp_2d_fields%dlon .ge. 0 .and.
     +                 kpp_2d_fields%dlon .le. 45) THEN
c     WRITE(6,*) 'Change slab depth at ',
c     +                    kpp_2d_fields%dlat,kpp_2d_fields%dlon
                     hmixn=0.1
                  ELSE
                     hmixn=kpp_const_fields%slab_depth
                  ENDIF
               ELSE
                  hmixn=kpp_const_fields%slab_depth
               ENDIF
               kmixn=1
            ENDIF

            call ocnint(kpp_2d_fields,kpp_const_fields,1,kmixn,Uo,Xo)
            IF (kpp_const_fields%L_SLAB) kmixn=1
            iter = iter + 1

c     check iteration for convergence
            tol = hmixtolfrac*kpp_2d_fields%hm(kmixn)
            if(kmixn.eq.NZP1) tol = hmixtolfrac*kpp_2d_fields%hm(NZ)
c     write(40,*) abs(hmixn(ipt)-hmixe(ipt))/tol
            if(abs(hmixn-hmixe).gt.tol)  then
c     Uncommeting the following the lines iconv=0 to IF (iconv ...)
c     will make the model do two consecutive tests for convergence of the
c     hmix (added by SJW 17 Jan 03). This did not work well in testing for
c     long timestep, high resolution (the model generally failed to satisfy the
c     convergence test on two consecutive iterations.
               iconv=0
            ELSE
               iconv=iconv+1
            ENDIF
c     write(40,*) iconv
c     write(40,*) ntime,iter,iconv,hmixe,hmixn,
c     &        abs(hmixn(ipt)-hmixe(ipt))/tol
c     write(40+ntime,*) iter,iconv,hmixn,abs(hmixn(ipt)-hmixe(ipt))/tol
            IF (iconv .lt. 3) THEN
               if (iter.lt.itermax) then
                  hmixe = hmixn
                  kmixe = kmixn
                  goto 45
               else
c     use shallower hmix
                  if(hmixn.gt.hmixe) then
                     hmixe = hmixn ! comment out for hmix data
                     kmixe = kmixn ! ..      ..  ..  hmix data
                     goto 45    ! ..      ..  ..  hmix data
                  endif
               endif
            endif
!            if( iter.gt.(itermax+1) ) then
!               write(nuout,1009) kpp_const_fields%ntime, ! comment out for hmix data
!     +              kpp_2d_fields%dlon,kpp_2d_fields%dlat,
!     +              hmixe,hmixn,
!     +              hmixn-hmixe,kmixn,iter
! 1009          format('  long iteration at',i6,' steps',/,
!     +              ' location=(',f7.2,',',f6.2,')',/,
!     +              '  hmixest=',f7.2,' hmixnew=',f7.2,' diff=',f6.1,
!     +              ' kmixn=',i3,' iteration=',i3)
!            endif
         ENDIF
c     Trap for profiles that are very different from original profile
c     or clearly erroneous, to detect rare instances of instability
c     in the semi-implicit integration.  Reset to original profile,
c     add some noise via changing Coriolis term slightly, and try
c     integration again.
c     NPK 16/5/2013
         kpp_2d_fields%comp_flag=.FALSE.
      !WRITE(6,*) 'End of timestep:'
      !WRITE(6,*) 'U=',kpp_2d_fields%U(:,1)
      !WRITE(6,*) 'V=',kpp_2d_fields%U(:,2)
      !WRITE(6,*) 'T=',kpp_2d_fields%X(:,1)
      !WRITE(6,*) 'S=',kpp_2d_fields%X(:,2)
         DO k=1,NZ
            ! For slab ocean, do not test vertical gradient of temperature
            ! as this is likely to be very large due to the deep ocean below.
            IF (kpp_const_fields%L_SLAB) THEN
              IF (ABS(kpp_2d_fields%U(k,1)).ge. 10 .or.
     +           ABS(kpp_2d_fields%U(k,2)).ge.10) THEN
                kpp_2d_fields%comp_flag=.TRUE.
                kpp_2d_fields%f=kpp_2d_fields%f*
     +	           (1.01+MOD(kpp_2d_fields%reset_flag,2.0)*(-0.02))
              ENDIF
            ELSE
              IF (ABS(kpp_2d_fields%U(k,1)).ge. 10 .or.
     +            ABS(kpp_2d_fields%U(k,2)).ge.10 .or.
     +            ABS(kpp_2d_fields%X(k,1)-kpp_2d_fields%X(k+1,1))
     +              .ge. 10) THEN
                 WRITE(6,*) k,kpp_2d_fields%U(k,1),kpp_2d_fields%U(k,2),
     +                kpp_2d_fields%X(k,1),kpp_2d_fields%X(k+1,1)
                kpp_2d_fields%comp_flag=.TRUE.
                kpp_2d_fields%f=kpp_2d_fields%f*
     +			       (1.01+MOD(kpp_2d_fields%reset_flag,2.0)*(-0.02))
              ENDIF
            ENDIF
         END DO
         IF (.NOT. kpp_2d_fields%comp_flag) THEN
            rmsd(:)=0.
            DO k=1,NZP1
               rmsd(1)=rmsd(1)+(kpp_2d_fields%U(k,1)-Uo(k,1))*
     +              (kpp_2d_fields%U(k,1)-Uo(k,1))*
     +              kpp_2d_fields%hm(k)/kpp_2d_fields%dm(NZ)
               rmsd(2)=rmsd(2)+(kpp_2d_fields%U(k,2)-Uo(k,2))*
     +              (kpp_2d_fields%U(k,2)-Uo(k,2))*
     +              kpp_2d_fields%hm(k)/kpp_2d_fields%dm(NZ)
               rmsd(3)=rmsd(3)+(kpp_2d_fields%X(k,1)-Xo(k,1))*
     +              (kpp_2d_fields%X(k,1)-Xo(k,1))*
     +              kpp_2d_fields%hm(k)/kpp_2d_fields%dm(NZ)
               rmsd(4)=rmsd(4)+(kpp_2d_fields%X(k,2)-Xo(k,2))*
     +              (kpp_2d_fields%X(k,2)-Xo(k,2))*
     +              kpp_2d_fields%hm(k)/kpp_2d_fields%dm(NZ)
            ENDDO
            DO k=1,4
               rmsd(k)=SQRT(rmsd(k))
               IF (rmsd(k).ge.rmsd_threshold(k)) THEN
                  !WRITE(6,*) k,rmsd(k),rmsd_threshold(k)
                  kpp_2d_fields%comp_flag=.TRUE.
                  kpp_2d_fields%f=kpp_2d_fields%f*1.01
               ENDIF
            ENDDO
         ENDIF
         kpp_2d_fields%reset_flag=kpp_2d_fields%reset_flag+1
         IF (kpp_2d_fields%reset_flag .gt. comp_iter_max) THEN
            WRITE(6,*) 'Failed to find a reasonable solution '//
     +           'in the semi-implicit integration after ',
     +           comp_iter_max,' iterations.'
            WRITE(6,*) 'At point lat = ',
     +           kpp_2d_fields%dlat,' lon =',kpp_2d_fields%dlon,': '
            !WRITE(6,*) 'U = ',kpp_2d_fields%U(:,1)
            !WRITE(6,*) 'V = ',kpp_2d_fields%U(:,2)
            !WRITE(6,*) 'T = ',kpp_2d_fields%X(:,1)
            !WRITE(6,*) 'S = ',kpp_2d_fields%X(:,2)
            !WRITE(6,*) 'hmix = ',hmixn,kmixn
         ENDIF
      ENDDO
c     End of trapping code.

c     Output  Results from permanent grid iterations to common.inc
c     Compute diagnostic fluxes for writing to dat file
      do k=1,NZ
         deltaz = 0.5*(kpp_2d_fields%hm(k)+kpp_2d_fields%hm(k+1))
         do n=1,NSCLR
            kpp_2d_fields%wX(k,n)=-kpp_2d_fields%difs(k)*
     +           ((kpp_2d_fields%X(k,n)-kpp_2d_fields%X(k+1,n))/deltaz-
     +           kpp_2d_fields%ghat(k)*kpp_2d_fields%wX(0,n))
         enddo
         if(kpp_const_fields%LDD)
     +        kpp_2d_fields%wX(k,1)= -kpp_2d_fields%dift(k)*
     +        ((kpp_2d_fields%X(k,1)-
     +        kpp_2d_fields%X(k+1,1))/deltaz-kpp_2d_fields%ghat(k)*
     +        kpp_2d_fields%wX(0,1))
         kpp_2d_fields%wX(k,nsp1)= kpp_const_fields%grav *
     +        (kpp_2d_fields%talpha(k)*kpp_2d_fields%wX(k,1) -
     +        kpp_2d_fields%sbeta(k) * kpp_2d_fields%wX(k,2))
         do n=1,NVEL
            kpp_2d_fields%wU(k,n)= -kpp_2d_fields%difm(k)*
     +           (kpp_2d_fields%U(k,n)-kpp_2d_fields%U(k+1,n))/deltaz
         enddo
      enddo

c     Compute energetics
      rhonot = 1026.
      Eflx = 0.5 * ( (Uo(1,1) + kpp_2d_fields%U(1,1)) *
     +     kpp_2d_fields%sflux(1,5,0) +
     &     (Uo(1,2) + kpp_2d_fields%U(1,2)) *
     +     kpp_2d_fields%sflux(2,5,0) )
      Esnk =-0.5*rhonot* ( (Uo(NZ,1) + kpp_2d_fields%U(NZ,1)) *
     &     kpp_2d_fields%wU(NZ,1) +
     &     (Uo(NZ,2) + kpp_2d_fields%U(NZ,2)) *
     &     kpp_2d_fields%wU(NZ,2) )
            Ptke = 0.0
c     use "amax1" to prevent "underflow" in single precision
            do 120 k=1,NZ-1
               Ptke = Ptke - 0.5*( amax1(kpp_2d_fields%wU(k,1),1.E-10)*
     &              (rhonot   * (Uo(k,1) + kpp_2d_fields%U(k,1)) -
     &              rhonot   * (Uo(k+1,1) + kpp_2d_fields%U(k+1,1)) ) +
     &              amax1(kpp_2d_fields%wU(k,2),1.E-10)*
     &              (rhonot   * (Uo(k,2) + kpp_2d_fields%U(k,2)) -
     &              rhonot   * (Uo(k+1,2) + kpp_2d_fields%U(k+1,2)) ) )
 120        continue
            Tmke = 0.0
            do 130 k=1,NZP1
               rmke(k) = 0.5 * rhonot * (kpp_2d_fields%U(k,1)**2 +
     +              kpp_2d_fields%U(k,2)**2) * kpp_2d_fields%hm(k)
               Tmke = Tmke + rmke(k)
 130        continue

c     check heat and salt budgets
c     call budget(Xo,kpp_2d_fields,kpp_const_fields)

c     Set new profiles
c            do k=1,NZP1         ! values at NZP1 only change for slab ocean
c               do n=1,NVEL
c                  kpp_2d_fields%U(k,n) = kpp_2d_fields%U(k,n)
c               enddo
c               do n=1,NSCLR
c                  kpp_2d_fields%X(k,n) = kpp_2d_fields%X(k,n)
c               enddo
c            enddo
c     Set correct surface values, and CP and rho profiles for new profiles
c     Get latest profiles
c     Removed to ensure that hmix,diff,cp,rho diagnostics are consistent with
c     those used in the final iteration of the timestep
c     (CP,RHO are not quite correct for updated values, but are the values
c     by the integration) (SJW 16/01/03)
c     Need to consider improving convergence test!!! (SJW 16/01/03)

c     The final call to vmix is removed to ensure that the diffusion and
c     boundary layer profiles in the diagnostics are the ones used to calculate
c     the fluxes, as it stands at the moment this means that the CP and rho are
c     also the values used in the timestepping not the values appropriate to the
c     S,T at the new time level.
c            call vmix(kpp_2d_fields%U,kpp_2d_fields%X,hmixe,kmixe)
c     write(40,*) time,iter, hmixn,hmixe,kmixn,kmixe
            kpp_2d_fields%hmix = hmixn
            kpp_2d_fields%kmix = kmixn
            kpp_2d_fields%uref = kpp_2d_fields%U(1,1)
            kpp_2d_fields%vref = kpp_2d_fields%U(1,2)
            kpp_2d_fields%Tref = kpp_2d_fields%X(1,1)
            IF (kpp_const_fields%L_SSref) THEN
               kpp_2d_fields%Ssurf=kpp_2d_fields%SSref
            ELSE
               kpp_2d_fields%Ssurf=kpp_2d_fields%X(1,2)+
     +              kpp_2d_fields%Sref
            ENDIF

c     Damping currents, Added LH (06/08/2013)

            IF (kpp_const_fields%L_DAMP_CURR) THEN
               dampu(:)=0.
               do k=1,NZP1
                  do l=1,NVEL
                     a=0.99*ABS(kpp_2d_fields%U(k,l))
                     b=kpp_2d_fields%U(k,l)**2/
     +                    (kpp_const_fields%dt_uvdamp*
     +                    (86400./kpp_const_fields%dto))
                     Ui=MIN(a,b)
c     LH (29/08/2013) Add Flags to check which Ui (a or b) is chosen,
c     dtuvdamp=360 (specified in namelist).
c     The flags for u and v can be requested as diagnostics dampu_flag,
c     dampv_flag (singout 11,12). Note that the value of the flag is equal to
c     the *fraction* of levels at that point where (U**2)/r .lt. alpha*ABS(U),
c     1.0=all Ui are (U**2)/r
                     IF (b .lt. a) THEN
                        dampU(l)=dampU(l)+1.0/REAL(NZP1)
                     ENDIF

c     Apply damping
                     kpp_2d_fields%U(k,l)= kpp_2d_fields%U(k,l) -
     +                    SIGN(Ui,kpp_2d_fields%U(k,l))
                  enddo
               enddo
               kpp_2d_fields%dampu_flag=dampU(1)
               kpp_2d_fields%dampv_flag=dampU(2)
            ENDIF

c     End of damping

c     Save variables for next timestep
            kpp_2d_fields%old = kpp_2d_fields%new
            kpp_2d_fields%new = 1 - kpp_2d_fields%old
            kpp_2d_fields%hmixd(kpp_2d_fields%new) =
     +           kpp_2d_fields%hmix
            do k=1,NZP1
               do l=1,NVEL
                  kpp_2d_fields%Us(k,l,kpp_2d_fields%new)=
     +                 kpp_2d_fields%U(k,l)
               enddo
               do l=1,NSCLR
                  kpp_2d_fields%Xs(k,l,kpp_2d_fields%new)=
     +                 kpp_2d_fields%X(k,l)
               enddo
            enddo
c     close(40+ntime)
            return
            end

***********************************************************************

      SUBROUTINE ocnint(kpp_2d_fields,kpp_const_fields,intri,kmixe,
     +     Uo,Xo)

c     Integrate the ocn model by backwards Euler(implicit)discretization
c     On input : Un,Xn are estimated profiles which are used
c                to estimate diffusivity profiles at new time.
c              : Updated diffusivities from Un Xn are in common
c     On output: Un,Xn are new profiles after integration.

c     Written  19 March 1991 - jan

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
c Input
      integer intri             ! index for tri.diag. coeff
c     +     nzi
c      real  Uo(nzi+1,NVEL), Xo(nzi+1,NSCLR), ! old profiles
c     +      z(nzi+1),h(nzi+1),d(0:nzi)
      REAL Uo(NZP1,NVEL),Xo(NZP1,NSCLR)
c Output
c      real  Un(nzi+1,NVEL), Xn(nzi+1,NSCLR)  ! new profiles
      TYPE(kpp_2d_type) :: kpp_2d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

c Common tridiagonal matrix factors (set in "init_ocn")
c      real tri(0:NZtmax,0:1,NGRID)    ! dt/dz/dz factors in trid. matrix
c      common/ trifac / tri
c Local
      real cu (NZtmax), ! upper coeff for (k-1) on k line of trid.matrix
     +     cc (NZtmax), ! central ...     (k  ) ..
     +     cl (NZtmax), ! lower .....     (k-1) ..
     +     rhs(NZtmax)  ! right-hand-side terms
      real diff(0:NZtmax),gcap(NZtmax),ntflx(0:NZtmax,NSCLR),
     +     weight_hek,ekvel_max,dist_hek
c More local variables to make implicit none
      integer kmixe,i,npd,imode,n,k
      real ftemp,ghatflux,sturflux
      integer adv_mode
      real adv_mag
c
c     Small value to use when applying flux corrections at depth:
c     Amount to add to temperatures if flux corrections change
c     sign of vertical temperature gradient.
c
c      real prev_t

c      COMMON /save_fcorr_withz/ fcorr_withz, tinc_fcorr
c      COMMON /save_fcorr_twod/ fcorr_twod
c      COMMON /save_sal/ sal_clim


c      WRITE(nuout,*) 'Entering ocnint:'
c      WRITE(nuout,*) '--' L_FCORR =',L_FCORR
c      WRITE(nuout,*) '--> L_FCORR_WITHZ =',L_FCORR_WITHZ
c      WRITE(nuout,*) '--> L_RELAX_SST =',L_RELAX_SST
c      WRITE(nuout,*) '--> L_RELAX_CALCONLY =',L_RELAX_CALCONLY

c ********************************************************************
c U and V solution of tridiagonal matrix
c                set f = 0 for equatorial application
      ftemp = kpp_2d_fields%f
c      ftemp = 2. * (twopi/86164.) * sin(0.5*twopi/360.)
c      ftemp = 2. * (twopi/86164.) * sin(2.5*twopi/360.)
c      ftemp=0.
c                               set coefficients of tridiagonal matrix
      DO k=0,NZtmax
         diff(k)=kpp_2d_fields%difm(k)
      ENDDO
      call tridcof(diff,NZ,intri,cu,cc,cl,kpp_2d_fields)
c                               U right hand side and solution
      rhs(1)= Uo(1,1) + kpp_const_fields%dto*
     +     ( ftemp*.5*(Uo(1,2)+kpp_2d_fields%U(1,2)) -
     +     kpp_2d_fields%wU(0,1)/kpp_2d_fields%hm(1))
c      rhs(1)=Uo(1,1)+kpp_const_fields%dto*(ftemp*Uo(1,2) -
c     +     kpp_2d_fields%wU(0,1)/kpp_2d_fields%hm(1))
      do i=2,NZ-1
         rhs(i)= Uo(i,1) + kpp_const_fields%dto*ftemp*.5*(Uo(i,2)+
     +        kpp_2d_fields%U(i,2))
c         rhs(i)=Uo(i,1)+kpp_const_fields%dto*ftemp*Uo(i,2)
      enddo
      i=NZ                      ! bottom
      rhs(i)= Uo(i,1) + kpp_const_fields%dto*ftemp*.5*(Uo(i,2)+
     +     kpp_2d_fields%U(i,2)) + kpp_2d_fields%tri(i,1,intri)*
     +     kpp_2d_fields%difm(i)*Uo(i+1,1)
c      rhs(i)=Uo(i,1)+kpp_const_fields%dto*ftemp*Uo(i,2)+
c     +     kpp_2d_fields%tri(i,1,intri)*
c     +     kpp_2d_fields%difm(i)*Uo(i+1,1)

      call tridmat(cu,cc,cl,rhs,Uo(:,1),NZ,kpp_2d_fields%U(:,1))
c     V rhs and solution
      rhs(1)= Uo(1,2) - kpp_const_fields%dto*
     +     ( ftemp*.5*(Uo(1,1)+kpp_2d_fields%U(1,1)) +
     +     kpp_2d_fields%wU(0,2)/kpp_2d_fields%hm(1))
c      rhs(1)=Uo(1,2)-kpp_const_fields%dto*
c     +     ( ftemp*Uo(1,1)+
c     +     kpp_2d_fields%wU(0,2)/kpp_2d_fields%hm(1))
      do i=2,NZ-1
         rhs(i)= Uo(i,2) - kpp_const_fields%dto*ftemp*.5*(Uo(i,1)+
     +        kpp_2d_fields%U(i,1))
c         rhs(i)=Uo(i,2)-kpp_const_fields%dto*ftemp*Uo(i,1)
      enddo
      i=NZ
      rhs(i)= Uo(i,2) - kpp_const_fields%dto*ftemp*.5*(Uo(i,1)+
     +     kpp_2d_fields%U(i,1)) + kpp_2d_fields%tri(i,1,intri)*
     +     kpp_2d_fields%difm(i)*Uo(i+1,2)
c      rhs(i)=Uo(i,2)-kpp_const_fields%dto*ftemp*Uo(i,1)+
c     +     kpp_2d_fields%tri(i,1,intri)*
c     +     kpp_2d_fields%difm(i)*Uo(i+1,2)

      npd = 1
      call tridmat(cu,cc,cl,rhs,Uo(:,2),NZ,kpp_2d_fields%U(:,2))

      IF (kpp_const_fields%L_RELAX_CURR .and.
     +    kpp_2d_fields%relax_curr .gt. 0.0) THEN
          kpp_2d_fields%U(:,1)=kpp_2d_fields%U(:,1)+
     +      (kpp_2d_fields%u_clim(:)-kpp_2d_fields%U(:,1))*
     +      kpp_2d_fields%relax_curr*kpp_const_fields%dto
          kpp_2d_fields%U(:,2)=kpp_2d_fields%U(:,2)+
     +      (kpp_2d_fields%v_clim(:)-kpp_2d_fields%U(:,2))*
     +      kpp_2d_fields%relax_curr*kpp_const_fields%dto
      ENDIF
c      f(ipt)= ftemp
c      WRITE(6,*) 'After computation of U, sst = ',kpp_2d_fields%X(1,1)

c *******************************************************************
c Scalar solutions of tridiagonal matrix
c     Temperature (different from other scalars because of ghat-term
c                  and double diffusion)
c     ghatflux = wX(0,1) - (1-SWDK(-hmixe,real(time)))
c    $                     * sflux(3,5,0) / rho(ipt,0) / CP(ipt,0)
c     ghatflux = wX(0,1) - (1-SWDK(-d(1) ,real(time)))
c    $                     * sflux(3,5,0) / rho(ipt,0) / CP(ipt,0)
      ghatflux = kpp_2d_fields%wX(0,1)
      sturflux = kpp_2d_fields%wX(0,1)
      diff(0)=kpp_2d_fields%dift(0)
      ntflx(0,1)=kpp_2d_fields%wXNT(0,1)
      DO k=1,NZtmax
         diff(k)=kpp_2d_fields%dift(k)
         gcap(k)=kpp_2d_fields%ghat(k)
         ntflx(k,1)=kpp_2d_fields%wXNT(k,1)
      ENDDO
      call tridcof(diff,NZ,intri,cu,cc,cl,kpp_2d_fields)
c      WRITE(6,*) 'Before tridrhs, ntflx = ',ntflx(:,1),'diff=',diff,
c     + 'gcap=',gcap,'sturflx=',sturflux,'intri=',intri
      call tridrhs(npd,kpp_2d_fields%hm,Xo(:,1),ntflx(:,1),diff,gcap,
     +     sturflux,ghatflux,kpp_const_fields%dto,NZ,intri,rhs,
     +     kpp_2d_fields)
c      WRITE(6,*) 'After RHS on temperature, sst = ',kpp_2d_fields%X(1,1)
c     +     ,'rhs=',rhs
c     modify rhs for advection
c      do imode=1,kpp_2d_fields%nmodeadv(1)
c         adv_mode=kpp_2d_fields%modeadv(imode,1)
c         adv_mag=kpp_2d_fields%advection(imode,1)
c         call rhsmod(1,adv_mode,adv_mag,
c     +        kpp_const_fields%dto,kmixe,
c     +        kpp_2d_fields%dm(kmixe),NZ,
c     +        rhs,kpp_2d_fields,kpp_const_fields)
c      enddo

c     Surface relaxation is incompatible with
c     flux corrections at depth (NPK 12/02/08).
      IF (kpp_const_fields%L_RELAX_SST .AND. .NOT.
     +     kpp_const_fields%L_FCORR_WITHZ .AND. .NOT.
     +     kpp_const_fields%L_FCORR) THEN
c
c     Relax the Mixed layer temperature back to SST0
c     By using a flux correction at the surface
c     Added by SJW (06/04)
c
         IF (kpp_2d_fields%relax_sst .GT. 1.e-10) THEN
            IF (.NOT. kpp_const_fields%L_RELAX_CALCONLY) THEN
!               WRITE(6,*) kpp_2d_fields%relax_sst,kpp_const_fields%dto,
!     +              kpp_2d_fields%SST0,kpp_2d_fields%dm(kmixe)
               rhs(1)=rhs(1)+
     +              kpp_const_fields%dto*kpp_2d_fields%relax_sst*
     +              (kpp_2d_fields%SST0-Xo(1,1))*
     +              kpp_2d_fields%dm(kmixe)/kpp_2d_fields%hm(1)
            ENDIF
            kpp_2d_fields%fcorr=kpp_2d_fields%relax_sst*
     +           (kpp_2d_fields%SST0-Xo(1,1))*
     +           kpp_2d_fields%dm(kmixe)*kpp_2d_fields%rho(1)*
     +           kpp_2d_fields%cp(1)
         ELSE
            kpp_2d_fields%fcorr=0.0
         ENDIF
      ENDIF

c     Relax the mixed-layer temperature by using a USER-SPECIFIED
c     flux correction at the surface!  Requires setting L_FCORR
c     in the namelist and related options for input file and update
c     frequency. Values are stored in fcorr_twod (two-dimensional flux correction)
c
c     Added by NPK (29/6/08)

      IF (kpp_const_fields%L_FCORR .AND. .NOT.
     +     kpp_const_fields%L_RELAX_SST .AND.
     &     .NOT. kpp_const_fields%L_FCORR_WITHZ) THEN
         rhs(1)=rhs(1)+
     +        kpp_const_fields%dto*kpp_2d_fields%fcorr_twod/
     +        (kpp_2d_fields%rho(1)*kpp_2d_fields%cp(1)*
     +        kpp_2d_fields%hm(1))
      ENDIF
c
c     Correct the temperature at each layer in the model by using a flux
c     correction applied directly to each layer.  Requires a
c     three-dimensional (x,y,z) input file of flux corrections.
c     Input flux corrections must be in units of W/m-3.
c
c     Added by NPK (12/2/08)
c
c     Surface relaxation is incompatible with
c     flux corrections at depth (NPK 12/02/08).
      kpp_2d_fields%tinc_fcorr(:)=0.
      IF (kpp_const_fields%L_FCORR_WITHZ .AND. .NOT.
     +     kpp_const_fields%L_FCORR) THEN
         DO k=1,NZP1
            kpp_2d_fields%tinc_fcorr(k) = kpp_const_fields%dto*
     +           kpp_2d_fields%fcorr_withz(k)/
     +           (kpp_2d_fields%rho(k)*kpp_2d_fields%cp(k))
         ENDDO
         kpp_2d_fields%X(NZP1,1) = Xo(NZP1,1)+
     +		 kpp_2d_fields%tinc_fcorr(NZP1)
      ENDIF

c     Relax the temperature at each layer in the model by computing
c     a flux correction at each layer.  Requires a three-dimensional
c     (x,y,z) input file of ocean temperatures via subroutine
c     read_ocean_temperatures.
      IF (kpp_const_fields%L_RELAX_OCNT) THEN
         DO k=1,NZP1
c     Store the relaxation term as tinc_fcorr so that, on output,
c     that field contains the actual correction applied in K/timestep.
            kpp_2d_fields%tinc_fcorr(k)=kpp_2d_fields%tinc_fcorr(k)+
     +           kpp_const_fields%dto*kpp_2d_fields%relax_ocnT*
     +           (kpp_2d_fields%ocnT_clim(k)-Xo(k,1))
         ENDDO
      ENDIF
      DO k=1,NZP1
         rhs(k) = rhs(k) + kpp_2d_fields%tinc_fcorr(k)
c     Modify the correction field so that, when output, it is in
c     the correct units to be input as a flux correction via
c     L_FCORR_WITHZ (see above).
         kpp_2d_fields%ocnTcorr(k)=kpp_2d_fields%tinc_fcorr(k)*
     +        kpp_2d_fields%rho(k)*kpp_2d_fields%cp(k)/
     +        kpp_const_fields%dto
      ENDDO

c     ---
c     Ekman pumping
c     Vertical advection of temperature and salinity following Lu et al.
c     (2017, Ocean Dynamics, doi:10.1007/s10236-016-1029-9).
c
c     1. Compute Ekman depth from surface wind stress
c     2. If depth is shallower than a user-prescribed maximum depth
c        (ekmax; max_ekman_depth in NAME_FORCING), then advect.
c     3. Advection by Ekman pumping is assumed to decay sinusoidally from the
c        Ekman depth to the surface and to a user-prescribed maximum depth,
c        (ekadv_max; max_ekadv_depth in NAME_FORCING).  Advection uses
c        a centred difference in the vertical.  Follows reference above.
c
c     Curl of the wind stress must be provided in surface forcing data.
c     Positive curl = upwelling.
c     NPK 24/7/17.
c     ---
      kpp_2d_fields%ekvel(:)=0.0
      kpp_2d_fields%ekadv(:,:)=0.0
      kpp_2d_fields%hekman = 0.4*SQRT((
     +     SQRT(kpp_2d_fields%sflux(1,5,0)**2 +
     +     kpp_2d_fields%sflux(2,5,0)**2))
     +     / (kpp_2d_fields%rho(1)*kpp_2d_fields%f**2)) + 1e-16
      IF (kpp_const_fields%L_EKMAN_PUMP .and.
     +     kpp_2d_fields%hekman .lt. kpp_const_fields%ekmax) THEN
         ekvel_max = 1/(kpp_2d_fields%rho(1)*kpp_2d_fields%f)*
     +     kpp_2d_fields%sflux(7,5,0)
         DO k=2,NZ
            IF (ABS(kpp_2d_fields%zm(k)) .lt.
     +           kpp_const_fields%ekadv_max) THEN
               dist_hek = kpp_2d_fields%hekman-
     +              ABS(kpp_2d_fields%zm(k))
               ! If below Ekman depth, compute distance to bottom of Ekman layer
               IF (ABS(kpp_2d_fields%zm(k))
     +              .ge.kpp_2d_fields%hekman) THEN
                  weight_hek = SIN((kpp_const_fields%ekadv_max-
     +                 kpp_2d_fields%hekman-
     +                 ABS(dist_hek))*3.14159/2.0/
     +                 (200-kpp_2d_fields%hekman))
               ELSE
               ! If above Ekman depth, compute distance to surface
                  weight_hek = SIN((kpp_2d_fields%hekman-dist_hek)/
     +                 kpp_2d_fields%hekman*3.14159/2.0)
               ENDIF
               kpp_2d_fields%ekvel(k)=weight_hek*ekvel_max
               kpp_2d_fields%ekadv(k,1) = kpp_2d_fields%ekvel(k) *
     +              kpp_const_fields%dtsec / kpp_2d_fields%hm(k) *
     +              (kpp_2d_fields%X(k+1,1)-kpp_2d_fields%X(k-1,1))/2.0
               rhs(k) = rhs(k) + kpp_2d_fields%ekadv(k,1)
            ENDIF
         ENDDO
      ENDIF

      call tridmat(cu,cc,cl,rhs,Xo(:,1),NZ,kpp_2d_fields%X(:,1))

c     Salinity and other scalars
      DO k=0,NZtmax
         diff(k)=kpp_2d_fields%difs(k)
      ENDDO
      call tridcof(diff,NZ,intri,cu,cc,cl,kpp_2d_fields)
      do 200 n=2,NSCLR
         DO k=0,NZtmax
            ntflx(k,n)=kpp_2d_fields%wXNT(k,n)
         ENDDO
         ghatflux = kpp_2d_fields%wX(0,n)
         sturflux = kpp_2d_fields%wX(0,n)
         call tridrhs(npd,kpp_2d_fields%hm,Xo(:,n),ntflx(:,n),
     >        diff,gcap,sturflux,ghatflux,kpp_const_fields%dto,NZ,intri,
     +        rhs,kpp_2d_fields)

c     modify rhs for advections
         do imode=1,kpp_2d_fields%nmodeadv(2)
            adv_mode=kpp_2d_fields%modeadv(imode,2)
            adv_mag=kpp_2d_fields%advection(imode,2)
            call rhsmod(2,adv_mode,adv_mag,
     +           kpp_const_fields%dto,kmixe,
     +           kpp_2d_fields%dm(kmixe),
     +           NZ,rhs,kpp_2d_fields,kpp_const_fields)
         enddo

         IF (kpp_const_fields%L_EKMAN_PUMP) THEN
            DO k=2,NZ
               kpp_2d_fields%ekadv(k,n) = kpp_2d_fields%ekvel(k) *
     +              kpp_const_fields%dtsec / kpp_2d_fields%hm(k) *
     +              (kpp_2d_fields%X(k+1,n)-kpp_2d_fields%X(k-1,n))/2.0
               rhs(k) = rhs(k) + kpp_2d_fields%ekadv(k,n)
            ENDDO
         ENDIF

c     -----Added by LH (28/05/2013) modified NPK (4/7/13)

         if (n .eq. 2) then
            kpp_2d_fields%sinc_fcorr(:)=0.
c     Surface salinity  relaxation is incompatible with
c     flux corrections at depth.
            IF (kpp_const_fields%L_SFCORR_WITHZ .AND. .NOT.
     +           kpp_const_fields%L_SFCORR) THEN
               DO k=1,NZP1
                  kpp_2d_fields%sinc_fcorr(k) = kpp_const_fields%dto*
     +                 kpp_2d_fields%sfcorr_withz(k)
               ENDDO
            ENDIF
c     Relax the salinity at each layer in the model by computing
c     a flux correction at each layer.  Requires a three-dimensional
c     (x,y,z) input file of salinity via subroutine
c     read_salinity.
            IF (kpp_const_fields%L_RELAX_SAL) THEN
               DO k=1,NZP1
c     Store the relaxation term as sinc_fcorr so that, on output,
c     that field contains the actual correction applied in psu/timestep.
                  kpp_2d_fields%sinc_fcorr(k)=
     +                 kpp_2d_fields%sinc_fcorr(k)+
     +                 kpp_const_fields%dto*kpp_2d_fields%relax_sal*
     +                 (kpp_2d_fields%sal_clim(k)-Xo(k,n))
               ENDDO
            ENDIF
            DO k=1,NZP1
               rhs(k)=rhs(k)+kpp_2d_fields%sinc_fcorr(k)
c     Modify the correction field so that, when output, it is in
c     the correct units to be input as a flux correction via
c     L_SFCORR_WITHZ (see above).
               kpp_2d_fields%scorr(k)=kpp_2d_fields%sinc_fcorr(k)/
     +              kpp_const_fields%dto
            ENDDO
         ENDIF

         call tridmat(cu,cc,cl,rhs,Xo(:,n),NZ,kpp_2d_fields%X(:,n))
 200  continue
!      IF (ABS(kpp_2d_fields%hmix) .gt. 1000)
!     +     WRITE(6,*) 'hmix = ',kpp_2d_fields%hmix
!      return
      end

************************************************************************

      SUBROUTINE tridcof(diff,nzi,ind,cu,cc,cl,kpp_2d_fields)

c     Compute coefficients for tridiagonal matrix (dimension=nzi).
c     Note: cu(1) = 0. and cl(nzi) = 0. are necessary conditions.
c-----
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
c Input
      TYPE(kpp_2d_type) :: kpp_2d_fields
      integer nzi,              ! dimension of field
     +     ind                  ! index for tri-coefficients: = kmixo for t-grid,
c                                           =     1 for p-grid.
      real diff(0:nzi) ! diffusivity profile on interfaces
c Output
      real cu(nzi),    ! upper coeff. for (k-1) on k line of trid.matrix
     +     cc(nzi),    ! central ...      (k  ) ..
     +     cl(nzi)     ! lower .....      (k-1) ..
c Common tridiagonal factors (set in "init_ocn<")
c      real tri(0:NZtmax,0:1,NGRID)    ! dt/dz/dz factors in trid. matrix
c      common/ trifac / tri
c more local variables (to make implicit none)
      integer i

c
c In the surface layer
      cu(1) = 0.
      cc(1) = 1. + kpp_2d_fields%tri(1,1,ind)*diff(1)   ! 1.+ dto/h(1)/dzb(1)*diff(1)
      cl(1) =    - kpp_2d_fields%tri(1,1,ind)*diff(1)   !   - dto/h(1)/dzb(1)*diff(1)
c Inside the domain
      do 10 i=2,nzi
      cu(i) =    - kpp_2d_fields%tri(i,0,ind)*diff(i-1)
      cc(i) = 1. + kpp_2d_fields%tri(i,1,ind)*diff(i)   +
     +     kpp_2d_fields%tri(i,0,ind)*diff(i-1)
      cl(i) =    - kpp_2d_fields%tri(i,1,ind)*diff(i)
 10   continue
c In the bottom layer
      cl(nzi)= 0.
      return
      end

***********************************************************************

      SUBROUTINE tridrhs(npd,h,yo,ntflux,diff,ghat,sturflux,ghatflux,
     +                   dto,nzi,ind,rhs,kpp_2d_fields)

c     Compute right hand side of tridiagonal matrix for scalar fields:
c     =  yo (old field)
c      + flux-divergence of ghat
c      + flux-divergence of non-turbulant fluxes
c     Note: surface layer needs +dto/h(1) * surfaceflux
c           bottom  ..... ..... +dto/h(nzi)*diff(nzi)/dzb(nzi)*yo(nzi+1)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include "kpp_3d_type.com"

c Input
      TYPE(kpp_2d_type) :: kpp_2d_fields
      real dto           ! timestep interval (seconds)
      integer nzi,       ! dimension of field
     +        ind        ! index for tri-coefficients:=kmixo for t-grid,
c                              =    1 for p-grid
      real h(nzi+1),     ! layer thickness
     +     yo(nzi+1),    ! old profile
     +     ntflux(0:nzi),! non-turbulent flux = wXNT(0:nzi,1:2)
     +     diff(0:nzi),  ! diffusivity profile on interfaces
     +     ghat(nzi),    ! ghat turbulent flux
     +     sturflux,     ! surface turbulent (kinematic) flux = wX(0,n)
     +     ghatflux      ! surface flux for ghat: includes solar flux
      integer npd        ! included in list by sjw for implicit none
c Output
      real rhs(nzi)      ! right hand side
c Common tridiagonal factors (set in "init_ocn")
c      real tri(0:NZtmax,0:1,NGRID)    ! dt/dz/dz factors in trid. matrix
c      common/ trifac / tri

c more local variables to make implicit none
      integer i
      real divflx
       divflx =  1.0 / float(npd)

c In the surface layer (dto/h(1)=tri(0,1,ind)
      rhs(1)= yo(1) + dto/h(1) *
     +        (   ghatflux*diff(1)*ghat(1)
     +       - sturflux*divflx
     +       + ntflux(1) - ntflux( 0 ) )

c  Inside the domain to npd
      if(npd.ge.2) then
         do    i=2,npd
            rhs(i)= yo(i) + dto/h(i) *
c     + (  float(i)  * divflx *  ghatflux *  diff(i)  *ghat(i)
c     +  - float(i-1)* divflx *  ghatflux *  diff(i-1)*ghat(i-1)
     +           (                        ghatflux *  diff(i)  *ghat(i)
     +           -                       ghatflux *  diff(i-1)*ghat(i-1)
     +           - sturflux * divflx
     +           + ntflux(i) - ntflux(i-1) )
         end do
      endif

c Inside the rest of the domain
      do 10 i=npd+1,nzi-1
         rhs(i)= yo(i) + dto/h(i) *
     +        ( ghatflux*(diff(i)*ghat(i) - diff(i-1)*ghat(i-1))
     +        +ntflux(i) - ntflux(i-1) )
 10   continue

c     In the bottom layer
      if(nzi.gt.1) then   ! not for slab ocean
      i=nzi
      rhs(i)= yo(i) + dto/h(i) *
     +              ( ghatflux*(diff(i)*ghat(i) - diff(i-1)*ghat(i-1))
     +               +ntflux(i) - ntflux(i-1) )
     +      + yo(i+1)*kpp_2d_fields%tri(i,1,ind)*diff(i)
      endif
      return
      end

c**********************************************************************
      subroutine rhsmod(jsclr,mode,A,
     +     dto,km,dm,nzi,rhs,kpp_2d_fields,
     +     kpp_const_fields)

c     Modify rhs to correct scalar, jsclr,
c     for advection according to mode
c mode = 1 : Steady upper layer horizontal advection
c        2 : Steady mixed layer horizontal advection to km-1
c        3 : Steady horizontal advection throughout the entire column
c        4 : Steady vertical advection (= deep horizontal) below 100m
c            to bottom
c            (Change: start below 100m, instead of at layer 16, and
c            do not advect into bottom layer, 7-1-93)
c        5 : Steady bottom diffusion
c        6 : Seasonal mixed layer horizontal advection to dm
c        7 : Seasonal thermocline horizontal advection to 1.5 dm
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
c#include "ocn_paras.com"
c#include "local_pt.com"
c Input
      integer nzi,              ! vertical dimension of field
     +     km,                  ! index of gridpoint just below h
     +     mode,                ! type of advection
     +     jsclr                ! scalar
c     real h(nzi+1),     ! layer thickness
c     +     z(nzi+1),     ! z grid levels (added as input on 7-1-93)
      real rhs(nzi)             ! right hand side from tridrhs
c      DOUBLE PRECISION time     ! time in days from jan 1 of any year.
      real dto,                 ! ocean time step
c     +     dpy,                 ! days per year (added as input on 7-1-93)
     +     dm,                  ! depth d(km+.5)
     +     A                    ! advection of heat(W/m2) or Salt(PSU m/s)

      TYPE(kpp_2d_type) :: kpp_2d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

c Output
c     real rhs(nzi)      ! modified right hand side

      real am,fact,delta,depth,dmax
      integer n,n1,n2,nzend

c Internal
      real f(12)     ! monthly partion of annual advection
      real xsA(21)   ! yearly excess of heat

c     data f/.1,.1,6*0.0,.1,.3,.4,.2/
      data f/.05,.05,5*0.0,.05,.15,.20,.30,.20/
      data xsA/48.26,21.73,29.02,56.59,19.94,15.96,18.28,
     &         40.52,37.06,29.83,29.47,15.77, 1.47,14.55,
     &          4.22,28.19,39.54,19.58,20.27,11.19,21.72/

      if(mode.le.0) return
c               find ocean year
c        iyr = 1 + idint((time-75.0)/dpy)
c        day = time - dpy * (idint(time/dpy))  ! 365.25))
c        month = 1 + int(12. * day / dpy)    ! 367.)
c        if(month.gt.12) then
c           write(nuerr,*) 'STOP rhsmod (ocn.f):'
c           write(nuerr,*) '     rounding error, month gt 12 =',month
c           stop 97
c        endif

c       Am = -12. * f(month) * (xsA(iyr) - 0.0 )       ! Annual
c        Am =  12. * f(month) * A                       ! Seasonal
       Am = A                                          ! Steady

      if(mode.eq.1) then
c                          correct upper layer advection
         if(jsclr.eq.1) fact = dto * Am /
     +        (kpp_2d_fields%rho(1)*kpp_2d_fields%cp(1))
         if(jsclr.eq.2) fact = dto * Am * 0.033
         rhs(1) = rhs(1)
     &        + fact / kpp_2d_fields%hm(1)

      else if(mode.eq.2) then
c     correct mixed layer advection
         delta = 0.0
         do n=1,km-1
            delta = delta + kpp_2d_fields%hm(n)
         enddo
         do 215 n=1,km-1
            if(jsclr.eq.1) fact = dto * Am /
     +           (kpp_2d_fields%rho(n)*kpp_2d_fields%cp(n))
            if(jsclr.eq.2) fact = dto * Am * 0.033
            rhs(n) = rhs(n)
     &           + fact  / delta
 215     continue

      else if (mode.eq.3) then
c     throughout whole water column
         delta = 0.0
         do n=1,nzi
            delta = delta + kpp_2d_fields%hm(n)
         enddo
         do 315 n=1,nzi
            if(jsclr.eq.1) fact = dto * Am /
     +           (kpp_2d_fields%rho(n)*kpp_2d_fields%cp(n))
            if(jsclr.eq.2) fact = dto * Am * 0.033
            rhs(n) = rhs(n)
     &           + fact / delta
 315     continue


      else if (mode.eq.4) then
c     vertical advection = deep horizontal
         nzend=nzi-1            ! nzend=nzi (change:7-1-93)
         n1=0                   ! n1=16     (change:7-1-93)
 401     n1=n1+1
         if(kpp_2d_fields%zm(n1).ge.-100.) goto 401
         delta = 0.0
         do n=n1,nzend
            delta = delta + kpp_2d_fields%hm(n)
         enddo
         do 415 n=n1,nzend
            if(jsclr.eq.1) fact = dto * Am /
     +           (kpp_2d_fields%rho(n)*kpp_2d_fields%cp(n))
            if(jsclr.eq.2) fact = dto * Am * 0.033
            rhs(n) = rhs(n)
     &           + fact / delta
 415     continue

      else if(mode.eq.5) then
c     correct bottom layer diffusion
         if(jsclr.eq.1) fact = dto * Am / (kpp_2d_fields%rho(nzi)*
     +        kpp_2d_fields%cp(nzi))
         if(jsclr.eq.2) fact = dto * Am * 0.033
         rhs(nzi) = rhs(nzi)
     &        + fact / kpp_2d_fields%hm(nzi)

      else

c     seasonal mixed layer or thermocline advection
c     find ocean year
c     iyr = 1 + idint((time-75.0)/dpy)
c     day = time - dpy * (idint(time/dpy))  ! 365.25))
c     month = 1 + int(12. * day / dpy)    ! 367.)
c     diag
c     if(month.gt.12) then
c     write(nuerr,*) 'STOP rhsmod (ocn.f):'
c     write(nuerr,*) '     rounding error, month gt 12 =',month
c     stop 97
c     endif
c     diag
c     Am = -12. * f(month) * (xsA(iyr) - 0.0 )       ! Annual
c     Am =  12. * f(month) * A                       ! Seasonal
c     Am = A                                          ! Steady

         if(mode.eq.6) then
c     mixed layer to dm
            n1 = 1
            depth = kpp_2d_fields%hm(1)
            dmax  = dm -  0.5 * (kpp_2d_fields%hm(km) +
     +           kpp_2d_fields%hm(km-1))
            delta = 0.0
            do 605 n =n1,nzi
               n2    = n
               delta = delta + kpp_2d_fields%hm(n)
               depth = depth + kpp_2d_fields%hm(n+1)
               if(depth.ge.dmax) go to 606
 605        continue
 606        continue

         else if (mode.eq.7) then
c     thermocline to 100m
            n1 = km - 1
            depth = dm - 0.5 * kpp_2d_fields%hm(km)
            dmax = 100.
            delta = 0.0
            do 705 n=n1,nzi
               n2 = n
               delta = delta + kpp_2d_fields%hm(n)
               depth = depth + kpp_2d_fields%hm(n+1)
               if(depth.ge.dmax) go to 706
 705        continue
 706        continue

         else
            write(nuerr,*) 'STOP in rhsmod (ocn.f):'
            write(nuerr,*) '      mode out of range, mode=',mode
            CALL MIXED_ABORT
         endif

c     Finish both 6 and 7 here
         do 615 n=n1,n2
            if(jsclr.eq.1) fact = dto * Am /
     +           (kpp_2d_fields%rho(n)*kpp_2d_fields%cp(n))
            if(jsclr.eq.2) fact = dto * Am * 0.033
            rhs(n) = rhs(n) + fact  / delta
 615     continue

      endif

      return
      end

***********************************************************************

      SUBROUTINE tridmat(cu,cc,cl,rhs,yo,nzi,yn)
c
c     Solve tridiagonal matrix for new vector yn, given right hand side
c     vector rhs. Note: yn(nzi+1) = yo(nzi+1).
c-----
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include "parameter.inc"
c Input
      integer nzi               ! dimension of matrix
      real cu (nzi),            ! upper coeff. for (k-1) on k line of tridmatrix
     +     cc (nzi),            ! central ...      (k  ) ..
     +     cl (nzi),            ! lower .....      (k-1) ..
     +     rhs(nzi),            ! right hand side
     +     yo(nzi+1),yni            ! old field
c     +     diff(0:nzi)
c Output
      real yn(nzi+1)    ! new field
c Local
      real gam(NZtmax), ! temporary array for tridiagonal solver
     +     bet          ! ...
c more local for implicit none
      integer i
c Solve tridiagonal matrix.
      bet   = cc(1)
c     yn(1) = (rhs(1) + tri(0,1,ind)*surflux) / bet    ! surface
      yn(1) =  rhs(1) / bet    ! surface
      do 21 i=2,nzi
         gam(i)= cl(i-1)/bet
         bet   = cc(i) - cu(i)*gam(i)
         if(bet.eq.0.) then
            write(nuerr,*)'* algorithm for solving tridiag matrix fails'
            write(nuerr,*)'* bet=',bet
            write(nuerr,*)'*i-1=',i-1,' cc=',cc(i-1),'cl=',cl(i-1)
            write(nuerr,*)'*i=',i,' cc=',cc(i),' cu=',cu(i),
     +           ' gam=',gam(i)
            CALL MIXED_ABORT
            bet=1.E-12
c     Pause 3
         endif
c     to avoid "Underflow" at single precision on the sun
         yn(i) =      (rhs(i)  - cu(i)  *yn(i-1)  )/bet
c     yni   =      (rhs(i)  - cu(i)  *yn(i-1)  )/bet
c     yn(i) = max( (rhs(i)  - cu(i)  *yn(i-1)  )/bet , 1.E-12 )
c     if(yni.lt.0.)
c     +        yn(i) = min( (rhs(i)  - cu(i)  *yn(i-1)  )/bet ,-1.E-12 )
 21   continue

c     yn(nzi)  = (rhs(nzi)- cu(nzi)*yn(nzi-1)
c    +                    + tri(nzi,1,ind)*diff(nzi)*yo(nzi+1) )/bet
c                                                      ! bottom
      do 22 i=nzi-1,1,-1
         yn(i)  = yn(i) - gam(i+1)*yn(i+1)
 22   continue
      yn(nzi+1) = yn(nzi+1)
c      yn(nzi+1) = yo(nzi+1)

      return
      end

*********************************************************************

      SUBROUTINE init_ocn(kpp_3d_fields,kpp_const_fields)
c
c     Initialize ocean model:
c     Set coefficients for tridiagonal matrix solver.
c     Compute hmix and diffusivity profiles for initial profile.
c     Prepare for first time step.
c
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include "kpp_3d_type.com"

c     Input
c      real  U(NPTS,NZP1,NVEL),X(NPTS,NZP1,NSCLR)
      TYPE(kpp_const_type) :: kpp_const_fields

c     Output
c      real dbloc(NPTS,nz),shsq(NPTS,nzp1),Rig(NPTS,nzp1)
      TYPE(kpp_3d_type) :: kpp_3d_fields

c     Local
      TYPE(kpp_2d_type) :: kpp_2d_fields
c      real U0(NZP1,NVEL),X0(NZP1,NSCLR)

c     Common Blocks
c      real tri(0:NZtmax,0:1,NGRID)! dt/dz/dz factors in trid. matrix
c      common/ trifac / tri
c      logical L_initflag
c      common/initflag/L_initflag

c      real hmixd(NPTS,0:1),     ! storage arrays for extrapolations
c     +     Us(NPTS,NZP1,NVEL ,0:1),    ! ..      ..     ..  ..
c     +     Xs(NPTS,NZP1,NSCLR,0:1)     ! ..      ..     ..  ..
c     integer old,new           ! extrapolation index for Us,Xs,hmixd
c     common/ saveUXh /
c     +     old,new,Us,Xs,hmixd

c Local variables
      real dzb(NZ)              ! diff. between grid-levels below z(j)
c     more local for implicit none
      integer k,kmix0,n,l,ipt
      real hmix0,deltaz
c
c     Compute factors for coefficients of tridiagonal matrix elements.
c     tri(0     ,1,.........) : dt/h(1) factor for rhs flux
c     tri(k=1:NZ,0,.........) : dt/h(k)/ {dzb(k-1)=z(k-1)-z(k)=dzabove}
c     tri(k=1:NZ,1,.........) : dt/h(k)/ {dzb(k  )=z(k)-z(k+1)=dzbelow}
c

      kpp_2d_fields%dm = kpp_const_fields%dm
      kpp_2d_fields%hm = kpp_const_fields%hm
      kpp_2d_fields%zm = kpp_const_fields%zm
      IF (kpp_const_fields%L_SLAB) THEN
         kpp_2d_fields%dm(1)=kpp_const_fields%slab_depth
         kpp_2d_fields%zm(1)=kpp_const_fields%slab_depth*(-0.5)
         kpp_2d_fields%hm(1)=kpp_const_fields%slab_depth
      ENDIF
      do 10 k=1,NZ
         dzb(k)     = kpp_2d_fields%zm(k) - kpp_2d_fields%zm(k+1)
 10   continue
      kpp_const_fields%tri(0,1,1) = kpp_const_fields%dto/
     +     kpp_2d_fields%hm(1)
      kpp_const_fields%tri(1,1,1) = kpp_const_fields%dto/
     +     kpp_2d_fields%hm(1)/dzb(1)
      do 20 k=2,NZ
         kpp_const_fields%tri(k,1,1) = kpp_const_fields%dto/
     +        kpp_2d_fields%hm(k)/dzb(k)
         kpp_const_fields%tri(k,0,1) = kpp_const_fields%dto/
     +        kpp_2d_fields%hm(k)/dzb(k-1)
 20   continue

      IF ( .NOT. kpp_const_fields%L_RESTART) THEN
c
c     Determine hmix for initial profile:
c#ifdef OPENMP
c!$OMP PARALLEL DEFAULT(PRIVATE) SHARED(kpp_3d_fields,kpp_const_fields)
c!$OMP DO SCHEDULE(dynamic)
c#endif
         DO ipt=1,npts
            CALL kpp_fields_3dto2d(kpp_3d_fields,ipt,kpp_2d_fields)
            IF (kpp_2d_fields%L_OCEAN) THEN
c               DO k=1,NZP1
c                  DO i=1,nvel
c                     U0(k,i)=kpp_3d_fields%U(ipt,k,i)
c                  ENDDO
c                  DO i=1,nsclr
c                     X0(k,i)=kpp_3d_fields%X(ipt,k,i)
c                  ENDDO
c               ENDDO
               kpp_2d_fields%L_INITFLAG=.TRUE.
               CALL vmix(kpp_2d_fields,kpp_const_fields,hmix0,kmix0)
               kpp_2d_fields%L_INITFLAG=.FALSE.
               IF (.NOT. kpp_const_fields%L_SLAB) THEN
                  kpp_2d_fields%hmix = hmix0
                  kpp_2d_fields%kmix = kmix0
               ELSE
                  kpp_2d_fields%hmix = kpp_const_fields%slab_depth
                  kpp_2d_fields%kmix = 1
               ENDIF
               kpp_2d_fields%Tref = kpp_2d_fields%X(1,1)
c Evaluate initial fluxes (to write to output data file)
               do k=1,NZ
                  deltaz = 0.5*(kpp_2d_fields%hm(k)+
     +                 kpp_2d_fields%hm(k+1))
                  do n=1,NSCLR
                     kpp_2d_fields%wX(k,n)=-kpp_2d_fields%difs(k)*
     +                    ((kpp_2d_fields%X(k,n)-
     +                    kpp_2d_fields%X(k+1,n))/deltaz-
     +                    kpp_2d_fields%ghat(k)*kpp_2d_fields%wX(0,n))
                  enddo
                  if(kpp_const_fields%LDD)
     +                 kpp_2d_fields%wX(k,1)=-kpp_2d_fields%dift(k)*
     +                 ((kpp_2d_fields%X(k,1)-
     +                 kpp_2d_fields%X(k+1,1))/deltaz-
     +                 kpp_2d_fields%ghat(k)*kpp_2d_fields%wX(0,1))
                  kpp_2d_fields%wX(k,nsp1)= kpp_const_fields%grav *
     +                 (kpp_2d_fields%talpha(k)*kpp_2d_fields%wX(k,1) -
     +                 kpp_2d_fields%sbeta(k) * kpp_2d_fields%wX(k,2))
                  do  n=1,NVEL
                     kpp_2d_fields%wU(k,n)= -kpp_2d_fields%difm(k)*
     +                    (kpp_2d_fields%U(k,n)-kpp_2d_fields%U(k+1,n))/
     +                    deltaz
                  enddo
               enddo

c     Prepare for first time step

c     indices for extrapolation
               kpp_2d_fields%old = 0
               kpp_2d_fields%new = 1
c     initialize array for extrapolating hmixd,Us,Xs
               kpp_2d_fields%hmixd(0) = kpp_2d_fields%hmix
               kpp_2d_fields%hmixd(1) = kpp_2d_fields%hmix
               do k=1,NZP1
                  do l=1,NVEL
                     kpp_2d_fields%Us(k,l,0)=
     +                    kpp_2d_fields%U(k,l)
                     kpp_2d_fields%Us(k,l,1)=
     +                    kpp_2d_fields%U(k,l)
                  enddo
                  do l=1,NSCLR
                     kpp_2d_fields%Xs(k,l,0)=kpp_2d_fields%X(k,l)
                     kpp_2d_fields%Xs(k,l,1)=kpp_2d_fields%X(k,l)
                  enddo
               enddo
c     IF (ipt .eq. 3849) THEN
c     WRITE(6,*) 'lat = ',kpp_2d_fields%dlat,
c     +                 'lon = ',kpp_2d_fields%dlon
c     WRITE(6,*) 'difm = ',kpp_2d_fields%difm
c     STOP
c     ENDIF
            ENDIF
            !WRITE(6,*) 'Initial hmix = ',kpp_2d_fields%hmix
            CALL kpp_fields_2dto3d(kpp_2d_fields,ipt,kpp_3d_fields)
         ENDDO
c#ifdef OPENMP
c!$OMP END DO
c!$OMP END PARALLEL
c#endif
      ENDIF
      return
      end

*********************************************************************
      subroutine swfrac( fact, z, kpp_2d_fields, swdk )
c     compute fraction of solar short-wave flux penetrating to specified
c     depth (times fact) due to exponential decay in  Jerlov water type
c     reference : two band solar absorption model of simpson and
c     paulson (1977)

      IMPLICIT NONE
#include "kpp_3d_type.com"
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      TYPE(kpp_2d_type) :: kpp_2d_fields
      integer nwtype
      parameter(nwtype=5) ! max number of different water types
c
c  model
c#include "local_pt.com"

c      integer imt         ! number of horizontal grid points

c  input
      real fact           ! scale  factor to apply to depth array
      real z         ! vertical height ( <0.) for desired sw
c                           fraction                                 (m)
      integer jwtype ! index for jerlov water type

c  output
      real swdk      !  short wave (radiation) fractional decay

c  local
      real  rfac(nwtype),a1(nwtype),a2(nwtype)
      real rmin,r1,r2
c      integer i
c      save  rfac,a1,a2,rmin
c
c     jerlov water type :  I       IA      IB      II      III
c                jwtype    1       2       3       4       5
c
      data rfac         /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
      data a1           /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
      data a2           / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /
      data rmin         / -80. /
c
c      do 100 i = ipt,ipt

c         r1      = MAX(z*fact/a1(jwtype), rmin)
c	 r2      = MAX(z*fact/a2(jwtype), rmin)
c         swdk =      rfac(jwtype)  * exp(r1)
c     $        + (1.-rfac(jwtype)) * exp(r2)

      r1 = MAX(z*fact/kpp_2d_fields%h1,rmin)
      r2 = MAX(z*fact/kpp_2d_fields%h2,rmin)     
      swdk = kpp_2d_fields%rfac*exp(r1)+(1-kpp_2d_fields%rfac)*exp(r2)

c  100 continue

      return
      end

      subroutine budget (X1,kpp_2d_fields,kpp_const_fields)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include "kpp_3d_type.com"

      type(kpp_2d_type) :: kpp_2d_fields
      type(kpp_const_type) :: kpp_const_fields
      real X1,X2
      dimension X1(nzp1,nsclr)

      real T1,T2,S1,S2
      real fltn
      real dhdt,dfdt,qtop,ftop,qbot,fbot
      real delt,fact,rhs,diff
      integer k

      T1 = 0.0
      T2 = 0.0
      S1 = 0.0
      S2 = 0.0
      fltn = 1. / float(nz)
      write(nuout,*)
      do 5 k=1,nz
      T1  =  T1 + X1(k,1)  * kpp_2d_fields%hm(k)
      T2  =  T2 + kpp_2d_fields%X(k,1)  * kpp_2d_fields%hm(k)
      S1  =  S1 + X1(k,2)  * kpp_2d_fields%hm(k)
      S2  =  S2 + kpp_2d_fields%X(k,2)  * kpp_2d_fields%hm(k)
      write(nuout,*) k,X1(k,1),kpp_2d_fields%X(k,1),X1(k,2),
     +     kpp_2d_fields%X(k,2)
  5   continue
      dhdt = (T2 - T1) * kpp_2d_fields%rho(0) * kpp_2d_fields%cp(0) /
     +     kpp_const_fields%dto
      dfdt = (S1 - S2) * kpp_2d_fields%rho(0) / kpp_2d_fields%Ssurf /
     +     kpp_const_fields%dto
      Qtop = kpp_2d_fields%sflux(3,5,0) + kpp_2d_fields%sflux(4,5,0)
      Ftop = kpp_2d_fields%sflux(6,5,0)
      Qbot = -kpp_2d_fields%rho(0) * kpp_2d_fields%cp(0) *
     +     (kpp_2d_fields%wX(nz,1) + kpp_2d_fields%wXNT(nz,1) )
      Fbot =  kpp_2d_fields%rhoh2o / kpp_2d_fields%Ssurf  *
     +     kpp_2d_fields%wX(nz,2)

      write(nuout,*) 'heat ',Qtop,dhdt,Qbot
      write(nuout,*) 'salt ',Ftop,dfdt,Fbot
      write(nuout,*) kpp_const_fields%dto,kpp_2d_fields%rho(0),
     +     kpp_2d_fields%CP(0),kpp_2d_fields%Ssurf,
     +     kpp_2d_fields%sflux(3,5,0),kpp_2d_fields%sflux(4,5,0)

      do 15 k=1,nz
         delt = (kpp_2d_fields%X(k,1)-X1(k,1))
         fact = kpp_const_fields%dto / kpp_2d_fields%hm(k)
      rhs  = fact * (kpp_2d_fields%wX(k,1)-kpp_2d_fields%wX(k-1,1) +
     +     kpp_2d_fields%wXNT(k,1)-kpp_2d_fields%wXNT(k-1,1) )
      diff = delt - rhs
      write(nuout,*) kpp_2d_fields%wX(k-1,1),kpp_2d_fields%wXNT(k-1,1),
     +     kpp_2d_fields%wX(k,1),kpp_2d_fields%wXNT(k,1)
      write(nuout,*) fact,delt,rhs,diff
 15   continue

      return
      end


c *********************************************************


      SUBROUTINE zint(nz,zed,Gout,mz,zm,Gint)

      dimension zed(0:nz+1), Gout(0:nz+1)
      dimension  zm(mz+1), Gint(0:mz+1)

      Gint(0) = Gout(0)
      k = 1
      do 35 j=1,mz+1
  45  continue
          if(zed(k).gt.zm(j)) then
            k = k+1
            go to 45
          else
            Gint(j) = Gout(k-1) + (Gout(k)-Gout(k-1)) *
     >                          (zm(j)-zed(k-1))/(zed(k)-zed(k-1))
          endif

 35   continue

      return
      end


************************************************************************
      SUBROUTINE Identity(nzi,cu,cc,cl)
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include "parameter.inc"
      integer nzi
      real cu(nzi),cc(nzi),cl(nzi)
c      real tri(0:NZtmax,0:1,NGRID)
c      common/ trifac / tri
      integer i

      do i=1,nzi
      cu(i) = 0.
      cc(i) = 1.
      cl(i) = 0.
      enddo
      return
      end



************************************************************************
      subroutine swfrac_opt( fact, kpp_2d_fields,
     +     kpp_const_fields)
c     compute fraction of solar short-wave flux penetrating to specified
c     depth (times fact) due to exponential decay in  Jerlov water type
c     reference : two band solar absorption model of simpson and
c     paulson (1977)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      integer nwtype
      parameter(nwtype=5) ! max number of different water types
c
c  model
! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
c#include "local_pt.com"
c#include "times.com"
c#include "vert_pgrid.com"

c      integer imt         ! number of horizontal grid points

c  input
      real fact                 ! scale  factor to apply to depth array
c     integer k                 ! index of vertical grid point
c     real z                    ! vertical height ( <0.) for desired sw
c                                    fraction (m)
      type(kpp_2d_type) :: kpp_2d_fields
      type(kpp_const_type) :: kpp_const_fields

c  output
c      real swdk                 !  short wave (radiation) fractional decay

c  local
c      real swfrac_save(NPTS,NZP1)
      real  rfac(nwtype),a1(nwtype),a2(nwtype)
      real rmin,r1,r2
      integer l
c      save  rfac,a1,a2,rmin
c      common /save_swfrac/swfrac_save
c
c     jerlov water type :  I       IA      IB      II      III
c                jwtype    1       2       3       4       5
c
      data rfac         /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
      data a1           /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
      data a2           / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /
      data rmin         / -80. /
c
c      IF ((kpp_const_fields%ntime .eq. 1) .and. (k .eq. 2)) THEN
      DO l=1,NZP1
c     do 100 i = ipt,ipt
c         r1      = MAX(kpp_2d_fields%zm(l)*
c     +        fact/a1(kpp_2d_fields%jerlov), rmin)
         r1 = MAX(kpp_2d_fields%zm(l)*fact/kpp_2d_fields%h1,rmin)
         
c         r2      = MAX(kpp_2d_fields%zm(l)*
c     +        fact/a2(kpp_2d_fields%jerlov), rmin)
         r2 = MAX(kpp_2d_fields%zm(l)*fact/kpp_2d_fields%h2,rmin)
         
c         kpp_2d_fields%swfrac(l) = rfac(kpp_2d_fields%jerlov)  *
c     +        exp(r1) + (1.-rfac(kpp_2d_fields%jerlov)) * exp(r2)
         kpp_2d_fields%swfrac(l) = kpp_2d_fields%rfac*exp(r1)+
     +        (1-kpp_2d_fields%rfac)*exp(r2)
         
c     100        continue
      ENDDO

c      DO i=ipt,ipt
c     swdk=kpp_2d_fields%swfrac(l)
c      ENDDO

      return
      end
