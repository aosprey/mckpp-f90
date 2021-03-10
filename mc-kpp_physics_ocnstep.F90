SUBROUTINE mckpp_physics_ocnstep(kpp_1d_fields,kpp_const_fields)
#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_1d_type,kpp_const_type
#else
  USE mckpp_data_types
#endif
  USE mckpp_parameters
  
  !-----------------------------------------------------------------------
  ! Note in this version:
  !   -  ADVECTIVE CORRECTIONS AVAILABLE
  !   - read forcing, but compute SW from okta model   (fread in atmrad.f)
  !   - use okta model for PAPA          (use cloudpapa in solar in bio.f)
  !   - Jerlov water type II                            (SWDK in fluxes.f)
  !   - albedo for ocean is set to 0.06, which is used for QSW from fcomp
  !        or fread when rad/conv is not running:
  !        albocn=0.06                             (init cnsts in input.f)
  !   - no net fresh water flux into ocean when forcing with state
  !        variables (laflx >= 1):
  !        sflux(7,2,jptr) = - sflux(6,2,jptr)        (atmflx in fluxes.f)
  !   - use psnow flux data to input observed SST, 
  !                                                    (fread in atmrad.f)
  !  ALSO :
  !         General ability to read large scale forcing from euc and doc LES
  ! ----------------------------------------------------------------------
  ! Originally ~/KPP/LARGE/ocn.f (from Bill Large)
  !  Modified by SJW to remove those bits which a specific to Bill's 
  ! work and just leave the bit we need
  ! Started 18/03/02
  !-----------------------------------------------------------------------
  
  !     Main driver for ocean module.
  !     Integration is performed only on the permanent grid
  !     Written   3 Mar 1991 - WGL
  !     Modified  5 Jun 1992 - jan : implicit scheme
  !              16 Nov      - jan : latest version
  !              16 Nov 1994 - wgl : new KPP codes no temporary grid
  
  IMPLICIT NONE
  
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  real hmixe,hmixn,tol
  real Uo(NZP1,NVEL),Xo(NZP1,NSCLR)
  real Ux(NZP1,NVEL),XX(NZP1,NSCLR) ! Additional variables to provide smoothing in the iteration.
  real Ui                   ! Ui used in damping (LH 8/08/2013)
  real dampU(NVEL)          ! dampU used to flag which Ui chosen in damping (LH 29/08/2013)
  real lambda               ! Factor to control smoothing
  integer iter,iconv        ! number of iterations
  integer kmixe,kmixn
  
  ! More Local Variables (to make implicit none)
  real deltaz,rhonot,a,b
  integer k,l,n

  ! These aren't used elsewhere in the code but maybe they should go in kpp_fields
  REAL :: eflx, esnk, Tmke, Ptke    
  REAL, ALLOCATABLE, DIMENSION(:) :: rmke(:)  
  
  ! Number of iterations for computational instability
  integer comp_iter_max
  real rmsd(4),rmsd_threshold(4)
  data comp_iter_max /10/
  ! Critical depth-integrated RMS difference between old and new profiles
  ! for repeating integration, for (/U,V,T,S/). Typical (stable) 
  ! values are O(10^-2) for U and V, O(10^-3) for T and O(10^-4) for S.
  ! NPK 17/5/13
  data rmsd_threshold /1,1,1,1/
  
  data lambda /0.5/

  ALLOCATE( rmke(nzp1) )
  
  Uo=kpp_1d_fields%U(:,:)
  Xo=kpp_1d_fields%X(:,:)
  kpp_1d_fields%comp_flag=.TRUE.
  kpp_1d_fields%reset_flag=0
  kpp_1d_fields%dampu_flag=0
  kpp_1d_fields%dampv_flag=0
  
  DO WHILE (kpp_1d_fields%comp_flag .and. kpp_1d_fields%reset_flag .le. comp_iter_max)
     ! Estimate new profiles by  extrapolation
     DO k=1,NZP1
        DO l=1,NVEL
!	   WRITE(6,*) 'k=',k,'l=',l,'new=',kpp_1d_fields%new,'old=',&
!         kpp_1d_fields%old,'Us(new)=',kpp_1d_fields%Us(k,l,kpp_1d_fields%new),&
!         'Us(old)=',kpp_1d_fields%Us(k,l,kpp_1d_fields%old)
           IF (kpp_1d_fields%old .lt. 0 .or. kpp_1d_fields%old .gt. 1) THEN
              WRITE(6,*) 'Dodgy value of old at k=',k,'l=',l,'old=',kpp_1d_fields%old
              kpp_1d_fields%old=kpp_1d_fields%new
           ENDIF
           IF (kpp_1d_fields%new .lt. 0 .or. kpp_1d_fields%new .gt. 1) THEN
              WRITE(6,*) 'Dodgy value of new at k=',k,'l=',l,'new=',kpp_1d_fields%new
              kpp_1d_fields%new=kpp_1d_fields%old
           ENDIF        
           kpp_1d_fields%U(k,l)=2.*kpp_1d_fields%Us(k,l,kpp_1d_fields%new)- &
                kpp_1d_fields%Us(k,l,kpp_1d_fields%old)
           Ux(k,l)=kpp_1d_fields%U(k,l)
        ENDDO
        DO l=1,NSCLR
           kpp_1d_fields%X(k,l)=2.*kpp_1d_fields%Xs(k,l,kpp_1d_fields%new)-&
                kpp_1d_fields%Xs(k,l,kpp_1d_fields%old)
           Xx(k,l)=kpp_1d_fields%X(k,l)
        ENDDO
     ENDDO
     
     ! Iteration loop for semi-implicit integration
     ! Reset iteration counter
     iter=0
     iconv=0
     
     ! This loop controls the number of compulsory iterations
     ! The original value was 2, using (the upper value of the loop +2)
     ! added by SJW (17 Jan 03) to try an alleviate some non-convergences
     DO iter=0,2
        DO k=1,NZP1
           DO l=1,NVEL
              kpp_1d_fields%U(k,l)=lambda*Ux(k,l)+(1-lambda)*kpp_1d_fields%U(k,l)
              Ux(k,l)=kpp_1d_fields%U(k,l)
           ENDDO
           DO l=1,NSCLR
              kpp_1d_fields%X(k,l)=lambda*Xx(k,l)+(1-lambda)*kpp_1d_fields%X(k,l)
              Xx(k,l)=kpp_1d_fields%X(k,l)
           ENDDO
        ENDDO
        call MCKPP_PHYSICS_VERTICALMIXING(kpp_1d_fields,kpp_const_fields,hmixe,kmixe)
        call MCKPP_PHYSICS_OCNINT(kpp_1d_fields,kpp_const_fields,1,kmixe,Uo,Xo)
     ENDDO
     ! The original code can be restored by reseting iter=1 and removing the  
     ! above loop  
     ! iter=1
         
     IF (kpp_const_fields%LKPP) THEN
45      continue
        DO k=1,NZP1
           DO l=1,NVEL
              kpp_1d_fields%U(k,l)=lambda*Ux(k,l)+(1-lambda)*kpp_1d_fields%U(k,l)
              Ux(k,l)=kpp_1d_fields%U(k,l)
           ENDDO
           DO l=1,NSCLR
              kpp_1d_fields%X(k,l)=lambda*Xx(k,l)+(1-lambda)*kpp_1d_fields%X(k,l)
              Xx(k,l)=kpp_1d_fields%X(k,l)
           ENDDO
        ENDDO
        call MCKPP_PHYSICS_VERTICALMIXING(kpp_1d_fields,kpp_const_fields,hmixn,kmixn)       
        call MCKPP_PHYSICS_OCNINT(kpp_1d_fields,kpp_const_fields,1,kmixn,Uo,Xo)
        iter = iter + 1
         
        ! check iteration for convergence
        tol = hmixtolfrac*kpp_const_fields%hm(kmixn)
        if(kmixn.eq.NZP1) tol = hmixtolfrac*kpp_const_fields%hm(NZ)
        if(abs(hmixn-hmixe).gt.tol)  then

           ! Uncommeting the following the lines iconv=0 to IF (iconv ...)
           ! will make the model do two consecutive tests for convergence of the 
           ! hmix (added by SJW 17 Jan 03). This did not work well in testing for
           ! long timestep, high resolution (the model generally failed to satisfy the 
           ! convergence test on two consecutive iterations.
           iconv=0
        ELSE
           iconv=iconv+1
        ENDIF
        IF (iconv .lt. 3) THEN
           if (iter.lt.itermax) then
              hmixe = hmixn
              kmixe = kmixn
              goto 45
           else
              ! use shallower hmix
              if(hmixn.gt.hmixe) then
                 hmixe = hmixn ! comment out for hmix data
                 kmixe = kmixn ! ..      ..  ..  hmix data
                 goto 45    ! ..      ..  ..  hmix data 
              endif
           endif
        endif
        if( iter.gt.(itermax+1) ) then 
           write(nuout,1009) kpp_const_fields%ntime,& ! comment out for hmix data
                kpp_1d_fields%dlon,kpp_1d_fields%dlat,hmixe,hmixn,hmixn-hmixe,kmixn,iter
1009       format('  long iteration at',i6,' steps',/,' location=(',f7.2,',',f6.2,')',/,&
                '  hmixest=',f7.2,' hmixnew=',f7.2,' diff=',f6.1,' kmixn=',i3,' iteration=',i3)
        endif
     ENDIF

     ! Trap for profiles that are very different from original profile
     ! or clearly erroneous, to detect rare instances of instability
     ! in the semi-implicit integration.  Reset to original profile,
     ! add some noise via changing Coriolis term slightly, and try
     ! integration again.
     ! NPK 16/5/2013
     kpp_1d_fields%comp_flag=.FALSE.
     DO k=1,NZ
        IF (ABS(kpp_1d_fields%U(k,1)).ge. 10 .or. ABS(kpp_1d_fields%U(k,2)).ge.10 .or. &
             ABS(kpp_1d_fields%X(k,1)-kpp_1d_fields%X(k+1,1)) .ge. 10) THEN 
           kpp_1d_fields%comp_flag=.TRUE.
           kpp_1d_fields%f=kpp_1d_fields%f*1.01
        ENDIF
     END DO
     IF (.NOT. kpp_1d_fields%comp_flag) THEN
        rmsd(:)=0.
        DO k=1,NZP1
           rmsd(1)=rmsd(1)+(kpp_1d_fields%U(k,1)-Uo(k,1))*(kpp_1d_fields%U(k,1)-Uo(k,1))*&
                kpp_const_fields%hm(k)/kpp_const_fields%dm(NZ)
           rmsd(2)=rmsd(2)+(kpp_1d_fields%U(k,2)-Uo(k,2))*(kpp_1d_fields%U(k,2)-Uo(k,2))*&
                kpp_const_fields%hm(k)/kpp_const_fields%dm(NZ)
           rmsd(3)=rmsd(3)+(kpp_1d_fields%X(k,1)-Xo(k,1))*(kpp_1d_fields%X(k,1)-Xo(k,1))*&
                kpp_const_fields%hm(k)/kpp_const_fields%dm(NZ)
           rmsd(4)=rmsd(4)+(kpp_1d_fields%X(k,2)-Xo(k,2))*(kpp_1d_fields%X(k,2)-Xo(k,2))*&
                kpp_const_fields%hm(k)/kpp_const_fields%dm(NZ)
        ENDDO
        DO k=1,4
           rmsd(k)=SQRT(rmsd(k))
           IF (rmsd(k).ge.rmsd_threshold(k)) THEN
              kpp_1d_fields%comp_flag=.TRUE.
              kpp_1d_fields%f=kpp_1d_fields%f*1.01
           ENDIF
        ENDDO
     ENDIF
     kpp_1d_fields%reset_flag=kpp_1d_fields%reset_flag+1
     IF (kpp_1d_fields%reset_flag .gt. comp_iter_max) THEN
        WRITE(6,*) 'Failed to find a reasonable solution in the semi-implicit integration after ',&
             comp_iter_max,' iterations.'
        WRITE(6,*) 'At point lat = ',kpp_1d_fields%dlat,' lon =',kpp_1d_fields%dlon,&
             ' ipt = ',kpp_1d_fields%point,':'
        !WRITE(6,*) 'U = ',kpp_1d_fields%U(:,1)
        !WRITE(6,*) 'V = ',kpp_1d_fields%U(:,2)
        !WRITE(6,*) 'T = ',kpp_1d_fields%X(:,1)
        !WRITE(6,*) 'S = ',kpp_1d_fields%X(:,2)            
     ENDIF
  ENDDO
  ! End of trapping code.
  
  ! Output  Results from permanent grid iterations to common.inc
  ! Compute diagnostic fluxes for writing to dat file
  do k=1,NZ
     deltaz = 0.5*(kpp_const_fields%hm(k)+kpp_const_fields%hm(k+1))
     do n=1,NSCLR
        kpp_1d_fields%wX(k,n)=-kpp_1d_fields%difs(k)*((kpp_1d_fields%X(k,n)-kpp_1d_fields%X(k+1,n))/deltaz-&
             kpp_1d_fields%ghat(k)*kpp_1d_fields%wX(0,n))
     enddo
     if(kpp_const_fields%LDD) &
          kpp_1d_fields%wX(k,1)= -kpp_1d_fields%dift(k)*((kpp_1d_fields%X(k,1)-&
          kpp_1d_fields%X(k+1,1))/deltaz-kpp_1d_fields%ghat(k)*kpp_1d_fields%wX(0,1))
     kpp_1d_fields%wX(k,nsp1)= kpp_const_fields%grav * (kpp_1d_fields%talpha(k)*kpp_1d_fields%wX(k,1) - &
          kpp_1d_fields%sbeta(k) * kpp_1d_fields%wX(k,2))
     do n=1,NVEL
        kpp_1d_fields%wU(k,n)= -kpp_1d_fields%difm(k)*(kpp_1d_fields%U(k,n)-kpp_1d_fields%U(k+1,n))/deltaz
     enddo
  enddo
  ! Compute energetics
  rhonot = 1026.
  Eflx = 0.5 * ( (Uo(1,1) + kpp_1d_fields%U(1,1)) * &
       kpp_1d_fields%sflux(1,5,0) + (Uo(1,2) + kpp_1d_fields%U(1,2)) * &
       kpp_1d_fields%sflux(2,5,0) )
  Esnk = -0.5*rhonot* ( (Uo(NZ,1) + kpp_1d_fields%U(NZ,1)) * kpp_1d_fields%wU(NZ,1) + &
       (Uo(NZ,2) + kpp_1d_fields%U(NZ,2)) * kpp_1d_fields%wU(NZ,2) )
  Ptke = 0.0
  ! use "amax1" to prevent "underflow" in single precision
  do k=1,NZ-1
     Ptke = Ptke - 0.5*( amax1(kpp_1d_fields%wU(k,1),1.E-10)* &
          (rhonot   * (Uo(k,1) + kpp_1d_fields%U(k,1)) - rhonot * &
          (Uo(k+1,1) + kpp_1d_fields%U(k+1,1)) ) + amax1(kpp_1d_fields%wU(k,2),1.E-10)*&
          (rhonot   * (Uo(k,2) + kpp_1d_fields%U(k,2)) - rhonot   * (Uo(k+1,2) + kpp_1d_fields%U(k+1,2)) ) )
  ENDDO
  Tmke = 0.0
  do k=1,NZP1
     rmke(k) = 0.5 * rhonot * (kpp_1d_fields%U(k,1)**2 +  kpp_1d_fields%U(k,2)**2) * kpp_const_fields%hm(k)
     Tmke = Tmke + rmke(k)
  ENDDO
            
!     check heat and salt budgets
!     call budget(Xo,kpp_1d_fields,kpp_const_fields)
            
!     Set new profiles
!            do k=1,NZP1         ! values at NZP1 only change for slab ocean
!               do n=1,NVEL
!                  kpp_1d_fields%U(k,n) = kpp_1d_fields%U(k,n)
!               enddo
!               do n=1,NSCLR
!                  kpp_1d_fields%X(k,n) = kpp_1d_fields%X(k,n)
!               enddo
!            enddo
!     Set correct surface values, and CP and rho profiles for new profiles
!     Get latest profiles
!     Removed to ensure that hmix,diff,cp,rho diagnostics are consistent with 
!     those used in the final iteration of the timestep 
!     (CP,RHO are not quite correct for updated values, but are the values
!     by the integration) (SJW 16/01/03)
!     Need to consider improving convergence test!!! (SJW 16/01/03)
                 
!     The final call to vmix is removed to ensure that the diffusion and 
!     boundary layer profiles in the diagnostics are the ones used to calculate
!     the fluxes, as it stands at the moment this means that the CP and rho are
!     also the values used in the timestepping not the values appropriate to the
!     S,T at the new time level.  
!     call vmix(Un,Xn,hmixe,kmixe)  
!     write(40,*) time,iter, hmixn,hmixe,kmixn,kmixe
  kpp_1d_fields%hmix = hmixn
  kpp_1d_fields%kmix = kmixn
  kpp_1d_fields%uref = kpp_1d_fields%U(1,1)
  kpp_1d_fields%vref = kpp_1d_fields%U(1,2)
  kpp_1d_fields%Tref = kpp_1d_fields%X(1,1)
  IF (kpp_const_fields%L_SSref) THEN
     kpp_1d_fields%Ssurf=kpp_1d_fields%SSref
  ELSE
     kpp_1d_fields%Ssurf=kpp_1d_fields%X(1,2)+kpp_1d_fields%Sref
  ENDIF
  
  ! Damping currents, Added LH (06/08/2013)
  IF (kpp_const_fields%L_DAMP_CURR) THEN 
     dampu(:)=0.
     do k=1,NZP1
        do l=1,NVEL
           a=0.99*ABS(kpp_1d_fields%U(k,l))
           b=kpp_1d_fields%U(k,l)**2/(kpp_const_fields%dt_uvdamp*(86400./kpp_const_fields%dto))
           Ui=MIN(a,b)
           ! LH (29/08/2013) Add Flags to check which Ui (a or b) is chosen, 
           ! dtuvdamp=360 (specified in namelist). 
           ! The flags for u and v can be requested as diagnostics dampu_flag, 
           ! dampv_flag (singout 11,12). Note that the value of the flag is equal to 
           ! the *fraction* of levels at that point where (U**2)/r .lt. alpha*ABS(U), 
           ! 1.0=all Ui are (U**2)/r
           IF (b .lt. a) THEN  
              dampU(l)=dampU(l)+1.0/REAL(NZP1)
           ENDIF
           
           ! Apply damping
           kpp_1d_fields%U(k,l)= kpp_1d_fields%U(k,l) - SIGN(Ui,kpp_1d_fields%U(k,l))
        enddo
     enddo
     kpp_1d_fields%dampu_flag=dampU(1)
     kpp_1d_fields%dampv_flag=dampU(2)
  ENDIF

  ! Save variables for next timestep
  kpp_1d_fields%old = kpp_1d_fields%new
  kpp_1d_fields%new = 1 - kpp_1d_fields%old
  kpp_1d_fields%hmixd(kpp_1d_fields%new) = kpp_1d_fields%hmix
  do k=1,NZP1
     do l=1,NVEL
        kpp_1d_fields%Us(k,l,kpp_1d_fields%new)=kpp_1d_fields%U(k,l)
     enddo
     do l=1,NSCLR
        kpp_1d_fields%Xs(k,l,kpp_1d_fields%new)=kpp_1d_fields%X(k,l)
     enddo
  enddo
  
  ! close(40+ntime)
  return
end SUBROUTINE mckpp_physics_ocnstep
