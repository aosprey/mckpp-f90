SUBROUTINE mckpp_physics_ocnint(kpp_1d_fields,kpp_const_fields,intri,kmixe,Uo,Xo)
#ifdef MCKPP_CAM3
  USE mckpp_parameters
  USE mckpp_types, only: kpp_1d_type,kpp_const_type
#else 
  USE mckpp_data_fields
#endif
  IMPLICIT NONE

  ! Integrate the ocn model by backwards Euler(implicit)discretization
  ! On input : Un,Xn are estimated profiles which are used
  !            to estimate diffusivity profiles at new time.
  !          : Updated diffusivities from Un Xn are in common
  ! On output: Un,Xn are new profiles after integration.

  ! Written  19 March 1991 - jan

  INTEGER,parameter :: nuout=6,nuerr=0
  
  ! Input
  integer intri             ! index for tri.diag. coeff
  REAL Uo(NZP1,NVEL),Xo(NZP1,NSCLR)
  
  ! Output
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  ! Local
  real cu (NZtmax), &! upper coeff for (k-1) on k line of trid.matrix
       cc (NZtmax), &! central ...     (k  ) ..
       cl (NZtmax), &! lower .....     (k-1) ..
       rhs(NZtmax)  ! right-hand-side terms
  real diff(0:NZtmax),gcap(NZtmax),ntflx(0:NZtmax,NSCLR)
  ! More local variables to make implicit none
  integer kmixe,i,npd,imode,n,k
  real ftemp,ghatflux,sturflux
  integer adv_mode
  real adv_mag

  ! ********************************************************************
  ! U and V solution of tridiagonal matrix
  !                set f = 0 for equatorial application
  ftemp = kpp_1d_fields%f

  DO k=0,NZtmax
     diff(k)=kpp_1d_fields%difm(k)
  ENDDO
  call MCKPP_PHYSICS_SOLVERS_TRIDCOF(diff,NZ,intri,cu,cc,cl,kpp_const_fields)
  
  ! U right hand side and solution
  rhs(1)= Uo(1,1) + kpp_const_fields%dto* ( ftemp*.5*(Uo(1,2)+kpp_1d_fields%U(1,2)) - &
       kpp_1d_fields%wU(0,1)/kpp_const_fields%hm(1)) 
  do i=2,NZ-1
     rhs(i)= Uo(i,1) + kpp_const_fields%dto*ftemp*.5*(Uo(i,2)+kpp_1d_fields%U(i,2))      
  enddo
  i=NZ                      ! bottom
  rhs(i)= Uo(i,1) + kpp_const_fields%dto*ftemp*.5*(Uo(i,2)+ kpp_1d_fields%U(i,2)) + &
       kpp_const_fields%tri(i,1,intri)*kpp_1d_fields%difm(i)*Uo(i+1,1)
  call MCKPP_PHYSICS_SOLVERS_TRIDMAT(cu,cc,cl,rhs,Uo(:,1),NZ,kpp_1d_fields%U(:,1))
  
  ! V rhs and solution
  rhs(1)= Uo(1,2) - kpp_const_fields%dto* ( ftemp*.5*(Uo(1,1)+kpp_1d_fields%U(1,1)) + &
       kpp_1d_fields%wU(0,2)/kpp_const_fields%hm(1))
  do i=2,NZ-1
     rhs(i)= Uo(i,2) - kpp_const_fields%dto*ftemp*.5*(Uo(i,1)+kpp_1d_fields%U(i,1)) 
  enddo
  i=NZ
  rhs(i)= Uo(i,2) - kpp_const_fields%dto*ftemp*.5*(Uo(i,1)+kpp_1d_fields%U(i,1)) + &
       kpp_const_fields%tri(i,1,intri)*kpp_1d_fields%difm(i)*Uo(i+1,2)
  npd = 1
  call MCKPP_PHYSICS_SOLVERS_TRIDMAT(cu,cc,cl,rhs,Uo(:,2),NZ,kpp_1d_fields%U(:,2))
 
! *******************************************************************
! Scalar solutions of tridiagonal matrix
!     Temperature (different from other scalars because of ghat-term
!                  and double diffusion)
!     ghatflux = wX(0,1) - (1-SWDK(-hmixe,real(time)))
!    $                     * sflux(3,5,0) / rho(ipt,0) / CP(ipt,0)
!     ghatflux = wX(0,1) - (1-SWDK(-d(1) ,real(time)))
!    $                     * sflux(3,5,0) / rho(ipt,0) / CP(ipt,0)
  ghatflux = kpp_1d_fields%wX(0,1)
  sturflux = kpp_1d_fields%wX(0,1)
  diff(0)=kpp_1d_fields%dift(0)
  ntflx(0,1)=kpp_1d_fields%wXNT(0,1)
  DO k=1,NZtmax
     diff(k)=kpp_1d_fields%dift(k)
     gcap(k)=kpp_1d_fields%ghat(k)
     ntflx(k,1)=kpp_1d_fields%wXNT(k,1)
  ENDDO
  call MCKPP_PHYSICS_SOLVERS_TRIDCOF(diff,NZ,intri,cu,cc,cl,kpp_const_fields)
  
  call MCKPP_PHYSICS_SOLVERS_TRIDRHS(npd,kpp_const_fields%hm,Xo(:,1),ntflx(:,1),diff,gcap,&
       sturflux,ghatflux,kpp_const_fields%dto,NZ,intri,rhs,kpp_const_fields)  
  
  ! Surface relaxation is incompatible with flux corrections at depth (NPK 12/02/08).
  IF (kpp_const_fields%L_RELAX_SST .AND. .NOT. kpp_const_fields%L_FCORR_WITHZ .AND. .NOT. &
       kpp_const_fields%L_FCORR) THEN
     
     ! Relax the Mixed layer temperature back to SST0
     ! By using a flux correction at the surface
     ! Added by SJW (06/04)
     IF (kpp_1d_fields%relax_sst .GT. 1.e-10) THEN
        IF (.NOT. kpp_const_fields%L_RELAX_CALCONLY) THEN
           rhs(1)=rhs(1)+kpp_const_fields%dto*kpp_1d_fields%relax_sst*&
                (kpp_1d_fields%SST0-Xo(1,1))*kpp_const_fields%dm(kmixe)/kpp_const_fields%hm(1)
        ENDIF
        kpp_1d_fields%fcorr=kpp_1d_fields%relax_sst*(kpp_1d_fields%SST0-Xo(1,1))*&
             kpp_const_fields%dm(kmixe)*kpp_1d_fields%rho(1)*kpp_1d_fields%cp(1)
     ELSE
        kpp_1d_fields%fcorr=0.0
     ENDIF
  ENDIF
      
  ! Relax the mixed-layer temperature by using a USER-SPECIFIED
  ! flux correction at the surface!  Requires setting L_FCORR
  ! in the namelist and related options for input file and update
  ! frequency. Values are stored in fcorr_twod (two-dimensional flux correction)
  ! Added by NPK (29/6/08)
  IF (kpp_const_fields%L_FCORR .AND. .NOT. kpp_const_fields%L_RELAX_SST .AND. &
       .NOT. kpp_const_fields%L_FCORR_WITHZ) THEN
     rhs(1)=rhs(1)+kpp_const_fields%dto*kpp_1d_fields%fcorr_twod/&
          (kpp_1d_fields%rho(1)*kpp_1d_fields%cp(1)*kpp_const_fields%hm(1))
  ENDIF
     
  ! Correct the temperature at each layer in the model by using a flux
  ! correction applied directly to each layer.  Requires a 
  ! three-dimensional (x,y,z) input file of flux corrections.
  ! Input flux corrections must be in units of W/m-3.
  ! Added by NPK (12/2/08)
   
  ! Surface relaxation is incompatible with flux corrections at depth (NPK 12/02/08).
  kpp_1d_fields%tinc_fcorr(:)=0.
  IF (kpp_const_fields%L_FCORR_WITHZ .AND. .NOT. kpp_const_fields%L_FCORR) THEN
     DO k=1,NZP1
        kpp_1d_fields%tinc_fcorr(k) = kpp_const_fields%dto*kpp_1d_fields%fcorr_withz(k)/&
             (kpp_1d_fields%rho(k)*kpp_1d_fields%cp(k))
     ENDDO
  ENDIF
  
  ! Relax the temperature at each layer in the model by computing a flux correction 
  ! at each layer.  Requires a three-dimensional (x,y,z) input file of ocean temperatures 
  ! via subroutine MCKPP_READ_TEMPERATURES_3D.
  IF (kpp_const_fields%L_RELAX_OCNT) THEN
     DO k=1,NZP1
        ! Store the relaxation term as tinc_fcorr so that, on output,
        ! that field contains the actual correction applied in K/timestep.
        kpp_1d_fields%tinc_fcorr(k)=kpp_1d_fields%tinc_fcorr(k)+&
             kpp_const_fields%dto*kpp_1d_fields%relax_ocnT*&
             (kpp_1d_fields%ocnT_clim(k)-Xo(k,1))   
     ENDDO
  ENDIF
  DO k=1,NZP1
     rhs(k) = rhs(k) + kpp_1d_fields%tinc_fcorr(k)
     ! Modify the correction field so that, when output, it is in
     ! the correct units to be input as a flux correction via
     ! L_FCORR_WITHZ (see above).         
     kpp_1d_fields%ocnTcorr(k)=kpp_1d_fields%tinc_fcorr(k)*kpp_1d_fields%rho(k)*&
          kpp_1d_fields%cp(k)/kpp_const_fields%dto
  ENDDO

  call MCKPP_PHYSICS_SOLVERS_TRIDMAT(cu,cc,cl,rhs,Xo(:,1),NZ,kpp_1d_fields%X(:,1))

  ! Salinity and other scalars
  DO k=0,NZtmax
     diff(k)=kpp_1d_fields%difs(k)
  ENDDO
  call MCKPP_PHYSICS_SOLVERS_TRIDCOF(diff,NZ,intri,cu,cc,cl,kpp_const_fields)
  do n=2,NSCLR
     DO k=0,NZtmax
        ntflx(k,n)=kpp_1d_fields%wXNT(k,n)
     ENDDO
     ghatflux = kpp_1d_fields%wX(0,n) 
     sturflux = kpp_1d_fields%wX(0,n)
     call MCKPP_PHYSICS_SOLVERS_TRIDRHS(npd,kpp_const_fields%hm,Xo(:,n),ntflx(:,n),diff,gcap,sturflux,&
          ghatflux,kpp_const_fields%dto,NZ,intri,rhs,kpp_const_fields)
     
     ! modify rhs for advections
     do imode=1,kpp_1d_fields%nmodeadv(2)
        adv_mode=kpp_1d_fields%modeadv(imode,2)
        adv_mag=kpp_1d_fields%advection(imode,2)
        call MCKPP_PHYSICS_SOLVERS_RHSMOD(2,adv_mode,adv_mag,kpp_const_fields%dto,kmixe,&
             kpp_const_fields%dm(kmixe),NZ,rhs,kpp_1d_fields,kpp_const_fields)
     enddo
     
     ! -----Added by LH (28/05/2013) modified NPK (4/7/13)   
     if (n .eq. 2) then
        kpp_1d_fields%sinc_fcorr(:)=0.
        ! Surface salinity  relaxation is incompatible with flux corrections at depth.
        IF (kpp_const_fields%L_SFCORR_WITHZ .AND. .NOT. kpp_const_fields%L_SFCORR) THEN
           DO k=1,NZP1
              kpp_1d_fields%sinc_fcorr(k) = kpp_const_fields%dto*kpp_1d_fields%sfcorr_withz(k)
           ENDDO
        ENDIF
        
        ! Relax the salinity at each layer in the model by computing a flux correction 
        ! at each layer.  Requires a three-dimensional(x,y,z) input file of salinity 
        ! via subroutine MCKPP_READ_SALINITY_3D.
        
        IF (kpp_const_fields%L_RELAX_SAL) THEN
           DO k=1,NZP1
              ! Store the relaxation term as sinc_fcorr so that, on output,
              ! that field contains the actual correction applied in psu/timestep.
              kpp_1d_fields%sinc_fcorr(k) = kpp_1d_fields%sinc_fcorr(k)+&
                   kpp_const_fields%dto*kpp_1d_fields%relax_sal*(kpp_1d_fields%sal_clim(k)-Xo(k,n))
           ENDDO
        ENDIF
        DO k=1,NZP1
           rhs(k)=rhs(k)+kpp_1d_fields%sinc_fcorr(k)
           ! Modify the correction field so that, when output, it is in
           ! the correct units to be input as a flux correction via
           ! L_SFCORR_WITHZ (see above).
           kpp_1d_fields%scorr(k)=kpp_1d_fields%sinc_fcorr(k)/kpp_const_fields%dto
        ENDDO
     ENDIF
     ! -----end of Added LH 28/05/2013   
         
     call MCKPP_PHYSICS_SOLVERS_TRIDMAT(cu,cc,cl,rhs,Xo(:,n),NZ,kpp_1d_fields%X(:,n))
  ENDDO

  RETURN
END SUBROUTINE mckpp_physics_ocnint
