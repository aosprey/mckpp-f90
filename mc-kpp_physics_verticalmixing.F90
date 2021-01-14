SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING(kpp_1d_fields,kpp_const_fields,hmixn,kmixn)
#ifdef MCKPP_CAM3
  USE mckpp_parameters
  USE mckpp_types, only: kpp_1d_type,kpp_const_type
#else 
  USE mckpp_data_types
#endif
  !  Interface between 1-d model and vertical mixing
  IMPLICIT NONE

! inputs 
  type(kpp_1d_type) :: kpp_1d_fields
  type(kpp_const_type) :: kpp_const_fields
  real B0,B0sol,ustar

! outputs 
  real hmixn                ! boundary layer depth (m)
  integer kmixn       
  real rhob
  
! local
  real dVsq(nzp1)    ! (velocity shear re sfc)^2      (m/s)^2
  real Ritop(nz)     ! numerator of bulk Richardson Number (m/s)^2
  real alphaDT(nz)   ! alpha * DT across interfaces
  real betaDS(nz)    ! beta  * DS across interfaces
  real epsilon,epsln
  real alpha,beta,exppr
  real sigma,sigma0
  real MCKPP_CPSW
  real tau
  real zref,wz,bref
  integer k,n,kl
  real del
  real dlimit,vlimit
  integer jerl(12)
  ! month  1   2   3   4   5   6   7   8   9   10  11  12
  data jerl / 2 , 2 , 2 , 3 , 3 , 3 , 4 , 4 , 4 , 4 , 3 , 2 /
  
  epsilon = 0.1
  epsln   = 1.e-20

  ! calculate density of fresh water and brine in surface layer
  alpha = 1.
  beta  = 1.
  exppr = 0.0
  sigma0=0
  sigma=0
  call MCKPP_ABK80(0.0,kpp_1d_fields%X(1,1),-kpp_const_fields%zm(1),alpha,beta,exppr,sigma0,sigma)
  kpp_1d_fields%rhoh2o = 1000. + sigma0
  call MCKPP_ABK80(kpp_const_fields%SICE,kpp_1d_fields%X(1,1),-kpp_const_fields%zm(1),alpha,beta,exppr,sigma0,sigma)
  rhob = 1000. + sigma0
 
  ! calculate temperature and salt contributions of buoyancy gradients  
  ! calculate buoyancy profile (m/s**2) on gridlevels
  do k=1,nzp1
     call MCKPP_ABK80(kpp_1d_fields%X(k,2)+kpp_1d_fields%Sref,kpp_1d_fields%X(k,1),-kpp_const_fields%zm(k),&
          alpha,beta,exppr,sigma0,sigma)
     kpp_1d_fields%rho(k)= 1000. + sigma0
     kpp_1d_fields%CP(k) = MCKPP_CPSW(kpp_1d_fields%X(k,2)+kpp_1d_fields%Sref,kpp_1d_fields%X(k,1),&
          -kpp_const_fields%zm(k))
     kpp_1d_fields%talpha(k) = alpha
     kpp_1d_fields%sbeta(k)  = beta
     kpp_1d_fields%buoy(k) = -kpp_const_fields%grav * sigma0 / 1000.
  ENDDO
  
  kpp_1d_fields%rho(0) = kpp_1d_fields%rho(1) 
  kpp_1d_fields%CP(0) = kpp_1d_fields%CP(1)   
  kpp_1d_fields%talpha(0) = kpp_1d_fields%talpha(1)   
  kpp_1d_fields%sbeta(0)  = kpp_1d_fields%sbeta(1)           

  ! Call to ntflx, put here to allow removal of diagnostic call to vmix
  ! and to ensure the most recent cp,rho used (consistent with other 
  ! surface fluxes?)
  call MCKPP_FLUXES_NTFLUX(kpp_1d_fields,kpp_const_fields)

  ! calculate kinematic surface momentum fluxes
  kpp_1d_fields%wU(0,1) = -kpp_1d_fields%sflux(1,5,0) / kpp_1d_fields%rho(0)
  kpp_1d_fields%wU(0,2) = -kpp_1d_fields%sflux(2,5,0) / kpp_1d_fields%rho(0)
  tau = sqrt( kpp_1d_fields%sflux(1,5,0)**2 + kpp_1d_fields%sflux(2,5,0)**2 ) + 1.e-16
  !  1.e-16 added to stop subsequent division by zero if tau=0.0
  ustar = sqrt( tau / kpp_1d_fields%rho(0) )
 
  ! total turbulent kinematic temperature flux (C m/s)
  kpp_1d_fields%wX(0,1)  = -kpp_1d_fields%sflux(4,5,0) / kpp_1d_fields%rho(0) / kpp_1d_fields%CP(0)
 
  ! total turbulent kinematic salinity flux (o/oo m/s)
  kpp_1d_fields%wX(0,2) = kpp_1d_fields%Ssurf*kpp_1d_fields%sflux(6,5,0)/&
       kpp_1d_fields%rhoh2o+(kpp_1d_fields%Ssurf-kpp_const_fields%SICE)*&
       kpp_1d_fields%sflux(5,5,0)/rhob

  ! calculate total kinematic surface buoyancy flux (m**2/s**3)
  B0 = -kpp_const_fields%grav*(kpp_1d_fields%talpha(0)*kpp_1d_fields%wX(0,1) - &
       kpp_1d_fields%sbeta(0)*kpp_1d_fields%wX(0,2) )
  kpp_1d_fields%wX(0,NSP1) =  - B0
  B0sol = kpp_const_fields%grav * kpp_1d_fields%talpha(0) * kpp_1d_fields%sflux(3,5,0) / &
       (kpp_1d_fields%rho(0) * kpp_1d_fields%CP(0))

  ! calculate temperature and salt contributions of buoyancy gradients on interfaces for double diffusion      
  do n = 1,nz
     alphaDT(n) =0.5 *(kpp_1d_fields%talpha(n)+kpp_1d_fields%talpha(n+1)) * &
          (kpp_1d_fields%X(n,1) - kpp_1d_fields%X(n+1,1))
     betaDS(n)  =0.5 *(kpp_1d_fields%sbeta(n) + kpp_1d_fields%sbeta(n+1)) * &
          (kpp_1d_fields%X(n,2) - kpp_1d_fields%X(n+1,2))
  ENDDO
      
  ! compute buoyancy and shear profiles
  DO n = 1,nz
     zref =  epsilon * kpp_const_fields%zm(n)
     ! compute reference buoyancy and velocity
     wz = AMAX1(kpp_const_fields%zm(1),zref) 
     kpp_1d_fields%uref  = kpp_1d_fields%U(1,1) * wz / zref
     kpp_1d_fields%vref  = kpp_1d_fields%U(1,2) * wz / zref
     bref  = kpp_1d_fields%buoy(1)* wz / zref
     do kl = 1,nz
        IF(zref.ge.kpp_const_fields%zm(kl)) go to 126
        wz = AMIN1(kpp_const_fields%zm(kl)-kpp_const_fields%zm(kl+1),kpp_const_fields%zm(kl)-zref) 
        del = 0.5 * wz / (kpp_const_fields%zm(kl) - kpp_const_fields%zm(kl+1))
        kpp_1d_fields%uref = kpp_1d_fields%uref - wz*( kpp_1d_fields%U(kl,1) + del * &
             (kpp_1d_fields%U(kl+1,1)- kpp_1d_fields%U(kl,1))) /zref
        kpp_1d_fields%vref = kpp_1d_fields%vref - wz*( kpp_1d_fields%U(kl,2) + del * &
             (kpp_1d_fields%U(kl+1,2)- kpp_1d_fields%U(kl,2))) /zref
        bref=bref -wz*(kpp_1d_fields%buoy(kl) + del * &
             (kpp_1d_fields%buoy(kl+1)-kpp_1d_fields%buoy(kl))) /zref
     ENDDO
126  continue
     Ritop(n) = (zref - kpp_const_fields%zm(n)) * (bref - kpp_1d_fields%buoy(n))
     ! NPK Additions (25/9/2008). Prevent Ritop from going negative.
     ! IF (Ritop(ipt,n) .lt. 0) Ritop(ipt,n) = epsln
     kpp_1d_fields%dbloc(n) = kpp_1d_fields%buoy(n) - kpp_1d_fields%buoy(n+1)
     dVsq(n)  = (kpp_1d_fields%Uref - kpp_1d_fields%U(n,1))**2 + (kpp_1d_fields%Vref - kpp_1d_fields%U(n,2))**2
     kpp_1d_fields%shsq(n)  = (kpp_1d_fields%U(n,1)-kpp_1d_fields%U(n+1,1))**2 + &
          (kpp_1d_fields%U(n,2)-kpp_1d_fields%U(n+1,2))**2         
  ENDDO
     
  CALL MCKPP_PHYSICS_VERTICALMIXING_KPPMIX(nz,nzp1,dVsq,ustar,B0,B0sol,alphaDT,betaDS,Ritop,hmixn,kmixn,kpp_1d_fields,kpp_const_fields)
        
  ! limit the bottom diffusity and viscosity
  ! zero diffusivities for no bottom flux option
  ! if(LNBFLX) then 
  !      dlimit = 0.0
  !      vlimit = 0.0
  !      do n=1,nsclr
  !      wxNT(ipt,nz,n) = 0.0
  !      enddo
  ! else
  dlimit = 0.00001
  vlimit = 0.0001
  ! endif
  do k=nz,nzp1
     kpp_1d_fields%difm(k) = vlimit
     kpp_1d_fields%difs(k) = dlimit
     kpp_1d_fields%dift(k) = dlimit
  enddo
  kpp_1d_fields%ghat(nz) = 0.0

  RETURN
END SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING
