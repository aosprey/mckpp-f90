#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif
SUBROUTINE mckpp_physics_solvers_tridcof(diff,nzi,ind,cu,cc,cl,kpp_const_fields)
#ifdef MCKPP_CAM3
  USE mckpp_parameters
  USE mckpp_types, only: kpp_const_type
#else 
  USE mckpp_data_types
#endif

  ! Compute coefficients for tridiagonal matrix (dimension=nzi).
  !     Note: cu(1) = 0. and cl(nzi) = 0. are necessary conditions.
  
  IMPLICIT NONE

  ! Input
  TYPE(kpp_const_type) :: kpp_const_fields
  integer nzi,&             ! dimension of field
       ind                  ! index for tri-coefficients: = kmixo for t-grid,
                            !                             =     1 for p-grid.
  real diff(0:nzi) ! diffusivity profile on interfaces
  
  ! Output
  real cu(nzi),&    ! upper coeff. for (k-1) on k line of trid.matrix
       cc(nzi),&    ! central ...      (k  ) ..
       cl(nzi)      ! lower .....      (k-1) ..
  integer i

  ! In the surface layer
  cu(1) = 0.
  cc(1) = 1. + kpp_const_fields%tri(1,1,ind)*diff(1)   ! 1.+ dto/h(1)/dzb(1)*diff(1)
  cl(1) =    - kpp_const_fields%tri(1,1,ind)*diff(1)   !   - dto/h(1)/dzb(1)*diff(1)

  ! Inside the domain
  do i=2,nzi
     cu(i) =    - kpp_const_fields%tri(i,0,ind)*diff(i-1)
     cc(i) = 1. + kpp_const_fields%tri(i,1,ind)*diff(i) + kpp_const_fields%tri(i,0,ind)*diff(i-1)
     cl(i) =    - kpp_const_fields%tri(i,1,ind)*diff(i)
  ENDDO

  !In the bottom layer
  cl(nzi)= 0.

  return
end SUBROUTINE mckpp_physics_solvers_tridcof

SUBROUTINE mckpp_physics_solvers_tridrhs(npd,h,yo,ntflux,diff,ghat,sturflux,ghatflux,&
     dto,nzi,ind,rhs,kpp_const_fields)

#ifdef MCKPP_CAM3
  USE mckpp_parameters
  USE mckpp_types, only: kpp_const_type
#else 
  USE mckpp_data_types
#endif

  ! Compute right hand side of tridiagonal matrix for scalar fields:
  !  =  yo (old field) 
  !     + flux-divergence of ghat
  !     + flux-divergence of non-turbulant fluxes
  ! Note: surface layer needs +dto/h(1) * surfaceflux
  ! bottom  ..... ..... +dto/h(nzi)*diff(nzi)/dzb(nzi)*yo(nzi+1)

  IMPLICIT NONE
  
  !  Input
  TYPE(kpp_const_type) :: kpp_const_fields
  real dto            ! timestep interval (seconds)
  integer nzi,&       ! dimension of field
           ind        ! index for tri-coefficients:=kmixo for t-grid,
                      !                            =1 for p-grid
  real h(nzi+1),&     ! layer thickness
       yo(nzi+1),&    ! old profile
       ntflux(0:nzi),&! non-turbulent flux = wXNT(0:nzi,1:2)
       diff(0:nzi),&  ! diffusivity profile on interfaces
       ghat(nzi),&    ! ghat turbulent flux   
       sturflux,&     ! surface turbulent (kinematic) flux = wX(0,n)
       ghatflux       ! surface flux for ghat: includes solar flux      
  integer npd         ! included in list by sjw for implicit none

  ! Output
  real rhs(nzi)      ! right hand side

  ! more local variables to make implicit none
  integer i
  real divflx
  divflx =  1.0 / float(npd)
  
  ! In the surface layer (dto/h(1)=tri(0,1,ind)
  rhs(1)= yo(1) + dto/h(1) * (ghatflux*diff(1)*ghat(1)- sturflux*divflx + &
       ntflux(1) - ntflux( 0 ) )
  
  ! Inside the domain to npd
  if(npd.ge.2) then
     do    i=2,npd
        rhs(i)= yo(i) + dto/h(i) * &
             (                        ghatflux *  diff(i)  *ghat(i) &
             -                       ghatflux *  diff(i-1)*ghat(i-1) &
             - sturflux * divflx &
             + ntflux(i) - ntflux(i-1) ) 
     end do
  endif

  ! Inside the rest of the domain
  do i=npd+1,nzi-1
     rhs(i)= yo(i) + dto/h(i) * ( ghatflux*(diff(i)*ghat(i) - diff(i-1)*ghat(i-1)) &
          +ntflux(i) - ntflux(i-1) )
  ENDDO
      
  !     In the bottom layer     
  if(nzi.gt.1) then   ! not for slab ocean
     i=nzi
     rhs(i)= yo(i) + dto/h(i) * ( ghatflux*(diff(i)*ghat(i) - diff(i-1)*ghat(i-1)) + &
          ntflux(i) - ntflux(i-1) ) + yo(i+1)*kpp_const_fields%tri(i,1,ind)*diff(i)
  endif
  
  RETURN
END SUBROUTINE mckpp_physics_solvers_tridrhs

SUBROUTINE mckpp_physics_solvers_tridmat(cu,cc,cl,rhs,yo,nzi,yn)

  USE mckpp_parameters

  ! Solve tridiagonal matrix for new vector yn, given right hand side
  ! vector rhs. Note: yn(nzi+1) = yo(nzi+1).
  IMPLICIT NONE

  ! Input
  integer nzi               ! dimension of matrix
  real cu (nzi),&            ! upper coeff. for (k-1) on k line of tridmatrix
       cc (nzi),&            ! central ...      (k  ) ..
       cl (nzi),&            ! lower .....      (k-1) ..
       rhs(nzi),&           ! right hand side
       yo(nzi+1),yni            ! old field

  ! Output
  real yn(nzi+1)    ! new field
  
  ! Local 
  real gam(NZtmax),& ! temporary array for tridiagonal solver
       bet           ! ...
  ! more local for implicit none
  integer i

  ! Solve tridiagonal matrix.
  bet   = cc(1)
  yn(1) =  rhs(1) / bet    ! surface
  DO i=2,nzi
     gam(i)= cl(i-1)/bet
     bet   = cc(i) - cu(i)*gam(i)
     if(bet.eq.0.) then
        write(nuerr,*)'* algorithm for solving tridiag matrix fails'
        write(nuerr,*)'* bet=',bet
        write(nuerr,*)'*i-1=',i-1,' cc=',cc(i-1),'cl=',cl(i-1)
        write(nuerr,*)'*i=',i,' cc=',cc(i),' cu=',cu(i),' gam=',gam(i)
        CALL MCKPP_ABORT
        bet=1.E-12
        !     Pause 3
     endif
     ! to avoid "Underflow" at single precision on the sun
     yn(i) =      (rhs(i)  - cu(i)  *yn(i-1)  )/bet     
  ENDDO

  do i=nzi-1,1,-1
     yn(i)  = yn(i) - gam(i+1)*yn(i+1)
  ENDDO
  yn(nzi+1) = yo(nzi+1)
  
  RETURN
END SUBROUTINE mckpp_physics_solvers_tridmat

subroutine mckpp_physics_solvers_rhsmod(jsclr,mode,A,dto,km,dm,nzi,rhs,kpp_1d_fields,kpp_const_fields)
  
#ifdef MCKPP_CAM3
  USE mckpp_parameters
  USE mckpp_types, only: kpp_1d_type,kpp_const_type
#else
  USE mckpp_data_types
#endif

!     Modify rhs to correct scalar, jsclr, 
!     for advection according to mode
! mode = 1 : Steady upper layer horizontal advection
!        2 : Steady mixed layer horizontal advection to km-1 
!        3 : Steady horizontal advection throughout the entire column
!        4 : Steady vertical advection (= deep horizontal) below 100m 
!            to bottom 
!            (Change: start below 100m, instead of at layer 16, and
!            do not advect into bottom layer, 7-1-93)
!        5 : Steady bottom diffusion
!        6 : Seasonal mixed layer horizontal advection to dm
!        7 : Seasonal thermocline horizontal advection to 1.5 dm
  IMPLICIT NONE

  ! Input
  integer nzi,&              ! vertical dimension of field
       km,&                  ! index of gridpoint just below h
       mode,&                ! type of advection
       jsclr                ! scalar
  real rhs(nzi)             ! right hand side from tridrhs
  real dto,&                 ! ocean time step
       dm,&                  ! depth d(km+.5)
       A                    ! advection of heat(W/m2) or Salt(PSU m/s)      
  
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  ! Output
  real am,fact,delta,depth,dmax
  integer n,n1,n2,nzend
  
  ! Internal
  real f(12)     ! monthly partion of annual advection
  real xsA(21)   ! yearly excess of heat

  !data f/.1,.1,6*0.0,.1,.3,.4,.2/
  data f/.05,.05,5*0.0,.05,.15,.20,.30,.20/
  data xsA/48.26,21.73,29.02,56.59,19.94,15.96,18.28,&
           40.52,37.06,29.83,29.47,15.77, 1.47,14.55,&
           4.22,28.19,39.54,19.58,20.27,11.19,21.72/

  if(mode.le.0) return 
!               find ocean year
!        iyr = 1 + idint((time-75.0)/dpy)
!        day = time - dpy * (idint(time/dpy))  ! 365.25))
!        month = 1 + int(12. * day / dpy)    ! 367.)
!        if(month.gt.12) then
!           write(nuerr,*) 'STOP rhsmod (ocn.f):'
!           write(nuerr,*) '     rounding error, month gt 12 =',month
!           stop 97
!        endif

!       Am = -12. * f(month) * (xsA(iyr) - 0.0 )       ! Annual
!        Am =  12. * f(month) * A                       ! Seasonal
  Am = A                                          ! Steady
  
  if(mode.eq.1) then
     ! correct upper layer advection
     if(jsclr.eq.1) fact = dto * Am / (kpp_1d_fields%rho(1)*kpp_1d_fields%cp(1))
     if(jsclr.eq.2) fact = dto * Am * 0.033
     rhs(1) = rhs(1) + fact / kpp_const_fields%hm(1)
     
  else if(mode.eq.2) then
     ! correct mixed layer advection
     delta = 0.0
     do n=1,km-1
        delta = delta + kpp_const_fields%hm(n)
     enddo
     do n=1,km-1
        if(jsclr.eq.1) fact = dto * Am / (kpp_1d_fields%rho(n)*kpp_1d_fields%cp(n)) 
        if(jsclr.eq.2) fact = dto * Am * 0.033
        rhs(n) = rhs(n) + fact  / delta        
     ENDDO
     
  else if (mode.eq.3) then
     ! throughout whole water column
     delta = 0.0
     do n=1,nzi
        delta = delta + kpp_const_fields%hm(n)
     enddo
     do n=1,nzi
        if(jsclr.eq.1) fact = dto * Am / (kpp_1d_fields%rho(n)*kpp_1d_fields%cp(n))
        if(jsclr.eq.2) fact = dto * Am * 0.033
        rhs(n) = rhs(n) + fact / delta 
     ENDDO
     
  else if (mode.eq.4) then
     ! vertical advection = deep horizontal
     nzend=nzi-1            ! nzend=nzi (change:7-1-93)
     n1=0                   ! n1=16     (change:7-1-93)
401  n1=n1+1
     if(kpp_const_fields%zm(n1).ge.-100.) goto 401
     delta = 0.0
     do n=n1,nzend
        delta = delta + kpp_const_fields%hm(n)
     enddo
     do n=n1,nzend
        if(jsclr.eq.1) fact = dto * Am / (kpp_1d_fields%rho(n)*kpp_1d_fields%cp(n))
        if(jsclr.eq.2) fact = dto * Am * 0.033
        rhs(n) = rhs(n) + fact / delta 
     ENDDO
     
  else if(mode.eq.5) then
     ! correct bottom layer diffusion
     if(jsclr.eq.1) fact = dto * Am / (kpp_1d_fields%rho(nzi)*kpp_1d_fields%cp(nzi))
     if(jsclr.eq.2) fact = dto * Am * 0.033
     rhs(nzi) = rhs(nzi) + fact / kpp_const_fields%hm(nzi)         
     
  else             
     !     seasonal mixed layer or thermocline advection  
     !     find ocean year
     !     iyr = 1 + idint((time-75.0)/dpy)
     !     day = time - dpy * (idint(time/dpy))  ! 365.25))
     !     month = 1 + int(12. * day / dpy)    ! 367.)
     !     diag
     !     if(month.gt.12) then
     !     write(nuerr,*) 'STOP rhsmod (ocn.f):'
     !     write(nuerr,*) '     rounding error, month gt 12 =',month
     !     stop 97
     !     endif  
     !     diag
     !     Am = -12. * f(month) * (xsA(iyr) - 0.0 )       ! Annual
     !     Am =  12. * f(month) * A                       ! Seasonal
     !     Am = A                                          ! Steady
     if(mode.eq.6) then
        ! mixed layer to dm
        n1 = 1
        depth = kpp_const_fields%hm(1)
        dmax  = dm -  0.5 * (kpp_const_fields%hm(km) + kpp_const_fields%hm(km-1))
        delta = 0.0
        do n =n1,nzi
           n2    = n
           delta = delta + kpp_const_fields%hm(n)
           depth = depth + kpp_const_fields%hm(n+1)
           if(depth.ge.dmax) go to 606
        ENDDO
606     continue
        
     else if (mode.eq.7) then
        !     thermocline to 100m   
        n1 = km - 1
        depth = dm - 0.5 * kpp_const_fields%hm(km) 
        dmax = 100.        
        delta = 0.0
        do n=n1,nzi
           n2 = n
           delta = delta + kpp_const_fields%hm(n)
           depth = depth + kpp_const_fields%hm(n+1)
           if(depth.ge.dmax) go to 706
        ENDDO
706     continue
        
     else
        write(nuerr,*) 'STOP in rhsmod (ocn.f):'
        write(nuerr,*) '      mode out of range, mode=',mode
        CALL MCKPP_ABORT
     endif
     
     !     Finish both 6 and 7 here
     do n=n1,n2  
        if(jsclr.eq.1) fact = dto * Am / (kpp_1d_fields%rho(n)*kpp_1d_fields%cp(n)) 
        if(jsclr.eq.2) fact = dto * Am * 0.033
        rhs(n) = rhs(n) + fact  / delta
     ENDDO

  endif
  
  return
end subroutine mckpp_physics_solvers_rhsmod
