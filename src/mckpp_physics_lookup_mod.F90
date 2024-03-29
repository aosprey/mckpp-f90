MODULE mckpp_physics_lookup_mod

  USE mckpp_data_fields, ONLY: kpp_const_type
  USE mckpp_log_messages, ONLY: mckpp_print

  IMPLICIT NONE

CONTAINS

! This is only called once at initialization 
SUBROUTINE mckpp_physics_lookup(kpp_const_fields)
  
  TYPE(kpp_const_type) :: kpp_const_fields
  real zmin,zmax,umin,umax,usta,zeta,zehat,epsln,&
       am,cm,c1,c2,zetam,as,cs,c3,zetas,cstar,deltau,deltaz
  integer i,j,ni,nj
  CHARACTER(LEN=20) :: routine = "MCKPP_PHYSICS_LOOKUP"

  CALL mckpp_print(routine, "") 

  ni    = 890     ! number of values for zehat
  nj    = 48      ! number of values for ustar      
  epsln = 1.e-20
  c1    = 5.0  
  zmin  = -4.e-7
  zmax  = 0.0     ! m3/s3
  umin  = 0.0
  umax  = 0.04    ! m/s
  am    = 1.257
  cm    = 8.380
  c2    = 16.0
  zetam = -0.2    !7-24-92
  as    = -28.86
  cs    = 98.96
  c3    = 16.0
  zetas = -1.0
  cstar = 5.0 
  
  deltaz = (zmax-zmin)/(ni+1) 
  deltau = (umax-umin)/(nj+1)
  
  DO i=0,ni+1
     zehat = deltaz*(i) + zmin
     DO j=0,nj+1
        usta = deltau*(j) + umin
        zeta = zehat/(usta**3+epsln)
        
        if(zehat.ge.0.) then
           kpp_const_fields%wmt(i,j) = kpp_const_fields%vonk*usta/(1.+c1*zeta)
           kpp_const_fields%wst(i,j) = kpp_const_fields%wmt(i,j)
        else
           if(zeta.gt.zetam) then
              kpp_const_fields%wmt(i,j) = kpp_const_fields%vonk * usta * (1.-c2*zeta)**(1./4.)
           else
              kpp_const_fields%wmt(i,j) = kpp_const_fields%vonk * (am*usta**3 - cm*zehat)**(1./3.)
           endif
           if(zeta.gt.zetas) then
              kpp_const_fields%wst(i,j) = kpp_const_fields%vonk * usta * (1.-c3*zeta)**(1./2.)
           else
              kpp_const_fields%wst(i,j) = kpp_const_fields%vonk * (as*usta**3 - cs*zehat)**(1./3.)
           endif
        endif
     ENDDO
  ENDDO
      
END SUBROUTINE mckpp_physics_lookup

END MODULE mckpp_physics_lookup_mod
