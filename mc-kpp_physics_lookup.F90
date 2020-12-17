SUBROUTINE mckpp_physics_lookup(kpp_const_fields)
#ifdef MCKPP_CAM3
  USE mckpp_parameters
  USE mckpp_types, only : kpp_const_type
#else 
  USE mckpp_data_fields
#endif
  IMPLICIT NONE

  TYPE(kpp_const_type) :: kpp_const_fields
  
  real zmin,zmax,umin,umax,usta,zeta,zehat,epsln,&
       am,cm,c1,c2,zetam,as,cs,c3,zetas,cstar,deltau,deltaz
  integer i,j,ni,nj
  parameter ( ni = 890,&     ! number of values for zehat
       nj = 48)             ! number of values for ustar      
  data epsln             /   1.e-20 /
  data c1                /   5.0    /
  data zmin,zmax  / -4.e-7, 0.0   / ! m3/s3
  data umin,umax  /  0.   , .04   / ! m/s
  data am,cm,c2,zetam    /   1.257  ,  8.380, 16.0, - 0.2 / !7-24-92
  data as,cs,c3,zetas    / -28.86   , 98.96 , 16.0, - 1.0 /
  data cstar             /    5.    /
  
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
      
  RETURN
END SUBROUTINE mckpp_physics_lookup
