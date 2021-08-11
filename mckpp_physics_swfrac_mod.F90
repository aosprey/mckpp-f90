MODULE mckpp_physics_swfrac_mod

CONTAINS 

subroutine MCKPP_PHYSICS_SWFRAC_OPT(fact, kpp_1d_fields, kpp_const_fields)

#ifdef MCKPP_CAM3
  USE mckpp_types, only : kpp_1d_type, kpp_const_type
#else
  USE mckpp_data_fields, ONLY: kpp_1d_type, kpp_const_type
#endif
  USE mckpp_parameters, ONLY: nzp1

  ! compute fraction of solar short-wave flux penetrating to specified
  ! depth (times fact) due to exponential decay in  Jerlov water type
  ! reference : two band solar absorption model of simpson and 
  ! paulson (1977)
      
  IMPLICIT NONE
  
  integer nwtype
  parameter(nwtype=5) ! max number of different water types 

  ! Input
  real fact                 ! scale  factor to apply to depth array
  type(kpp_1d_type) :: kpp_1d_fields
  type(kpp_const_type) :: kpp_const_fields
  
  real  rfac(nwtype),a1(nwtype),a2(nwtype)
  real rmin,r1,r2
  integer l
  
  ! jerlov water type :  I       IA      IB      II      III
  !            jwtype    1       2       3       4       5
  !
  data rfac         /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
  data a1           /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
  data a2           / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /
  data rmin         / -80. /
  
  DO l=1,NZP1
     r1      = MAX(kpp_const_fields%zm(l)*fact/a1(kpp_1d_fields%jerlov), rmin)
     r2      = MAX(kpp_const_fields%zm(l)*fact/a2(kpp_1d_fields%jerlov), rmin)
     kpp_1d_fields%swfrac(l) = rfac(kpp_1d_fields%jerlov)  * &
          exp(r1) + (1.-rfac(kpp_1d_fields%jerlov)) * exp(r2) 
  ENDDO
  
end subroutine MCKPP_PHYSICS_SWFRAC_OPT


subroutine MCKPP_PHYSICS_SWFRAC(fact, z, jwtype, swdk)

! compute fraction of solar short-wave flux penetrating to specified
! depth (times fact) due to exponential decay in  Jerlov water type
! reference : two band solar absorption model of simpson and paulson (1977)

  IMPLICIT NONE
  
  integer nwtype
  parameter(nwtype=5) ! max number of different water types 

  ! input
  real fact      ! scale factor to apply to depth array
  real z         ! vertical height ( <0.) for desired sw 
  !                           fraction                                 (m)
  integer jwtype ! index for jerlov water type

  !  output
  real swdk      !  short wave (radiation) fractional decay

  !  local
  real  rfac(nwtype),a1(nwtype),a2(nwtype)
  real rmin,r1,r2 

  ! jerlov water type :  I       IA      IB      II      III
  !            jwtype    1       2       3       4       5
  data rfac         /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
  data a1           /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
  data a2           / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /
  data rmin         / -80. /
  
        
  r1   = MAX(z*fact/a1(jwtype), rmin)
  r2   = MAX(z*fact/a2(jwtype), rmin)
  swdk = rfac(jwtype)  * exp(r1) + (1.-rfac(jwtype)) * exp(r2)

END SUBROUTINE MCKPP_PHYSICS_SWFRAC

END MODULE mckpp_physics_swfrac_mod
