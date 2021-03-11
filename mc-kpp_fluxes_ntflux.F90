SUBROUTINE mckpp_fluxes_ntflux(kpp_1d_fields,kpp_const_fields)
#ifdef MCKPP_CAM3
  USE mckpp_parameters
  USE mckpp_types, only: kpp_1d_type,kpp_const_type
#else 
  USE mckpp_data_types
#endif 

  IMPLICIT NONE

  INTEGER k
  REAL MCKPP_FLUXES_SWDK
  EXTERNAL MCKPP_FLUXES_SWDK
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  !WRITE(6,*) 'MCKPP_FLUXES_NTFLUX at time = ',kpp_const_fields%ntime
  IF (kpp_const_fields%ntime .le. 1) THEN
     DO k=0,NZ
        kpp_1d_fields%swdk_opt(k)=MCKPP_FLUXES_SWDK(-kpp_const_fields%dm(k),kpp_1d_fields%jerlov)
     ENDDO
  ENDIF
  DO k=0,NZ
     IF (kpp_const_fields%ntime .ge. 1) THEN 
        kpp_1d_fields%wXNT(k,1)=-kpp_1d_fields%sflux(3,5,0)*kpp_1d_fields%swdk_opt(k)&
             /(kpp_1d_fields%rho(0)*kpp_1d_fields%CP(0))
     ENDIF
  ENDDO

  RETURN
END SUBROUTINE mckpp_fluxes_ntflux

REAL FUNCTION MCKPP_FLUXES_SWDK(z,jerlov)

  parameter(max=5)
  real Rfac(max),a1(max),a2(max)
!         types =  I       IA      IB      II      III
!             j =  1       2       3       4       5
  data Rfac /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
  data a1   /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
  data a2   / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /
  
  j = jerlov
  MCKPP_FLUXES_SWDK = Rfac(j) * dexp(dble(z/a1(j))) + (1.0-Rfac(j)) * dexp(dble(z/a2(j)))

  return
end FUNCTION MCKPP_FLUXES_SWDK

