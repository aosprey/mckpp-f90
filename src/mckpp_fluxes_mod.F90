MODULE mckpp_fluxes_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields, kpp_1d_type, kpp_const_type
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: npts, nz
  USE mckpp_read_fluxes_mod, ONLY: mckpp_read_fluxes
  USE mckpp_types_transfer, ONLY: mckpp_fields_3dto1d, mckpp_fields_1dto3d
  
  IMPLICIT NONE

CONTAINS 

SUBROUTINE MCKPP_FLUXES()

  TYPE(kpp_1d_type) :: kpp_1d_fields
  REAL(8), DIMENSION(NPTS) :: taux, tauy, swf, lwf, lhf, shf, rain, snow
  INTEGER :: ipt

  IF (.NOT. kpp_const_fields%L_FLUXDATA) THEN
!#ifdef OPENMP
!!$OMP PARALLEL DEFAULT(shared) PRIVATE(ipt)
!!$OMP DO SCHEDULE(static)
!#endif
    DO ipt=1,npts
      taux(ipt)=0.01
      tauy(ipt)=0.0
      swf(ipt)=200.0
      lwf(ipt)=0.0
      lhf(ipt)=-150.0
      shf(ipt)=0.0
      rain(ipt)=6e-5
      snow(ipt)=0.0
    ENDDO   
!#ifdef OPENMP
!!$OMP END DO
!!$OMP END PARALLEL
!#endif
  ELSE
    CALL mckpp_read_fluxes(taux,tauy,swf,lwf,lhf,shf,rain,snow)
  ENDIF
      
!#ifdef OPENMP
!!$OMP PARALLEL DEFAULT(shared) PRIVATE(ipt,kpp_1d_fields)
!!$OMP DO SCHEDULE(dynamic)
!#endif
  DO ipt=1,npts         
    IF (kpp_3d_fields%L_OCEAN(ipt)) THEN 
      IF ((taux(ipt) .EQ. 0.0) .AND. (tauy(ipt) .EQ. 0.0)) THEN
        taux(ipt)=1.e-10
      ENDIF
      IF (.NOT. kpp_const_fields%L_REST) THEN
        kpp_3d_fields%sflux(ipt,1,5,0)=taux(ipt)
        kpp_3d_fields%sflux(ipt,2,5,0)=tauy(ipt)
        kpp_3d_fields%sflux(ipt,3,5,0)=swf(ipt)            
        kpp_3d_fields%sflux(ipt,4,5,0)=lwf(ipt)+lhf(ipt)+shf(ipt)-snow(ipt)*kpp_const_fields%FLSN
        kpp_3d_fields%sflux(ipt,5,5,0)=1e-10 ! Melting of sea-ice = 0.0
        kpp_3d_fields%sflux(ipt,6,5,0)=(rain(ipt)+snow(ipt)+(lhf(ipt)/kpp_const_fields%EL))
      ELSE
        kpp_3d_fields%sflux(ipt,1,5,0)=1.e-10
        kpp_3d_fields%sflux(ipt,2,5,0)=0.00
        kpp_3d_fields%sflux(ipt,3,5,0)=300.00
        kpp_3d_fields%sflux(ipt,4,5,0)=-300.00
        kpp_3d_fields%sflux(ipt,5,5,0)=0.00
        kpp_3d_fields%sflux(ipt,6,5,0)=0.00
      ENDIF
      CALL mckpp_fields_3dto1d(kpp_3d_fields,ipt,kpp_1d_fields)
      CALL mckpp_fluxes_ntflux(kpp_1d_fields,kpp_const_fields)
      CALL mckpp_fields_1dto3d(kpp_1d_fields,ipt,kpp_3d_fields)
    ENDIF
  ENDDO
!#ifdef OPENMP
!!$OMP END DO
!!$OMP END PARALLEL
!#endif
  
END SUBROUTINE MCKPP_FLUXES


SUBROUTINE mckpp_fluxes_ntflux(kpp_1d_fields,kpp_const_fields)
  
  INTEGER k
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  CHARACTER(LEN=31) :: routine = "MCKPP_FLUXES_NTFLUX"
  CHARACTER(LEN=max_message_len) :: message
  
  ! WRITE(message,*) "At time = ", kpp_const_fields%ntime
  ! CALL mckpp_print(routine, message)
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

END SUBROUTINE mckpp_fluxes_ntflux


REAL FUNCTION MCKPP_FLUXES_SWDK(z,jerlov)
  REAL :: z
  INTEGER :: j, jerlov
  
  INTEGER, PARAMETER :: max=5
  real Rfac(max),a1(max),a2(max)
!         types =  I       IA      IB      II      III
!             j =  1       2       3       4       5
  Rfac = (/  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /)
  a1 = (/  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /)
  a2 = (/ 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /)
  
  j = jerlov
  MCKPP_FLUXES_SWDK = Rfac(j) * dexp(dble(z/a1(j))) + (1.0-Rfac(j)) * dexp(dble(z/a2(j)))

END FUNCTION MCKPP_FLUXES_SWDK


END MODULE mckpp_fluxes_mod
