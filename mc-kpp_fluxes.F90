SUBROUTINE MCKPP_FLUXES(kpp_3d_fields,kpp_const_fields)

  IMPLICIT NONE

#include <mc-kpp_3d_type.com>
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  REAL(8), DIMENSION(NPTS) :: taux, tauy, swf, lwf, lhf, shf, rain, snow
  INTEGER :: ipt

  IF (.NOT. kpp_const_fields%L_FLUXDATA) THEN
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(shared) PRIVATE(ipt)
!$OMP DO SCHEDULE(static)
#endif
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
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
  ELSE
    CALL mckpp_read_fluxes(taux,tauy,swf,lwf,lhf,shf,rain,snow,&
             kpp_3d_fields,kpp_const_fields)
  ENDIF
      
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(shared) PRIVATE(ipt,kpp_1d_fields)
!$OMP DO SCHEDULE(dynamic)
#endif
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
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
  
END SUBROUTINE MCKPP_FLUXES
