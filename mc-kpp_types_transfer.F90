SUBROUTINE mckpp_fields_3dto1d(kpp_fields_3d,point,kpp_fields_1d)
#ifdef MCKPP_CAM3
  USE mckpp_parameters
  USE mckpp_types, only: kpp_3d_type,kpp_1d_type
#else 
  USE mckpp_data_fields
#endif
  IMPLICIT NONE

!     Accepts a 3D variable of the KPP derived type.
!     Returns a 1D variable of the KPP derived type, extraced from the 3D variable
!     at a specified point.

  TYPE(kpp_3d_type),intent(in)  :: kpp_fields_3d
  TYPE(kpp_1d_type),intent(out) :: kpp_fields_1d
  INTEGER,intent(in) :: point
  INTEGER :: i,j,k
  REAL :: temp
  LOGICAL :: logical_temp

  CALL mckpp_allocate_1d_fields(kpp_fields_1d)

  DO i=1,NZP1         
     DO j=1,NVEL
        temp=kpp_fields_3d%U(point,i,j)
        kpp_fields_1d%U(i,j)=temp
        temp=kpp_fields_3d%U_init(point,i,j) ! Not updated within physics
        kpp_fields_1d%U_init(i,j)=temp
        DO k=0,1
           temp=kpp_fields_3d%Us(point,i,j,k)
           kpp_fields_1d%Us(i,j,k)=temp
        ENDDO
     ENDDO
     DO j=1,NSCLR
        temp=kpp_fields_3d%X(point,i,j)
        kpp_fields_1d%X(i,j)=temp
        DO k=0,1
           temp=kpp_fields_3d%Xs(point,i,j,k)
           kpp_fields_1d%Xs(i,j,k)=temp
        ENDDO
     ENDDO
     temp=kpp_fields_3d%Rig(point,i)
     kpp_fields_1d%Rig(i)=temp
     temp=kpp_fields_3d%Shsq(point,i)
     kpp_fields_1d%Shsq(i)=temp
     temp=kpp_fields_3d%swfrac(point,i)
     kpp_fields_1d%swfrac(i)=temp
     temp=kpp_fields_3d%tinc_fcorr(point,i)
     kpp_fields_1d%tinc_fcorr(i)=temp
     temp=kpp_fields_3d%sinc_fcorr(point,i)
     kpp_fields_1d%sinc_fcorr(i)=temp
     temp=kpp_fields_3d%fcorr_withz(point,i)
     kpp_fields_1d%fcorr_withz(i)=temp      
     temp=kpp_fields_3d%sfcorr_withz(point,i)
     kpp_fields_1d%sfcorr_withz(i)=temp
     temp=kpp_fields_3d%sal_clim(point,i) !Not updated within physics
     kpp_fields_1d%sal_clim(i)=temp
     temp=kpp_fields_3d%ocnT_clim(point,i) !Not updated within physics
     kpp_fields_1d%ocnT_clim(i)=temp
     temp=kpp_fields_3d%ocnTcorr(point,i)
     kpp_fields_1d%ocnTcorr(i)=temp
     IF (i .le. NZ) THEN
        temp=kpp_fields_3d%dbloc(point,i)
        kpp_fields_1d%dbloc(i)=temp
     ENDIF
  ENDDO

  DO i=0,1
     temp=kpp_fields_3d%hmixd(point,i)
     kpp_fields_1d%hmixd(i)=temp
  ENDDO
  DO i=0,NZP1tmax
     temp=kpp_fields_3d%rho(point,i)
     kpp_fields_1d%rho(i)=temp
     temp=kpp_fields_3d%cp(point,i)
     kpp_fields_1d%cp(i)=temp         
     IF (i .gt. 0) THEN
        temp=kpp_fields_3d%buoy(point,i)
        kpp_fields_1d%buoy(i)=temp
     ENDIF
     IF (i .le. NZtmax) THEN
        temp=kpp_fields_3d%difm(point,i)
        kpp_fields_1d%difm(i)=temp
        temp=kpp_fields_3d%difs(point,i)
        kpp_fields_1d%difs(i)=temp
        temp=kpp_fields_3d%dift(point,i)
        kpp_fields_1d%dift(i)=temp
        DO j=1,NVP1
           temp=kpp_fields_3d%wU(point,i,j)
           kpp_fields_1d%wU(i,j)=temp
        ENDDO
        DO j=1,NSP1
           temp=kpp_fields_3d%wX(point,i,j)
           kpp_fields_1d%wX(i,j)=temp
        ENDDO
        DO j=1,NSCLR
           temp=kpp_fields_3d%wXNT(point,i,j)
           kpp_fields_1d%wXNT(i,j)=temp
        ENDDO
        IF (i .gt. 0) THEN
           temp=kpp_fields_3d%ghat(point,i)
           kpp_fields_1d%ghat(i)=temp
        ENDIF
     ENDIF
     IF (i .le. NZ) THEN
        temp=kpp_fields_3d%swdk_opt(point,i)
        kpp_fields_1d%swdk_opt(i)=temp
     ENDIF
  ENDDO
  
  DO i=1,2
     k=kpp_fields_3d%nmodeadv(point,i) ! Not updated within physics
     kpp_fields_1d%nmodeadv(i)=k
     DO j=1,maxmodeadv
        k=kpp_fields_3d%modeadv(point,j,i) !Not updated within physics
        kpp_fields_1d%modeadv(j,i)=k
        temp=kpp_fields_3d%advection(point,j,i) !Not updated within physics
        kpp_fields_1d%advection(j,i)=temp
     ENDDO
  ENDDO
  
  DO i=1,NSFLXS
     DO j=1,5
        DO k=0,NJDT
           temp=kpp_fields_3d%sflux(point,i,j,k)  !Not updated within physics
           kpp_fields_1d%sflux(i,j,k)=temp
        ENDDO
     ENDDO
  ENDDO
  
  temp=kpp_fields_3d%ocdepth(point) !Not updated within physics
  kpp_fields_1d%ocdepth=temp      
  logical_temp=kpp_fields_3d%L_OCEAN(point) !Not updated within physics
  kpp_fields_1d%L_OCEAN=logical_temp
  logical_temp=kpp_fields_3d%L_INITFLAG(point)
  kpp_fields_1d%L_INITFLAG=logical_temp
  temp=kpp_fields_3d%f(point) !Not updated within physics
  kpp_fields_1d%f=temp
  
  temp=kpp_fields_3d%relax_sst(point) !Not updated within physics
  kpp_fields_1d%relax_sst=temp
  temp=kpp_fields_3d%fcorr(point)
  kpp_fields_1d%fcorr=temp
  temp=kpp_fields_3d%fcorr_twod(point)
  kpp_fields_1d%fcorr_twod=temp
  temp=kpp_fields_3d%SST0(point)      !Not updated within physics
  kpp_fields_1d%SST0=temp
  temp=kpp_fields_3d%relax_sal(point) !Not updated within physics
  kpp_fields_1d%relax_sal=temp
  temp=kpp_fields_3d%relax_ocnt(point) !Not updated within physics
  kpp_fields_1d%relax_ocnt=temp
  
  temp=kpp_fields_3d%hmix(point)
  kpp_fields_1d%hmix=temp
  temp=kpp_fields_3d%kmix(point)
  kpp_fields_1d%kmix=temp
  temp=kpp_fields_3d%Tref(point)
  kpp_fields_1d%Tref=temp
  temp=kpp_fields_3d%uref(point)
  kpp_fields_1d%uref=temp
  temp=kpp_fields_3d%vref(point)
  kpp_fields_1d%vref=temp
  temp=kpp_fields_3d%Ssurf(point)
  kpp_fields_1d%Ssurf=temp
  temp=kpp_fields_3d%Sref(point) !Not updated within physics
  kpp_fields_1d%Sref=temp
  temp=kpp_fields_3d%SSref(point) !Not updated within physics
  kpp_fields_1d%SSref=temp
  
  i=kpp_fields_3d%old(point)
  kpp_fields_1d%old=i
  i=kpp_fields_3d%new(point)
  kpp_fields_1d%new=i
  i=kpp_fields_3d%jerlov(point) ! Not updated within physics
  kpp_fields_1d%jerlov=i
  
  temp=kpp_fields_3d%dlat(point) !Not updated within physics
  kpp_fields_1d%dlat=temp          
  temp=kpp_fields_3d%dlon(point) !Not updated within physics
  kpp_fields_1d%dlon=temp
  temp=kpp_fields_3d%cplwght(point)
  kpp_fields_1d%cplwght=temp
  
  kpp_fields_1d%point=point
  
  RETURN
END SUBROUTINE mckpp_fields_3dto1d

SUBROUTINE mckpp_fields_1dto3d(kpp_fields_1d,point,kpp_fields_3d)
#ifdef MCKPP_CAM3
  USE mckpp_parameters
  USE mckpp_types, only: kpp_1d_type,kpp_3d_type
#else 
  USE mckpp_data_fields
#endif
  IMPLICIT NONE
  
  ! Accepts a 1D and a 3D variable of the KPP derived type.
  ! Returns the 3D variable, updated at a specified point with the 
  ! values from the 1D variable.

  TYPE(kpp_3d_type),intent(inout)  :: kpp_fields_3d
  TYPE(kpp_1d_type),intent(in) :: kpp_fields_1d
  INTEGER,intent(in) :: point
  INTEGER :: i,j,k
  REAL :: temp
  LOGICAL :: logical_temp
  
  DO i=1,NZP1
     DO j=1,NVEL
        temp=kpp_fields_1d%U(i,j)
        kpp_fields_3d%U(point,i,j)=temp
        DO k=0,1
           temp=kpp_fields_1d%Us(i,j,k)
           kpp_fields_3d%Us(point,i,j,k)=temp
        ENDDO
     ENDDO
     DO j=1,NSCLR
        temp=kpp_fields_1d%X(i,j)
        kpp_fields_3d%X(point,i,j)=temp
        DO k=0,1
           temp=kpp_fields_1d%Xs(i,j,k)
           kpp_fields_3d%Xs(point,i,j,k)=temp
        ENDDO
     ENDDO
     temp=kpp_fields_1d%Rig(i)
     kpp_fields_3d%Rig(point,i)=temp
     temp=kpp_fields_1d%Shsq(i)
     kpp_fields_3d%Shsq(point,i)=temp
     temp=kpp_fields_1d%swfrac(i)
     kpp_fields_3d%swfrac(point,i)=temp
     temp=kpp_fields_1d%tinc_fcorr(i)
     kpp_fields_3d%tinc_fcorr(point,i)=temp
     temp=kpp_fields_1d%sinc_fcorr(i)
     kpp_fields_3d%sinc_fcorr(point,i)=temp
     temp=kpp_fields_1d%fcorr_withz(i)
     kpp_fields_3d%fcorr_withz(point,i)=temp
     temp=kpp_fields_1d%scorr(i)
     kpp_fields_3d%scorr(point,i)=temp
     temp=kpp_fields_1d%sfcorr_withz(i)
     kpp_fields_3d%sfcorr_withz(point,i)=temp
     temp=kpp_fields_1d%ocnTcorr(i)
     kpp_fields_3d%ocnTcorr(point,i)=temp
     IF (i.le.NZ) THEN 
        temp=kpp_fields_1d%dbloc(i)
        kpp_fields_3d%dbloc(point,i)=temp
     ENDIF
  ENDDO
  DO i=0,1
     temp=kpp_fields_1d%hmixd(i)
     kpp_fields_3d%hmixd(point,i)=temp
  ENDDO
  DO i=0,NZP1tmax
     temp=kpp_fields_1d%rho(i)
     kpp_fields_3d%rho(point,i)=temp
     temp=kpp_fields_1d%cp(i)
     kpp_fields_3d%cp(point,i)=temp         
     IF (i .gt. 0) THEN
        temp=kpp_fields_1d%buoy(i)
        kpp_fields_3d%buoy(point,i)=temp
     ENDIF
     IF (i .le. NZtmax) THEN
        temp=kpp_fields_1d%difm(i)
        kpp_fields_3d%difm(point,i)=temp
        temp=kpp_fields_1d%difs(i)
        kpp_fields_3d%difs(point,i)=temp
        temp=kpp_fields_1d%dift(i)
        kpp_fields_3d%dift(point,i)=temp
        DO j=1,NVP1
           temp=kpp_fields_1d%wU(i,j)
           kpp_fields_3d%wU(point,i,j)=temp
        ENDDO
        DO j=1,NSP1
           temp=kpp_fields_1d%wX(i,j)
           kpp_fields_3d%wX(point,i,j)=temp
        ENDDO
        DO j=1,NSCLR
           temp=kpp_fields_1d%wXNT(i,j)
           kpp_fields_3d%wXNT(point,i,j)=temp
        ENDDO
        IF (i .gt. 0) THEN 
           temp=kpp_fields_1d%ghat(i)
           kpp_fields_3d%ghat(point,i)=temp
        ENDIF
     ENDIF
     IF (i .le. NZ) THEN
        temp=kpp_fields_1d%swdk_opt(i)
        kpp_fields_3d%swdk_opt(point,i)=temp
     ENDIF
  ENDDO
  
  logical_temp=kpp_fields_1d%L_INITFLAG
  kpp_fields_3d%L_INITFLAG(point)=logical_temp
  temp=kpp_fields_1d%freeze_flag
  kpp_fields_3d%freeze_flag(point)=temp
  
  temp=kpp_fields_1d%fcorr
  kpp_fields_3d%fcorr(point)=temp
  temp=kpp_fields_1d%fcorr_twod
  kpp_fields_3d%fcorr_twod(point)=temp
  
  temp=kpp_fields_1d%hmix
  kpp_fields_3d%hmix(point)=temp
  temp=kpp_fields_1d%kmix
  kpp_fields_3d%kmix(point)=temp
  temp=kpp_fields_1d%Tref
  kpp_fields_3d%Tref(point)=temp
  temp=kpp_fields_1d%uref
  kpp_fields_3d%uref(point)=temp
  temp=kpp_fields_1d%vref
  kpp_fields_3d%vref(point)=temp
  temp=kpp_fields_1d%Ssurf
  kpp_fields_3d%Ssurf(point)=temp

  i=kpp_fields_1d%old
  kpp_fields_3d%old(point)=i
  i=kpp_fields_1d%new
  kpp_fields_3d%new(point)=i
  
  temp=kpp_fields_1d%reset_flag
  kpp_fields_3d%reset_flag(point)=temp
  
  temp=kpp_fields_1d%dampu_flag
  kpp_fields_3d%dampu_flag(point)=temp
  temp=kpp_fields_1d%dampv_flag
  kpp_fields_3d%dampv_flag(point)=temp
  
  RETURN
END SUBROUTINE mckpp_fields_1dto3d
