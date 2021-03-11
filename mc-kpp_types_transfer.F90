
SUBROUTINE mckpp_fields_3dto1d(kpp_3d_fields,point,kpp_1d_fields)
#ifdef MCKPP_CAM3
  USE mckpp_parameters
  USE mckpp_types, only: kpp_3d_type,kpp_1d_type
#else 
  USE mckpp_data_types
#endif
  IMPLICIT NONE

!     Accepts a 3D variable of the KPP derived type.
!     Returns a 1D variable of the KPP derived type, extraced from the 3D variable
!     at a specified point.

  TYPE(kpp_3d_type),intent(in)  :: kpp_3d_fields
  TYPE(kpp_1d_type),intent(out) :: kpp_1d_fields
  INTEGER,intent(in) :: point
  INTEGER :: i,j,k
  REAL :: temp
  LOGICAL :: logical_temp

  CALL mckpp_allocate_1d_fields(kpp_1d_fields)

  DO i=1,NZP1         
     DO j=1,NVEL
        temp=kpp_3d_fields%U(point,i,j)
        kpp_1d_fields%U(i,j)=temp
        temp=kpp_3d_fields%U_init(point,i,j) ! Not updated within physics
        kpp_1d_fields%U_init(i,j)=temp
        DO k=0,1
           temp=kpp_3d_fields%Us(point,i,j,k)
           kpp_1d_fields%Us(i,j,k)=temp
        ENDDO
     ENDDO
     DO j=1,NSCLR
        temp=kpp_3d_fields%X(point,i,j)
        kpp_1d_fields%X(i,j)=temp
        DO k=0,1
           temp=kpp_3d_fields%Xs(point,i,j,k)
           kpp_1d_fields%Xs(i,j,k)=temp
        ENDDO
     ENDDO
     temp=kpp_3d_fields%Rig(point,i)
     kpp_1d_fields%Rig(i)=temp
     temp=kpp_3d_fields%Shsq(point,i)
     kpp_1d_fields%Shsq(i)=temp
     temp=kpp_3d_fields%swfrac(point,i)
     kpp_1d_fields%swfrac(i)=temp
     temp=kpp_3d_fields%tinc_fcorr(point,i)
     kpp_1d_fields%tinc_fcorr(i)=temp
     temp=kpp_3d_fields%sinc_fcorr(point,i)
     kpp_1d_fields%sinc_fcorr(i)=temp
     temp=kpp_3d_fields%fcorr_withz(point,i)
     kpp_1d_fields%fcorr_withz(i)=temp      
     temp=kpp_3d_fields%scorr(point,i)
     kpp_1d_fields%scorr(i)=temp
     temp=kpp_3d_fields%sfcorr_withz(point,i)
     kpp_1d_fields%sfcorr_withz(i)=temp
     temp=kpp_3d_fields%sal_clim(point,i) !Not updated within physics
     kpp_1d_fields%sal_clim(i)=temp
     temp=kpp_3d_fields%ocnT_clim(point,i) !Not updated within physics
     kpp_1d_fields%ocnT_clim(i)=temp
     temp=kpp_3d_fields%ocnTcorr(point,i)
     kpp_1d_fields%ocnTcorr(i)=temp
     IF (i .le. NZ) THEN
        temp=kpp_3d_fields%dbloc(point,i)
        kpp_1d_fields%dbloc(i)=temp
     ENDIF
  ENDDO

  DO i=0,1
     temp=kpp_3d_fields%hmixd(point,i)
     kpp_1d_fields%hmixd(i)=temp
  ENDDO
  DO i=0,NZP1tmax
     temp=kpp_3d_fields%rho(point,i)
     kpp_1d_fields%rho(i)=temp
     temp=kpp_3d_fields%cp(point,i)
     kpp_1d_fields%cp(i)=temp         
     IF (i .gt. 0) THEN
        temp=kpp_3d_fields%buoy(point,i)
        kpp_1d_fields%buoy(i)=temp
     ENDIF
     IF (i .le. NZtmax) THEN
        temp=kpp_3d_fields%difm(point,i)
        kpp_1d_fields%difm(i)=temp
        temp=kpp_3d_fields%difs(point,i)
        kpp_1d_fields%difs(i)=temp
        temp=kpp_3d_fields%dift(point,i)
        kpp_1d_fields%dift(i)=temp
        DO j=1,NVP1
           temp=kpp_3d_fields%wU(point,i,j)
           kpp_1d_fields%wU(i,j)=temp
        ENDDO
        DO j=1,NSP1
           temp=kpp_3d_fields%wX(point,i,j)
           kpp_1d_fields%wX(i,j)=temp
        ENDDO
        DO j=1,NSCLR
           temp=kpp_3d_fields%wXNT(point,i,j)
           kpp_1d_fields%wXNT(i,j)=temp
        ENDDO
        IF (i .gt. 0) THEN
           temp=kpp_3d_fields%ghat(point,i)
           kpp_1d_fields%ghat(i)=temp
        ENDIF
     ENDIF
     IF (i .le. NZ) THEN
        temp=kpp_3d_fields%swdk_opt(point,i)
        kpp_1d_fields%swdk_opt(i)=temp
     ENDIF
  ENDDO
  
  DO i=1,2
     k=kpp_3d_fields%nmodeadv(point,i) ! Not updated within physics
     kpp_1d_fields%nmodeadv(i)=k
     DO j=1,maxmodeadv
        k=kpp_3d_fields%modeadv(point,j,i) !Not updated within physics
        kpp_1d_fields%modeadv(j,i)=k
        temp=kpp_3d_fields%advection(point,j,i) !Not updated within physics
        kpp_1d_fields%advection(j,i)=temp
     ENDDO
  ENDDO
  
  DO i=1,NSFLXS
     DO j=1,5
        DO k=0,NJDT
           temp=kpp_3d_fields%sflux(point,i,j,k)  !Not updated within physics
           kpp_1d_fields%sflux(i,j,k)=temp
        ENDDO
     ENDDO
  ENDDO
  
  temp=kpp_3d_fields%ocdepth(point) !Not updated within physics
  kpp_1d_fields%ocdepth=temp      
  logical_temp=kpp_3d_fields%L_OCEAN(point) !Not updated within physics
  kpp_1d_fields%L_OCEAN=logical_temp
  logical_temp=kpp_3d_fields%L_INITFLAG(point)
  kpp_1d_fields%L_INITFLAG=logical_temp
  temp=kpp_3d_fields%f(point) !Not updated within physics
  kpp_1d_fields%f=temp
  temp=kpp_3d_fields%freeze_flag(point)
  kpp_1d_fields%freeze_flag=temp
  
  temp=kpp_3d_fields%relax_sst(point) !Not updated within physics
  kpp_1d_fields%relax_sst=temp
  temp=kpp_3d_fields%fcorr(point)
  kpp_1d_fields%fcorr=temp
  temp=kpp_3d_fields%fcorr_twod(point)
  kpp_1d_fields%fcorr_twod=temp
  temp=kpp_3d_fields%SST0(point)      !Not updated within physics
  kpp_1d_fields%SST0=temp
  temp=kpp_3d_fields%relax_sal(point) !Not updated within physics
  kpp_1d_fields%relax_sal=temp
  temp=kpp_3d_fields%relax_ocnt(point) !Not updated within physics
  kpp_1d_fields%relax_ocnt=temp
  
  temp=kpp_3d_fields%hmix(point)
  kpp_1d_fields%hmix=temp
  temp=kpp_3d_fields%kmix(point)
  kpp_1d_fields%kmix=temp
  temp=kpp_3d_fields%Tref(point)
  kpp_1d_fields%Tref=temp
  temp=kpp_3d_fields%uref(point)
  kpp_1d_fields%uref=temp
  temp=kpp_3d_fields%vref(point)
  kpp_1d_fields%vref=temp
  temp=kpp_3d_fields%Ssurf(point)
  kpp_1d_fields%Ssurf=temp
  temp=kpp_3d_fields%Sref(point) !Not updated within physics
  kpp_1d_fields%Sref=temp
  temp=kpp_3d_fields%SSref(point) !Not updated within physics
  kpp_1d_fields%SSref=temp
  
  i=kpp_3d_fields%old(point)
  kpp_1d_fields%old=i
  i=kpp_3d_fields%new(point)
  kpp_1d_fields%new=i
  i=kpp_3d_fields%jerlov(point) ! Not updated within physics
  kpp_1d_fields%jerlov=i
  
  temp=kpp_3d_fields%dlat(point) !Not updated within physics
  kpp_1d_fields%dlat=temp          
  temp=kpp_3d_fields%dlon(point) !Not updated within physics
  kpp_1d_fields%dlon=temp
  temp=kpp_3d_fields%cplwght(point) !Not updated within physics
  kpp_1d_fields%cplwght=temp
  
  kpp_1d_fields%point=point
  
  RETURN
END SUBROUTINE mckpp_fields_3dto1d

SUBROUTINE mckpp_fields_1dto3d(kpp_1d_fields,point,kpp_3d_fields)
#ifdef MCKPP_CAM3
  USE mckpp_parameters
  USE mckpp_types, only: kpp_1d_type,kpp_3d_type
#else 
  USE mckpp_data_types
#endif
  IMPLICIT NONE
  
  ! Accepts a 1D and a 3D variable of the KPP derived type.
  ! Returns the 3D variable, updated at a specified point with the 
  ! values from the 1D variable.

  TYPE(kpp_3d_type),intent(inout)  :: kpp_3d_fields
  TYPE(kpp_1d_type),intent(in) :: kpp_1d_fields
  INTEGER,intent(in) :: point
  INTEGER :: i,j,k
  REAL :: temp
  LOGICAL :: logical_temp
  
  DO i=1,NZP1
     DO j=1,NVEL
        temp=kpp_1d_fields%U(i,j)
        kpp_3d_fields%U(point,i,j)=temp
        DO k=0,1
           temp=kpp_1d_fields%Us(i,j,k)
           kpp_3d_fields%Us(point,i,j,k)=temp
        ENDDO
     ENDDO
     DO j=1,NSCLR
        temp=kpp_1d_fields%X(i,j)
        kpp_3d_fields%X(point,i,j)=temp
        DO k=0,1
           temp=kpp_1d_fields%Xs(i,j,k)
           kpp_3d_fields%Xs(point,i,j,k)=temp
        ENDDO
     ENDDO
     temp=kpp_1d_fields%Rig(i)
     kpp_3d_fields%Rig(point,i)=temp
     temp=kpp_1d_fields%Shsq(i)
     kpp_3d_fields%Shsq(point,i)=temp
     temp=kpp_1d_fields%swfrac(i)
     kpp_3d_fields%swfrac(point,i)=temp
     temp=kpp_1d_fields%tinc_fcorr(i)
     kpp_3d_fields%tinc_fcorr(point,i)=temp
     temp=kpp_1d_fields%sinc_fcorr(i)
     kpp_3d_fields%sinc_fcorr(point,i)=temp
     temp=kpp_1d_fields%fcorr_withz(i)
     kpp_3d_fields%fcorr_withz(point,i)=temp
     temp=kpp_1d_fields%scorr(i)
     kpp_3d_fields%scorr(point,i)=temp
     temp=kpp_1d_fields%sfcorr_withz(i)
     kpp_3d_fields%sfcorr_withz(point,i)=temp
     temp=kpp_1d_fields%ocnTcorr(i)
     kpp_3d_fields%ocnTcorr(point,i)=temp
     IF (i.le.NZ) THEN 
        temp=kpp_1d_fields%dbloc(i)
        kpp_3d_fields%dbloc(point,i)=temp
     ENDIF
  ENDDO
  DO i=0,1
     temp=kpp_1d_fields%hmixd(i)
     kpp_3d_fields%hmixd(point,i)=temp
  ENDDO
  DO i=0,NZP1tmax
     temp=kpp_1d_fields%rho(i)
     kpp_3d_fields%rho(point,i)=temp
     temp=kpp_1d_fields%cp(i)
     kpp_3d_fields%cp(point,i)=temp         
     IF (i .gt. 0) THEN
        temp=kpp_1d_fields%buoy(i)
        kpp_3d_fields%buoy(point,i)=temp
     ENDIF
     IF (i .le. NZtmax) THEN
        temp=kpp_1d_fields%difm(i)
        kpp_3d_fields%difm(point,i)=temp
        temp=kpp_1d_fields%difs(i)
        kpp_3d_fields%difs(point,i)=temp
        temp=kpp_1d_fields%dift(i)
        kpp_3d_fields%dift(point,i)=temp
        DO j=1,NVP1
           temp=kpp_1d_fields%wU(i,j)
           kpp_3d_fields%wU(point,i,j)=temp
        ENDDO
        DO j=1,NSP1
           temp=kpp_1d_fields%wX(i,j)
           kpp_3d_fields%wX(point,i,j)=temp
        ENDDO
        DO j=1,NSCLR
           temp=kpp_1d_fields%wXNT(i,j)
           kpp_3d_fields%wXNT(point,i,j)=temp
        ENDDO
        IF (i .gt. 0) THEN 
           temp=kpp_1d_fields%ghat(i)
           kpp_3d_fields%ghat(point,i)=temp
        ENDIF
     ENDIF
     IF (i .le. NZ) THEN
        temp=kpp_1d_fields%swdk_opt(i)
        kpp_3d_fields%swdk_opt(point,i)=temp
     ENDIF
  ENDDO
  
  logical_temp=kpp_1d_fields%L_INITFLAG
  kpp_3d_fields%L_INITFLAG(point)=logical_temp
  temp=kpp_1d_fields%freeze_flag
  kpp_3d_fields%freeze_flag(point)=temp
  
  temp=kpp_1d_fields%fcorr
  kpp_3d_fields%fcorr(point)=temp
  temp=kpp_1d_fields%fcorr_twod
  kpp_3d_fields%fcorr_twod(point)=temp
  
  temp=kpp_1d_fields%hmix
  kpp_3d_fields%hmix(point)=temp
  temp=kpp_1d_fields%kmix
  kpp_3d_fields%kmix(point)=temp
  temp=kpp_1d_fields%Tref
  kpp_3d_fields%Tref(point)=temp
  temp=kpp_1d_fields%uref
  kpp_3d_fields%uref(point)=temp
  temp=kpp_1d_fields%vref
  kpp_3d_fields%vref(point)=temp
  temp=kpp_1d_fields%Ssurf
  kpp_3d_fields%Ssurf(point)=temp

  i=kpp_1d_fields%old
  kpp_3d_fields%old(point)=i
  i=kpp_1d_fields%new
  kpp_3d_fields%new(point)=i
  
  temp=kpp_1d_fields%reset_flag
  kpp_3d_fields%reset_flag(point)=temp
  
  temp=kpp_1d_fields%dampu_flag
  kpp_3d_fields%dampu_flag(point)=temp
  temp=kpp_1d_fields%dampv_flag
  kpp_3d_fields%dampv_flag(point)=temp
  
  RETURN
END SUBROUTINE mckpp_fields_1dto3d
