      SUBROUTINE fluxes(kpp_3d_fields,kpp_const_fields,kpp_timer)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
      
! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
#include "landsea.com"
#ifdef COUPLE
#ifdef OASIS3
#include "kpp_oasis3.inc"
#endif /*OASIS3*/
#endif /*COUPLE*/
      include 'ocn_paras.com'
      include 'flx_paras.com'
      include 'flx_in.com'
      include 'local_pt.com'
#include "timocn.com"
#include "initialcon.com"
      include 'times.com'
      include 'couple.com'

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_2d_type) :: kpp_2d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      TYPE(kpp_timer_type) :: kpp_timer

      REAL taux(NPTS),tauy(NPTS),
     $     swf(NPTS),lwf(NPTS),lhf(NPTS),shf(NPTS),
     $     rain(NPTS),snow(NPTS),curl_tau(NPTS)
      REAL SST_in(NX_GLOBE,NY_GLOBE,1),ICE_in(NX_GLOBE,NY_GLOBE,1),
     +     icedepth_in(NX_GLOBE,NY_GLOBE,1),vsf_in(NX_GLOBE,NY_GLOBE),
     +     snowdepth_in(NX_GLOBE,NY_GLOBE,1),usf_in(NX_GLOBE,NY_GLOBE)
      CHARACTER(LEN=19) :: trans_timer_name
      COMMON /save_sstin/ SST_in,ICE_in,icedepth_in,snowdepth_in,
     +     usf_in,vsf_in
      INTEGER ix,jy,ipt_globe
      
#ifdef OPENMP
      INTEGER tid,OMP_GET_THREAD_NUM      
#endif

#ifdef COUPLE
#ifdef OASIS2
      call coupled_flux(swf,lwf,rain,ntime-1)
      call coupled_stress(taux,tauy,ntime-1)
#else
#ifdef CFS
      CALL read_gfs_forcing(swf,lwf,rain,taux,tauy)
#else
#ifdef OASIS3
!     Normal coupling - no writing to or reading from netCDF files
      CALL mpi1_oasis3_input(swf,lwf,rain,taux,tauy,curl_tau,
     +     kpp_3d_fields,kpp_const_fields)

      ! HadGEM3 passes zeros at the first timestep for a new run (i.e., NRUN)
      ! Thus, if this is NOT a restart run, we need to provide a file
      ! of fluxes for KPP for the first coupling timestep.           
      ! NPK 15/10/09, revised 6/11/09 to specify .NOT. L_RESTART
      ! as HadGEM3 does pass good fields for a restart run (i.e., CRUN)
      IF (kpp_const_fields%ntime .EQ. 1 .AND. .NOT. L_RESTART) THEN          
         WRITE(6,*) 'KPP: Reading fluxes from file ',initflux_file
         CALL init_flxdata(initflux_file,kpp_const_fields)
         CALL read_fluxes(taux,tauy,swf,lwf,lhf,shf,rain,snow,curl_tau,
     +        kpp_3d_fields,kpp_const_fields)
! Convert to variables expected for a coupled model
         DO ipt=1,npts
            IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
               lwf(ipt)=lwf(ipt)+lhf(ipt)+shf(ipt)-snow(ipt)*FLSN
            ENDIF
         ENDDO
      ENDIF
#endif /*OASIS3*/
#endif /*CFS*/
#endif /*OASIS2*/
! All coupled models do this step      
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(shared) 
!$OMP& PRIVATE(ipt,kpp_2d_fields,tid,trans_timer_name)
      tid=OMP_GET_THREAD_NUM()
      WRITE(trans_timer_name,'(A17,I2)') 'KPP 3D/2D thread ',tid
#else
      WRITE(trans_timer_name,'(A19)') 'KPP 3D/2D thread 01'
#endif
      CALL kpp_timer_time(kpp_timer,trans_timer_name,1)
#ifdef OPENMP
!$OMP DO SCHEDULE(dynamic)
#endif
      DO ipt=1,npts
c        WRITE(6,*) 'ipt=',ipt,'lwf=',lwf(ipt)
         IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
            IF ((taux(ipt) .EQ. 0.0) .AND. (tauy(ipt) .EQ. 0.0)) THEN
               taux(ipt)=1.e-10
            ENDIF
            kpp_3d_fields%sflux(ipt,1,5,0)=taux(ipt)
            kpp_3d_fields%sflux(ipt,2,5,0)=tauy(ipt)
            kpp_3d_fields%sflux(ipt,3,5,0)=swf(ipt)               
            kpp_3d_fields%sflux(ipt,4,5,0)=lwf(ipt) 
            kpp_3d_fields%sflux(ipt,5,5,0)=0.0 ! Melting of sea-ice = 0.0               
            kpp_3d_fields%sflux(ipt,6,5,0)=rain(ipt) ! assuming rain = P-E
            IF (kpp_const_fields%L_EKMAN_PUMP) THEN
               kpp_3d_fields%sflux(ipt,7,5,0)=curl_tau(ipt)
            ELSE
               kpp_3d_fields%sflux(ipt,7,5,0)=0.0
            ENDIF
            CALL kpp_fields_3dto2d(kpp_3d_fields,ipt,kpp_2d_fields) 
            
            call ntflx(kpp_2d_fields,kpp_const_fields)
            CALL kpp_fields_2dto3d(kpp_2d_fields,ipt,kpp_3d_fields)
         ENDIF
c         WRITE(6,*) 'After fluxes, ipt=',ipt,'sst=',
c     +        kpp_3d_fields%X(ipt,1,1)
      ENDDO
#ifdef OPENMP
!$OMP END DO
#endif
      CALL kpp_timer_time(kpp_timer,trans_timer_name,0)
#ifdef OPENMP
!$OMP END PARALLEL
#endif
#else  /* NOT COUPLED */
      ! Get fluxes for the forced case
      IF (.NOT. L_FLUXDATA) THEN
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
         call read_fluxes(taux,tauy,swf,lwf,lhf,shf,rain,snow,curl_tau,
     +        kpp_3d_fields,kpp_const_fields)
      ENDIF
      
c      WRITE(6,*) 'L_REST=',L_REST
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(shared) PRIVATE(ipt,kpp_2d_fields)
!$OMP DO SCHEDULE(dynamic)
#endif
      DO ipt=1,npts         
         IF (kpp_3d_fields%L_OCEAN(ipt)) THEN 
            IF ((taux(ipt) .EQ. 0.0) .AND. (tauy(ipt) .EQ. 0.0)) THEN
               taux(ipt)=1.e-10
            ENDIF
            IF (.NOT. L_REST) THEN
               kpp_3d_fields%sflux(ipt,1,5,0)=taux(ipt)
               kpp_3d_fields%sflux(ipt,2,5,0)=tauy(ipt)
               kpp_3d_fields%sflux(ipt,3,5,0)=swf(ipt)            
               kpp_3d_fields%sflux(ipt,4,5,0)=lwf(ipt)+lhf(ipt)+shf(ipt)
     +              -snow(ipt)*FLSN
               kpp_3d_fields%sflux(ipt,6,5,0)=(rain(ipt)+snow(ipt)+
     +              (lhf(ipt)/EL))
               kpp_3d_fields%sflux(ipt,5,5,0)=1e-10 ! Melting of sea-ice = 0.0
               IF (kpp_const_fields%L_EKMAN_PUMP) THEN
                  kpp_3d_fields%sflux(ipt,7,5,0)=curl_tau(ipt)
               ELSE
                  kpp_3d_fields%sflux(ipt,7,5,0)=0.0
               ENDIF
            ELSE
               kpp_3d_fields%sflux(ipt,1,5,0)=1.e-10
               kpp_3d_fields%sflux(ipt,2,5,0)=0.00
               kpp_3d_fields%sflux(ipt,3,5,0)=300.00
               kpp_3d_fields%sflux(ipt,4,5,0)=-300.00
               kpp_3d_fields%sflux(ipt,5,5,0)=0.00
               kpp_3d_fields%sflux(ipt,6,5,0)=0.00
               kpp_3d_fields%sflux(ipt,7,5,0)=0.0
            ENDIF
            CALL kpp_fields_3dto2d(kpp_3d_fields,ipt,kpp_2d_fields)
            call ntflx(kpp_2d_fields,kpp_const_fields)
            CALL kpp_fields_2dto3d(kpp_2d_fields,ipt,kpp_3d_fields)
         ENDIF
      ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
#endif /*COUPLE*/

!     If L_FCORR_NSOL is enabled, modify the non-solar heat flux by
!     a factor proportional to the SST bias against a user-specified
!     climatology.
      IF (kpp_const_fields%L_FCORR_NSOL) THEN
         DO ix=ifirst,ilast
            DO jy=jfirst,jlast
               ipt=(jy-jfirst)*nx+(ix-ifirst)+1
               ipt_globe=(jy-1)*NX_GLOBE+ix
               IF (kpp_3d_fields%L_OCEAN(ipt) 
#ifdef COUPLE
     &              .and. kpp_3d_fields%cplwght(ipt_globe) .gt. 0) THEN
#else
     &              ) THEN
#endif
                  kpp_3d_fields%fcorr_nsol(ipt) = 
     &                 kpp_3d_fields%fcorr_nsol_coeff(ipt)*
     &                 (kpp_3d_fields%X(ipt,1,1)-SST_in(ix,jy,1))
                  kpp_3d_fields%sflux(ipt,4,5,0) = 
     &                 kpp_3d_fields%sflux(ipt,4,5,0) + 
     &                 kpp_3d_fields%fcorr_nsol(ipt)
               ELSE
                  kpp_3d_fields%fcorr_nsol(ipt)=1e20
               ENDIF               
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END

************************************************************************
      SUBROUTINE ntflx(kpp_2d_fields,kpp_const_fields)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
      
      INTEGER k
      REAL SWDK, dm(0:NZ)
c      REAL SWDK_OPT(NPTS,0:NZ)
      EXTERNAL SWDK
      TYPE(kpp_2d_type) :: kpp_2d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

c      COMMON /SWDK_SAVE/SWDK_OPT

c      WRITE(6,*) 'ntime = ',kpp_const_fields%ntime
      IF (kpp_const_fields%ntime .le. 1) THEN
	 IF (kpp_const_fields%L_SLAB) THEN
	    kpp_2d_fields%swdk_opt(0)=1
	    kpp_2d_fields%swdk_opt(1)=0
            IF (kpp_const_fields%L_COLUMBIA_LAND .and.
     +           kpp_2d_fields%dlat .ge. -30 .and.
     +           kpp_2d_fields%dlat .le. 30 .and.
     +           kpp_2d_fields%dlon .ge. 0 .and.
     +           kpp_2d_fields%dlon .le. 45) THEN
               dm(1) = 0.1
            ELSE
               dm(1) = kpp_const_fields%slab_depth     
            ENDIF
         ELSE
            dm = kpp_const_fields%dm
         ENDIF
         IF (.NOT. kpp_const_fields%L_VARY_OPT) THEN
            DO k=0,NZ
               kpp_2d_fields%swdk_opt(k)=swdk(-dm(k),
     +              kpp_2d_fields%jerlov) 
            ENDDO
         ENDIF
      ENDIF
      IF (kpp_const_fields%ntime .ge. 1) THEN
         DO k=0,NZ
            kpp_2d_fields%wXNT(k,1)=-kpp_2d_fields%sflux(3,5,0)
     +           *kpp_2d_fields%swdk_opt(k)
     &           /(kpp_2d_fields%rho(0)*kpp_2d_fields%CP(0))
         ENDDO
      ENDIF
c     WRITE(6,*) 'Computed wXNT, sflux=',kpp_2d_fields%sflux(3,5,0),
c     + 'swdk_opt=',kpp_2d_fields%swdk_opt,'rho=',
c     +    kpp_2d_fields%rho(0),'cp=',kpp_2d_fields%CP(0),'jerlov=',
c     +    kpp_2d_fields%jerlov

c     DO k=0,NZ
c        wXNT(ipt,k,1)=-sflux(ipt,3,5,0)*swdk(-dm(k))
c    &        /(rho(ipt,0)*CP(ipt,0))
c     ENDDO

      RETURN
      END

*******************************************************************

      REAL FUNCTION SWDK(z,jerlov)
#include "parameter.inc"
c     include 'proc_pars.com'
c     include 'local_pt.com'

      parameter(max=5)
      real Rfac(max),a1(max),a2(max)
c         types =  I       IA      IB      II      III
c             j =  1       2       3       4       5
      data Rfac /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
      data a1   /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
      data a2   / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /

         j = jerlov

c      write(nuout,*) 'time=',ftime,' mon=',mon,' j=',j

      SWDK =        Rfac(j)  * dexp(dble(z/a1(j))) 
     >       + (1.0-Rfac(j)) * dexp(dble(z/a2(j)))

      return
      end
**************************************************************

      SUBROUTINE update_optical(kpp_3d_fields,kpp_const_fields)
      IMPLICIT NONE
#include "kpp_3d_type.com"

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      INTEGER :: ipt,z,max
      parameter(max=5)
      real Rfac(max),a1(max),a2(max)
c         types =  I       IA      IB      II      III
c             j =  1       2       3       4       5
      data Rfac /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
      data a1   /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
      data a2   / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /
                
      ! Recompute double exponential on demand
      IF (kpp_const_fields%L_VARY_OPT) THEN
         DO ipt=1,NPTS
            DO z=0,NZ
               kpp_3d_fields%swdk_opt(ipt,z) = kpp_3d_fields%Rfac(ipt)*
     +            dexp(dble(-kpp_const_fields%dm(z)/
     +            kpp_3d_fields%h1(ipt)))+(1.0-kpp_3d_fields%Rfac(ipt))*
     +            dexp(dble(-kpp_const_fields%dm(z)/
     +            kpp_3d_fields%h2(ipt)))
            ENDDO
         ENDDO
      ELSE
         DO ipt=1,NPTS
            kpp_3d_fields%rfac(ipt)=Rfac(kpp_3d_fields%jerlov(ipt))
            kpp_3d_fields%h1(ipt)=a1(kpp_3d_fields%jerlov(ipt))
            kpp_3d_fields%h2(ipt)=a2(kpp_3d_fields%jerlov(ipt))
         ENDDO
      ENDIF
         
      RETURN
      END
            
      
