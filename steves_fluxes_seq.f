      SUBROUTINE fluxes

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
      
#ifdef COUPLE
#ifdef CFS
      include 'parameter.cfs_coupled.inc'
#else
      include 'parameter.oasis2.inc'
#endif
#else
#ifdef CFS
      include 'parameter.cfs_forced.inc'
#else
      include 'parameter.forced.inc'
#endif
#endif
      include 'flx_sfc.com'
      include 'ocn_paras.com'
      include 'flx_paras.com'
      include 'flx_in.com'
      include 'local_pt.com'
      include 'times.com'
      include 'couple.com'

      REAL taux(NPTS),tauy(NPTS),
     $     swf(NPTS),lwf(NPTS),lhf(NPTS),shf(NPTS),
     $     rain(NPTS),snow(NPTS)

      IF (.NOT. L_COUPLE) THEN

         IF (.NOT. L_FLUXDATA) THEN
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
         ELSE
            call read_fluxes(taux,tauy,swf,lwf,lhf,shf,rain,snow)
         ENDIF
            
         DO ipt=1,npts
            IF ((taux(ipt) .EQ. 0.0) .AND. (tauy(ipt) .EQ. 0.0)) THEN
               taux(ipt)=1.e-10
            ENDIF
            sflux(ipt,1,5,0)=taux(ipt)
            sflux(ipt,2,5,0)=tauy(ipt)
            
            sflux(ipt,3,5,0)=swf(ipt)
            
            sflux(ipt,4,5,0)=lwf(ipt)+lhf(ipt)+shf(ipt)-snow(ipt)*FLSN
            
            sflux(ipt,5,5,0)=0.0 ! Melting of sea-ice = 0.0
            
            sflux(ipt,6,5,0)=(rain(ipt)+snow(ipt)+(lhf(ipt)/EL))
            
            call ntflx
            
         ENDDO
      ELSE
         call coupled_flux(swf,lwf,rain,ntime+ndtocn-1)
         call coupled_stress(taux,tauy,ntime+ndtocn-1)
         DO ipt=1,npts
            IF ((taux(ipt) .EQ. 0.0) .AND. (tauy(ipt) .EQ. 0.0)) THEN
               taux(ipt)=1.e-10
            ENDIF
            sflux(ipt,1,5,0)=taux(ipt)
            sflux(ipt,2,5,0)=tauy(ipt)
            
            sflux(ipt,3,5,0)=swf(ipt)
            
            sflux(ipt,4,5,0)=lwf(ipt) 
            sflux(ipt,5,5,0)=0.0 ! Melting of sea-ice = 0.0
            
            sflux(ipt,6,5,0)=rain(ipt)  ! assuming rain = P-E          
            call ntflx
         ENDDO
      ENDIF
         
      

      RETURN
      END

************************************************************************
      SUBROUTINE ntflx

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#ifdef COUPLE
#ifdef CFS
      include 'parameter.cfs_coupled.inc'
#else
      include 'parameter.oasis2.inc'
#endif
#else
#ifdef CFS
      include 'parameter.cfs_forced.inc'
#else
      include 'parameter.forced.inc'
#endif
#endif
      include 'times.com'
      include 'flx_sfc.com'
      include 'flx_profs.com'
      include 'vert_pgrid.com'
      include 'ocn_paras.com'
      include 'local_pt.com'

      INTEGER k
      REAL SWDK
      REAL SWDK_OPT(NPTS,0:NZ)
      EXTERNAL SWDK

      COMMON /SWDK_SAVE/SWDK_OPT

      IF (ntime .eq. 1) THEN
         DO k=0,NZ
            swdk_opt(ipt,k)=swdk(-dm(k))
         ENDDO
      ENDIF
      DO k=0,NZ
         wXNT(ipt,k,1)=-sflux(ipt,3,5,0)*swdk_opt(ipt,k)
     &        /(rho(ipt,0)*CP(ipt,0))
      ENDDO

c     DO k=0,NZ
c        wXNT(ipt,k,1)=-sflux(ipt,3,5,0)*swdk(-dm(k))
c    &        /(rho(ipt,0)*CP(ipt,0))
c     ENDDO

      RETURN
      END

*******************************************************************

      REAL FUNCTION SWDK(z)
#ifdef COUPLE
#ifdef CFS
      include 'parameter.cfs_coupled.inc'
#else
      include 'parameter.oasis2.inc'
#endif
#else
#ifdef CFS
      include 'parameter.cfs_forced.inc'
#else
      include 'parameter.forced.inc'
#endif
#endif
      include 'proc_pars.com'
      include 'local_pt.com'

      parameter(max=5)
      real Rfac(max),a1(max),a2(max)
c         types =  I       IA      IB      II      III
c             j =  1       2       3       4       5
      data Rfac /  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /
      data a1   /  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /
      data a2   / 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /

         j = jerlov(ipt)

c      write(nuout,*) 'time=',ftime,' mon=',mon,' j=',j

      SWDK =        Rfac(j)  * dexp(dble(z/a1(j))) 
     >       + (1.0-Rfac(j)) * dexp(dble(z/a2(j)))

      return
      end
**************************************************************
