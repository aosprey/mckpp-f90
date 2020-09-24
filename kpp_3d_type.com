#include <parameter.inc>
#include <kpp_timer.com>
	TYPE kpp_3D_type

	real :: U(NPTS,NZP1,NVEL),
     +     X(NPTS,NZP1,NSCLR),
     +     Rig(NPTS,NZP1),
     +     dbloc(NPTS,NZ),
     +     Shsq(NPTS,NZP1),
     +     hmixd(NPTS,0:1),
     +     Us(NPTS,NZP1,NVEL,0:1),
     +     Xs(NPTS,NZP1,NSCLR,0:1),
     +     rho(NPTS,0:NZP1tmax),
     +     cp(NPTS,0:NZP1tmax),
     +     buoy(NPTS,NZP1tmax),
     +     rhoh2o(npts),
     +     ocdepth(npts),f(npts),
     +     swfrac(npts,NZP1),
     +     swdk_opt(npts,0:NZ),
     +     difm(npts,0:NZtmax),
     +     difs(npts,0:NZtmax),
     +     dift(npts,0:NZtmax),
     +     wU(npts,0:NZTmax,NVP1),
     +     wX(npts,0:NZTmax,NSP1),
     +     wXNT(npts,0:NZTmax,NSCLR),
     +     ghat(npts,NZTmax),
     +     relax_sst(npts),fcorr(npts),
     +     cplwght(npts_globe),
     +     SST0(npts),fcorr_twod(npts),
     +     sfcorr_twod(npts),
     +     tinc_fcorr(npts,NZP1),
     +     sinc_fcorr(npts,NZP1),
     +     fcorr_withz(npts,NZP1),sfcorr(npts),
     +     sfcorr_withz(npts,NZP1),
     +     advection(npts,maxmodeadv,2),
     +     relax_sal(npts),scorr(npts,NZP1),
     +     relax_ocnT(npts),ocnTcorr(npts,NZP1),
     +     sal_clim(npts,NZP1),ocnT_clim(npts,NZP1),
     +     hmix(npts),kmix(npts),Tref(npts),
     +     uref(npts),vref(npts),Ssurf(npts),
     +     Sref(npts),SSref(npts),
     +     sflux(npts,NSFLXS,5,0:NJDT),
     +     dlat(npts),dlon(npts),
     +     freeze_flag(npts),reset_flag(npts),
     +     dampu_flag(npts),dampv_flag(npts),
     +     U_init(NPTS,NZP1,NVEL)
      logical :: l_ocean(npts),l_initflag(npts)
      integer :: old(NPTS),old_pt,
     +     new(NPTS),new_pt,
     +     jerlov(NPTS),jerlov_pt,
     +     nmodeadv(NPTS,2),
     +     modeadv(NPTS,maxmodeadv,2)
      ENDTYPE kpp_3D_type
      
      TYPE kpp_2D_type
      real :: U(NZP1,NVEL),
     +     X(NZP1,NSCLR),
     +     Rig(NZP1),
     +     dbloc(NZ),
     +     Shsq(NZP1),
     +     hmixd(0:1),
     +     Us(NZP1,NVEL,0:1),
     +     Xs(NZP1,NSCLR,0:1),
     +     rho(0:NZP1tmax),
     +     cp(0:NZP1tmax),
     +     buoy(NZP1tmax),
     +     rhoh2o,ocdepth,f,
     +     swfrac(NZP1),
     +     swdk_opt(0:NZ),
     +     difm(0:NZtmax),
     +     difs(0:NZtmax),
     +     dift(0:NZtmax),
     +     wU(0:NZTmax,NVP1),
     +     wX(0:NZTmax,NSP1),
     +     wXNT(0:NZTmax,NSCLR),
     +     ghat(NZTmax),
     +     relax_sst,fcorr,SST0,
     +     fcorr_twod,
     +     sfcorr_twod,sfcorr,
     +     tinc_fcorr(NZP1),
     +     sinc_fcorr(NZP1),
     +     fcorr_withz(NZP1),
     +     sfcorr_withz(NZP1),
     +     advection(maxmodeadv,2),
     +     relax_sal,scorr(NZP1),relax_ocnT,ocnTcorr(NZP1),
     +     sal_clim(NZP1),ocnT_clim(NZP1),
     +     hmix,kmix,Tref,
     +     uref,vref,Ssurf,
     +     Sref,SSref,
     +     sflux(NSFLXS,5,0:NJDT),
     +     dlat,dlon,
     +     talpha(0:NZP1tmax),sbeta(0:NZP1tmax), ! Not needed outside physics
     +     reset_flag,dampu_flag,dampv_flag
      integer :: old,new,jerlov,
     +     nmodeadv(2),modeadv(maxmodeadv,2)
      logical :: l_ocean,l_initflag,comp_flag
      ENDTYPE kpp_2D_type

      TYPE kpp_const_type
      real :: zm(nzp1),
     +     hm(nzp1),
     +     dm(0:nz),
     +     spd,dpy,twopi,onepi,
     +     grav,vonk,TK0,sbc,epsw,
     +     albocn,sice,EL,SL,FL,FLSN,dto,time,
     +     startt,finalt,dtsec,
     +     iso_thresh
c      real wmt(0:891,0:49)      ! lookup table for wm
c      real wst(0:891,0:49)      ! lookup table for ws
	real, allocatable :: wmt(:,:),wst(:,:),tri(:,:,:)
      integer :: ntime,iso_bot,dt_uvdamp
      logical :: LKPP,LRI,LDD,LICE,LBIO,
     +     LTGRID,LNBFLX,LRHS,L_SSref,
     +     L_RELAX_SST,
     +     L_RELAX_CALCONLY,L_FCORR,
     +     L_FCORR_WITHZ,L_RESTART,
     +     L_SFCORR,L_SFCORR_WITHZ,
     +     L_RELAX_SAL,L_RELAX_OCNT
      ENDTYPE kpp_const_type

      TYPE kpp_timer_type
#ifdef OPENMP
      REAL,dimension(timer_max_timers) :: 
     +     timer_elapsed_time,timer_start_time
#else
      REAL,dimension(timer_max_timers) :: timer_elapsed_time,
     +     timer_start_time
#endif
      LOGICAL,dimension(timer_max_timers) :: timer_running
      CHARACTER(LEN=30),dimension(timer_max_timers) :: timer_all_names
      
      INTEGER :: timer_number_allocated
      ENDTYPE
