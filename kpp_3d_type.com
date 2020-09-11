#include "parameter.inc"
#include "kpp_timer.com"
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
     +     freeze_flag(npts),
     +     dampu_flag(npts),dampv_flag(npts),
     +     U_init(NPTS,NZP1,NVEL),fcorr_nsol_coeff(NPTS),
     +     fcorr_nsol(NPTS),reset_flag(npts),hekman(NPTS),
     +	   tinc_ekadv(NPTS,NZP1),sinc_ekadv(NPTS,NZP1),
     +     runoff_incr(NPTS),u_clim(NPTS,NZP1),
     +     v_clim(NPTS,NZP1),relax_curr(npts),
     +	   clim_sst(NX_GLOBE,NY_GLOBE),anom_sst(NX_GLOBE,NY_GLOBE),
     +     clim_ice(NX_GLOBE,NY_GLOBE),anom_ice(NX_GLOBE,NY_GLOBE),
     +	   Rfac(NPTS),h1(NPTS),h2(NPTS),sst_lag(NPTS),
     +     sst_lag_tmp(NPTS)
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
     +     hmix,kmix,Tref,hekman,
     +     uref,vref,Ssurf,
     +     Sref,SSref,
     +     sflux(NSFLXS,5,0:NJDT),
     +     dlat,dlon,
     +     talpha(0:NZP1tmax),sbeta(0:NZP1tmax), ! Not needed outside physics
     +     dampu_flag,dampv_flag,reset_flag,ekvel(NZP1),ekadv(NZP1,2),
     +	   u_clim(NZP1),v_clim(NZP1),relax_curr,rfac,h1,h2,
     +     dm(0:nz),hm(nzp1),zm(nzp1)
      real :: tri(0:NZtmax,0:1,NGRID)
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
     +     iso_thresh,slab_depth,barrier_dT,barrier_subdepth
c      real wmt(0:891,0:49)      ! lookup table for wm
c      real wst(0:891,0:49)      ! lookup table for ws
      real, allocatable :: wmt(:,:),wst(:,:),tri(:,:,:)
      integer :: ntime,iso_bot,dt_uvdamp,ekmax,ekadv_max,sst_lag_len,
     +     sst_smooth_ifirst,sst_smooth_ilast,sst_smooth_jfirst,
     +     sst_smooth_jlast,sst_smooth_blend,barrier_ifirst,
     +     barrier_ilast,barrier_jfirst,barrier_jlast
      logical :: LKPP,LRI,LDD,LICE,LBIO,
     +     LTGRID,LNBFLX,LRHS,L_SSref,
     +     L_RELAX_SST,
     +     L_RELAX_CALCONLY,L_FCORR,
     +     L_FCORR_WITHZ,L_RESTART,
     +     L_SFCORR,L_SFCORR_WITHZ,
     +     L_RELAX_SAL,L_RELAX_OCNT,L_DIST_RUNOFF,
     +     L_DAMP_CURR,L_SLAB,L_COLUMBIA_LAND,L_FCORR_NSOL,
     +	   L_EKMAN_PUMP,L_RELAX_CURR,L_VARY_OPT,L_SST_LAG_FUDGE,
     +     L_SST_LAG,L_SST_SMOOTH,L_SST_SMOOTH_X,L_SST_SMOOTH_Y,
     +     L_SST_SMOOTH_ANOM,L_SST_ANOM_FUDGE,L_BARRIER_REMOVE,
     +     L_BARRIER_SALISO,L_BARRIER_SALVAVG,L_NO_EGTP
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
