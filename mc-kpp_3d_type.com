#include <parameter.inc>
#include <kpp_timer.com>
 
  INTEGER,parameter :: N_VAROUTS=23,N_SINGOUTS=12,N_ZPROFS_MAX=5
  INTEGER,parameter :: NVEC_MEAN=N_VAROUTS,NVEC_RANGE=N_VAROUTS,&
       NSCLR_MEAN=N_SINGOUTS,NSCLR_RANGE=N_SINGOUTS

  TYPE kpp_3D_type
     real :: U(NPTS,NZP1,NVEL),&
          X(NPTS,NZP1,NSCLR),&
          Rig(NPTS,NZP1),&
          dbloc(NPTS,NZ),&
          Shsq(NPTS,NZP1),&
          hmixd(NPTS,0:1),&
          Us(NPTS,NZP1,NVEL,0:1),&
          Xs(NPTS,NZP1,NSCLR,0:1),&
          rho(NPTS,0:NZP1tmax),&
          cp(NPTS,0:NZP1tmax),&
          buoy(NPTS,NZP1tmax),&
          rhoh2o(npts),&
          ocdepth(npts),f(npts),&
          swfrac(npts,NZP1),&
          swdk_opt(npts,0:NZ),&
          difm(npts,0:NZtmax),&
          difs(npts,0:NZtmax),&
          dift(npts,0:NZtmax),&
          wU(npts,0:NZTmax,NVP1),&
          wX(npts,0:NZTmax,NSP1),&
          wXNT(npts,0:NZTmax,NSCLR),&
          ghat(npts,NZTmax),&
          relax_sst(npts),fcorr(npts),&
          cplwght(npts_globe),&
          SST0(npts),fcorr_twod(npts),&
          sfcorr_twod(npts),&
          tinc_fcorr(npts,NZP1),&
          sinc_fcorr(npts,NZP1),&
          fcorr_withz(npts,NZP1),sfcorr(npts),&
          sfcorr_withz(npts,NZP1),&
          advection(npts,maxmodeadv,2),&
          relax_sal(npts),scorr(npts,NZP1),&
          relax_ocnT(npts),ocnTcorr(npts,NZP1),&
          sal_clim(npts,NZP1),ocnT_clim(npts,NZP1),&
          hmix(npts),kmix(npts),Tref(npts),&
          uref(npts),vref(npts),Ssurf(npts),&
          Sref(npts),SSref(npts),&
          sflux(npts,NSFLXS,5,0:NJDT),&
          dlat(npts),dlon(npts),&
          freeze_flag(npts),reset_flag(npts),&
          dampu_flag(npts),dampv_flag(npts),&
          U_init(NPTS,NZP1,NVEL),&
          bottom_temp(NPTS),&
          taux(NPTS),tauy(NPTS),swf(NPTS),lwf(NPTS),lhf(NPTS),shf(NPTS),&
          rain(NPTS),snow(NPTS),&
#ifdef MCKPP_COUPLE
          sst(NX_GLOBE,NY_GLOBE),iceconc(NX_GLOBE,NY_GLOBE),usf(NX_GLOBE,NY_GLOBE),vsf(NX_GLOBE,NY_GLOBE),&
          icedepth(NX_GLOBE,NY_GLOBE),snowdepth(NX_GLOBE,NY_GLOBE)
#else
          sst(NX,NY),iceconc(NX,NY),usf(NX,NY),vsf(NX,NY),icedepth(NX,NY),snowdepth(NX,NY)
#endif
     logical :: l_ocean(npts),l_initflag(npts)
     integer :: old(NPTS),old_pt,new(NPTS),new_pt,&
          jerlov(NPTS),jerlov_pt,nmodeadv(NPTS,2),&
          modeadv(NPTS,maxmodeadv,2)     
     real,allocatable :: VEC_mean(:,:,:),VEC_range(:,:,:,:),&
          SCLR_mean(:,:),SCLR_range(:,:,:)
  ENDTYPE kpp_3D_type
  
  TYPE kpp_1D_type
     real :: U(NZP1,NVEL),U_init(NZP1,NVEL),&
          X(NZP1,NSCLR),&
          Rig(NZP1),&
          dbloc(NZ),&
          Shsq(NZP1),&
          hmixd(0:1),&
          Us(NZP1,NVEL,0:1),&
          Xs(NZP1,NSCLR,0:1),&
          rho(0:NZP1tmax),&
          cp(0:NZP1tmax),&
          buoy(NZP1tmax),&
          rhoh2o,ocdepth,f,&
          swfrac(NZP1),&
          swdk_opt(0:NZ),&
          difm(0:NZtmax),&
          difs(0:NZtmax),&
          dift(0:NZtmax),&
          wU(0:NZTmax,NVP1),&
          wX(0:NZTmax,NSP1),&
          wXNT(0:NZTmax,NSCLR),&
          ghat(NZTmax),&
          cplwght,&
          relax_sst,fcorr,SST0,&
          fcorr_twod,&
          sfcorr_twod,sfcorr,&
          tinc_fcorr(NZP1),&
          sinc_fcorr(NZP1),&
          fcorr_withz(NZP1),&
          sfcorr_withz(NZP1),&
          advection(maxmodeadv,2),&
          relax_sal,scorr(NZP1),relax_ocnT,ocnTcorr(NZP1),&
          sal_clim(NZP1),ocnT_clim(NZP1),&
          hmix,kmix,Tref,uref,vref,Ssurf,&
          Sref,SSref,sflux(NSFLXS,5,0:NJDT),&
          dlat,dlon,&
          talpha(0:NZP1tmax),sbeta(0:NZP1tmax),& ! Not needed outside physics
          reset_flag,dampu_flag,dampv_flag,freeze_flag
     integer :: old,new,jerlov,&
          nmodeadv(2),modeadv(maxmodeadv,2),point
     logical :: l_ocean,l_initflag,comp_flag
  ENDTYPE kpp_1D_type
  
  TYPE kpp_const_type
     real :: zm(nzp1),&
          hm(nzp1),&
          dm(0:nz),&
          spd,dpy,twopi,onepi,&
          grav,vonk,TK0,sbc,epsw,&
          albocn,sice,EL,SL,FL,FLSN,dto,time,&
          startt,finalt,dtsec,&
          iso_thresh,dmax,& 
          relax_sst_in(ny),relax_sal_in(ny),relax_ocnt_in(ny),& 
          alat,alon,delta_lat,delta_lon,dscale
     real, allocatable :: wmt(:,:),wst(:,:),tri(:,:,:)
     integer :: nstart,nend,ndtocn,ntime,iso_bot,dt_uvdamp
     logical :: LKPP,LRI,LDD,LICE,LBIO,&
          LTGRID,LNBFLX,LRHS,L_SSref,&
          L_RELAX_SST,&
          L_RELAX_CALCONLY,L_FCORR,&
          L_FCORR_WITHZ,L_RESTART,&
          L_SFCORR,L_SFCORR_WITHZ,&
          L_RELAX_SAL,L_RELAX_OCNT
     CHARACTER*50 :: forcing_file, sst_file, ice_file, fcorr_file, sfcorr_file,&
          ocnT_file, sal_file, bottom_file, advect_file, landsea_file, vgrid_file, &
          cplwght_file, paras_file, initdata_file
     CHARACTER*17 :: restart_outfile, restart_infile
     INTEGER :: ncid_out,mean_ncid_out,min_ncid_out,max_ncid_out,&
          ndt_per_file,ndt_per_restart,flx_ncid,flx_timein_id,flx_varin_id(7),&
          ndtupdsst,climsst_period,ndtupdice,climice_period,ndtupdcurr,&
          ndtupdfcorr,ndtupdsfcorr,fcorr_period,sfcorr_period,&
          ndtupdocnt,ocnt_period,ndtupdsal,sal_period,ndt_interp_sal,ndt_interp_ocnt,&
          ndtupdbottom,bottom_temp_period,&
          ifirst,ilast,jfirst,jlast,day_out
     INTEGER,dimension(N_VAROUTS) :: &
          ndt_varout_inst,ndt_varout_range,ndt_varout_mean,&
          zprof_varout_inst,zprof_varout_mean,zprof_varout_range,&
          varid_vec_inst,varid_vec_mean,varid_vec_range,&
          ntout_vec_inst,ntout_vec_mean,ntout_vec_range
     INTEGER,dimension(N_SINGOUTS) :: &
          ndt_singout_inst,ndt_singout_range,ndt_singout_mean,&
          varid_sing_inst,varid_sing_mean,varid_sing_range,&
          ntout_sing_inst,ntout_sing_mean,ntout_sing_range
     LOGICAL :: L_OUTPUT_MEAN,L_OUTPUT_INST,L_OUTPUT_RANGE,&
          L_RESTARTW,L_CLIMSST,L_UPD_CLIMSST,L_PERIODIC_CLIMSST,&
          L_CLIMICE,L_UPD_CLIMICE,L_PERIODIC_CLIMICE,L_CLIM_ICE_DEPTH,L_CLIM_SNOW_ON_ICE,&
          L_BAD_ICE_DEPTH,L_CLIMCURR,L_UPD_CLIMCURR,L_PERIODIC_CLIMCURR,&
          L_UPD_FCORR,L_PERIODIC_FCORR,L_UPD_SFCORR,L_PERIODIC_SFCORR,&
          L_UPD_OCNT,L_PERIODIC_OCNT,L_INTERP_OCNT,L_UPD_SAL,L_PERIODIC_SAL,L_INTERP_SAL,&
          L_VARY_BOTTOM_TEMP,L_UPD_BOTTOM_TEMP,L_PERIODIC_BOTTOM_TEMP,&
          L_OUTKELVIN,L_COUPLE_CURRENTS,L_FLUXDATA,L_REST,L_ADVECT,&
	  L_REGGRID,L_LANDSEA,L_VGRID_FILE,L_STRETCHGRID,L_CPLWGHT,L_JERLOV,& 
          L_INITDATA,L_INTERPINIT,L_NO_ISOTHERM,L_NO_FREEZE,L_DAMP_CURR
     LOGICAL,dimension(NZP1,0:N_ZPROFS_MAX) :: zprofs_mask
     INTEGER,dimension(NZP1,0:N_ZPROFS_MAX) :: zprofs
     INTEGER,dimension(0:N_ZPROFS_MAX) :: zprofs_nvalid   
     REAL*4 :: dtout,flx_first_timein
     
  ENDTYPE kpp_const_type
  
  TYPE kpp_timer_type
#ifdef OPENMP
     REAL,dimension(timer_max_timers) ::& 
          timer_elapsed_time,timer_start_time
#else
     REAL,dimension(timer_max_timers) :: timer_elapsed_time,&
          timer_start_time
#endif
     LOGICAL,dimension(timer_max_timers) :: timer_running
     CHARACTER(LEN=30),dimension(timer_max_timers) :: timer_all_names
     
     INTEGER :: timer_number_allocated
  ENDTYPE kpp_timer_type
  
