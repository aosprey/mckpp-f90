#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
MODULE mckpp_types
  
  USE shr_kind_mod,only: r8 => shr_kind_r8, r4=>shr_kind_r4
  USE ppgrid,      only: pcols,begchunk,endchunk
  USE phys_grid,   only: get_ncols_p
  USE mckpp_parameters, ONLY: nz, nzp1, nztmax, nzp1tmax, nvel, nvp1, nsclr, nsp1, maxmodeadv, njdt, nsflxs
  
  IMPLICIT NONE
  
  PUBLIC

  PUBLIC initialize_mckpp
  PUBLIC kpp_global_type
  PUBLIC kpp_3d_type
  PUBLIC kpp_1d_type
  PUBLIC kpp_const_type

  TYPE kpp_global_type
     real(r8) :: U(PLON,PLAT,NZP1,NVEL),&
          X(PLON,PLAT,NZP1,NVEL),&
          ocdepth(PLON,PLAT),&
          fcorr_twod(PLON,PLAT),&
          sfcorr_twod(PLON,PLAT),&
          sal_clim(PLON,PLAT,NZP1),&
          ocnT_clim(PLON,PLAT,NZP1),&
          bottom_temp(PLON,PLAT),&
          sst(PLON,PLAT),&
          iceconc(PLON,PLAT),&
          usf(PLON,PLAT),&
          vsf(PLON,PLAT),&
          icedepth(PLON,PLAT),&
          snowdepth(PLON,PLAT),dlat,dlon,&
          relax(PLON,PLAT),latitude(PLAT),longitude(PLON)
  END TYPE kpp_global_type

  TYPE kpp_3d_type
     real(r8) :: U(PCOLS,NZP1,NVEL),&
          X(PCOLS,NZP1,NSCLR),&
          Rig(PCOLS,NZP1),&
          dbloc(PCOLS,NZ),&
          Shsq(PCOLS,NZP1),&
          hmixd(PCOLS,0:1),&
          Us(PCOLS,NZP1,NVEL,0:1),&
          Xs(PCOLS,NZP1,NSCLR,0:1),&
          rho(PCOLS,0:NZP1tmax),&
          cp(PCOLS,0:NZP1tmax),&
          buoy(PCOLS,NZP1tmax),&
          rhoh2o(PCOLS),&
          ocdepth(PCOLS),f(PCOLS),&
          swfrac(PCOLS,NZP1),&
          swdk_opt(PCOLS,0:NZ),&
          difm(PCOLS,0:NZtmax),&
          difs(PCOLS,0:NZtmax),&
          dift(PCOLS,0:NZtmax),&
          wU(PCOLS,0:NZTmax,NVP1),&
          wX(PCOLS,0:NZTmax,NSP1),&
          wXNT(PCOLS,0:NZTmax,NSCLR),&
          ghat(PCOLS,NZTmax),&
          relax_sst(PCOLS),fcorr(PCOLS),&
          cplwght(PCOLS),&
          SST0(PCOLS),fcorr_twod(PCOLS),&
          sfcorr_twod(PCOLS),&
          tinc_fcorr(PCOLS,NZP1),&
          sinc_fcorr(PCOLS,NZP1),&
          fcorr_withz(PCOLS,NZP1),sfcorr(PCOLS),&
          sfcorr_withz(PCOLS,NZP1),&
          advection(PCOLS,maxmodeadv,2),&
          relax_sal(PCOLS),scorr(PCOLS,NZP1),&
          relax_ocnT(PCOLS),ocnTcorr(PCOLS,NZP1),&
          sal_clim(PCOLS,NZP1),ocnT_clim(PCOLS,NZP1),&
          hmix(PCOLS),kmix(PCOLS),Tref(PCOLS),&
          uref(PCOLS),vref(PCOLS),Ssurf(PCOLS),&
          Sref(PCOLS),SSref(PCOLS),&
          sflux(PCOLS,NSFLXS,5,0:NJDT),&
          dlat(PCOLS),dlon(PCOLS),&
          freeze_flag(PCOLS),reset_flag(PCOLS),&
          dampu_flag(PCOLS),dampv_flag(PCOLS),&
          U_init(PCOLS,NZP1,NVEL),&
          bottom_temp(PCOLS),&
          taux(PCOLS),tauy(PCOLS),swf(PCOLS),lwf(PCOLS),lhf(PCOLS),shf(PCOLS),&
          rain(PCOLS),snow(PCOLS),landfrac(PCOLS),&
          sst(PCOLS),iceconc(PCOLS),usf(PCOLS),vsf(PCOLS),icedepth(PCOLS),snowdepth(PCOLS), &
	  sflux_cpl(PCOLS,NSFLXS),sst_cpl(PCOLS)
     LOGICAL :: L_OCEAN(PCOLS),L_INITFLAG(PCOLS)
     INTEGER :: old(PCOLS),old_pt,new(PCOLS),new_pt,&
          jerlov(PCOLS),jerlov_pt,nmodeadv(PCOLS,2),&
          modeadv(PCOLS,maxmodeadv,2)          
  END TYPE kpp_3d_type
  
  TYPE kpp_1d_type
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
          reset_flag,dampu_flag,dampv_flag,freeze_flag,cplwght
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
          iso_thresh,alat,alon,delta_lat,delta_lon,DMAX,dscale
     real, allocatable :: wmt(:,:),wst(:,:),tri(:,:,:)
     integer :: nstart,nend,ndtocn,ntime,iso_bot,dt_uvdamp
     logical :: LKPP,LRI,LDD,LICE,LBIO,&
          LTGRID,LNBFLX,LRHS,L_SSref,&
          L_RELAX_SST,&
          L_RELAX_CALCONLY,L_FCORR,&
          L_FCORR_WITHZ,L_RESTART,&
          L_SFCORR,L_SFCORR_WITHZ,&
          L_RELAX_SAL,L_RELAX_OCNT,L_REGGRID,L_LANDSEA,L_JERLOV
     CHARACTER*50 :: forcing_file, sst_file, ice_file, fcorr_file, sfcorr_file,&
          ocnT_file, sal_file, bottom_file, landsea_file, vgrid_file, cplwght_file, &
          advect_file, paras_file, initdata_file
     CHARACTER*17 :: restart_outfile, restart_infile
     INTEGER :: ncid_out,mean_ncid_out,min_ncid_out,max_ncid_out,&
          ndt_per_file,ndt_per_restart,flx_ncid,flx_timein_id,flx_varin_id(7),&
          ndtupdsst,climsst_period,ndtupdice,climice_period,ndtupdcurr,&
          ndtupdfcorr,ndtupdsfcorr,fcorr_period,sfcorr_period,&
          ndtupdocnt,ocnt_period,ndtupdsal,sal_period,ndt_interp_sal,ndt_interp_ocnt,&
          ndtupdbottom,bottom_temp_period,&
          ifirst,ilast,jfirst,jlast,day_out
     LOGICAL :: L_OUTPUT_MEAN,L_OUTPUT_INST,L_OUTPUT_RANGE,&
          L_RESTARTW,L_CLIMSST,L_UPD_CLIMSST,L_PERIODIC_CLIMSST,&
          L_CLIMICE,L_UPD_CLIMICE,L_PERIODIC_CLIMICE,L_CLIM_ICE_DEPTH,L_CLIM_SNOW_ON_ICE,&
          L_BAD_ICE_DEPTH,L_CLIMCURR,L_UPD_CLIMCURR,L_PERIODIC_CLIMCURR,&
          L_UPD_FCORR,L_PERIODIC_FCORR,L_UPD_SFCORR,L_PERIODIC_SFCORR,&
          L_UPD_OCNT,L_PERIODIC_OCNT,L_INTERP_OCNT,L_UPD_SAL,L_PERIODIC_SAL,L_INTERP_SAL,&
          L_VARY_BOTTOM_TEMP,L_UPD_BOTTOM_TEMP,L_PERIODIC_BOTTOM_TEMP,&
          L_OUTKELVIN,L_COUPLE_CURRENTS,L_FLUXDATA,L_REST,L_VGRID_FILE,L_STRETCHGRID,&
          L_CPLWGHT,L_ADVECT,L_INITDATA,L_INTERPINIT,L_NO_ISOTHERM,L_NO_FREEZE,L_DAMP_CURR
#ifdef MCKPP_CAM3
     INTEGER :: relax_sst_in(NY_GLOBE), relax_ocnT_in(NY_GLOBE), relax_sal_in(NY_GLOBE)
#else
     INTEGER :: relax_sst_in(NY), relax_ocnT_in(NY), relax_sal_in(NY)
#endif
     REAL*4 :: dtout,flx_first_timein     
  ENDTYPE kpp_const_type

  TYPE(kpp_global_type) :: kpp_global_fields
  TYPE(kpp_3d_type),allocatable :: kpp_3d_fields(:)
  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

CONTAINS
  
  SUBROUTINE initialize_mckpp
    IMPLICIT NONE

    allocate(kpp_3d_fields(begchunk:endchunk))

    RETURN
  END SUBROUTINE initialize_mckpp
  
END MODULE mckpp_types
#else 
SUBROUTINE dummy 
END SUBROUTINE dummy 
#endif 
