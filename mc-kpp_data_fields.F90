MODULE mckpp_data_fields

  USE mckpp_parameters

  IMPLICIT NONE 

  TYPE kpp_3D_type

     REAL, ALLOCATABLE :: & 
          U(:,:,:), &
          X(:,:,:), &
          Rig(:,:), &
          dbloc(:,:), &
          Shsq(:,:), &
          hmixd(:,:), &
          Us(:,:,:,:), &
          Xs(:,:,:,:), &
          rho(:,:), &
          cp(:,:), &
          buoy(:,:), &
          rhoh2o(:), &
          ocdepth(:), & 
          f(:), &
          swfrac(:,:), &
          swdk_opt(:,:), &
          difm(:,:), &
          difs(:,:), &
          dift(:,:), &
          wU(:,:,:), &
          wX(:,:,:), &
          wXNT(:,:,:), &
          ghat(:,:), &
          relax_sst(:), & 
          fcorr(:), &
          cplwght(:), &
          SST0(:), & 
          fcorr_twod(:), &
          sfcorr_twod(:), &
          tinc_fcorr(:,:), &
          sinc_fcorr(:,:), &
          fcorr_withz(:,:), & 
          sfcorr(:), &
          sfcorr_withz(:,:), &
          advection(:,:,:), &
          relax_sal(:), & 
          scorr(:,:), &
          relax_ocnT(:), & 
          ocnTcorr(:,:), &
          sal_clim(:,:), & 
          ocnT_clim(:,:), &
          hmix(:), & 
          kmix(:), & 
          Tref(:), &
          uref(:), & 
          vref(:), & 
          Ssurf(:), &
          Sref(:), & 
          SSref(:), &
          sflux(:,:,:,:), &
          dlat(:), & 
          dlon(:), &
          freeze_flag(:), & 
          reset_flag(:), &
          dampu_flag(:), & 
          dampv_flag(:), &
          U_init(:,:,:), &
          bottom_temp(:), &
          taux(:), & 
          tauy(:), & 
          swf(:), & 
          lwf(:), & 
          lhf(:), & 
          shf(:), &
          rain(:), & 
          snow(:), &
          sst(:,:), &
          iceconc(:,:), & 
          usf(:,:), & 
          vsf(:,:), &
          icedepth(:,:), & 
          snowdepth(:,:)

     LOGICAL, ALLOCATABLE :: & 
          l_ocean(:), & 
          l_initflag(:)

     INTEGER :: & 
          old_pt, & 
          new_pt, &
          jerlov_pt

     INTEGER, ALLOCATABLE :: & 
          old(:), &
          new(:), & 
          jerlov(:), & 
          nmodeadv(:,:), &
          modeadv(:,:,:)     

  ENDTYPE kpp_3D_type


  TYPE kpp_1D_type

     REAL :: & 
          rhoh2o, &
          ocdepth, &
          f, &
          cplwght, &
          relax_sst, & 
          fcorr, & 
          SST0, &
          fcorr_twod, &
          sfcorr_twod, & 
          sfcorr, &
          dlat, & 
          dlon, &
          relax_sal, & 
          relax_ocnT, & 
          hmix, & 
          kmix, & 
          Tref, & 
          uref, & 
          vref, & 
          Ssurf, &
          Sref, & 
          SSref,& 
          reset_flag, & 
          dampu_flag, & 
          dampv_flag, & 
          freeze_flag

     REAL, ALLOCATABLE :: & 
          U(:,:), &
          U_init(:,:), &
          X(:,:), &
          Rig(:), &
          dbloc(:), &
          Shsq(:), &
          hmixd(:), &
          Us(:,:,:), &
          Xs(:,:,:), &
          rho(:), &
          cp(:), &
          buoy(:), &
          swfrac(:), &
          swdk_opt(:), &
          difm(:), &
          difs(:), &
          dift(:), &
          wU(:,:), &
          wX(:,:), &
          wXNT(:,:), &
          ghat(:), &
          tinc_fcorr(:), &
          sinc_fcorr(:), &
          fcorr_withz(:), &
          sfcorr_withz(:), &
          advection(:,:), &
          scorr(:), & 
          ocnTcorr(:),&
          sal_clim(:), & 
          ocnT_clim(:), &
          sflux(:,:,:), &
          talpha(:), & 
          sbeta(:)

     INTEGER :: & 
          old, & 
          new, & 
          jerlov, &
          nmodeadv(2), & 
          point

     INTEGER, ALLOCATABLE :: & 
          modeadv(:,:)

     LOGICAL :: & 
          l_ocean, & 
          l_initflag, & 
          comp_flag

  ENDTYPE kpp_1D_type


  TYPE kpp_const_type

     REAL*4 :: & 
          dtout, & 
          flx_first_timein

     REAL :: & 
          spd, & 
          dpy, & 
          twopi, & 
          onepi, &
          grav, & 
          vonk, & 
          TK0, & 
          sbc, & 
          epsw, &
          albocn, & 
          sice, & 
          EL, & 
          SL, & 
          FL,  & 
          FLSN, & 
          dto, & 
          time, &
          startt, & 
          finalt, & 
          dtsec, &
          iso_thresh, & 
          dmax, & 
          alat, &
          alon, & 
          delta_lat, & 
          delta_lon, & 
          dscale

     REAL, ALLOCATABLE :: & 
          zm(:), &
          hm(:), &
          dm(:), &
          relax_sst_in(:), & 
          relax_sal_in(:), & 
          relax_ocnt_in(:), & 
          wmt(:,:), & 
          wst(:,:), & 
          tri(:,:,:)

     INTEGER :: & 
          nstart, & 
          nend, & 
          ndtocn, & 
          ntime, & 
          iso_bot, & 
          dt_uvdamp, & 
          ndt_per_restart, & 
          flx_ncid, & 
          flx_timein_id, & 
          flx_varin_id(7), &
          ndtupdsst, & 
          climsst_period, & 
          ndtupdice, & 
          climice_period, & 
          ndtupdcurr, &
          ndtupdfcorr, & 
          ndtupdsfcorr, & 
          fcorr_period, & 
          sfcorr_period, &
          ndtupdocnt, & 
          ocnt_period, & 
          ndtupdsal, & 
          sal_period, & 
          ndt_interp_sal, & 
          ndt_interp_ocnt, &
          ndtupdbottom, & 
          bottom_temp_period, &
          ifirst, & 
          ilast, & 
          jfirst, & 
          jlast, & 
          day_out

     LOGICAL :: & 
          LKPP, & 
          LRI, & 
          LDD, & 
          LICE, & 
          LBIO, &
          LTGRID, & 
          LNBFLX, & 
          LRHS, & 
          L_SSref, &
          L_RELAX_SST, &
          L_RELAX_CALCONLY, & 
          L_FCORR, &
          L_FCORR_WITHZ, & 
          L_RESTART, &
          L_SFCORR, & 
          L_SFCORR_WITHZ, &
          L_RELAX_SAL, & 
          L_RELAX_OCNT, & 
          L_RESTARTW, & 
          L_CLIMSST, & 
          L_UPD_CLIMSST, & 
          L_PERIODIC_CLIMSST, &
          L_CLIMICE, & 
          L_UPD_CLIMICE, & 
          L_PERIODIC_CLIMICE, & 
          L_CLIM_ICE_DEPTH, & 
          L_CLIM_SNOW_ON_ICE, &
          L_BAD_ICE_DEPTH, & 
          L_CLIMCURR, & 
          L_UPD_CLIMCURR, & 
          L_PERIODIC_CLIMCURR, &
          L_UPD_FCORR, & 
          L_PERIODIC_FCORR, & 
          L_UPD_SFCORR, & 
          L_PERIODIC_SFCORR, &
          L_UPD_OCNT, & 
          L_PERIODIC_OCNT, & 
          L_INTERP_OCNT, & 
          L_UPD_SAL, & 
          L_PERIODIC_SAL, & 
          L_INTERP_SAL, &
          L_VARY_BOTTOM_TEMP, & 
          L_UPD_BOTTOM_TEMP, & 
          L_PERIODIC_BOTTOM_TEMP, &
          L_OUTKELVIN, & 
          L_COUPLE_CURRENTS, & 
          L_FLUXDATA, & 
          L_REST, & 
          L_ADVECT, &
          L_REGGRID, & 
          L_LANDSEA, & 
          L_VGRID_FILE, & 
          L_STRETCHGRID, & 
          L_CPLWGHT, & 
          L_JERLOV, & 
          L_INITDATA, & 
          L_INTERPINIT, & 
          L_NO_ISOTHERM, & 
          L_NO_FREEZE, & 
          L_DAMP_CURR

     CHARACTER(LEN=max_nc_filename_len) :: & 
          forcing_file, & 
          sst_file, & 
          ice_file, & 
          fcorr_file, & 
          sfcorr_file, &
          ocnT_file, & 
          sal_file, & 
          bottom_file, & 
          advect_file, & 
          landsea_file, & 
          vgrid_file, &
          cplwght_file, & 
          paras_file, & 
          initdata_file

     CHARACTER(LEN=max_restart_filename_len) :: & 
          restart_outfile, & 
          restart_infile
  
  ENDTYPE kpp_const_type

END MODULE mckpp_data_fields
