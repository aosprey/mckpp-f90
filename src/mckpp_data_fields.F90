MODULE mckpp_data_fields

  USE mckpp_mpi_control, ONLY : npts_local
  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, max_restart_filename_len
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
        dlat_all(:), & 
        dlon(:), &
        dlon_all(:), &
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
        l_initflag(:), &
        run_physics(:)

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
        num_timesteps, & 
        iso_bot, & 
        dt_uvdamp, & 
        ndt_per_restart, & 
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
        L_DAMP_CURR, &
        L_COUPLE

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

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

CONTAINS

  SUBROUTINE mckpp_allocate_3d_fields()

    ! Most fields are on local subdomain (npts_local)
    ALLOCATE( kpp_3d_fields%u(npts_local,nzp1,nvel) ) 
    ALLOCATE( kpp_3d_fields%x(npts_local,nzp1,nsclr) ) 
    ALLOCATE( kpp_3d_fields%rig(npts_local,nzp1) )
    ALLOCATE( kpp_3d_fields%dbloc(npts_local,nz) )
    ALLOCATE( kpp_3d_fields%shsq(npts_local,nzp1) )
    ALLOCATE( kpp_3d_fields%hmixd(npts_local,0:1) )
    ALLOCATE( kpp_3d_fields%us(npts_local,nzp1,nvel,0:1) )
    ALLOCATE( kpp_3d_fields%xs(npts_local,nzp1,nsclr,0:1) )
    ALLOCATE( kpp_3d_fields%rho(npts_local,0:nzp1tmax) )
    ALLOCATE( kpp_3d_fields%cp(npts_local,0:nzp1tmax) )
    ALLOCATE( kpp_3d_fields%buoy(npts_local,nzp1tmax) )
    ALLOCATE( kpp_3d_fields%rhoh2o(npts_local) )
    ALLOCATE( kpp_3d_fields%ocdepth(npts_local) )
    ALLOCATE( kpp_3d_fields%f(npts_local) )
    ALLOCATE( kpp_3d_fields%swfrac(npts_local,nzp1) )
    ALLOCATE( kpp_3d_fields%swdk_opt(npts_local,0:nz) )
    ALLOCATE( kpp_3d_fields%difm(npts_local,0:nztmax) )
    ALLOCATE( kpp_3d_fields%difs(npts_local,0:nztmax) )
    ALLOCATE( kpp_3d_fields%dift(npts_local,0:nztmax) )
    ALLOCATE( kpp_3d_fields%wu(npts_local,0:nztmax,nvp1) )
    ALLOCATE( kpp_3d_fields%wx(npts_local,0:nztmax,nsp1) )
    ALLOCATE( kpp_3d_fields%wxnt(npts_local,0:nztmax,nsclr) )
    ALLOCATE( kpp_3d_fields%ghat(npts_local,nztmax) )
    ALLOCATE( kpp_3d_fields%relax_sst(npts_local) ) 
    ALLOCATE( kpp_3d_fields%fcorr(npts_local) )
    ALLOCATE( kpp_3d_fields%sst0(npts_local) )
    ALLOCATE( kpp_3d_fields%fcorr_twod(npts_local) )
    ALLOCATE( kpp_3d_fields%sfcorr_twod(npts_local) )
    ALLOCATE( kpp_3d_fields%tinc_fcorr(npts_local,nzp1) )
    ALLOCATE( kpp_3d_fields%sinc_fcorr(npts_local,nzp1) )
    ALLOCATE( kpp_3d_fields%fcorr_withz(npts_local,nzp1) ) 
    ALLOCATE( kpp_3d_fields%sfcorr(npts_local) )
    ALLOCATE( kpp_3d_fields%sfcorr_withz(npts_local,nzp1) )
    ALLOCATE( kpp_3d_fields%advection(npts_local,maxmodeadv,2) )
    ALLOCATE( kpp_3d_fields%relax_sal(npts_local) ) 
    ALLOCATE( kpp_3d_fields%scorr(npts_local,nzp1) )
    ALLOCATE( kpp_3d_fields%relax_ocnt(npts_local) ) 
    ALLOCATE( kpp_3d_fields%ocntcorr(npts_local,nzp1) )
    ALLOCATE( kpp_3d_fields%sal_clim(npts_local,nzp1) ) 
    ALLOCATE( kpp_3d_fields%ocnt_clim(npts_local,nzp1) )
    ALLOCATE( kpp_3d_fields%hmix(npts_local) )
    ALLOCATE( kpp_3d_fields%kmix(npts_local) )
    ALLOCATE( kpp_3d_fields%tref(npts_local) )
    ALLOCATE( kpp_3d_fields%uref(npts_local) )
    ALLOCATE( kpp_3d_fields%vref(npts_local) )
    ALLOCATE( kpp_3d_fields%ssurf(npts_local) )
    ALLOCATE( kpp_3d_fields%sref(npts_local) )
    ALLOCATE( kpp_3d_fields%ssref(npts_local) )
    ALLOCATE( kpp_3d_fields%sflux(npts_local,nsflxs,5,0:njdt) )
    ALLOCATE( kpp_3d_fields%dlat(npts_local) )
    ALLOCATE( kpp_3d_fields%dlon(npts_local) )
    ALLOCATE( kpp_3d_fields%freeze_flag(npts_local) )
    ALLOCATE( kpp_3d_fields%reset_flag(npts_local) )
    ALLOCATE( kpp_3d_fields%dampu_flag(npts_local) )
    ALLOCATE( kpp_3d_fields%dampv_flag(npts_local) )
    ALLOCATE( kpp_3d_fields%u_init(npts_local,nzp1,nvel) )
    ALLOCATE( kpp_3d_fields%bottom_temp(npts_local) )
    ALLOCATE( kpp_3d_fields%taux(npts_local) )
    ALLOCATE( kpp_3d_fields%tauy(npts_local) )
    ALLOCATE( kpp_3d_fields%swf(npts_local) )
    ALLOCATE( kpp_3d_fields%lwf(npts_local) )
    ALLOCATE( kpp_3d_fields%lhf(npts_local) ) 
    ALLOCATE( kpp_3d_fields%shf(npts_local) )
    ALLOCATE( kpp_3d_fields%rain(npts_local) )
    ALLOCATE( kpp_3d_fields%snow(npts_local) )
    ALLOCATE( kpp_3d_fields%l_ocean(npts_local) ) 
    ALLOCATE( kpp_3d_fields%l_initflag(npts_local) )
    ALLOCATE( kpp_3d_fields%run_physics(npts_local) ) 
    ALLOCATE( kpp_3d_fields%old(npts_local) ) 
    ALLOCATE( kpp_3d_fields%new(npts_local) )
    ALLOCATE( kpp_3d_fields%jerlov(npts_local) )
    ALLOCATE( kpp_3d_fields%nmodeadv(npts_local,2) )
    ALLOCATE( kpp_3d_fields%modeadv(npts_local,maxmodeadv,2) )
    ALLOCATE( kpp_3d_fields%cplwght(npts_local) )

    ! These fields are on the full domain (npts) or
    ! driving domain (npts_globe)
    ! May only be valid on root pe
    ! Do we need to store these globally? 
    IF ( kpp_const_fields%L_COUPLE ) THEN 
      ALLOCATE( kpp_3d_fields%sst(nx_globe,ny_globe) )
      ALLOCATE( kpp_3d_fields%iceconc(nx_globe,ny_globe) )
      ALLOCATE( kpp_3d_fields%usf(nx_globe,ny_globe) )
      ALLOCATE( kpp_3d_fields%vsf(nx_globe,ny_globe) )
      ALLOCATE( kpp_3d_fields%icedepth(nx_globe,ny_globe) )
      ALLOCATE( kpp_3d_fields%snowdepth(nx_globe,ny_globe) )
    ELSE
      ALLOCATE( kpp_3d_fields%sst(nx,ny) )
      ALLOCATE( kpp_3d_fields%iceconc(nx,ny) )
      ALLOCATE( kpp_3d_fields%usf(nx,ny) )
      ALLOCATE( kpp_3d_fields%vsf(nx,ny) )
      ALLOCATE( kpp_3d_fields%icedepth(nx,ny) )
      ALLOCATE( kpp_3d_fields%snowdepth(nx,ny) )
    END IF 
    ALLOCATE( kpp_3d_fields%dlat_all(npts) )
    ALLOCATE( kpp_3d_fields%dlon_all(npts) )

  END SUBROUTINE mckpp_allocate_3d_fields


  SUBROUTINE mckpp_allocate_1d_fields(kpp_1d_fields) 

    TYPE(kpp_1d_type) :: kpp_1d_fields 

    ALLOCATE( kpp_1d_fields%u(nzp1,nvel) )
    ALLOCATE( kpp_1d_fields%u_init(nzp1,nvel) )
    ALLOCATE( kpp_1d_fields%x(nzp1,nsclr) )
    ALLOCATE( kpp_1d_fields%rig(nzp1) )
    ALLOCATE( kpp_1d_fields%dbloc(nz) )
    ALLOCATE( kpp_1d_fields%shsq(nzp1) )
    ALLOCATE( kpp_1d_fields%hmixd(0:1) )
    ALLOCATE( kpp_1d_fields%us(nzp1,nvel,0:1) )
    ALLOCATE( kpp_1d_fields%xs(nzp1,nsclr,0:1) )
    ALLOCATE( kpp_1d_fields%rho(0:nzp1tmax) )
    ALLOCATE( kpp_1d_fields%cp(0:nzp1tmax) )
    ALLOCATE( kpp_1d_fields%buoy(nzp1tmax) )
    ALLOCATE( kpp_1d_fields%swfrac(nzp1) )
    ALLOCATE( kpp_1d_fields%swdk_opt(0:nz) )
    ALLOCATE( kpp_1d_fields%difm(0:nztmax) )
    ALLOCATE( kpp_1d_fields%difs(0:nztmax) )
    ALLOCATE( kpp_1d_fields%dift(0:nztmax) )
    ALLOCATE( kpp_1d_fields%wu(0:nztmax,nvp1) )
    ALLOCATE( kpp_1d_fields%wx(0:nztmax,nsp1) )
    ALLOCATE( kpp_1d_fields%wxnt(0:nztmax,nsclr) )
    ALLOCATE( kpp_1d_fields%ghat(nztmax) )
    ALLOCATE( kpp_1d_fields%tinc_fcorr(nzp1) )
    ALLOCATE( kpp_1d_fields%sinc_fcorr(nzp1) )
    ALLOCATE( kpp_1d_fields%fcorr_withz(nzp1) )
    ALLOCATE( kpp_1d_fields%sfcorr_withz(nzp1) )
    ALLOCATE( kpp_1d_fields%advection(maxmodeadv,2) )
    ALLOCATE( kpp_1d_fields%scorr(nzp1) ) 
    ALLOCATE( kpp_1d_fields%ocntcorr(nzp1) )
    ALLOCATE( kpp_1d_fields%sal_clim(nzp1) ) 
    ALLOCATE( kpp_1d_fields%ocnt_clim(nzp1) )
    ALLOCATE( kpp_1d_fields%sflux(nsflxs,5,0:njdt) )
    ALLOCATE( kpp_1d_fields%talpha(0:nzp1tmax) )
    ALLOCATE( kpp_1d_fields%sbeta(0:nzp1tmax) ) 
    ALLOCATE( kpp_1d_fields%modeadv(maxmodeadv,2) )

  END SUBROUTINE mckpp_allocate_1d_fields


  SUBROUTINE mckpp_allocate_const_fields()
    
    ALLOCATE( kpp_const_fields%zm(nzp1) ) 
    ALLOCATE( kpp_const_fields%hm(nzp1) )
    ALLOCATE( kpp_const_fields%dm(0:nz) )
    ALLOCATE( kpp_const_fields%relax_sst_in(ny) )
    ALLOCATE( kpp_const_fields%relax_ocnt_in(ny) )  
    ALLOCATE( kpp_const_fields%relax_sal_in(ny) )
    ALLOCATE(kpp_const_fields%wmt(0:891,0:49))
    ALLOCATE(kpp_const_fields%wst(0:891,0:49))
    ALLOCATE(kpp_const_fields%tri(0:NZtmax,0:1,NGRID))
  
  END SUBROUTINE mckpp_allocate_const_fields

END MODULE mckpp_data_fields
