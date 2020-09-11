      LOGICAL L_COUPLE,L_CPLWGHT,L_CLIMICE,L_OUTKELVIN,
     +     L_CLIM_ICE_DEPTH, L_PERIODIC_CLIMICE,
     +     L_CLIM_SNOW_ON_ICE,L_UPD_CLIMICE,L_COUPLE_CURRENTS,
     +     L_BAD_ICE_DEPTH,L_DIST_RUNOFF,L_SST_LAG_FUDGE,L_SST_LAG,
     +     L_SST_SMOOTH, L_SST_SMOOTH_Y, L_SST_SMOOTH_X,
     +     L_SST_SMOOTH_ANOM,L_SST_ANOM_FUDGE
      INTEGER ifirst,ilast,jfirst,jlast,ncid_cplwght,ndtupdice,
     +     climice_period,sst_lag_len,sst_smooth_ifirst,
     +     sst_smooth_ilast,sst_smooth_jfirst,sst_smooth_jlast,
     +     sst_smooth_blend
      CHARACTER*200 sstin_file,cplwght_file,icein_file,
     +     currin_file, initflux_file    
      common /couple/ ifirst,ilast,jfirst,jlast,sstin_file,
     +     icein_file,L_COUPLE,L_CLIMICE,L_CLIM_ICE_DEPTH,
     +     L_CLIM_SNOW_ON_ICE,ndtupdice,
     +     cplwght_file,ncid_cplwght,L_OUTKELVIN,
     +	   L_UPD_CLIMICE,L_PERIODIC_CLIMICE,L_COUPLE_CURRENTS,
     +     currin_file,climice_period,L_CPLWGHT,
     +     L_BAD_ICE_DEPTH,L_DIST_RUNOFF,initflux_file,sst_lag_len,
     +     L_SST_LAG_FUDGE,L_SST_LAG,L_SST_SMOOTH,L_SST_SMOOTH_Y,
     +     L_SST_SMOOTH_X,sst_smooth_ifirst,sst_smooth_ilast,
     +     sst_smooth_jfirst,sst_smooth_jlast,sst_smooth_blend,
     +     L_SST_SMOOTH_ANOM,L_SST_ANOM_FUDGE
