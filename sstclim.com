      LOGICAL L_UPDCLIM, L_CLIMSST, L_UPD_CLIMSST, L_PERIODIC_CLIMSST
      INTEGER ndtupdsst, ifirst_sst, jfirst_sst, climsst_period
      common /sst_clim/ L_CLIMSST, L_UPD_CLIMSST, L_PERIODIC_CLIMSST, 
     +  ndtupdsst, ifirst_sst, jfirst_sst, climsst_period
