      LOGICAL L_RELAX_SAL, L_UPD_SAL, L_PERIODIC_SAL,L_INTERP_SAL
      INTEGER ndtupdsal, sal_period, ndt_interp_sal
      REAL relax_sal_in(NY)
      CHARACTER*200 sal_file
      COMMON /sal_in/ L_RELAX_SAL, L_UPD_SAL, ndtupdsal, sal_file,
     +     L_PERIODIC_SAL, sal_period, relax_sal_in,L_INTERP_SAL,
     +     ndt_interp_sal

      LOGICAL L_RELAX_OCNT,L_UPD_OCNT,L_PERIODIC_OCNT,L_INTERP_OCNT
      INTEGER ndtupdocnt, ocnt_period, ndt_interp_ocnt
      REAL relax_ocnT_in(NY)
      CHARACTER*200 ocnT_file
      COMMON /ocnT_in/ L_RELAX_OCNT,L_UPD_OCNT,ndtupdocnt,ocnt_period,
     +     L_PERIODIC_OCNT,ocnt_file,relax_ocnT_in,L_INTERP_OCNT,
     +     ndt_interp_ocnt

      LOGICAL L_RELAX_CURR,L_UPD_CURR,L_PERIODIC_CURR
      INTEGER ndtupdcurr, curr_period
      REAL relax_curr_in(NY)
      CHARACTER*200 curr_file
      COMMON /curr_in/ L_RELAX_CURR,L_UPD_CURR,L_PERIODIC_CURR,
     +     ndtupdcurr,curr_period,curr_file,relax_curr_in

      LOGICAL L_DAMP_CURR
      INTEGER dtuvdamp
      common /curr_clim/ L_DAMP_CURR, dtuvdamp

      LOGICAL L_RELAX_FILE
      CHARACTER*200 relax_file
      common /relax_file_in/ L_RELAX_FILE,relax_file
