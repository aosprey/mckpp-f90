      LOGICAL L_FCORR_WITHZ, L_FCORR, L_UPD_FCORR, L_PERIODIC_FCORR, 
     &     L_NO_FREEZE, L_NO_ISOTHERM, L_FCORR_NSOL, L_FCORR_NSOL_FILE,
     &     L_INTERP_FCORR
      INTEGER ndtupdfcorr, fcorr_period, isotherm_bottom,
     +     ndt_interp_fcorr
      REAL isotherm_threshold, fcorr_nsol_coeff
      CHARACTER*200 fcorrin_file, fcorr_nsol_file
      common /fcorr_in/ L_FCORR_WITHZ,L_FCORR,
     &     ndtupdfcorr, fcorrin_file,L_UPD_FCORR, L_PERIODIC_FCORR,
     &     fcorr_period, L_NO_FREEZE,
     &     L_NO_ISOTHERM, isotherm_bottom, isotherm_threshold,
     &     L_FCORR_NSOL, L_FCORR_NSOL_FILE, fcorr_nsol_file,
     &     fcorr_nsol_coeff, L_INTERP_FCORR, ndt_interp_fcorr
