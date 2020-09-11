      LOGICAL L_SFCORR_WITHZ, L_SFCORR, L_UPD_SFCORR, L_PERIODIC_SFCORR,
     +     L_INTERP_SFCORR
      INTEGER ndtupdsfcorr, sfcorr_period, ndt_interp_sfcorr
      CHARACTER*200 sfcorrin_file
      common /sfcorr_in/ L_SFCORR_WITHZ,L_SFCORR,
     &     ndtupdsfcorr, sfcorrin_file,L_UPD_SFCORR, L_PERIODIC_SFCORR,
     &     sfcorr_period, L_INTERP_SFCORR, ndt_interp_sfcorr

