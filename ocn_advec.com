      LOGICAL L_ADVECT
      CHARACTER*40 advect_file
      INTEGER ncid_advec
      common /ocn_advec/ L_ADVECT,advect_file,ncid_advec

      LOGICAL L_RELAX_SST,L_RELAX_CALCONLY
#ifdef MCKPP_CAM3
      REAL relax_sst_in(ny_globe)
#else
      REAL relax_sst_in(ny)
#endif

      common /ocn_relax/ relax_sst_in,L_RELAX_SST,L_RELAX_CALCONLY
