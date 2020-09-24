      LOGICAL L_RELAX_SAL, L_UPD_SAL, L_PERIODIC_SAL,L_INTERP_SAL
      INTEGER ndtupdsal, sal_period, ndt_interp_sal
#ifdef MCKPP_CAM3
      REAL relax_sal_in(NY_GLOBE)
#else
      REAL relax_sal_in(NY)
#endif
      CHARACTER*40 sal_file
      COMMON /sal_in/ L_RELAX_SAL, L_UPD_SAL, ndtupdsal, sal_file,L_PERIODIC_SAL, &
           sal_period, relax_sal_in,L_INTERP_SAL,ndt_interp_sal

      LOGICAL L_RELAX_OCNT,L_UPD_OCNT,L_PERIODIC_OCNT,L_INTERP_OCNT
      INTEGER ndtupdocnt, ocnt_period, ndt_interp_ocnt
#ifdef MCKPP_CAM3
      REAL relax_ocnT_in(NY_GLOBE)
#else
      REAL relax_ocnT_in(NY)
#endif
      CHARACTER*40 ocnT_file
      COMMON /ocnT_in/ L_RELAX_OCNT,L_UPD_OCNT,ndtupdocnt,ocnt_period,&
           L_PERIODIC_OCNT,ocnt_file,relax_ocnT_in,L_INTERP_OCNT,ndt_interp_ocnt
      
