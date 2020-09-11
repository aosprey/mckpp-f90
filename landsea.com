      LOGICAL L_LANDSEA
      INTEGER ncid_landsea
      CHARACTER*200 landsea_file
      REAL ocdepth(NPTS)
#ifdef COUPLE
#ifdef CFS
      REAL*4 global_longitudes(NX_GLOBE),global_latitudes(NY_GLOBE)
      REAL*4 global_lsm(NX_GLOBE,NY_GLOBE)
      COMMON /landsea/ ocdepth,L_LANDSEA,L_OCEAN,landsea_file,
     +     ncid_landsea,global_longitudes,global_latitudes,global_lsm
#else
      COMMON /landsea/ ocdepth,L_LANDSEA,landsea_file,
     +     ncid_landsea        
#endif /*CFS*/
#else
      COMMON /landsea/ ocdepth,L_LANDSEA,landsea_file,
     +     ncid_landsea
#endif /*COUPLE*/
