      LOGICAL L_LANDSEA
      INTEGER ncid_landsea
      CHARACTER*40 landsea_file
#ifdef MCKPP_COUPLE
#ifdef CFS
      REAL*4 global_longitudes(NX_GLOBE),global_latitudes(NY_GLOBE)
      REAL*4 global_lsm(NX_GLOBE,NY_GLOBE)
      COMMON /landsea/ ocdepth,L_LANDSEA,L_OCEAN,landsea_file,&
           ncid_landsea,global_longitudes,global_latitudes,global_lsm
#else
      COMMON /landsea/ L_LANDSEA,landsea_file
#endif /*CFS*/
#else
      COMMON /landsea/ L_LANDSEA,landsea_file,ncid_landsea
#endif /*COUPLE*/
