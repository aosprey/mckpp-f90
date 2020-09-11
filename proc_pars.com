c      INTEGER jerlov(npts)
c      common/ proc_pars / jerlov
      
      CHARACTER*200 paras_file
      INTEGER ncid_paras,max_ekman_depth,max_ekadv_depth,ndtupdopt,
     $     opt_period
      LOGICAL L_JERLOV,L_VARY_OPT,L_PERIODIC_OPT
      common/ l_proc_pars / paras_file,ncid_paras,
     $     l_jerlov,max_ekman_depth,max_ekadv_depth,ndtupdopt,
     $     L_VARY_OPT,L_PERIODIC_OPT,opt_period
      
     
