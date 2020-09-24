      CHARACTER*40 paras_file
      INTEGER ncid_paras
      LOGICAL L_JERLOV
      common/ l_proc_pars / paras_file,ncid_paras,l_jerlov
      
