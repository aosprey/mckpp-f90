      LOGICAL L_INITDATA,L_INTERPINIT,L_RESTART
      CHARACTER*50 initdata_file
      CHARACTER*17 restart_infile
      common /initdata/ L_INITDATA,L_INTERPINIT,L_RESTART,initdata_file,restart_infile
      
