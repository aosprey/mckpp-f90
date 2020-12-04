! Diagnostic output now fully handled by XIOS. Define diagnostics and output processing in XML file. 
! This just controls restart writing now. 
  LOGICAL L_RESTARTW
  INTEGER ndt_per_restart
  CHARACTER*17 restart_outfile
            
  common /output/ restart_outfile,l_restartw,ndt_per_restart
