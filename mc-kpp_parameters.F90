MODULE mckpp_parameters

! This can be read in at runtime if we have dynamic arrays. 
#include <parameter.inc>

  INTEGER, PARAMETER :: max_nc_filename_len = 30 
  INTEGER, PARAMETER :: max_restart_filename_len = 17 

  INTEGER, PARAMETER :: nuout = 6, nuerr = 5

END MODULE mckpp_parameters
