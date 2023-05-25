MODULE mckpp_abort_mod

CONTAINS 

SUBROUTINE MCKPP_ABORT()
  
  IMPLICIT NONE
  
#ifdef MCKPP_COUPLE
#ifdef OASIS2
  ! Support for stopping in OASIS2
  ! NPK April 2009 - R2
  CALL halte('KPP : MIXED_ABORT called')
#endif
  ! Support for stopping in OASIS3
  ! NPK 2/11/09 - R3
#ifdef OASIS3
  CALL mpi1_oasis3_terminate()
#endif
#ifdef CFS
  ! Support for stopping with the CFS coupler
  ! Unsure how to stop the model for the GFS - Just stop?
  ! NPK June 2009 - R2     
  STOP
#endif /*CFS*/
#else
  STOP
#endif /*MCKPP_COUPLE*/
  
END SUBROUTINE MCKPP_ABORT

END MODULE mckpp_abort_mod
