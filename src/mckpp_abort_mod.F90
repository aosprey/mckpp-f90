MODULE mckpp_abort_mod
  
  IMPLICIT NONE

CONTAINS 

SUBROUTINE MCKPP_ABORT()
  
#ifdef OASIS3
  ! Support for stopping in OASIS3
  ! NPK 2/11/09 - R3
  CALL mpi1_oasis3_terminate()
#else
  STOP
#endif
  
END SUBROUTINE MCKPP_ABORT

END MODULE mckpp_abort_mod
