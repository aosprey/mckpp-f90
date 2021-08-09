MODULE mckpp_time_control

  ! Routines to work with timestep counter, current time, and
  ! deriving validity times for I/O routines.
  ! Time variables all stored in kpp_const_fields 

USE mckpp_data_fields, ONLY: kpp_const_fields
  
IMPLICIT NONE

PUBLIC
  
CONTAINS

  ! Update ntime and time in kpp_const_fields
  SUBROUTINE mckpp_update_time(ntime)

    INTEGER, INTENT(IN) :: ntime ! new timestep

    kpp_const_fields%ntime = ntime
    kpp_const_fields%time = mckpp_get_time(ntime)

  END SUBROUTINE mckpp_update_time
  

  ! Return time for a given timestep
  REAL FUNCTION mckpp_get_time(ntime)

    INTEGER, INTENT(IN) :: ntime ! timestep

    mckpp_get_time = kpp_const_fields%startt + (ntime-1) * &
       kpp_const_fields%dto / kpp_const_fields%spd
  
  END FUNCTION mckpp_get_time  


END MODULE mckpp_time_control
