#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE mckpp_output_control
  USE mckpp_types, only: kpp_const_fields
  USE pmgrid, only: masterproc
#else
SUBROUTINE mckpp_output_control(kpp_3d_fields,kpp_const_fields,kpp_timer)
#endif

  USE mckpp_xios_control

  IMPLICIT NONE

#ifndef MCKPP_CAM3
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  TYPE(kpp_timer_type) :: kpp_timer
#endif

  CALL mckpp_xios_output_control(kpp_3d_fields,kpp_const_fields) 
  
  RETURN
END SUBROUTINE mckpp_output_control
