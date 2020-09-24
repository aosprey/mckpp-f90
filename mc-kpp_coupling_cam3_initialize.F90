#include <misc.h>
#include <params.h>

SUBROUTINE mckpp_coupling_cam3_initialize(kpp_global_fields,kpp_const_fields)
  
  USE shr_kind_mod, only: r8=>shr_kind_r8
  USE mckpp_types,  only: kpp_global_type,kpp_const_type

  TYPE(kpp_global_type) :: kpp_global_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  

END SUBROUTINE mckpp_coupling_cam3_initialize
