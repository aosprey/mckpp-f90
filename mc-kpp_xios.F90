SUBROUTINE mckpp_initialize_xios(kpp_3d_fields, kpp_const_fields)

USE mpi
USE xios 

IMPLICIT NONE 

#include <mc-kpp_3d_type.com>

  INTEGER :: ierr, comm 
  TYPE(xios_duration) :: dtime
  TYPE(xios_context) :: ctx_hdl
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER, PARAMETER :: my_nx=nx, my_ny=ny, my_nzp1=nzp1
  REAL, DIMENSION(my_nx) :: lons
  REAL, DIMENSION(my_ny) :: lats
  REAL, DIMENSION(my_nzp1) :: levs
  LOGICAL, DIMENSION(my_nx*my_ny) :: mask


  ! This should go elsewhere but leave here for now. 
  CALL mpi_init(ierr) 

  CALL xios_initialize("client", return_comm=comm)
  CALL xios_context_initialize("kpp",comm)

  CALL xios_get_handle("kpp",ctx_hdl)
  CALL xios_set_current_context(ctx_hdl)

  CALL xios_define_calendar(type="Gregorian", &
                            start_date=xios_date(2015, 11, 01, 00, 00, 00), &
                            time_origin=xios_date(2015, 11, 01, 00, 00, 00))

  lons = kpp_3d_fields%dlon(1:NX)
  lats = kpp_3d_fields%dlat(1::NX)
  levs = kpp_const_fields%zm
  mask = kpp_3d_fields%L_OCEAN

  CALL xios_set_domain_attr("domain_kpp",type="rectilinear", & 
                            ni_glo=nx, nj_glo=ny, & 
                            data_dim=1, &  
                            lonvalue_1d=lons, latvalue_1d=lats, & 
                            mask_1d=mask)
  CALL xios_set_axis_attr("axis_kpp",n_glo=nzp1 ,value=levs) 


  dtime%second=kpp_const_fields%dto
  CALL xios_set_timestep(timestep=dtime)

  CALL xios_close_context_definition()

END SUBROUTINE mckpp_initialize_xios


SUBROUTINE mckpp_output_xios(kpp_3d_fields, kpp_const_fields)

USE xios 

#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  CALL xios_update_calendar(kpp_const_fields%ntime)

  CALL xios_send_field("T", kpp_3d_fields%X(:,:,1)) 

END SUBROUTINE mckpp_output_xios


SUBROUTINE mckpp_finalize_xios() 

USE mpi
USE xios 

  CALL xios_context_finalize()
  CALL xios_finalize()

  CALL mpi_finalize(ierr)


END SUBROUTINE mckpp_finalize_xios
