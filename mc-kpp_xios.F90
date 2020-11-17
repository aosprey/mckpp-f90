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
  INTEGER, PARAMETER :: my_nx=NX, my_ny=NY, my_nzp1=NZP1
  REAL, DIMENSION(my_nx) :: lons
  REAL, DIMENSION(my_ny) :: lats
  REAL, DIMENSION(my_nzp1) :: levs
  LOGICAL, DIMENSION(my_nx*my_ny) :: mask


  ! This should go elsewhere but leave here for now. 
  CALL mpi_init(ierr) 

  CALL xios_initialize("client", return_comm=comm)
  CALL xios_context_initialize("kpp", comm)

  CALL xios_get_handle("kpp", ctx_hdl)
  CALL xios_set_current_context(ctx_hdl)

  CALL xios_define_calendar(type="Gregorian", &
                            start_date=xios_date(2015, 11, 01, 00, 00, 00), &
                            time_origin=xios_date(2015, 11, 01, 00, 00, 00))

  lons = kpp_3d_fields%dlon(1:NX)
  lats = kpp_3d_fields%dlat(1::NX)
  levs = kpp_const_fields%zm
  mask = kpp_3d_fields%L_OCEAN

  CALL xios_set_domain_attr("domain_kpp", type="rectilinear", & 
                            ni_glo=NX, nj_glo=NY, & 
                            data_dim=1, &  
                            lonvalue_1d=lons, latvalue_1d=lats, & 
                            mask_1d=mask)
  CALL xios_set_axis_attr("axis_kpp", n_glo=NZP1, value=levs) 


  dtime%second=kpp_const_fields%dto
  CALL xios_set_timestep(timestep=dtime)

  CALL xios_close_context_definition()

END SUBROUTINE mckpp_initialize_xios


SUBROUTINE mckpp_output_xios(kpp_3d_fields, kpp_const_fields)

USE xios 

#include <mc-kpp_3d_type.com>

  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER, PARAMETER :: my_nx=NX, my_ny=NY, my_nzp1=NZP1
  REAL, DIMENSION(my_nx*my_ny,nzp1) :: temp_2d
  REAL, DIMENSION(my_nx*my_ny) :: temp_1d
  INTEGER :: k, ix, iy

  CALL xios_update_calendar(kpp_const_fields%ntime)

  !!! Depth-varying diagnostics 

  ! Zonal current
  CALL xios_send_field("u", kpp_3d_fields%U(:,:,1))

  ! Meridional current 
  CALL xios_send_field("v", kpp_3d_fields%U(:,:,2))

  ! Temperature 
  CALL xios_send_field("T", kpp_3d_fields%X(:,:,1)) 

  ! Salinity
  DO k=1,NZP1
    temp_2d(:,k) = kpp_3d_fields%X(:,k,2)+kpp_3d_fields%Sref(:)
  END DO
  CALL xios_send_field("S", temp_2d)

  ! Buoyancy
  CALL xios_send_field("B", kpp_3d_fields%buoy(:,1:NZP1))

  ! Turbulent zonal velocity flux
  CALL xios_send_field("wu", kpp_3d_fields%wU(:,0:NZ,1))

  ! Turbulent meridional velocity flux
  CALL xios_send_field("wv", kpp_3d_fields%wU(:,0:NZ,2))

  ! Turbulent temperature flux
  CALL xios_send_field("wT", kpp_3d_fields%wX(:,0:NZ,1))

  ! Turbulent salinity flux
  CALL xios_send_field("wS", kpp_3d_fields%wX(:,0:NZ,2))

  ! Turbulent buoyancy flux
  CALL xios_send_field("wB", kpp_3d_fields%wX(:,0:NZ,NSP1))

  ! Non-turbulent temperature flux
  CALL xios_send_field("wTnt", kpp_3d_fields%wXNT(:,0:NZ,1))

  ! Diffusion coefficient for momentum
  temp_2d(:,1)=0.0
  temp_2d(:,2:NZP1)=kpp_3d_fields%difm(:,1:NZ)
  CALL xios_send_field("difm", temp_2d)

  ! Diffusion coefficient for temperature
  temp_2d(:,1)=0.0
  temp_2d(:,2:NZP1)=kpp_3d_fields%dift(:,1:NZ)
  CALL xios_send_field("dift", temp_2d)

  ! Diffusion coefficient for salinity
  temp_2d(:,1)=0.0
  temp_2d(:,2:NZP1)=kpp_3d_fields%difs(:,1:NZ)
  CALL xios_send_field("difs", temp_2d)

  ! Density
  CALL xios_send_field("rho", kpp_3d_fields%rho(:,1:NZP1))

  ! Specific heat capacity
  CALL xios_send_field("cp", kpp_3d_fields%cp(:,1:NZP1))

  ! Salinity correction 
  CALL xios_send_field("scorr", kpp_3d_fields%scorr)

  ! Local Richardson number
  CALL xios_send_field("Rig", kpp_3d_fields%Rig)

  ! Local delta buoyancy
  temp_2d(:,1:NZ)=kpp_3d_fields%dbloc(:,1:NZ)
  temp_2d(:,NZP1)=0.0
  CALL xios_send_field("dbloc", temp_2d)

  ! Local shear-squared terms
  CALL xios_send_field("Shsq", kpp_3d_fields%Shsq)

  ! Temperature increment 
  CALL xios_send_field("tinc_fcorr", kpp_3d_fields%Tinc_fcorr)

  ! Temperature correction
  CALL xios_send_field("fcorr_z", kpp_3d_fields%ocnTcorr)

  ! Salinity increment
  CALL xios_send_field("sinc_fcorr", kpp_3d_fields%Sinc_fcorr)

  !!! Single-level diagnostics 

  ! Mixed layer depth  
  CALL xios_send_field("hmix", kpp_3d_fields%hmix) 

  ! Flux correction 
  CALL xios_send_field("fcorr", kpp_3d_fields%fcorr)

  ! Zonal wind stress input
  CALL xios_send_field("taux_in", kpp_3d_fields%sflux(:,1,5,0))

  ! Meridional wind stress input
  CALL xios_send_field("tauy_in", kpp_3d_fields%sflux(:,2,5,0))

  ! Shortwave radiation input 
  CALL xios_send_field("solar_in", kpp_3d_fields%sflux(:,3,5,0))

  ! Non-shortwave radiation input 
  CALL xios_send_field("nsolar_in", kpp_3d_fields%sflux(:,4,5,0))

  ! Net freshwater input 
  CALL xios_send_field("PminusE_in", kpp_3d_fields%sflux(:,6,5,0))

  ! Coupling weight 
  DO ix=kpp_const_fields%ifirst,kpp_const_fields%ilast
    DO iy=kpp_const_fields%jfirst,kpp_const_fields%jlast
      ipt=(iy-1)*NX_GLOBE+ix
      temp_1d((iy-kpp_const_fields%jfirst)*NX+ix-kpp_const_fields%ifirst+1)=kpp_3d_fields%cplwght(ipt)
    ENDDO
  ENDDO
  CALL xios_send_field("cplwght", temp_1d)

  ! Fraction of points below freezing 
  CALL xios_send_field("freeze_flag", kpp_3d_fields%freeze_flag)

  ! Flag for isothermal detection 
  CALL xios_send_field("comp_flag", kpp_3d_fields%reset_flag)

  ! Fraction of points with ui~u**2
  CALL xios_send_field("dampu_flag", kpp_3d_fields%dampu_flag)

  ! Fraction of points with  vi~v**2
  CALL xios_send_field("dampv_flag", kpp_3d_fields%dampv_flag)

END SUBROUTINE mckpp_output_xios


SUBROUTINE mckpp_finalize_xios() 

USE mpi
USE xios 

  CALL xios_context_finalize()
  CALL xios_finalize()

  CALL mpi_finalize(ierr)


END SUBROUTINE mckpp_finalize_xios
