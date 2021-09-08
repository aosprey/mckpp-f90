! XIOS code for reading/writing restarts and writing diagnostics.
! Used by mckpp_xios_control module which contains higher-level routines used by KPP model. 
MODULE mckpp_xios_io

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_log_messages, ONLY: mckpp_print_error, max_message_len, update_context
  USE mckpp_mpi_control, ONLY: comm
  USE mckpp_parameters, ONLY: nx, ny, nx_globe, ny_globe, npts, nz, nzp1, nsp1
  USE mckpp_time_control, ONLY: ntime, time

  USE xios 

  IMPLICIT NONE

  PUBLIC :: mckpp_xios_diagnostic_definition, mckpp_xios_restart_definition, &
       mckpp_xios_diagnostic_output, mckpp_xios_write_restart, mckpp_xios_read_restart

  PRIVATE

  REAL, ALLOCATABLE, DIMENSION(:) :: lons
  REAL, ALLOCATABLE, DIMENSION(:) :: lats
  REAL, ALLOCATABLE, DIMENSION(:) :: levs
  LOGICAL, ALLOCATABLE, DIMENSION(:) :: mask
  TYPE(xios_duration) :: dtime
  TYPE(xios_date) :: start_date
  
  TYPE(xios_context) :: ctx_hdl_diags

CONTAINS 

   !!! Routines for diagnostic output. 
  
  SUBROUTINE mckpp_xios_diagnostic_definition()

    CALL xios_context_initialize("kpp", comm)
    CALL xios_get_handle("kpp", ctx_hdl_diags)
    CALL xios_set_current_context(ctx_hdl_diags)

    CALL xios_define_calendar(type="Gregorian")

    CALL mckpp_xios_set_dimensions() 

    CALL xios_set_start_date(start_date) 
    CALL xios_set_time_origin(start_date) 
    CALL xios_set_timestep(timestep=dtime) 
    CALL xios_set_domain_attr("domain_kpp", type="rectilinear", & 
        data_dim=1, mask_1d=mask, &  
        ni_glo=NX, nj_glo=NY, & 
        lonvalue_1d=lons, latvalue_1d=lats)
    CALL xios_set_axis_attr("levels_kpp", n_glo=NZP1, value=levs) 

    CALL xios_close_context_definition()

  END SUBROUTINE mckpp_xios_diagnostic_definition


  SUBROUTINE mckpp_xios_set_dimensions() 

    dtime%second = kpp_const_fields%dto
    ! Work out date from days counter, and add 1 as 1st Jan is day 0. 
    start_date = xios_date(0000,01,01,00,00,00)+xios_day*(kpp_const_fields%startt+1)

    ALLOCATE( lons(nx), lats(ny), levs(nzp1), mask(npts) ) 
    lons = kpp_3d_fields%dlon(1:NX)
    lats = kpp_3d_fields%dlat(1::NX)
    levs = kpp_const_fields%zm
    mask = kpp_3d_fields%L_OCEAN

  END SUBROUTINE mckpp_xios_set_dimensions
  

  SUBROUTINE mckpp_xios_diagnostic_output()

    REAL, ALLOCATABLE, DIMENSION(:,:) :: temp_2d
    REAL, ALLOCATABLE, DIMENSION(:) :: temp_1d
    INTEGER :: k, ix, iy, ipt

    CALL xios_update_calendar(ntime)

    !!! Depth-varying diagnostics 
    ALLOCATE( temp_2d(npts, nzp1) ) 

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
    ALLOCATE( temp_1d(npts) ) 

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

  END SUBROUTINE mckpp_xios_diagnostic_output


  !!! Restart input/output. File and field definition all done in code. 

  SUBROUTINE mckpp_xios_write_restart(restart_time) 

    REAL, INTENT(IN) :: restart_time

    CHARACTER(LEN=9) :: restart_time_str
    CHARACTER(LEN=17) :: restart_filename
    CHARACTER(LEN=21) :: context_name
    TYPE(xios_context) :: ctx_hdl_restart
    CHARACTER(LEN=24) :: routine = "MCKPP_XIOS_WRITE_RESTART"

    WRITE(restart_time_str,'(F9.3)') restart_time
    restart_filename = "restart_"//TRIM(ADJUSTL(restart_time_str))
    context_name = "ctx_restart_"//TRIM(ADJUSTL(restart_time_str))

    ! Define a new context for this restart file 
    CALL xios_context_initialize(context_name, comm)
    CALL xios_get_handle(context_name, ctx_hdl_restart)
    CALL xios_set_current_context(ctx_hdl_restart)

    ! Define file  
    CALL mckpp_xios_restart_definition(routine, restart_filename, "write") 

    ! Write fields 
    CALL mckpp_xios_restart_output() 

    ! Now close context and switch back to diagnostic one 
    CALL xios_context_finalize()
    CALL xios_set_current_context(ctx_hdl_diags)

  END SUBROUTINE mckpp_xios_write_restart

  
  SUBROUTINE mckpp_xios_read_restart(restart_filename) 

    CHARACTER(LEN=17), INTENT(IN) :: restart_filename
    
    CHARACTER(LEN=16) :: context_name = "ctx_restart_read"
    TYPE(xios_context) :: ctx_hdl_restart
    CHARACTER(LEN=23) :: routine = "MCKPP_XIOS_READ_RESTART"


    ! Define a new context for this restart file 
    CALL xios_context_initialize(context_name, comm)
    CALL xios_get_handle(context_name, ctx_hdl_restart)
    CALL xios_set_current_context(ctx_hdl_restart)

    ! Define file  
    CALL mckpp_xios_restart_definition(routine, restart_filename, "read") 

    ! Read fields 
    CALL mckpp_xios_restart_input() 

    ! Now close context
    CALL xios_context_finalize()

  END SUBROUTINE mckpp_xios_read_restart


  ! We can use the same defintion code for reading or writing restart, with a couple of small changes
  SUBROUTINE mckpp_xios_restart_definition(calling_routine, filename, mode) 

    CHARACTER(*), INTENT(IN) :: calling_routine, filename, mode

    LOGICAL :: read_mode
    TYPE(xios_filegroup) :: filedefn_hdl 
    TYPE(xios_file) :: file_hdl
    TYPE(xios_axisgroup) :: axisdefn_hdl
    TYPE(xios_axis) :: axis_hdl
    TYPE(xios_gridgroup) :: griddefn_hdl
    TYPE(xios_grid) :: grid_hdl
    TYPE(xios_domaingroup) :: domaindefn_hdl
    TYPE(xios_domain) :: domain_hdl
    TYPE(xios_scalargroup) :: scalardefn_hdl
    TYPE(xios_scalar) :: scalar_hdl
    TYPE(xios_field) :: field_hdl   
    TYPE(xios_date) :: restart_date
    CHARACTER(LEN=max_message_len) :: context, message
    CHARACTER(LEN=29) :: routine = "MCKPP_XIOS_RESTART_DEFINITION"

    ! Check if read or write mode
    IF (mode .EQ. "read") THEN
      read_mode = .TRUE.
    ELSE IF (mode .EQ. "write") THEN
      read_mode = .FALSE.
    ELSE
      context = update_context(calling_routine, routine) 
      WRITE(message,*) "Routine called with mode = ", mode, ", but must be read or write"
      CALL mckpp_abort(context, message)
    END IF  

    ! Define calendar and ts
    CALL xios_define_calendar(type="Gregorian") 
    restart_date = xios_date(0000,01,01,00,00,00)+xios_day*(time+1)
    CALL xios_set_start_date(restart_date) 
    CALL xios_set_time_origin(restart_date) 
    dtime%second = kpp_const_fields%dto
    CALL xios_set_timestep(timestep=dtime) 

    ! Define axes and grids - no landsea mask is applied
    CALL xios_get_handle("domain_definition", domaindefn_hdl) 
    CALL xios_add_child(domaindefn_hdl, domain_hdl, "domain_kpp_nomask")
    CALL xios_set_domain_attr("domain_kpp_nomask", type="rectilinear", & 
        data_dim=1, ni_glo=NX, nj_glo=NY, lon_name="longitude", lat_name="latitude")
    IF (read_mode) THEN
      CALL xios_set_domain_attr("domain_kpp_nomask", ni=nx, nj=ny, ibegin=0, jbegin=0)
    ELSE
      CALL xios_set_domain_attr("domain_kpp_nomask", lonvalue_1d=lons, latvalue_1d=lats)
    END IF 

    CALL xios_get_handle("axis_definition", axisdefn_hdl) 

    CALL xios_add_child(axisdefn_hdl, axis_hdl, "levels_kpp") 
    CALL xios_set_axis_attr("levels_kpp", n_glo=NZP1, value=levs, name="z") 

    CALL xios_add_child(axisdefn_hdl, axis_hdl, "intcnt_kpp") 
    CALL xios_set_axis_attr("intcnt_kpp", n_glo=2, value=(/1.,2./))

    CALL xios_get_handle("grid_definition", griddefn_hdl) 

    CALL xios_add_child(griddefn_hdl, grid_hdl, "grid_kpp_3d_nomask") 
    CALL xios_add_child(grid_hdl, domain_hdl, "domain_kpp_nomask") 
    CALL xios_add_child(grid_hdl, axis_hdl, "levels_kpp") 

    CALL xios_add_child(griddefn_hdl, grid_hdl, "grid_kpp_2d_nomask") 
    CALL xios_add_child(grid_hdl, domain_hdl, "domain_kpp_nomask") 

    CALL xios_add_child(griddefn_hdl, grid_hdl, "grid_kpp_2d_nomask_intcnt") 
    CALL xios_add_child(grid_hdl, domain_hdl, "domain_kpp_nomask") 
    CALL xios_add_child(grid_hdl, axis_hdl, "intcnt_kpp") 

    CALL xios_add_child(griddefn_hdl, grid_hdl, "grid_kpp_3d_nomask_intcnt") 
    CALL xios_add_child(grid_hdl, domain_hdl, "domain_kpp_nomask") 
    CALL xios_add_child(grid_hdl, axis_hdl, "levels_kpp") 
    CALL xios_add_child(grid_hdl, axis_hdl, "intcnt_kpp") 

    CALL xios_get_handle("scalar_definition", scalardefn_hdl) 
    CALL xios_add_child(scalardefn_hdl, scalar_hdl, "grid_scalar") 

    ! Define restart file
    CALL xios_get_handle("file_definition", filedefn_hdl)
    CALL xios_add_child(filedefn_hdl, file_hdl, "restart") 
    CALL xios_set_file_attr("restart", name=filename, output_freq=xios_timestep, & 
        type="one_file", par_access="collective")
    IF (read_mode) THEN
      CALL xios_set_file_attr("restart", mode="read")
    END IF

    ! Define variables to write to restart
    CALL xios_add_child(file_hdl, field_hdl, "time")
    CALL xios_set_attr(field_hdl, name="time", long_name="Time", unit="days", &  
        scalar_ref="grid_scalar", operation="instant") 

    CALL mckpp_xios_restart_define_field(file_hdl, "uvel", "Zonal velocity", "m/s", "grid_kpp_3d_nomask") 
    CALL mckpp_xios_restart_define_field(file_hdl, "vvel", "Meridional velocity", "m/s", "grid_kpp_3d_nomask") 
    CALL mckpp_xios_restart_define_field(file_hdl, "T", "Temperature", "degC", "grid_kpp_3d_nomask")
    CALL mckpp_xios_restart_define_field(file_hdl, "S", "Salinity", "o/oo", "grid_kpp_3d_nomask") 
    CALL mckpp_xios_restart_define_field(file_hdl, "CP", "Specific heat capacity", "J/kg/K", "grid_kpp_3d_nomask") 
    CALL mckpp_xios_restart_define_field(file_hdl, "rho", "Density", "kg/m^3", "grid_kpp_3d_nomask") 
    CALL mckpp_xios_restart_define_field(file_hdl, "hmix", "Mixed layer depth", "m", "grid_kpp_2d_nomask") 
    CALL mckpp_xios_restart_define_field(file_hdl, "kmix", "Mixed layer depth", "m", "grid_kpp_2d_nomask") 
    CALL mckpp_xios_restart_define_field(file_hdl, "Sref", "Reference salinity", "o/oo", "grid_kpp_2d_nomask") 
    CALL mckpp_xios_restart_define_field(file_hdl, "SSref", "Reference surface salinity", "o/oo", "grid_kpp_2d_nomask") 
    CALL mckpp_xios_restart_define_field(file_hdl, "Ssurf", "Surafce salinity", "o/oo", "grid_kpp_2d_nomask") 
    CALL mckpp_xios_restart_define_field(file_hdl, "Tref", "Reference temperature", "degC", "grid_kpp_2d_nomask") 
    CALL mckpp_xios_restart_define_field(file_hdl, "old", "Integration counter", "unitless", "grid_kpp_2d_nomask") 
    CALL mckpp_xios_restart_define_field(file_hdl, "new", "Integration counter", "unitless", "grid_kpp_2d_nomask") 
    CALL mckpp_xios_restart_define_field(file_hdl, "Us", "Zonal velocity in integration", "m/s", "grid_kpp_3d_nomask_intcnt") 
    CALL mckpp_xios_restart_define_field(file_hdl, "Vs", "Meridional velocity in integration", "m/s", "grid_kpp_3d_nomask_intcnt")
    CALL mckpp_xios_restart_define_field(file_hdl, "Ts", "Teemperature in integration", "degC", "grid_kpp_3d_nomask_intcnt")
    CALL mckpp_xios_restart_define_field(file_hdl, "Ss", "Salinity in integration", "o/oo", "grid_kpp_3d_nomask_intcnt")
    CALL mckpp_xios_restart_define_field(file_hdl, "hmixd", "Mixed layer depth in integration", "m", "grid_kpp_2d_nomask_intcnt")

    CALL xios_close_context_definition()

  END SUBROUTINE mckpp_xios_restart_definition


  SUBROUTINE mckpp_xios_restart_define_field(file_hdl, name, long_name, unit, grid) 

    TYPE(xios_file), INTENT(IN) :: file_hdl 
    CHARACTER(*), INTENT(IN) :: name, long_name, unit, grid
    TYPE(xios_field) :: field_hdl   

    CALL xios_add_child(file_hdl, field_hdl, name)
    CALL xios_set_attr(field_hdl, name=name, long_name=long_name, & 
        unit=unit, grid_ref=grid, operation="instant") 

  END SUBROUTINE mckpp_xios_restart_define_field


  SUBROUTINE mckpp_xios_restart_output()

    ! Set correct time for validity of restart fields 
    ! (end of this timestep = start of next timestep)
    CALL xios_update_calendar(ntime+1)

    CALL xios_send_field("time", time+kpp_const_fields%dto/kpp_const_fields%spd) 
    CALL xios_send_field("uvel", kpp_3d_fields%U(:,:,1)) 
    CALL xios_send_field("vvel", kpp_3d_fields%U(:,:,2))
    CALL xios_send_field("T", kpp_3d_fields%X(:,:,1)) 
    CALL xios_send_field("S", kpp_3d_fields%X(:,:,2))
    CALL xios_send_field("CP", kpp_3d_fields%cp(:,1:NZP1))
    CALL xios_send_field("rho", kpp_3d_fields%rho(:,1:NZP1))
    CALL xios_send_field("hmix", kpp_3d_fields%hmix)
    CALL xios_send_field("kmix", kpp_3d_fields%kmix)
    CALL xios_send_field("Sref", kpp_3d_fields%Sref)
    CALL xios_send_field("SSref", kpp_3d_fields%SSref)
    CALL xios_send_field("Ssurf", kpp_3d_fields%Ssurf)
    CALL xios_send_field("Tref", kpp_3d_fields%Tref)
    CALL xios_send_field("old", REAL(kpp_3d_fields%old))
    CALL xios_send_field("new", REAL(kpp_3d_fields%new))
    CALL xios_send_field("Us", kpp_3d_fields%Us(:,:,1,0:1))
    CALL xios_send_field("Vs", kpp_3d_fields%Us(:,:,2,0:1))
    CALL xios_send_field("Ts", kpp_3d_fields%Xs(:,:,1,0:1))
    CALL xios_send_field("Ss", kpp_3d_fields%Xs(:,:,2,0:1))
    CALL xios_send_field("hmixd", kpp_3d_fields%hmixd(:,0:1))

  END SUBROUTINE mckpp_xios_restart_output

  
  SUBROUTINE mckpp_xios_restart_input()

    REAL, DIMENSION(npts) :: tmp

    ! Assume only reading restart once at start of run
    CALL xios_update_calendar(0)

    CALL xios_recv_field("uvel", kpp_3d_fields%U(:,:,1)) 
    CALL xios_recv_field("vvel", kpp_3d_fields%U(:,:,2))
    CALL xios_recv_field("T", kpp_3d_fields%X(:,:,1)) 
    CALL xios_recv_field("S", kpp_3d_fields%X(:,:,2))
    CALL xios_recv_field("CP", kpp_3d_fields%cp(:,1:NZP1))
    CALL xios_recv_field("rho", kpp_3d_fields%rho(:,1:NZP1))
    CALL xios_recv_field("hmix", kpp_3d_fields%hmix)
    CALL xios_recv_field("kmix", kpp_3d_fields%kmix)
    CALL xios_recv_field("Sref", kpp_3d_fields%Sref)
    CALL xios_recv_field("SSref", kpp_3d_fields%SSref)
    CALL xios_recv_field("Ssurf", kpp_3d_fields%Ssurf)
    CALL xios_recv_field("Tref", kpp_3d_fields%Tref)
    CALL xios_recv_field("old", tmp)
    kpp_3d_fields%old = INT(tmp)
    CALL xios_recv_field("new", tmp)
    kpp_3d_fields%new = INT(tmp)
    CALL xios_recv_field("Us", kpp_3d_fields%Us(:,:,1,0:1))
    CALL xios_recv_field("Vs", kpp_3d_fields%Us(:,:,2,0:1))
    CALL xios_recv_field("Ts", kpp_3d_fields%Xs(:,:,1,0:1))
    CALL xios_recv_field("Ss", kpp_3d_fields%Xs(:,:,2,0:1))
    CALL xios_recv_field("hmixd", kpp_3d_fields%hmixd(:,0:1))
    
  END SUBROUTINE mckpp_xios_restart_input 

  
END MODULE mckpp_xios_io
