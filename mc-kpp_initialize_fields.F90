#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif

SUBROUTINE MCKPP_INITIALIZE_FIELDS()
  
#ifdef MCKPP_CAM3
  USE shr_kind_mod, only : r8=>shr_kind_r8
  USE mckpp_types, only : kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE phys_grid, only : get_ncols_p,get_rlat_all_p,get_rlon_all_p
  USE ppgrid, only : pcols,begchunk,endchunk
  USE pmgrid, only : masterproc
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields, &
      mckpp_allocate_3d_fields
#endif
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_warning, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, npts

  IMPLICIT NONE

#ifdef MCKPP_CAM3
  INTEGER :: icol,ncol,ichnk
  REAL(r8) :: clat1(pcols),clon1(pcols)
#endif
  INTEGER :: iy,ix,ipt
  CHARACTER(LEN=23) :: routine = "MCKPP_INITIALIZE_FIELDS"
  CHARACTER(LEN=max_message_len) :: message 

#ifndef MCKPP_CAM3
  CALL mckpp_allocate_3d_fields()
#endif

  ! Set initial values for flags in kpp_3d_fields, which otherwise
  ! might never be set if points are not coupled.
  kpp_3d_fields%dampu_flag(:)=0.
  kpp_3d_fields%dampv_flag(:)=0.
  kpp_3d_fields%freeze_flag(:)=0.
  kpp_3d_fields%reset_flag(:)=0.
   
  ! Initialize cplwght. 
  kpp_3d_fields%cplwght(:)=0.
  
  ! Initialize latitude and longitude areas and the land/sea mask
#ifdef MCKPP_CAM3
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     CALL GET_RLAT_ALL_P(ichnk,ncol,clat1)
     CALL GET_RLON_ALL_P(ichnk,ncol,clon1)
     kpp_3d_fields(ichnk)%dlat(:)=clat1*360./kpp_const_fields%twopi
     kpp_3d_fields(ichnk)%dlon(:)=clon1*360./kpp_const_fields%twopi    
  ENDDO
  IF (kpp_const_fields%L_LANDSEA) THEN
     CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_LANDSEA")
     CALL MCKPP_INITIALIZE_LANDSEA()
     CALL mckpp_print(routine, "Returned from MCKPP_INITIALIZE_LANDSEA")
  ENDIF
#elif CFS
  IF (L_LANDSEA) CALL read_landsea_global()
#else
  IF (kpp_const_fields%L_LANDSEA) THEN
     CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_LANDSEA")
     kpp_3d_fields%dlat(1)=kpp_const_fields%alat
     kpp_3d_fields%dlon(1)=kpp_const_fields%alon
     CALL MCKPP_INITIALIZE_LANDSEA()
  ELSEIF (kpp_const_fields%L_REGGRID) THEN
     DO iy=1,ny
        DO ix=1,nx
           ipt=(iy-1)*nx+ix
           kpp_3d_fields%dlat(ipt)=kpp_const_fields%alat+(iy-1)*kpp_const_fields%delta_lat
           kpp_3d_fields%dlon(ipt)=kpp_const_fields%alon+(ix-1)*kpp_const_fields%delta_lon
           kpp_3d_fields%ocdepth(ipt)=-10000.
           kpp_3d_fields%L_OCEAN(ipt)=.TRUE.
        ENDDO
     ENDDO
   ELSEIF (.NOT. kpp_const_fields%L_REGGRID .AND. .NOT. kpp_const_fields%L_LANDSEA) THEN
     message = "If you set L_REGGRID=.FALSE., you must specify a land-sea mask file from which" &
         // " to read the locations of the gridpoints in the horizontal."
     CALL mckpp_warning(routine, message)
  ENDIF
#endif

  ! Initialize the vertical grid
  CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_GEOGRAPHY")
  CALL MCKPP_INITIALIZE_GEOGRAPHY()
  CALL mckpp_print(routine, "Returned from MCKPP_INITIALIZE_GEOGRAPHY")

  ! Initialize coupling weights   
#if (defined OASIS2 || defined OASIS3 || defined CFS)
  CALL MCKPP_INITIALIZE_COUPLINGWEIGHT()
#elif (defined MCKPP_CAM3)
  CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_COUPLINGWEIGHT") 
  CALL MCKPP_INITIALIZE_COUPLINGWEIGHT()
  CALL mckpp_print(routine, "Returned from MCKPP_INITIALIZE_COUPLINGWEIGHT")
  ! Initialize coupling-period mean fluxes and SST fields
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%sflux_cpl(:,:)=0
     kpp_3d_fields(ichnk)%sst_cpl(:)=0
  ENDDO
#else
  IF (kpp_const_fields%L_CPLWGHT) CALL MCKPP_INITIALIZE_COUPLINGWEIGHT()
#endif

  ! Initialize advection options
  IF (kpp_const_fields%L_ADVECT) THEN
     CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_ADVECTION")
     CALL MCKPP_INITIALIZE_ADVECTION()
     CALL mckpp_print(routine, "Returned from MCKPP_INITIALIZE_ADVECTION")
  ELSE
#ifdef MCKPP_CAM3
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%nmodeadv(1:ncol,:)=0
      ENDDO
#else
     DO ipt=1,npts
        kpp_3d_fields%nmodeadv(ipt,1)=0
        kpp_3d_fields%nmodeadv(ipt,2)=0
     ENDDO
#endif 
     CALL mckpp_print(routine, "No advection has been specified.") 
   ENDIF

  ! Initialize relaxation of SST, temperature and/or salinity
  IF (kpp_const_fields%L_RELAX_SST .OR. kpp_const_fields%L_RELAX_SAL &
       .OR. kpp_const_fields%L_RELAX_OCNT) THEN     
     CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_RELAXATION")
     CALL MCKPP_INITIALIZE_RELAXATION()
     CALL mckpp_print(routine, "Returned from MCKPP_INITIALIZE_RELAXATION")
  ENDIF

  ! Initialize water type for optical properties of seawater
  CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_OPTICS")
  CALL MCKPP_INITIALIZE_OPTICS()
  CALL mckpp_print(routine, "Returned from MCKPP_INITIALIZE_OPTICS")

  ! Initialize ocean profiles
  IF (kpp_const_fields%L_RESTART) THEN     
     CALL mckpp_print(routine, "Calling MCKPP_RESTART_IO_READ")
     ! Still needs scattering code
     CALL MCKPP_RESTART_IO_READ_NETCDF()
     CALL mckpp_print(routine, "Returned from MCKPP_RESTART_IO_READ")
  ELSE
     CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_OCEAN_PROFILES")
     CALL MCKPP_INITIALIZE_OCEAN_PROFILES()
     CALL mckpp_print(routine, "Returned from MCKPP_INITIALIZE_OCEAN_PROFILES")
  ENDIF
  
  ! Initialize boundary conditions
  IF (kpp_const_fields%L_CLIMSST) THEN     
     CALL mckpp_print(routine, "Calling MCKPP_READ_SST")
     CALL MCKPP_READ_SST()
     CALL mckpp_print(routine, "Returned from MCKPP_READ_SST")
  ENDIF

  IF (kpp_const_fields%L_CLIMICE) THEN    
     CALL mckpp_print(routine, "Calling MCKPP_READ_ICE")
     CALL MCKPP_READ_ICE()
     CALL mckpp_print(routine, "Returned from MCKPP_READ_ICE")
  ENDIF
!!$    !IF (L_CLIMCURR) CALL read_surface_currents(kpp_3d_fields,kpp_const_fields)
  
  IF (kpp_const_fields%L_FCORR_WITHZ) THEN     
     CALL mckpp_print(routine, "Calling MCKPP_READ_FCORR_3D")
     CALL MCKPP_READ_FCORR_3D()
     CALL mckpp_print(routine, "Returned from MCKPP_READ_FCORR_3D")
  ELSEIF (kpp_const_fields%L_FCORR) THEN      
     CALL mckpp_print(routine, "Calling MCKPP_READ_FCORR_2D")
     CALL MCKPP_READ_FCORR_2D()
     CALL mckpp_print(routine, "Returned from MCKPP_READ_FCORR_2D")
  ENDIF

  IF (kpp_const_fields%L_SFCORR_WITHZ) THEN   
     CALL mckpp_print(routine, "Calling MCKPP_READ_SFCORR_3D")
     CALL MCKPP_READ_SFCORR_3D()
     CALL mckpp_print(routine, "Returned from MCKPP_READ_SFCORR_3D")
  ELSEIF (kpp_const_fields%L_SFCORR) THEN
     CALL mckpp_print(routine, "Calling MCKPP_READ_SFCORR_2D")
     CALL MCKPP_READ_SFCORR_2D()
     CALL mckpp_print(routine, "Returned from MCKPP_READ_SFCORR_2D")
  ENDIF

  IF (kpp_const_fields%L_VARY_BOTTOM_TEMP) THEN     
     CALL mckpp_print(routine, "Calling MCKPP_READ_TEMPERATURES_BOTTOM")
     CALL MCKPP_READ_TEMPERATURES_BOTTOM()
     CALL mckpp_print(routine, "Returned from MCKPP_READ_TEMPERATURES_BOTTOM")
  ENDIF

  !!! THESE SHOULD NOT BE CALLED IF DOING L_INTERP_OCNT / L_INTERP_SAL !!!

  IF (kpp_const_fields%L_RELAX_OCNT .AND. .NOT. kpp_const_fields%L_INTERP_OCNT) THEN 
     CALL mckpp_print(routine, "Calling MCKPP_READ_TEMPERATURES_3D")
     CALL MCKPP_READ_TEMPERATURES_3D()
     CALL mckpp_print(routine, "Returned from MCKPP_READ_TEMPERATURES_3D")
  ELSEIF (kpp_const_fields%L_RELAX_OCNT .AND. kpp_const_fields%L_INTERP_OCNT) THEN     
     CALL mckpp_print(routine, "Calling MCKPP_BOUNDARY_INTERPOLATE_TEMP")
     CALL MCKPP_BOUNDARY_INTERPOLATE_TEMP()
     CALL mckpp_print(routine, "Returned from MCKPP_BOUNDARY_INTERPOLATE_TEMP")
  ENDIF

  IF (kpp_const_fields%L_RELAX_SAL .AND. .NOT. kpp_const_fields%L_INTERP_SAL) THEN
     CALL mckpp_print(routine, "Calling MCKPP_READ_SALINITY_3D")
     CALL MCKPP_READ_SALINITY_3D()
     CALL mckpp_print(routine, "Returned from MCKPP_READ_SALINITY_3D")
  ELSEIF (kpp_const_fields%L_RELAX_SAL .AND. kpp_const_fields%L_INTERP_SAL) THEN
     CALL MCKPP_BOUNDARY_INTERPOLATE_SAL()
  ENDIF

  CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_FLUXES_VARIABLES")
  CALL MCKPP_INITIALIZE_FLUXES_VARIABLES()
  CALL mckpp_print(routine, "Returned from MCKPP_INITIALIZE_FLUXES_VARIABLES")

  ! Isothermal detection routine requires 3D ocean temperature and salinity fields
  IF (kpp_const_fields%L_NO_ISOTHERM .AND. .NOT. kpp_const_fields%L_RELAX_SAL &
       .AND. .NOT. kpp_const_fields%L_RELAX_OCNT) THEN         
     CALL mckpp_print(routine, " Calling MCKPP_READ_TEMPERATURES_3D")
     CALL MCKPP_READ_TEMPERATURES_3D()
     CALL mckpp_print(routine, "Returned from MCKPP_READ_TEMPERATURES_3D")

     CALL mckpp_print(routine, "Calling MCKPP_READ_SALINITY_3D")
     CALL MCKPP_READ_SALINITY_3D()
     CALL mckpp_print(routine, "Returned from MCKPP_SALINITY_3D")
  ENDIF
  
  ! L_INTERP_OCNT and L_INTERP_SAL imply L_PERIODIC_OCNT and L_PERIODIC_SAL,
  ! respectively, to deal with times before the first time in the input file.
  IF (kpp_const_fields%L_INTERP_OCNT) kpp_const_fields%L_PERIODIC_OCNT=.TRUE.
  IF (kpp_const_fields%L_INTERP_SAL) kpp_const_fields%L_PERIODIC_SAL=.TRUE.

#ifndef MCKPP_COUPLE
  IF (kpp_const_fields%L_FLUXDATA) THEN
     CALL MCKPP_INITIALIZE_FLUXES_FILE()       
  ENDIF
#endif  

  kpp_const_fields%ntime=0
  ! CALL mckpp_print(routine, "Calling MCKPP_PHYSICS_LOOKUP")
  CALL MCKPP_PHYSICS_LOOKUP(kpp_const_fields)
  ! CALL mckpp_print(routine, "Returned from MCKPP_PHYSICS_LOOKUP")

  CALL mckpp_print(routine, "Calling MCKPP_INITIALIZE_OCEAN_MODEL")
  CALL MCKPP_INITIALIZE_OCEAN_MODEL()
  CALL mckpp_print(routine, "Returned from MCKPP_INITIALIZE_OCEAN_MODEL")

END SUBROUTINE MCKPP_INITIALIZE_FIELDS
