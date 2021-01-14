#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE MCKPP_INITIALIZE_FIELDS
  USE shr_kind_mod, only : r8=>shr_kind_r8
  USE mckpp_parameters
  USE mckpp_types, only : kpp_global_fields,kpp_3d_fields,kpp_const_fields
  USE phys_grid, only : get_ncols_p,get_rlat_all_p,get_rlon_all_p
  USE ppgrid, only : pcols,begchunk,endchunk
  USE pmgrid, only : masterproc
#else
SUBROUTINE MCKPP_INITIALIZE_FIELDS(kpp_3d_fields,kpp_const_fields)
  USE mckpp_data_types
#endif /*MCKPP_CAM3*/

  IMPLICIT NONE

#ifdef MCKPP_CAM3
  INTEGER :: icol,ncol,ichnk
  REAL(r8) :: clat1(pcols),clon1(pcols)
#else
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
#endif
  INTEGER :: iy,ix,ipt

#ifndef MCKPP_CAM3
  CALL mckpp_allocate_3d_fields(kpp_3d_fields)
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
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_INITIALIZE_LANDSEA'
     CALL MCKPP_INITIALIZE_LANDSEA
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_INITIALIZE_LANDSEA'
  ENDIF
#elif CFS
  IF (L_LANDSEA) CALL read_landsea_global
#else
  IF (kpp_const_fields%L_LANDSEA) THEN
     WRITE(6,*) "MCKPP_INITIALIZE_FIELDS: Calling MCKPP_INITIALIZE_LANDSEA"
     kpp_3d_fields%dlat(1)=kpp_const_fields%alat
     kpp_3d_fields%dlon(1)=kpp_const_fields%alon
     CALL MCKPP_INITIALIZE_LANDSEA(kpp_3d_fields,kpp_const_fields)
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
     WRITE(nuerr,*) 'KPP : If you set L_REGGRID=.FALSE., you must',&
          ' specify a land-sea mask file from which to read',&
          ' the locations of the gridpoints in the horizontal.'
  ENDIF
#endif

  ! Initialize the vertical grid
#ifdef MCKPP_CAM3
  IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_INITIALIZE_GEOGRAPHY'
  CALL MCKPP_INITIALIZE_GEOGRAPHY()
  IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_INITIALIZE_GEOGRAPHY'
#else
  CALL MCKPP_INITIALIZE_GEOGRAPHY(kpp_3d_fields,kpp_const_fields)
#endif  

  ! Initialize coupling weights   
#if (defined OASIS2 || defined OASIS3 || defined CFS)
  CALL MCKPP_INITIALIZE_COUPLINGWEIGHT(kpp_3d_fields,kpp_const_fields)
#elif (defined MCKPP_CAM3)
  IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_INITIALIZE_COUPLINGWEIGHT'
  CALL MCKPP_INITIALIZE_COUPLINGWEIGHT()
  IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_INITIALIZE_COUPLINGWEIGHT'
  ! Initialize coupling-period mean fluxes and SST fields
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%sflux_cpl(:,:)=0
     kpp_3d_fields(ichnk)%sst_cpl(:)=0
  ENDDO
#else
  IF (kpp_const_fields%L_CPLWGHT) CALL MCKPP_INITIALIZE_COUPLINGWEIGHT(kpp_3d_fields,kpp_const_fields)
#endif

  ! Initialize advection options
  IF (kpp_const_fields%L_ADVECT) THEN
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_INITIALIZE_ADVECTION'
     CALL MCKPP_INITIALIZE_ADVECTION()
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_INITIALIZE_ADVECTION'
  ELSE
     DO ichnk=begchunk,endchunk
        ncol=get_ncols_p(ichnk)
        kpp_3d_fields(ichnk)%nmodeadv(1:ncol,:)=0
     ENDDO
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: No advection has been specified.'
  ENDIF
#else
     CALL MCKPP_INITIALIZE_ADVECTION(kpp_3d_fields,kpp_const_fields)
  ELSE
     DO ipt=1,npts
        kpp_3d_fields%nmodeadv(ipt,1)=0
        kpp_3d_fields%nmodeadv(ipt,2)=0
     ENDDO
     write(6,*) 'KPP : No advection has been specified'
  ENDIF
#endif

  ! Initialize relaxation of SST, temperature and/or salinity
  IF (kpp_const_fields%L_RELAX_SST .OR. kpp_const_fields%L_RELAX_SAL &
       .OR. kpp_const_fields%L_RELAX_OCNT) THEN     
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_INITIALIZE_RELAXATION'
     CALL MCKPP_INITIALIZE_RELAXATION()
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_INITIALIZE_RELAXATION'
#else
     CALL MCKPP_INITIALIZE_RELAXATION(kpp_3d_fields,kpp_const_fields)
#endif     
  ENDIF

  ! Initialize water type for optical properties of seawater
#ifdef MCKPP_CAM3
  IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_INITIALIZE_OPTICS'
  CALL MCKPP_INITIALIZE_OPTICS
  IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_INITIALIZE_OPTICS' 
#else
  CALL MCKPP_INITIALIZE_OPTICS(kpp_3d_fields,kpp_const_fields)
#endif

  ! Initialize ocean profiles
  IF (kpp_const_fields%L_RESTART) THEN     
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_RESTART_IO_READ'
     ! Still needs scattering code
     CALL MCKPP_RESTART_IO_READ_NETCDF
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_RESTART_IO_READ'
#else
     CALL MCKPP_RESTART_IO_READ_NETCDF(kpp_3d_fields,kpp_const_fields)
#endif
  ELSE
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_INITIALIZE_OCEAN_PROFILES'
     CALL MCKPP_INITIALIZE_OCEAN_PROFILES
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_INITIALIZE_OCEAN_PROFILES'
#else
     CALL MCKPP_INITIALIZE_OCEAN_PROFILES(kpp_3d_fields,kpp_const_fields)
#endif  
  ENDIF
  
  ! Initialize boundary conditions
  IF (kpp_const_fields%L_CLIMSST) THEN     
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_READ_SST'
     CALL MCKPP_READ_SST()
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_READ_SST'
#else
     CALL MCKPP_READ_SST(kpp_3d_fields,kpp_const_fields)
#endif     
  ENDIF

  IF (kpp_const_fields%L_CLIMICE) THEN    
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_READ_ICE'
     CALL MCKPP_READ_ICE()
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_READ_ICE'
#else
     CALL MCKPP_READ_ICE(kpp_3d_fields,kpp_const_fields)
#endif    
  ENDIF
!!$    !IF (L_CLIMCURR) CALL read_surface_currents(kpp_3d_fields,kpp_const_fields)
  
  IF (kpp_const_fields%L_FCORR_WITHZ) THEN     
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_READ_FCORR_3D'
     CALL MCKPP_READ_FCORR_3D()
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_READ_FCORR_3D'
#else
     CALL MCKPP_READ_FCORR_3D(kpp_3d_fields,kpp_const_fields)
#endif     
  ELSEIF (kpp_const_fields%L_FCORR) THEN      
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_READ_FCORR_2D'
     CALL MCKPP_READ_FCORR_2D()
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_READ_FCORR_2D'
#else
     CALL MCKPP_READ_FCORR_2D(kpp_3d_fields,kpp_const_fields)
#endif  
  ENDIF

  IF (kpp_const_fields%L_SFCORR_WITHZ) THEN   
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_READ_SFCORR_3D'
     CALL MCKPP_READ_SFCORR_3D()
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_READ_SFCORR_3D'
#else
     CALL MCKPP_READ_SFCORR_3D(kpp_3d_fields,kpp_const_fields)
#endif   
  ELSEIF (kpp_const_fields%L_SFCORR) THEN
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_READ_SFCORR_2D'
     CALL MCKPP_READ_SFCORR_2D()
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_READ_SFCORR_2D'
#else
     CALL MCKPP_READ_SFCORR_2D(kpp_3d_fields,kpp_const_fields)
#endif
  ENDIF

  IF (kpp_const_fields%L_VARY_BOTTOM_TEMP) THEN     
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_READ_TEMPERATURES_BOTTOM'
     CALL MCKPP_READ_TEMPERATURES_BOTTOM()
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_READ_TEMPERATURES_BOTTOM'
#else
     CALL MCKPP_READ_TEMPERATURES_BOTTOM(kpp_3d_fields,kpp_const_fields)
#endif 
  ENDIF

  !!! THESE SHOULD NOT BE CALLED IF DOING L_INTERP_OCNT / L_INTERP_SAL !!!

  IF (kpp_const_fields%L_RELAX_OCNT .AND. .NOT. kpp_const_fields%L_INTERP_OCNT) THEN 
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_READ_TEMPERATURES_3D'
     CALL MCKPP_READ_TEMPERATURES_3D()
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_READ_TEMPERATURES_3D'
#else
     CALL MCKPP_READ_TEMPERATURES_3D(kpp_3d_fields,kpp_const_fields)
#endif     
  ELSEIF (kpp_const_fields%L_RELAX_OCNT .AND. kpp_const_fields%L_INTERP_OCNT) THEN     
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_BOUNDARY_INTERPOLATE_TEMP'
     CALL MCKPP_BOUNDARY_INTERPOLATE_TEMP
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_BOUNDARY_INTERPOLATE_TEMP'
#else
     CALL MCKPP_BOUNDARY_INTERPOLATE_TEMP(kpp_3d_fields,kpp_const_fields)
#endif   
  ENDIF

  IF (kpp_const_fields%L_RELAX_SAL .AND. .NOT. kpp_const_fields%L_INTERP_SAL) THEN
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_READ_SALINITY_3D' 
     CALL MCKPP_READ_SALINITY_3D()
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_READ_SALINITY_3D'  
#else
     CALL MCKPP_READ_SALINITY_3D(kpp_3d_fields,kpp_const_fields)
#endif   
  ELSEIF (kpp_const_fields%L_RELAX_SAL .AND. kpp_const_fields%L_INTERP_SAL) THEN
#ifdef MCKPP_CAM3     
     CALL MCKPP_BOUNDARY_INTERPOLATE_SAL
#else
     CALL MCKPP_BOUNDARY_INTERPOLATE_SAL(kpp_3d_fields,kpp_const_fields)
#endif
  ENDIF

#ifdef MCKPP_CAM3
  IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_INITIALIZE_FLUXES_VARIABLES'
  CALL MCKPP_INITIALIZE_FLUXES_VARIABLES
  IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_INITIALIZE_FLUXES_VARIABLES'
#else
  CALL MCKPP_INITIALIZE_FLUXES_VARIABLES(kpp_3d_fields)
#endif

  ! Isothermal detection routine requires 3D ocean temperature and salinity fields
  IF (kpp_const_fields%L_NO_ISOTHERM .AND. .NOT. kpp_const_fields%L_RELAX_SAL &
       .AND. .NOT. kpp_const_fields%L_RELAX_OCNT) THEN         
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_READ_TEMPERATURES_3D'
     CALL MCKPP_READ_TEMPERATURES_3D
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_READ_TEMPERATURES_3D'
#else
     CALL MCKPP_READ_TEMPERATURES_3D(kpp_3d_fields,kpp_const_fields)
#endif
#ifdef MCKPP_CAM3
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_READ_SALINITY_3D'
     CALL MCKPP_READ_SALINITY_3D
     IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_SALINITY_3D'
#else
     CALL MCKPP_READ_SALINITY_3D(kpp_3d_fields,kpp_const_fields)
#endif
  ENDIF
  
  ! L_INTERP_OCNT and L_INTERP_SAL imply L_PERIODIC_OCNT and L_PERIODIC_SAL,
  ! respectively, to deal with times before the first time in the input file.
  IF (kpp_const_fields%L_INTERP_OCNT) kpp_const_fields%L_PERIODIC_OCNT=.TRUE.
  IF (kpp_const_fields%L_INTERP_SAL) kpp_const_fields%L_PERIODIC_SAL=.TRUE.

#ifndef MCKPP_COUPLE
  IF (kpp_const_fields%L_FLUXDATA) THEN
     CALL MCKPP_INITIALIZE_FLUXES_FILE(kpp_const_fields)       
  ENDIF
#endif  

  kpp_const_fields%ntime=0
  !WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_PHYSICS_LOOKUP'
  CALL MCKPP_PHYSICS_LOOKUP(kpp_const_fields)
  !WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_PHYSICS_LOOKUP'

#ifdef MCKPP_CAM3
  IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Calling MCKPP_INITIALIZE_OCEAN_MODEL'
  CALL MCKPP_INITIALIZE_OCEAN_MODEL
  IF (masterproc) WRITE(6,*) 'MCKPP_INITIALIZE_FIELDS: Returned from MCKPP_INITIALIZE_OCEAN_MODEL'
#else
  CALL MCKPP_INITIALIZE_OCEAN_MODEL(kpp_3d_fields,kpp_const_fields)
#endif

  RETURN
END SUBROUTINE MCKPP_INITIALIZE_FIELDS
