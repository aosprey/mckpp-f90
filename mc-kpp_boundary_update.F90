SUBROUTINE MCKPP_BOUNDARY_UPDATE()

#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_const_fields
#else
  USE mckpp_data_fields, ONLY: kpp_const_fields
#endif
  USE mckpp_parameters, ONLY: nuout

  ! Update all boundary conditions that are read from netCDF files,
  ! except for surface fluxes, which are handled in MCKPP_UPDATE_FLUXES or coupling routines.

  IMPLICIT NONE

  ! Update SST
  IF (kpp_const_fields%L_UPD_CLIMSST .AND. &
       MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndtupdsst) .EQ. 0) THEN
     !WRITE(6,*) 'MCKPP_BOUNDARY_UPDATE: Calling MCKPP_READ_SST'
     CALL MCKPP_READ_SST()
     !WRITE(6,*) 'MCKPP_BOUNDARY_UPDATE: Returned from MCKPP_READ_SST'
     !WRITE(nuout,*) 'KPP: Called read_sstin, ntime =',kpp_const_fields%ntime

     ! SST0 specifies the temperature to which to relax the mixed-layer temperature.
     ! The use of the relaxation is controlled by L_RELAX_SST.
     ! Changing L_RELAX_SST .OR. L_COUPLE to just L_RELAX_SST - NPK 2/11/09 - R3
     IF (kpp_const_fields%L_RELAX_SST) THEN
        !WRITE(6,*) 'MCKPP_BOUNDARY_UPDATE: Calling MCKPP_PHYSICS_OVERRIDES_SST0'
        CALL MCKPP_PHYSICS_OVERRIDES_SST0()
        !WRITE(nuout,*) 'KPP: Called upd_sst0, ntime =',kpp_const_fields%ntime               
     ENDIF
  ENDIF
  
  ! Update sea ice  
  IF (kpp_const_fields%L_UPD_CLIMICE .AND. MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndtupdice) .EQ. 0) THEN
     !WRITE(6,*) 'MCKPP_BOUNDARY_UPDATE: Calling MCKPP_READ_ICE'
     CALL MCKPP_READ_ICE()
     !WRITE(6,*) 'MCKPP_BOUNDARY_UPDATE: Returned from MCKPP_READ_ICE'
     !WRITE(nuout,*) 'KPP: Called read_icein, ntime =',kpp_const_fields%ntime
  ENDIF
  
  ! Update surface currents
  !IF(kpp_const_fields%L_UPD_CLIMCURR .AND. MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndtupdcurr) .EQ. 0) THEN
     !CALL read_surface_currents(kpp_3d_fields,kpp_const_fields)
     !WRITE(nuout,*) 'KPP: Called read_surface_currents, ntime =',kpp_const_fields%ntime
  !ENDIF
  
  ! Update heat corrections
  IF (kpp_const_fields%L_UPD_FCORR .AND. MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndtupdfcorr) .EQ. 0) THEN
     IF (kpp_const_fields%L_FCORR_WITHZ) THEN
        CALL MCKPP_READ_FCORR_3D()
        !WRITE(nuout,*) 'KPP: Called read_fcorrwithz, ntime =',kpp_const_fields%ntime
     ELSEIF (kpp_const_fields%L_FCORR) THEN
        CALL MCKPP_READ_FCORR_2D()
        !WRITE(nuout,*) 'KPP: Called read_fcorr, ntime =',kpp_const_fields%ntime
     ENDIF
  ENDIF

  ! Update salt corrections
  IF (kpp_const_fields%L_UPD_SFCORR .AND. MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndtupdsfcorr) .EQ. 0) THEN
     IF (kpp_const_fields%L_SFCORR_WITHZ) THEN
        CALL MCKPP_READ_SFCORR_3D()
        !WRITE(nuout,*) 'KPP: Called read_sfcorrwithz, ntime =',kpp_const_fields%ntime
     ELSEIF (kpp_const_fields%L_SFCORR) THEN
        CALL MCKPP_READ_SFCORR_2D()
        !WRITE(nuout,*) 'KPP: Called read_sfcorr, ntime =',kpp_const_fields%ntime
     ENDIF
  ENDIF
  
  ! Update bottom temperatures
  IF (kpp_const_fields%L_UPD_BOTTOM_TEMP .AND. MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndtupdbottom) .EQ. 0) THEN            
     CALL MCKPP_READ_TEMPERATURES_BOTTOM()
     !WRITE(nuout,*) 'KPP: Called read_bottom_temp, ntime =',kpp_const_fields%ntime
  ENDIF

  ! Update reference salinity profiles
  IF (kpp_const_fields%L_UPD_SAL .AND. MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndtupdsal) .EQ. 0 .AND. &
       .NOT. kpp_const_fields%L_INTERP_SAL) THEN
     CALL MCKPP_READ_SALINITY_3D()
     !WRITE(nuout,*) 'KPP: Called read_salinity, ntime =',kpp_const_fields%ntime
  ELSEIF (kpp_const_fields%L_UPD_SAL .AND. kpp_const_fields%L_INTERP_SAL .AND. &
       MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndt_interp_sal).EQ.0) THEN
     CALL MCKPP_BOUNDARY_INTERPOLATE_SAL()
     !WRITE(nuout,*) 'KPP: Interpolated ocean salinity, ntime =',kpp_const_fields%ntime
  ENDIF

  ! Update reference temperature profiles
  IF (kpp_const_fields%L_UPD_OCNT .AND. MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndtupdocnt) .EQ. 0 .AND. &
       .NOT. kpp_const_fields%L_INTERP_OCNT) THEN
     CALL MCKPP_READ_TEMPERATURES_3D()
     !WRITE(nuout,*) 'KPP: Called read_ocean_temperatures, ntime =',kpp_const_fields%ntime
  ELSEIF (kpp_const_fields%L_UPD_OCNT .AND. kpp_const_fields%L_INTERP_OCNT .AND. &
       MOD(kpp_const_fields%ntime-1,kpp_const_fields%ndt_interp_ocnt).EQ.0) THEN
     CALL MCKPP_BOUNDARY_INTERPOLATE_TEMP()
     !WRITE(nuout,*) 'KPP: Interpolated ocean temperatures, ntime =',kpp_const_fields%ntime
  ENDIF

END SUBROUTINE MCKPP_BOUNDARY_UPDATE
