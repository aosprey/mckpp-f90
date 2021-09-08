MODULE mckpp_namelists

  USE mckpp_netcdf_read, ONLY: max_nc_filename_len, max_restart_filename_len
  USE mckpp_parameters

  IMPLICIT NONE 

  ! Variables read in from namelists. These used to be held in common blocks. 
  ! This module is just used for the initialization, as the useful   
  ! variables are copied into kpp_const_type derived type.

  ! name_advec
  LOGICAL :: L_ADVECT, L_RELAX_SST, L_RELAX_CALCONLY, L_RELAX_SAL, L_RELAX_OCNT
  CHARACTER(max_restart_filename_len) :: advect_file
  REAL, ALLOCATABLE, DIMENSION(:) :: relax_sst_in, relax_sal_in, relax_ocnT_in
  INTEGER :: ncid_advec 
  NAMELIST/name_advec/ L_ADVECT, advect_file, L_RELAX_SST, &
      relax_sst_in ,relax_sal_in, L_RELAX_CALCONLY, L_RELAX_SAL, &
      L_RELAX_OCNT, relax_ocnt_in

  ! name_constants
  REAL :: grav, vonk, sbc, twopi, onepi, TK0, spd, dpy, epsw, albocn, sice, EL, SL, FL, &
      FLSN
  NAMELIST/name_constants/ grav, vonk, sbc, twopi, onepi, TK0, spd, dpy, &
      epsw, albocn, sice, EL, SL, FL, FLSN

  ! name_coupling
  LOGICAL :: L_COUPLE, L_CLIMSST, L_UPD_CLIMSST, L_CPLWGHT, L_CLIMICE, L_UPD_CLIMICE, &
      L_CLIM_ICE_DEPTH, L_CLIM_SNOW_ON_ICE, L_OUTKELVIN, L_COUPLE_CURRENTS, L_CLIMCURR, &
      L_UPD_CLIMCURR, L_PERIODIC_CLIMICE, L_PERIODIC_CLIMSST, L_BAD_ICE_DEPTH      
  INTEGER :: ifirst, ilast, jfirst, jlast, ndtupdsst, ndtupdice, ndtupdcurr, &
      climsst_period, climice_period
  CHARACTER(LEN=max_nc_filename_len) :: sstin_file, cplwght_file, icein_file, currin_file
  NAMELIST/name_couple/ L_COUPLE, ifirst, ilast, jfirst, jlast, L_CLIMSST, sstin_file, &
      L_UPD_CLIMSST, ndtupdsst, L_CPLWGHT, cplwght_file, icein_file, L_CLIMICE, &
      L_UPD_CLIMICE, ndtupdice, L_CLIM_ICE_DEPTH, L_CLIM_SNOW_ON_ICE, L_OUTKELVIN, &
      L_COUPLE_CURRENTS, currin_file, L_CLIMCURR, L_UPD_CLIMCURR, ndtupdcurr, & 
      L_PERIODIC_CLIMICE, L_PERIODIC_CLIMSST, climsst_period, climice_period, &
      L_BAD_ICE_DEPTH

  ! name_domain
  LOGICAL :: L_STRETCHGRID, L_REGGRID, L_VGRID_FILE
  REAL :: dmax, alon, alat, delta_lat, delta_lon, dscale
  CHARACTER(LEN=max_nc_filename_len) :: vgrid_file
  NAMELIST/name_domain/ DMAX, alon, alat, delta_lat, delta_lon, L_STRETCHGRID, dscale, &
      L_REGGRID, L_VGRID_FILE, vgrid_file

  ! name_forcing
  LOGICAL :: L_FLUXDATA, L_FCORR_WITHZ, L_VARY_BOTTOM_TEMP, L_FCORR, L_UPD_FCORR, &
      L_UPD_BOTTOM_TEMP, L_REST, L_PERIODIC_FCORR, L_PERIODIC_BOTTOM_TEMP, &
      L_SFCORR_WITHZ, L_SFCORR, L_UPD_SFCORR, L_PERIODIC_SFCORR, L_UPD_SAL, &
      L_PERIODIC_SAL, L_UPD_OCNT, L_PERIODIC_OCNT, L_NO_FREEZE, L_NO_ISOTHERM, &
      L_DAMP_CURR, L_INTERP_OCNT, L_INTERP_SAL
  INTEGER :: ndtupdfcorr, ndtupdbottom, fcorr_period, ndtupdsfcorr, sfcorr_period, &
      bottom_temp_period, sal_period, ndtupdsal, ocnt_period, ndtupdocnt, &
      isotherm_bottom, dtuvdamp, ndt_interp_ocnt, ndt_interp_sal
  REAL :: isotherm_threshold
  CHARACTER (LEN=max_nc_filename_len) :: forcing_file, fcorrin_file, bottomin_file, &
      sfcorrin_file, sal_file, ocnt_file
  NAMELIST/name_forcing/ L_FLUXDATA, forcing_file, L_FCORR_WITHZ, fcorrin_file, &
      ndtupdfcorr, L_VARY_BOTTOM_TEMP, ndtupdbottom,  bottomin_file, L_FCORR, &
      L_UPD_FCORR, L_UPD_BOTTOM_TEMP, L_REST, L_PERIODIC_FCORR, L_PERIODIC_BOTTOM_TEMP, &
      fcorr_period, L_SFCORR_WITHZ, sfcorrin_file, ndtupdsfcorr, L_SFCORR, L_UPD_SFCORR, &
      L_PERIODIC_SFCORR, sfcorr_period, bottom_temp_period, sal_file, L_UPD_SAL, &
      L_PERIODIC_SAL, sal_period, ndtupdsal, ocnt_file, L_UPD_OCNT, L_PERIODIC_OCNT, &
      ocnt_period, ndtupdocnt, L_NO_FREEZE, L_NO_ISOTHERM, isotherm_bottom, &
      isotherm_threshold, L_DAMP_CURR, dtuvdamp, L_INTERP_OCNT, ndt_interp_ocnt, &
      L_INTERP_SAL, ndt_interp_sal

  ! name_landsea
  LOGICAL :: L_LANDSEA
  CHARACTER(LEN=max_nc_filename_len) :: landsea_file
  NAMELIST/name_landsea/ L_LANDSEA, landsea_file

  ! name_output
  LOGICAL :: L_RESTARTW
  INTEGER :: ndt_per_restart
  CHARACTER(max_restart_filename_len) ::  restart_outfile
  NAMELIST/name_output/ L_RESTARTW, restart_outfile, ndt_per_restart
  
  ! name_parameters
  ! These variables are defined in mckpp_parameters as they are used throughout the code. 
  NAMELIST/name_parameters/ nz, ndim, nx, ny, nvel, nsclr, nsb, itermax, hmixtolfrac, & 
      ngrid, nzl, nzu, nzdivmax, nztmax, igridmax, nsflxs, njdt, ndharm, maxmodeadv, &
      mr, nx_globe, ny_globe

  ! name_paras
  CHARACTER(LEN=max_nc_filename_len) :: paras_file
  LOGICAL :: L_JERLOV
  NAMELIST/name_paras/ paras_file, L_JERLOV

  ! name_procswit
  LOGICAL LKPP , LRI, LDD, LICE, LBIO, LNBFLX, LTGRID, LRHS, L_SSref
  NAMELIST/name_procswit/ LKPP, LRI, LDD, LICE, LBIO, LNBFLX, LTGRID, LRHS, L_SSref

  ! name_start
  LOGICAL :: L_INITDATA, L_INTERPINIT, L_RESTART
  CHARACTER(max_nc_filename_len) :: initdata_file
  CHARACTER(max_restart_filename_len) :: restart_infile
  NAMELIST/name_start/ L_INITDATA, initdata_file, L_INTERPINIT, L_RESTART, restart_infile

  ! name_times
  REAL :: dtsec, startt, finalt
  INTEGER :: ndtocn, nyear
  NAMELIST/name_times/ dtsec, startt, finalt, ndtocn, nyear

END MODULE mckpp_namelists
