MODULE mckpp_initialize_constants_mod

CONTAINS

SUBROUTINE mckpp_initialize_constants(kpp_const_fields)

  ! This should be called after *ALL* constants have been read in
  ! in mckkp_initialize_namelists.  It should be called
  ! only once.

#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_const_type
#else 
  USE mckpp_data_fields, ONLY: kpp_const_type
#endif 
  USE mckpp_namelists
  USE mckpp_parameters
  
  IMPLICIT NONE

  TYPE(kpp_const_type),intent(inout) :: kpp_const_fields
  
  kpp_const_fields%spd=spd
  kpp_const_fields%dpy=dpy
  kpp_const_fields%twopi=twopi
  kpp_const_fields%onepi=onepi
  kpp_const_fields%grav=grav
  kpp_const_fields%vonk=vonk
  kpp_const_fields%TK0=TK0
  kpp_const_fields%sbc=sbc
  kpp_const_fields%epsw=epsw
  kpp_const_fields%albocn=albocn
  kpp_const_fields%sice=sice
  kpp_const_fields%EL=EL
  kpp_const_fields%SL=SL
  kpp_const_fields%FL=FL
  kpp_const_fields%FLSN=FLSN
  
  kpp_const_fields%LKPP=LKPP
  kpp_const_fields%LRI=LRI
  kpp_const_fields%LDD=LDD
  kpp_const_fields%LICE=LICE
  kpp_const_fields%LBIO=LBIO
  kpp_const_fields%LTGRID=LTGRID
  kpp_const_fields%LNBFLX=LNBFLX
  kpp_const_fields%LRHS=LRHS
  kpp_const_fields%L_SSref=L_SSref
  kpp_const_fields%L_RELAX_SST=L_RELAX_SST
  kpp_const_fields%L_RELAX_CALCONLY=L_RELAX_CALCONLY
  kpp_const_fields%L_FCORR=L_FCORR
  kpp_const_fields%L_SFCORR=L_SFCORR
  kpp_const_fields%L_FCORR_WITHZ=L_FCORR_WITHZ
  kpp_const_fields%L_SFCORR_WITHZ=L_SFCORR_WITHZ
  kpp_const_fields%L_RESTART=L_RESTART
  kpp_const_fields%L_RESTARTW=L_RESTARTW
  kpp_const_fields%L_RELAX_SAL=L_RELAX_SAL
  kpp_const_fields%L_RELAX_OCNT=L_RELAX_OCNT
  
  kpp_const_fields%ndt_per_restart=ndt_per_restart
  
  kpp_const_fields%sst_file=sstin_file
  kpp_const_fields%ndtupdsst=ndtupdsst
  kpp_const_fields%L_PERIODIC_CLIMSST=L_PERIODIC_CLIMSST
  kpp_const_fields%L_CLIMSST=L_CLIMSST
  kpp_const_fields%L_UPD_CLIMSST=L_UPD_CLIMSST
  kpp_const_fields%climsst_period=climsst_period
  
  kpp_const_fields%ice_file=icein_file
  kpp_const_fields%L_CLIMICE=L_CLIMICE
  kpp_const_fields%L_CLIM_ICE_DEPTH=L_CLIM_ICE_DEPTH
  kpp_const_fields%L_CLIM_SNOW_ON_ICE=L_CLIM_SNOW_ON_ICE
  kpp_const_fields%L_UPD_CLIMICE=L_UPD_CLIMICE
  kpp_const_fields%L_PERIODIC_CLIMICE=L_PERIODIC_CLIMICE
  kpp_const_fields%climice_period=climice_period
  kpp_const_fields%L_BAD_ICE_DEPTH=L_BAD_ICE_DEPTH
  kpp_const_fields%ndtupdice=ndtupdice

  kpp_const_fields%L_CLIMCURR=L_CLIMCURR
  kpp_const_fields%L_UPD_CLIMCURR=L_UPD_CLIMCURR
  kpp_const_fields%ndtupdcurr=ndtupdcurr
  
  kpp_const_fields%fcorr_file=fcorrin_file
  kpp_const_fields%L_UPD_FCORR=L_UPD_FCORR
  kpp_const_fields%ndtupdfcorr=ndtupdfcorr
  kpp_const_fields%fcorr_period=fcorr_period
  kpp_const_fields%L_PERIODIC_FCORR=L_PERIODIC_FCORR
  
  kpp_const_fields%sfcorr_file=sfcorrin_file
  kpp_const_fields%L_UPD_SFCORR=L_UPD_SFCORR
  kpp_const_fields%ndtupdsfcorr=ndtupdsfcorr
  kpp_const_fields%sfcorr_period=sfcorr_period
  kpp_const_fields%L_PERIODIC_SFCORR=L_PERIODIC_SFCORR
  
  kpp_const_fields%ocnT_file=ocnT_file
  kpp_const_fields%L_UPD_OCNT=L_UPD_OCNT      
  kpp_const_fields%ndtupdocnt=ndtupdocnt
  kpp_const_fields%L_PERIODIC_OCNT=L_PERIODIC_OCNT
  kpp_const_fields%L_INTERP_OCNT=L_INTERP_OCNT
  kpp_const_fields%ocnT_period=ocnT_period
  kpp_const_fields%ndt_interp_ocnt=ndt_interp_ocnt

  kpp_const_fields%sal_file=sal_file
  kpp_const_fields%L_UPD_SAL=L_UPD_SAL      
  kpp_const_fields%ndtupdsal=ndtupdsal
  kpp_const_fields%L_PERIODIC_SAL=L_PERIODIC_SAL
  kpp_const_fields%L_INTERP_SAL=L_INTERP_SAL
  kpp_const_fields%sal_period=sal_period
  kpp_const_fields%ndt_interp_sal=ndt_interp_sal
  
  kpp_const_fields%bottom_file=bottomin_file
  kpp_const_fields%L_VARY_BOTTOM_TEMP=L_VARY_BOTTOM_TEMP
  kpp_const_fields%L_UPD_BOTTOM_TEMP=L_UPD_BOTTOM_TEMP
  kpp_const_fields%ndtupdbottom=ndtupdbottom
  kpp_const_fields%L_PERIODIC_BOTTOM_TEMP=L_PERIODIC_BOTTOM_TEMP
  kpp_const_fields%bottom_temp_period=bottom_temp_period
  
  kpp_const_fields%L_OUTKELVIN=L_OUTKELVIN
  kpp_const_fields%ifirst=ifirst
  kpp_const_fields%jfirst=jfirst
  kpp_const_fields%ilast=ilast
  kpp_const_fields%jlast=jlast
  kpp_const_fields%L_COUPLE_CURRENTS=L_COUPLE_CURRENTS

  kpp_const_fields%L_FLUXDATA=L_FLUXDATA
  kpp_const_fields%forcing_file=forcing_file
  kpp_const_fields%L_REST=L_REST
  kpp_const_fields%L_REGGRID=L_REGGRID
  kpp_const_fields%L_LANDSEA=L_LANDSEA
  kpp_const_fields%landsea_file=landsea_file
  kpp_const_fields%L_VGRID_FILE=L_VGRID_FILE
  kpp_const_fields%vgrid_file=vgrid_file
  kpp_const_fields%L_STRETCHGRID=L_STRETCHGRID
  kpp_const_fields%DMAX=DMAX
  kpp_const_fields%L_CPLWGHT=L_CPLWGHT
  kpp_const_fields%cplwght_file=cplwght_file
  kpp_const_fields%L_ADVECT=L_ADVECT
  kpp_const_fields%advect_file=advect_file
  kpp_const_fields%L_JERLOV=L_JERLOV
  kpp_const_fields%paras_file=paras_file

  kpp_const_fields%relax_sst_in=relax_sst_in
  kpp_const_fields%relax_sal_in=relax_sal_in
  kpp_const_fields%relax_ocnT_in=relax_ocnT_in

  kpp_const_fields%restart_infile=restart_infile
  kpp_const_fields%restart_outfile=restart_outfile

  kpp_const_fields%L_INITDATA=L_INITDATA
  kpp_const_fields%initdata_file=initdata_file
  kpp_const_fields%L_INTERPINIT=L_INTERPINIT

  kpp_const_fields%L_NO_ISOTHERM=L_NO_ISOTHERM
  kpp_const_fields%L_NO_FREEZE=L_NO_FREEZE
  kpp_const_fields%L_DAMP_CURR=L_DAMP_CURR

END SUBROUTINE mckpp_initialize_constants

END MODULE mckpp_initialize_constants_mod
