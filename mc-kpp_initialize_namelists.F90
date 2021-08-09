#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
#endif

SUBROUTINE MCKPP_INITIALIZE_NAMELIST()
  
#ifdef MCKPP_CAM3  
  USE mckpp_types, only: kpp_const_fields
#else
  USE mckpp_data_fields, ONLY: kpp_const_fields, mckpp_allocate_const_fields
#endif
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, max_message_len
  USE mckpp_namelists
  USE mckpp_parameters 

  IMPLICIT NONE
  
  ! Local variables    
  INTEGER :: i,j,k,l,ipt,ix,iy
  CHARACTER(LEN=25) :: routine = "MCKPP_INITIALIZE_NAMELIST"
  CHARACTER(LEN=max_message_len) :: message
 
  ! This is a bug fix for the IBM XLF compiler, which otherwise complains
  ! about "incorrect characters" in the namelist.  If you are using the
  ! IBM compiler, you need to pass -WF,-DXLF_OLDNAME when you compile MC-KPP.
  ! NPK 5/6/09 - R2
#ifdef XLF_OLDNAME
  CALL SETRTEOPTS("namelist=old")
#endif
  
  ! Open the namelist
  OPEN(75,FILE='3D_ocn.nml')  

  ! Read parameters namelist first

  ! Set some defaults
  ndim = 1 
  nvel = 2
  nsclr = 2
  nsb = 1
  itermax = 200
  hmixtolfrac = 0.1
  nzl = 1 
  nzu = 2 
  nzdivmax = 8
  igridmax = 5
  nsflxs = 9
  njdt = 1
  ndharm = 5
  maxmodeadv = 6
  mr = 100

  ! These need to be defined in namelist 
  nx = 0.0 
  ny = 0.0 
  nz = 0.0 
  ngrid = 0.0
  nx_globe = 0.0
  ny_globe = 0.0
  
  READ(75,NAME_PARAMETERS) 
  CALL mckpp_print(routine, "Read Namelist PARAMETERS") 
  IF ( (nx .LE. 0) .OR. (ny .LE. 0) .OR. (nz .LE. 0) ) THEN 
    CALL mckpp_print_error(routine, "You must specify values of nx, ny and nz in the namelist")
    CALL MCKPP_ABORT()
  END IF 
  IF (ngrid .LE. 0) THEN 
    CALL mckpp_print_error(routine, "You must specify a value of ngrid in the namelist")
    CALL MCKPP_ABORT()
  END IF 
  IF (nztmax .LE. 0) THEN 
    CALL mckpp_print_error(routine, "You must specify a value of nztmax in the namelist")
    CALL MCKPP_ABORT()
  END IF 
  IF ( (nx_globe .LE. 0) .OR. (ny_globe .LE. 0) ) THEN 
    CALL mckpp_print_error(routine, "You must specify a value of nx_globe and ny_globe in the namelist")
    CALL MCKPP_ABORT()
  END IF 

  nzm1 = nz -1 
  nzp1 = nz+1
  npts = nx * ny 
  nvp1 = nvel + 1 
  nsp1 = nsclr + 1 
  nzp1tmax = nztmax + 1 
  nsflxsm1 = nsflxs - 1
  nsflxsp2 = nsflxs + 2
  mrp1 = mr + 1 
  npts_globe = nx_globe * ny_globe 

  WRITE(message,*) "nzm1, nzp1, npts, nvp1, nsp1, nzp1tmax, nsflxsm1, nsflxsp1, mrp1, npts_globe = "
  CALL mckpp_print(routine, message)
  WRITE(message,*) nzm1, nzp1, npts, nvp1, nsp1, nzp1tmax, nsflxsm1, nsflxsp2, mrp1, npts_globe
  CALL mckpp_print(routine, message)

#ifndef MCKPP_CAM3
  CALL mckpp_allocate_const_fields() 
#endif 
  allocate(kpp_const_fields%wmt(0:891,0:49))
  allocate(kpp_const_fields%wst(0:891,0:49))
  allocate(kpp_const_fields%tri(0:NZtmax,0:1,NGRID))

  ! Initialse and read the constants name list
  spd=86400.                ! secs/day
  dpy=360.                  ! days/year
  twopi=8*atan(1.)          ! 2pi
  onepi=twopi/2.            ! pi
  grav=9.816                ! gravity
  vonk=0.4                  ! Von Karman's constant
  TK0=273.15                ! Kelvin of 0degC
  sbc=5.67e-8               ! Stefan Boltzmann Constant
  epsw=1.0                  ! cor.fac for departure of H2O from B.body
  albocn=0.06               ! albedo for seawater
  sice=4.0                  ! salinity of ice(?)
  EL=2.50e6                 ! Latent heat of evap. at 0C (or constant)
  SL=2512200.               ! Latent heat of evap for ice
  FL=334000.                ! Latent heat of fusion for ice
  FLSN=FL                   ! Latent heat of fusion for snow
  READ(75,NAME_CONSTANTS)
  CALL mckpp_print(routine, "Read Namelist CONSTANTS")
  
  ! Initialize and read the processes namelist
  LKPP=.TRUE.
  LRI=.TRUE.
  LDD=.FALSE.
  LICE=.FALSE.
  LBIO=.FALSE.
  LTGRID=.FALSE.
  LNBFLX=.FALSE.
  LRHS=.FALSE.
  L_SSref=.TRUE.
  READ(75,NAME_PROCSWIT)
  CALL mckpp_print(routine, "Read Namelist PROCSWIT")
  
  ! Initilalize and read the location name list
  DMAX=0.0
  alat=0.0
  alon=0.0
  delta_lat=2.5
  delta_lon=3.75
  dscale=0.0
  L_STRETCHGRID=.FALSE.
  L_REGGRID=.TRUE.
  L_VGRID_FILE=.FALSE.
  READ(75,NAME_DOMAIN)
  IF (DMAX .LE. 0.0) THEN 
     CALL mckpp_print_error(routine, "You must specify a depth for the domain")
     CALL MCKPP_ABORT()
  ENDIF
  IF ((L_STRETCHGRID) .AND. (dscale .EQ. 0.0)) THEN
     CALL mckpp_print_error(routine, "You cannot have dscale=0 for stretched grids") 
     CALL MCKPP_ABORT()
  ENDIF
  CALL mckpp_print(routine, "Read Namelist DOMAIN")
  kpp_const_fields%alat=alat
  kpp_const_fields%alon=alon
  kpp_const_fields%delta_lat=delta_lat
  kpp_const_fields%delta_lon=delta_lon
  kpp_const_fields%dscale=dscale
  
  ! Initialize and read the landsea name list
  L_LANDSEA=.FALSE.
  READ(75,NAME_LANDSEA)
  CALL mckpp_print(routine, "Read Namelist LANDSEA")
  
  ! Initialize and read the start name list
  L_INITDATA= .TRUE.
  L_INTERPINIT= .TRUE.
  L_RESTART= .FALSE.
  WRITE(restart_infile,*) 'fort.30'
  READ(75,NAME_START) 
  CALL mckpp_print(routine, "Read Namelist START")
  
  ! Initialize and read the times namelist
  ndtocn=1
  dtsec=0.0
  startt=-999.999
  finalt=-999.999
  READ(75,NAME_TIMES) 
  IF ((dtsec .LE. 0.0) .OR. (startt .LT. 0.0) .OR. (finalt .LT. 0.0)) THEN 
     CALL mckpp_print_error(routine, "You must specify values of dtsec,startt,finalt in the namelist")
     CALL MCKPP_ABORT()
  ENDIF
  kpp_const_fields%ndtocn=ndtocn
  kpp_const_fields%spd=spd
  kpp_const_fields%dtsec=dtsec
  kpp_const_fields%startt=startt*kpp_const_fields%spd
  kpp_const_fields%finalt=finalt*kpp_const_fields%spd
  kpp_const_fields%dto=kpp_const_fields%dtsec/float(kpp_const_fields%ndtocn)
  kpp_const_fields%nend=int((kpp_const_fields%finalt-kpp_const_fields%startt)/kpp_const_fields%dtsec)
  kpp_const_fields%nstart=nint(kpp_const_fields%startt)/kpp_const_fields%dto
  kpp_const_fields%num_timesteps=kpp_const_fields%nend*kpp_const_fields%ndtocn
  IF (float(kpp_const_fields%num_timesteps) .NE. &
       (kpp_const_fields%finalt-kpp_const_fields%startt)/kpp_const_fields%dto) THEN
     CALL mckpp_print_error(routine, "The integration length is not a multiple of the ocean timestep")
     WRITE(message, *) "dto = ", kpp_const_fields%dto, ", finalt = ", kpp_const_fields%finalt, &
        ", startt = ", kpp_const_fields%startt
     CALL mckpp_print_error(routine, message) 
     CALL MCKPP_ABORT()
  ENDIF
  kpp_const_fields%startt=kpp_const_fields%startt/kpp_const_fields%spd
  kpp_const_fields%finalt=kpp_const_fields%finalt/kpp_const_fields%spd
  kpp_const_fields%time=kpp_const_fields%startt
  CALL mckpp_print(routine, "Read Namelist TIMES") 
  
  ! Initialize and read the couple namelist
#ifdef MCKPP_COUPLE
  L_COUPLE=.TRUE.
#else
  L_COUPLE=.FALSE.
#endif
  L_COUPLE_CURRENTS=.FALSE.
  L_OUTKELVIN=.FALSE.
  L_UPD_CLIMSST=.FALSE.
  L_UPD_CLIMICE=.FALSE.
  L_CLIMICE=.FALSE.
  L_CLIMSST=.FALSE.
  L_CLIMCURR=.FALSE. 
  L_BAD_ICE_DEPTH=.FALSE.
  ifirst=1
  ilast=nx
  jfirst=1
  jlast=ny
  READ(75,NAME_COUPLE)
  CALL mckpp_print(routine, "Read Namelist COUPLE")
  
  ! Initialize and read the advection namelist
  L_ADVECT=.FALSE.
  L_RELAX_SST=.FALSE.
  L_RELAX_CALCONLY=.FALSE.
#ifdef MCKPP_CAM3
  ALLOCATE( relax_sst_in(ny_globe) )
  ALLOCATE( relax_sal_in(ny_globe) )
  ALLOCATE( relax_ocnt_in(ny_globe) )
#else
  ALLOCATE( relax_sst_in(ny) )
  ALLOCATE( relax_sal_in(ny) )
  ALLOCATE( relax_ocnt_in(ny) )
#endif
  DO iy=1,ny
     relax_sst_in(iy)=0.0
     relax_sal_in(iy)=0.0
  ENDDO
  READ(75,NAME_ADVEC)
  
  ! Initialize and read the paras namelist
  paras_file='3D_ocnparas.nc'
  L_JERLOV=.TRUE.
  READ(75,NAME_PARAS)

  ! Initialize and read the forcing namelist
  L_FLUXDATA=.FALSE.
  L_FCORR_WITHZ=.FALSE.
  L_FCORR=.FALSE.
  L_UPD_FCORR=.FALSE.
  L_SFCORR_WITHZ=.FALSE.
  L_SFCORR=.FALSE.
  L_UPD_SFCORR=.FALSE.
  L_UPD_SAL=.FALSE.
  L_VARY_BOTTOM_TEMP=.FALSE.
  L_UPD_BOTTOM_TEMP=.FALSE.
  L_REST=.FALSE.
  L_NO_FREEZE=.FALSE.
  L_NO_ISOTHERM=.FALSE.
  L_DAMP_CURR=.FALSE.
  forcing_file='1D_ocean_forcing.nc'
  ocnT_file='none'
  READ(75,NAME_FORCING)
  CALL mckpp_print(routine, "Read Namelist FORCING")  
  IF (L_FCORR_WITHZ .AND. L_FCORR) THEN
     WRITE(message, *) "L_FCORR and L_FCORR_WITHZ are mutually exclusive. Choose one or neither."
     CALL mckpp_print_error(routine, message)
     CALL MCKPP_ABORT()
  ENDIF
  IF (L_SFCORR_WITHZ .AND. L_SFCORR) THEN
     WRITE(message, *) "L_SFCORR and L_SFCORR_WITHZ are mutually exclusive. Choose one or neither."
     CALL mckpp_print_error(routine, message)
     CALL MCKPP_ABORT()
  ENDIF
  IF (L_FCORR_WITHZ .AND. L_RELAX_SST) THEN
     WRITE(message, *) "L_FCORR_WITHZ and L_RELAX_SST are mutually exclusive. Choose one or neither."
     CALL mckpp_print_error(routine, message)
     CALL MCKPP_ABORT()
  ENDIF
  IF (L_NO_ISOTHERM .AND. (ocnT_file .eq. 'none' .or.&
      sal_file .eq. 'none')) THEN
     WRITE(message, *) "If you specify L_NO_ISOTHERM for reseting of isothermal points, " &
        // "you must specify files from which to read climatological ocean temperature " &
        // "(ocnT_file) and salinity (sal_file)."
     CALL mckpp_print_error(routine, message)
     CALL MCKPP_ABORT()
  ELSEIF (L_NO_ISOTHERM) THEN
     kpp_const_fields%iso_bot=isotherm_bottom
     kpp_const_fields%iso_thresh=isotherm_threshold
  ENDIF
  IF (L_DAMP_CURR) THEN
     kpp_const_fields%dt_uvdamp=dtuvdamp
  ENDIF
  
  L_RESTARTW=.TRUE.      
  kpp_const_fields%ndt_per_restart=kpp_const_fields%nend*kpp_const_fields%ndtocn    
    READ(75,NAME_OUTPUT)
  CALL mckpp_print(routine, "Read Namelist OUTPUT") 
  
  ! Call routine to copy constants and logicals needed for ocean
  ! physics into the kpp_const_fields derived type.  Added for 
  ! compatability with OpenMP DEFAULT(private). NPK 8/2/13
  CALL mckpp_initialize_constants(kpp_const_fields)    
  
END SUBROUTINE MCKPP_INITIALIZE_NAMELIST
