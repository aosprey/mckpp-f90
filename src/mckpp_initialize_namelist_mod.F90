MODULE mckpp_initialize_namelist_mod

  USE mckpp_data_fields, ONLY: kpp_const_fields, mckpp_allocate_const_fields
  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_initialize_constants_mod, ONLY: mckpp_initialize_constants
  USE mckpp_log_messages, ONLY: mckpp_print, mckpp_print_error, & 
        max_message_len, update_context, l_debug
  USE mckpp_mpi_control, ONLY: root, l_root, comm
  USE mckpp_namelists
  USE mckpp_parameters
  USE mpi

  IMPLICIT NONE

  PUBLIC :: mckpp_initialize_namelist

  PRIVATE

CONTAINS

  SUBROUTINE mckpp_initialize_namelist()

    INTEGER :: nml_unit = 37
    CHARACTER(LEN=10) :: nml_file = "3D_ocn.nml"

    ! Only single pe reads namelist data, then sends data to others
    IF (l_root) THEN
      OPEN(nml_unit, file=nml_file)
    END IF

    CALL read_parameters_namelist(nml_unit)
    CALL read_constants_namelist(nml_unit)
    CALL read_procswit_namelist(nml_unit)
    CALL read_domain_namelist(nml_unit)
    CALL read_landsea_namelist(nml_unit)
    CALL read_start_namelist(nml_unit)
    CALL read_times_namelist(nml_unit)
    CALL read_couple_namelist(nml_unit)
    CALL read_advec_namelist(nml_unit)   
    CALL read_paras_namelist(nml_unit)
    CALL read_forcing_namelist(nml_unit)
    CALL read_output_namelist(nml_unit)

    IF (l_root) THEN 
      CLOSE(nml_unit)
    END IF

    ! Set up const_fields based on namelist data
    CALL mckpp_allocate_const_fields() 
    CALL mckpp_initialize_constants(kpp_const_fields)

  END SUBROUTINE mckpp_initialize_namelist


  SUBROUTINE read_parameters_namelist(nml_unit)

    INTEGER, INTENT(IN) :: nml_unit
    CHARACTER(LEN=24) :: routine = "READ_PARAMETERS_NAMELIST"
    CHARACTER(LEN=max_message_len) :: message
 
    ! Structure to hold parameters namelist values
    TYPE name_type
      SEQUENCE
      REAL(KIND=8) :: hmixtolfrac
      INTEGER :: nz, ndim, nx, ny, nvel, nsclr, nsb, itermax, &
          ngrid, nzl, nzu, nzdivmax, nztmax, nzp1tmax, igridmax, &
          nsflxs, ndharm, maxmodeadv, mr, nx_globe, ny_globe
    END TYPE name_type
    TYPE(name_type) my_name

    INTEGER :: count, name_type_mpi, ierr
    INTEGER, DIMENSION(2) :: types, lens

    ! Define MPI type to hold namelist data
    count = 2
    types(1) = MPI_REAL8
    types(2) = MPI_INT
    lens(1) = 1
    lens(2) = 21
    CALL make_mpi_type(routine, count, types, lens, name_type_mpi)

    IF (l_root) THEN
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

      READ(nml_unit,NAME_PARAMETERS) 
      CALL mckpp_print(routine, "Read Namelist PARAMETERS")
      
      IF ( (nx .LE. 0) .OR. (ny .LE. 0) .OR. (nz .LE. 0) ) THEN
        CALL mckpp_abort(routine, &
            "You must specify values of nx, ny and nz in the namelist")
      END IF
      IF (ngrid .LE. 0) THEN 
        CALL mckpp_abort(routine, &
            "You must specify a value of ngrid in the namelist")
      END IF
      IF (nztmax .LE. 0) THEN 
        CALL mckpp_abort(routine, &
            "You must specify a value of nztmax in the namelist")
      END IF
      IF ( (nx_globe .LE. 0) .OR. (ny_globe .LE. 0) ) THEN 
        CALL mckpp_abort(routine, &
            "You must specify a value of nx_globe and ny_globe in the namelist")
      END IF

      ! Pack into derived type
      my_name % hmixtolfrac = hmixtolfrac
      my_name % nz = nz
      my_name % ndim = ndim
      my_name % nx = nx
      my_name % ny = ny
      my_name % nvel = nvel
      my_name % nsclr = nsclr
      my_name % nsb = nsb
      my_name % itermax = itermax
      my_name % ngrid = ngrid
      my_name % nzl = nzl
      my_name % nzu = nzu
      my_name % nzdivmax = nzdivmax
      my_name % nztmax = nztmax
      my_name % nzp1tmax = nzp1tmax
      my_name % igridmax = igridmax
      my_name % nsflxs = nsflxs
      my_name % ndharm = ndharm
      my_name % maxmodeadv = maxmodeadv
      my_name % mr = mr
      my_name % nx_globe = nx_globe
      my_name % ny_globe = ny_globe     
    ENDIF

    ! Send data
    CALL mpi_bcast(my_name, 1, name_type_mpi, root, comm, ierr)

    ! Free memory from MPI type 
    CALL mpi_type_free(name_type_mpi, ierr)

    ! Unpack data
    IF (.NOT. l_root) THEN
      hmixtolfrac = my_name % hmixtolfrac
      nz = my_name % nz
      ndim = my_name % ndim
      nx = my_name % nx
      ny = my_name % ny
      nvel = my_name % nvel
      nsclr = my_name % nsclr
      nsb = my_name % nsb
      itermax = my_name % itermax
      ngrid = my_name % ngrid
      nzl = my_name % nzl
      nzu = my_name % nzu
      nzdivmax = my_name % nzdivmax
      nztmax = my_name % nztmax
      nzp1tmax = my_name % nzp1tmax
      igridmax = my_name % igridmax
      nsflxs = my_name % nsflxs
      ndharm = my_name % ndharm
      maxmodeadv = my_name % maxmodeadv
      mr = my_name % mr
      nx_globe = my_name % nx_globe
      ny_globe = my_name % ny_globe     
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

    IF (l_debug) THEN 
      WRITE(message,*) "hmixtolfrac, nz, ndim, nx, ny, nvel, nsclr, nsb = ", & 
                        hmixtolfrac, nz, ndim, nx, ny, nvel, nsclr, nsb 
      CALL mckpp_print(routine, message)
      WRITE(message,*) "itermax, ngrid, nzl, nzu, nzdivmax, nztmax = ", &  
                        itermax, ngrid, nzl, nzu, nzdivmax, nztmax 
      CALL mckpp_print(routine, message)
      WRITE(message,*) "nzp1tmax, igridmax, nsflxs, ndharm, maxmodeadv, mr = ", &
                        nzp1tmax, igridmax, nsflxs, ndharm, maxmodeadv, mr
      CALL mckpp_print(routine, message)
      WRITE(message,*) "nx_globe, ny_globe, nzm1, nzp1, npts, nvp1, nsp1 = ", & 
                        nx_globe, ny_globe, nzm1, nzp1, npts, nvp1, nsp1
      CALL mckpp_print(routine, message)
      WRITE(message,*) "nzp1tmax, nsflxsm1, nsflxsp2, mrp1, npts_globe = ", &
                        nzp1tmax, nsflxsm1, nsflxsp2, mrp1, npts_globe
      CALL mckpp_print(routine, message)
   END IF 
      
  END SUBROUTINE read_parameters_namelist


  SUBROUTINE read_constants_namelist(nml_unit)

    INTEGER, INTENT(IN) :: nml_unit
    CHARACTER(LEN=23) :: routine = "READ_CONSTANTS_NAMELIST"
    CHARACTER(LEN=max_message_len) :: message
 
    TYPE name_type
      SEQUENCE
      REAL(KIND=8) :: grav, vonk, sbc, twopi, onepi, TK0, spd, dpy, epsw, & 
                      albocn, EL, SL, FL, FLSN
    END TYPE name_type
    TYPE(name_type) my_name

    INTEGER :: count, name_type_mpi, ierr
    INTEGER, DIMENSION(1) :: types, lens

    ! Define MPI type
    count = 1
    types(1) = MPI_REAL8
    lens(1) = 14
    CALL make_mpi_type(routine, count, types, lens, name_type_mpi)

    ! Initialize, read nml and pack into buffer
    IF (l_root) THEN 
      spd = 86400.                ! secs/day
      dpy = 360.                  ! days/year
      twopi = 8*atan(1.)          ! 2pi
      onepi = twopi/2.            ! pi
      grav = 9.816                ! gravity
      vonk = 0.4                  ! Von Karman's constant
      TK0 = 273.15                ! Kelvin of 0degC
      sbc = 5.67e-8               ! Stefan Boltzmann Constant
      epsw = 1.0                  ! cor.fac for departure of H2O from B.body
      albocn = 0.06               ! albedo for seawater
      sice = 4.0                  ! salinity of ice(?)
      EL = 2.50e6                 ! Latent heat of evap. at 0C (or constant)
      SL = 2512200.               ! Latent heat of evap for ice
      FL = 334000.                ! Latent heat of fusion for ice
      FLSN = FL                   ! Latent heat of fusion for snow

      READ(nml_unit,NAME_CONSTANTS)
      CALL mckpp_print(routine, "Read Namelist CONSTANTS")

      my_name % grav = grav
      my_name % vonk = vonk
      my_name % sbc = sbc
      my_name % twopi = twopi
      my_name % onepi = onepi
      my_name % TK0 = TK0
      my_name % spd = spd
      my_name % dpy = dpy
      my_name % epsw = epsw
      my_name % albocn = albocn
      my_name % EL = EL
      my_name % SL = SL
      my_name % FL = FL
      my_name % FLSN = FLSN
    END IF

    ! Send data
    CALL mpi_bcast(my_name, 1, name_type_mpi, root, comm, ierr)  
    CALL mpi_type_free(name_type_mpi, ierr)

    ! Unpack data
    IF (.NOT. l_root) THEN
      grav = my_name % grav
      vonk = my_name % vonk
      sbc = my_name % sbc
      twopi = my_name % twopi
      onepi = my_name % onepi
      TK0 = my_name % TK0
      spd = my_name % spd
      dpy = my_name % dpy
      epsw = my_name % epsw
      albocn = my_name % albocn
      EL = my_name % EL
      SL = my_name % SL
      FL = my_name % FL
      FLSN = my_name % FLSN
    END IF

    IF (l_debug) THEN 
      WRITE(message, *) "grav, vonk, sbc, twopi, onepi = ",  & 
                         grav, vonk, sbc, twopi, onepi
      CALL mckpp_print(routine, message)
      WRITE(message, *) "TK0, spd, dpy, epsw, albocn = ", & 
                         TK0, spd, dpy, epsw, albocn
      CALL mckpp_print(routine, message)
      WRITE(message, *) "EL, SL, FL, FLSN = ", & 
                         EL, SL, FL, FLSN
      CALL mckpp_print(routine, message)
    END IF 
      
  END SUBROUTINE read_constants_namelist


  SUBROUTINE read_procswit_namelist(nml_unit)

    INTEGER, INTENT(IN) :: nml_unit 
    CHARACTER(LEN=22) :: routine = "READ_PROCSWIT_NAMELIST"
    CHARACTER(LEN=max_message_len) :: message

    TYPE name_type
      SEQUENCE
      LOGICAL :: LKPP, LRI, LDD, LICE, LBIO, LNBFLX, LTGRID, LRHS, L_SSref
    END TYPE name_type
    TYPE(name_type) my_name

    INTEGER :: count, name_type_mpi, ierr
    INTEGER, DIMENSION(1) :: types, lens

    ! Define MPI type
    count = 1
    types(1) = MPI_LOGICAL
    lens(1) = 14
    CALL make_mpi_type(routine, count, types, lens, name_type_mpi)

    ! Initialize, read nml and pack into buffer
    IF (l_root) THEN 
      LKPP = .TRUE.
      LRI = .TRUE.
      LDD = .FALSE.
      LICE = .FALSE.
      LBIO = .FALSE.
      LTGRID = .FALSE.
      LNBFLX = .FALSE.
      LRHS = .FALSE.
      L_SSref = .TRUE.

      READ(nml_unit,NAME_PROCSWIT)
      CALL mckpp_print(routine, "Read Namelist PROCSWIT")

      my_name % LKPP = LKPP
      my_name % LRI = LRI
      my_name % LDD = LDD
      my_name % LICE = LICE
      my_name % LBIO = LBIO
      my_name % LNBFLX = LNBFLX
      my_name % LTGRID = LTGRID
      my_name % LRHS = LRHS
      my_name % L_SSref = L_SSref
    END IF

    ! Send data
    CALL mpi_bcast(my_name, 1, name_type_mpi, root, comm, ierr)
    CALL mpi_type_free(name_type_mpi, ierr)

    ! Unpack data
    IF (.NOT. l_root) THEN
      LKPP = my_name % LKPP
      LRI = my_name % LRI
      LDD = my_name % LDD
      LICE = my_name % LICE
      LBIO = my_name % LBIO
      LNBFLX = my_name % LNBFLX
      LTGRID = my_name % LTGRID
      LRHS = my_name % LRHS
      L_SSref = my_name % L_SSref
    END IF

    IF (l_debug) THEN 
      WRITE(message, *) "LKPP, LRI,LDD ,LICE, LBIO, LNBFLX, LTGRID = ", & 
                         LKPP, LRI,LDD ,LICE, LBIO, LNBFLX, LTGRID
      CALL mckpp_print(routine, message)
      WRITE(message, *) "LRHS, L_SSref = ", LRHS, L_SSref
      CALL mckpp_print(routine, message)
    END IF 
      
  END SUBROUTINE read_procswit_namelist


  SUBROUTINE read_domain_namelist(nml_unit)

    INTEGER, INTENT(IN) :: nml_unit 
    CHARACTER(LEN=20) :: routine = "READ_DOMAIN_NAMELIST"
    CHARACTER(LEN=max_message_len) :: message
 
    TYPE name_type
      SEQUENCE
      LOGICAL :: L_REGGRID, L_STRETCHGRID, L_VGRID_FILE
      REAL :: dmax, alat, alon, delta_lat, delta_lon, dscale
      CHARACTER(LEN=max_nc_filename_len) :: vgrid_file = ""
    END TYPE name_type
    TYPE(name_type) my_name

    INTEGER :: count, name_type_mpi, ierr
    INTEGER, DIMENSION(3) :: types, lens

    ! Define MPI type
    count = 3
    types(1) = MPI_LOGICAL
    types(2) = MPI_REAL8
    types(3) = MPI_CHARACTER
    lens(1) = 3
    lens(2) = 6
    lens(3) = max_nc_filename_len
    CALL make_mpi_type(routine, count, types, lens, name_type_mpi)

    ! Initialize, read nml and pack into buffer
    IF (l_root) THEN
      DMAX = 0.0
      alat = 0.0
      alon = 0.0
      delta_lat = 2.5
      delta_lon = 3.75
      dscale = 0.0
      L_STRETCHGRID = .FALSE.
      L_REGGRID = .TRUE.
      L_VGRID_FILE = .FALSE.

      READ(nml_unit,NAME_DOMAIN)
      CALL mckpp_print(routine, "Read Namelist DOMAIN")
      
      IF (DMAX .LE. 0.0) THEN 
        CALL mckpp_abort(routine, "You must specify a depth for the domain")
      END IF
      IF ((L_STRETCHGRID) .AND. (dscale .EQ. 0.0)) THEN
        CALL mckpp_abort(routine, & 
                         "You cannot have dscale=0 for stretched grids") 
      END IF
 
      my_name % L_REGGRID = L_REGGRID
      my_name % L_STRETCHGRID = L_STRETCHGRID
      my_name % L_VGRID_FILE = L_VGRID_FILE
      my_name % dmax = dmax
      my_name % alat = alat
      my_name % alon = alon
      my_name % delta_lat = delta_lat
      my_name % delta_lon = delta_lon
      my_name % dscale = dscale
      my_name % vgrid_file = vgrid_file
    ENDIF

    ! Send data
    CALL mpi_bcast(my_name, 1, name_type_mpi, root, comm, ierr)
    CALL mpi_type_free(name_type_mpi, ierr)

    ! Unpack data
    IF (.NOT. l_root) THEN
      L_REGGRID = my_name % L_REGGRID
      L_STRETCHGRID = my_name % L_STRETCHGRID
      L_VGRID_FILE = my_name % L_VGRID_FILE
      dmax = my_name % dmax
      alat = my_name % alat
      alon = my_name % alon
      delta_lat = my_name % delta_lat
      delta_lon = my_name % delta_lon
      dscale = my_name % dscale
      vgrid_file = my_name % vgrid_file
    END IF

    IF (l_debug) THEN 
      WRITE(message,*) "DMAX, alon, alat, delta_lat, delta_lon = ", &
                        DMAX, alon, alat, delta_lat, delta_lon
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "L_STRETCHGRID, dscale, L_REGGRID, L_VGRID_FILE", &
                        L_STRETCHGRID, dscale, L_REGGRID, L_VGRID_FILE
      CALL mckpp_print(routine, message) 
      WRITE(message, *) "vgrid_file = ", TRIM(ADJUSTL(vgrid_file))
      CALL mckpp_print(routine, message) 
    END IF 
      
  END SUBROUTINE read_domain_namelist


  SUBROUTINE read_landsea_namelist(nml_unit)

    INTEGER, INTENT(IN) :: nml_unit
    CHARACTER(LEN=21) :: routine = "READ_LANDSEA_NAMELIST"
    CHARACTER(LEN=max_message_len) :: message

    TYPE name_type
      SEQUENCE
      LOGICAL :: L_LANDSEA
      CHARACTER(LEN=max_nc_filename_len) :: landsea_file = ""
    END TYPE name_type
    TYPE(name_type) my_name

    INTEGER :: count, name_type_mpi, ierr
    INTEGER, DIMENSION(2) :: types, lens

    ! Define MPI type
    count = 2
    types(1) = MPI_LOGICAL
    types(2) = MPI_CHARACTER
    lens(1) = 1
    lens(2) = max_nc_filename_len
    CALL make_mpi_type(routine, count, types, lens, name_type_mpi)

    ! Initialize, read nml and pack into buffer
    IF (l_root) THEN
      L_LANDSEA = .FALSE.
      
      READ(nml_unit,NAME_LANDSEA)
      CALL mckpp_print(routine, "Read Namelist LANDSEA")

      my_name % L_LANDSEA = L_LANDSEA
      my_name % landsea_file = landsea_file
    END IF

    ! Send data  
    CALL mpi_bcast(my_name, 1, name_type_mpi, root, comm, ierr)
    CALL mpi_type_free(name_type_mpi, ierr)

    ! Unpack data
    IF (.NOT. l_root) THEN
      L_LANDSEA = my_name % L_LANDSEA
      landsea_file = my_name % landsea_file
    END IF

    IF (l_debug) THEN 
      WRITE(message,*) "L_LANDSEA, landsea_file = ", &
                        L_LANDSEA, TRIM(landsea_file)
      CALL mckpp_print(routine, message)
    END IF 

  END SUBROUTINE read_landsea_namelist


  SUBROUTINE read_start_namelist(nml_unit)

    INTEGER, INTENT(IN) :: nml_unit
    CHARACTER(LEN=19) :: routine = "READ_START_NAMELIST"
    CHARACTER(LEN=max_message_len) :: message

    TYPE name_type
      SEQUENCE
      LOGICAL :: L_INITDATA, L_INTERPINIT, L_RESTART
      CHARACTER(max_nc_filename_len) :: initdata_file = ""
      CHARACTER(max_restart_filename_len) :: restart_infile = ""
    END TYPE name_type
    TYPE(name_type) my_name

    INTEGER :: count, name_type_mpi, ierr
    INTEGER, DIMENSION(2) :: types, lens

    ! Define MPI type
    count = 2
    types(1) = MPI_LOGICAL
    types(2) = MPI_CHARACTER
    lens(1) = 3
    lens(2) = max_nc_filename_len + max_restart_filename_len
    CALL make_mpi_type(routine, count, types, lens, name_type_mpi)

    ! Initialize, read nml and pack into buffer
    IF (l_root) THEN     
      L_INITDATA = .TRUE.
      L_INTERPINIT = .TRUE.
      L_RESTART = .FALSE.

      READ(nml_unit,NAME_START) 
      CALL mckpp_print(routine, "Read Namelist START")

      my_name % L_INITDATA = L_INITDATA
      my_name % L_INTERPINIT = L_INTERPINIT
      my_name % L_RESTART = L_RESTART
      my_name % initdata_file = initdata_file
      my_name % restart_infile = restart_infile
    ENDIF

    ! Send data
    CALL mpi_bcast(my_name, 1, name_type_mpi, root, comm, ierr)
    CALL mpi_type_free(name_type_mpi, ierr)

    ! Unpack data
    IF (.NOT. l_root) THEN
      L_INITDATA = my_name % L_INITDATA
      L_INTERPINIT = my_name % L_INTERPINIT
      L_RESTART = my_name % L_RESTART
      initdata_file = my_name % initdata_file
      restart_infile = my_name % restart_infile
    END IF

    IF (l_debug) THEN 
      WRITE(message,*) "L_INITDATA, initdata_file, L_INTERPINIT = ", & 
                        L_INITDATA, TRIM(initdata_file), L_INTERPINIT
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "L_RESTART, restart_infile = ", & 
                        L_RESTART, TRIM(restart_infile)
      CALL mckpp_print(routine, message)
    END IF 
 
  END SUBROUTINE read_start_namelist


  SUBROUTINE read_times_namelist(nml_unit)

    INTEGER, INTENT(IN) :: nml_unit
    CHARACTER(LEN=19) :: routine = "READ_TIMES_NAMELIST"
    CHARACTER(LEN=max_message_len) :: message

    TYPE name_type
      SEQUENCE
      REAL :: dtsec, startt, finalt
      INTEGER ndtocn, nyear
    END TYPE name_type
    TYPE(name_type) my_name

    INTEGER :: count, name_type_mpi, ierr
    INTEGER, DIMENSION(2) :: types, lens

    ! Define MPI type
    count = 2
    types(1) = MPI_REAL8
    types(2) = MPI_INTEGER
    lens(1) = 3
    lens(2) = 2
    CALL make_mpi_type(routine, count, types, lens, name_type_mpi)

    ! Initialize, read nml and pack into buffer
    IF (l_root) THEN
      ndtocn = 1
      dtsec = 0.0
      startt = -999.999
      finalt = -999.999

      READ(nml_unit,NAME_TIMES) 
      CALL mckpp_print(routine, "Read Namelist TIMES")
      
      IF ((dtsec .LE. 0.0) .OR. (startt .LT. 0.0) .OR. (finalt .LT. 0.0)) THEN 
        CALL mckpp_abort(routine, &
        "You must specify values of dtsec,startt,finalt in the namelist")
      ENDIF

      my_name % dtsec = dtsec
      my_name % startt = startt
      my_name % finalt = finalt
      my_name % ndtocn = ndtocn
      my_name % nyear = nyear
    END IF

    ! Send data
    CALL mpi_bcast(my_name, 1, name_type_mpi, root, comm, ierr)
    CALL mpi_type_free(name_type_mpi, ierr)

    ! Unpack data
    IF (.NOT. l_root) THEN
      dtsec = my_name % dtsec
      startt = my_name % startt
      finalt = my_name % finalt
      ndtocn = my_name % ndtocn
      nyear = my_name % nyear
    END IF

    IF (l_debug) THEN 
      WRITE(message,*) "dtsec, startt, finalt, ndtocn, nyear = ", &
                        dtsec, startt, finalt, ndtocn, nyear
      CALL mckpp_print(routine, message)
    END IF 

  END SUBROUTINE read_times_namelist


  SUBROUTINE read_couple_namelist(nml_unit)

    INTEGER, INTENT(IN) :: nml_unit
    CHARACTER(LEN=20) :: routine = "READ_COUPLE_NAMELIST"
    CHARACTER(LEN=max_message_len) :: message

    TYPE name_type
      SEQUENCE
      LOGICAL :: L_COUPLE, L_CLIMSST, L_UPD_CLIMSST, L_CPLWGHT, L_CLIMICE, & 
                 L_UPD_CLIMICE, L_CLIM_ICE_DEPTH, L_CLIM_SNOW_ON_ICE, & 
                 L_OUTKELVIN, L_COUPLE_CURRENTS, L_CLIMCURR, L_UPD_CLIMCURR, & 
                 L_PERIODIC_CLIMICE, L_PERIODIC_CLIMSST, L_BAD_ICE_DEPTH
      INTEGER :: ifirst, ilast, jfirst, jlast, ndtupdsst, ndtupdice, & 
                 ndtupdcurr, climsst_period, climice_period 
      CHARACTER(LEN=max_nc_filename_len) :: sstin_file = "", & 
          cplwght_file = "", icein_file = "", currin_file = ""
    END TYPE name_type
    TYPE(name_type) my_name

    INTEGER :: count, name_type_mpi, ierr
    INTEGER, DIMENSION(3) :: types, lens

    ! Define MPI type
    count = 3
    types(1) = MPI_LOGICAL
    types(2) = MPI_INTEGER
    types(3) = MPI_CHARACTER
    lens(1) = 15
    lens(2) = 9
    lens(3) = max_nc_filename_len * 4
    CALL make_mpi_type(routine, count, types, lens, name_type_mpi)

    ! Initialize, read nml and pack into buffer
    IF (l_root) THEN
#if defined MCKPP_COUPLE
      L_COUPLE = .TRUE.
#else
      L_COUPLE = .FALSE.
#endif
      L_COUPLE_CURRENTS = .FALSE.
      L_OUTKELVIN = .FALSE.
      L_UPD_CLIMSST = .FALSE.
      L_UPD_CLIMICE = .FALSE.
      L_CLIMICE = .FALSE.
      L_CLIMSST = .FALSE.
      L_CLIMCURR = .FALSE. 
      L_BAD_ICE_DEPTH = .FALSE.
      ifirst = 1
      ilast = nx
      jfirst = 1
      jlast = ny
       
      READ(nml_unit,NAME_COUPLE)
      CALL mckpp_print(routine, "Read Namelist COUPLE")

      my_name % L_COUPLE = L_COUPLE
      my_name % L_CLIMSST = L_CLIMSST
      my_name % L_UPD_CLIMSST = L_UPD_CLIMSST
      my_name % L_CPLWGHT = L_CPLWGHT
      my_name % L_CLIMICE = L_CLIMICE
      my_name % L_UPD_CLIMICE = L_UPD_CLIMICE
      my_name % L_CLIM_ICE_DEPTH = L_CLIM_ICE_DEPTH
      my_name % L_CLIM_SNOW_ON_ICE = L_CLIM_SNOW_ON_ICE
      my_name % L_OUTKELVIN = L_OUTKELVIN
      my_name % L_COUPLE_CURRENTS = L_COUPLE_CURRENTS
      my_name % L_CLIMCURR = L_CLIMCURR
      my_name % L_UPD_CLIMCURR = L_UPD_CLIMCURR
      my_name % L_PERIODIC_CLIMICE = L_PERIODIC_CLIMICE
      my_name % L_PERIODIC_CLIMSST = L_PERIODIC_CLIMSST
      my_name % L_BAD_ICE_DEPTH = L_BAD_ICE_DEPTH
      my_name % ifirst = ifirst
      my_name % ilast = ilast
      my_name % jfirst = jfirst
      my_name % jlast = jlast
      my_name % ndtupdsst = ndtupdsst
      my_name % ndtupdice = ndtupdice
      my_name % ndtupdcurr = ndtupdcurr
      my_name % climsst_period = climsst_period
      my_name % climice_period = climice_period
      my_name % sstin_file = sstin_file
      my_name % cplwght_file = cplwght_file
      my_name % icein_file = icein_file
      my_name % currin_file = currin_file 
    ENDIF

    ! Send data
    CALL mpi_bcast(my_name, 1, name_type_mpi, root, comm, ierr)
    CALL mpi_type_free(name_type_mpi, ierr)

    ! Unpack data
    IF (.NOT. l_root) THEN
      L_COUPLE = my_name % L_COUPLE
      L_CLIMSST = my_name % L_CLIMSST
      L_UPD_CLIMSST= my_name % L_UPD_CLIMSST
      L_CPLWGHT = my_name % L_CPLWGHT
      L_CLIMICE = my_name % L_CLIMICE
      L_UPD_CLIMICE = my_name % L_UPD_CLIMICE
      L_CLIM_ICE_DEPTH = my_name % L_CLIM_ICE_DEPTH
      L_CLIM_SNOW_ON_ICE = my_name % L_CLIM_SNOW_ON_ICE
      L_OUTKELVIN = my_name % L_OUTKELVIN
      L_COUPLE_CURRENTS = my_name % L_COUPLE_CURRENTS
      L_CLIMCURR = my_name % L_CLIMCURR
      L_UPD_CLIMCURR = my_name % L_UPD_CLIMCURR
      L_PERIODIC_CLIMICE = my_name % L_PERIODIC_CLIMICE
      L_PERIODIC_CLIMSST = my_name % L_PERIODIC_CLIMSST
      L_BAD_ICE_DEPTH = my_name % L_BAD_ICE_DEPTH
      ifirst = my_name % ifirst
      ilast = my_name % ilast
      jfirst = my_name % jfirst
      jlast = my_name % jlast
      ndtupdsst = my_name % ndtupdsst
      ndtupdice = my_name % ndtupdice
      ndtupdcurr = my_name % ndtupdcurr
      climsst_period = my_name % climsst_period
      climice_period = my_name % climice_period
      sstin_file = my_name % sstin_file
      cplwght_file = my_name % cplwght_file
      icein_file = my_name % icein_file
      currin_file = my_name % currin_file
    END IF

    IF (l_debug) THEN 
      WRITE(message,*) "L_COUPLE, ifirst, ilast, jfirst, jlast, L_CLIMSST = ", &
                        L_COUPLE, ifirst, ilast, jfirst, jlast, L_CLIMSST
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "sstin_file, L_UPD_CLIMSST, ndtupdsst, L_CPLWGHT = ", & 
                        TRIM(sstin_file), L_UPD_CLIMSST, ndtupdsst, L_CPLWGHT
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "cplwght_file = ", TRIM(ADJUSTL(cplwght_file))
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "icein_file = ", TRIM(ADJUSTL(icein_file))
      CALL mckpp_print(routine, message)
      WRITE(message,*) "L_CLIMICE, L_UPD_CLIMICE, ndtupdice = ", &
                        L_CLIMICE, L_UPD_CLIMICE, ndtupdice
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "L_CLIM_ICE_DEPTH, L_CLIM_SNOW_ON_ICE, L_OUTKELVIN = ", &
                      L_CLIM_ICE_DEPTH, L_CLIM_SNOW_ON_ICE, L_OUTKELVIN
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "currin_file = ", TRIM(currin_file)
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "L_COUPLE_CURRENTS, L_CLIMCURR, L_UPD_CLIMCURR = ", & 
                        L_COUPLE_CURRENTS, L_CLIMCURR, L_UPD_CLIMCURR
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "ndtupdcurr, L_PERIODIC_CLIMICE, L_PERIODIC_CLIMSST = ", &
                        ndtupdcurr, L_PERIODIC_CLIMICE, L_PERIODIC_CLIMSST
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "climsst_period, climice_period, L_BAD_ICE_DEPTH = ", & 
                        climsst_period, climice_period, L_BAD_ICE_DEPTH
      CALL mckpp_print(routine, message)
    END IF 

  END SUBROUTINE read_couple_namelist


  SUBROUTINE read_advec_namelist(nml_unit)

    INTEGER, INTENT(IN) :: nml_unit
    CHARACTER(LEN=19) :: routine = "READ_ADVEC_NAMELIST"
    CHARACTER(LEN=max_message_len) :: message

    TYPE name_type
      SEQUENCE
      LOGICAL :: L_ADVECT, L_RELAX_SST, L_RELAX_CALCONLY, L_RELAX_SAL, & 
                 L_RELAX_OCNT
      CHARACTER(max_restart_filename_len) :: advect_file = ""
    END TYPE name_type
    TYPE(name_type) my_name

    INTEGER :: count, name_type_mpi, ierr
    INTEGER, DIMENSION(2) :: types, lens

    ! Define MPI type
    count = 2
    types(1) = MPI_LOGICAL
    types(2) = MPI_CHARACTER
    lens(1) = 5
    lens(2) = max_restart_filename_len
    CALL make_mpi_type(routine, count, types, lens, name_type_mpi)

    ALLOCATE( relax_sst_in(ny) )
    ALLOCATE( relax_sal_in(ny) )
    ALLOCATE( relax_ocnt_in(ny) )

    ! Initialize, read nml and pack into buffer
    IF (l_root) THEN
      L_ADVECT = .FALSE.
      L_RELAX_SST = .FALSE.
      L_RELAX_CALCONLY = .FALSE.
      relax_sst_in = 0.0
      relax_sal_in = 0.0
      relax_ocnt_in = 0.0

      READ(nml_unit,NAME_ADVEC)
      CALL mckpp_print(routine, "Read Namelist ADVEC")

      my_name % L_ADVECT = L_ADVECT
      my_name % L_RELAX_SST = L_RELAX_SST
      my_name % L_RELAX_CALCONLY = L_RELAX_CALCONLY
      my_name % L_RELAX_SAL = L_RELAX_SAL
      my_name % L_RELAX_OCNT = L_RELAX_OCNT
      my_name % advect_file = advect_file
    ENDIF

    ! Send data
    CALL mpi_bcast(my_name, 1, name_type_mpi, root, comm, ierr)
    CALL mpi_type_free(name_type_mpi, ierr)

    ! Need to send allocatable arrays separately
    CALL mpi_bcast(relax_sst_in, ny, MPI_REAL8, root, comm, ierr)
    CALL mpi_bcast(relax_sal_in, ny, MPI_REAL8, root, comm, ierr)
    CALL mpi_bcast(relax_ocnt_in, ny, MPI_REAL8, root, comm, ierr)

    ! Unpack data
    IF (.NOT. l_root) THEN
      L_ADVECT = my_name % L_ADVECT
      L_RELAX_SST = my_name % L_RELAX_SST
      L_RELAX_CALCONLY = my_name % L_RELAX_CALCONLY
      L_RELAX_SAL = my_name % L_RELAX_SAL
      L_RELAX_OCNT = my_name % L_RELAX_OCNT
      advect_file = my_name % advect_file
    END IF

    IF (l_debug) THEN 
      WRITE(message,*) "L_ADVECT, advect_file, L_RELAX_SST = ", &
                        L_ADVECT, TRIM(advect_file), L_RELAX_SST
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "L_RELAX_CALCONLY, L_RELAX_SAL, L_RELAX_OCNT = ", &
                        L_RELAX_CALCONLY, L_RELAX_SAL, L_RELAX_OCNT
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "relax_sst_in = ", relax_sst_in
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "relax_sal_in = ", relax_sal_in
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "relax_ocnt_in = ", relax_ocnt_in
      CALL mckpp_print(routine, message) 
    END IF 
      
  END SUBROUTINE read_advec_namelist


  SUBROUTINE read_paras_namelist(nml_unit)

    INTEGER, INTENT(IN) :: nml_unit
    CHARACTER(LEN=19) :: routine = "READ_PARAS_NAMELIST"
    CHARACTER(LEN=max_message_len) :: message
 
    TYPE name_type
      SEQUENCE
      LOGICAL :: L_JERLOV
      CHARACTER(LEN=max_nc_filename_len) :: paras_file = ""
    END TYPE name_type
    TYPE(name_type) my_name

    INTEGER :: count, name_type_mpi, ierr
    INTEGER, DIMENSION(2) :: types, lens

    ! Define MPI type
    count = 2
    types(1) = MPI_LOGICAL
    types(2) = MPI_CHARACTER
    lens(1) = 1
    lens(2) = max_nc_filename_len
    CALL make_mpi_type(routine, count, types, lens, name_type_mpi)

    ! Initialize, read nml and pack into buffer
    IF (l_root) THEN 
      L_JERLOV = .TRUE.
 
      READ(nml_unit,NAME_PARAS)
      CALL mckpp_print(routine, "Read Namelist PARAS")

      my_name % L_JERLOV = L_JERLOV
      my_name % paras_file = paras_file
    END IF

    ! Send data  
    CALL mpi_bcast(my_name, 1, name_type_mpi, root, comm, ierr)
    CALL mpi_type_free(name_type_mpi, ierr)

    ! Unpack data
    IF (.NOT. l_root) THEN
      L_JERLOV = my_name % L_JERLOV
      paras_file = my_name % paras_file
    END IF

    IF (l_debug) THEN 
      WRITE(message,*) "paras_file, L_JERLOV = ", TRIM(paras_file), L_JERLOV
      CALL mckpp_print(routine, message) 
    END IF 
      
  END SUBROUTINE read_paras_namelist


  SUBROUTINE read_forcing_namelist(nml_unit)

    INTEGER, INTENT(IN) :: nml_unit
    CHARACTER(LEN=21) :: routine = "READ_FORCING_NAMELIST"
    CHARACTER(LEN=max_message_len) :: message
 
    TYPE name_type
      SEQUENCE
      LOGICAL :: L_FCORR_WITHZ, L_FCORR, L_UPD_FCORR, L_PERIODIC_FCORR, & 
                 L_NO_FREEZE, L_NO_ISOTHERM, L_FLUXDATA, L_REST, L_UPD_SAL, & 
                 L_PERIODIC_SAL, L_INTERP_SAL, L_UPD_OCNT, L_PERIODIC_OCNT, & 
                 L_INTERP_OCNT, L_SFCORR_WITHZ, L_SFCORR, L_UPD_SFCORR, & 
                 L_PERIODIC_SFCORR, L_DAMP_CURR, L_VARY_BOTTOM_TEMP, &
                 L_UPD_BOTTOM_TEMP, L_PERIODIC_BOTTOM_TEMP
      INTEGER :: ndtupdfcorr, fcorr_period, isotherm_bottom, ndtupdsal, & 
                 sal_period, ndt_interp_sal, ndtupdocnt, ocnt_period, & 
                 ndt_interp_ocnt, ndtupdsfcorr, sfcorr_period, ndtupdbottom, & 
                 bottom_temp_period, dtuvdamp
      REAL :: isotherm_threshold
      CHARACTER (LEN=max_nc_filename_len) :: fcorrin_file = "", & 
          forcing_file = "", sal_file = "", ocnT_file = "", & 
          sfcorrin_file = "", bottomin_file = ""
    END TYPE name_type
    TYPE(name_type) my_name

    INTEGER :: count, name_type_mpi, ierr
    INTEGER, DIMENSION(4) :: types, lens

    ! Define MPI type
    count = 4
    types(1) = MPI_LOGICAL
    types(2) = MPI_INTEGER
    types(3) = MPI_REAL8
    types(4) = MPI_CHARACTER
    lens(1) = 22
    lens(2) = 14
    lens(3) = 1 
    lens(4) = max_nc_filename_len * 6
    CALL make_mpi_type(routine, count, types, lens, name_type_mpi)

    ! Initialize, read nml and pack into buffer
    IF (l_root) THEN
      L_FLUXDATA = .FALSE.
      L_FCORR_WITHZ = .FALSE.
      L_FCORR = .FALSE.
      L_UPD_FCORR = .FALSE.
      L_SFCORR_WITHZ = .FALSE.
      L_SFCORR = .FALSE.
      L_UPD_SFCORR = .FALSE.
      L_UPD_SAL = .FALSE.
      L_VARY_BOTTOM_TEMP = .FALSE.
      L_UPD_BOTTOM_TEMP = .FALSE.
      L_REST = .FALSE.
      L_NO_FREEZE = .FALSE.
      L_NO_ISOTHERM = .FALSE.
      L_DAMP_CURR = .FALSE.

      READ(nml_unit,NAME_FORCING)  
      CALL mckpp_print(routine, "Read Namelist FORCING")

      IF (L_FCORR_WITHZ .AND. L_FCORR) THEN
        WRITE(message, *) "L_FCORR and L_FCORR_WITHZ are mutually ", & 
                          "exclusive. Choose one or neither."
        CALL mckpp_abort(routine, message)
      ENDIF
      IF (L_SFCORR_WITHZ .AND. L_SFCORR) THEN
        WRITE(message, *) "L_SFCORR and L_SFCORR_WITHZ are mutually ", & 
                          "exclusive. Choose one or neither."
        CALL mckpp_abort(routine, message)
      ENDIF
      IF (L_FCORR_WITHZ .AND. L_RELAX_SST) THEN
        WRITE(message, *) "L_FCORR_WITHZ and L_RELAX_SST are mutually ", & 
                          "exclusive. Choose one or neither."
        CALL mckpp_abort(routine, message)
      ENDIF
      IF (L_NO_ISOTHERM .AND. (ocnT_file .eq. 'none' .or.&
          sal_file .eq. 'none')) THEN
        WRITE(message, *) "If you specify L_NO_ISOTHERM for reseting of ", & 
                          "isothermal points, you must specify files from ", & 
                          "which to read climatological ocean temperature ", &
                          "(ocnT_file) and salinity (sal_file)."
        CALL mckpp_abort(routine, message)
      ENDIF

      my_name % L_FCORR_WITHZ = L_FCORR_WITHZ
      my_name % L_FCORR = L_FCORR
      my_name % L_UPD_FCORR = L_UPD_FCORR
      my_name % L_PERIODIC_FCORR = L_PERIODIC_FCORR
      my_name % L_NO_FREEZE = L_NO_FREEZE
      my_name % L_NO_ISOTHERM = L_NO_ISOTHERM
      my_name % L_FLUXDATA = L_FLUXDATA
      my_name % L_REST = L_REST
      my_name % L_UPD_SAL = L_UPD_SAL
      my_name % L_PERIODIC_SAL = L_PERIODIC_SAL
      my_name % L_INTERP_SAL = L_INTERP_SAL
      my_name % L_UPD_OCNT = L_UPD_OCNT
      my_name % L_PERIODIC_OCNT = L_PERIODIC_OCNT
      my_name % L_INTERP_OCNT = L_INTERP_OCNT
      my_name % L_SFCORR_WITHZ = L_SFCORR_WITHZ
      my_name % L_SFCORR = L_SFCORR
      my_name % L_UPD_SFCORR = L_UPD_SFCORR
      my_name % L_PERIODIC_SFCORR = L_PERIODIC_SFCORR
      my_name % L_DAMP_CURR = L_DAMP_CURR
      my_name % L_VARY_BOTTOM_TEMP = L_VARY_BOTTOM_TEMP
      my_name % L_UPD_BOTTOM_TEMP = L_UPD_BOTTOM_TEMP
      my_name % L_PERIODIC_BOTTOM_TEMP = L_PERIODIC_BOTTOM_TEMP
      my_name % ndtupdfcorr = ndtupdfcorr
      my_name % fcorr_period = fcorr_period
      my_name % isotherm_bottom = isotherm_bottom
      my_name % ndtupdsal = ndtupdsal
      my_name % sal_period = sal_period
      my_name % ndt_interp_sal = ndt_interp_sal
      my_name % ndtupdocnt = ndtupdocnt
      my_name % ocnt_period = ocnt_period
      my_name % ndt_interp_ocnt = ndt_interp_ocnt
      my_name % ndtupdsfcorr = ndtupdsfcorr
      my_name % sfcorr_period = sfcorr_period
      my_name % ndtupdbottom = ndtupdbottom
      my_name % bottom_temp_period = bottom_temp_period
      my_name % dtuvdamp = dtuvdamp
      my_name % isotherm_threshold = isotherm_threshold
      my_name % fcorrin_file = fcorrin_file
      my_name % forcing_file = forcing_file
      my_name % sal_file = sal_file
      my_name % ocnT_file = ocnT_file
      my_name % sfcorrin_file = sfcorrin_file
      my_name % bottomin_file = bottomin_file
    END IF

    ! Send data  
    CALL mpi_bcast(my_name, 1, name_type_mpi, root, comm, ierr)
    CALL mpi_type_free(name_type_mpi, ierr)

    ! Unpack data
    IF (.NOT. l_root) THEN
      L_FCORR_WITHZ = my_name % L_FCORR_WITHZ
      L_FCORR = my_name % L_FCORR
      L_UPD_FCORR = my_name % L_UPD_FCORR
      L_PERIODIC_FCORR = my_name % L_PERIODIC_FCORR
      L_NO_FREEZE = my_name % L_NO_FREEZE
      L_NO_ISOTHERM = my_name % L_NO_ISOTHERM
      L_FLUXDATA = my_name % L_FLUXDATA
      L_REST = my_name % L_REST
      L_UPD_SAL = my_name % L_UPD_SAL
      L_PERIODIC_SAL = my_name % L_PERIODIC_SAL
      L_INTERP_SAL = my_name % L_INTERP_SAL
      L_UPD_OCNT = my_name % L_UPD_OCNT
      L_PERIODIC_OCNT = my_name % L_PERIODIC_OCNT
      L_INTERP_OCNT = my_name % L_INTERP_OCNT
      L_SFCORR_WITHZ = my_name % L_SFCORR_WITHZ
      L_SFCORR = my_name % L_SFCORR
      L_UPD_SFCORR = my_name % L_UPD_SFCORR
      L_PERIODIC_SFCORR = my_name % L_PERIODIC_SFCORR
      L_DAMP_CURR = my_name % L_DAMP_CURR
      L_VARY_BOTTOM_TEMP = my_name % L_VARY_BOTTOM_TEMP
      L_UPD_BOTTOM_TEMP = my_name % L_UPD_BOTTOM_TEMP
      L_PERIODIC_BOTTOM_TEMP = my_name % L_PERIODIC_BOTTOM_TEMP    
      ndtupdfcorr = my_name % ndtupdfcorr
      fcorr_period = my_name % fcorr_period 
      isotherm_bottom = my_name % isotherm_bottom
      ndtupdsal = my_name % ndtupdsal 
      sal_period = my_name % sal_period
      ndt_interp_sal = my_name % ndt_interp_sal
      ndtupdocnt = my_name % ndtupdocnt
      ocnt_period = my_name % ocnt_period  
      ndt_interp_ocnt = my_name % ndt_interp_ocnt     
      ndtupdsfcorr = my_name % ndtupdsfcorr  
      sfcorr_period = my_name % sfcorr_period
      ndtupdbottom = my_name % ndtupdbottom    
      bottom_temp_period = my_name % bottom_temp_period 
      dtuvdamp = my_name % dtuvdamp    
      isotherm_threshold = my_name % isotherm_threshold
      fcorrin_file = my_name % fcorrin_file    
      forcing_file = my_name % forcing_file  
      sal_file = my_name % sal_file    
      ocnT_file = my_name % ocnT_file
      sfcorrin_file = my_name % sfcorrin_file    
      bottomin_file = my_name % bottomin_file    
    END IF

    IF (l_debug) THEN 
      WRITE(message,*) "L_FCORR_WITHZ, L_FCORR, L_UPD_FCORR = ", & 
                        L_FCORR_WITHZ, L_FCORR, L_UPD_FCORR
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "L_PERIODIC_FCORR, L_NO_FREEZE, L_NO_ISOTHERM = ", & 
                        L_PERIODIC_FCORR, L_NO_FREEZE, L_NO_ISOTHERM
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "L_FLUXDATA, L_REST, L_UPD_SAL, L_PERIODIC_SAL = ", & 
                        L_FLUXDATA, L_REST, L_UPD_SAL, L_PERIODIC_SAL
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "L_INTERP_SAL, L_UPD_OCNT, L_PERIODIC_OCNT = ", & 
                        L_INTERP_SAL, L_UPD_OCNT, L_PERIODIC_OCNT
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "L_INTERP_OCNT, L_SFCORR_WITHZ, L_SFCORR = ", & 
                        L_INTERP_OCNT, L_SFCORR_WITHZ, L_SFCORR
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "L_UPD_SFCORR, L_PERIODIC_SFCORR, L_DAMP_CURR = ", & 
                        L_UPD_SFCORR, L_PERIODIC_SFCORR, L_DAMP_CURR
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "L_VARY_BOTTOM_TEMP, L_UPD_BOTTOM_TEMP = ", & 
                        L_VARY_BOTTOM_TEMP, L_UPD_BOTTOM_TEMP
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "L_PERIODIC_BOTTOM_TEMP = ", & 
                        L_PERIODIC_BOTTOM_TEMP
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "ndtupdfcorr, fcorr_period, isotherm_bottom = ", &
                        ndtupdfcorr, fcorr_period, isotherm_bottom
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "ndtupdsal, sal_period, ndt_interp_sal, ndtupdocnt = ", &
                        ndtupdsal, sal_period, ndt_interp_sal, ndtupdocnt
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "ocnt_period, ndt_interp_ocnt, isotherm_threshold = ", &
                        ocnt_period, ndt_interp_ocnt, isotherm_threshold
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "isotherm_threshold = ", isotherm_threshold
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "fcorrin_file, forcing_file, sal_file = ", &
                        TRIM(fcorrin_file), TRIM(forcing_file), TRIM(sal_file)
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "ocnT_file, sfcorrin_file, bottomin_file = ", &
                       TRIM(ocnT_file), TRIM(sfcorrin_file), TRIM(bottomin_file)
      CALL mckpp_print(routine, message)
    END IF 

  END SUBROUTINE read_forcing_namelist


  SUBROUTINE read_output_namelist(nml_unit)

    INTEGER, INTENT(IN) :: nml_unit
    CHARACTER(LEN=20) :: routine = "READ_OUTPUT_NAMELIST"
    CHARACTER(LEN=max_message_len) :: message
 
    TYPE name_type
      SEQUENCE
      LOGICAL :: L_RESTARTW
      INTEGER :: ndt_per_restart
      CHARACTER(max_restart_filename_len) :: restart_outfile = ""
    END TYPE name_type
    TYPE(name_type) my_name

    INTEGER :: count, name_type_mpi, ierr
    INTEGER, DIMENSION(3) :: types, lens

    ! Define MPI type
    count = 3
    types(1) = MPI_LOGICAL
    types(2) = MPI_INTEGER
    types(3) = MPI_CHARACTER
    lens(1) = 1
    lens(2) = 2
    lens(3) = max_nc_filename_len
    CALL make_mpi_type(routine, count, types, lens, name_type_mpi)

    ! Initialize, read nml and pack into buffer
    IF (l_root) THEN
      L_RESTARTW = .TRUE.
      
      READ(nml_unit,NAME_OUTPUT)
      CALL mckpp_print(routine, "Read Namelist OUTPUT") 

      my_name % L_RESTARTW = L_RESTARTW
      my_name % ndt_per_restart = ndt_per_restart
      my_name % restart_outfile = restart_outfile
    END IF

    ! Send data  
    CALL mpi_bcast(my_name, 1, name_type_mpi, root, comm, ierr)
    CALL mpi_type_free(name_type_mpi, ierr)

    ! Unpack data
    IF (.NOT. l_root) THEN
      L_RESTARTW = my_name % L_RESTARTW
      ndt_per_restart = my_name % ndt_per_restart
      restart_outfile = my_name % restart_outfile
    END IF

    IF (l_debug) THEN 
      WRITE(message,*) "L_RESTARTW, restart_outfile, ndt_per_restart = ", &
                        L_RESTARTW, restart_outfile, ndt_per_restart
      CALL mckpp_print(routine, message) 
    END IF 
      
  END SUBROUTINE read_output_namelist


  SUBROUTINE make_mpi_type(calling_routine, count, types, blocklens, & 
                           new_mpi_type)

    CHARACTER(LEN=*), INTENT(IN) ::  calling_routine
    INTEGER, INTENT(IN) :: count
    INTEGER, DIMENSION(count), INTENT(IN) :: types, blocklens
    INTEGER, INTENT(OUT) :: new_mpi_type

    INTEGER :: i, ierr
    INTEGER(KIND=MPI_ADDRESS_KIND), DIMENSION(count) :: displs
    INTEGER(KIND=MPI_ADDRESS_KIND) :: extent, lb 
    CHARACTER(LEN=13) :: routine = "MAKE_MPI_TYPE"
    CHARACTER(LEN=max_message_len) :: context, message

    context = update_context(calling_routine, routine) 
    IF (count .LT. 1) THEN
      CALL mckpp_abort(context, & 
                       "Error in initialize_namelists creating MPI type.")
    END IF

    ! Calculate displacements
    displs(1) = 0
    IF (count .GT. 1) THEN 
      DO i = 2, count
        CALL mpi_type_get_extent(types(i-1), lb, extent, ierr)
        displs(i) = displs(i-1) + extent * blocklens(i-1) 
      END DO
    END IF

    ! Create new type
    CALL mpi_type_create_struct(count, blocklens, displs, types, & 
                                new_mpi_type, ierr)
    CALL mpi_type_commit(new_mpi_type, ierr)

  END SUBROUTINE make_mpi_type

END MODULE mckpp_initialize_namelist_mod
