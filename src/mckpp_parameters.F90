MODULE mckpp_parameters

  INTEGER, PARAMETER :: & 
       max_nc_filename_len = 50, &
       max_restart_filename_len = 50

  ! These are read in from parameters namelist
  INTEGER :: & 

       ! Main (permanent) grid
       nz, &           ! number of layers      
       nzm1, &         ! = nz - 1 
       nzp1, &         ! = nz + 1 (number of grid points) 
       ndim, &         ! dimension of the model (=1) 
       nx, &           ! 
       ny, &           ! 
       npts, &         ! = nx * ny 
       nvel, &         ! number of velocity components, i.e. 2 for U and V
       nsclr, &        ! number of scalars, i.e. T, S, and additional scalars
       nvp1, &         ! = nvel + 1
       nsp1, &         ! = nsclr + 1 
       nsb, &          ! number of biological scalars, used in "biocommon.inc"
       itermax         ! maximum number of hmix iterations 
                       !  (on main or temporary grids).
   REAL :: &
       hmixtolfrac     ! convergence tolerance for hmix(new)-hmix(old) 
                       !  iteration in ocnstep: fraction of layer thickness
                       !  hm(kmix)
   INTEGER :: & 

       ! Temporary grid 
       ngrid, &        ! number of grids = permanent grid + temporary grids,
                       !  if only p-grid used: ngrid = 1, 
                       !  if t-grid used: ngrid = nz
       nzl, &          ! refinement interval: it is defined on permanent-grid
       nzu, &          !  from (kmix+nzu) to (kmix-nzl)
       nzdivmax, &     ! maximum number of fine(temporary) grid intervals
                       !  per permanent grid interval. It is needed to dimension 
                       !  the necessary arrays; the number to be used is being 
                       !  read in as nzdiv.
       nztmax, &       ! maximum number of layers on temporary grid,
                       !  the number to be used in this run is read in as nzt. 
       nzp1tmax, &     ! = nztmax + 1 
       igridmax, &     ! maximum number of new grid refinements in ocntgrid

       ! Fluxes and forcing 
       nsflxs, &       ! number of fluxes: sflux(nsflxs,5,0:njdt)
       njdt, &         ! number of older flux values used to extrapolate new fluxes
       nsflxsm1, &     ! = nsflxs - 1 
       nsflxsp2, &     ! = nsflxs + 2 
       ndharm, &       ! maximum number of harmonics used to specify forcing

       ! Ocean advection
       maxmodeadv, &   ! maximum number of different modes for advection

       ! Richardson mixing
       mr, &           ! dimension of "F of Ri" in ri_mix routine       
       mrp1, &         ! = mr + 1 

       ! Parameters for regional coupling 
       nx_globe, &  
       ny_globe, &   
       npts_globe      ! = nx_globe * ny_globe 

END MODULE mckpp_parameters
