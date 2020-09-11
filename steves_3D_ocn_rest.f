      PROGRAM ocn_model_3D
**************************************************************************
* 3D version of the 
* 1D ocean model using the kpp mixing scheme of Large et al, with his 
* interface between the ocean and kpp scheme
* calls his subroutines
* init_ocn : to initialize the ocean model
* ocn_step : to update the model
*
* Also uses his parameter namelists and common blocks
* There maybe a lot of unnecessary variables and parameters here, but until
* we get something working it seems foolish to try and strip too much out.
* Written April 2002
*
*  Steve Woolnough
**************************************************************************

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
      
#ifdef COUPLE
#ifdef CFS
      include 'parameter.cfs_coupled.inc'
#else
      include 'parameter.oasis2.inc'
#endif
#else
#ifdef CFS
      include 'parameter.cfs_forced.inc'
#else
      include 'parameter.forced.inc'
#endif
#endif
      include 'param.h'
      include 'constants.com'
      include 'times.com'
      include 'timocn.com'

      include 'flx_sfc.com'
      include 'flx_profs.com'
      include 'vert_pgrid.com'
      include 'ocn_paras.com'
      include 'kprof_out.com'
      include 'dble_diff.com'
      include 'local_pt.com'
      include 'output.com'
      include 'couple.com'
      
* Local variables
      REAL X(NPTS,NZP1,NSCLR),U(NPTS,NZP1,NVEL) ! Scalar and wind fields
      
      INTEGER k,nflx,nstep
* Initialize the model 
* Setup the constants, read the namelists, setup the initial conditions
      call initialize(U,X)
      
      write(nuout,*) 'Model initialized with, ',NZ,' levels'
      write(nuout,*) 'and ',npts,' points'
c      write(nuout,*) 'Initial conditions are given by'
c      write(nuout,*) ' k     z(m)    U(m/s)   V(m/s)    T(C)   S(o/oo)'
c      DO k=1,NZP1
c         print('(i3,2x,f8.3,2(2x,f7.5,),2x,f7.4,2x,f7.4)'),
c     &        k,zm(k),U(k,1),U(k,2),X(k,1),X(k,2)+Sref
c      ENDDO

      IF (L_COUPLE) THEN
         IF ((im .NE. NX_GLOBE) .OR. (jm .NE. NY_GLOBE)) THEN
***       Test that the global ocean size is the same in both definitions
***       lazy but probably easier than rewriting the bits of code
             call halte('im,jm not equal to NX_GLOBE,NY_GLOBE')
         ENDIF
         call inicmo(nend*ndtocn,ndtocn,int(dto))
      ENDIF

c     DO nflx=1,nend
      write(nuout,*) ' Starting the ocean integration' !CSJW
      DO ntime=1,nend*ndtocn
          
c     time=startt+ntime*dto/spd
         IF (MOD(ntime-1,ndtocn) .EQ. 0) THEN
            call fluxes
         ENDIF
c         DO nstep=1,ndtocn
c            ntime=(nflx-1)*ndtocn+nstep
            time=startt+ntime*dto/spd
            DO ipt=1,npts
               call ocnstep(U,X)
            ENDDO
            IF (MOD(ntime,ndtout) .EQ. 0) THEN
c            IF ( abs(float(ntime/ndtout)-float(ntime)/float(ndtout)) 
c     &           .LT. 1.e-3) THEN
               IF (time .GT. float(day_out)) THEN
                  day_out=day_out+5
                  write(output_file(flen+2:flen+5),'(i4.4)'),day_out
                  call init_output
               ENDIF
               call output_inst(U,X)
            ENDIF  
            IF ((MOD(ntime,ndtocn) .EQ. 0) .AND. L_COUPLE ) THEN
               call coupled_out (X,ntime,.FALSE.)
            ENDIF
c         ENDDO
         
      ENDDO
      
c      call output_close
      
      write(nuout,*) 'Successful termination of the model integration'
   
      IF (L_RESTARTW) THEN 
         CALL WRITE_RESTART(U,X)
      ENDIF

      IF (L_COUPLE) THEN
         call coupled_out (X,ntime,.TRUE.)
      ENDIF

      END

      SUBROUTINE initialize(U,X)
************************************************************************
*     Subroutine to initialize the model, some of the output is passed
*     through the common blocks
************************************************************************
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#ifdef COUPLE
#ifdef CFS
      include 'parameter.cfs_coupled.inc'
#else
      include 'parameter.oasis2.inc'
#endif
#else
#ifdef CFS
      include 'parameter.cfs_forced.inc'
#else
      include 'parameter.forced.inc'
#endif
#endif
      include 'constants.com'
      include 'times.com'
      include 'timocn.com'
      include 'location.com'
      include 'vert_pgrid.com'
      include 'proc_swit.com'
      include 'proc_pars.com'
      include 'initialcon.com'
      include 'ocn_advec.com'
      include 'ocn_state.com'
      include 'ocn_paras.com'
      include 'ice_paras.com'
      include 'flx_paras.com'
      include 'flx_in.com'
      include 'output.com'
      include 'couple.com'
      
*     Outputs
      REAL U(NPTS,NZP1,NVEL),   ! On output contains
     $     X(NPTS,NZP1,NSCLR)   ! initial U,X fields.               
*
* Local Variablies, including some read in from name lists
      REAL Ubot(NVEL),Xbot(NSCLR) ! U,X at the bottom of the domain
      REAL dscale ! (neg)lambda parameter for defining the stretch
      REAL alat,alon,delta_lat,delta_lon                          

      INTEGER k,l,ipt,ix,iy

      CHARACTER*40 forcing_file
      
      NAMELIST/NAME_CONSTANTS/grav,vonk,sbc,twopi,onepi,TK0,spd,dpy,
     &     epsw,albocn,EL,SL,FL,FLSN
      NAMELIST/NAME_PROCSWIT/LKPP,LRI,LDD,LICE,
     &     LBIO,LNBFLX,LTGRID,LRHS,L_SSref
      NAMELIST/NAME_DOMAIN/DMAX,alon,alat,delta_lat,delta_lon,
     &     lstretchgrid,dscale
      NAMELIST/NAME_START/ L_INITDATA,initdata_file,L_INTERPINIT,
     &     L_RESTART
      NAMELIST/NAME_TIMES/ dtsec,startt,finalt,ndtocn 
      NAMELIST/NAME_ADVEC/ L_ADVECT,advect_file,L_RELAX,relax_in
      NAMELIST/NAME_PARAS/ paras_file,L_JERLOV
      NAMELIST/NAME_OUTPUT/ L_VAROUT,L_SINGOUT,L_AVES,output_file,
     &     ndtout,L_RESTARTW
      NAMELIST/NAME_FORCING/ L_FLUXDATA,forcing_file
      NAMELIST/NAME_COUPLE/ L_COUPLE,ifirst,ilast,jfirst,jlast,
     &     sstin_file

      OPEN(75,FILE='3D_ocn.nml')

* Initialse and read the constants name list
      spd=86400.                ! secs/day
      dpy=360.                  ! days/year
      twopi=8*atan(1.)          ! 2pi
      onepi=twopi/2             ! pi
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
      write(nuout,*) 'Read Namelist CONSTANTS'

* Initialize and read the processes namelist
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
      write(nuout,*) 'Read Namelist PROCSWIT'
      

* Initilalize and read the location name list
      DMAX=0.0
      alat=0.0
      alon=0.0
      delta_lat=2.5
      delta_lon=3.75
      dscale=0.0
      lstretchgrid=.FALSE.
      READ(75,NAME_DOMAIN)
      IF (DMAX .LE. 0.0) THEN 
         write(nuerr,*) 'You must specify a depth for the domain'
         CALL MIXED_ABORT
      ENDIF
      IF ((lstretchgrid) .AND. (dscale .EQ. 0.0)) THEN
         write(nuerr,*) "You can't have dscale=0 for stretched grids"
         CALL MIXED_ABORT
      ENDIF
      DO iy=1,ny
         DO ix=1,nx
            ipt=(iy-1)*nx+ix
            dlat(ipt)=alat+(iy-1)*delta_lat
            dlon(ipt)=alon+(ix-1)*delta_lon
         ENDDO
      ENDDO
      write(nuout,*) 'Read Namelist DOMAIN'

      CALL init_env(lstretchgrid,dscale)

* Initialize and read the start name list
      L_INITDATA= .TRUE.
      L_INTERPINIT= .TRUE.
      L_RESTART= .FALSE.
      READ(75,NAME_START) 
      write(nuout,*) 'Read Namelist START'

* Initialize and read the times namelist
      ndtocn=1
      dtsec=0.0
      startt=-999.999
      finalt=-999.999
      READ(75,NAME_TIMES) 
      IF ( (dtsec .LE. 0.0) .OR. (startt .LT. 0.0)
     &     .OR. (finalt .LT. 0.0)) THEN 
         write(nuerr,*) 'You must specify values of dtsec,startt,finalt'
         CALL MIXED_ABORT
      ENDIF
      startt=startt*spd
      finalt=finalt*spd
      dto=dtsec/float(ndtocn)
      nend=int((finalt-startt)/dtsec)
      nstart=nint(startt)/dto
      IF (float(nend*ndtocn) .NE. (finalt-startt)/dto ) THEN
         write(nuerr,*) 'The integration length is not '
     &        //'a multiple of the ocean timestep' 
         write(nuerr,*) 'dto=',dto
         write(nuerr,*) 'finalt=',finalt
         write(nuerr,*) 'startt=',startt
         CALL MIXED_ABORT
      ENDIF
      startt=startt/spd
      finalt=finalt/spd
      time=startt
      write(nuout,*) 'Read Namelist TIMES'

*Initialize and read the couple namelist
      L_COUPLE=.FALSE.
      ifirst=1
      ilast=nx
      jfirst=1
      jfirst=ny
      READ(75,NAME_COUPLE)
      write(nuout,*) 'Read Namelist COUPLE'
      IF (L_COUPLE) THEN
         call read_sstin
      ENDIF

* Initialize and read the advection namelist
      L_ADVECT=.FALSE.
      L_RELAX=.FALSE.
      DO iy=1,ny
        relax_in(iy)=0.0
      ENDDO
      READ(75,NAME_ADVEC)
      IF (L_ADVECT) THEN 
         CALL init_advect
      ELSE
         DO ipt=1,npts
            nmodeadv(ipt,1)=0
            nmodeadv(ipt,2)=0
         ENDDO
         write(nuout,*) 'No advection has been specified'
      ENDIF
      write(nuout,*) 'Read Namelist ADVEC'
      write(nuout,*) 'relax_in'
      IF (L_RELAX) THEN
         call init_relax
      ENDIF


* Initialize and read the paras namelist
      paras_file='3D_ocnparas.nc'
      L_JERLOV=.TRUE.
      READ(75,NAME_PARAS)
      call init_paras
      write(nuout,*) 'Read Namelist PARAS'

* Initialize and read the forcing namelist
      L_FLUXDATA=.FALSE.
      forcing_file='1D_ocean_forcing.nc'
      READ(75,NAME_FORCING)
      write(nuout,*) 'Read Namelist FORCING'


      IF (L_RESTART) THEN
         CALL READ_RESTART(U,X)
      ELSE
         CALL init_flds(U,X)
         write(nuout,*) 'Fields Initialized'
      ENDIF
      CALL init_flx
      IF ((L_FLUXDATA) .AND. .NOT. L_COUPLE) THEN
         CALL init_flxdata(forcing_file)
      ENDIF

* Initialize and read the output name list
      DO l=1,N_VAROUTS
         L_VAROUT(l)=.TRUE.
      ENDDO
      DO l=1,N_SINGOUTS
         L_SINGOUT(l)=.TRUE.
      ENDDO
      L_AVES=.FALSE.
      L_RESTARTW=.TRUE.
      ndtout=1
      output_file='KPPocean'
      READ(75,NAME_OUTPUT)
      write(nuout,*) 'Read Namelist OUTPUT'
      
      flen=INDEX(output_file,' ')-1
      day_out=int(startt)+5
      write(output_file(flen+1:flen+1),'(a)') '_'
      write(output_file(flen+2:flen+5),'(i4.4)'),day_out
      write(output_file(flen+6:flen+8),'(3A)') '.nc'

      dtout=ndtout*dto/spd
      call init_output

      
      call init_ocn(U,X)

* Write out the data from the initial condition

      IF ( .NOT. L_RESTART) THEN
         call output_inst(U,X)
      ENDIF
     
      CLOSE(75)

      RETURN
      END

      SUBROUTINE WRITE_RESTART(U,X)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#ifdef COUPLE
#ifdef CFS
      include 'parameter.cfs_coupled.inc'
#else
      include 'parameter.oasis2.inc'
#endif
#else
#ifdef CFS
      include 'parameter.cfs_forced.inc'
#else
      include 'parameter.forced.inc'
#endif
#endif
      include 'times.com'
      include 'ocn_paras.com'
      include 'ocn_state.com'
      include 'kprof_out.com'
      
*     Inputs
      REAL U(NPTS,NZP1,NVEL),   ! On output contains
     $     X(NPTS,NZP1,NSCLR)   ! initial U,X fields.               

c
c Local Common Blocks
c
      real hmixd(NPTS,0:1),     ! storage arrays for extrapolations
     +     Us(NPTS,NZP1,NVEL ,0:1), ! ..      ..     ..  ..
     +     Xs(NPTS,NZP1,NSCLR,0:1)  ! ..      ..     ..  ..
      integer old(NPTS),new(NPTS) ! extrapolation index for Us,Xs,hmixd
      common/ saveUXh /
     +     old,new,Us,Xs,hmixd
      

      WRITE(31) time,U,X,CP,rho,hmix,kmix,Sref,SSref,Ssurf,Tref,
     &     old,new,Us,Xs,hmixd
      CLOSE(31)

      RETURN
      END

      SUBROUTINE READ_RESTART(U,X)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#ifdef COUPLE
#ifdef CFS
      include 'parameter.cfs_coupled.inc'
#else
      include 'parameter.oasis2.inc'
#endif
#else
#ifdef CFS
      include 'parameter.cfs_forced.inc'
#else
      include 'parameter.forced.inc'
#endif
#endif
      include 'times.com'
      include 'ocn_paras.com'
      include 'ocn_state.com'
      include 'kprof_out.com'
      
*     Inputs
      REAL U(NPTS,NZP1,NVEL),   ! On output contains
     $     X(NPTS,NZP1,NSCLR)   ! initial U,X fields.               

c
c Local Common Blocks
c
      real hmixd(NPTS,0:1),     ! storage arrays for extrapolations
     +     Us(NPTS,NZP1,NVEL ,0:1), ! ..      ..     ..  ..
     +     Xs(NPTS,NZP1,NSCLR,0:1)  ! ..      ..     ..  ..
      integer old(NPTS),new(NPTS) ! extrapolation index for Us,Xs,hmixd
      common/ saveUXh /
     +     old,new,Us,Xs,hmixd
      

      READ(30) time,U,X,CP,rho,hmix,kmix,Sref,SSref,Ssurf,Tref,
     &     old,new,Us,Xs,hmixd
      CLOSE(30)

c     IF (abs(time-startt) .GT. 1.e-4) THEN 
c        write(nuerr,*) 'Start time doesn''t match the restart record'
c        CALL MIXED_ABORT
c     ENDIF

      RETURN
      END

      SUBROUTINE MIXED_ABORT

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

      include 'couple.com'

      IF (L_COUPLE) THEN
         call halte(' ')
      ELSE
         CALL MIXED_ABORT
      ENDIF

      END

      SUBROUTINE init_relax

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#ifdef COUPLE
#ifdef CFS
      include 'parameter.cfs_coupled.inc'
#else
      include 'parameter.oasis2.inc'
#endif
#else
#ifdef CFS
      include 'parameter.cfs_forced.inc'
#else
      include 'parameter.forced.inc'
#endif
#endif
      include 'constants.com'
      include 'ocn_advec.com'
      include 'couple.com'

      INTEGER ix,iy,ipoint

      REAL sst_in(NX_GLOBE,NY_GLOBE),ice_in(NX_GLOBE,NY_GLOBE)
      REAL usf_in(NX_GLOBE,NY_GLOBE),vsf_in(NX_GLOBE,NY_GLOBE)

      COMMON /save_sstin/ sst_in,ice_in,usf_in,vsf_in

      DO iy=1,ny
         IF (relax_in(iy) .EQ. 0.0) THEN
            DO ix=1,nx
               ipoint=(iy-1)*nx+ix
               relax(ipoint)=0.0
            ENDDO
         ELSE
            DO ix=1,nx
               ipoint=(iy-1)*nx+ix
               relax(ipoint)=1./(relax_in(iy)*spd)
            ENDDO
         ENDIF
      ENDDO

      DO iy=1,ny
         DO ix=1,nx
            ipoint=(iy-1)*nx+ix
            SST0(ipoint)=SST_in(ix+ifirst-1,iy+jfirst-1)
            fcorr(ipoint)=0.0
         ENDDO
      ENDDO
      
      RETURN
      END
