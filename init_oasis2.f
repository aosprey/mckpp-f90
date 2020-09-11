#ifdef COUPLE
#ifdef OASIS2
      SUBROUTINE mpi1_cpl_init(kastp,kexch,kstep,il_commlocal)
C
C     Initialize the MPI environment for the model.
C     NPK 13/11/07 and modified 15/11/07
C
C     Call the Oasis CLIM_Init.F routine to receive a local
C     MPI communicator and transmit timestep information.
C
C     Attach a (LARGE) MPI buffer for sending and receiving messages.
C
      INTEGER*4,PARAMETER :: nuout=6,nuerr=0,buffer_size=10000000
c     
#INCLUDE "clim.h"
#INCLUDE "mpiclim.h"
#INCLUDE "oasis.h"
#INCLUDE "param.h"
#INCLUDE "param_cou.h"
#INCLUDE "inc_cpl.h"
#INCLUDE <mpif.h>
c
      INTEGER il_commlocal,kastp,kstep,kexch,inmods,ierr
      INTEGER il_size,il_rank
      CHARACTER*6 clmodnam
      CHARACTER*3 cljobnam
      LOGICAL llflag
      REAL*8 send_buffer(buffer_size)
c
      WRITE(clmodnam,'(A6)') 'KPPoce'
      WRITE(cljobnam,'(A3)') 'BRU'
c
      inmods=2
c
      WRITE(nuout,*) 'KPP:Now calling CLIM_Init'
      WRITE(nuout,*) 'KPP:Job name is ',cljobnam
      WRITE(nuout,*) 'KPP:Model name is ',clmodnam
      WRITE(nuout,*) 'KPP:Number of models is ',inmods
      CALL CLIM_Init(cljobnam,clmodnam,inmods,7,kastp,kexch,kstep,5,3600
     &     ,3600,il_commlocal)
      WRITE(nuout,*) 'KPP:Finished call to CLIM_Init'
c
c     il_commlocal is the local communicator, not an error code.
c     I am not sure why it is being used as an error code, particularly
c     since it can easily be negative. Comment out this test.
c
c     NPK 12/12/08
c
c      IF (il_commlocal .lt. CLIM_Ok) THEN
c         WRITE(nuout,*) 'KPP:inicmo : pb init clim'
c         CALL halte('STOP in inicmo')
c         WRITE(nuout,*) 'KPP:error code is = ', il_commlocal
c      ELSE
         WRITE(nuout,*) 'KPP:inicmo : init clim ok'
c      ENDIF
c
      WRITE(nuout,*)'KPP:Now calling MPI_BUFFER_ATTACH'
C
C     Must attach a buffer of buffer_size*8 for REAL*8 variables
C
      CALL MPI_BUFFER_ATTACH(send_buffer,buffer_size*8,ierr)
      WRITE(nuout,*)'KPP:Done calling MPI_BUFFER_ATTACH,',
     &     'buffer_size=',buffer_size
c
      RETURN
      END

      SUBROUTINE mpi1_cpl_stepi(cloasis,istep, ifcpl, idt)
c
      IMPLICIT NONE
c
c     Wrapper routine for CLIM_Stepi when using MPI1 interface
c     NPK 16/11/07
c
      INTEGER istep,ifcpl,idt,info
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
      CHARACTER*5 cloasis
c
      CALL CLIM_Stepi(cloasis,istep,ifcpl,idt,info)
c
c     Note that I have had to change the test on info here.  It used
c     to be:
c
c     IF (info .NE. clim_ok)
c
c     which is incorrect, as info is not an output code: it is the number
c     of processors involved in the coupling.  Therefore I am changing
c     this to be a test on info=0, which would mean that no processors
c     matched the name of the model passed to CLIM_Stepi (i.e., cloasis)
c
      IF (info .LE. 0) THEN
         WRITE ( UNIT = nuerr, FMT = *)
     $        ' warning : problem in getting step info ',
     $        'from oasis '
         WRITE (UNIT = nuerr, FMT = *)
     $        ' =======   error code number = ', info
      ELSE
         WRITE (UNIT = nuout, FMT = *)
     $        ' got step information from oasis '
         WRITE ( nuout, *) ' number of tstep in oasis ', istep
         WRITE ( nuout, *) ' exchange frequency in oasis ', ifcpl
         WRITE ( nuout, *) ' length of tstep in oasis ', idt           
      ENDIF
c
      RETURN
      END

c
c**** *inicmo*  - initialize coupled mode communication for ocean
c                 and exchange some initial information with Oasis
c
c     input:
c     -----
c       kitro  : total number of timesteps in oceanic model
c       kexco  : frequency of exchange for the fluxes (in time steps)
c       kstpo  : timestep value (in seconds)
c
c     -----------------------------------------------------------
c
      SUBROUTINE inicmo(kastp,kexch,kstep)
c
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
#INCLUDE "param.h"
c
      INTEGER kastp, kexch, kstep
      INTEGER iparal(3)
      INTEGER ifcpl, idt, info, imxtag, istep, jf, il_commlocal
c
#INCLUDE "param_cou.h"
#INCLUDE "inc_cpl.h"
      CHARACTER*3 cljobnam      ! experiment name
      CHARACTER*6 clmodnam      ! model name
      CHARACTER*5 cloasis       ! coupler name (Oasis)
      INTEGER imess(4)
      INTEGER getpid            ! system functions
c
#INCLUDE "clim.h"
#INCLUDE "mpiclim.h"
c
#INCLUDE "oasis.h"         ! contains the name of communication technique.
                                ! Here cchan=CLIM only is possible.
c
C     -----------------------------------------------------------
C
C*    1. Initializations
C        ---------------
C      
      WRITE(nuout,*) ' '
      WRITE(nuout,*) ' '
      WRITE(nuout,*) ' ROUTINE INICMO'
      WRITE(nuout,*) ' **************'
      WRITE(nuout,*) ' '
      WRITE(nuout,*) ' '
c
c     Define the model name
c
      clmodnam = 'KPPoce'       ! as $NBMODEL in namcouple
c
c     Define the coupler name
c
      cloasis = 'Oasis'        !  as in coupler
c
c
c     Define symbolic name for fields exchanged from ocean to coupler,
c         must be the same as (1) of the field  definition in namcouple:
c
      cl_writ(2)='SOSSTSST'
      cl_writ(1)='SOICECOV'
      cl_writ(3)='SOZOCRTX'
      cl_writ(4)='SOMECRTY'
c
c     Define files name for fields exchanged from ocean to coupler,
c         must be the same as (6) of the field  definition in namcouple:
c
      cl_f_writ(1)='sstocean'
      cl_f_writ(2)='sstocean'
      cl_f_writ(3)='sstocean'
      cl_f_writ(4)='sstocean'
c
c     Define symbolic name for fields exchanged from coupler to ocean
c         must be the same as (2) of the field  definition in namcouple:
c
      cl_read(1)='SONSHLDO'     ! non solar heat flux (positif vers l'ocean)
      cl_read(2)='SOSHFLDO'     ! solar flux
      cl_read(3)='SOWAFLDO'     ! wter flux
      cl_read(4)='SOZOTAUX'     ! first zonal wind stress
      cl_read(5)='SOMETAUY'     ! first meridien wind stress
c
c     Define files name for fields exchanged from coupler to ocean
c         must be the same as (7) of the field  definition in namcouple:
c
      cl_f_read(1)='flxocean'
      cl_f_read(2)='flxocean'
      cl_f_read(3)='flxocean'
      cl_f_read(4)='flxocean'
      cl_f_read(5)='flxocean'
c
c     Define the number of processors involved in the coupling for
c     Oasis (=1) and each model (as last two INTEGER on $CHATYPE line
c     in the namcouple); they will be stored in a COMMON in mpiclim.h
c     (used for CLIM/MPI2 only)
      mpi_nproc(0)=1
      mpi_nproc(1)=1
      mpi_nproc(2)=1 
c
c
c     Define infos for sending to oasis
c
      imess(1) = kastp
      imess(2) = kexch
      imess(3) = kstep
      imess(4) = getpid()
c
c     Initialization and exchange of initial info in the CLIM technique
c
      IF (cchan.eq.'CLIM') THEN
c  
c         Define the experiment name :
c
          cljobnam = 'BRU'      ! as $JOBNAM in namcouple
c
c         Define the number of processors used by each model as in
c         $CHATYPE line of namcouple (used for CLIM/MPI2 only)
          mpi_totproc(1)=1
          mpi_totproc(2)=1
c
c         Define names of each model as in $NBMODEL line of namcouple
c         (used for CLIM/MPI2 only)
          cmpi_modnam(1)='toyatm'
          cmpi_modnam(2)='mixed'
c
c         Start the coupling 
          CALL MPI1_CPL_INIT(kastp,kexch,kstep,il_commlocal)
c         (see lib/clim/src/CLIM_Init for the definition of input parameters)
c
c          WRITE(nuout,*) ' inicmo: job name is ',cljobnam
c          WRITE(nuout,*) ' inicmo: model name is ',clmodnam
c          WRITE(nuout,*) ' inicmo: now calling CLIM_Init, kastp=',kastp
c          CALL CLIM_Init(cljobnam,clmodnam,2,7,kastp,kexch,kstep,
c     *                 5, 180, 180, info )
c
c          WRITE(nuout,*) ' inicmo: passed called to CLIM_Init'
c          IF (info.ne.clim_ok) THEN
c              WRITE ( nuerr, *) ' inicmo : pb init clim '
c              WRITE ( nuerr, *) ' error code is = ', info
c              CALL halte('STOP in inicmo')
c            ELSE
c              WRITE(nuout,*) 'inicmo : init clim ok '
c          ENDIF
c
c         For each coupling field, association of a port to its symbolic name
c
c         -Define the parallel decomposition associated to the port of each
c          field; here no decomposition for all ports.

          iparal ( clim_strategy ) = clim_serial
          iparal ( clim_length   ) = imjm
          iparal ( clim_offset   ) = 0
c
c         -Loop on total number of coupler-to-ocean fields
c         (see lib/clim/src/CLIM_Define for the definition of input parameters)

          DO jf=1, jpflda2o1+jpflda2o2
            CALL CLIM_Define (cl_read(jf), clim_in , clim_double, iparal
     $          , info )  
          END DO 
c
c         -Loop on total number of ocean-to-coupler fields 
c         (see lib/clim/src/CLIM_Define for the definition of input parameters)
          DO jf=1, jpfldo2a
            CALL CLIM_Define (cl_writ(jf), clim_out , clim_double,
     $          iparal, info )   
          END DO 
          WRITE(nuout,*) 'inicmo : clim_define ok '
c
c         -Join a pvm group, wait for other programs and broadcast usefull 
c          informations to Oasis and to the atmosphere 
c          (see lib/clim/src/CLIM_Start)

          CALL CLIM_Start ( imxtag, info )
          IF (info.ne.clim_ok) THEN
              WRITE ( nuerr, *) 'inicmo : pb start clim '
              WRITE ( nuerr, *) ' error code is = ', info
              CALL halte('stop in inicmo')
            ELSE
              WRITE ( nuout, *)  'inicmo : start clim ok '
          ENDIF
c
c         -Get initial information from Oasis
c          (see lib/clim/src/CLIM_Stepi)
          CALL MPI1_CPL_STEPI(cloasis,istep,ifcpl,idt)
c
c          CALL CLIM_Stepi (cloasis, istep, ifcpl, idt, info)
c          IF (info .NE. CLIM_Ok) THEN
c              WRITE ( UNIT = nuerr, FMT = *)
c     $            ' warning : problem in getting step info ',
c     $            'from oasis '
c              WRITE (UNIT = nuerr, FMT = *)
c     $            ' =======   error code number = ', info
c            ELSE
c              WRITE (UNIT = nuout, FMT = *)
c     $            ' got step information from oasis '
c          ENDIF
c          WRITE ( nuout, *) ' number of tstep in oasis ', istep
c          WRITE ( nuout, *) ' exchange frequency in oasis ', ifcpl
c          WRITE ( nuout, *) ' length of tstep in oasis ', idt
      ENDIF
c
      RETURN
      END
#endif /*OASIS2*/
#endif /*COUPLE*/
