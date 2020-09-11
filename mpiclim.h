C
C -- mpiclim.h  26-10-99   Version 2.4   Author: Jean Latour (F.S.E.)
C    *********
C@
C@  Contents : variables related to MPI-2 message passing
C@  --------
C@
C@ -- mpi_totproc: number of processors on which to launch each model
C@
C@ -- mpi_nproc: number of processors involved in the coupling for
C@               each model
C@ -- cmpi_modnam: models name
C     -----------------------------------------------------------------
C
      INTEGER*4 mpi_totproc(1:CLIM_MaxMod-1),mpi_nproc(0:CLIM_MaxMod-1)
C
      CHARACTER*6 cmpi_modnam(1:CLIM_MaxMod-1)
C
      common/CLIM_mpiclim/mpi_totproc, mpi_nproc, cmpi_modnam 
C
C     -----------------------------------------------------------------
