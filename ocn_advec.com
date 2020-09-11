c      INTEGER nmodeadv(npts,2),modeadv(npts,maxmodeadv,2)
c      REAL advection(npts,maxmodeadv,2)
      LOGICAL L_ADVECT
      CHARACTER*200 advect_file
      INTEGER ncid_advec
      common/ ocn advec / 
     &     L_ADVECT,advect_file,ncid_advec

      LOGICAL L_RELAX_SST,L_RELAX_CALCONLY,L_RELAX_INIT
c      REAL relax_sst(npts),SST0(npts),fcorr(npts)
      REAL relax_sst_in(ny)

      common/ ocn_relax / relax_sst_in, 
     &   L_RELAX_SST,L_RELAX_CALCONLY,L_RELAX_INIT
