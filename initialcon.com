c      REAL hu(NPTS),SSU(NPTS,NVEL),hx(NPTS),SSX(NPTS,NSCLR)
c      INTEGER mu(NPTS),mx(NPTS)
c      CHARACTER*40 initprofile
c      common/ initialcon/ hu,SSU,hx,mu,mx,SSX,
c     +                    initprofile

      LOGICAL L_INITDATA,L_INTERPINIT,L_RESTART,L_PERSIST_SST,
     +   L_PERSIST_SST_ANOM, L_PERSIST_ICE, L_PERSIST_ICE_ANOM
      CHARACTER*200 initdata_file
      CHARACTER*200 restart_infile
      common/ initdata / L_INITDATA,L_INTERPINIT,L_RESTART,
     +   L_PERSIST_SST,L_PERSIST_SST_ANOM,initdata_file,restart_infile,
     +   L_PERSIST_ICE,L_PERSIST_ICE_ANOM
      
