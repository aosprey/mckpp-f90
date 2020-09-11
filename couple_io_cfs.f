      SUBROUTINE READ_GFS_FORCING(solar,non_solar,PminusE,
     +     u_stress,v_stress)
c
c     SUBROUTINE read_gfs_forcing supports coupling to the NCEP GFS.
c     This subroutine reads the flux data output by the GFS
c     and converts it into the format expected by KPP.
c
c     Written by Nicholas Klingaman, 19-23 March 2009.
c     Updated by Nicholas Klingaman, 25 May 2009.
c
c     This procedure currently supports coupling only once per day.
c     The GFS outputs two surface-forcing files per day, called
c     flxcplgrb12 and flxcplgrb24.  This subroutine takes the mean
c     of the fields in those files and passes the results to KPP.
c     
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
      
#ifdef COUPLE
#ifdef CFS
#include "parameter.cfs_coupled.inc"
#else
#include "parameter.oasis2.inc"
#endif
#else
#ifdef CFS
#include "parameter.cfs_forced.inc"
#else
#include "parameter.forced.inc"
#endif
#endif
#include "times.com"
#include "couple.com"
#include "location.com"

      REAL solar(NPTS),non_solar(NPTS),PminusE(NPTS),
     +     u_stress(NPTS),v_stress(NPTS),precip(NPTS),evap(NPTS),
     +     efcnv
      REAL*4 work_twod(NX_GLOBE,NY_GLOBE),work_oned(NPTS)

      INTEGER i,ierr,atmos_unit
      PARAMETER(atmos_unit=11)
      CHARACTER*11 filename

c     Zero out arrays to initialize.
      solar=0.
      non_solar=0.
      PminusE=0.
      u_stress=0.
      v_stress=0.
      evap=0.
      precip=0.

c     Constant to convert latent heat to evaporation
      efcnv=1./(2.5E6)

c     Read in the files from the atmospheric model and create
c     a daily mean (by taking mean of two 12-hourly flux files)

      DO i=1,2
         IF (i .eq. 1) THEN
            WRITE(filename,'(A11)') 'flxcplgrb12'
         ELSEIF (i .eq. 2) THEN
            WRITE(filename,'(A11)') 'flxcplgrb24'
         ENDIF
         WRITE(nuout,*) 'Now attempting to open file ',filename
         CALL baopenr(atmos_unit,filename,ierr)
         IF (ierr .ne. 0) THEN
            WRITE(nuout,*) 'Read_GFS_Forcing: Error opening file ',
     +           filename,' ierr=',ierr
            CALL mixed_abort
         ELSE
            WRITE(nuout,*) 'Read_GFS_Forcing: Successfully opened file '
     +           ,filename
         ENDIF
         
c     Note that KPP has the following sign conventions:
c     Solar flux: positive downward (positive for warming)
c     Non-solar flux: positive downward (positive for warming)

c     Get zonal momentum flux (N/m^2)
         WRITE(6,*) 'Read_GFS_Forcing: Getting zonal momentum flux ...'
         CALL CFS_GETFLX(atmos_unit,124,1,0,work_twod,NX_GLOBE,NY_GLOBE)         
         CALL twod_global_oned_regional(work_twod,work_oned)
         u_stress=u_stress+work_oned*0.5
c     Get meridional momentum flux (N/m^2)
         WRITE(6,*) 'Read_GFS_Forcing: Getting meridional momentum ',
     +        'flux ...'
         CALL CFS_GETFLX(atmos_unit,125,1,0,work_twod,NX_GLOBE,NY_GLOBE)
         CALL twod_global_oned_regional(work_twod,work_oned)
         v_stress=v_stress+work_oned*0.5
c     Get downward solar radiation (W/m^2)
         WRITE(6,*) 'Read_GFS_Forcing: Getting downward solar ',
     +        'radiation ...'
         CALL CFS_GETFLX(atmos_unit,204,1,0,work_twod,NX_GLOBE,NY_GLOBE)
         CALL twod_global_oned_regional(work_twod,work_oned)
         solar=solar+work_oned*0.5
c     Get upward solar radiation (W/m^2)
         WRITE(6,*) 'Read_GFS_Forcing: Getting upward solar ',
     +        'radiation ...'
         CALL CFS_GETFLX(atmos_unit,211,1,0,work_twod,NX_GLOBE,NY_GLOBE)
         CALL twod_global_oned_regional(work_twod,work_oned)
         solar=solar-work_oned*0.5
c     Get downward longwave flux (W/m^2) <- Note sign convention as above
         WRITE(6,*) 'Read_GFS_Forcing: Getting downward longwave ',
     +        'flux ...'
         CALL CFS_GETFLX(atmos_unit,205,1,0,work_twod,NX_GLOBE,NY_GLOBE)
         CALL twod_global_oned_regional(work_twod,work_oned)
         non_solar=non_solar+work_oned*0.5
c     Get upward longwave flux (W/m^2) <- Note sign convention as above
         WRITE(6,*) 'Read_GFS_Forcing: Getting upward longwave ',
     +        'flux ...'
         CALL CFS_GETFLX(atmos_unit,212,1,0,work_twod,NX_GLOBE,NY_GLOBE)
         CALL twod_global_oned_regional(work_twod,work_oned)
         non_solar=non_solar-work_oned*0.5
c     Get precipitation (mm/s)
         WRITE(6,*) 'Read_GFS_Forcing: Getting precipitation ...'
         CALL CFS_GETFLX(atmos_unit,59,1,0,work_twod,NX_GLOBE,NY_GLOBE)
         CALL twod_global_oned_regional(work_twod,work_oned)
         precip=precip+work_oned*0.5
c     Get latent heat flux (positive upward) <- Note sign convention as above
         WRITE(6,*) 'Read_GFS_Forcing: Getting latent heat flux ...'
         CALL CFS_GETFLX(atmos_unit,121,1,0,work_twod,NX_GLOBE,NY_GLOBE)
         CALL twod_global_oned_regional(work_twod,work_oned)
         non_solar=non_solar-work_oned*0.5
         evap=evap+work_oned*0.5*efcnv
c     Take precipitation minus evaporation for return to KPP
         PminusE=PminusE+precip-evap         
c     Get sensible heat flux (positive upward) <- Note sign convention as above
         WRITE(6,*) 'Read_GFS_Forcing: Getting sensible heat flux ...'
         CALL CFS_GETFLX(atmos_unit,122,1,0,work_twod,NX_GLOBE,NY_GLOBE)
         CALL twod_global_oned_regional(work_twod,work_oned)
         non_solar=non_solar-work_oned*0.5     

         WRITE(6,*) 'Read_GFS_Forcing: Now closing input file ',
     +        filename,' ...'
         CALL baclose(atmos_unit,ierr)
      ENDDO               

      RETURN
      END

      SUBROUTINE CFS_WRITE_GRIB_SSTS(kpp_3d_fields)
c 
c     SUBROUTINE write_grib_ssts supports coupling to the NCEP GFS. 
c     This subroutine blends the KPP SSTs with the climatological
c     or observed SSTs prescribed via the L_CLIMSST and sstin_file
c     namelist options.  The resulting SST field is written to 
c     a GRIB file that can be read by the GFS.
c
c     Note that this subroutine writes the *global* SST field
c     to the GRIB file.
c
c     The blending can be controlled via the netCDF file of coupling
c     weights (use the L_CPLWGHT and cplwght_file options).
c
c     Written by Nicholas Klingaman 22-25 May 2009.
c
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER(nuout=6,nuerr=0)
     
#include "kpp_3d_type.com"
#include "location.com"
#include "couple.com"
#include "constants.com"
#include "landsea.cfs_coupled.com"

c     kpp_sst is a one-dimensional array of the model SST on the model grid.
      REAL SST_in(NX_GLOBE,NY_GLOBE,1)
c     Local variables
c     1X1 degree grids (output)
      INTEGER NX_GLOBE_1X1, NY_GLOBE_1X1
      PARAMETER(NX_GLOBE_1X1=360,NY_GLOBE_1X1=180)
      REAL*4 global_longitudes_1X1(NX_GLOBE_1X1),
     +     global_latitudes_1X1(NY_GLOBE_1X1),
     +     SST_out(NX_GLOBE,NY_GLOBE),
     +     SST_out_1X1(NX_GLOBE_1X1,NY_GLOBE_1X1),
     +     global_latitudes_StoN(NY_GLOBE),
     +     half_globe_temp(NY_GLOBE/2)
c     iyrcur, idatecur and imocur are the variables read in
c     from the "sstparm" namelist file - do not adjust these names.
      INTEGER iyrcur, idaycur, imocur, current_year_of_century, 
     +     next_year_of_century,century,iyrnxt, idaynxt, imonxt
c     kpds and kgds contain parameters for the GRIB file metadata.
      INTEGER kpds(200), kgds(200),i,ix,j,jy,ipoint,ipoint_globe,error
c     Name for GRIB file
      CHARACTER(LEN=20) kpp_sst_grib_anal,kpp_sst_grib_fcst
      DATA kpp_sst_grib_anal/'kpp_sst_00000000.grb'/
      DATA kpp_sst_grib_fcst/'kpp_sst_00000000.grb'/
      LOGICAL*1 grib_logical(NX_GLOBE_1X1,NY_GLOBE_1X1)
      REAL*4 bad
      PARAMETER(bad=-999.9)

      COMMON /save_sstin/ SST_in
      
      NAMELIST/datecur/iyrcur,idaycur,imocur
      NAMELIST/datenxt/iyrnxt,idaynxt,imonxt

c     Information about the current and next date, month and year in the GFS
c     is contained in the file "sstparm".  KPP does not require this
c     information, but we do need it for the GRIB output file.
      OPEN(22,FILE='sstparm')
c     Information about the current date, month and year is in the
c     "datecur" namelist.
      READ(22,datecur)
c     Information about "tomorrow's" date, month and year is in the
c     "datenxt" namelist.
      READ(22,datenxt)
      CLOSE(22)
      
      current_year_of_century=MOD(iyrcur-1,100)+1
      next_year_of_century=MOD(iyrnxt-1,100)+1
      century=(iyrcur+99)/100

c     Set up 1X1 grids
      DO i=1,NX_GLOBE_1X1
         global_longitudes_1X1(i)=i-0.5
      ENDDO
      DO j=1,NY_GLOBE_1X1
         global_latitudes_1X1(j)=j-90.5
      ENDDO

      DO ix=1,NX_GLOBE
         DO jy=1,NY_GLOBE
            ipoint_globe = (jy-1)*NX_GLOBE+ix
            IF (kpp_3d_fields%cplwght(ipoint_globe) .LT. -1e-10) THEN              
c     Point is outside the coupling domain; set to SST climatology
               SST_out(ix,jy) = SST_in(ix,jy,1)
            ELSE
c     Point is inside the coupling domain; set to weighted value
               ipoint=(jy-jfirst)*NX+(ix-ifirst)+1
               SST_out(ix,jy) = kpp_3d_fields(ipoint,1,1)*
     &              kpp_3d_fields%cplwght(ipoint_globe)+SST_in(ix,jy,1)
     &              *(1-kpp_3d_fields%cplwght(ipoint_globe))
            ENDIF            
         ENDDO
      ENDDO
      IF (L_OUTKELVIN) SST_out = SST_out+TK0

      WRITE(nuout,*) 'Write_GRIB_SSTs: Produced global SST field'
      WRITE(nuout,*) 'Write_GRIB_SSTs: ... at KPP resolution'

c      WRITE(nuout,*) 'global_lsm=',global_lsm

c     Fill land points in KPP SST field with missing value (bad).
      DO j=1,NY_GLOBE
         DO i=1,NX_GLOBE
            IF (global_lsm(i,j) .ge. 0.5) SST_out(i,j)=bad
         ENDDO
      ENDDO

c     Fill land points with adjacent ocean values
      CALL CFS_FILL_LAND(SST_out,NX_GLOBE,NY_GLOBE,bad)

c     Interpolation routines require data to be in South->North,
c     but KPP is North->South.  Reverse latitude dimension.
      CALL CFS_REVERSE_LATITUDE_INPLACE(SST_out,global_latitudes,
     +     NX_GLOBE,NY_GLOBE)
      
c     Interpolate KPP SSTs to 1x1 grid required for GFS
      CALL CFS_KPP_TO_1X1(SST_out,NX_GLOBE,NY_GLOBE,
     +     global_longitudes,global_latitudes,SST_out_1X1,
     +     NX_GLOBE_1X1,NY_GLOBE_1X1,global_longitudes_1X1,
     +     global_latitudes_1X1,bad)

      WRITE(nuout,*) 'Write_GRIB_SSTs: ... and at 1x1 deg resolution'

c     Output must be in North->South, but interpolated SSTs are
c     South->North. Reverse latitude dimension (again).
      CALL CFS_REVERSE_LATITUDE_INPLACE(SST_out_1X1,
     +     global_latitudes_1X1,NX_GLOBE_1X1,NY_GLOBE_1X1)

c     Construct output filenames for analysis file ...
      WRITE(kpp_sst_grib_anal(9:12),'(i4.4)') iyrcur
      WRITE(kpp_sst_grib_anal(13:14),'(i2.2)') imocur
      WRITE(kpp_sst_grib_anal(15:16),'(i2.2)') idaycur
c     ... and forecast file.
      WRITE(kpp_sst_grib_fcst(9:12),'(i4.4)') iyrnxt
      WRITE(kpp_sst_grib_fcst(13:14),'(i2.2)') imonxt
      WRITE(kpp_sst_grib_fcst(15:16),'(i2.2)') idaynxt

c     Write SSTs to the analysis (today) file.
      WRITE(6,*) 'Write_GRIB_SSTs: Now outputting to file ',
     +     kpp_sst_grib_fcst  
      CALL BAOPENWT(60,kpp_sst_grib_anal,error)
      IF (error.ne.0) THEN
         WRITE(nuerr,*) 'Write_GRIB_SSTs: Error opening the GRIB SST ',
     +        'output file ',kpp_sst_grib_anal
         CALL MIXED_ABORT
      ENDIF

c     Fill the KPDS and KGDS arrays for the analysis SST file
      CALL CFS_FILL_KPDS_KGDS(global_longitudes_1X1,NX_GLOBE_1X1,
     +     global_latitudes_1X1,NY_GLOBE_1X1,century,
     +     current_year_of_century,imocur,idaycur,0,kpds,kgds)
      WRITE(6,*) 'Write_GRIB_SSTs: Filled the KPDS and KGDS arrays ',
     +     'for the analysis file'
c
c     NOTE: Do *not* substitute NPTS_GLOBE for NX_GLOBE*NY_GLOBE here!
c     NPTS_GLOBE is declared as INTEGER*8 and using it here will make
c     PUTGB either hang or crash with a segmentation fault.
c     NPK 25/05/09.
c
      CALL PUTGB(60,NX_GLOBE_1X1*NY_GLOBE_1X1,kpds,kgds,grib_logical,
     +     SST_out_1X1,error)
      IF (error.ne.0) THEN
         WRITE(nuerr,*) 'Write_GRIB_SSTs: Error writing to the GRIB SST'
     +        ,' output file ',kpp_sst_grib_anal
         CALL MIXED_ABORT
      ENDIF
      CALL BACLOSE(60,error)
      IF (error.ne.0) THEN
         WRITE(nuerr,*) 'Write_GRIB_SSTs: Error closing the GRIB SST ',
     +        'output file ',kpp_sst_grib_anal
         CALL MIXED_ABORT
      ENDIF

c     Now write the same SSTs to the forecast (tomorrow) file.
      WRITE(6,*) 'Write_GRIB_SSTs: Now outputting to file ',
     +     kpp_sst_grib_fcst
      CALL BAOPENWT(61,kpp_sst_grib_fcst,error)
      IF (error.ne.0) THEN
         WRITE(nuerr,*) 'Write_GRIB_SSTs: Error opening the GRIB SST ',
     +        'output file ',kpp_sst_grib_fcst
         CALL MIXED_ABORT
      ENDIF

c     Fill the KPDS and KGDS arrays for the forecast SST file
      CALL CFS_FILL_KPDS_KGDS(global_longitudes_1X1,NX_GLOBE_1X1,
     +     global_latitudes_1X1,NY_GLOBE_1X1,century,
     +     next_year_of_century,imonxt,idaynxt,0,kpds,kgds)
      WRITE(6,*) 'Write_GRIB_SSTs: Filled the KPDS and KGDS arrays ',
     +     'for the forecast GRIB file'

      
      CALL PUTGB(61,NX_GLOBE_1X1*NY_GLOBE_1X1,kpds,kgds,grib_logical,
     +     SST_out_1X1,error)
      IF (error.ne. 0) THEN
         WRITE(nuerr,*) 'Write_GRIB_SSTs: Error writing to the GRIB SST'
     +        ,' output file ',kpp_sst_grib_fcst
         CALL MIXED_ABORT
      ENDIF
        
      CALL BACLOSE(61,error)
      IF (error.ne.0) THEN
         WRITE(nuerr,*) 'Write_GRIB_SSTs: Error closing the GRIB SST ',
     +        'output file ',kpp_sst_grib_fcst
         CALL MIXED_ABORT
      ENDIF

      RETURN
      END

      SUBROUTINE twod_global_oned_regional(twod_global,oned_regional)
      IMPLICIT NONE
c
c     Transforms a two-dimensional global array (e.g., fluxes read
c     in from the atmospheric model) to a one-dimensional regional array
c     (e.g., those same fluxes on the domain of KPP)
c
c     Written by Nicholas Klingaman 19-23 March 2009.
c
#ifdef COUPLE
#ifdef CFS
#include "parameter.cfs_coupled.inc"
#else
#include "parameter.oasis2.inc"
#endif
#else
#ifdef CFS
#include "parameter.cfs_forced.inc"
#else
#include "parameter.forced.inc"
#endif
#endif
#include "couple.com"

      REAL*4 twod_global(NX_GLOBE,NY_GLOBE)
      REAL*4 oned_regional(NPTS)
      INTEGER i,j,ipoint
      
      DO j=jfirst,jlast
         DO i=ifirst,ilast
            ipoint=(j-jfirst)*nx+(i-ifirst)+1
            oned_regional(ipoint)=twod_global(i,j)
         ENDDO
      ENDDO

      END

c     SUBROUTINE getflx - copied from CFS coupA2O.f (coupler routine)
c     NPK made modifications as noted.
      SUBROUTINE CFS_GETFLX(nfi,kpds5,kpds6,kpds7,wrkT62,imT62,jmT62)
      
      integer JPDS(200),JGDS(200)
      integer KPDS(200),KGDS(200)
      integer IJDIM62,IRET
      real*4  tmp(imT62,jmT62)
      real*4  wrkT62(imT62,jmT62)
      logical*1 lbms(imT62,jmT62)

      IJDIM62=imT62*jmT62
      JPDS(:)=-1
      JPDS(5)=kpds5
      JPDS(6)=kpds6
      JPDS(7)=kpds7
      JGDS(:)=-1
      N=-1
      IRET=0
      DO i=1,200
         kpds(i)=0
         kgds(i)=0
      ENDDO

      CALL GETGB(nfi,0,IJDIM62,N,JPDS,JGDS,
     *     NDATA,KSKP,KPDS,KGDS,LBMS,tmp,IRET)
c     NPK added the following error trap
      IF (IRET .NE. 0) THEN
         WRITE(6,*) 'GETGB returned an error code of IRET=',IRET
      ENDIF

      lat1=kgds(4)              ! beginning latitude in grib

c     NPK has changed this from lat1 .gt. 0 to lat .lt. 0, because
c     the KPP grid runs from North to South.  Thus, we need reverse
c     the latitude dimension only if the grid is running from South
c     to North.
      if(lat1.lt.0) then
         do j=jmT62,1,-1
            jj=jmT62-j+1
            do i=1,imT62
               wrkT62(i,j)=tmp(i,jj)
            enddo
         enddo
      else
         wrkT62=tmp
      endif
      
      return
      end

c     SUBROUTINE FillLand - copied from CFS coupO2A.f (coupler routine)
c     Not modified
      SUBROUTINE CFS_FILL_LAND(f1X1,im1X1,jm1X1,bad)
      real*4 f1X1(im1X1,jm1X1)
      real*4 zf1X1(jm1X1)
      real*4 f1X1tmp1(im1X1,jm1X1)
      real*4 f1X1tmp2(im1X1,jm1X1)

      do j=1,jm1X1
      do i=1,im1X1
       f1X1tmp1(i,j)=f1X1(i,j)
      enddo
      enddo

      do nloop=1,5  ! loops to fill bad points with adjacent values
       do j=1,jm1X1
       do i=1,im1X1
        f1X1tmp2(i,j)=f1X1tmp1(i,j)
       enddo
       enddo
       do j=1,jm1X1
       do i=1,im1X1
        if(f1X1tmp1(i,j).eq.bad) then
         fe=1.0
         fw=1.0
         fn=1.0
         fs=1.0
         iw=i-1
         ie=i+1
         js=j-1
         jn=j+1
         if(iw.lt.1) iw=im1X1
         if(ie.gt.im1X1) iw=1
         if(js.lt.1) js=1
         if(jn.gt.jm1X1) jn=jm1X1
         if(f1X1tmp1(ie,j).eq.bad) fe=0.0
         if(f1X1tmp1(iw,j).eq.bad) fw=0.0
         if(f1X1tmp1(i,js).eq.bad) fs=0.0
         if(f1X1tmp1(i,jn).eq.bad) fn=0.0
         ftot=fe+fw+fs+fn
         if(ftot.gt.0.5) then
          f1X1tmp1(i,j)=(f1X1tmp1(ie,j)*fe+f1X1tmp1(iw,j)*fw+
     1                   f1X1tmp1(i,js)*fs+f1X1tmp1(i,jn)*fn)/ftot
         endif
        endif
       enddo
       enddo
      enddo

c fill bad ponits with zonal average
c     do j=1,jm1X1
c      zavg=0.0
c      c=0.0
c      do i=1,im1X1
c       if(f1X1tmp1(i,j).ne.bad) then
c        zavg=zavg+f1X1tmp1(i,j)
c        c=c+1.0
c       endif
c      enddo
c      if(c.gt.3.0) then
c       zavg=zavg/c
c       do i=1,im1X1
c        if(f1X1tmp1(i,j).eq.bad) f1X1tmp1(i,j)=zavg
c       enddo
c      endif
c     enddo

      do j=1,jm1X1
      do i=1,im1X1
       f1X1(i,j)=f1X1tmp1(i,j)
      enddo
      enddo

      return
      end
c
c     SUBROUTINE SSTMOM3to1x1 - copied from CFS coupO2A.f (coupler routine)
c     NPK modifications as noted
c
      subroutine CFS_KPP_TO_1X1(datin,INMAX,JNMAX,XIN,YIN
     1              ,datout,IOMAX,JOMAX,XOUT,YOUT,BAD)
c
c     Added definition of REAL*4 to avoid the compiler turning these into
c     REAL*8 as that is the default. NPK 8/7/09.
c
      REAL*4 datin,XIN,YIN,datout,XOUT,YOUT,XIN2,DATIN2,BAD
      INTEGER INMAX,JNMAX,IOMAX,JOMAX
      dimension YIN(JNMAX),XIN(INMAX),XOUT(IOMAX),YOUT(JOMAX)
      dimension datin(INMAX,JNMAX),datout(IOMAX,JOMAX)
      dimension XIN2(INMAX+2),datin2(INMAX+2,JNMAX)
 
       dxin=xin(2)-xin(1)
       do j=1,jnmax
        do i=2,inmax+1
         datin2(i,j)=datin(i-1,j)
        enddo
        datin2(1,j)=datin2(inmax+1,j)
        datin2(inmax+2,j)=datin2(2,j)
       enddo
       do i=2,inmax+1
        xin2(i)=xin(i-1)
       enddo
       xin2(1)=xin(1)-dxin
       xin2(inmax+2)=xin(inmax)+dxin
       call cfs_intp2d(datin2,inmax+2,jnmax,xin2,yin
     1            ,datout,iomax,jomax,xout,yout,-1,BAD)
       return
       end

c     SUBROUTINE INTP2D - copied from CFS coupO2A.f (coupler routine)
c     NPK modifications as noted
c
c     Added MFLG argument to match list of arguments in 
c     call from CFS_KPP_TO_1X1.  This doesn't appear to be used, but
c     otherwise the routine thinks that BAD=-1 (see call to CFS_INTP2D above).
c     NPK 8/7/09.
      SUBROUTINE CFS_INTP2D(A,IMX,IMY,XA,YA,B,JMX,JMY,XB,YB,MFLG,BAD)
      INTEGER MAX,MAXP
      PARAMETER(MAX=721,MAXP=721*361)
c
c     Commented out the common blocks as they are not used. NPK 8/7/09.
c      COMMON /COMASK/ DEFLT,MASK(maxp)
c      COMMON /COMWRK/ IERR,JERR,DXP(MAX),DYP(MAX),DXM(MAX),DYM(MAX),
c     1                IPTR(MAX),JPTR(MAX)
c
c     Added definition of REAL*4 to prevent the compiler turning these
c     into REAL*8 as that is the default. NPK 8/7/09.
      REAL*4 A,B,XA,YA,XB,YB,MFLG,BAD,WGT,DXP,DYP,DXM,DYM,
     +     D1,D2,D3,D4,DD
      INTEGER IERR,JERR,IPTR,JPTR
c     Defined these dimensions, which were previously assumed-size.
c     NPK 8/7/09.
      DIMENSION A(IMX,IMY),B(JMX,JMY),XA(IMX),YA(IMY),XB(JMX),YB(JMY),
     +     WGT(4),DXP(MAX),DYP(MAX),DXM(MAX),DYM(MAX),IPTR(MAX),
     +     JPTR(MAX)

      DO I=1,MAX
         IPTR(I)=0
         JPTR(I)=0
      ENDDO

      CALL CFS_SETPTR(XA,IMX,XB,JMX,IPTR,DXP,DXM,IERR)
      CALL CFS_SETPTR(YA,IMY,YB,JMY,JPTR,DYP,DYM,JERR)
c
      DO I = 1,4
         WGT(I) = 1.0
      ENDDO
      DO J = 1,JMY
         DO I = 1,JMX
            B(I,J) = BAD
         ENDDO
      ENDDO
      DO 20 J = 1,JMY
        JM = JPTR(J)
        IF (JM.LT.0) GOTO 20
        JP = JM + 1
      DO 10 I = 1,JMX
        IM = IPTR(I)
        IF (IM.LT.0) GOTO 10
        IP = IM + 1
        D1 = DXM(I)*DYM(J)*WGT(1)
        D2 = DXM(I)*DYP(J)*WGT(2)
        D3 = DXP(I)*DYM(J)*WGT(3)
        D4 = DXP(I)*DYP(J)*WGT(4)
        if (a(ip,jp) .EQ. BAD) d1=0.0
        if (a(ip,jm) .EQ. BAD) d2=0.0
        if (a(im,jp) .EQ. BAD) d3=0.0
        if (a(im,jm) .EQ. BAD) d4=0.0
        DD = D1 + D2 +D3 + D4
        IF (DD.EQ.0.0) GOTO 10
        B(I,J) = (D4*A(IM,JM)+D3*A(IM,JP)+D2*A(IP,JM)+D1*A(IP,JP))/DD
   10 CONTINUE
   20 CONTINUE
      RETURN
      END
c     
c     SUBROUTINE SETPTR - copied from CFS coupO2A.f (coupler routine)
c     NPK modifications as noted
c     Note that arrays must be passed in ordered from lowest to highest,
c     so latitude must run from south to north
c
      SUBROUTINE CFS_SETPTR(X,M,Y,N,IPTR,DP,DM,IERR)
c     Added definition of REAL*4 to prevent the compiler turning these
c     into REAL*8 as that is the default. NPK 8/7/09.
      REAL*4 X,Y,DP,DM
      INTEGER IPTR,IERR
      DIMENSION X(*),Y(*),IPTR(*),DP(*),DM(*)
      IERR = 0
      DO 10 J = 1,N
        YL = Y(J)
        IF (YL.LT.X(1)) THEN
          IERR = IERR + 1
          IPTR(J) = -1
          GOTO 10
        ELSEIF (YL.GT.X(M)) THEN
          IERR = IERR + 1
          IPTR(J) = -1
          GOTO 10
        ENDIF
        DO 20 I = 1,M-1
          IF (YL.GT.X(I+1)) GOTO 20
          DM(J) = YL-X(I)
          DP(J) = X(I+1)-YL
          IPTR(J) = I
          GOTO 10
  20    CONTINUE
  10  CONTINUE
      RETURN
      END

      SUBROUTINE CFS_FILL_KPDS_KGDS(longitude,nlon,latitude,nlat,
     +     century,year,month,date,hour,kpds,kgds)
      IMPLICIT NONE

      INTEGER nlon,nlat,century,year,month,date,hour,i,
     +     kpds(200),kgds(200)
      REAL*4 longitude(nlon),latitude(nlat)

c     Fill the KPDS array
      DO i=1,200
         kpds(i)=0
         kgds(i)=0
      ENDDO
      kpds(1)  = 7               ! Identification of the center
      kpds(2)  = 0               ! Generating process id number
      kpds(3)  = 255             ! Grid definition
      kpds(4)  = 128             ! GDS/BMS flag (right adj copy of octet 8)
      kpds(5)  = 11              ! Parameter identifier
      kpds(6)  = 1               ! Type of level
      kpds(7)  = 0               ! Height/pressure of level
      kpds(8)  = year            ! Year of the current century
      kpds(9)  = month           ! Month of the current year
      kpds(10) = date            ! Day of the current year
      kpds(11) = 0               ! Hour of the current day <- might need to change if sub-daily coupling?
      kpds(12) = 0               ! Minute of the current hour
      kpds(13) = 1               ! Indicator of forecast time unit
      kpds(14) = 0               ! Time range 1
      kpds(15) = 0               ! Time range 2
      kpds(16) = 1               ! Time range flag (should be 1 if number of years in averaging period <= 1)
      kpds(17) = 0               ! Number included in average
      kpds(18) = 1               ! Version number of grib specification
      kpds(19) = 3               ! Version number of parameter table
      kpds(20) = 0               ! Number of points missing from average/accumulation
      kpds(21) = century         ! Century of reference time of data
      kpds(22) = 2               ! Units decimal scale factor (?)
      kpds(23) = 0               ! Subcenter number
      kpds(24) = 0               ! For NMC ensemble products, additive of 128 (forecast error field), 64 (bias corrected), 32 (smoothed)
      kpds(25) = 0               ! Not used
c     Note that it is necessary to complete KPDS only to position 25;
c     the remainder of the parameters are either reserved or filled by
c     the GRIB routines themselves
c
c     Fill the KGDS array
      kgds(1)  = 0                                 ! Data representation type (0 for regular grid, 4 for Gaussian)
      kgds(2)  = nlon                              ! Number of points on a latitude circle (number of points in longitude)
      kgds(3)  = nlat                              ! Number of points on a longitude meridian (number of points in latitude)
      kgds(4)  = NINT(latitude(1)*1000)            ! First latitude point (millidegrees)
      kgds(5)  = NINT(longitude(1)*1000)           ! First longitude point (millidegrees)     
      kgds(6)  = 128                               ! Resolution flag (right adj copy of octet 17)
      kgds(7)  = NINT(latitude(nlat)
     +     *1000)                                  ! Latitude of last point (millidegrees)
      IF (longitude(nlon) .gt. 180) THEN
         kgds(8) = NINT((longitude(nlon)
     +        -360)*1000)                          ! Longitude of last point (millidegrees)
      ELSE
         kgds(8) = NINT(longitude(nlon)
     +        *1000)                               ! Longitude of last point (millidegrees)
      ENDIF
      kgds(9)  = NINT(ABS((longitude(2)-
     +     longitude(1)))*1000)                    ! Increment in longitude (millidegrees)
      kgds(10) = NINT(ABS((latitude(2)-
     +     latitude(1)))*1000)                     ! Increment in latitude
      kgds(11) = 0                                 ! Scanning mode flags
      kgds(19) = 0                                 ! Number of vertical coordinate parameters
      kgds(20) = 255                               ! Octet number of list of vertical coordinate parameters or octet number of the 
                                                   ! list of numbers of points in each row, or 255 if neither are present
      kgds(21) = 0                                 ! Not used
      kgds(22) = 0                                 ! Not used
      DO i=12,18
         kgds(i) = 0                               ! Positions 12-18 not used
      ENDDO

      WRITE(6,*) 'KPDS = ',kpds
      WRITE(6,*)
      WRITE(6,*) 'KGDS = ',kgds

      RETURN
      END

      SUBROUTINE CFS_REVERSE_LATITUDE_INPLACE(x,latitude,nx,ny)
      IMPLICIT NONE
c
c     Reverses the second dimension of a two-dimensional array (x)
c     and its associated one-dimensional coordinate array (latitude).
c
c     Note that this reversal is done IN PLACE (i.e., the input array
c     is overwritten).
c
c     Nicholas P. Klingaman 9/7/09
c
      INTEGER NX,NY,i,j
      REAL*4 x(NX,NY),latitude(NY),half_globe_temp(NY/2)
     
      DO j=1,NY/2
         half_globe_temp(j)=latitude(j)
      ENDDO
      DO j=NY,NY/2,-1
         latitude(NY-j+1)=latitude(j)
      ENDDO
      DO j=1,NY/2
         latitude(NY-j+1)=half_globe_temp(j)
      ENDDO
      DO i=1,NX
         DO j=1,NY/2
            half_globe_temp(j)=x(i,j)
         ENDDO
         DO j=NY,NY/2,-1
            x(i,NY-j+1)=x(i,j)
         ENDDO
         DO j=1,NY/2
            x(i,NY-j+1)=half_globe_temp(j)
         ENDDO
      ENDDO

      RETURN
      END

      
