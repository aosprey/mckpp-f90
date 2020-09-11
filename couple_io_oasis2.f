#ifdef COUPLE
#ifdef OASIS2
      SUBROUTINE coupled_flux(solar,non_solar,PminusE,kt)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

c     Include files for KPP
#include "parameter.inc"
#include "times.com"
#include "couple.com"
#include "location.com"

c     Include files for OASIS2
#include "clim.h"
#include "oasis.h"
#include "param_cou.h"
#include "inc_cpl.h"

      REAL solar(NPTS),non_solar(npts),PminusE(npts)

      REAL*8 solar_in(NX_GLOBE,NY_GLOBE),
     +     non_solar_in(NX_GLOBE,NY_GLOBE),
     +     PminusE_in(NX_GLOBE,NY_GLOBE)
     
      INTEGER kt,info,jf,ipoint,ix,jy

      IF (cchan .EQ. 'CLIM') THEN
         DO jf=1,jpflda2o1
            IF (jf .eq. 1) CALL CLIM_Import (cl_read(jf), kt, 
     +           non_solar_in,info)
            IF (jf .eq. 2) CALL CLIM_Import (cl_read(jf), kt, 
     +           solar_in,info)
            IF (jf .eq. 3) CALL CLIM_Import (cl_read(jf), kt, 
     +           PminusE_in,info)
            IF (info .NE. CLIM_Ok) THEN
               WRITE(nuerr,*)'Error in reading ',cl_read(jf), jf
               WRITE(nuerr,*)'Coupler timestep is ',kt
               WRITE(nuerr,*)'CLIM error code is ',info
               CALL halte ('STOP in coupled_flux')
            ENDIF
         ENDDO
      ENDIF    

      DO jy=jfirst,jlast
         DO ix=ifirst,ilast
c     
c     Debugging WRITE statements - uncomment if required
c     NPK 2/11/09 - R3 (now commented by default)
c
c            IF (ix .eq. ifirst+10 .and. jy .eq. jfirst+10) THEN
c               WRITE(nuout,*)'At ',ifirst+10,',',jfirst+10,' solar_in=',
c     +              solar_in(ix,jy),' non_solar_in=',non_solar_in(ix,jy)
c     +              ,' PminusE_in=',PminusE_in(ix,jy)
c            ENDIF
c
c     Transform global fields (as received from atmosphere)
c     to regional fields (as required for KPP)
c     
            ipoint=(jy-jfirst)*nx+(ix-ifirst)+1
            solar(ipoint)=solar_in(ix,jy)
            non_solar(ipoint)=non_solar_in(ix,jy)
            PminusE(ipoint)=PminusE_in(ix,jy)
         ENDDO
      ENDDO
      
      RETURN
      END

      SUBROUTINE coupled_stress(ustress,vstress,kt)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

c     Include files for KPP
#include "parameter.inc"      
#include "times.com"
#include "couple.com"
#include "location.com"

      REAL ustress(NPTS),vstress(npts)
      INTEGER kt,info
      REAL ustress_in(NX_GLOBE,NY_GLOBE),vstress_in(NX_GLOBE,NY_GLOBE)
      
c     Include files for OASIS2
#include "clim.h"
#include "oasis.h"
#include "param_cou.h"
#include "inc_cpl.h"

      INTEGER jf,ipoint,ix,jy,jf_in

      IF (cchan .EQ. 'CLIM') THEN
         DO jf=1,jpflda2o2
            jf_in=jf+jpflda2o1
            IF (jf .eq. 1) CALL CLIM_Import (cl_read(jf_in), kt, 
     +           ustress_in,info)
            IF (jf .eq. 2) CALL CLIM_Import (cl_read(jf_in), kt, 
     +           vstress_in,info)
            IF (info .NE. CLIM_Ok) THEN
               WRITE(nuerr,*)'Error in reading ',cl_read(jf_in), jf
               WRITE(nuerr,*)'Coupler timestep is ',kt
               WRITE(nuerr,*)'CLIM error code is ',info
               CALL halte ('STOP in coupled_stress')
            ENDIF
         ENDDO
      ENDIF
c
c     Transform global fields (as received from atmosphere)
c     to regional fields (as required for KPP)
c
      DO jy=jfirst,jlast
         DO ix=ifirst,ilast
            ipoint=(jy-jfirst)*nx+(ix-ifirst)+1
            ustress(ipoint)=ustress_in(ix,jy)
            vstress(ipoint)=vstress_in(ix,jy)
         ENDDO
      ENDDO
      
      RETURN
      END

      SUBROUTINE coupled_out(X,kt,L_LAST)

c***********************************************************
c     Routine to output the fields to the coupler
c     Has to deal with regional coupling
c     The way it is written one can only do regional coupling in 
c     the longitude direction IF (ifirst .NE. 1)
c***********************************************************
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

c     Include files for KPP
#include "parameter.inc"
#include "constants.com"
#include "times.com"
#include "couple.com"
#include "location.com"
#include "ocn_advec.com"

      REAL X(NPTS,NZP1,NSCLR)
      INTEGER kt

c     Include files for OASIS2
#include "clim.h"
#include "oasis.h"
#include "param_cou.h"
#include "inc_cpl.h"

      INTEGER ix,jy,ipoint,ipoint_globe
      REAL ICE_in(NX_GLOBE,NY_GLOBE,1)
      REAL usf_in(NX_GLOBE,NY_GLOBE)
      REAL vsf_in(NX_GLOBE,NY_GLOBE)
      REAL SST(NX_GLOBE,NY_GLOBE)
      REAL ICE(NX_GLOBE,NY_GLOBE)
      REAL SST_in(NX_GLOBE,NY_GLOBE,1)

      CHARACTER*8 file_name(jpmaxfld)
      INTEGER max_file
      INTEGER file_unit_max,file_unit(jpmaxfld)
      INTEGER file_unit_field(jpmaxfld)
      
      INTEGER icstep, info, jn, jf, ierror
      LOGICAL L_TROUVE,L_LAST
            
      COMMON /save_sstin/ SST_in,ICE_in,usf_in,vsf_in
c
c     Use the coupling weight to modify the SSTs before passing
c     back to the coupler.
c     N.B. : We want to do this whether or not the user has specified
c     coupling weights.  Cplwght is automatically set IF (.NOT. L_CPLWGHT)
c     NPK April 2008 - R1
c      
      WRITE(nuout,*) ' Entering loop to weight SSTs'
      DO ix=1,NX_GLOBE
         DO jy=1,NY_GLOBE
            ipoint_globe = (jy-1)*NX_GLOBE+ix
            IF (cplwght(ipoint_globe) .LT. -1e-10) THEN              
c     Point is outside the coupling domain; set to SST climatology
               SST(ix,jy) = SST_in(ix,jy,1)
            ELSE
c     Point is inside the coupling domain; set to weighted value
               ipoint=(jy-jfirst)*nx+(ix-ifirst)+1
               SST(ix,jy) = X(ipoint,1,1)*cplwght(ipoint_globe) +
     &              SST_in(ix,jy,1)*(1-cplwght(ipoint_globe))
            ENDIF
            ice(ix,jy) = ice_in(ix,jy,1)
         ENDDO
      ENDDO
      IF (L_OUTKELVIN) SST = SST+TK0
c         
c     This bit is taken almost directly from stpmco in the toymodel
c      
      WRITE(nuout,*) ' '
      WRITE(nuout,*) 'stpcmo: sending fields to CPL with kt= ', kt
      WRITE(nuout,*) ' '

      IF (L_LAST) THEN 
c
c     -WRITE fields to binary files for coupler restart at last time step
c
c         -initialisation
c
          max_file=1
          file_unit_max=99
c         -keeps first file name
          file_name(max_file)=cl_f_writ(max_file)
c         -keeps first file unit
          file_unit(max_file)=file_unit_max
c         -decrements file unit maximum
          file_unit_max=file_unit_max-1
c         -keeps file unit for field
          file_unit_field(1)=file_unit(max_file)
c
c         -different files names counter
c
          
          DO jf= 2, jpfldo2a
             L_TROUVE=.false.
             DO jn= 1, max_file
                IF (.not. L_TROUVE) THEN
                   IF (cl_f_writ(jf).EQ.file_name(jn)) THEN
c     keep file unit for field
                      file_unit_field(jf)=file_unit(jn)
                      L_TROUVE=.true.
                   END IF 
                END IF 
             END DO 
             IF (.not. L_TROUVE) then
c     -increment the number of different files
                max_file=max_file+1
c     -keep file name
                file_name(max_file)=cl_f_writ(jf)
c     -keep file unit for file
                file_unit(max_file)=file_unit_max
c     -keep file unit for field
                file_unit_field(jf)=file_unit(max_file)
c     -decrement unit maximum number from 99 to 98, ...
                file_unit_max=file_unit_max-1
             END IF 
          END DO 
c     
          DO jn=1, max_file 
             OPEN (file_unit(jn), FILE=file_name(jn), 
     &            FORM='UNFORMATTED')
          END DO 
c     
c     WRITE fields to files           
c
          WRITE(nuout,*) 'KPP: Writing to output files in couple_io.f'
          CALL FLUSH(nuout)
          DO jf=1, jpfldo2a
             IF (jf.eq.2)
     $            CALL locwrite(cl_writ(jf),sst, NX_GLOBE*NY_GLOBE,
     $            file_unit_field(jf), ierror, nuout) 
             IF (jf.eq.1) 
     $            CALL locwrite(cl_writ(jf),ice, NX_GLOBE*NY_GLOBE,
     $            file_unit_field(jf), ierror, nuout)
             IF (jf.eq.3) 
     $            CALL locwrite(cl_writ(jf),usf_in, NX_GLOBE*NY_GLOBE,
     $            file_unit_field(jf), ierror, nuout)
             IF (jf.eq.4) 
     $            CALL locwrite(cl_writ(jf),vsf_in, NX_GLOBE*NY_GLOBE,
     $            file_unit_field(jf), ierror, nuout)
          END DO 
C     
C     -simulate a FLUSH
C
          WRITE(nuout,*) 'KPP: Closing output files in couple_io.f'
          CALL FLUSH(nuout)
          DO jn=1, max_file 
             CLOSE (file_unit(jn))
          END DO 
c     
c     
c     
          IF(cchan.eq.'CLIM') THEN 
c
c         -inform PVM daemon, I have finished
c
             WRITE(nuout,*) 'KPP: Calling CLIM_Quit'
             CALL FLUSH(nuout)
             CALL CLIM_Quit (CLIM_StopPvm, info)
             IF (info .NE. CLIM_Ok) THEN
                WRITE (6,*) 
     $               'An error occured while leaving CLIM. Error = ',
     $               info
             ENDIF
             
          END IF 
          RETURN             
       END IF 
c     
       IF(cchan.eq.'CLIM') THEN
C     
C     -Give oceanic fields to Oasis
C     
          DO jn=1, jpfldo2a
             
             IF (jn.eq.2) THEN
                CALL CLIM_Export(cl_writ(jn), kt, sst, info)
             ENDIF
             IF (jn.eq.1) THEN
                CALL CLIM_Export(cl_writ(jn), kt, ice, info)
             ENDIF
             IF (jn.eq.3) THEN
                CALL CLIM_Export(cl_writ(jn), kt, usf_in, info)
             ENDIF            
             IF (jn.eq.4) THEN
                CALL CLIM_Export(cl_writ(jn), kt, vsf_in, info)
             ENDIF
c             
             IF (info .NE. CLIM_Ok) THEN
                WRITE (nuerr,*) 'STEP : Pb giving ',cl_writ(jn), ':',jn
                WRITE (nuerr,*) ' at timestep = ', icstep,'kt = ',kt
                WRITE (nuerr,*) 'Clim error code is = ',info
                CALL halte('STOP in coupled_out ')
             ENDIF
          END DO 
       ENDIF
       
       RETURN
       END
      

      SUBROUTINE locread ( cdfldn, pfield, kdimax, knulre, kflgre
     $                   , kout)
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 0 *
C               * -------------     ------- *
C               *****************************
C
C**** *locread*  - Read binary field on unit knulre
C
C     Purpose:
C     -------
C     Find string cdfldn on unit knulre and read array pfield
C
C**   Interface:
C     ---------
C       *CALL*  *locread (cdfldn, pfield, kdimax, knulre, kflgre, kout)*
C
C     Input:
C     -----
C                cdfldn : character string locator
C                kdimax : dimension of field to be read 
C                knulre : logical unit to be read 
C                kout   : logical unit to write messages
C
C     Output:
C     ------
C                pfield : field array (real 1D)
C                kflgre : error status flag
C
C     Reference:
C     ---------
C     See OASIS manual (1995) 
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       2.0       L. Terray      95/09/01  created
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C

      IMPLICIT CHARACTER(c)
      IMPLICIT LOGICAL(l)
C
C* ---------------------------- Argument declarations -------------------
C
      REAL pfield(kdimax)
      CHARACTER*8 cdfldn
C
C* ---------------------------- Local declarations ----------------------
C
      CHARACTER*8 clecfl
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Initialization
C        --------------
C
      WRITE (UNIT = kout,FMT = 1001) knulre
C
C* Formats
C
 1001 FORMAT('Locread : Read binary file connected to unit = ',I3)
C
C     2. Find field in file
C        ------------------
C
      REWIND knulre
 200  CONTINUE
C* Find string
      READ (UNIT = knulre, ERR = 200, END = 210) clecfl
      IF (clecfl .NE. cdfldn) GO TO  200
C* Read associated field
      READ (UNIT = knulre, ERR = 210, END = 210) pfield
C* Reading done and ok
      kflgre = 0
      GO TO 220
C* Problem in reading
 210  kflgre = 1
 220  CONTINUE
C
C
C*    3. End of routine
C        --------------
C
      WRITE (UNIT = kout,FMT = *) 'Locread : done'
c     CALL FLUSH (kout) ! removed to compile on HCPX
      RETURN
      END

      SUBROUTINE locwrite (cdfldn, pfield, kdimax, knulre, kflgre, kout)
      IMPLICIT none
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL 0 *
C               * -------------     ------- *
C               *****************************
C
C**** *locwrite*  - Write binary field on unit knulre
C
C     Purpose:
C     -------
C     Write string cdfldn and array pfield on unit knulre
C
C**   Interface:
C     ---------
C       *CALL*  *locwrite (cdfldn, pfield, kdimax, knulre, kflgre, kout)*
C
C     Input:
C     -----
C                cdfldn : character string locator
C                kdimax : dimension of field to be written 
C                knulre : logical unit to be written
C                pfield : field array (real 1D) 
C                kout   : logical unit to write messages
C
C     Output:
C     ------
C                kflgre : error status flag
C
C     Workspace:
C     ---------
C     None
C
C     Externals:
C     ---------
C     None
C
C     Reference:
C     ---------
C     See OASIS manual (1995) 
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       2.0       L. Terray      95/09/01  created
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C
C* ---------------------------- Argument declarations -------------------
C
      INTEGER kdimax, knulre, kflgre, kout
      REAL pfield(kdimax)
      CHARACTER*8 cdfldn
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Initialization
C        --------------
C
      WRITE (UNIT = kout,FMT = 1001) knulre
C
C* Formats
C
 1001 FORMAT(5X,' Write binary file connected to unit = ',I3)
C
C     2. Find field in file
C        ------------------
C
C* Write string
      WRITE (UNIT = knulre, ERR = 210) cdfldn
C* Write associated field
      WRITE (UNIT = knulre, ERR = 210) pfield
C* Writing done and ok
      kflgre = 0
      GO TO 220
C* Problem in Writing
 210  kflgre = 1
 220  CONTINUE
C
C
C*    3. End of routine
C        --------------
C
      WRITE (UNIT = kout,FMT = *) 'Locwrite : done'
c     CALL FLUSH (kout) ! removed to complile on HPCX
      RETURN
      END



      SUBROUTINE halte (cdtext)
C****
C               *****************************
C               * OASIS ROUTINE  -  LEVEL C *
C               * -------------     ------- *
C               *****************************
C
C**** *halte*  - Abort the program
C
C     Purpose:
C     -------
C     Print an error message to standard output and abort the coupler
C
C**   Interface:
C     ---------
C       *CALL*  *halte (cdtext)*
C
C     Input:
C     -----
C                cdtext   : character string to be printed
C
C     Output:
C     ------
C     None
C
C     Workspace:
C     ---------
C     None
C
C     Externals:
C     ---------
C     None
C
C     Reference:
C     ---------
C     See OASIS 2.2 manual (1997) 
C
C     History:
C     -------
C       Version   Programmer     Date      Description
C       -------   ----------     ----      -----------  
C       2.2       S. Valcke      97/11/18  created
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C* ---------------------------- Argument declarations ----------------------
C
      CHARACTER*(*) cdtext
C
C* ---------------------------- Local declarations ----------------------
C
      CHARACTER cpbase
      CHARACTER*10 cprpt, cpdots
      CHARACTER*69 cline
      PARAMETER ( cpbase = '-' )
      PARAMETER ( cprpt = '* ===>>> :' )
      PARAMETER ( cpdots = '  ------  ' )
C
C* ---------------------------- Poema verses ----------------------------
C
C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C
C*    1. Print text 
C        ----------
C
          cline = ' '
          ilen = len(cdtext)
          DO 110 jl = 1, ilen
            cline(jl:jl) = cpbase
 110      CONTINUE
          nulou=6
          WRITE(UNIT = nulou,FMT='(/,A,X,A)') cpdots, cline
          WRITE(UNIT = nulou,FMT='(/,A,X,A,/)') cprpt, cdtext
          WRITE(UNIT = nulou,FMT='(A,X,A,/)')cpdots, cline
C
C
C*    2. FLUSH the coupler output
C        ------------------------
C
c     CALL FLUSH (nulou) ! removed to compile on HCPx
C
C
C*    3. Abort the coupler
C        -----------------
C
      CALL abort
      RETURN
      END

      SUBROUTINE TWO_TO_ONE(twod,oned)

      IMPLICIT NONE

#include "parameter.inc"
  
      REAL twod(NX_GLOBE,NY_GLOBE),oned(NX_GLOBE*NY_GLOBE)
      INTEGER*4 i,j,point

      DO i=1,NX_GLOBE
         DO j=1,NY_GLOBE
            point = (i-1)*NY_GLOBE+j
            oned(point) = twod(i,j)
         END DO
      END DO
         
      RETURN
      END
#endif /*OASIS2*/
#endif /*COUPLE*/
