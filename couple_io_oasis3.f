#ifdef COUPLE
#ifdef OASIS3
      SUBROUTINE mpi1_oasis3_input(solar,non_solar,PminusE,ustress,
     +     vstress,curl_tau,kpp_3d_fields,kpp_const_fields)
c
c     Facilitate the exchange of coupled fields via the OASIS3 coupler.
c     NPK 18/09/09 for the OASIS3 toy model - completed 28/09/09 - R3
c
c     Use OASIS3 (PRISM) modules - the compiler must be able to find
c     these via the Makefile.
c
#ifdef OASIS3_MCT
      USE mod_prism
#else
      USE mod_kinds_model
      USE mod_prism_proto
      USE mod_prism_def_partition_proto
      USE mod_prism_put_proto
      USE mod_prism_get_proto
      USE mod_prism_grids_writing
#endif
c
      IMPLICIT NONE
c
c     Standard include files
c
#include <netcdf.inc>
! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
#include "kpp_oasis3.inc"
#include "times.com"
#include "constants.com"
#include "couple.com"
      TYPE(kpp_const_type) :: kpp_const_fields
      TYPE(kpp_3d_type) :: kpp_3d_fields
c
c     Output variables on the KPP regional grid - returned to
c     the calling routine (usually <fluxes>).
c
      REAL solar(NPTS)
      REAL non_solar(NPTS)
      REAL PminusE(NPTS)
      REAL ustress(NPTS)
      REAL vstress(NPTS)
      REAL curl_tau(NPTS)
c
c     Local variables
c     
c     Note: "time_in_seconds" must be defined as INTEGER to
c     be consistent with the definition in OASIS3.
c
c     Note: the "temporary" variable must be defined
c     as the OASIS3 type "ip_realwp_p".  OASIS3 must be compiled
c     with default double precision (-fdefault-real-8 or similar).
c
#ifdef TOYCLIM
      REAL(KIND=ip_realwp_p) temporary(NPTS_GLOBE)
#else
      REAL(KIND=ip_realwp_p) temporary(NX_GLOBE,NY_GLOBE)
#endif
      REAL rain(NPTS),evap(NPTS),runoff(NPTS),runoff_mean,
     +     weights(NPTS)
      INTEGER i,ix,j,jy,ipt,ierror,npts_ocean,my_jpfldin
      INTEGER(KIND=4) time_in_seconds
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
c
c     Note - time in seconds for receive should be the end of the
c     previous timestep, not the end of this timestep.  Subtract 1
c     timestep.
c
      time_in_seconds=NINT((kpp_const_fields%ntime-1)*
     +     kpp_const_fields%dto)
      WRITE(nuout,*) 'KPP: Time for coupling in is ',
     +     time_in_seconds,' seconds'
c      
c     Get the coupled fields from the OASIS coupler.  Note that you
c     can discard any fields you do not want by simply not defining
c     a CASE for them.
c

      ! Initialise all coupled inputs to zero.
      solar(:)=0.0
      non_solar(:)=0.0
      PminusE(:)=0.0
      ustress(:)=0.0
      vstress(:)=0.0
      curl_tau(:)=0.0

      ! If not passing river runoff, need to reduce number of input
      ! fields by one
      IF (.NOT. kpp_const_fields%L_DIST_RUNOFF) THEN
         my_jpfldin=jpfldin-1
      ELSE
         my_jpfldin=jpfldin
      ENDIF
      
      DO i=1,my_jpfldin
         CALL prism_get_proto(il_var_id_in(i),
     +        time_in_seconds,temporary,ierror)
!         WRITE(6,*) 'KPP: For field number ',i,' called ',cl_read(i),
!     +        ' received ierror = ',ierror
!         WRITE(6,*) 'KPP: For field number ',i,' called ',cl_read(i),
!     +        ' unweighted global sum = ',SUM(temporary)
         IF (ierror.NE.PRISM_Ok .and. ierror .LT. PRISM_Recvd) THEN
            WRITE(il_mparout) 'KPP: Received error from ',
     +           'PRISM_Get_Proto =',ierror,' receiving variable ',
     +           cl_read(i),' at model time ',time_in_seconds,' sec.'
            WRITE(il_mparout) 'KPP: Aborting coupled integration ...'
            CALL prism_abort_proto(il_comp_id,'couple_io_oasis3.f',
     +           'get')
         ELSE
            SELECT CASE (cl_read(i))
c
c     For each field that we are coupling, use <TWOD_GLOBAL_ONED_REGIONAL>
c     to transform that field from the global atmospheric grid to the
c     KPP regional grid.
c
            CASE('HEATFLUX')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,non_solar)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,non_solar)
#endif
            CASE('SOLAR')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,solar)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,solar)
#endif
	    CASE('PEN_SOL')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,solar)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,solar)
#endif

            CASE('TAUX')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,ustress)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,ustress)
#endif
            CASE('TAUY')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,vstress)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,vstress)
#endif
            CASE ('TRAIN')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,rain)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,rain)
#endif
            CASE ('EVAP2D')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,evap)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,evap)
#endif
            CASE ('RUNOFF')
#ifdef TOYCLIM
               CALL ONED_GLOBAL_ONED_REGIONAL(temporary,runoff)
#else
               CALL TWOD_GLOBAL_ONED_REGIONAL(temporary,runoff)
#endif

            CASE DEFAULT
               WRITE(il_mparout,*) 'KPP: Discarding field ',
     +              cl_read(i)
            END SELECT
         ENDIF
      ENDDO
      IF (kpp_const_fields%L_DIST_RUNOFF) THEN
         weights(:) = 0.
         DO i=1,NPTS
            IF (kpp_3d_fields%L_OCEAN(i))
     +           weights(i) = COS(kpp_3d_fields%dlat(i)*3.14159/180.)     
         ENDDO
         weights = weights/SUM(weights)
         runoff_mean = SUM(runoff*weights)
!         WRITE(6,*) 'KPP: Domain-mean freshwater flux river '//
!     +        ' runoff at time ',time_in_seconds,' seconds = ',
!     +        runoff_mean,' mm/s'
         kpp_3d_fields%runoff_incr = runoff_mean
         runoff_mean = SUM(kpp_3d_fields%runoff_incr*weights)
!         WRITE(6,*) 'KPP: Domain-mean increment to P minus E '//
!     +        ' at time ',time_in_seconds,' seconds = ',
!     +        runoff_mean,' mm/s'
      ELSE
         kpp_3d_fields%runoff_incr(:)=0
         runoff_mean=0
      ENDIF
      DO ix=ifirst,ilast
         DO jy=jfirst,jlast
            ipt=(jy-jfirst)*NX+(ix-ifirst)+1
            IF (kpp_3d_fields%L_OCEAN(ipt) .and.
     +           kpp_const_fields%L_NO_EGTP .and. 
     +           ix .ge. kpp_const_fields%barrier_ifirst .and. 
     +           ix .le. kpp_const_fields%barrier_ilast .and.
     +           jy .ge. kpp_const_fields%barrier_jfirst .and.
     +           jy .le. kpp_const_fields%barrier_jlast) THEN
!     Remove any flux of freshwater out of the ocean
               IF (rain(ipt)+kpp_3d_fields%runoff_incr(ipt)
     +              .lt. evap(ipt)) 
     +              evap(ipt) = rain(ipt)+kpp_3d_fields%runoff_incr(ipt)
            ENDIF
         ENDDO
      ENDDO
      DO ipt=1,NPTS
         PminusE(ipt)=rain(ipt)+
     +        kpp_3d_fields%runoff_incr(ipt)-evap(ipt)
      ENDDO

!     Set curl of wind stress to zero for now.  Need to add
!     code for computing curl from input taux and tauy.
      curl_tau(:)=0.0
      

c      IF (READ_FROM_NETCDF)
c     +     CALL mpi1_oasis3_read_netcdf(solar,non_solar,PminusE,
c     +     ustress,vstress)
c      IF (WRITE_TO_NETCDF) 
c     +     CALL mpi1_oasis3_write_netcdf(solar,non_solar,PminusE,
c     +     ustress,vstress)         
c
c     End
c     
      RETURN      
      END SUBROUTINE mpi1_oasis3_input
      
      SUBROUTINE mpi1_oasis3_output(kpp_3d_fields,kpp_const_fields)
c     
c     Facilitate the exchange of coupled fields via the OASIS3 coupler.
c     NPK 18/09/09 for the OASIS3 toy model - completed 28/09/09 - R3
c
c     Use OASIS3 (PRISM) modules - the compiler must be able to find
c     these via the Makefile.
c
#ifdef OASIS3_MCT
      USE mod_prism
#else
      USE mod_kinds_model
      USE mod_prism_proto
      USE mod_prism_def_partition_proto
      USE mod_prism_put_proto
      USE mod_prism_get_proto
      USE mod_prism_grids_writing
#endif
c      
      IMPLICIT NONE
c
c     Standard include files
c
#include <netcdf.inc>
! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
#include "kpp_oasis3.inc"
#include "times.com"
#include "couple.com"
c#include "location.com"
#include "constants.com"
#include "currclim.com"
c#include "initialcon.com"

      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
c
c     Pass in the regional SST from the calling routine
c     (usually <MAIN>).
c
      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
c
c     "SST_in" and "ICE_in" are common-block variables from
c     the netCDF input routines.
c
      REAL SST_in(NX_GLOBE,NY_GLOBE,1)
      REAL ICE_in(NX_GLOBE,NY_GLOBE,1)
      REAL ICEDEPTH_in(NX_GLOBE,NY_GLOBE,1)
      REAL SNOWDEPTH_in(NX_GLOBE,NY_GLOBE,1)
      REAL usf_in(NX_GLOBE,NY_GLOBE)
      REAL vsf_in(NX_GLOBE,NY_GLOBE)
c
c     Global SST and ICE variables that will be exported to OASIS
c
      REAL SST(NPTS_GLOBE)
      REAL ICE(NPTS_GLOBE)
      REAL SNOWDEPTH(NPTS_GLOBE)
      REAL ICEDEPTH(NPTS_GLOBE)
      REAL SURF_CURR_X(NPTS_GLOBE)
      REAL SURF_CURR_Y(NPTS_GLOBE)
c
c     Temporary variable to allow the use of a SELECT CASE block.
c     Note that this must be defined as "ip_realwp_p" and OASIS3 must
c     be compiled with default REAL*8 (-fdefault-real-8 or similar).
c
#ifdef TOYCLIM
      REAL(KIND=ip_realwp_p) temporary(NX_GLOBE*NY_GLOBE)
#else
      REAL(KIND=ip_realwp_p) temporary(NX_GLOBE,NY_GLOBE)
#endif
c
c     Other local variables
c     Note: "time_in_seconds" must be INTEGER to agree with the
c     definition in OASIS3.
c
c     Temporary variable to hold smoothed SST and SST anomaly
      REAL, allocatable :: SST_smooth(:,:), SST_anom(:,:)
c     Temporary variable to hold SST anomaly
      INTEGER i,ix,jy,ipoint_globe,ipoint,ierror
      INTEGER(KIND=4) time_in_seconds
c
c     COMMON block for SST_in and ICE_in
c
      COMMON /save_sstin/ SST_in,ICE_in,icedepth_in,snowdepth_in,
     +     usf_in,vsf_in
c
c     Note: the time in the send needs to be the time at the end
c     of this timestep.
c     
      time_in_seconds=NINT(kpp_const_fields%ntime*kpp_const_fields%dto)
!      WRITE(nuout,*) 'KPP: Time for coupling out is ',
!     +     time_in_seconds,' seconds'
c
c     Use the coupling weight to modify the SSTs before passing
c     back to the coupler.
c     N.B. : We want to do this whether or not the user has specified
c     coupling weights.  Cplwght is automatically set IF (.NOT. L_CPLWGHT)
c      
!      WRITE(nuout,*) 'KPP: Creating coupled output fields'
      DO ix=1,NX_GLOBE
         DO jy=1,NY_GLOBE
            ipoint_globe = (jy-1)*NX_GLOBE+ix
            IF (kpp_3d_fields%cplwght(ipoint_globe) .LT. -1e-10) THEN              
c     Point is outside the coupling domain; set to SST climatology
               SST(ipoint_globe) = SST_in(ix,jy,1)
c               IF (.NOT. L_CLIMCURR) THEN
               SURF_CURR_X(ipoint_globe)=0.
               SURF_CURR_Y(ipoint_globe)=0.
c               ELSE
c                 SURF_CURR_X(ipoint_globe)=usf_in(ix,jy)
c                 SURF_CURR_Y(ipoint_globe)=vsf_in(ix,jy)
c               ENDIF
            ELSE
c     Point is inside the coupling domain; set to weighted value
               ipoint=(jy-jfirst)*nx+(ix-ifirst)+1             
      	       IF (kpp_const_fields%L_SST_LAG) THEN
                  SST(ipoint_globe) = kpp_3d_fields%sst_lag(ipoint)*
     +                 kpp_3d_fields%cplwght(ipoint_globe)+
     +                 SST_in(ix,jy,1)*
     +                 (1-kpp_3d_fields%cplwght(ipoint_globe))
  	       ELSE
                  SST(ipoint_globe) = kpp_3d_fields%X(ipoint,1,1)*
     +                 kpp_3d_fields%cplwght(ipoint_globe)+
     +                 SST_in(ix,jy,1)*
     +                 (1-kpp_3d_fields%cplwght(ipoint_globe))
	       ENDIF
               IF (L_COUPLE_CURRENTS) THEN
                  SURF_CURR_X(ipoint_globe)=kpp_3d_fields%U(ipoint,1,1)*
     +                 kpp_3d_fields%cplwght(ipoint_globe)
                  SURF_CURR_Y(ipoint_globe)=kpp_3d_fields%U(ipoint,1,2)*
     +                 kpp_3d_fields%cplwght(ipoint_globe)
c               ELSEIF (L_COUPLE_CURRENTS .AND. L_CLIMCURR) THEN
c                  SURF_CURR_X(ipoint_globe)=usf_in(ix,jy)
c                  SURF_CURR_Y(ipoint_globe)=vsf_in(ix,jy)
               ENDIF
            ENDIF
            ice(ipoint_globe)=ice_in(ix,jy,1)
c
c     If the user does not provide a climatological ice depth,
c     then set the ice depth to be 2 metres.  Note that we set
c     the depth to be 2 * the ice concentration ** 2 because HadGEM3-A
c     divides the provided depth by the ice concentration to obtain
c     the mean depth over the ice-covered portion of the gridbox, 
c     assuming that the ocean model provides the mean depth over the
c     the entire gridbox.
c
c     The previous, erroneous behaviour of setting our icedepth
c     variable to 2m can be restored by setting L_BAD_ICE_DEPTH.
c            
c     NPK updated 25/6/14.
c
            IF (.NOT. L_CLIM_ICE_DEPTH) THEN
               IF (ice(ipoint_globe) .GE. 0) THEN
                  IF (.NOT. L_BAD_ICE_DEPTH) THEN
                     ipoint=(jy-jfirst)*nx+(ix-ifirst)+1
                     icedepth(ipoint_globe)=2.*ice(ipoint_globe)**2                     
                     IF (icedepth(ipoint_globe) .gt. 2)
     +                    icedepth(ipoint_globe) = 2.00                        
                  ELSE
                     icedepth(ipoint_globe)=2.
                  ENDIF
               ELSE
                  icedepth(ipoint_globe)=0.00
               ENDIF
            ELSE
               icedepth(ipoint_globe)=icedepth_in(ix,jy,1)
            ENDIF
c
c     If there are no data available for the amount of snow on the ice,
c     then assume that there isn't any.  This is, again, based on the
c     default behavior in HadGEM3-A when running with AMIP-2 sea ice.
c
c     NPK 16/12/09 - R3
c
            IF (.NOT. L_CLIM_SNOW_ON_ICE) THEN
	       snowdepth(ipoint_globe)=0.00
    	       IF (kpp_const_fields%L_SST_LAG_FUDGE) THEN
c     
c     Send a lagged SST through the coupler using the snowdepth on ice field.
c     NPK 08/03/19
                  IF (kpp_3d_fields%cplwght(ipoint_globe) 
     +			.LT. 1e-10 .or. .not. kpp_3d_fields%L_OCEAN(ipoint)) THEN
                     snowdepth(ipoint_globe)=0.0
                  ELSE
                     ipoint=(jy-jfirst)*nx+(ix-ifirst)+1
                     snowdepth(ipoint_globe)=TK0+
     +                    kpp_3d_fields%sst_lag(ipoint)*
     +                    kpp_3d_fields%cplwght(ipoint_globe)+
     +                    SST_in(ix,jy,1)*
     +                    (1.0-kpp_3d_fields%cplwght(ipoint_globe))
                     snowdepth(ipoint_globe)=
     +  		  snowdepth(ipoint_globe)/1.0e6
		  ENDIF
               ENDIF
            ELSE
               snowdepth(ipoint_globe)=snowdepth_in(ix,jy,1)
	    ENDIF	 
         ENDDO
      ENDDO
      !IF (L_SST_LAG_FUDGE) WRITE(6,*) 'Lagged SST: ',snowdepth 
       
!     WRITE(il_mparout,*) 'KPP: Finished creating coupled output fields'
!     WRITE(nuout,*) 'KPP: Finished creating coupled output fields'

      DO i=1,jpfldout
c     
c     Export each field to OASIS.  Use a SELECT CASE block to avoid
c     repeated bits of code.  Use the "temporary" variable to transfer
c     the SST and ICE fields to the OASIS "ip_realwp_p" TYPE.
c     
         SELECT CASE (cl_writ(i))
         CASE('OCN_SST')
#ifdef TOYCLIM
            temporary=SST
#else               
            CALL ONED_GLOBAL_TWOD_GLOBAL(SST,temporary)           
            IF (kpp_const_fields%L_SST_SMOOTH_ANOM) THEN
               allocate(SST_anom(NX_GLOBE,NY_GLOBE))
! Take SST anomaly from climatological SST provided
               SST_anom = temporary - SST_in(:,:,1)
            ELSE              
            ENDIF
            IF (kpp_const_fields%L_SST_SMOOTH .or.
     +           kpp_const_fields%L_SST_SMOOTH_ANOM) THEN
               allocate(SST_smooth(NX_GLOBE,NY_GLOBE))
               IF (kpp_const_fields%L_SST_SMOOTH_ANOM) THEN
                  WRITE(6,*) 'KPP: Smoothing SST anomaly'
                  CALL smooth_sst_out(SST_anom,kpp_3d_fields,
     +                 kpp_const_fields,SST_smooth)               
               ELSE IF(kpp_const_fields%L_SST_SMOOTH) THEN
                  WRITE(6,*) 'KPP: Smoothing SST'
                  CALL smooth_sst_out(temporary,kpp_3d_fields,
     +                 kpp_const_fields,SST_smooth)
               ENDIF
               IF (kpp_const_fields%L_SST_SMOOTH_ANOM) THEN                  
                  !WRITE(6,*) 'KPP: Smoothed SST anomaly = ',SST_smooth                  
                  ! Add smoothed anomaly (sst_smooth) to climatology (SST_in)
                  ! Make sure to also remove unsmoothed anomaly!
                  !WRITE(6,*) 'SST_in = ',SST_in(:,160,1)
                  !WRITE(6,*) 'SST_smooth = ',SST_smooth(:,160)
                  !WRITE(6,*) 'SST_anom = ',SST_anom(:,160)                  
                  IF (kpp_const_fields%L_SST_ANOM_FUDGE) THEN
                     ! DO NOT overwrite the SST passed to the coupler ("temporary")
                     SST_smooth = SST_in(:,:,1) + SST_smooth ! Use this for overwriting SST
                  ELSE
                     temporary = SST_in(:,:,1) + SST_smooth
                  ENDIF
                  deallocate(SST_anom)
               ELSE IF (kpp_const_fields%L_SST_SMOOTH) THEN               
                  temporary=SST_smooth
               ENDIF
               !deallocate(SST_smooth)
            ENDIF
            IF (L_OUTKELVIN) temporary=temporary+TK0
            !WRITE(6,*) 'temporary = ',temporary(:,160)
#endif
         CASE('OFRZN01')
#ifdef TOYCLIM
            temporary=ICE
#else
            CALL ONED_GLOBAL_TWOD_GLOBAL(ICE,temporary)
#endif
         CASE ('OSNWTN01')
#ifdef TOYCLIM
            temporary=SNOWDEPTH
#else          
            IF (kpp_const_fields%L_SST_ANOM_FUDGE) THEN
               DO ix=1,NX_GLOBE
                  DO jy=1,NY_GLOBE
                     ipoint_globe=jy*NX_GLOBE+ix
                     IF (ix .ge. ifirst .and. ix .le. ilast .and.
     +                    jy .ge. jfirst .and. jy .le. jlast) THEN                        
                        ipoint=(jy-jfirst)*nx+(ix-ifirst)+1
                        IF (kpp_3d_fields%cplwght(ipoint_globe).LT.1e-10 
     +                  .or. .not. kpp_3d_fields%L_OCEAN(ipoint)) THEN 
                           snowdepth(ipoint_globe)=0.0
                        ELSE
                           snowdepth(ipoint_globe)=TK0+
     +                          SST_smooth(ix,jy)*
     +                          kpp_3d_fields%cplwght(ipoint_globe)+
     +                          SST_in(ix,jy,1)*(1.0-
     +                          kpp_3d_fields%cplwght(ipoint_globe))
                           snowdepth(ipoint_globe)=
     +                          snowdepth(ipoint_globe)/1.0e6
                        ENDIF
                     ELSE
                        snowdepth(ipoint_globe)=0.0
                     ENDIF
                  ENDDO
               ENDDO
            ELSE
               snowdepth(:)=0.0
	    ENDIF
            CALL ONED_GLOBAL_TWOD_GLOBAL(SNOWDEPTH,temporary)
#endif            
         CASE ('OHICN01')
#ifdef TOYCLIM
            temporary=ICEDEPTH
#else
            CALL ONED_GLOBAL_TWOD_GLOBAL(ICEDEPTH,temporary)
#endif
         CASE ('SUNOCEAN')
#ifdef TOYCLIM
            temporary=SURF_CURR_X
#else
            IF (L_COUPLE_CURRENTS) THEN
               CALL ONED_GLOBAL_TWOD_GLOBAL(SURF_CURR_X,temporary)
            ELSE
               temporary=0.
            ENDIF
#endif
         CASE ('SVNOCEAN')
#ifdef TOYCLIM
            temporary=SURF_CURR_Y
#else
            IF (L_COUPLE_CURRENTS) THEN
               CALL ONED_GLOBAL_TWOD_GLOBAL(SURF_CURR_Y,temporary)
            ELSE
               temporary=0.
            ENDIF
#endif
         CASE DEFAULT
c
c     The user should never see this - if they do, it means that the model
c     is exporting a field that is not defined here (and likely not defined
c     in the "namcouple" file either).
c
            WRITE(nuout,*) 'KPP: Unexpected CASE DEFAULT for i=',i
         END SELECT
!         WRITE(nuout,*) 'KPP: Calling PRISM_Put_Proto ',
!     +        'for variable ',cl_writ(i)         
         CALL prism_put_proto(il_var_id_out(i),
     +        time_in_seconds,temporary,ierror)
         IF (ierror.NE.PRISM_Ok.and.ierror.LT.PRISM_Sent) THEN
            WRITE(nuout,*) 'KPP: Received error from ',
     +           'PRISM_Put_Proto =',ierror,' sending variable ',
     +           cl_writ(i),' at model time ',time_in_seconds,' sec.'
            WRITE(nuout,*) 'KPP: Aborting coupled integration ...'
            CALL prism_abort_proto(il_comp_id,'couple_io_oasis3.f',
     +           'send')
!         ELSE
!            WRITE(nuout,*) 'KPP: Successfully called ',
!     +           'PRISM_Put_Proto for variable ',cl_writ(i)
         ENDIF
      ENDDO
c
c     End
c
      RETURN
      END SUBROUTINE mpi1_oasis3_output

      SUBROUTINE mpi1_oasis3_terminate
c     
c     Allows KPP to terminate its role in the coupled integration
c     Essentially a wrapper around <prism_terminate_proto>
c
c     NPK 28/09/09 - R3
c     
c     Use OASIS3 (PRISM) modules - the compiler must be able to find
c     these via the Makefile.
c      
#ifdef OASIS3_MCT
      USE mod_prism
#else
      USE mod_kinds_model 
      USE mod_prism_proto
#endif
c
      IMPLICIT NONE
c
#include "kpp_oasis3.inc"
c
c     Local variables
c
      INTEGER ierror
c
c     Terminate the integration
c
      WRITE(6,*) 'KPP : Calling prism_terminate_proto(ierror)'
      CALL prism_terminate_proto(ierror)
      WRITE(6,*) 'KPP : Called prism_terminate_proto(ierror)'
      IF (ierror .NE. PRISM_Ok) THEN
         WRITE(il_mparout,*) 'KPP: Received error from ',
     +        'PRISM_Terminate_Proto =',ierror,
     +        ' when terminating model.'
      ENDIF
      WRITE(il_mparout,*) '---- End of the KPP integration ----'
      CLOSE(il_mparout)
      STOP
c     
c     End
c
      RETURN
      END SUBROUTINE mpi1_oasis3_terminate

      SUBROUTINE mpi1_oasis3_read_netcdf(solar,non_solar,PminusE,
     +     ustress,vstress)
      IMPLICIT NONE

#include "parameter.inc"
#include <netcdf.inc>
     
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
      REAL solar(NPTS),non_solar(NPTS),PminusE(NPTS),
     +     ustress(NPTS),vstress(NPTS)
      REAL*4 temp_in(NPTS)
      INTEGER ncid,status,dimid,varids(6)

      status=NF_OPEN('kpp_atmos_fields.nc',NF_NOWRITE,ncid)
      IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'KPP : Reading coupling fields from input file'
      
      status=NF_INQ_VARID(ncid,'solar',varids(1))
      status=NF_INQ_VARID(ncid,'non_solar',varids(2))
      status=NF_INQ_VARID(ncid,'PminusE',varids(3))
      status=NF_INQ_VARID(ncid,'ustress',varids(4))
      status=NF_INQ_VARID(ncid,'vstress',varids(5))

      status=NF_GET_VARA_REAL(ncid,varids(1),(/1/),(/NPTS/),
     +     temp_in)
      solar=temp_in
      status=NF_GET_VARA_REAL(ncid,varids(2),(/1/),(/NPTS/),
     +     temp_in)
      non_solar=temp_in
      status=NF_GET_VARA_REAL(ncid,varids(3),(/1/),(/NPTS/),
     +     temp_in)
      PminusE=temp_in
      status=NF_GET_VARA_REAL(ncid,varids(4),(/1/),(/NPTS/),
     +     temp_in)
      ustress=temp_in
      status=NF_GET_VARA_REAL(ncid,varids(5),(/1/),(/NPTS/),
     +     temp_in)
      vstress=temp_in

      status=NF_CLOSE(ncid)

      RETURN
      END SUBROUTINE mpi1_oasis3_read_netcdf      

      SUBROUTINE mpi1_oasis3_write_netcdf(solar,non_solar,PminusE,
     +     ustress,vstress)
      IMPLICIT NONE

#include "parameter.inc"
#include <netcdf.inc>
      
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
      REAL solar(NPTS),non_solar(NPTS),PminusE(NPTS),
     +     ustress(NPTS),vstress(NPTS)
      REAL*4 temp_out(NPTS)
      INTEGER ncid,status,dimid,varids(6)

      status=NF_CREATE('kpp_atmos_fields.nc',NF_CLOBBER,ncid)
      IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'KPP : Coupling fields output file created'
      
      status=NF_DEF_DIM(ncid,'point',NPTS,dimid)
      status=NF_DEF_VAR(ncid,'solar',NF_FLOAT,1,(/dimid/),varids(1))
      status=NF_DEF_VAR(ncid,'non_solar',NF_FLOAT,1,(/dimid/),
     +     varids(2))
      status=NF_DEF_VAR(ncid,'PminusE',NF_FLOAT,1,(/dimid/),varids(3))
      status=NF_DEF_VAR(ncid,'ustress',NF_FLOAT,1,(/dimid/),varids(4))
      status=NF_DEF_VAR(ncid,'vstress',NF_FLOAT,1,(/dimid/),varids(5))
      
      status=NF_ENDDEF(ncid)
      
      temp_out=solar
      status=NF_PUT_VARA_REAL(ncid,varids(1),temp_out)
      temp_out=non_solar
      status=NF_PUT_VARA_REAL(ncid,varids(2),temp_out)
      temp_out=PminusE
      status=NF_PUT_VARA_REAL(ncid,varids(3),temp_out)
      temp_out=ustress
      status=NF_PUT_VARA_REAL(ncid,varids(4),temp_out)
      temp_out=vstress
      status=NF_PUT_VARA_REAL(ncid,varids(5),temp_out)
      
      status=NF_CLOSE(ncid)

      RETURN
      END SUBROUTINE mpi1_oasis3_write_netcdf
            
      SUBROUTINE ONED_GLOBAL_ONED_REGIONAL(global,regional)
c
c     Transforms a one-dimensional global field to a one-dimensional
c     regional field.
c     NPK 19/9/09 - R3
c     
      IMPLICIT NONE

#include "parameter.inc"
#include "couple.com"
      
      REAL global(NPTS_GLOBE),regional(NPTS)
      INTEGER*4 ix,jy,ipoint,ipoint_globe

      DO jy=jfirst,jlast
         DO ix=ifirst,ilast
            ipoint_globe=jy*NX_GLOBE+ix
            ipoint=(jy-jfirst)*nx+(ix-ifirst)+1
            regional(ipoint)=global(ipoint_globe)
         ENDDO
      ENDDO
         
      RETURN
      END

      SUBROUTINE TWOD_GLOBAL_ONED_REGIONAL(global,regional)
c
c     Transforms a two-dimensional global field (e.g., one received
c     from HadGEM3-A) to a one-dimensional regional field (as required
c     for KPP).
c     NPK 28/9/09 - R3
c
      IMPLICIT NONE
#include "parameter.inc"
#include "couple.com"
      
      REAL global(NX_GLOBE,NY_GLOBE),regional(NPTS)
      INTEGER*4 ix,jy,ipoint,ipoint_globe

      DO jy=jfirst,jlast
         DO ix=ifirst,ilast
            ipoint=(jy-jfirst)*NX+(ix-ifirst)+1
            regional(ipoint)=global(ix,jy)
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE ONED_GLOBAL_TWOD_GLOBAL(oned,twod)
      IMPLICIT NONE
c
c     Transforms a one-dimensional global field (as produced by
c     the combination of a KPP regional field and a global climatology)
c     to a two-dimensional global field (as required by HadGEM3-A).
c     NPK 28/9/09 - R3
c
#include "parameter.inc"

      REAL oned(NX_GLOBE*NY_GLOBE),twod(NX_GLOBE,NY_GLOBE)
      INTEGER*4 ix,jy,ipoint_globe

      DO jy=1,NY_GLOBE
         DO ix=1,NX_GLOBE
            ipoint_globe=(jy-1)*NX_GLOBE+ix
            twod(ix,jy)=oned(ipoint_globe)
         ENDDO
      ENDDO

      RETURN
      END      

      SUBROUTINE smooth_sst_out(sst_in,kpp_3d_fields,kpp_const_fields,
     +     sst_out)
      IMPLICIT NONE

#include "kpp_3d_type.com"

      TYPE(kpp_const_type), intent(in) :: kpp_const_fields
      TYPE(kpp_3d_type),intent(in) :: kpp_3d_fields
      REAL sst_in(NX_GLOBE,NY_GLOBE), sst_out(NX_GLOBE,NY_GLOBE)
      REAL,allocatable :: sst_smooth(:,:)
      REAL weight
      INTEGER ix,my_ix,jy,my_jy,ifirst,ilast,jfirst,jlast,blend,
     +     ipoint_globe,my_npts

      ifirst = kpp_const_fields%sst_smooth_ifirst
      ilast = kpp_const_fields%sst_smooth_ilast
      jfirst = kpp_const_fields%sst_smooth_jfirst
      jlast = kpp_const_fields%sst_smooth_jlast
      blend = kpp_const_fields%sst_smooth_blend
      
      sst_out = sst_in
      allocate(sst_smooth(ifirst:ilast,jfirst:jlast))
      sst_smooth(:,:) = 0
      IF (kpp_const_fields%L_SST_SMOOTH_X .and. .not. 
     +     kpp_const_fields%L_SST_SMOOTH_Y) THEN
!     Smooth in X.  Need a separate smoothed value at each Y point (mean over all X).
         WRITE(6,*) 'KPP: Smooth SST in X between ',jfirst,' and ',jlast
     +        ,' and between ',ifirst,' and ',ilast         
         DO jy=jfirst,jlast
            my_npts = 0
            DO ix=ifirst,ilast
               ipoint_globe = (jy-1)*NX_GLOBE+ix
               IF (kpp_3d_fields%cplwght(ipoint_globe) .gt. 0) THEN
                  sst_smooth(ifirst,jy) = sst_smooth(ifirst,jy) + 
     +                 sst_in(ix,jy)
                  my_npts = my_npts+1
               ENDIF
            ENDDO
            IF (my_npts .gt. 0) THEN 
               sst_smooth(:,jy) = sst_smooth(ifirst,jy)/FLOAT(my_npts)
            ELSE
               sst_smooth(:,jy) = -999.0
            ENDIF
         ENDDO
         IF (blend .gt. 0) 
     +        CALL smooth_blended_sst(sst_in,sst_smooth,sst_out,
     +        kpp_3d_fields,kpp_const_fields,ifirst,ilast,jfirst,jlast)  
      ELSE IF (kpp_const_fields%L_SST_SMOOTH_Y .and. .not. 
     +        kpp_const_fields%L_SST_SMOOTH_X) THEN
!     Smooth in Y. Need a separate smoothed value at each X point (mean over all Y).
         WRITE(6,*) 'KPP: Smooth SST in Y between ',jfirst,' and ',jlast
     +        ,' and between ',ifirst,' and ',ilast
         DO ix=ifirst,ilast
            sst_smooth(ix,:) = 0
            my_npts=0
            DO jy=jfirst,jlast
               ipoint_globe = (jy-1)*NX_GLOBE+ix
               IF (kpp_3d_fields%cplwght(ipoint_globe) .gt. 0) THEN          
                  sst_smooth(ix,jfirst) = sst_smooth(ix,jfirst) + 
     +                 sst_in(ix,jy)
                  my_npts = my_npts+1
               ENDIF
            ENDDO
            IF (my_npts .gt. 0) THEN
               sst_smooth(ix,:) = sst_smooth(ix,jfirst)/FLOAT(my_npts)
            ELSE
               sst_smooth(ix,:) = -999.0
            ENDIF
         ENDDO         
         IF (blend .gt. 0)
     +        CALL smooth_blended_sst(sst_in,sst_smooth,sst_out,
     +        kpp_3d_fields,kpp_const_fields,ifirst,ilast,jfirst,jlast)
      ELSE IF (kpp_const_fields%L_SST_SMOOTH_Y .and. 
     +        kpp_const_fields%L_SST_SMOOTH_X) THEN         
         WRITE(6,*) 'KPP: Smooth SST in X and Y between ',jfirst,' and '
     +        ,jlast,' and between ',ifirst,' and ',ilast
!     Smooth in both X and Y. Need one value.
         my_npts=0         
         DO ix=ifirst,ilast
            DO jy=jfirst,jlast
               ipoint_globe = (jy-1)*NX_GLOBE+ix
               IF (kpp_3d_fields%cplwght(ipoint_globe) .gt. 0) THEN
                  sst_smooth(ifirst,jfirst) = sst_smooth(ifirst,jfirst)+
     +                 sst_in(ix,jy)
                  my_npts = my_npts+1
               ENDIF
            ENDDO            
         ENDDO
         IF (my_npts .gt. 0) THEN
            sst_smooth(:,:) = sst_smooth(ifirst,jfirst)/FLOAT(my_npts)
         ELSE
            sst_smooth(:,:) = -999.0
         ENDIF
         IF (blend .gt. 0) 
     +        CALL smooth_blended_sst(sst_in,sst_smooth,sst_out,
     +        kpp_3d_fields,kpp_const_fields,ifirst,ilast,jfirst,jlast)    
      ENDIF
      deallocate(sst_smooth)
            
      RETURN
      END SUBROUTINE smooth_sst_out

      SUBROUTINE smooth_blended_sst(sst_in,sst_smooth,sst_out,
     +     kpp_3d_fields,kpp_const_fields,ifirst,ilast,jfirst,jlast)
      IMPLICIT NONE

#include "kpp_3d_type.com"

      REAL, intent(in) :: sst_in(NX_GLOBE,NY_GLOBE)
      REAL, intent(in) :: sst_smooth(ifirst:ilast,jfirst:jlast)
      REAL, intent(out) :: sst_out(NX_GLOBE,NY_GLOBE)
      REAL :: weight(NX_GLOBE,NY_GLOBE),sst_tmp
      TYPE(kpp_3d_type),intent(in) :: kpp_3d_fields
      TYPE(kpp_const_type), intent(in) :: kpp_const_fields
      INTEGER :: blend, ix,jy,my_ix,my_jy,ipoint_globe
      INTEGER, intent(in) :: ifirst,ilast,jfirst,jlast

      blend = kpp_const_fields%sst_smooth_blend

      sst_out = sst_in
      
      weight(:,:)=0.0
      weight(ifirst:ilast,jfirst:jlast) = 1.0
      
! West boundary
      DO ix=ifirst-blend,ifirst
         IF (ix .lt. 1) THEN
            my_ix = ix + NX_GLOBE
         ELSE IF (ix .gt. NX_GLOBE) THEN
            my_ix = ix - NX_GLOBE
         ELSE
            my_ix = ix
         ENDIF
         DO jy=jfirst,jlast
            weight(ix,jy) = MAX(weight(ix,jy),
     +           1.0-ABS(ix-ifirst)/FLOAT(blend))
         ENDDO
      ENDDO
! East boundary
      DO ix=ilast,ilast+blend
         IF (ix .lt. 1) THEN
            my_ix = ix + NX_GLOBE
         ELSE IF (ix .gt. NX_GLOBE) THEN
            my_ix = ix - NX_GLOBE
         ELSE
            my_ix = ix
         ENDIF
         DO jy=jfirst,jlast
            weight(ix,jy) = MAX(weight(ix,jy),
     +           1.0-ABS(ix-ilast)/FLOAT(blend))
         ENDDO
      ENDDO
! South boundary
      DO jy=jfirst-blend,jfirst
         IF (jy .lt. 1 .or. jy .gt. NY_GLOBE) THEN
            WRITE(0,*) 'KPP: Blending region for smoothed SST is', 
     +           ' outside global domain. Abort !'
            CALL MIXED_ABORT
         ELSE
            DO ix=ifirst,ilast
               weight(ix,jy) = MAX(weight(ix,jy),
     +              1.0-ABS(jy-jfirst)/FLOAT(blend))
            ENDDO
         ENDIF               
      ENDDO
! North boundary
      DO jy=jlast,jlast+blend
         IF (jy .lt. 1 .or. jy .gt. NY_GLOBE) THEN
            WRITE(0,*) 'KPP: Blending region for smoothed SST is', 
     +           ' outside global domain. Abort !'
            CALL MIXED_ABORT
         ELSE
            DO ix=ifirst,ilast
               weight(ix,jy) = MAX(weight(ix,jy),
     +              1.0-ABS(jy-jlast)/FLOAT(blend))
            ENDDO
         ENDIF               
      ENDDO

      DO ix=ifirst-blend,ifirst
! Southwest corner         
         DO jy=jfirst-blend,jfirst
            weight(ix,jy) = MIN(weight(ifirst,jy),
     +           weight(ix,jfirst))
         ENDDO
! Northwest corner         
         DO jy=jlast,jlast+blend
            weight(ix,jy) = MIN(weight(ifirst,jy),
     +           weight(ix,jlast))
         ENDDO
      ENDDO
      DO ix=ilast,ilast+blend
! Southeast corner
         DO jy=jfirst-blend,jfirst
            weight(ix,jy) = MIN(weight(ilast,jy),
     +           weight(ix,jfirst))
         ENDDO
! Northeast corner
         DO jy=jlast,jlast+blend
            weight(ix,jy) = MIN(weight(ilast,jy),
     +           weight(ix,jlast))
         ENDDO
      ENDDO

      DO ix=1,NX_GLOBE
         DO jy=1,NY_GLOBE
            IF (weight(ix,jy) .eq. 0) THEN
               sst_out(ix,jy) = sst_in(ix,jy)
            ELSE IF (weight(ix,jy) .eq. 1) THEN
               sst_out(ix,jy) = sst_smooth(ix,jy)
            ELSE IF (weight(ix,jy) .gt. 0 .and.
     +              weight(ix,jy) .lt. 1) THEN
               my_ix = MAX(ifirst,ix)
               my_ix = MIN(my_ix,ilast)
               my_jy = MAX(jfirst,jy)
               my_jy = MIN(my_jy,jlast)
!               WRITE(6,*) 'KPP: Smooth SST at ',ix,',',jy,
!     +              ' using ',my_ix,',',my_jy 
               sst_tmp = sst_smooth(my_ix,my_jy)
               ipoint_globe = (jy-1)*NX_GLOBE+ix
               IF (kpp_3d_fields%cplwght(ipoint_globe) .gt. 0 .and.
     +              sst_tmp .gt. -100 .and. sst_tmp .lt. 1000) THEN
!     WRITE(6,*) 'KPP: ',sst_tmp
                  sst_out(ix,jy) = weight(ix,jy)*sst_tmp +
     +                 (1.0-weight(ix,jy))*sst_in(ix,jy)
!                  WRITE(6,*) 'KPP: ',sst_tmp,sst_in(ix,jy),
!     +                 sst_out(ix,jy),weight(ix,jy)
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE smooth_blended_sst
      
      
      
#endif /*OASIS3*/
#endif /*COUPLE*/
