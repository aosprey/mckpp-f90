#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>

SUBROUTINE MCKPP_COUPLING_CAM3_STEP(srf_state,srfflx)

  USE shr_kind_mod, only: r8=>shr_kind_r8
  USE ppgrid,       only: pcols,begchunk,endchunk
  USE comsrf,       only: surface_state,srfflx_parm,landfrac,asdirocn,aldirocn,asdifocn,aldifocn,tsocn
  USE phys_grid,    only: get_ncols_p,get_rlat_all_p,get_rlon_all_p
  USE time_manager, only: get_step_size,get_curr_calday,get_nstep,is_end_curr_day
  USE physconst,    only: latvap
  USE mckpp_types,  only: kpp_3d_fields,kpp_1d_fields,kpp_const_fields
  USE pmgrid,       only: masterproc
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len

  IMPLICIT NONE

  TYPE(surface_state),intent(inout) :: srf_state(begchunk:endchunk)
  TYPE(srfflx_parm),intent(inout) :: srfflx(begchunk:endchunk)

#include <comctl.h>

  REAL(r8) :: sstk(pcols),sst(pcols),ltheat(pcols),shflx(pcols),lhflx(pcols),&
       lwup(pcols),tref(pcols),taux(pcols),tauy(pcols)
  INTEGER :: ncol,ichnk,icol,j,k,indx(pcols),my_npts
  
  CHARACTER(LEN=24) :: routine = "MCKPP_COUPLING_CAM3_STEP"
  CHARACTER(LEN=max_message_len) :: message
    
  ! Increment model time
  kpp_const_fields%ntime=kpp_const_fields%ntime+1
  kpp_const_fields%time=kpp_const_fields%startt+(kpp_const_fields%ntime-1)*&
       kpp_const_fields%dto/kpp_const_fields%spd
  
  WRITE(message,*) 'Running for timestep ',kpp_const_fields%ntime,' time = ',kpp_const_fields%time
  CALL mckpp_print(routine, message)

  ! Update MC-KPP boundary conditions if not first timestep
  IF (kpp_const_fields%ntime .ne. 1) THEN 
     ! CALL mckpp_print(routine,'Calling MCKPP_BOUNDARY_UPDATE')
     CALL MCKPP_BOUNDARY_UPDATE
     ! CALL mckpp_print(routine, 'Returned from MCKPP_BOUNDARY_UPDATE')
  ENDIF


  ! Get CAM3 fluxes
  ! Removed OPENMP loop here because it looked dodgy.

  ! CALL mckpp_print(routine, 'Getting fluxes from CAM atmosphere')
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)

     my_npts=0
     DO icol=1,ncol
        IF (landfrac(icol,ichnk) .lt. 1) THEN !IF (kpp_3d_fields(ichnk)%L_OCEAN(icol) .and. kpp_3d_fields(ichnk)%cplwght(icol) .gt. 0.0) THEN
           my_npts=my_npts+1
           indx(my_npts) = icol
        ELSE
           sst(icol)=0.
        ENDIF
     ENDDO
     
     ltheat(1:ncol)=latvap

     !kpp_3d_fields(ichnk)%sflux(:,:,5,0)=0.0

     IF (my_npts .gt. 0) THEN
       ! WRITE(message,*) 'KPP SST = ',kpp_3d_fields(ichnk)%X(1:ncol,1,1)
       ! CALL mckpp_print(routine, message)
        
        IF (kpp_const_fields%ntime .ne. 1) THEN
           ! If not first timestep, use mean value over coupling period
           sstk(1:ncol) = tsocn(1:ncol,ichnk)+kpp_const_fields%TK0
           ! WRITE(message,*) sstk(1:ncol)
           ! CALL mckpp_print(routine, message)
           ! sstk(1:ncol) = kpp_3d_fields(ichnk)%sst_cpl(1:ncol)+kpp_const_fields%TK0
        ELSE 
           ! If first timestep, use initial value
           sstk(1:ncol) = kpp_3d_fields(ichnk)%X(1:ncol,1,1)+kpp_const_fields%TK0
        ENDIF
        ! Delete when means over coupling period work correctly
        ! sstk(1:ncol) = kpp_3d_fields(ichnk)%X(1:ncol,1,1)+kpp_const_fields%TK0

        CALL FLXOCE(indx,my_npts,srf_state(ichnk)%pbot,srf_state(ichnk)%ubot,&
             srf_state(ichnk)%vbot,srf_state(ichnk)%tbot,srf_state(ichnk)%qbot,&
             srf_state(ichnk)%thbot,srf_state(ichnk)%zbot,sstk,ltheat,&
             srfflx(ichnk)%shf,srfflx(ichnk)%lhf,srfflx(ichnk)%wsx,srfflx(ichnk)%wsy,&
             srfflx(ichnk)%lwup,srfflx(ichnk)%tref)
        
        DO j=1,my_npts
           icol=indx(j)
                      
           ! Store surface fluxes in temporary variables to allow for
           ! coupling timesteps longer than the model timestep
           
           ! Zonal wind stress
           kpp_3d_fields(ichnk)%sflux_cpl(icol,1) = srfflx(ichnk)%wsx(icol)/FLOAT(kpp_const_fields%ndtocn) + &
                kpp_3d_fields(ichnk)%sflux_cpl(icol,1) 
           
           ! Meridional wind stress
           kpp_3d_fields(ichnk)%sflux_cpl(icol,2) = srfflx(ichnk)%wsy(icol)/FLOAT(kpp_const_fields%ndtocn) + &
                kpp_3d_fields(ichnk)%sflux_cpl(icol,2)
           
           ! Net downward shortwave flux
           kpp_3d_fields(ichnk)%sflux_cpl(icol,3) = (srf_state(ichnk)%sols(icol)*(1.0_r8-asdirocn(icol,ichnk))+&
                srf_state(ichnk)%solsd(icol)*(1.0_r8-asdifocn(icol,ichnk))+&
                srf_state(ichnk)%soll(icol)*(1.0_r8-aldirocn(icol,ichnk))+&
                srf_state(ichnk)%solld(icol)*(1.0_r8-aldifocn(icol,ichnk))) / FLOAT(kpp_const_fields%ndtocn) + &
                kpp_3d_fields(ichnk)%sflux_cpl(icol,3)
           
           ! Net downward longwave flux
           kpp_3d_fields(ichnk)%sflux_cpl(icol,4) = (srf_state(ichnk)%flwds(icol)-&
                srfflx(ichnk)%lwup(icol)-srfflx(ichnk)%lhf(icol)-srfflx(ichnk)%shf(icol)) / &
                FLOAT(kpp_const_fields%ndtocn) + kpp_3d_fields(ichnk)%sflux_cpl(icol,4)
           
           ! Melting of sea ice = 0.0
           kpp_3d_fields(ichnk)%sflux_cpl(icol,5) = 0.0_r8
           
           ! Precipitation minus evaporation
           kpp_3d_fields(ichnk)%sflux_cpl(icol,6) = ((srf_state(ichnk)%precl(icol)+srf_state(ichnk)%precc(icol)+&
                srf_state(ichnk)%precsc(icol)+srf_state(ichnk)%precsl(icol))*1000.-&
                srfflx(ichnk)%lhf(icol)/kpp_const_fields%EL) / FLOAT(kpp_const_fields%ndtocn) + &
                kpp_3d_fields(ichnk)%sflux_cpl(icol,6)
           
           IF (kpp_const_fields%ntime .ne. 1) THEN
              ! Update values in MC-KPP only if this is a coupling timestep
              IF (MOD(kpp_const_fields%ntime,kpp_const_fields%ndtocn) .eq. 0) THEN
                 kpp_3d_fields(ichnk)%sflux(icol,:,5,0)=kpp_3d_fields(ichnk)%sflux_cpl(icol,:)
                 ! Reset means to zero
                 kpp_3d_fields(ichnk)%sflux_cpl(icol,:)=0.
              ENDIF
           ELSE        
              ! Use instantaneous values for the first timestep
              
              ! Zonal wind stress
              kpp_3d_fields(ichnk)%sflux(icol,1,5,0)=srfflx(ichnk)%wsx(icol)
              ! Meridional wind stress
              kpp_3d_fields(ichnk)%sflux(icol,2,5,0)=srfflx(ichnk)%wsy(icol)             
              ! Net downward shortwave flux, using albedos from CAM (positive down)
              kpp_3d_fields(ichnk)%sflux(icol,3,5,0)=srf_state(ichnk)%sols(icol)*(1.0_r8-asdirocn(icol,ichnk))+&
                   srf_state(ichnk)%solsd(icol)*(1.0_r8-asdifocn(icol,ichnk))+&
                   srf_state(ichnk)%soll(icol)*(1.0_r8-aldirocn(icol,ichnk))+&
                   srf_state(ichnk)%solld(icol)*(1.0_r8-aldifocn(icol,ichnk)) 
              ! Net downward longwave flux
              kpp_3d_fields(ichnk)%sflux(icol,4,5,0)=srf_state(ichnk)%flwds(icol)-&
                   srfflx(ichnk)%lwup(icol)-srfflx(ichnk)%lhf(icol)-srfflx(ichnk)%shf(icol)
              ! Melting of sea ice = 0.0
              kpp_3d_fields(ichnk)%sflux(icol,5,5,0)=0.0_r8

              ! Precipitation minus evaporation
              kpp_3d_fields(ichnk)%sflux(icol,6,5,0)=(srf_state(ichnk)%precl(icol)+srf_state(ichnk)%precc(icol)+&
                   srf_state(ichnk)%precsc(icol)+srf_state(ichnk)%precsl(icol))*1000.-&
                   srfflx(ichnk)%lhf(icol)/kpp_const_fields%EL
           ENDIF
           
        ENDDO        
        !STOP
     ENDIF
  ENDDO
  ! CALL mckpp_print(routine, 'Done getting fluxes from CAM atmosphere')

  ! Call MC-KPP physics
  ! CALL mckpp_print(routine, 'Calling MCKPP_PHYSICS_DRIVER')
  CALL MCKPP_PHYSICS_DRIVER
  ! CALL mckpp_print(routine, 'Returned from MCKPP_PHYSICS_DRIVER')

  ! Write MC-KPP output if necessary
  ! CALL mckpp_print(routine, 'Calling MCKPP_OUTPUT_CONTROL')
  CALL MCKPP_OUTPUT_CONTROL
  ! CALL mckpp_print(routine, 'Returned from MCKPP_OUTPUT_CONTROL')

  ! Write MC-KPP checkpoint files if necessary
  CALL MCKPP_RESTART_CONTROL

  ! Output SST and ice back to CAM
  ! CALL mckpp_print(routine, 'Sending SST and ice back to CAM atmosphere')
  CALL MCKPP_COUPLING_CAM3_OUTPUT(srfflx)
  WRITE(message,*) 'Finished coupling for timestep ',kpp_const_fields%ntime
  CALL mckpp_print(routine, message)

  RETURN
END SUBROUTINE MCKPP_COUPLING_CAM3_STEP

SUBROUTINE MCKPP_COUPLING_CAM3_OUTPUT(srfflx)
  USE shr_kind_mod, only: r8=>shr_kind_r8
  USE mckpp_types, only: kpp_3d_fields,kpp_const_fields
  USE comsrf, only: tsocn, icefrac, sicthk, ocnfrac, landfrac, aice, frzmlt, lwupocn, srfflx_parm, Focn
  USE ice_constants, only: Tffresh, tfrez, rhow, cp_ocn
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid, only: get_ncols_p
  USE physconst, only: latvap
  USE constituents, only: pcnst,pnats
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len

  IMPLICIT NONE

  TYPE(srfflx_parm) :: srfflx(begchunk:endchunk)
  INTEGER :: j,k,ichnk,ncol,icol,my_npts,indx(PCOLS)
  REAL(r8) :: sst_out(PCOLS),& ! SST in deg C internally, converted to K before output as Tsocn
       ltheat(PCOLS),blended_sst(PCOLS)

  CHARACTER(LEN=26) :: routine = "MCKPP_COUPLING_CAM3_OUTPUT"
  CHARACTER(LEN=max_message_len) :: message
    
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     
     ltheat(1:ncol)=latvap
     my_npts=0
     DO icol=1,ncol
        IF (landfrac(icol,ichnk) .lt. 1) THEN 
           my_npts=my_npts+1
           indx(my_npts) = icol
           Focn(icol,ichnk)=0.
           frzmlt(icol,ichnk)=0.
        ELSE
           sst_out(icol)=0.
           Focn(icol,ichnk)=0.
           frzmlt(icol,ichnk)=0.
        ENDIF
     ENDDO

     IF (my_npts .gt. 0) THEN
        DO j=1,my_npts
           icol=indx(j)

           IF (kpp_3d_fields(ichnk)%cplwght(icol) .gt. 0 .and. kpp_3d_fields(ichnk)%cplwght(icol) .lt. 1) THEN
              blended_sst(icol)=kpp_3d_fields(ichnk)%X(icol,1,1)*kpp_3d_fields(ichnk)%cplwght(icol)+&
                   (kpp_3d_fields(ichnk)%sst(icol))*(1-kpp_3d_fields(ichnk)%cplwght(icol))
           ELSE IF (kpp_3d_fields(ichnk)%cplwght(icol) .eq. 1) THEN
              blended_sst(icol)=kpp_3d_fields(ichnk)%X(icol,1,1)
           ELSE IF (kpp_3d_fields(ichnk)%cplwght(icol) .eq. 0) THEN
              blended_sst(icol)=kpp_3d_fields(ichnk)%sst(icol)
           ELSE
              WRITE(message,*) 'Invalid value of coupling weight at ichnk = ',ichnk, &
                  'icol = ',icol,' coupling weight = ',kpp_3d_fields(ichnk)%cplwght(icol)
              CALL mckpp_print_warning(routine, message) 
           ENDIF
           
           ! Put MC-KPP ice fraction to CAM ice fraction
           icefrac(icol,ichnk)=kpp_3d_fields(ichnk)%iceconc(icol)     
           ! Make arbritrary ice thickness?
           IF (icefrac(icol,ichnk) .gt. 0.001) THEN
              sicthk(icol,ichnk)=2.0
           ELSE
              sicthk(icol,ichnk)=0.0
              icefrac(icol,ichnk)=0.0
           ENDIF

           !frzmlt(ichnk,icol)=0.
           !frzmlt(ichnk,icol) = (Tfrez-sst(icol))*(rhow*cp_ocn*kpp_3d_fields(ichnk)%hmix(icol))/kpp_const_fields%dtsec
           blended_sst(icol)=MAX(blended_sst(icol),Tfrez)
                      
           ! Store SST in a temporary variable to allow for a coupling timestep
           ! longer than the model timestep
           kpp_3d_fields(ichnk)%sst_cpl(icol) = blended_sst(icol) / FLOAT(kpp_const_fields%ndtocn) + &
                kpp_3d_fields(ichnk)%sst_cpl(icol)
           IF (kpp_const_fields%ntime .ne. 1) THEN
              IF (MOD(kpp_const_fields%ntime,kpp_const_fields%ndtocn) .eq. 0) THEN
                 ! If this is a coupling timestep, use mean SST in output variable
                 sst_out(icol)=kpp_3d_fields(ichnk)%sst_cpl(icol)
                 ! Reset coupling SST to zero
                 kpp_3d_fields(ichnk)%sst_cpl(icol)=0.
              ELSE
                 ! If not a coupling timestep, use previous value of SST
                 sst_out(icol)=tsocn(icol,ichnk)
              ENDIF
           ELSE
              ! Use instantaneous value for the first timestep
              sst_out(icol)=blended_sst(icol)
           ENDIF

           ! Set constituent water flux to be equal to -1*evap
           srfflx(ichnk)%cflx(icol,1)=srfflx(ichnk)%lhf(icol)/ltheat(icol)
           ! Set non-water constituent fluxes to zero
           DO k=2,pcnst+pnats
              srfflx(ichnk)%cflx(icol,k)=0.
           ENDDO
        ENDDO
        !lwupocn(1:ncol,ichnk)=kpp_3d_fields(ichnk)%sflux(1:ncol,4,5,0)
        lwupocn(1:ncol,ichnk)=srfflx(ichnk)%lwup(1:ncol)
     ENDIF

     ! Update SST in CAM
     srfflx(ichnk)%ts(1:ncol)=sst_out(1:ncol)+kpp_const_fields%TK0 ! Store in Kelvin    
     tsocn(1:ncol,ichnk)=sst_out(1:ncol)                           ! Store in Celsius
        
     ! Code copied from inidat.F90
     where (icefrac(:ncol,ichnk) + landfrac(:ncol,ichnk) > 1.0)
        icefrac(:ncol,ichnk) = 1.-landfrac(:ncol,ichnk)
     end where
     where (landfrac(:ncol,ichnk) < 1.0)
        aice(:ncol,ichnk) = icefrac(:ncol,ichnk)/(1. - landfrac(:ncol,ichnk))
     elsewhere
        aice(:ncol,ichnk) = 0.
     end where
     ocnfrac(:ncol,ichnk) = 1.-landfrac(:ncol,ichnk)-icefrac(:ncol,ichnk)
     

  ENDDO
  
  RETURN
END SUBROUTINE MCKPP_COUPLING_CAM3_OUTPUT
#else
SUBROUTINE MCKPP_COUPLING_CAM3_OUTPUT
END SUBROUTINE MCKPP_COUPLING_CAM3_OUTPUT
#endif
