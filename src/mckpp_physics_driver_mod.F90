MODULE mckpp_physics_driver_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_1d_type, kpp_const_fields
  USE mckpp_timer, ONLY: mckpp_define_new_timer, mckpp_start_timer, mckpp_stop_timer
  USE mckpp_parameters
  USE mckpp_physics_ocnstep_mod, ONLY: mckpp_physics_ocnstep
  USE mckpp_physics_overrides, ONLY: mckpp_physics_overrides_bottomtemp, mckpp_physics_overrides_check_profile
  USE mckpp_types_transfer, ONLY: mckpp_fields_3dto1d, mckpp_fields_1dto3d
  
  IMPLICIT NONE

CONTAINS

SUBROUTINE mckpp_physics_driver()
  
  ! Local
  TYPE(kpp_1d_type) :: kpp_1d_fields
  INTEGER :: ipt, index
#ifdef OPENMP
  INTEGER :: tid, OMP_GET_THREAD_NUM
#endif
  CHARACTER(LEN=21) phys_timer_name
  CHARACTER(LEN=19) trans_timer_name
  
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP SHARED(kpp_3d_fields, kpp_const_fields) &
!$OMP SHARED(nz, nzp1, nx, ny, npts, nvel, nsclr, nvp1, nsp1, itermax) &
!$OMP SHARED(hmixtolfrac, nztmax, nzp1tmax, nsflxs, njdt, maxmodeadv) &
!$OMP PRIVATE(trans_timer_name, phys_timer_name, tid, index, kpp_1d_fields)
  tid=OMP_GET_THREAD_NUM()
  WRITE(trans_timer_name,'(A17,I2)') 'KPP 3D/1D thread ',tid
  WRITE(phys_timer_name,'(A19,I2)') 'KPP Physics thread ',tid
  IF (kpp_const_fields%ntime .EQ. 1) THEN
!$OMP CRITICAL
    CALL mckpp_define_new_timer(trans_timer_name, index)
    CALL mckpp_define_new_timer(phys_timer_name, index)
!$OMP END CRITICAL
  END IF 
!$OMP DO SCHEDULE(dynamic)
#else
  WRITE(trans_timer_name,'(A19)') 'KPP 3D/1D thread 01'
  WRITE(phys_timer_name,'(A21)') 'KPP Physics thread 01'
#endif /*OPENMP*/
  DO ipt=1,npts
    IF (kpp_3d_fields%run_physics(ipt)) THEN

        CALL mckpp_start_timer(trans_timer_name)
        CALL mckpp_fields_3dto1d(kpp_3d_fields,ipt,kpp_1d_fields)
        CALL mckpp_stop_timer(trans_timer_name) 

        CALL mckpp_start_timer(phys_timer_name)
        CALL mckpp_physics_ocnstep(kpp_1d_fields,kpp_const_fields)
        CALL mckpp_physics_overrides_check_profile(kpp_1d_fields,kpp_const_fields)
        CALL mckpp_stop_timer(phys_timer_name)

        CALL mckpp_start_timer(trans_timer_name)
        CALL mckpp_fields_1dto3d(kpp_1d_fields,ipt,kpp_3d_fields)
        CALL mckpp_stop_timer(trans_timer_name)
        
     ENDIF
   ENDDO
!$OMP END DO
!$OMP END PARALLEL
  
  IF (kpp_const_fields%L_VARY_BOTTOM_TEMP) THEN
     CALL mckpp_start_timer("Update ancillaries")
     CALL MCKPP_PHYSICS_OVERRIDES_BOTTOMTEMP()
     CALL mckpp_stop_timer("Update ancillaries")
  ENDIF

END SUBROUTINE mckpp_physics_driver

END MODULE mckpp_physics_driver_mod
