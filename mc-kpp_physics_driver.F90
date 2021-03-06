#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE mckpp_physics_driver
  USE mckpp_parameters
  USE mckpp_types, only: kpp_3d_fields,kpp_1d_type,kpp_const_fields
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid,only: get_ncols_p
#else
SUBROUTINE mckpp_physics_driver(kpp_3d_fields,kpp_const_fields)
  USE mckpp_data_types
  USE mckpp_timer
#endif
  IMPLICIT NONE

#ifdef MCKPP_CAM3
  INTEGER :: icol,ncol,ichnk
#else
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
#endif
  
  ! Local
  TYPE(kpp_1d_type) :: kpp_1d_fields
  INTEGER :: ipt, index
#ifdef OPENMP
  INTEGER :: tid,OMP_GET_THREAD_NUM
#endif
  CHARACTER(LEN=21) phys_timer_name
  CHARACTER(LEN=19) trans_timer_name
  
#ifdef MCKPP_CAM3
  !WRITE(6,*) 'Before ocnstep, U = ',kpp_3d_fields(begchunk)%U(1:ncol,1,1)
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     DO icol=1,ncol
        IF (kpp_3d_fields(ichnk)%landfrac(icol) .lt. 1.0 .and. &
             kpp_3d_fields(ichnk)%cplwght(icol) .gt. 0.0) THEN
           CALL MCKPP_FIELDS_3Dto1D(kpp_3d_fields(ichnk),icol,kpp_1d_fields)
           CALL MCKPP_PHYSICS_OCNSTEP(kpp_1d_fields,kpp_const_fields)
           CALL MCKPP_PHYSICS_OVERRIDES_CHECK_PROFILE(kpp_1d_fields,kpp_const_fields)
           CALL MCKPP_FIELDS_1Dto3D(kpp_1d_fields,icol,kpp_3d_fields(ichnk))
        ENDIF
     ENDDO
  ENDDO   
  !WRITE(6,*) 'After ocnstep, U = ',kpp_3d_fields(begchunk)%U(1:ncol,1,1)
#else
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
#ifdef MCKPP_COUPLE
        IF (kpp_3d_fields%L_OCEAN(ipt) .and. kpp_3d_fields%cplwght(ipt) .gt. 0.0) THEN
#else       
        IF (kpp_3d_fields%L_OCEAN(ipt)) THEN
#endif
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
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif /*OPENMP*/
#endif /*MCKPP_CAM3*/
  
  IF (kpp_const_fields%L_VARY_BOTTOM_TEMP) THEN
     CALL mckpp_start_timer("Update ancillaries")
#ifdef MCKPP_CAM3
     CALL MCKPP_PHYSICS_OVERRIDES_BOTTOMTEMP
#else
     CALL MCKPP_PHYSICS_OVERRIDES_BOTTOMTEMP(kpp_3d_fields,kpp_const_fields)
#endif
     CALL mckpp_stop_timer("Update ancillaries")
  ENDIF

  RETURN
END SUBROUTINE mckpp_physics_driver
