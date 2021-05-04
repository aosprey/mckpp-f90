SUBROUTINE mckpp_physics_overrides_bottomtemp()

  ! Written by NPK 10/4/08

#ifdef MCKPP_CAM3
  USE shr_kind_mod, only: r8=>shr_kind_r8
  USE mckpp_types, only: kpp_3d_fields,kpp_const_fields
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid, only: get_ncols_p
#else 
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
#endif
  USE mckpp_parameters, ONLY: npts, nzp1
  
  IMPLICIT NONE
  
#ifdef MCKPP_CAM3
  INTEGER :: ichnk,ncol,icol
#endif  
  INTEGER ipt,z
  
#ifdef MCKPP_CAM3
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%tinc_fcorr(1:ncol,NZP1)=kpp_3d_fields(ichnk)%bottom_temp(1:ncol)-kpp_3D_fields(ichnk)%X(1:ncol,NZP1,1)
     kpp_3d_fields(ichnk)%ocnTcorr(1:ncol,NZP1)  =kpp_3d_fields(ichnk)%tinc_fcorr(1:ncol,NZP1)*&
          kpp_3d_fields(ichnk)%rho(1:ncol,NZP1)*kpp_3d_fields(ichnk)%cp(1:ncol,NZP1)/&
          kpp_const_fields%dto
     kpp_3d_fields(ichnk)%X(1:ncol,NZP1,1) = kpp_3d_fields(ichnk)%bottom_temp(1:ncol)
  ENDDO
#else
  DO ipt=1,npts
     kpp_3d_fields%tinc_fcorr(ipt,NZP1)=kpp_3d_fields%bottom_temp(ipt)-kpp_3D_fields%X(ipt,NZP1,1)
     kpp_3d_fields%ocnTcorr(ipt,NZP1)=kpp_3d_fields%tinc_fcorr(ipt,NZP1)*&
          kpp_3d_fields%rho(ipt,NZP1)*kpp_3d_fields%cp(ipt,NZP1)/&
          kpp_const_fields%dto
     kpp_3d_fields%X(ipt,NZP1,1) = kpp_3d_fields%bottom_temp(ipt)
  ENDDO
#endif
  
END SUBROUTINE mckpp_physics_overrides_bottomtemp


SUBROUTINE mckpp_physics_overrides_sst0()
  
  ! Written by NPK 27/8/07

#ifdef MCKPP_CAM3
  USE shr_kind_mod,only: r8=>shr_kind_r8
  USE mckpp_types, only: kpp_3d_fields
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid, only: get_ncols_p
#else
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
#endif
  USE mckpp_parameters, ONLY: nx, ny

  IMPLICIT NONE

#ifdef MCKPP_CAM3
  INTEGER :: ichnk,ncol,icol
#endif

  INTEGER ix,iy,ipoint

#ifdef MCKPP_CAM3
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%SST0(1:ncol) = kpp_3d_fields(ichnk)%sst(1:ncol)
  ENDDO
#else
  DO iy=1,ny
     DO ix=1,nx
        ipoint=(iy-1)*nx+ix
        kpp_3d_fields%SST0(ipoint)=kpp_3d_fields%SST(ix+kpp_const_fields%ifirst-1,iy+kpp_const_fields%jfirst-1)
     ENDDO
  ENDDO
#endif
  
END SUBROUTINE mckpp_physics_overrides_sst0


SUBROUTINE mckpp_physics_overrides_check_profile(kpp_1d_fields,kpp_const_fields)

#ifdef MCKPP_CAM3
  USE mckpp_types, only: kpp_1d_type,kpp_const_type
#else 
  USE mckpp_data_fields, ONLY: kpp_1d_type,kpp_const_type
#endif
  USE mckpp_log_messages, ONLY: mckpp_print_warning, max_message_len
  USE mckpp_parameters, ONLY: nzp1

  IMPLICIT NONE

  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields

  INTEGER :: z,j
  REAL :: dz_total,dtdz_total,dz
  
  CHARACTER(LEN=36) :: routine = "MCKPP_PHSYICS_OVERRIDE_CHECK_PROFILE"
  CHARACTER(LEN=max_message_len) :: message

  ! If the integration has failed because of unrealistic values in T, S, U or V
  ! or very high RMS difference between the old and new profiles, then reset
  ! T and S to climatology (if available) and U and V to the initial profiles.
  ! NPK 17/5/13.
  IF (kpp_1d_fields%comp_flag .and. kpp_const_fields%ocnT_file .ne. 'none' .and. &
       kpp_const_fields%sal_file .ne. 'none') THEN
     CALL mckpp_print_warning(routine, "Resetting point to climatology.")
     kpp_1d_fields%X(:,1)=kpp_1d_fields%ocnT_clim(:)
     ! WRITE(message,*) 'T = ',kpp_1d_fields%ocnT_clim(:)
     ! CALL mckpp_print_warning(routine, message)
     kpp_1d_fields%X(:,2)=kpp_1d_fields%sal_clim(:)
     ! WRITE(message,*) 'S = ',kpp_1d_fields%sal_clim(:)
     ! CALL mckpp_print_warning(routine, message)
     kpp_1d_fields%U=kpp_1d_fields%U_init(:,:)
     ! WRITE(message,*) 'U = ',kpp_1d_fields%U_init(:,1)
     ! CALL mckpp_print_warning(routine, message)
     ! WRITE(message,*) 'V = ',kpp_1d_fields%U_init(:,2)
     ! CALL mckpp_print_warning(routine, message)
     kpp_1d_fields%reset_flag=999
  ELSE IF (kpp_1d_fields%comp_flag) THEN
     WRITE(message,*) 'Cannot reset point to T,S climatology as either ocean temperature or salinity data '//&
          'not provided.  Will reset currents to initial conditions and keep going.'
     CALL mckpp_print_warning(routine, message)
     kpp_1d_fields%U=kpp_1d_fields%U_init(:,:)
     kpp_1d_fields%reset_flag=999
  ENDIF
  
  ! Check whether the temperature at any (x,z) point is less than the
  ! threshold for sea ice (-1.8C).  If it is, reset it to -1.8C and
  ! set a flag.  The flag can be requested as a diagnostic (singout 9).
  ! Note that the value of the flag is equal to the *fraction* of levels
  ! at that point that were < -1.8C.
  IF (kpp_1d_fields%L_OCEAN .and. kpp_const_fields%L_NO_FREEZE) THEN
     DO z=1,NZP1
        IF (kpp_1d_fields%X(z,1) .lt. -1.8) THEN
           kpp_1d_fields%tinc_fcorr(z)=kpp_1d_fields%tinc_fcorr(z)+&
                (-1.8-kpp_1d_fields%X(z,1))
           kpp_1d_fields%X(z,1)=-1.8
           kpp_1d_fields%freeze_flag=kpp_1d_fields%freeze_flag+1.0/REAL(NZP1)
        ENDIF
     ENDDO
  ENDIF
        
  ! Check whether the temperature difference between the surface
  ! and a user-specified level (presumably deep) is less than a user-specified 
  ! threshold (presumably small).  If so, reset the temperature and salinity
  ! profiles to climatological values.  Added to prevent spurious very
  ! deep mixing that creates unrealistic isothermal (and isohaline) layers.
  ! NPK 15/5/2013 for R4.
  IF (kpp_1d_fields%L_OCEAN .and. kpp_const_fields%L_NO_ISOTHERM) THEN
     dtdz_total=0.
     dz_total=0.
     DO j=2,kpp_const_fields%iso_bot
        dz=kpp_const_fields%zm(j)-kpp_const_fields%zm(j-1)
        dtdz_total=dtdz_total+ABS((kpp_1d_fields%X(j,1)-&
             kpp_1d_fields%X(j-1,1)))*dz
        dz_total=dz_total+dz
     ENDDO
     dtdz_total=dtdz_total/dz_total
     
     ! If resetting to climatology because of isothermal layer (rather than because of 
     ! computational instability trap in ocn.f), then set reset_flag to a negative
     ! value (-1*number of interations in of semi-implicit integration in ocn.f).
     IF (ABS(dtdz_total).lt.kpp_const_fields%iso_thresh) THEN
        kpp_1d_fields%X(:,1)=kpp_1d_fields%ocnT_clim(:)
        kpp_1d_fields%X(:,2)=kpp_1d_fields%sal_clim(:)
        kpp_1d_fields%reset_flag=(-1.)*kpp_1d_fields%reset_flag
     ENDIF
  ELSE
     kpp_1d_fields%reset_flag=0
  ENDIF
  
END SUBROUTINE mckpp_physics_overrides_check_profile

