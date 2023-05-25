MODULE mckpp_initialize_relaxtion_mod

CONTAINS

SUBROUTINE MCKPP_INITIALIZE_RELAXATION()
  
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
  USE mckpp_parameters, ONLY: nx, ny, ny_globe
  USE mckpp_physics_overrides, ONLY: mckpp_physics_overrides_sst0
  

! Re-write logic to allow for relaxing either SST or
! salinity - NPK 24/08/11

  IMPLICIT NONE

  INTEGER ix,iy,ipoint,my_ny
  CHARACTER(LEN=27) :: routine = "MCKPP_INITIALIZE_RELAXATION"
 
!  REAL sst_in(NX_GLOBE,NY_GLOBE,1)
!  COMMON /save_sstin/ sst_in

  my_ny=NY

  IF (kpp_const_fields%L_RELAX_SST) THEN 
     DO iy=1,my_ny
        IF (kpp_const_fields%relax_sst_in(iy) .EQ. 0.0) THEN
           DO ix=1,nx
              ipoint=(iy-1)*nx+ix
              kpp_3d_fields%relax_sst(ipoint)=0.0
           ENDDO
        ELSE
           DO ix=1,nx
              ipoint=(iy-1)*nx+ix
              kpp_3d_fields%relax_sst(ipoint)=1./(kpp_const_fields%relax_sst_in(iy)*kpp_const_fields%spd)
           ENDDO
        ENDIF
     ENDDO
  ENDIF
  IF (kpp_const_fields%L_RELAX_SAL) THEN
     DO iy=1,my_ny
        IF (kpp_const_fields%relax_sal_in(iy) .EQ. 0.0) THEN
           DO ix=1,nx
              ipoint=(iy-1)*nx+ix
              kpp_3d_fields%relax_sal(ipoint)=0.0
           ENDDO
        ELSE
           DO ix=1,nx
              ipoint=(iy-1)*nx+ix
              kpp_3d_fields%relax_sal(ipoint)=1./(kpp_const_fields%relax_sal_in(iy)*kpp_const_fields%spd)
           ENDDO
        ENDIF
     ENDDO
  ENDIF
  IF (kpp_const_fields%L_RELAX_OCNT) THEN
     DO iy=1,my_ny
        IF (kpp_const_fields%relax_ocnt_in(iy) .EQ. 0.0) THEN
           DO ix=1,nx
              ipoint=(iy-1)*nx+ix
              kpp_3d_fields%relax_ocnT(ipoint)=0.0
           ENDDO
        ELSE
           DO ix=1,nx
              ipoint=(iy-1)*nx+ix
              kpp_3d_fields%relax_ocnT(ipoint)=1./(kpp_const_fields%relax_ocnT_in(iy)*kpp_const_fields%spd)
           ENDDO
        ENDIF
     ENDDO
  ENDIF
  CALL MCKPP_PHYSICS_OVERRIDES_SST0()

! Do we need these initialization statements?
  DO iy=1,ny
     DO ix=1,nx
        ipoint=(iy-1)*nx+ix
        kpp_3d_fields%fcorr(ipoint)=0.0
        kpp_3d_fields%scorr(ipoint,:)=0.0
     ENDDO
  ENDDO
  
  CALL mckpp_print(routine, "Calculated SST0, fcorr and scorr")

END SUBROUTINE mckpp_initialize_relaxation

END MODULE mckpp_initialize_relaxtion_mod
