MODULE mckpp_initialize_relaxtion_mod

CONTAINS

  SUBROUTINE mckpp_initialize_relaxation()

    USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields
    USE mckpp_log_messages, ONLY: mckpp_print, max_message_len
    USE mckpp_mpi_control, ONLY: inds_global, npts_local
    USE mckpp_physics_overrides, ONLY: mckpp_physics_overrides_sst0

    ! Re-write logic to allow for relaxing either SST or
    ! salinity - NPK 24/08/11

    IMPLICIT NONE

    INTEGER i, j, ipt
    CHARACTER(LEN=27) :: routine = "MCKPP_INITIALIZE_RELAXATION"

    CALL mckpp_print(routine, "") 

    IF ( kpp_const_fields%l_relax_sst ) THEN 

      DO ipt = 1, npts_local 
        j = inds_global(ipt,2) 

        IF ( kpp_const_fields%relax_sst_in(j) .EQ. 0.0 ) THEN
          kpp_3d_fields%relax_sst(ipt) = 0.0
        ELSE 
          kpp_3d_fields%relax_sst(ipt) = 1. / & 
              ( kpp_const_fields%relax_sst_in(j) * kpp_const_fields%spd )
        END IF

      END DO

    END IF

    IF ( kpp_const_fields%l_relax_sal ) THEN

      DO ipt = 1, npts_local 
        j = inds_global(ipt,2) 
     
        IF ( kpp_const_fields%relax_sal_in(j) .EQ. 0.0 ) THEN
          kpp_3d_fields%relax_sal(ipt) = 0.0
        ELSE
          kpp_3d_fields%relax_sal(ipt) = 1. / & 
              ( kpp_const_fields%relax_sal_in(j) * kpp_const_fields%spd )
        ENDIF

      ENDDO

    ENDIF
 
    IF ( kpp_const_fields%l_relax_ocnt ) THEN

      DO ipt = 1, npts_local 
        j = inds_global(ipt,2) 

        IF ( kpp_const_fields%relax_ocnt_in(j) .EQ. 0.0 ) THEN
          kpp_3d_fields%relax_ocnt(ipt) = 0.0
        ELSE
          kpp_3d_fields%relax_ocnt(ipt) = 1. / & 
              ( kpp_const_fields%relax_ocnT_in(j) * kpp_const_fields%spd )
         ENDIF

      ENDDO

    ENDIF

    CALL mckpp_physics_overrides_sst0()

    kpp_3d_fields%fcorr(:) = 0.0
    kpp_3d_fields%scorr(:,:) = 0.0

    CALL mckpp_print(routine, "Calculated SST0, fcorr and scorr")

  END SUBROUTINE mckpp_initialize_relaxation


END MODULE mckpp_initialize_relaxtion_mod
