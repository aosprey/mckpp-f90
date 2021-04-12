! Control routines to be called from main KPP code. 
MODULE mckpp_xios_control

  USE mckpp_data_fields, ONLY: kpp_3d_fields,kpp_const_fields
  USE mckpp_parameters, ONLY: nuout
  USE mckpp_timer, ONLY: mckpp_start_timer, mckpp_stop_timer
  USE mckpp_xios_io, ONLY: xios_comm, ctx_hdl_diags, &
      mckpp_xios_diagnostic_definition, mckpp_xios_diagnostic_output, &
      mckpp_xios_restart_definition, mckpp_xios_write_restart

  USE mpi 
  USE xios 

  IMPLICIT NONE 

CONTAINS 

  ! Startup MPI and XIOS, then initialize diagnostic output 
  SUBROUTINE mckpp_initialize_output()

    INTEGER :: ierr

    CALL mpi_init(ierr)
    CALL xios_initialize("client", return_comm=xios_comm)    
    CALL mckpp_xios_diagnostic_definition()

  END SUBROUTINE mckpp_initialize_output


  ! Shutdown MPI and XIOS
  SUBROUTINE mckpp_finalize_output() 

    INTEGER :: ierr

    CALL xios_context_finalize()
    CALL xios_finalize()
    CALL mpi_finalize(ierr)

  END SUBROUTINE mckpp_finalize_output


  ! Diagnostic output
  SUBROUTINE mckpp_output_control() 

    ! Send diags to XIOS at every ts 
    CALL mckpp_xios_diagnostic_output() 

  END SUBROUTINE mckpp_output_control


  ! Restart output
  SUBROUTINE mckpp_restart_control() 

    REAL :: restart_time

    ! Check if restart timestep 
    ! - always write restart at end of run 
    CALL mckpp_start_timer("Write restart output") 
    IF (MOD(kpp_const_fields%ntime,kpp_const_fields%ndt_per_restart) .EQ. 0 .OR. &
        kpp_const_fields%ntime .EQ. kpp_const_fields%nend*kpp_const_fields%ndtocn) THEN

      ! Set correct time for validity of restart file
      ! (end of this timestep = start of next timestep)
      restart_time=kpp_const_fields%time+kpp_const_fields%dto/kpp_const_fields%spd

      WRITE(nuout,*) "Writing restart at time ", restart_time    
      CALL mckpp_xios_write_restart(restart_time) 
    END IF
    CALL mckpp_stop_timer("Write restart output") 

  END SUBROUTINE mckpp_restart_control

END MODULE mckpp_xios_control
