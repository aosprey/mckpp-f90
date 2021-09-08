MODULE mckpp_mpi_control 

  USE mckpp_log_messages, ONLY: mckpp_initialize_logs, mckpp_finalize_logs, &
       mckpp_print, max_message_len 
  USE mckpp_parameters, ONLY : npts
  USE mpi
  USE xios

  IMPLICIT NONE 

  INTEGER :: comm, rank, nproc, npts_local, offset_global
  INTEGER :: root_proc = 0
  LOGICAL :: l_root_proc
  INTEGER, DIMENSION(:), ALLOCATABLE :: npts_local_all, offset_global_all

CONTAINS

  ! Initialize MPI, XIOS and log files. 
  ! Get communicator from XIOS 
  SUBROUTINE mckpp_initialize_mpi() 

    INTEGER :: ierror 
    CHARACTER(LEN=20) :: routine = "MCKPP_INITIALIZE_MPI"
    CHARACTER(LEN=max_message_len) :: message

    CALL mpi_init(ierror) 
    CALL xios_initialize("client", return_comm=comm)

    CALL mpi_comm_size(comm, nproc, ierror)
    CALL mpi_comm_rank(comm, rank, ierror)

    CALL mckpp_initialize_logs(rank)
    WRITE(message,*) "Rank ", rank, " of ", nproc
    CALL mckpp_print(routine, message)
    
    IF (rank .EQ. root_proc) THEN
      l_root_proc = .TRUE.
    ELSE
      l_root_proc = .FALSE.
    END IF

  END SUBROUTINE mckpp_initialize_mpi


  ! Shutdown XIOS, log files and MPI
  SUBROUTINE mckpp_finalize_mpi

    INTEGER :: ierr
    
    CALL xios_context_finalize()
    CALL xios_finalize()
    CALL mckpp_finalize_logs()
    CALL mpi_finalize(ierr)

  END SUBROUTINE mckpp_finalize_mpi
  

  ! 1d decomposition 
  ! Extra columns assigned from highest rank.
  ! Call after mckpp_initialize_namelist (so npts is set)
  SUBROUTINE mckpp_decompose_domain()

    INTEGER :: n, min, rem, index
    CHARACTER(LEN=22) :: routine = "MCKPP_DECOMPOSE_DOMAIN"
    CHARACTER(LEN=max_message_len) :: message

    ALLOCATE( npts_local_all(0:nproc-1) )
    ALLOCATE( offset_global_all(0:nproc-1) )

    min = npts / nproc
    rem = npts - (min*nproc)

    index = 0
    DO n = 0, nproc-1
      IF ( (nproc - n) .LE. rem) THEN
        npts_local_all(n) = min + rem
      ELSE
        npts_local_all(n) = min
      END IF

      offset_global_all(n) = index
      index = index + npts_local_all(n)

      IF (rank .EQ. 0) THEN
        WRITE(message,*) n, npts_local_all(n), offset_global_all(n)
        CALL mckpp_print(routine, message) 
      END IF
    END DO

    npts_local = npts_local_all(rank)
    offset_global = offset_global_all(rank)

  END SUBROUTINE mckpp_decompose_domain


  SUBROUTINE mckpp_scatter_field(global, local, root)

    REAL, DIMENSION(:), INTENT(IN) :: global
    REAL, DIMENSION(:), INTENT(OUT) :: local
    INTEGER, INTENT(IN) :: root

    INTEGER :: ierr

    CALL MPI_scatterv(global, npts_local_all, offset_global_all, MPI_DOUBLE_PRECISION, &
        local, npts_local, MPI_DOUBLE_PRECISION, root, comm, ierr)  

  END SUBROUTINE mckpp_scatter_field


  SUBROUTINE mckpp_broadcast_field(field, count, root)

    REAL, DIMENSION(:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: count, root

    INTEGER :: ierr

    CALL MPI_bcast(field, count, MPI_DOUBLE_PRECISION, root, comm, ierr) 

  END SUBROUTINE mckpp_broadcast_field


END MODULE mckpp_mpi_control
