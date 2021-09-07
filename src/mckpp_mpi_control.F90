MODULE mckpp_mpi_control 

  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len 
  USE mpi
  USE xios

  IMPLICIT NONE 

  INTEGER :: comm, rank, nproc, root_proc
  INTEGER :: npts_local, offset_global
  INTEGER, DIMENSION(:), ALLOCATABLE :: npts_local_all, offset_global_all

CONTAINS

  ! Initialize MPI and XIOS. 
  ! Get communicator from XIOS 
  SUBROUTINE mckpp_initialize_mpi() 

    INTEGER :: ierror 
    CHARACTER(LEN=20) :: routine = "MCKPP_INITIALIZE_MPI"
    CHARACTER(LEN=max_message_len) :: message

    CALL mpi_init(ierror) 
    CALL xios_initialize("client", return_comm=comm)

    CALL mpi_comm_size(comm, nproc, ierror)
    CALL mpi_comm_rank(comm, rank, ierror)
    WRITE(message,*) "Rank ", rank, " of ", nproc
    CALL mckpp_print(routine, message) 

    IF (rank .EQ. 0) THEN
      root_proc = .TRUE.
    ELSE
      root_proc = .FALSE.
    END IF

  END SUBROUTINE mckpp_initialize_mpi


  SUBROUTINE mckpp_finalize_mpi

    INTEGER :: ierr

    CALL mpi_finalize(ierr)

  END SUBROUTINE mckpp_finalize_mpi
  

  ! 1d decomposition 
  ! Extra columns assigned from highest rank.
  SUBROUTINE mckpp_decompose_domain()

    USE mckpp_parameters, ONLY : npts, nuout

    INTEGER :: n, min, rem, index

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
        WRITE(nuout,*) n, npts_local_all(n), offset_global_all(n)
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
