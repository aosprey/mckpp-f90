MODULE mckpp_mpi_control 

  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_log_messages, ONLY: mckpp_initialize_logs, mckpp_finalize_logs, &
        mckpp_print, mckpp_print_error, max_message_len 
  USE mckpp_parameters, ONLY : npts, nz
  USE mpi
  USE xios

  IMPLICIT NONE

  INTEGER :: comm, rank, nproc, npts_local, offset_global, start_global, & 
             end_global, subdomain_type
  INTEGER :: root_proc = 0
  LOGICAL :: l_root_proc
  INTEGER, DIMENSION(:), ALLOCATABLE :: npts_local_all, offset_global_all

  INTERFACE mckpp_scatter_field
    MODULE PROCEDURE mckpp_scatter_field_1d, mckpp_scatter_field_2d
  END INTERFACE mckpp_scatter_field 

CONTAINS

  ! Initialize MPI, XIOS and log files. 
  ! Get communicator from XIOS 
  SUBROUTINE mckpp_initialize_mpi() 

    INTEGER :: ierr
    CHARACTER(LEN=20) :: routine = "MCKPP_INITIALIZE_MPI"
    CHARACTER(LEN=max_message_len) :: message

    CALL mpi_init(ierr) 
    CALL xios_initialize("client", return_comm=comm)

    CALL mpi_comm_size(comm, nproc, ierr)
    CALL mpi_comm_rank(comm, rank, ierr)

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
  ! nproc needs to divide npts exactly.
  ! Call after mckpp_initialize_namelist (so npts is set)
  SUBROUTINE mckpp_decompose_domain()

    INTEGER, DIMENSION(2) :: global_sizes, local_sizes, starts
    INTEGER :: tmp_type, dbl_size, ierr
    INTEGER(kind=MPI_ADDRESS_KIND) :: extent, start 
    CHARACTER(LEN=22) :: routine = "MCKPP_DECOMPOSE_DOMAIN"
    CHARACTER(LEN=max_message_len) :: message

    ! Check even distribution 
    IF ( MOD(npts, nproc) .NE. 0 ) THEN
      WRITE(message, *) "nproc, npts = ", nproc, npts
      CALL mckpp_print_error(routine, message)
      CALL mckpp_abort(routine, & 
        "Number of MPI processes must divide number of grid pts.")
    END IF

    ! Decompose
    npts_local = npts / nproc
   
    ! Work out global indices
    offset_global = (rank-1)*npts_local
    start_global = offset_global + 1 
    end_global = offset_global + npts_local 

    IF (l_root) THEN
      WRITE(message,*) "nproc, npts_local = ", nproc, npts_local 
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "offset_global, start_global, end_global = ", & 
                        offset_global, start_global, end_global
      CALL mckpp_print(routine, message) 
    END IF

    ! Setup derived type for sub-domains that we can use for scatters
    global_sizes = (/ npts, nz /)
    local_sizes = (/ npts_local, nz /)
    starts = (/ 0,0 /)
    CALL MPI_type_create_subarray( 2, global_sizes, local_sizes, starts, &
                                   MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, & 
                                   tmp_type, ierr )
    
    CALL MPI_type_size( MPI_DOUBLE_PRECISION, dbl_size, ierr )
    start = 0
    extent = dbl_size*npts_local

    CALL MPI_type_create_resized( tmp_type, start, extent, subdomain_type, ierr )
    CALL MPI_type_commit(subdomain_type, ierr)

  END SUBROUTINE mckpp_decompose_domain


  ! Scatter 1d array (npts) 
  SUBROUTINE mckpp_scatter_field_1d(global, local, root)

    REAL, DIMENSION(:), INTENT(IN) :: global
    REAL, DIMENSION(:), INTENT(OUT) :: local
    INTEGER, INTENT(IN) :: root

    INTEGER :: ierr

    CALL MPI_scatter( global, npts_local, MPI_DOUBLE_PRECISION, &
                      local, npts_local, MPI_DOUBLE_PRECISION, & 
                      root, comm, ierr )  

  END SUBROUTINE mckpp_scatter_field_1d

  
  ! Scatter 2d array (npts, nz) 
  SUBROUTINE mckpp_scatter_field_2d(global, local, root)

    REAL, DIMENSION(:,:), INTENT(IN) :: global
    REAL, DIMENSION(:,:), INTENT(OUT) :: local
    INTEGER, INTENT(IN) :: root

    INTEGER :: ierr

    CALL MPI_scatter( global, 1, subdomain_type, &
                      local, npts_local*nz, MPI_DOUBLE_PRECISION, & 
                      root, comm, ierr)  

  END SUBROUTINE mckpp_scatter_field_2d
  

  SUBROUTINE mckpp_broadcast_field(field, count, root)

    REAL, DIMENSION(:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: count, root

    INTEGER :: ierr

    CALL MPI_bcast( field, count, MPI_DOUBLE_PRECISION, root, comm, ierr ) 

  END SUBROUTINE mckpp_broadcast_field


  ! 1d decomposition
  ! - Where we suppose an unequal divsion. 
  ! - Extra columns assigned from highest rank.
  ! Not used currently 
  ! Call after mckpp_initialize_namelist (so npts is set)
  SUBROUTINE mckpp_decompose_domain_uneven()

    INTEGER :: n, min, rem, index
    CHARACTER(LEN=22) :: routine = "MCKPP_DECOMPOSE_DOMAIN_UNEVEN"
    CHARACTER(LEN=max_message_len) :: message

    IF (l_root) CALL mckpp_print(routine, "proc, offset_global, npts_local")
    
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

      IF (l_root_proc) THEN
        WRITE(message,*) n, offset_global_all(n), npts_local_all(n)
        CALL mckpp_print(routine, message) 
      END IF
    END DO

    npts_local = npts_local_all(rank)
    offset_global = offset_global_all(rank)

    start_global = offset_global + 1
    end_global = offset_global + npts_local
    
    WRITE(message,*) "offset_global, npts_local, start_global, end_global = ", &
                      offset_global, npts_local, start_global, end_global
    CALL mckpp_print(routine, message)

  END SUBROUTINE mckpp_decompose_domain_uneven

   
  ! Scatter 1d array dimensioned on npts
  ! - When we have an uneven distribution
  SUBROUTINE mckpp_scatter_field_uneven(global, local, root)

    REAL, DIMENSION(:), INTENT(IN) :: global
    REAL, DIMENSION(:), INTENT(OUT) :: local
    INTEGER, INTENT(IN) :: root

    INTEGER :: ierr

    CALL MPI_scatterv(global, npts_local_all, offset_global_all, & 
                      MPI_DOUBLE_PRECISION, &  
                      local, npts_local, MPI_DOUBLE_PRECISION, & 
                      root, comm, ierr)

  END SUBROUTINE mckpp_scatter_field_uneven


END MODULE mckpp_mpi_control
