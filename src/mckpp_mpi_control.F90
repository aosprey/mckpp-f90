MODULE mckpp_mpi_control 

  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_log_messages, ONLY: mckpp_initialize_logs, mckpp_finalize_logs, &
        mckpp_print, mckpp_print_error, max_message_len, l_debug
  USE mckpp_parameters, ONLY : nx, npts, nzp1

  USE mpi  
  USE xios

  IMPLICIT NONE

  INTEGER :: comm, rank, nproc
  INTEGER :: root = 0
  LOGICAL :: l_root

  INTEGER :: npts_local, &                     ! local data size 
             offset_global, &                  ! pos in global 1d array of npts
             start_global, end_global          ! inds in global 1d array of npts
  
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: inds_global ! i,j inds in global 2d 
                                                      ! array of nx, ny

  INTEGER :: ni, nj, ibegin, jbegin, data_ni, data_ibegin   ! XIOS domain defn
  
  ! Only reqd for 2d domain decomposition - not currently used 
  INTEGER, DIMENSION(:), ALLOCATABLE :: npts_local_all, offset_global_all

  INTEGER, PRIVATE :: subdomain_type


  INTERFACE mckpp_broadcast_field 
    MODULE PROCEDURE mckpp_broadcast_field_real, mckpp_broadcast_field_logical
  END INTERFACE mckpp_broadcast_field 

  INTERFACE mckpp_scatter_field
    MODULE PROCEDURE mckpp_scatter_field_real_1d, mckpp_scatter_field_real_2d, &
      mckpp_scatter_field_int_1d
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
    
    IF (rank .EQ. root) THEN
      l_root = .TRUE.
    ELSE
      l_root = .FALSE.
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

    INTEGER :: ind, ipt
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
   
    ! Work out index in global 1d (npts) array 
    offset_global = rank*npts_local
    start_global = offset_global + 1 
    end_global = offset_global + npts_local 

    ! Work out i, j indices in 2d (nx,ny) array
    ALLOCATE( inds_global(npts_local,2) ) 

    ipt = 1 
    DO ind = start_global, end_global 
      inds_global(ipt,1) = MOD( (ind-1), nx ) + 1
      inds_global(ipt,2) = ( (ind-1) / nx ) + 1 
      ipt = ipt+1
    END DO 

    IF (l_debug) THEN 
      WRITE(message,*) "nproc, npts_local = ", nproc, npts_local 
      CALL mckpp_print(routine, message) 
      WRITE(message,*) "offset_global, start_global, end_global = ", & 
                        offset_global, start_global, end_global
      CALL mckpp_print(routine, message) 
    END IF 
      
    ! Set up derived type required for scatter 
    CALL mckpp_setup_types() 
   
    ! Set up variables to define local domain to XIOS 
    CALL mckpp_define_xios_domain()

  END SUBROUTINE mckpp_decompose_domain


  ! Setup derived type for sub-domains that we can use for scatters
  ! This is only required for 2d (npts, nz) fields. 
  SUBROUTINE mckpp_setup_types() 

    INTEGER, DIMENSION(2) :: global_sizes, local_sizes, starts
    INTEGER :: ierr, tmp_type, dbl_size
    INTEGER(kind=MPI_ADDRESS_KIND) :: extent, start 

    global_sizes = (/ npts, nzp1 /)
    local_sizes = (/ npts_local, nzp1 /)
    starts = (/ 0,0 /)

    CALL MPI_type_create_subarray( 2, global_sizes, local_sizes, starts, &
                                   MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, & 
                                   tmp_type, ierr )
    
    CALL MPI_type_size( MPI_DOUBLE_PRECISION, dbl_size, ierr )
    start = 0
    extent = dbl_size*npts_local

    CALL MPI_type_create_resized( tmp_type, start, extent, & 
                                  subdomain_type, ierr )
    CALL MPI_type_commit( subdomain_type, ierr )

  END SUBROUTINE mckpp_setup_types


  ! Define local domain for XIOS 
  SUBROUTINE mckpp_define_xios_domain()

    CHARACTER(LEN=24) :: routine = "MCKPP_DEFINE_XIOS_DOMAIN"
    CHARACTER(LEN=max_message_len) :: message

    ! Work out local domain in the way XIOS expects. 
    ! This is complicated for a 1d decmop, as we need to define a 2d ni*nj
    ! domain which covers all local points. 

    nj = inds_global(npts_local, 2) - inds_global(1, 2) + 1 
    jbegin = inds_global(1,2) - 1

    ! Two cases: 
    ! - The subdomain is a part of a row
    IF (nj .EQ. 1) THEN 
      ni = npts_local
      ibegin = inds_global(1,1) - 1
      data_ibegin = 0

    ! - The subdomain spans multiple rows
    ELSE 
      ni = nx 
      ibegin = 0 
      data_ibegin = inds_global(1,1) - 1
    END IF

    IF (l_debug) THEN
      WRITE(message, *) "ni, nj, ibegin, jbegin, data_ni, data_ibegin = ", & 
                         ni, nj, ibegin, jbegin, data_ni, data_ibegin
    CALL mckpp_print(message, routine) 
    END IF 

  END SUBROUTINE mckpp_define_xios_domain


  ! Scatter 1d real array (npts) 
  SUBROUTINE mckpp_scatter_field_real_1d(global, local, root)

    REAL, DIMENSION(npts), INTENT(IN) :: global
    REAL, DIMENSION(npts_local), INTENT(OUT) :: local
    INTEGER, INTENT(IN) :: root

    INTEGER :: ierr

    CALL MPI_scatter( global, npts_local, MPI_DOUBLE_PRECISION, &
                      local, npts_local, MPI_DOUBLE_PRECISION, & 
                      root, comm, ierr )  

  END SUBROUTINE mckpp_scatter_field_real_1d


  ! Scatter 1d integer array (npts) 
  SUBROUTINE mckpp_scatter_field_int_1d(global, local, root)

    INTEGER, DIMENSION(npts), INTENT(IN) :: global
    INTEGER, DIMENSION(npts_local), INTENT(OUT) :: local
    INTEGER, INTENT(IN) :: root
    INTEGER :: ierr

    CALL MPI_scatter( global, npts_local, MPI_INTEGER, &
                      local, npts_local, MPI_INTEGER, & 
                      root, comm, ierr )  

  END SUBROUTINE mckpp_scatter_field_int_1d

  
  ! Scatter 2d array (npts, nzp1) 
  ! - needs to be these exact sizes 
  SUBROUTINE mckpp_scatter_field_real_2d(global, local, root)

    REAL, DIMENSION(npts,nzp1), INTENT(IN) :: global
    REAL, DIMENSION(npts_local,nzp1), INTENT(OUT) :: local
    INTEGER, INTENT(IN) :: root
    INTEGER :: ierr

    CALL MPI_scatter( global, 1, subdomain_type, &
                      local, npts_local*nzp1, MPI_DOUBLE_PRECISION, & 
                      root, comm, ierr)  

  END SUBROUTINE mckpp_scatter_field_real_2d
  

  SUBROUTINE mckpp_broadcast_field_real(field, count, root)

    REAL, DIMENSION(:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: count, root
    INTEGER :: ierr

    CALL MPI_bcast( field, count, MPI_DOUBLE_PRECISION, root, comm, ierr ) 

  END SUBROUTINE mckpp_broadcast_field_real 


  SUBROUTINE mckpp_broadcast_field_logical(field, count, root) 

    LOGICAL, DIMENSION(:), INTENT(INOUT) :: field
    INTEGER, INTENT(IN) :: count, root
    INTEGER :: ierr

    CALL MPI_bcast( field, count, MPI_LOGICAL, root, comm, ierr ) 

  END SUBROUTINE mckpp_broadcast_field_logical


  ! 1d decomposition
  ! - Where we suppose an unequal divsion. 
  ! - Extra columns assigned from highest rank.
  ! Not used currently 
  ! Call after mckpp_initialize_namelist (so npts is set)
  SUBROUTINE mckpp_decompose_domain_uneven()

    INTEGER :: n, min, rem, index
    CHARACTER(LEN=22) :: routine = "MCKPP_DECOMPOSE_DOMAIN_UNEVEN"
    CHARACTER(LEN=max_message_len) :: message

    IF (l_debug .AND. l_root) THEN
      CALL mckpp_print(routine, "proc, offset_global, npts_local")
    END IF 
      
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

      IF (l_debug .AND. l_root) THEN
        WRITE(message,*) n, offset_global_all(n), npts_local_all(n)
        CALL mckpp_print(routine, message) 
      END IF
    END DO

    npts_local = npts_local_all(rank)
    offset_global = offset_global_all(rank)

    start_global = offset_global + 1
    end_global = offset_global + npts_local

    IF (l_debug) THEN 
      WRITE(message,*) "offset_global, npts_local, start_global, end_global = ", &
                        offset_global, npts_local, start_global, end_global
      CALL mckpp_print(routine, message)
    END IF

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
