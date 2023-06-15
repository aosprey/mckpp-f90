MODULE mckpp_initialize_ocean

  USE mckpp_abort_mod, ONLY: mckpp_abort
  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields, kpp_1d_type
  USE mckpp_log_messages, ONLY: mckpp_print
  USE mckpp_mpi_control, ONLY: npts_local
  USE mckpp_parameters
  USE mckpp_types_transfer, ONLY: mckpp_fields_3dto1d, mckpp_fields_1dto3d
  USE mckpp_physics_verticalmixing_mod, ONLY: mckpp_physics_verticalmixing
  
  IMPLICIT NONE

CONTAINS

! Initialize ocean model:
! Set coefficients for tridiagonal matrix solver.
! Compute hmix and diffusivity profiles for initial profile.
! Prepare for first time step.
SUBROUTINE MCKPP_INITIALIZE_OCEAN_MODEL()
    
  ! Local
  TYPE(kpp_1d_type) :: kpp_1d_fields
  real dzb(NZ)              ! diff. between grid-levels below z(j)
  integer k,kmix0,n,l,ipt
  real hmix0,deltaz
  CHARACTER(LEN=28) :: routine = "MCKPP_INITIALIZE_OCEAN_MODEL"

  CALL mckpp_print(routine, "")

  ! Compute factors for coefficients of tridiagonal matrix elements.
  ! tri(0     ,1,.........) : dt/h(1) factor for rhs flux
  ! tri(k=1:NZ,0,.........) : dt/h(k)/ {dzb(k-1)=z(k-1)-z(k)=dzabove}
  ! tri(k=1:NZ,1,.........) : dt/h(k)/ {dzb(k  )=z(k)-z(k+1)=dzbelow}

  DO k=1,NZ
     dzb(k) = kpp_const_fields%zm(k) - kpp_const_fields%zm(k+1)
  ENDDO
  
  kpp_const_fields%tri(0,1,1) = kpp_const_fields%dto/kpp_const_fields%hm(1)
  kpp_const_fields%tri(1,1,1) = kpp_const_fields%dto/kpp_const_fields%hm(1)/dzb(1)
  DO k=2,NZ
     kpp_const_fields%tri(k,1,1) = kpp_const_fields%dto/kpp_const_fields%hm(k)/dzb(k)
     kpp_const_fields%tri(k,0,1) = kpp_const_fields%dto/kpp_const_fields%hm(k)/dzb(k-1)
  ENDDO
  
  IF ( .NOT. kpp_const_fields%L_RESTART) THEN    
     ! Determine hmix for initial profile:

!$OMP PARALLEL DEFAULT(none) &
!$OMP SHARED(kpp_3d_fields, kpp_const_fields) &
!$OMP SHARED(nz, nzp1, nx, ny, npts_local, nvel, nsclr, nvp1, nsp1, itermax) &
!$OMP SHARED(hmixtolfrac, nztmax, nzp1tmax, nsflxs, njdt, maxmodeadv) &
!$OMP PRIVATE(ipt, k, l, deltaz, kpp_1d_fields, hmix0, kmix0)
!$OMP DO SCHEDULE(dynamic)
     DO ipt=1,npts_local
        IF (kpp_3d_fields%run_physics(ipt)) THEN

          CALL mckpp_fields_3dto1d(kpp_3d_fields,ipt,kpp_1d_fields)

           kpp_1d_fields%L_INITFLAG=.TRUE.
           CALL MCKPP_PHYSICS_VERTICALMIXING(kpp_1d_fields,kpp_const_fields,hmix0,kmix0)
           kpp_1d_fields%L_INITFLAG=.FALSE.
           kpp_1d_fields%hmix = hmix0
           kpp_1d_fields%kmix = kmix0
           kpp_1d_fields%Tref = kpp_1d_fields%X(1,1)
           ! Evaluate initial fluxes (to write to output data file)
           DO k=1,NZ
              deltaz = 0.5*(kpp_const_fields%hm(k)+kpp_const_fields%hm(k+1))
              DO n=1,NSCLR
                 kpp_1d_fields%wX(k,n)=-kpp_1d_fields%difs(k)*&
                      ((kpp_1d_fields%X(k,n)-kpp_1d_fields%X(k+1,n))/deltaz-&
                      kpp_1d_fields%ghat(k)*kpp_1d_fields%wX(0,n))
              ENDDO
              IF(kpp_const_fields%LDD) kpp_1d_fields%wX(k,1)=-kpp_1d_fields%dift(k)*&
                   ((kpp_1d_fields%X(k,1)-kpp_1d_fields%X(k+1,1))/deltaz-kpp_1d_fields%ghat(k)*kpp_1d_fields%wX(0,1))
              kpp_1d_fields%wX(k,nsp1)= kpp_const_fields%grav * (kpp_1d_fields%talpha(k)*kpp_1d_fields%wX(k,1) - &
                   kpp_1d_fields%sbeta(k) * kpp_1d_fields%wX(k,2))
              DO  n=1,NVEL
                 kpp_1d_fields%wU(k,n)= -kpp_1d_fields%difm(k)*&
                      (kpp_1d_fields%U(k,n)-kpp_1d_fields%U(k+1,n))/deltaz
              ENDDO
           ENDDO
           
           ! Prepare for first time step
           
           ! indices for extrapolation
           kpp_1d_fields%old = 0
           kpp_1d_fields%new = 1               
           ! initialize array for extrapolating hmixd,Us,Xs
           kpp_1d_fields%hmixd(0) = kpp_1d_fields%hmix
           kpp_1d_fields%hmixd(1) = kpp_1d_fields%hmix
           DO k=1,NZP1
              DO l=1,NVEL
                 kpp_1d_fields%Us(k,l,0)=kpp_1d_fields%U(k,l)
                 kpp_1d_fields%Us(k,l,1)=kpp_1d_fields%U(k,l)
              ENDDO
              DO l=1,NSCLR
                 kpp_1d_fields%Xs(k,l,0)=kpp_1d_fields%X(k,l)
                 kpp_1d_fields%Xs(k,l,1)=kpp_1d_fields%X(k,l)
              ENDDO
            ENDDO
            
          CALL mckpp_fields_1dto3d(kpp_1d_fields,ipt,kpp_3d_fields)

        ENDIF
     ENDDO
!$OMP END DO
!$OMP END PARALLEL

  ENDIF

END SUBROUTINE mckpp_initialize_ocean_model

END MODULE mckpp_initialize_ocean
