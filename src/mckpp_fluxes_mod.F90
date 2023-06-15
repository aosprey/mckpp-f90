MODULE mckpp_fluxes_mod

  USE mckpp_data_fields, ONLY: kpp_3d_fields, kpp_const_fields, kpp_1d_type, & 
        kpp_const_type
  USE mckpp_log_messages, ONLY: mckpp_print, max_message_len, nupe
  USE mckpp_mpi_control, ONLY: l_root, npts_local
  USE mckpp_parameters, ONLY: npts, nz, nsflxs
  USE mckpp_read_fluxes_mod, ONLY: mckpp_initialize_fluxes_file, & 
        mckpp_read_fluxes
  USE mckpp_time_control, ONLY: ntime
  USE mckpp_types_transfer, ONLY: mckpp_fields_3dto1d, mckpp_fields_1dto3d
  
  IMPLICIT NONE

CONTAINS 


  ! Set up parameters for calculating fluxes and initialize fluxes.
  ! intermediate values computed every ndtld
  SUBROUTINE mckpp_initialize_fluxes()

    ! extrapolation parameters for fluxes
    ! Initialize flux arrays
    kpp_3d_fields%wu = 0.0
    kpp_3d_fields%wx = 0.0
    kpp_3d_fields%wxnt = 0.0
    kpp_3d_fields%sflux = 0.0
    kpp_3d_fields%sflux(:,:,5,0) = 1e-20

    IF ( kpp_const_fields%l_fluxdata ) THEN 
      IF ( l_root ) CALL mckpp_initialize_fluxes_file()
    END IF  

  END SUBROUTINE mckpp_initialize_fluxes


  SUBROUTINE mckpp_fluxes()

    TYPE(kpp_1d_type) :: kpp_1d_fields
    REAL(8), DIMENSION(npts) :: taux, tauy, swf, lwf, lhf, shf, rain, snow
    INTEGER :: ipt

    IF ( .NOT. kpp_const_fields%l_fluxdata ) THEN
      taux = 0.01
      tauy = 0.0
      swf = 200.0
      lwf = 0.0
      lhf = -150.0
      shf = 0.0
      rain = 6e-5
      snow = 0.0   
    ELSE
      CALL mckpp_read_fluxes(taux, tauy, swf, lwf, lhf, shf, rain, snow)
    END IF

    WRITE(nupe, *) "taux = ", taux
    WRITE(nupe, *) "tauy = ", tauy
    WRITE(nupe, *) "swf = ", swf
    WRITE(nupe, *) "lwf = ", lwf
    WRITE(nupe, *) "lhf = ", lhf
    WRITE(nupe, *) "shf = ", shf 
    WRITE(nupe, *) "rain = ", rain 
    WRITE(nupe, *) "snow = ", snow

      
!!$OMP PARALLEL DEFAULT(shared) PRIVATE(ipt,kpp_1d_fields)
!!$OMP DO SCHEDULE(dynamic)
    DO ipt = 1, npts_local        
      IF ( kpp_3d_fields%l_ocean(ipt) ) THEN 

        IF ( (taux(ipt) .EQ. 0.0) .AND. (tauy(ipt) .EQ. 0.0) ) & 
          taux(ipt) = 1.e-10

        IF ( .NOT. kpp_const_fields%l_rest ) THEN
          kpp_3d_fields%sflux(ipt,1,5,0) = taux(ipt)
          kpp_3d_fields%sflux(ipt,2,5,0) = tauy(ipt)
          kpp_3d_fields%sflux(ipt,3,5,0) = swf(ipt)            
          kpp_3d_fields%sflux(ipt,4,5,0) = lwf(ipt) + lhf(ipt) + shf(ipt) - & 
                                           snow(ipt) * kpp_const_fields%flsn
          kpp_3d_fields%sflux(ipt,5,5,0) = 1e-10 ! Melting of sea-ice = 0.0
          kpp_3d_fields%sflux(ipt,6,5,0) = rain(ipt) + snow(ipt) + & 
                                           (lhf(ipt) / kpp_const_fields%el)
        ELSE
          kpp_3d_fields%sflux(ipt,1,5,0) = 1.e-10
          kpp_3d_fields%sflux(ipt,2,5,0) = 0.00
          kpp_3d_fields%sflux(ipt,3,5,0) = 300.00
          kpp_3d_fields%sflux(ipt,4,5,0) = -300.00
          kpp_3d_fields%sflux(ipt,5,5,0) = 0.00
          kpp_3d_fields%sflux(ipt,6,5,0) = 0.00
        END IF

        CALL mckpp_fields_3dto1d(kpp_3d_fields, ipt, kpp_1d_fields)
        CALL mckpp_fluxes_ntflux(kpp_1d_fields, kpp_const_fields)
        CALL mckpp_fields_1dto3d(kpp_1d_fields, ipt, kpp_3d_fields)

      END IF
    END DO
!!$OMP END DO
!!$OMP END PARALLEL

    WRITE(nupe,*) "kpp_3d_fields%wxnt(:,1,1) = ", kpp_3d_fields%wxnt(:,1,1)
  
  END SUBROUTINE mckpp_fluxes


  ! Non-turbulent fluxes
  SUBROUTINE mckpp_fluxes_ntflux(kpp_1d_fields, kpp_const_fields)
  
    INTEGER :: k
    TYPE(kpp_1d_type) :: kpp_1d_fields
    TYPE(kpp_const_type) :: kpp_const_fields
    CHARACTER(LEN=31) :: routine = "MCKPP_FLUXES_NTFLUX"
    CHARACTER(LEN=max_message_len) :: message
  
    ! WRITE(message,*) "At time = ", ntime
    ! CALL mckpp_print(routine, message)
    IF ( ntime .LE. 1 ) THEN
      DO k = 0, nz
        kpp_1d_fields%swdk_opt(k) = mckpp_fluxes_swdk( & 
          -kpp_const_fields%dm(k), kpp_1d_fields%jerlov )
      END DO
    END IF

    IF ( ntime .GE. 1 ) THEN 
      DO k = 0, nz
         kpp_1d_fields%wxnt(k,1) = &
           -kpp_1d_fields%sflux(3,5,0) * kpp_1d_fields%swdk_opt(k) / & 
           ( kpp_1d_fields%rho(0) * kpp_1d_fields%cp(0) )
      END DO 
    END IF

  END SUBROUTINE mckpp_fluxes_ntflux


  REAL FUNCTION mckpp_fluxes_swdk(z, j)
  
    REAL :: z
    INTEGER :: j ! jerlov
    INTEGER, PARAMETER :: max=5
    REAL, DIMENSION(max) :: rfac, a1, a2

!      types =  I       IA      IB      II      III
!          j =  1       2       3       4       5
    rfac = (/  0.58 ,  0.62 ,  0.67 ,  0.77 ,  0.78 /)
    a1 =   (/  0.35 ,  0.6  ,  1.0  ,  1.5  ,  1.4  /)
    a2 =   (/ 23.0  , 20.0  , 17.0  , 14.0  ,  7.9  /)
  
    mckpp_fluxes_swdk = rfac(j) * EXP( z / a1(j) ) + ( 1.0 - Rfac(j) ) * & 
                        EXP( z / a2(j) )

  END FUNCTION mckpp_fluxes_swdk

END MODULE mckpp_fluxes_mod
