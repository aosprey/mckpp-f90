MODULE mckpp_physics_verticalmixing_rimix_mod

  USE mckpp_data_fields, ONLY: kpp_1d_type, kpp_const_type
  USE mckpp_physics_verticalmixing_z121_mod, ONLY: mckpp_physics_verticalmixing_z121
 
  IMPLICIT NONE     

CONTAINS

! compute interior viscosity diffusivity coefficients due to
! shear instability (dependent on a local richardson number)
! and due to background internal wave activity.
SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_RIMIX(km,kmp1,kpp_1d_fields,kpp_const_fields)

  TYPE(kpp_1d_type) :: kpp_1d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  
  integer km,kmp1           ! number of vertical levels

  ! local variables
  real Rigg
  real fri,fcon             ! function of Rig
  real ratio
  real epsln,Riinfty,Ricon,difm0,difs0,difmiw,difsiw,difmcon,difscon,c1,c0
  INTEGER ki,mRi,j
      
  epsln   = 1.e-16   ! a small number          
  Riinfty = 0.8      ! LMD default was = 0.7
  Ricon   = -0.2     ! note: exp was repl by multiplication
  difm0   = 0.005    ! max visc due to shear instability
  difs0   = 0.005    ! max diff ..  .. ..    ..
  difmiw  = 0.0001   ! background/internal waves visc(m^2/s)
  difsiw  = 0.00001  ! ..         ..       ..    diff(m^2/s)
  difmcon = 0.0000   ! max visc for convection  (m^2/s)
  difscon = 0.0000   ! max diff for convection  (m^2/s)
  c1      = 1.0
  c0      = 0.0
  mRi     = 1        ! number of vertical smoothing passes
     
  !     compute interior gradient Ri at all interfaces, except surface     
  !-----------------------------------------------------------------------
  !     compute interior gradient Ri at all interfaces ki=1,km, (not surface)
  !     use visc(imt,ki=1,km) as temporary storage to be smoothed
  !     use dift(imt,ki=1,km) as temporary storage of unsmoothed Ri
  !     use difs(imt,ki=1,km) as dummy in smoothing call
      
  do ki = 1, km
     kpp_1d_fields%Rig(ki)  = kpp_1d_fields%dbloc(ki) * (kpp_const_fields%zm(ki)-kpp_const_fields%zm(ki+1))/&
          (kpp_1d_fields%Shsq(ki) + epsln)
     kpp_1d_fields%dift(ki) = kpp_1d_fields%Rig(ki)
     kpp_1d_fields%difm(ki) = kpp_1d_fields%dift(ki)                        
  ENDDO
  
  !-----------------------------------------------------------------------
  !     vertically smooth Ri mRi times
  do j = 1,mRi
     call MCKPP_PHYSICS_VERTICALMIXING_Z121(kmp1,c0,Riinfty,kpp_1d_fields%difm,kpp_1d_fields%difs)
  enddo

  !-----------------------------------------------------------------------
  !                           after smoothing loop
  DO ki = 1, km
     ! evaluate f of unsmooth Ri (fri) for convection        store in fcon
     ! evaluate f of   smooth Ri (fri) for shear instability store in fri
 
     Rigg  = AMAX1( kpp_1d_fields%dift(ki) , Ricon )
     ratio = AMIN1( (Ricon-Rigg)/Ricon , c1 )
     fcon  = (c1 - ratio*ratio)
     fcon  = fcon * fcon * fcon
     
     Rigg  = AMAX1( kpp_1d_fields%difm(ki) , c0 )
     ratio = AMIN1( Rigg/Riinfty , c1 )
     fri   = (c1 - ratio*ratio)
     fri   = fri * fri * fri
         	     
     ! ************************   Overwrite with Gent's PP **********
     !           fcon = 0.0
     !           Rigg  = AMAX1( dift(i,ki) , c0 )
     !           fri   = c1 / (c1 + 10. * Rigg )
     !           difm0 = 0.1 * fri
     !           difs0 = 0.1 * fri * fri

     !  ************************   Overwrite with original PP
     !           fcon = 0.0
     !           Rigg  = AMAX1( dift(i,ki) , c0 )
     !           fri   = c1 / (c1 +  5. * Rigg )
     !           difm0 = 0.01 * fri
     !           difs0 = (difmiw + fri * difm0)

     ! ----------------------------------------------------------------------
     !            evaluate diffusivities and viscosity
     ! mixing due to internal waves, and shear and static instability
 
     kpp_1d_fields%difm(ki) = (difmiw + fcon * difmcon + fri * difm0)
     kpp_1d_fields%difs(ki) = (difsiw + fcon * difscon + fri * difs0)
     kpp_1d_fields%dift(ki) = kpp_1d_fields%difs(ki)
  END DO

  ! ------------------------------------------------------------------------
  !         set surface values to 0.0
  
  kpp_1d_fields%difm(0) = c0
  kpp_1d_fields%dift(0) = c0
  kpp_1d_fields%difs(0) = c0      
  
end SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_RIMIX

END MODULE mckpp_physics_verticalmixing_rimix_mod
