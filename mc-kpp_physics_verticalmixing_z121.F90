SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_Z121 (kmp1,vlo,vhi,V,w)
  IMPLICIT NONE

  ! Necessary for IMPLICIT NONE (NPK 11/2/13)
  INTEGER kmp1
  REAL vlo,vhi,tmp,wait
  INTEGER k,km
  
  ! Apply 121 smoothing in k to 2-d array V(i,k=1,km)
  ! top (0) value is used as a dummy
  ! bottom (kmp1) value is set to input value from above.
 
  real V(0:kmp1)  ! 2-D array to be smoothed in kmp1 direction
  real w(0:kmp1)  ! 2-D array of internal weights to be computed
  
  km  = kmp1 - 1
  
  w(0)    =   0.0         
  w(kmp1) =   0.0             
  V(0)    =   0.0      
  V(kmp1) =   0.0
  
  do k=1,km  
     if((V(k).lt.vlo).or.(V(k).gt.vhi)) then
        w(k) = 0.0
        !     w(i,k) = 1.0
     else 
        w(k) = 1.0
     endif
     
  enddo
  
  do k=1,km
     tmp    = V(k)
     V(k) = w(k-1)*V(0)+2.*V(k)+w(k+1)*V(k+1)
     wait   = w(k-1) + 2.0 + w(k+1)
     V(k) = V(k) / wait             
     V(0) = tmp
  enddo
  
  RETURN
END SUBROUTINE MCKPP_PHYSICS_VERTICALMIXING_Z121
