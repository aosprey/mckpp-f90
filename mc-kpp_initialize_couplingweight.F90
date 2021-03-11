#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE MCKPP_INITIALIZE_COUPLINGWEIGHT
  USE shr_kind_mod,only: r8=>shr_kind_r8, r4=>shr_kind_r4
  USE pmgrid, only: masterproc
  USE ppgrid, only: begchunk,endchunk,pcols
  USE phys_grid, only: scatter_field_to_chunk,get_ncols_p
  USE mckpp_parameters
  USE mckpp_types,only: kpp_global_fields,kpp_3d_fields,kpp_const_fields
#else
SUBROUTINE mckpp_initialize_couplingweight(kpp_3d_fields,kpp_const_fields)
  USE mckpp_data_types
#endif  

  IMPLICIT NONE

#ifdef MCKPP_CAM3
  REAL(r8) :: cplwght(PLON,PLAT),cplwght_chunk(PCOLS,begchunk:endchunk)
  INTEGER :: icol,ncol,ichnk,lat_varid,lon_varid
  REAL(r4) :: latitude_in(PLAT),longitude_in(PLON)
#else
  TYPE(kpp_3d_type) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
#endif

  INTEGER start(2),count(2)
  INTEGER ix,jy,ipoint,cplwght_varid,status,ncid_cplwght
  INTEGER ipoint_globe

#include <netcdf.inc>
!#include <couple.com>
  
  REAL*4 ixx, jyy, cplwght_in(NX_GLOBE,NY_GLOBE)
  
!  If L_CPLWGHT has been set, then we will use the
!  NetCDF file to set values of cplwght over the
!  entire globe.
!  Otherwise, we will set the values ourselves, based
!  on the coupling region.
!  NPK 10/9/07 - R1
  
#ifdef MCKPP_CAM3
  IF (masterproc) THEN
#endif
  status=NF_OPEN(kpp_const_fields%cplwght_file,0,ncid_cplwght)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  start(1) = 1
  start(2) = 1
  count(1) = NX_GLOBE
  count(2) = NY_GLOBE
  WRITE(6,*) 'MCKPP_INITIALIZE_COUPLINGWEIGHT: Reading coupling weight (alpha)'
  status=NF_INQ_VARID(ncid_cplwght,'alpha',cplwght_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VARA_REAL(ncid_cplwght,cplwght_varid,start,count,cplwght_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
#ifdef MCKPP_CAM3

  ! Use coupling weight to get global latitude and longitude grid
  status=NF_INQ_VARID(ncid_cplwght,'latitude',lat_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VARA_REAL(ncid_cplwght,lat_varid,1,PLAT,latitude_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)  
  
  status=NF_INQ_VARID(ncid_cplwght,'longitude',lon_varid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  status=NF_GET_VARA_REAL(ncid_cplwght,lon_varid,1,PLON,longitude_in)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  kpp_global_fields%latitude(:)=latitude_in
  kpp_global_fields%longitude(:)=longitude_in 

  status=NF_CLOSE(ncid_cplwght)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ENDIF ! End of masterproc section
  cplwght(:,:)=cplwght_in(:,:)
  CALL scatter_field_to_chunk(1,1,1,PLON,cplwght,cplwght_chunk(1,begchunk))
  DO ichnk=begchunk,endchunk
     ncol=get_ncols_p(ichnk)
     kpp_3d_fields(ichnk)%cplwght(1:ncol)=cplwght_chunk(1:ncol,ichnk)        
  ENDDO
#else
  DO ix=1,NX_GLOBE
     DO jy=1,NY_GLOBE
        ipoint_globe=(jy-1)*NX_GLOBE+ix
        kpp_3d_fields%cplwght(ipoint_globe)=cplwght_in(ix,jy)
     ENDDO
  ENDDO  
  status=NF_CLOSE(ncid_cplwght)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
#endif

!  ELSE
!     kpp_3d_fields%cplwght(:) = 2.
!  ENDIF

  ! Do we even need this code anymore if reading a global coupling weight?
!!$  DO ix=1,NX_GLOBE
!!$     DO jy=1,NY_GLOBE
!!$        ipoint_globe=(jy-1)*NX_GLOBE+ix            
!!$        ixx=MIN(ix-ifirst,ilast-ix)
!!$        jyy=MIN(jy-jfirst,jlast-jy)
!!$        IF (ixx .GE. 0 .AND. jyy .GE. 0) THEN
!!$           ! Point is inside coupling domain.  
!!$           ! Set cplwght equal to one (if not already set from NetCDF file) 
!!$           ! to obtain model SSTs.
!!$           kpp_3d_fields%cplwght(ipoint_globe) = MIN(kpp_3d_fields%cplwght(ipoint_globe),1.)
!!$           ipoint=NX*jyy+ixx
!!$           IF (kpp_3d_fields%L_OCEAN(ipoint) .and. kpp_3d_fields%ocdepth(ipoint) .gt. 100) THEN
!!$              kpp_3d_fields%L_OCEAN(ipoint)=.FALSE.
!!$              kpp_3d_fields%cplwght(ipoint_globe)=0
!!$              WRITE(6,*) 'Overwriting coupling mask at'
!!$              WRITE(6,*) 'ixx=',ixx,'jyy=',jyy,'ipoint_globe=',ipoint_globe,'ipoint=',ipoint,'cplwght=',&
!!$                   kpp_3d_fields%cplwght(ipoint_globe),'ocdepth=',kpp_3d_fields%ocdepth(ipoint),&
!!$                   'L_OCEAN=',kpp_3d_fields%L_OCEAN(ipoint)
!!$           ENDIF
!!$        ELSE
!!$           ! Point is outside coupling domain.
!!$           ! Set cplwght equal to a negative value to obtain
!!$           ! climatological SSTs (or persisted SSTs, IF (.NOT. L_UPDCLIM))
!!$           kpp_3d_fields%cplwght(ipoint_globe) = -1.
!!$        ENDIF
!!$     ENDDO
!!$  ENDDO
  
  RETURN  
END SUBROUTINE mckpp_initialize_couplingweight
