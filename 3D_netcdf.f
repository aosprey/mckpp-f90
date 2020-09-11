      SUBROUTINE MY_NCDF_DEF_DIM (
     &     ncid,dimid,dimlen,varid,name,units,delta,long_name)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include <netcdf.inc>
      
      INTEGER ncid,dimid,varid,dimlen
      REAL*4 delta
      CHARACTER*(*) name,units,long_name

      INTEGER status,unit_len,long_len

      unit_len=len(units)
      long_len=len(long_name)
      
      status=NF_DEF_DIM(ncid,name,dimlen,dimid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_DEF_VAR(ncid,name,NF_FLOAT,1,dimid,varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_PUT_ATT_TEXT (ncid,varid,'units',unit_len,units)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_PUT_ATT_REAL (ncid,varid,'spacing',NF_FLOAT,1,delta)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_PUT_ATT_TEXT (ncid,varid,'long_name',
     &     long_len,long_name)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      
      RETURN
      END

      SUBROUTINE MY_NCDF_DEF_VAR (
     &     ncid,varid,ndims,dims,name,units,long_name)


      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
#include <netcdf.inc>


      REAL*4 valid_min,valid_max,miss_val,fill_val
      PARAMETER (valid_max=1.e19,valid_min=-valid_max,
     &     miss_val=1.e20,fill_val=miss_val)

      INTEGER ncid,varid,ndims
      INTEGER dims(ndims)
      CHARACTER*(*) name,units,long_name
      
      INTEGER status,unit_len,long_len

      unit_len=len(units)
      long_len=len(long_name)

      status=NF_DEF_VAR(ncid,name,nf_float,ndims,dims,varid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_PUT_ATT_TEXT (ncid,varid,'long_name',
     &     long_len,long_name)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_PUT_ATT_TEXT (ncid,varid,'units',unit_len,units)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_PUT_ATT_REAL (ncid,varid,'valid_min',nf_float,
     &     1,valid_min)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_PUT_ATT_REAL (ncid,varid,'valid_max',nf_float,
     &     1,valid_max)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_PUT_ATT_REAL (ncid,varid,'missing_value',nf_float,
     &     1,miss_val)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      status=NF_PUT_ATT_REAL (ncid,varid,'_FillValue',nf_float,
     &     1,fill_val)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      RETURN
      END

      SUBROUTINE HANDLE_ERR(status)
      
      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
#include <netcdf.inc>    

      integer status
      CHARACTER*80 err_message

      err_message=NF_STRERROR(status)
      write(nuerr,*) err_message
      
      CALL MIXED_ABORT
      END

