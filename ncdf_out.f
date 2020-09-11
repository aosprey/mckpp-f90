      SUBROUTINE init_output(filename,ncid,kpp_3d_fields,
     +     kpp_const_fields,dt_vec,dt_sing,vec_varids,sing_varids,
     +     zprof_vec,extra_time,L_MEAN,L_INST)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include <netcdf.inc>
! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
c#include "landsea.com"
#include "vert_pgrid.com"
#include "output.com"
#include "times.com"
#include "initialcon.com"
c#include "location.com"

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields
      INTEGER i,j,k,m,ix,iy,ncid,status,extra_time
      CHARACTER*11 varname(N_VAROUTS),singname(N_SINGOUTS)
      CHARACTER*50 longname(N_VAROUTS),singlong(N_SINGOUTS)
      CHARACTER(LEN=*) filename
      CHARACTER*15 units(N_VAROUTS),singunits(N_SINGOUTS)
      CHARACTER*6 type
      CHARACTER*7 time_name
      CHARACTER*4 z_name
      REAL*4, allocatable :: time_out(:)
      INTEGER :: vec_varids(N_VAROUTS),sing_varids(N_SINGOUTS),
     &     dt_vec(N_VAROUTS),dt_sing(N_VAROUTS),n_tdims,n_zdims,n_ddims,
     &     vert_dimid(N_VAROUTS),dims(4),zflag(N_VAROUTS), ! 1 for z-levels, 0 for d-levels
     &     time_dimids(N_VAROUTS+N_SINGOUTS), ! Contains a value for each diagnostic
     &     time_varids(N_VAROUTS+N_SINGOUTS), ! Contains only unique time variables
     &     dt_timeids(N_VAROUTS+N_SINGOUTS),
     &     zprof_vec(N_VAROUTS),zprof_zids(N_VAROUTS),
     &     zprof_dids(N_VAROUTS),z_dimids(N_VAROUTS),
     &     z_varids(N_VAROUTS),d_varids(N_VAROUTS)
      REAL*4 ZOUT(NZP1),ALON(NX),ALAT(NY)
      REAL*4 delta
      LOGICAL tdim_found,zdim_found,L_MEAN,L_INST ! If L_MEAN, time values are output as midpoints
                                                  ! If L_INST, an extra value is added for initial conditions
                                                  !     if .NOT. L_RESTART.
      DATA varname/'u','v','T','S','B',
     &     'wu','wv','wT','wS','wB',
     &     'wTnt',
     &     'difm','dift','difs',
     &     'rho','cp','scorr','Rig','dbloc','Shsq','tinc_fcorr',
     &     'fcorr_z','sinc_fcorr','tinc_ekadv','sinc_ekadv'/
      DATA zflag/1,1,1,1,1,
     &     0,0,0,0,0,
     &     0,
     &     0,0,0,
     &     1,1,1,1,1,1,1,
     &     1,1,1,1/
      DATA longname/
     &     'Zonal velocity',
     &     'Meridional velocity',
     &     'Temperature',
     &     'Salinity',
     &     'Buoyancy',
     &     'Turbulent Zonal Velocity Flux',
     &     'Turbulent Meridional Velocity Flux',
     &     'Turbulent Temperature Flux',
     &     'Turbulent Salinity Flux',
     &     'Turbulent Buoyancy Flux',
     &     'Non-Turbulent Temperature Flux',
     &     'Diffusion Coefficient (Momentum)',
     &     'Diffusion Coefficient (Temperature)',
     &     'Diffusion Coefficient (Salinity)',
     &     'Density',
     &     'Specific Heat Capacity',
     &     'Salinity correction (with depth)',
     &     'Local Richardson Number in kpp.f',
     &     'Local delta buoyancy in kpp.f',
     &     'Local shear-squared term in kpp.f',
     &     'Temperature increment flux correction',
     &     'Heat correction as flux (dT/dt*rho*cp)',
     &     'Salinity increment flux correction',
     &     'Temperature increment from Ekman pumping',
     &     'Salinity increment from Ekman pumping'/
      DATA units/
     &     'm/s',
     &     'm/s',
     &     'degC',
     &     'o/oo',
     &     'm/s^2',
     &     'm^2/s^2',
     &     'm^2/s^2',
     &     'degC m/s',
     &     'o/oo m/s',
     &     'm^2/s^3',
     &     'degC m/s',
     &     'm^2/s',
     &     'm^2/s',
     &     'm^2/s',
     &     'kg/m^3',
     &     'J/kg/K',
     &     'o/oo/s',
     &     'unitless',
     &     'm/s^2',
     &     'm^2/s^2',
     &     'K/timestep',
     &     'W/m^3',
     &     'o/oo/timestep',
     &     'K/timestep',
     &     'o/oo/timestep'/
      DATA singname/
     &     'hmix',
     &     'fcorr',
     &     'taux_in',
     &     'tauy_in',
     &     'solar_in',
     &     'nsolar_in',
     &     'PminusE_in',
     &     'cplwght',
     &     'freeze_flag',
     &     'comp_flag',
     &     'dampu_flag',
     &     'dampv_flag',
     &     'fcorr_nsol',
     &     'hekman',
     &     'runoff_incr'/
#ifdef COUPLE
#ifdef OASIS2
      DATA type /'OASIS2'/
#else
#ifdef OASIS3
      DATA type /'OASIS3'/
#else
#ifdef GFS
      DATA type /' GFS  '/
#endif /*GFS*/
#endif /*OASIS3*/
#endif /*OASIS2*/
#else
      DATA type /'netCDF'/
#endif /*COUPLE*/

      DATA singlong/
     &     'Mixed Layer Depth',
     &     'Flux Correction',
     &     'Zonal wind stress from ',
     &     'Meridional wind stress from ',
     &     'Solar from ',
     &     'Non-solar from ',
     &     'P minus E from ',
     &     'Coupling weight',
     &     'Fraction of levels below freezing',
     &     'Number of integrations (<0 = isothermal reset)',
     &     'Fraction of levels with ui~u**2',
     &     'Fraction of levels with vi~v**2',
     &     'Non-solar heat flux restoring term',
     &     'Depth of Ekman layer',
     &     'Increment to freshwater flux from river runoff'/
      DATA singunits/
     &     'm',
     &     'W/m^2',
     &     'N/m^2',
     &     'N/m^2',
     &     'W/m^2',
     &     'W/m^2',
     &     'mm/s',
     &     'none',
     &     'fraction',
     &     'unitless',
     &     'fraction',
     &     'fraction',
     &     'W/m^2',
     &     'm',
     &     'kg/m^2/s' /

      nout=1
      singlong(3:7)=singlong(3:7)//type

      DO ix=1,nx
         alon(ix)=kpp_3d_fields%dlon(ix)
      ENDDO
      DO iy=1,ny
         alat(iy)=kpp_3d_fields%dlat((iy-1)*nx+1)
      ENDDO

      status=NF_CREATE(filename, nf_clobber, ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      WRITE(nuout,*) 'Output file ',filename,' created successfully.'

      delta=0.0
      IF (NX .GT. 1 ) delta=alon(2)-alon(1)
      CALL MY_NCDF_DEF_DIM (
     &     ncid,londim,nx,lon_id,'longitude','deg',delta,' ')
      delta=0.0
      IF (NY .GT. 1) delta=alat(2)-alat(1)
      CALL MY_NCDF_DEF_DIM (
     &     ncid,latdim,ny,lat_id,'latitude','deg',delta,' ')
      CALL MY_NCDF_DEF_DIM (
     &     ncid,hdim,NZP1,h_id,'h','m',0.0,'Layer Thickness')
      delta=dtout

      n_tdims=0
      n_ddims=0
      n_zdims=0
      IF (dt_vec(1).gt.0) THEN
         n_tdims=1
	      dt_timeids(1)=dt_vec(1)
         WRITE(6,*) 'dt_vec time_1'
         CALL MY_NCDF_DEF_DIM(ncid,time_dimids(1),NF_UNLIMITED,
     +        time_varids(1),'time_1','days',delta,' ')
         j=NZP1
         IF (zprof_vec(1) .gt. 0) j=zprofs_nvalid(zprof_vec(1))
!         WRITE(6,*) j,NZP1,zprofs_nvalid(j)
         IF (zflag(1).eq.0) THEN
            n_ddims=1
            CALL MY_NCDF_DEF_DIM(ncid,z_dimids(1),
     +           j,d_varids(1),'d_1','m',0.0,'Depth of interfaces')
            zprof_dids(1)=zprof_vec(1)
         ELSEIF (zflag(1).eq.1) THEN
            n_zdims=1
            CALL MY_NCDF_DEF_DIM(ncid,z_dimids(1),
     +           j,z_varids(1),'z_1','m',0.0,'')
            zprof_zids(1)=zprof_vec(1)
         ELSE
            WRITE(6,*) 'Incorrect value of zflag for '//longname(k)
            CALL MIXED_ABORT
         ENDIF
      ENDIF
      WRITE(6,*) 'Doing VAROUTS'
      DO j=2,N_VAROUTS
         WRITE(6,*) j,dt_vec(j)
         IF (dt_vec(j) .gt. 0) THEN
            tdim_found=.FALSE.
            zdim_found=.FALSE.
            DO k=1,j-1
               IF (dt_vec(j) .eq. dt_vec(k)) THEN
                  tdim_found=.TRUE.
                  time_dimids(j)=time_dimids(k)
               ENDIF
               IF (zprof_vec(j).eq.zprof_vec(k) .and.
     +              zflag(j) .eq. zflag(k) .and. dt_vec(k).gt.0) THEN
                  zdim_found=.TRUE.
                  z_dimids(j)=z_dimids(k)
               ENDIF
            ENDDO
            IF (.NOT.tdim_found) THEN
               n_tdims=n_tdims+1
               dt_timeids(n_tdims)=dt_vec(j)
               WRITE(time_name,'(A5,I0)') 'time_',n_tdims
               WRITE(6,*) time_name
               IF (n_tdims .eq. 1) THEN
                  CALL MY_NCDF_DEF_DIM(ncid,time_dimids(j),NF_UNLIMITED,
     +                 time_varids(1),'time_1','days',delta,' ')
                  WRITE(6,*) 'Made time_1'
               ELSE
                  CALL MY_NCDF_DEF_DIM (ncid,time_dimids(j),
     +                 ndt_per_file/dt_vec(j)+extra_time,
     &                 time_varids(n_tdims),time_name,'days',delta,' ')
               ENDIF
            ENDIF
            IF (.NOT.zdim_found) THEN
               m=NZP1
               IF (zprof_vec(j).gt.0) m=zprofs_nvalid(zprof_vec(j))
               IF (zflag(k).eq.0) THEN
                  n_ddims=n_ddims+1
                  WRITE(z_name,'(A2,I0)') 'd_',n_ddims
                  WRITE(6,*) z_name,m
                  CALL MY_NCDF_DEF_DIM(ncid,z_dimids(j),
     +                 m,d_varids(n_ddims),z_name,'m',0.0,
     +                 'Depth of interfaces')
                  zprof_dids(n_ddims)=zprof_vec(j)
               ELSEIF (zflag(k).eq.1) THEN
                  n_zdims=n_zdims+1               
                  WRITE(z_name,'(A2,I0)') 'z_',n_zdims
                  WRITE(6,*) z_name,m
                  CALL MY_NCDF_DEF_DIM(ncid,z_dimids(j),
     +                 m,z_varids(n_zdims),z_name,'m',0.0,' ')
                  zprof_zids(n_zdims)=zprof_vec(j)
               ELSE
                  WRITE(6,*)'Incorrect value of zflag for '//longname(j)
                  CALL MIXED_ABORT
               ENDIF
            ENDIF
         ELSE
            time_dimids(j)=0
            z_dimids(j)=0
         ENDIF
      ENDDO
      WRITE(6,*) 'Finished VAROUTS'

      WRITE(6,*) 'Doing SINGOUTS'
      DO j=1,N_SINGOUTS
         IF (dt_sing(j) .gt. 0) THEN
            tdim_found=.FALSE.
            DO k=1,N_VAROUTS
               IF (dt_sing(j) .eq. dt_vec(k)) THEN
                  tdim_found=.TRUE.
                  time_dimids(j+N_VAROUTS)=time_dimids(k)
               ENDIF
            ENDDO
            IF (j .gt. 1 .and. .NOT. tdim_found) THEN
               DO k=1,j-1
                  IF (dt_sing(j) .eq. dt_sing(k)) THEN
                     tdim_found=.TRUE.
                     time_dimids(j+N_VAROUTS)=time_dimids(k+N_VAROUTS)
                  ENDIF
               ENDDO
            ENDIF
            IF (.NOT. tdim_found) THEN
               n_tdims=n_tdims+1
               dt_timeids(n_tdims)=dt_sing(j)
               WRITE(time_name,'(A5,I0)') 'time_',n_tdims
               WRITE(6,*) time_name
               IF (n_tdims .eq. 1) THEN
                  CALL MY_NCDF_DEF_DIM(ncid,time_dimids(j+N_VAROUTS),
     +                 NF_UNLIMITED,time_varids(1),'time_1','days',
     +                 delta,' ')
                  WRITE(6,*) 'Made time_1 in singouts'
               ELSE
                  CALL MY_NCDF_DEF_DIM(ncid,time_dimids(j+N_VAROUTS),
     +                 ndt_per_file/dt_sing(j)+extra_time,
     +                 time_varids(n_tdims),time_name,'days',delta,' ')
               ENDIF
            ENDIF
         ELSE
            time_dimids(j+N_VAROUTS)=0
         ENDIF
      ENDDO

      dims(1)=londim
      dims(2)=latdim
      DO k=1,N_VAROUTS
         dims(3)=z_dimids(k)
         dims(4)=time_dimids(k)
         IF (dt_vec(k) .gt. 0)
     &        CALL MY_NCDF_DEF_VAR (
     &        ncid,vec_varids(k),4,dims,varname(k),units(k),
     &        longname(k))
      ENDDO

      DO k=1,N_SINGOUTS
         dims(3)=time_dimids(N_VAROUTS+k)
         IF (dt_sing(k) .gt. 0)
     &        CALL MY_NCDF_DEF_VAR (
     &        ncid,sing_varids(k),3,dims,singname(k),singunits(k),
     &        singlong(k))
      ENDDO

      status=NF_ENDDEF(ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      status=NF_PUT_VAR_REAL(ncid,lon_id,alon)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      status=NF_PUT_VAR_REAL(ncid,lat_id,alat)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      DO j=1,n_zdims
         m=1
         DO k=1,NZP1
            IF (zprofs_mask(k,zprof_zids(j))) THEN
               zout(m)=kpp_const_fields%zm(k)
               m=m+1
            ENDIF
         ENDDO
         status=NF_PUT_VAR_REAL(ncid,z_varids(j),ZOUT)
         IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
      ENDDO
      DO j=1,n_ddims
         m=1
         DO k=1,NZP1
            IF (zprofs_mask(k,zprof_dids(j))) THEN
               zout(m)=kpp_const_fields%dm(k-1)
               m=m+1
            ENDIF
         ENDDO
         status=NF_PUT_VAR_REAL(ncid,d_varids(j),ZOUT)
         IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
      ENDDO

      DO k=1,NZP1
         ZOUT(k)=kpp_const_fields%hm(k)
      ENDDO
      status=NF_PUT_VAR_REAL(ncid,h_id,ZOUT)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      DO k=1,n_tdims
         allocate(time_out(ndt_per_file/dt_timeids(k)+extra_time))
         DO i=1-extra_time,ndt_per_file/dt_timeids(k)
            time_out(i+extra_time)=i*kpp_const_fields%dto/
     +           kpp_const_fields%spd*dt_timeids(k)+
     +	         kpp_const_fields%time
            IF (L_MEAN)
     +           time_out(i+extra_time)=time_out(i+extra_time)-
     +           0.5*kpp_const_fields%dto/kpp_const_fields%spd*
     +           dt_timeids(k)
         ENDDO
!         WRITE(6,*) 'Writing time_out=',time_out(:)
!         WRITE(6,*) 'ndt_per_file=',ndt_per_file,'time_varids(k)=',
!     +        time_varids(k)
         status=NF_PUT_VARA_REAL(ncid,time_varids(k),
     +        (/1/),(/ndt_per_file/dt_timeids(k)+extra_time/),
     +        time_out)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         deallocate(time_out)
      ENDDO

      write(nuout,*) 'NCDF output file initialized with name ',
     +      filename

      status=NF_CLOSE(ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

!     Reset extra_time parameter if set to one, as this
!     should be used only on the first timestep of a new simulation
!     to leave enough space to output the initial condition.
      IF (extra_time .eq. 1) extra_time = 0

      RETURN
      END

      SUBROUTINE output_inst(kpp_3d_fields,kpp_const_fields,diag_num,
     +     varid,zprof,ntout)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include <netcdf.inc>
! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
#include "output.com"

#include "times.com"
#include "ocn_advec.com"
#include "landsea.com"
#include "relax_3d.com"
#include "couple.com"
#include "fcorr_in.com"
#include "sfcorr_in.com"

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

      REAL*4, allocatable :: varout(:,:,:), singout(:,:),
     +     temp_2d(:,:),temp_1d(:),temp_zprof(:,:)
      REAL*4 TOUT

      INTEGER start(4),count(4),status
      INTEGER i,j,k,ivar,diag_num,varid,ntout
      INTEGER ix,iy,ipt,zprof

      count(1)=NX
      count(2)=NY
      count(3)=zprofs_nvalid(zprof)
      count(4)=1

      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=ntout

!      write(nuout,*) 'Writing output for diagnostic ',diag_num,
!     + 'to start=',start,'count=',count,'zprof=',zprof

c     NPK 25/2/08.
c     Stop opening and closing output files with each flush of output.
c     Just do it once at the beginning/end of the simulation.
c      call output_open
c      TOUT=kpp_const_fields%time
c      status=NF_PUT_VAR1_REAL(ncid_out,time_id,nout,TOUT)


      IF (diag_num .le. N_VAROUTS) THEN
         allocate(varout(NX,NY,zprofs_nvalid(zprof)))
         allocate(temp_2d(NPTS,NZP1))
         SELECT CASE (diag_num)
         CASE (1)
            temp_2d(:,:)=kpp_3d_fields%U(:,:,1)
         CASE (2)
            temp_2d(:,:)=kpp_3d_fields%U(:,:,2)
         CASE (3)
            temp_2d(:,:)=kpp_3d_fields%X(:,:,1)
         CASE (4)
            DO k=1,NZP1
               temp_2d(:,k)=kpp_3d_fields%X(:,k,2)+
     +              kpp_3d_fields%Sref(:)
            ENDDO
         CASE(5)
            temp_2d(:,:)=kpp_3d_fields%buoy(:,1:NZP1)
         CASE(6)
            temp_2d(:,:)=kpp_3d_fields%wU(:,0:NZ,1)
         CASE(7)
            temp_2d(:,:)=kpp_3d_fields%wU(:,0:NZ,2)
         CASE(8)
            temp_2d(:,:)=kpp_3d_fields%wX(:,0:NZ,1)
         CASE(9)
            temp_2d(:,:)=kpp_3d_fields%wX(:,0:NZ,2)
         CASE(10)
            temp_2d(:,:)=kpp_3d_fields%wX(:,0:NZ,NSP1)
         CASE(11)
            temp_2d(:,:)=kpp_3d_fields%wXNT(:,0:NZ,1)
         CASE(12)
            temp_2d(:,1)=0.0
            temp_2d(:,2:NZP1)=kpp_3d_fields%difm(:,1:NZ)
         CASE(13)
            temp_2d(:,1)=0.0
            temp_2d(:,2:NZP1)=kpp_3d_fields%dift(:,1:NZ)
         CASE(14)
            temp_2d(:,1)=0.0
            temp_2d(:,2:NZP1)=kpp_3d_fields%difs(:,1:NZ)
         CASE(15)
            temp_2d(:,:)=kpp_3d_fields%rho(:,1:NZP1)
         CASE(16)
            temp_2d(:,:)=kpp_3d_fields%cp(:,1:NZP1)
         CASE(17)
            temp_2d(:,:)=kpp_3d_fields%scorr(:,:)
         CASE(18)
            temp_2d(:,:)=kpp_3d_fields%Rig(:,:)
         CASE(19)
            temp_2d(:,1:NZ)=kpp_3d_fields%dbloc(:,1:NZ)
            temp_2d(:,NZP1)=0.0
         CASE(20)
            temp_2d(:,:)=kpp_3d_fields%Shsq(:,:)
         CASE(21)
            temp_2d(:,:)=kpp_3d_fields%Tinc_fcorr(:,:)
         CASE(22)
            temp_2d(:,:)=kpp_3d_fields%ocnTcorr(:,:)
         CASE(23)
            temp_2d(:,:)=kpp_3d_fields%Sinc_fcorr(:,:)
         CASE(24)
            temp_2d(:,:)=kpp_3d_fields%Tinc_ekadv(:,:)
         CASE(25)
            temp_2d(:,:)=kpp_3d_fields%Sinc_ekadv(:,:)
         CASE DEFAULT
            WRITE(6,*) 'You need to add more outputs in OUTPUT_INST'
         END SELECT

         IF (zprof .gt. 0) THEN
            j=1
            allocate(temp_zprof(NPTS,zprofs_nvalid(zprof)))
            DO i=1,NZP1
               IF (zprofs_mask(i,zprof)) THEN
                  temp_zprof(:,j)=temp_2d(:,i)
                  j=j+1
               ENDIF
            ENDDO
            CALL REFORMAT_MASK_OUTPUT_2D(temp_zprof,zprofs_nvalid(zprof)
     +           ,kpp_3d_fields%L_OCEAN,missval,varout)
         ELSE
            CALL REFORMAT_MASK_OUTPUT_2D(temp_2d,NZP1,
     +           kpp_3d_fields%L_OCEAN,missval,varout)
         ENDIF

         status=NF_PUT_VARA_REAL(
     &        ncid_out,varid,start,count,VAROUT)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      ELSE
         allocate(singout(NX,NY))
         allocate(temp_1d(NPTS))
         start(3)=start(4)
         count(3)=1
         SELECT CASE (diag_num-N_VAROUTS)
         CASE (1)
            temp_1d(:)=kpp_3d_fields%hmix(:)
         CASE (2)
            temp_1d(:)=kpp_3d_fields%fcorr(:)
         CASE (3)
            temp_1d(:)=kpp_3d_fields%sflux(:,1,5,0)
         CASE (4)
            temp_1d(:)=kpp_3d_fields%sflux(:,2,5,0)
         CASE (5)
            temp_1d(:)=kpp_3d_fields%sflux(:,3,5,0)
         CASE (6)
            temp_1d(:)=kpp_3d_fields%sflux(:,4,5,0)
         CASE (7)
            temp_1d(:)=kpp_3d_fields%sflux(:,6,5,0)
         CASE (8)
            DO ix=ifirst,ilast
               DO iy=jfirst,jlast
                  ipt=(iy-1)*NX_GLOBE+ix
                  temp_1d((iy-jfirst)*NX+ix-ifirst+1)=
     +                 kpp_3d_fields%cplwght(ipt)
               ENDDO
            ENDDO
         CASE (9)
            temp_1d(:)=kpp_3d_fields%freeze_flag(:)
         CASE (10)
            temp_1d(:)=kpp_3d_fields%reset_flag(:)
         CASE (11)
            temp_1d(:)=kpp_3d_fields%dampu_flag(:)
         CASE (12)
            temp_1d(:)=kpp_3d_fields%dampv_flag(:)
         CASE (13)
            temp_1d(:)=kpp_3d_fields%fcorr_nsol(:)
         CASE (14)
            temp_1d(:)=kpp_3d_fields%hekman(:)
         CASE (15)
            temp_1d(:)=kpp_3d_fields%runoff_incr(:)
         CASE DEFAULT
            WRITE(6,*) 'You need to add more outputs in '//
     +           'OUTPUT_INST'
         END SELECT

         CALL REFORMAT_MASK_OUTPUT_1D(temp_1d,kpp_3d_fields%L_OCEAN,
     +        missval,singout)

         status=NF_PUT_VARA_REAL(
     &        ncid_out,varid,start,count,SINGOUT)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
      ENDIF
      start(3)=1

      ntout=ntout+1
c     NPK 25/2/08
c     Stop opening and closing output files with each flush of output.
c     Just do it once at the beginning/end of the simulation.
c      call output_close
!      WRITE(nuout,*) ' Output successfully written'
      RETURN
      END

      SUBROUTINE reformat_mask_output_1d(oned_in,mask,missval,
     +     twod_out)
      IMPLICIT NONE
#include "parameter.inc"

      REAL*4,intent(in) :: oned_in(NPTS),missval
      REAL*4,intent(out) :: twod_out(NX,NY)
      LOGICAL,intent(in) :: mask(NPTS)
      INTEGER :: i,j,ipt

      DO i=1,NX
         DO j=1,NY
            ipt=(j-1)*NX+i
            IF (mask(ipt)) THEN
               twod_out(i,j)=oned_in(ipt)
            ELSE
               twod_out(i,j)=missval
            ENDIF
         ENDDO
      ENDDO
      RETURN
      END SUBROUTINE reformat_mask_output_1d

      SUBROUTINE reformat_mask_output_2d(twod_in,nz_in,mask,missval,
     +     threed_out)

      IMPLICIT NONE
#include "parameter.inc"

      INTEGER,intent(in) :: nz_in
      REAL*4,intent(in) :: twod_in(NPTS,nz_in),missval
      REAL*4,intent(out) :: threed_out(NX,NY,nz_in)
      LOGICAL,intent(in) :: mask(NPTS)
      INTEGER :: i,j,ipt

      DO i=1,NX
         DO j=1,NY
            ipt=(j-1)*NX+i
            IF (mask(ipt)) THEN
               threed_out(i,j,:)=twod_in(ipt,:)
            ELSE
               threed_out(i,j,:)=missval
            ENDIF
         ENDDO
      ENDDO

      RETURN
      END SUBROUTINE reformat_mask_output_2d

      SUBROUTINE write_means(kpp_3d_fields,kpp_const_fields,
     +     VEC_mean,SCLR_mean,diag_num,varid,zprof,ntout)
      IMPLICIT NONE

      INTEGER nuout,nuerr
      PARAMETER(nuout=6,nuerr=0)

#include <netcdf.inc>
! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
#include "output.com"
#include "landsea.com"
#include "times.com"
#include "timocn.com"
#include "ocn_advec.com"
#include "couple.com"
#include "constants.com"

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

      REAL :: VEC_mean(NPTS,NZP1,NVEC_MEAN),
     + SCLR_mean(NPTS,NSCLR_MEAN)
      REAL*4,allocatable :: VAROUT(:,:,:), SINGOUT(:,:),
     +     temp_2d(:,:),temp_zprof(:,:)
      REAL*4 TOUT

      INTEGER i,j,ivar,ipt,ix,iy,start(4),count(4),k,status,
     +     diag_num,varid,ntout,mean_num,zprof


c      TOUT=kpp_const_fields%time-(ndtout_mean*kpp_const_fields%dtsec/
c     +     (ndtocn*kpp_const_fields%spd))*0.5
c      status=NF_PUT_VAR1_REAL(mean_ncid_out,time_id,nout_mean,TOUT)
c      IF (status.NE.NF_NOERR) CALL HANDLE_ERR(status)
c      write(nuout,*) 'Writing mean output at timestep ',ntime+nstart,
c     +     ' Time=',TOUT,' time_id = ',time_id

      count(1)=NX
      count(2)=NY
      count(3)=zprofs_nvalid(zprof)
      count(4)=1

      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=ntout

      IF (diag_num .le. N_VAROUTS) THEN
         allocate(VAROUT(NX,NY,NZP1))
         allocate(temp_2d(NPTS,NZP1))
         mean_num=1
         DO i=1,diag_num-1
            IF (ndt_varout_mean(i) .gt. 0)
     +           mean_num=mean_num+1
         ENDDO
         SELECT CASE (diag_num)
         CASE (4)
            DO k=1,NZP1
               temp_2d(:,k)=VEC_mean(:,k,mean_num)+
     +              kpp_3d_fields%Sref(:)
            ENDDO
         CASE (12,13,14)
            temp_2d(:,1)=0.
            temp_2d(:,2:NZP1)=VEC_mean(:,2:NZP1,mean_num)
         CASE (18,19)
            temp_2d(:,1:NZ)=VEC_mean(:,1:NZ,mean_num)
            temp_2d(:,NZP1)=0.
         CASE DEFAULT
            temp_2d(:,:)=VEC_mean(:,:,mean_num)
         END SELECT
!         WRITE(6,*) 'In WRITE_MEANS for diag_num=',diag_num,'zprof=',
!     +        zprof
c         WRITE(6,*) 'Calling reformat_mask_output_2d for i=',mean_num

         IF (zprof .gt. 0) THEN
            j=1
            allocate(temp_zprof(NPTS,zprofs_nvalid(zprof)))
            DO i=1,NZP1
               IF (zprofs_mask(i,zprof)) THEN
                  temp_zprof(:,j)=temp_2d(:,i)
                  j=j+1
               ENDIF
            ENDDO
            CALL REFORMAT_MASK_OUTPUT_2D(temp_zprof,zprofs_nvalid(zprof)
     +           ,kpp_3d_fields%L_OCEAN,missval,varout)
         ELSE
            CALL REFORMAT_MASK_OUTPUT_2D(temp_2d,NZP1,
     +           kpp_3d_fields%L_OCEAN,missval,varout)
         ENDIF

!         WRITE(6,*) 'Writing to ncid=',mean_ncid_out,' varid=',
!     +        varid,' with start =',start,' and count =',
!     +        count
         status=NF_PUT_VARA_REAL(
     &        mean_ncid_out,varid,start,count,VAROUT)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         VEC_mean(:,:,mean_num)=0.
         i=i+1
      ELSE
         allocate(SINGOUT(NX,NY))
         start(3)=start(4)
         count(3)=1
         mean_num=1
         DO i=1,(diag_num-N_VAROUTS-1)
            IF (ndt_singout_mean(i) .gt. 0)
     +           mean_num=mean_num+1
         ENDDO
         SINGOUT(:,:) = missval
         DO ix=1,nx
            DO iy=1,ny
               ipt=(iy-1)*nx+ix
               IF (kpp_3d_fields%L_OCEAN(ipt))
     &              SINGOUT(ix,iy)=SCLR_mean(ipt,mean_num)
            ENDDO
         ENDDO
!         WRITE(nuout,*) 'In write_means for singout, i=',mean_num,
!     &        'diag_num=',diag_num
         status=NF_PUT_VARA_REAL(
     &        mean_ncid_out,varid,start,count,SINGOUT)
         IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
         SCLR_mean(:,mean_num)=0.
         i=i+1
      ENDIF
      start(3)=1

c     Increment counter for time dimension of NetCDF file
      ntout=ntout+1

      RETURN
      END

      SUBROUTINE output_range(kpp_3d_fields,kpp_const_fields,
     +     VEC_range,SCLR_range,diag_num,varid,zprof,ntout)
      IMPLICIT NONE

      INTEGER nuout,nuerr
      PARAMETER(nuout=6,nuerr=0)

#include <netcdf.inc>
! Automatically includes parameter.inc!
#include "kpp_3d_type.com"
#include "output.com"
#include "landsea.com"
#include "times.com"
#include "timocn.com"
#include "ocn_advec.com"
#include "couple.com"
#include "constants.com"

      TYPE(kpp_3d_type) :: kpp_3d_fields
      TYPE(kpp_const_type) :: kpp_const_fields

      REAL :: VEC_range(NPTS,NZP1,NVEC_RANGE,2),
     +     SCLR_range(NPTS,NSCLR_RANGE,2)
      REAL*4,allocatable :: VAROUT(:,:,:), SINGOUT(:,:),
     +     temp_2d(:,:),temp_zprof(:,:)
      REAL*4 TOUT

      INTEGER i,j,ivar,ipt,ix,iy,start(4),count(4),k,status,
     +     diag_num,varid,ntout,range_num,zprof

      count(1)=NX
      count(2)=NY
      count(3)=zprofs_nvalid(zprof)
      count(4)=1

      start(1)=1
      start(2)=1
      start(3)=1
      start(4)=ntout

      allocate(SINGOUT(NX,NY))
      allocate(VAROUT(NX,NY,NZP1))
      allocate(temp_2d(NPTS,NZP1))
      DO j=1,2
         IF (diag_num .le. N_VAROUTS) THEN
            range_num=1
            DO i=1,diag_num-1
               IF (ndt_varout_range(i) .gt. 0)
     +              range_num=range_num+1
            ENDDO
            SELECT CASE (diag_num)
            CASE (4)
               DO k=1,NZP1
                  temp_2d(:,k)=VEC_range(:,k,range_num,j)+
     +                 kpp_3d_fields%Sref(:)
               ENDDO
            CASE (12,13,14)
               temp_2d(:,1)=0.
               temp_2d(:,2:NZP1)=VEC_range(:,2:NZP1,range_num,j)
            CASE (18,19)
               temp_2d(:,1:NZ)=VEC_range(:,1:NZ,range_num,j)
               temp_2d(:,NZP1)=0.
            CASE DEFAULT
               temp_2d(:,:)=VEC_range(:,:,range_num,j)
            END SELECT
c            WRITE(6,*) 'Calling reformat_mask_output_2d diag=',diag_num
            IF (zprof .gt. 0) THEN
               k=1
               allocate(temp_zprof(NPTS,zprofs_nvalid(zprof)))
               DO i=1,NZP1
                  IF (zprofs_mask(i,zprof)) THEN
                     temp_zprof(:,k)=temp_2d(:,i)
                     k=k+1
                  ENDIF
               ENDDO
               CALL REFORMAT_MASK_OUTPUT_2D(temp_zprof,
     +              zprofs_nvalid(zprof),kpp_3d_fields%L_OCEAN,
     +              missval,varout)
               deallocate(temp_zprof)
            ELSE
               CALL REFORMAT_MASK_OUTPUT_2D(temp_2d,NZP1,
     +              kpp_3d_fields%L_OCEAN,missval,varout)
            ENDIF
            SELECT CASE (j)
            CASE (1)
!               WRITE(6,*) 'Writing to ncid=',min_ncid_out,' varid=',
!     +              varid,' with start =',start,' and count =',
!     +              count
               status=NF_PUT_VARA_REAL(
     &              min_ncid_out,varid,start,count,VAROUT)
               IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
               VEC_range(:,:,range_num,j)=2e20
            CASE (2)
!               WRITE(6,*) 'Writing to ncid=',max_ncid_out,' varid=',
!     +              varid,' with start =',start,' and count =',
!     +              count
               status=NF_PUT_VARA_REAL(
     &              max_ncid_out,varid,start,count,VAROUT)
               IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
               VEC_range(:,:,range_num,j)=-2e20
            END SELECT
c     i=i+1
         ELSE
            start(3)=start(4)
            count(3)=1
            range_num=1
            DO i=1,(diag_num-N_VAROUTS-1)
               IF (ndt_singout_range(i) .gt. 0)
     +              range_num=range_num+1
            ENDDO
            SINGOUT(:,:) = missval
            DO ix=1,nx
               DO iy=1,ny
                  ipt=(iy-1)*nx+ix
                  IF (kpp_3d_fields%L_OCEAN(ipt))
     &                 SINGOUT(ix,iy)=SCLR_range(ipt,range_num,j)
               ENDDO
            ENDDO
            SELECT CASE (j)
            CASE (1)
               status=NF_PUT_VARA_REAL(
     &              min_ncid_out,varid,start,count,SINGOUT)
               IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
               SCLR_range(:,range_num,j)=2e20
            CASE (2)
               status=NF_PUT_VARA_REAL(
     &              max_ncid_out,varid,start,count,SINGOUT)
               IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)
               SCLR_range(:,range_num,j)=-2e20
            END SELECT
c     i=i+1
         ENDIF
      ENDDO
      start(3)=1

c     Increment counter for time dimension of NetCDF file
      ntout=ntout+1

      RETURN
      END

      SUBROUTINE mean_output(kpp_3d_fields,VEC_mean,SCLR_mean)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
! Automatically includes parameter.inc
#include "kpp_3d_type.com"
#include "output.com"
#include "times.com"
#include "landsea.com"
#include "relax_3d.com"
#include "couple.com"
#include "ocn_advec.com"
#include "fcorr_in.com"
#include "sfcorr_in.com"

      REAL,intent(inout) :: VEC_mean(NPTS,NZP1,NVEC_MEAN),
     +  SCLR_mean(NPTS,NSCLR_MEAN)
      REAL :: field(NPTS,NZP1),vec(NPTS)

      TYPE(kpp_3d_type) :: kpp_3d_fields
      INTEGER i,j,k,ivar,ix,iy,ipt,ipt_globe,upper_limit,lower_limit

!      allocate(field(NPTS,NZP1))
!      allocate(vec(NPTS))

      WRITE(6,*) 'In mean_output ',N_VAROUTS
      i=1
      DO ivar=1,N_VAROUTS         
!         WRITE(6,*) ivar
         IF (ndt_varout_mean(ivar) .gt. 0) THEN
c            WRITE(6,*) 'Computing means for ivar = ',ivar,'i=',i
c            WRITE(6,*) 'ndt_varout_mean(ivar)=',ndt_varout_mean(ivar)
            SELECT CASE (ivar)
            CASE(1)
               field(:,:)=kpp_3d_fields%U(:,:,1)
            CASE(2)
               field(:,:)=kpp_3d_fields%U(:,:,2)
            CASE(3)
               field(:,:)=kpp_3d_fields%X(:,:,1)
            CASE(4)
               field(:,:)=kpp_3d_fields%X(:,:,2)
            CASE(5)
               field(:,:)=kpp_3d_fields%buoy(:,1:NZP1)
            CASE(6)
               field(:,:)=kpp_3d_fields%wu(:,0:NZ,1)
            CASE(7)
               field(:,:)=kpp_3d_fields%wu(:,0:NZ,2)
            CASE(8)
               field(:,:)=kpp_3d_fields%wx(:,0:NZ,1)
            CASE(9)
               field(:,:)=kpp_3d_fields%wx(:,0:NZ,2)
            CASE(10)
               field(:,:)=kpp_3d_fields%wx(:,0:NZ,NSP1)
            CASE(11)
               field(:,:)=kpp_3d_fields%wXNT(:,0:NZ,1)
            CASE(12)
               field(:,:)=kpp_3d_fields%difm(:,1:NZP1)
            CASE(13)
               field(:,:)=kpp_3d_fields%dift(:,1:NZP1)
            CASE(14)
               field(:,:)=kpp_3d_fields%difs(:,1:NZP1)
            CASE(15)
               field(:,:)=kpp_3d_fields%rho(:,1:NZP1)
            CASE(16)
               field(:,:)=kpp_3d_fields%cp(:,1:NZP1)
            CASE(17)
               field(:,:)=kpp_3d_fields%scorr(:,:)
            CASE(18)
               field(:,:)=kpp_3d_fields%Rig(:,:)
            CASE(19)
               field(:,1:NZ)=kpp_3d_fields%dbloc(:,1:NZ)
               field(:,NZP1)=0.
            CASE(20)
               field(:,:)=kpp_3d_fields%Shsq(:,:)
            CASE(21)
               field(:,:)=kpp_3d_fields%Tinc_fcorr(:,:)
            CASE(22)
               field(:,:)=kpp_3d_fields%ocnTcorr(:,:)
            CASE(23)
               field(:,:)=kpp_3d_fields%sinc_fcorr(:,:)
            CASE(24)
               field(:,:)=kpp_3d_fields%Tinc_ekadv(:,:)
            CASE(25)
               field(:,:)=kpp_3d_fields%Sinc_ekadv(:,:)
            END SELECT
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields,ndt_varout_mean)
!$OMP& SHARED(VEC_mean,i,field,ivar)
!$OMP DO SCHEDULE(static)
#endif
            DO j=1,NPTS
               IF (kpp_3d_fields%L_OCEAN(j)) THEN
                  DO k=1,NZP1
                     VEC_mean(j,k,i)=field(j,k) / ndt_varout_mean(ivar)
     +                    + VEC_mean(j,k,i)
                  ENDDO
               ENDIF
            ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            i=i+1
         ENDIF
      ENDDO
      i=1
      DO ivar=1,N_SINGOUTS
!         WRITE(6,*) 'Means with ivar=',ivar,'i=',i
         IF (ndt_singout_mean(ivar) .gt. 0) THEN
            SELECT CASE (ivar)
            CASE(1)
               vec(:)=kpp_3d_fields%hmix(:)
            CASE(2)
               vec(:)=kpp_3d_fields%fcorr(:)
            CASE(3)
               vec(:)=kpp_3d_fields%sflux(:,1,5,0)
            CASE(4)
               vec(:)=kpp_3d_fields%sflux(:,2,5,0)
            CASE(5)
               vec(:)=kpp_3d_fields%sflux(:,3,5,0)
            CASE(6)
               vec(:)=kpp_3d_fields%sflux(:,4,5,0)
            CASE(7)
               vec(:)=kpp_3d_fields%sflux(:,6,5,0)
            CASE(8)
               DO ix=ifirst,ilast
                  DO iy=jfirst,jlast
                     ipt_globe=(iy-1)*NX_GLOBE+ix
                     ipt=(iy-jfirst)*nx+(ix-ifirst+1)
                     vec(ipt)=kpp_3d_fields%cplwght(ipt_globe)
                  ENDDO
               ENDDO
            CASE(9)
               vec(:)=kpp_3d_fields%freeze_flag(:)
            CASE(10)
               vec(:)=kpp_3d_fields%reset_flag(:)
            CASE(11)
               vec(:)=kpp_3d_fields%dampu_flag(:)
            CASE(12)
               vec(:)=kpp_3d_fields%dampv_flag(:)
            CASE(13)
               vec(:)=kpp_3d_fields%fcorr_nsol(:)
            CASE(14)
               vec(:)=kpp_3d_fields%hekman(:)
            CASE(15)
               vec(:)=kpp_3d_fields%runoff_incr(:)
            END SELECT
            DO j=1,NPTS
               IF (kpp_3d_fields%L_OCEAN(j))
     +              SCLR_mean(j,i)=vec(j)/ndt_singout_mean(ivar)+
     +              SCLR_mean(j,i)
            ENDDO
            i=i+1
         ENDIF
      ENDDO
!      deallocate(field)
!      deallocate(vec)

      RETURN
      END

      SUBROUTINE range_output(kpp_3d_fields,VEC_range,SCLR_range)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)
! Automatically includes parameter.inc
#include "kpp_3d_type.com"
#include "output.com"
#include "times.com"
#include "landsea.com"
#include "relax_3d.com"
#include "couple.com"
#include "ocn_advec.com"
#include "fcorr_in.com"
#include "sfcorr_in.com"

      REAL,intent(inout) :: VEC_range(NPTS,NZP1,NVEC_RANGE,2),
     +  SCLR_range(NPTS,NSCLR_RANGE,2)
      REAL :: field(NPTS,NZP1),vec(NPTS)

      TYPE(kpp_3d_type) :: kpp_3d_fields
      INTEGER i,j,k,ivar,ix,iy,ipt,ipt_globe

!      allocate(field(NPTS,NZP1))
!      allocate(vec(NPTS))

      i=1
      DO ivar=1,N_VAROUTS
         IF (ndt_varout_range(ivar) .gt. 0) THEN
            SELECT CASE (ivar)
            CASE(1)
               field(:,:)=kpp_3d_fields%U(:,:,1)
            CASE(2)
               field(:,:)=kpp_3d_fields%U(:,:,2)
            CASE(3)
               field(:,:)=kpp_3d_fields%X(:,:,1)
            CASE(4)
               field(:,:)=kpp_3d_fields%X(:,:,2)
            CASE(5)
               field(:,:)=kpp_3d_fields%buoy(:,1:NZP1)
            CASE(6)
               field(:,:)=kpp_3d_fields%wu(:,0:NZ,1)
            CASE(7)
               field(:,:)=kpp_3d_fields%wu(:,0:NZ,2)
            CASE(8)
               field(:,:)=kpp_3d_fields%wx(:,0:NZ,1)
            CASE(9)
               field(:,:)=kpp_3d_fields%wx(:,0:NZ,2)
            CASE(10)
               field(:,:)=kpp_3d_fields%wx(:,0:NZ,NSP1)
            CASE(11)
               field(:,:)=kpp_3d_fields%wXNT(:,0:NZ,1)
            CASE(12)
               field(:,:)=kpp_3d_fields%difm(:,1:NZP1)
            CASE(13)
               field(:,:)=kpp_3d_fields%dift(:,1:NZP1)
            CASE(14)
               field(:,:)=kpp_3d_fields%difs(:,1:NZP1)
            CASE(15)
               field(:,:)=kpp_3d_fields%rho(:,1:NZP1)
            CASE(16)
               field(:,:)=kpp_3d_fields%cp(:,1:NZP1)
            CASE(17)
               field(:,:)=kpp_3d_fields%scorr(:,:)
            CASE(18)
               field(:,:)=kpp_3d_fields%Rig(:,:)
            CASE(19)
               field(:,1:NZ)=kpp_3d_fields%dbloc(:,1:NZ)
               field(:,NZP1)=0.
            CASE(20)
               field(:,:)=kpp_3d_fields%Shsq(:,:)
            CASE(21)
               field(:,:)=kpp_3d_fields%Tinc_fcorr(:,:)
            CASE(22)
               field(:,:)=kpp_3d_fields%ocnTcorr(:,:)
            CASE(23)
               field(:,:)=kpp_3d_fields%sinc_fcorr(:,:)
            CASE(24)
               field(:,:)=kpp_3d_fields%tinc_fcorr(:,:)
            CASE(25)
               field(:,:)=kpp_3d_fields%sinc_fcorr(:,:)
            END SELECT
#ifdef OPENMP
!$OMP PARALLEL DEFAULT(private) SHARED(kpp_3d_fields)
!$OMP& SHARED(VEC_range,i,field)
!$OMP DO SCHEDULE(static)
#endif
            DO j=1,NPTS
               IF (kpp_3d_fields%L_OCEAN(j)) THEN
                  DO k=1,NZP1
                     IF (field(j,k) .lt. VEC_range(j,k,i,1))
     +                    VEC_range(j,k,i,1)=field(j,k)
                     IF (field(j,k) .gt. VEC_range(j,k,i,2))
     +                    VEC_range(j,k,i,2)=field(j,k)
                  ENDDO
               ENDIF
            ENDDO
#ifdef OPENMP
!$OMP END DO
!$OMP END PARALLEL
#endif
            i=i+1
         ENDIF
      ENDDO
      i=1
      DO ivar=1,N_SINGOUTS
!         WRITE(6,*) 'Means with ivar=',ivar,'i=',i
         IF (ndt_singout_range(ivar) .gt. 0) THEN
            SELECT CASE (ivar)
            CASE(1)
               vec(:)=kpp_3d_fields%hmix(:)
            CASE(2)
               vec(:)=kpp_3d_fields%fcorr(:)
            CASE(3)
               vec(:)=kpp_3d_fields%sflux(:,1,5,0)
            CASE(4)
               vec(:)=kpp_3d_fields%sflux(:,2,5,0)
            CASE(5)
               vec(:)=kpp_3d_fields%sflux(:,3,5,0)
            CASE(6)
               vec(:)=kpp_3d_fields%sflux(:,4,5,0)
            CASE(7)
               vec(:)=kpp_3d_fields%sflux(:,6,5,0)
            CASE(8)
               DO ix=ifirst,ilast
                  DO iy=jfirst,jlast
                     ipt_globe=(iy-1)*NX_GLOBE+ix
                     ipt=(iy-jfirst)*nx+(ix-ifirst+1)
                     vec(ipt)=kpp_3d_fields%cplwght(ipt_globe)
                  ENDDO
               ENDDO
            CASE(9)
               vec(:)=kpp_3d_fields%freeze_flag(:)
            CASE(10)
               vec(:)=kpp_3d_fields%reset_flag(:)
            CASE(11)
               vec(:)=kpp_3d_fields%dampu_flag(:)
            CASE(12)
               vec(:)=kpp_3d_fields%dampv_flag(:)
            CASE(13)
               vec(:)=kpp_3d_fields%fcorr_nsol(:)
            CASE(14)
               vec(:)=kpp_3d_fields%hekman(:)
            CASE(15)
               vec(:)=kpp_3d_fields%runoff_incr(:)
            END SELECT
            DO j=1,NPTS
               IF (kpp_3d_fields%L_OCEAN(j)) THEN
                  IF (vec(j) .lt. SCLR_range(j,i,1))
     +                 SCLR_range(j,i,1)=vec(j)
                  IF (vec(j) .gt. SCLR_range(j,i,2))
     +                 SCLR_range(j,i,2)=vec(j)
               ENDIF
            ENDDO
         ENDIF
         i=i+1
      ENDDO
!      deallocate(field)
!      deallocate(vec)

      RETURN
      END

      SUBROUTINE output_close(ncid)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include <netcdf.inc>
#include "kpp_3d_type.com"
#include "output.com"

      INTEGER status
      INTEGER, intent(in) :: ncid

      status=NF_CLOSE(ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      RETURN
      END

      SUBROUTINE output_open(file,ncid)

      IMPLICIT NONE
      INTEGER nuout,nuerr
      PARAMETER (nuout=6,nuerr=0)

#include <netcdf.inc>
#include "kpp_3d_type.com"
#include "output.com"

      INTEGER status
      INTEGER,intent(out) :: ncid
      CHARACTER(LEN=*),intent(in) :: file

      status=NF_OPEN(file,NF_WRITE,ncid)
      IF (status .NE. NF_NOERR) CALL HANDLE_ERR(status)

      RETURN
      END
