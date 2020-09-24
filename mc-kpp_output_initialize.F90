#ifdef MCKPP_CAM3
#include <misc.h>
#include <params.h>
SUBROUTINE mckpp_output_initialize(filename,output_type)
  USE mckpp_types, only : kpp_global_fields,kpp_const_fields
#else
SUBROUTINE mckpp_output_initialize(filename,kpp_3d_fields,kpp_const_fields,output_type)
#endif
  
  IMPLICIT NONE
  INTEGER,parameter :: nuout=6,nuerr=0
#include <netcdf.inc>

#ifdef MCKPP_CAM3
#include <parameter.inc>
  INTEGER, parameter :: my_nx=NX_GLOBE, my_ny=NY_GLOBE
#else
  ! Automatically includes parameter.inc!
#include <mc-kpp_3d_type.com>
  TYPE(kpp_3d_type),intent(in) :: kpp_3d_fields
  TYPE(kpp_const_type) :: kpp_const_fields
  INTEGER, parameter :: my_nx=NX, my_ny=NY
#endif

  INTEGER i,j,k,m,ix,iy,ncid,status,extra_time
  CHARACTER*11 varname(N_VAROUTS),singname(N_SINGOUTS)
  CHARACTER*50 longname(N_VAROUTS),singlong(N_SINGOUTS),filename
  CHARACTER*15 units(N_VAROUTS),singunits(N_SINGOUTS)      
  CHARACTER*6 type
  CHARACTER*7 time_name
  CHARACTER*4 z_name,output_type
  REAL*4, allocatable :: time_out(:)
  INTEGER :: vec_varids(N_VAROUTS),sing_varids(N_SINGOUTS),&
       dt_vec(N_VAROUTS),dt_sing(N_SINGOUTS),n_tdims,n_zdims,n_ddims,&
       vert_dimid(N_VAROUTS),dims(4),zflag(N_VAROUTS),& ! 1 for z-levels, 0 for d-levels      
       time_dimids(N_VAROUTS+N_SINGOUTS),& ! Contains a value for each diagnostic
       time_varids(N_VAROUTS+N_SINGOUTS),& ! Contains only unique time variables
       dt_timeids(N_VAROUTS+N_SINGOUTS),&
       zprof_vec(N_VAROUTS),zprof_zids(N_VAROUTS),&
       zprof_dids(N_VAROUTS),z_dimids(N_VAROUTS),&
       z_varids(N_VAROUTS),d_varids(N_VAROUTS)
  REAL*4 ZOUT(NZP1),ALON(my_nx),ALAT(my_ny)
  REAL*4 delta
  LOGICAL tdim_found,zdim_found  
  INTEGER lon_id,lat_id,z_id,h_id,d_id,time_id,&
       londim,latdim,zdim,ddim,hdim,timdim
  
  DATA varname/'u','v','T','S','B',&
       'wu','wv','wT','wS','wB',&
       'wTnt',&
       'difm','dift','difs',&
       'rho','cp','scorr','Rig','dbloc','Shsq','tinc_fcorr',&
       'fcorr_z','sinc_fcorr'/
  DATA zflag/1,1,1,1,1,&
       0,0,0,0,0,&
       0,&
       0,0,0,&
       1,1,1,1,1,1,1,&
       1,1/
  DATA longname/&
       'Zonal velocity',&
       'Meridional velocity',&
       'Temperature',&
       'Salinity',&
       'Buoyancy',&
       'Turbulent Zonal Velocity Flux',&
       'Turbulent Meridional Velocity Flux',&
       'Turbulent Temperature Flux',&
       'Turbulent Salinity Flux',&
       'Turbulent Buoyancy Flux',&
       'Non-Turbulent Temperature Flux',&
       'Diffusion Coefficient (Momentum)',&
       'Diffusion Coefficient (Temperature)',&
       'Diffusion Coefficient (Salinity)',&
       'Density',&
       'Specific Heat Capacity',&
       'Salinity correction (with depth)',&
       'Local Richardson Number in kpp.f',&
       'Local delta buoyancy in kpp.f',&
       'Local shear-squared term in kpp.f',&
       'Temperature increment flux correction',&
       'Heat correction as flux (dT/dt*rho*cp)',&
       'Salinity increment flux correction'/
  DATA units/&
       'm/s',&
       'm/s',&
       'degC',&
       'o/oo',&
       'm/s^2',&
       'm^2/s^2',&
       'm^2/s^2',&
       'degC m/s',&
       'o/oo m/s',&
       'm^2/s^3',&
       'degC m/s',&
       'm^2/s',&
       'm^2/s',&
       'm^2/s',&
       'kg/m^3',&
       'J/kg/K',&
       'o/oo/s',&
       'unitless',&
       'm/s^2',&
       'm^2/s^2',&
       'K/timestep',&
       'W/m^3',&
       'o/oo/timestep'/
  DATA singname/&
       'hmix',&
       'fcorr',&
       'taux_in',&
       'tauy_in',&
       'solar_in',&
       'nsolar_in',&
       'PminusE_in',&
       'cplwght',& 
       'freeze_flag',&
       'comp_flag',&
       'dampu_flag',&
       'dampv_flag'/
#if defined(MCKPP_COUPLE) && defined(MCKPP_CAM3)
  DATA type /'CAM3  '/
#elif defined(MCKPP_COUPLE) && defined(OASIS2)
  DATA type /'OASIS2'/
#elif defined(MCKPP_COUPLE) && defined(OASIS3)
  DATA type /'OASIS3'/
#elif defined(MCKPP_COUPLE) && defined(GFS)
  DATA type /' GFS  '/
#else
  DATA type /'netCDF'/
#endif
  
  DATA singlong/&
       'Mixed Layer Depth',&
       'Flux Correction',&
       'Zonal wind stress from ',&
       'Meridional wind stress from ',&
       'Solar from ',&
       'Non-solar from ',&
       'P minus E from ',&
       'Coupling weight',&
       'Fraction of levels below freezing',&
       'Number of integrations (<0 = isothermal reset)',&
       'Fraction of levels with ui~u**2',&
       'Fraction of levels with vi~v**2'/
  
  DATA singunits/&
       'm',&
       'W/m^2',&
       'N/m^2',&
       'N/m^2',&
       'W/m^2',&
       'W/m^2',&
       'mm/s',&
       'none',&
       'fraction',&
       'unitless',&
       'fraction',&
       'fraction'/
  
  SELECT CASE(output_type)
  CASE('inst')        
     dt_vec=kpp_const_fields%ndt_varout_inst
     dt_sing=kpp_const_fields%ndt_singout_inst
     zprof_vec=kpp_const_fields%zprof_varout_inst
  CASE('mean')
     dt_vec=kpp_const_fields%ndt_varout_mean
     dt_sing=kpp_const_fields%ndt_singout_mean
     zprof_vec=kpp_const_fields%zprof_varout_mean
  CASE('minx')
     dt_vec=kpp_const_fields%ndt_varout_range
     dt_sing=kpp_const_fields%ndt_singout_range
     zprof_vec=kpp_const_fields%zprof_varout_range
  CASE('maxx')
     dt_vec=kpp_const_fields%ndt_varout_range
     dt_sing=kpp_const_fields%ndt_singout_range
     zprof_vec=kpp_const_fields%zprof_varout_range
  CASE DEFAULT
     WRITE(6,*) 'Unknown type of file specified in mc-kpp_init_output.f90'
     CALL MCKPP_ABORT 
  END SELECT

  singlong(3:7)=singlong(3:7)//type
  
  DO ix=1,my_nx
#ifdef MCKPP_CAM3
     alon(ix)=kpp_global_fields%longitude(ix)
#else
     alon(ix)=kpp_3d_fields%dlon(ix)
#endif
  ENDDO
  DO iy=1,my_ny
#ifdef MCKPP_CAM3
     alat(iy)=kpp_global_fields%latitude(iy)
#else
     alat(iy)=kpp_3d_fields%dlat((iy-1)*my_nx+1)
#endif
  ENDDO
  
  status=NF_CREATE(filename, nf_clobber, ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  WRITE(nuout,*) 'MCKPP_OUTPUT_INITIALIZE: Output file ',filename,' created successfully: ',ncid
  
  delta=0.0
  IF (my_nx .GT. 1 ) delta=alon(2)-alon(1)
  CALL MCKPP_NCDF_DEF_DIM (ncid,londim,my_nx,lon_id,'longitude','deg',delta,' ')
  delta=0.0
  IF (my_ny .GT. 1) delta=alat(2)-alat(1)
  CALL MCKPP_NCDF_DEF_DIM (ncid,latdim,my_ny,lat_id,'latitude','deg',delta,' ')
  CALL MCKPP_NCDF_DEF_DIM (ncid,hdim,NZP1,h_id,'h','m',0.0,'Layer Thickness')  
  delta=kpp_const_fields%dtout
  
  n_tdims=0
  n_ddims=0
  n_zdims=0
  IF (dt_vec(1).gt.0) THEN
     n_tdims=1
     dt_timeids(1)=dt_vec(1)
     CALL MCKPP_NCDF_DEF_DIM(ncid,time_dimids(1),NF_UNLIMITED,time_varids(1),'time_1','days',delta,' ')
     j=NZP1         
     IF (zprof_vec(1) .gt. 0) j=kpp_const_fields%zprofs_nvalid(zprof_vec(1))
     IF (zflag(1).eq.0) THEN
        n_ddims=1          
        CALL MCKPP_NCDF_DEF_DIM(ncid,z_dimids(1),j,d_varids(1),'d_1','m',0.0,'Depth of interfaces')
        zprof_dids(1)=zprof_vec(1)
     ELSEIF (zflag(1).eq.1) THEN
        n_zdims=1
        CALL MCKPP_NCDF_DEF_DIM(ncid,z_dimids(1),j,z_varids(1),'z_1','m',0.0,'')
        zprof_zids(1)=zprof_vec(1)
     ELSE
        WRITE(6,*) 'Incorrect value of zflag for '//longname(1)
        CALL MCKPP_ABORT               
     ENDIF
  ENDIF
  IF (output_type .EQ. 'inst' .AND. .NOT. kpp_const_fields%L_RESTART) THEN
     ! Add an extra point to time dimensions for initial condition
     extra_time=1
  ELSE
     extra_time=0
  ENDIF
  
  DO j=2,N_VAROUTS
     IF (dt_vec(j) .gt. 0) THEN
        tdim_found=.FALSE.
        zdim_found=.FALSE.
        DO k=1,j-1 
           IF (dt_vec(j) .eq. dt_vec(k)) THEN
              tdim_found=.TRUE.
              time_dimids(j)=time_dimids(k)
           ENDIF
           IF (zprof_vec(j).eq.zprof_vec(k) .and. zflag(j) .eq. zflag(k) &
                .and. dt_vec(k).gt.0) THEN
              zdim_found=.TRUE.
              z_dimids(j)=z_dimids(k)
           ENDIF
        ENDDO
        IF (.NOT.tdim_found) THEN 
           n_tdims=n_tdims+1
           dt_timeids(n_tdims)=dt_vec(j)
           WRITE(time_name,'(A5,I0)') 'time_',n_tdims          
           IF (n_tdims .eq. 1) THEN 
              CALL MCKPP_NCDF_DEF_DIM(ncid,time_dimids(j),NF_UNLIMITED,&
                   time_varids(1),'time_1','days',delta,' ')
           ELSE
              CALL MCKPP_NCDF_DEF_DIM (ncid,time_dimids(j),&
                   kpp_const_fields%ndt_per_file/dt_vec(j)+extra_time,&
                   time_varids(n_tdims),time_name,'days',delta,' ')
           ENDIF
        ENDIF
        IF (.NOT.zdim_found) THEN
           m=NZP1
           IF (zprof_vec(j).gt.0) m=kpp_const_fields%zprofs_nvalid(zprof_vec(j))           
           IF (zflag(k).eq.0) THEN
              n_ddims=n_ddims+1
              WRITE(z_name,'(A2,I0)') 'd_',n_ddims
              CALL MCKPP_NCDF_DEF_DIM(ncid,z_dimids(j),m,d_varids(n_ddims),z_name,'m',0.0,&
                   'Depth of interfaces')
              zprof_dids(n_ddims)=zprof_vec(j)
           ELSEIF (zflag(k).eq.1) THEN
              n_zdims=n_zdims+1
              WRITE(z_name,'(A2,I0)') 'z_',n_zdims
              CALL MCKPP_NCDF_DEF_DIM(ncid,z_dimids(j),m,z_varids(n_zdims),z_name,'m',0.0,' ')
              zprof_zids(n_zdims)=zprof_vec(j)
           ELSE
              WRITE(6,*)'Incorrect value of zflag for '//longname(j)
              CALL MCKPP_ABORT
           ENDIF
        ENDIF
     ELSE
        time_dimids(j)=0
        z_dimids(j)=0
     ENDIF
  ENDDO

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
           IF (n_tdims .eq. 1) THEN              
              CALL MCKPP_NCDF_DEF_DIM(ncid,time_dimids(j+N_VAROUTS),&
                   NF_UNLIMITED,time_varids(1),'time_1','days',delta,' ')
           ELSE
              CALL MCKPP_NCDF_DEF_DIM(ncid,time_dimids(j+N_VAROUTS),&
                   kpp_const_fields%ndt_per_file/dt_sing(j)+extra_time,&
                   time_varids(n_tdims),time_name,'days',delta,' ')
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
     IF (dt_vec(k) .gt. 0) THEN
        CALL MCKPP_NCDF_DEF_VAR(ncid,vec_varids(k),4,dims,varname(k),units(k),longname(k))
        WRITE(6,*) 'Defined variable ',varname(k),' with dims ',dims
     ENDIF
  ENDDO
  
  DO k=1,N_SINGOUTS
     dims(3)=time_dimids(N_VAROUTS+k)
     IF (dt_sing(k) .gt. 0) &
          CALL MCKPP_NCDF_DEF_VAR(ncid,sing_varids(k),3,dims,singname(k),singunits(k),singlong(k))
  ENDDO
  
  status=NF_ENDDEF(ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  status=NF_PUT_VAR_REAL(ncid,lon_id,alon)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  status=NF_PUT_VAR_REAL(ncid,lat_id,alat)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  DO j=1,n_zdims
     m=1
     DO k=1,NZP1
        IF (kpp_const_fields%zprofs_mask(k,zprof_zids(j))) THEN
           zout(m)=kpp_const_fields%zm(k)
           m=m+1
        ENDIF
     ENDDO
     status=NF_PUT_VAR_REAL(ncid,z_varids(j),ZOUT)
     IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ENDDO
  DO j=1,n_ddims
     m=1
     DO k=1,NZP1
        IF (kpp_const_fields%zprofs_mask(k,zprof_dids(j))) THEN
           zout(m)=kpp_const_fields%dm(k-1)
           m=m+1
        ENDIF
     ENDDO
     status=NF_PUT_VAR_REAL(ncid,d_varids(j),ZOUT)
     IF (status.NE.NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  ENDDO
  
  DO k=1,NZP1
     ZOUT(k)=kpp_const_fields%hm(k)
  ENDDO
  status=NF_PUT_VAR_REAL(ncid,h_id,ZOUT)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  DO k=1,n_tdims         
     allocate(time_out(kpp_const_fields%ndt_per_file/dt_timeids(k)+extra_time))
     DO i=1-extra_time,kpp_const_fields%ndt_per_file/dt_timeids(k)
        time_out(i+extra_time)=i*kpp_const_fields%dto/&
             kpp_const_fields%spd*dt_timeids(k)+kpp_const_fields%time
        IF (output_type .eq. 'mean') &
             time_out(i+extra_time)=time_out(i+extra_time)-&
             0.5*kpp_const_fields%dto/kpp_const_fields%spd*dt_timeids(k)
     ENDDO
     status=NF_PUT_VARA_REAL(ncid,time_varids(k),&
          (/1/),(/kpp_const_fields%ndt_per_file/dt_timeids(k)+extra_time/),time_out)
     IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
     deallocate(time_out)
  ENDDO
  
  !write(nuout,*) 'NCDF output file initialized with name ',filename
  
  status=NF_CLOSE(ncid)
  IF (status .NE. NF_NOERR) CALL MCKPP_HANDLE_ERR(status)
  
  SELECT CASE(output_type)
  CASE ('inst')
     kpp_const_fields%ncid_out=ncid
     kpp_const_fields%varid_vec_inst=vec_varids
     kpp_const_fields%varid_sing_inst=sing_varids
  CASE ('mean')            
     kpp_const_fields%mean_ncid_out=ncid
     kpp_const_fields%varid_vec_mean=vec_varids
     kpp_const_fields%varid_sing_mean=sing_varids
  CASE ('minx')
     kpp_const_fields%min_ncid_out=ncid
     kpp_const_fields%varid_vec_range=vec_varids
     kpp_const_fields%varid_sing_range=sing_varids
  CASE ('maxx')
     kpp_const_fields%max_ncid_out=ncid
     kpp_const_fields%varid_vec_range=vec_varids
     kpp_const_fields%varid_sing_range=sing_varids
  END SELECT
  
  RETURN
END SUBROUTINE mckpp_output_initialize

