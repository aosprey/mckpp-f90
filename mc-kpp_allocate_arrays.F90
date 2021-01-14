SUBROUTINE mckpp_allocate_3d_fields(kpp_3d_fields)

  USE mckpp_data_types
  USE mckpp_parameters 

  IMPLICIT NONE 

  TYPE(kpp_3d_type) :: kpp_3d_fields 

  ALLOCATE( kpp_3d_fields%u(npts,nzp1,nvel) ) 
  ALLOCATE( kpp_3d_fields%x(npts,nzp1,nsclr) ) 
  ALLOCATE( kpp_3d_fields%rig(npts,nzp1) )
  ALLOCATE( kpp_3d_fields%dbloc(npts,nz) )
  ALLOCATE( kpp_3d_fields%shsq(npts,nzp1) )
  ALLOCATE( kpp_3d_fields%hmixd(npts,0:1) )
  ALLOCATE( kpp_3d_fields%us(npts,nzp1,nvel,0:1) )
  ALLOCATE( kpp_3d_fields%xs(npts,nzp1,nsclr,0:1) )
  ALLOCATE( kpp_3d_fields%rho(npts,0:nzp1tmax) )
  ALLOCATE( kpp_3d_fields%cp(npts,0:nzp1tmax) )
  ALLOCATE( kpp_3d_fields%buoy(npts,nzp1tmax) )
  ALLOCATE( kpp_3d_fields%rhoh2o(npts) )
  ALLOCATE( kpp_3d_fields%ocdepth(npts) )
  ALLOCATE( kpp_3d_fields%f(npts) )
  ALLOCATE( kpp_3d_fields%swfrac(npts,nzp1) )
  ALLOCATE( kpp_3d_fields%swdk_opt(npts,0:nz) )
  ALLOCATE( kpp_3d_fields%difm(npts,0:nztmax) )
  ALLOCATE( kpp_3d_fields%difs(npts,0:nztmax) )
  ALLOCATE( kpp_3d_fields%dift(npts,0:nztmax) )
  ALLOCATE( kpp_3d_fields%wu(npts,0:nztmax,nvp1) )
  ALLOCATE( kpp_3d_fields%wx(npts,0:nztmax,nsp1) )
  ALLOCATE( kpp_3d_fields%wxnt(npts,0:nztmax,nsclr) )
  ALLOCATE( kpp_3d_fields%ghat(npts,nztmax) )
  ALLOCATE( kpp_3d_fields%relax_sst(npts) ) 
  ALLOCATE( kpp_3d_fields%fcorr(npts) )
  ALLOCATE( kpp_3d_fields%cplwght(npts_globe) )
  ALLOCATE( kpp_3d_fields%sst0(npts) )
  ALLOCATE( kpp_3d_fields%fcorr_twod(npts) )
  ALLOCATE( kpp_3d_fields%sfcorr_twod(npts) )
  ALLOCATE( kpp_3d_fields%tinc_fcorr(npts,nzp1) )
  ALLOCATE( kpp_3d_fields%sinc_fcorr(npts,nzp1) )
  ALLOCATE( kpp_3d_fields%fcorr_withz(npts,nzp1) ) 
  ALLOCATE( kpp_3d_fields%sfcorr(npts) )
  ALLOCATE( kpp_3d_fields%sfcorr_withz(npts,nzp1) )
  ALLOCATE( kpp_3d_fields%advection(npts,maxmodeadv,2) )
  ALLOCATE( kpp_3d_fields%relax_sal(npts) ) 
  ALLOCATE( kpp_3d_fields%scorr(npts,nzp1) )
  ALLOCATE( kpp_3d_fields%relax_ocnt(npts) ) 
  ALLOCATE( kpp_3d_fields%ocntcorr(npts,nzp1) )
  ALLOCATE( kpp_3d_fields%sal_clim(npts,nzp1) ) 
  ALLOCATE( kpp_3d_fields%ocnt_clim(npts,nzp1) )
  ALLOCATE( kpp_3d_fields%hmix(npts) )
  ALLOCATE( kpp_3d_fields%kmix(npts) )
  ALLOCATE( kpp_3d_fields%tref(npts) )
  ALLOCATE( kpp_3d_fields%uref(npts) )
  ALLOCATE( kpp_3d_fields%vref(npts) )
  ALLOCATE( kpp_3d_fields%ssurf(npts) )
  ALLOCATE( kpp_3d_fields%sref(npts) )
  ALLOCATE( kpp_3d_fields%ssref(npts) )
  ALLOCATE( kpp_3d_fields%sflux(npts,nsflxs,5,0:njdt) )
  ALLOCATE( kpp_3d_fields%dlat(npts) )
  ALLOCATE( kpp_3d_fields%dlon(npts) )
  ALLOCATE( kpp_3d_fields%freeze_flag(npts) )
  ALLOCATE( kpp_3d_fields%reset_flag(npts) )
  ALLOCATE( kpp_3d_fields%dampu_flag(npts) )
  ALLOCATE( kpp_3d_fields%dampv_flag(npts) )
  ALLOCATE( kpp_3d_fields%u_init(npts,nzp1,nvel) )
  ALLOCATE( kpp_3d_fields%bottom_temp(npts) )
  ALLOCATE( kpp_3d_fields%taux(npts) )
  ALLOCATE( kpp_3d_fields%tauy(npts) )
  ALLOCATE( kpp_3d_fields%swf(npts) )
  ALLOCATE( kpp_3d_fields%lwf(npts) )
  ALLOCATE( kpp_3d_fields%lhf(npts) ) 
  ALLOCATE( kpp_3d_fields%shf(npts) )
  ALLOCATE( kpp_3d_fields%rain(npts) )
  ALLOCATE( kpp_3d_fields%snow(npts) )
#ifdef mckpp_couple
  ALLOCATE( kpp_3d_fields%sst(nx_globe,ny_globe) )
  ALLOCATE( kpp_3d_fields%iceconc(nx_globe,ny_globe) )
  ALLOCATE( kpp_3d_fields%usf(nx_globe,ny_globe) )
  ALLOCATE( kpp_3d_fields%vsf(nx_globe,ny_globe) )
  ALLOCATE( kpp_3d_fields%icedepth(nx_globe,ny_globe) )
  ALLOCATE( kpp_3d_fields%snowdepth(nx_globe,ny_globe) )
#else
  ALLOCATE( kpp_3d_fields%sst(nx,ny) )
  ALLOCATE( kpp_3d_fields%iceconc(nx,ny) )
  ALLOCATE( kpp_3d_fields%usf(nx,ny) )
  ALLOCATE( kpp_3d_fields%vsf(nx,ny) )
  ALLOCATE( kpp_3d_fields%icedepth(nx,ny) )
  ALLOCATE( kpp_3d_fields%snowdepth(nx,ny) )
#endif
  ALLOCATE( kpp_3d_fields%l_ocean(npts) ) 
  ALLOCATE( kpp_3d_fields%l_initflag(npts) )
  ALLOCATE( kpp_3d_fields%old(npts) ) 
  ALLOCATE( kpp_3d_fields%new(npts) )
  ALLOCATE( kpp_3d_fields%jerlov(npts) )
  ALLOCATE( kpp_3d_fields%nmodeadv(npts,2) )
  ALLOCATE( kpp_3d_fields%modeadv(npts,maxmodeadv,2) )

END SUBROUTINE mckpp_allocate_3d_fields


SUBROUTINE mckpp_allocate_1d_fields(kpp_1d_fields) 

  USE mckpp_data_types
  USE mckpp_parameters 

  IMPLICIT NONE 

  TYPE(kpp_1d_type) :: kpp_1d_fields 

  ALLOCATE( kpp_1d_fields%u(nzp1,nvel) )
  ALLOCATE( kpp_1d_fields%u_init(nzp1,nvel) )
  ALLOCATE( kpp_1d_fields%x(nzp1,nsclr) )
  ALLOCATE( kpp_1d_fields%rig(nzp1) )
  ALLOCATE( kpp_1d_fields%dbloc(nz) )
  ALLOCATE( kpp_1d_fields%shsq(nzp1) )
  ALLOCATE( kpp_1d_fields%hmixd(0:1) )
  ALLOCATE( kpp_1d_fields%us(nzp1,nvel,0:1) )
  ALLOCATE( kpp_1d_fields%xs(nzp1,nsclr,0:1) )
  ALLOCATE( kpp_1d_fields%rho(0:nzp1tmax) )
  ALLOCATE( kpp_1d_fields%cp(0:nzp1tmax) )
  ALLOCATE( kpp_1d_fields%buoy(nzp1tmax) )
  ALLOCATE( kpp_1d_fields%swfrac(nzp1) )
  ALLOCATE( kpp_1d_fields%swdk_opt(0:nz) )
  ALLOCATE( kpp_1d_fields%difm(0:nztmax) )
  ALLOCATE( kpp_1d_fields%difs(0:nztmax) )
  ALLOCATE( kpp_1d_fields%dift(0:nztmax) )
  ALLOCATE( kpp_1d_fields%wu(0:nztmax,nvp1) )
  ALLOCATE( kpp_1d_fields%wx(0:nztmax,nsp1) )
  ALLOCATE( kpp_1d_fields%wxnt(0:nztmax,nsclr) )
  ALLOCATE( kpp_1d_fields%ghat(nztmax) )
  ALLOCATE( kpp_1d_fields%tinc_fcorr(nzp1) )
  ALLOCATE( kpp_1d_fields%sinc_fcorr(nzp1) )
  ALLOCATE( kpp_1d_fields%fcorr_withz(nzp1) )
  ALLOCATE( kpp_1d_fields%sfcorr_withz(nzp1) )
  ALLOCATE( kpp_1d_fields%advection(maxmodeadv,2) )
  ALLOCATE( kpp_1d_fields%scorr(nzp1) ) 
  ALLOCATE( kpp_1d_fields%ocntcorr(nzp1) )
  ALLOCATE( kpp_1d_fields%sal_clim(nzp1) ) 
  ALLOCATE( kpp_1d_fields%ocnt_clim(nzp1) )
  ALLOCATE( kpp_1d_fields%sflux(nsflxs,5,0:njdt) )
  ALLOCATE( kpp_1d_fields%talpha(0:nzp1tmax) )
  ALLOCATE( kpp_1d_fields%sbeta(0:nzp1tmax) ) 
  ALLOCATE( kpp_1d_fields%modeadv(maxmodeadv,2) )

END SUBROUTINE mckpp_allocate_1d_fields


SUBROUTINE mckpp_allocate_const_fields(kpp_const_fields) 


  USE mckpp_data_types
  USE mckpp_parameters 

  IMPLICIT NONE 

  TYPE(kpp_const_type) :: kpp_const_fields 

  ALLOCATE( kpp_const_fields%zm(nzp1) ) 
  ALLOCATE( kpp_const_fields%hm(nzp1) )
  ALLOCATE( kpp_const_fields%dm(0:nz) )
  ALLOCATE( kpp_const_fields%relax_sst_in(ny) )
  ALLOCATE( kpp_const_fields%relax_ocnt_in(ny) )  
  ALLOCATE( kpp_const_fields%relax_sal_in(ny) ) 

END SUBROUTINE mckpp_allocate_const_fields 

