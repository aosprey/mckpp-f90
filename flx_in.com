      LOGICAL L_FLUXDATA,L_UPD_FLUXDATA,L_REST

      INTEGER ncid_flx
      INTEGER timein_id,varin_id(8)
      REAL*4 first_timein
      common /flx_in/ L_FLUXDATA,ncid_flx,timein_id,varin_id,
     &     first_timein,L_REST,L_UPD_FLUXDATA

