      LOGICAL L_FLUXDATA,L_REST

      INTEGER ncid_flx
      INTEGER timein_id,varin_id(7)
      REAL*4 first_timein
      CHARACTER*50 :: forcing_file          
      common /flx_in/ L_FLUXDATA,ncid_flx,timein_id,varin_id,first_timein,L_REST,&
	      forcing_file

