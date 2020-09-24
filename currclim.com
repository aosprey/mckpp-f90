	LOGICAL L_CLIMCURR, L_UPD_CLIMCURR, L_DAMP_CURR
	INTEGER ndtupdcurr, dtuvdamp
	common /curr_clim/ L_CLIMCURR, L_UPD_CLIMCURR, ndtupdcurr,L_DAMP_CURR, dtuvdamp
