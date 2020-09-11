	LOGICAL L_DAMP_CURR
	INTEGER dtuvdamp
	common /curr_clim/ L_DAMP_CURR, dtuvdamp
