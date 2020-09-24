       LOGICAL L_VARY_BOTTOM_TEMP, L_UPD_BOTTOM_TEMP, L_PERIODIC_BOTTOM_TEMP
       INTEGER ndtupdbottom, bottom_temp_period
       CHARACTER*40 bottomin_file
       COMMON /bottomclim/ L_VARY_BOTTOM_TEMP,ndtupdbottom,&
            bottomin_file,L_UPD_BOTTOM_TEMP,L_PERIODIC_BOTTOM_TEMP,bottom_temp_period
