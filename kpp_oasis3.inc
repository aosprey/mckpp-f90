      INTEGER jpfldin,jpfldout            ! Number of fields sent and received, 
                                          ! respectively, by this model
#ifdef UM78
      PARAMETER(jpfldin=11,jpfldout=6)
#endif
#ifdef UM85 
      PARAMETER(jpfldin=9,jpfldout=6)
#endif
      INTEGER il_var_id_out(jpfldout)     ! IDs for the outgoing fields
      INTEGER il_var_id_in(jpfldin)       ! IDs for the incoming fields
      CHARACTER*8 cl_writ(jpfldout)       ! Names of the fields sent by this model
      CHARACTER*8 cl_read(jpfldin)        ! Names of the fields received by this model
      INTEGER il_mparout
      INTEGER il_comp_id

      COMMON /kpp_oasis3/ il_var_id_out,il_var_id_in,cl_writ,cl_read,
     +     il_mparout,il_comp_id
