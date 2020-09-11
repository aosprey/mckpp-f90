c      REAL dlat(npts),dlon(npts),rlat(npts),rlon(npts),f(npts)
      LOGICAL L_REGGRID
      common/ location  / L_REGGRID

