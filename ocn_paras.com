      REAL rhob,epsw  
      LOGICAL L_SSRef
      common/ ocn_paras / rhob,epsw,L_SSref

