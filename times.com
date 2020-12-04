      REAL dtsec,time,startt,finalt
      INTEGER ntime,nstart,nend,nyear
      common/ times     / dtsec,time,startt,finalt,ntime,nstart,nend,nyear

