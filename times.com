      REAL dtsec,time,startt,finalt
      INTEGER ntime,nstart,nend
      common/ times     / dtsec,time,startt,finalt,ntime,nstart,nend

