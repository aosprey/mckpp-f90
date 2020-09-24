      LOGICAL LKPP , LRI, LDD, LICE, LBIO, LNBFLX, LTGRID, LRHS,&
           lsaveaverages,L_STRETCHGRID,lrepeatdat,&
           lradtq,lfluxSSTdat,&
           lLWupSSTdat,lrhdat,lclddat 
      
      common /proc_swit/ LKPP , LRI, LDD, LICE, LBIO, LNBFLX, LTGRID,LRHS,&
           lsaveaverages,L_STRETCHGRID,lrepeatdat,lradtq,lfluxSSTdat,&
           lLWupSSTdat,lrhdat,lclddat
     
