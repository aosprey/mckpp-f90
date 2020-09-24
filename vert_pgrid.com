      REAL DMAX
      LOGICAL L_VGRID_FILE
      CHARACTER(LEN=50) vgrid_file
      common /vertpgrid/ DMAX,L_VGRID_FILE,vgrid_file

