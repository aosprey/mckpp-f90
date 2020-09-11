      REAL DMAX, slab_depth
      LOGICAL L_VGRID_FILE, L_SLAB, L_COLUMBIA_LAND
      CHARACTER(LEN=200) vgrid_file
      common/ vert pgrid/ DMAX,
     +     L_VGRID_FILE,vgrid_file,slab_depth,L_SLAB,
     +	   L_COLUMBIA_LAND

