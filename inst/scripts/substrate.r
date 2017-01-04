
  ## NOTE:: substrate size is really only relevant for SSE/snowcrab domain right now as no
  ##        other data source has been found/identified
  ##        but working at the size of canada.east.highres for compatibility with bathymetry
  ## TODO:: add data collected by snow crab survey and any others for that matter

  p = bio.substrate::substrate.parameters()

  if ( basedata.redo ) {
    substrate.db ( DS="substrate.initial.redo" ) # bring in Kostelev's data ... stored as a SpatialGridDataFrame
		substrate.db ( DS="lonlat.highres.redo" ) # in future .. additional data would be added here
  }

  p = bio.substrate::substrate.parameters() # reset to defaults
  p$lbm_local_modelengine = "krige" 
  p$storage.backend="bigmemory.ram"  # filebacked metods are still too slow ..
  p = bio.substrate::substrate.parameters( p=p, DS="lbm" )
  # p$clusters = rep("localhost",  detectCores() )
  p = lbm( p=p, DATA='substrate.db( p=p, DS="lbm.inputs" )' )
   
  substrate.db( p=p, DS="lbm.finalize.redo" )
  # B = substrate.db( p=p, DS="lbm.finalize" )

 
  # as the interpolation process is so expensive, regrid based off the above run
  # if you want more, will need to add to the list and modify the selection criteria
  # this requires "raster" (it is possible to use fields and be a bit faster but this is simpler for now)
  p$new.grids = c( "canada.east.highres", "canada.east", "SSE", "snowcrab", "SSE.mpa" )
  substrate.db( p=p, DS="complete.redo" )

  # test outputs/ access methods
  # plot( substrate.db( p, DS="complete", return.format="brick" )$substrate ) # raster brick
  # spplot( substrate.db( p, DS="complete", return.format="sp" ), "substrate" ) # spatial points/grid data frame




