
  ## NOTE:: substrate size is really only relevant for SSE/snowcrab domain right now as no
  ##        other data source has been found/identified
  ##        but working at the size of canada.east.highres for compatibility with bathymetry
  ## TODO:: add data collected by snow crab survey and any others for that matter

  p = bio.substrate::substrate.parameters( DS="bio.substrate" )

  if ( basedata.redo ) {
    substrate.db ( DS="substrate.initial.redo" ) # bring in Kostelev's data ... stored as a SpatialGridDataFrame
		substrate.db ( DS="lonlat.highres.redo" ) # in future .. additional data would be added here
  }

  p = bio.substrate::substrate.parameters() # reset to defaults
  p$hivemod_local_modelengine = "kernel.density" 
  p$storage.backend="bigmemory.ram"  # filebacked metods are still too slow ..
  p = bio.substrate::substrate.parameters( p=p, DS="hivemod" )

  p$clusters = rep("localhost",  detectCores() )
  DATA = 'substrate.db( p=p, DS="substrate.hivemod" )'
  p = hivemod( p=p, DATA=DATA )
   
  substrate.db( p=p, DS="substrate.hivemod.finalize.redo" )

  B = substrate.db( p=p, DS="substrate.hivemod.finalize" )

  ### -----------------------------------------------------------------
  # as the interpolation process is so expensive, regrid based off the above run
  # if you want more, will need to add to the list and modify the selection criteria
  p$grids.new = c( "canada.east.highres", "canada.east", "SSE", "snowcrab", "SSE.mpa" )
  substrate.db( p=p, DS="complete.redo" )

  # test outputs/ access methods
  # plot( substrate.db( p, DS="complete", return.format="brick" )$substrate ) # raster brick
  # spplot( substrate.db( p, DS="complete", return.format="sp" ), "substrate" ) # spatial points/grid data frame




