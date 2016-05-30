
# process Substrate information using SPDE /RINLA .. no GMT dependency
  
  ## NOTE:: substrate size is really only relevant for SSE/snowcrab domain right now as no 
  ##        other data source has been found/identified
  ##        but working at the size of canada.east.highres for compatibility with bathymetry 
  ##        .. might change this in future as it is also expensive in time .. but really only done once in a while, sooo...
  ## TODO:: add data collected by snow crab survey and any others for that matter 
   
  p = list( project.name = "substrate" )
  p$project.root = project.datadirectory( p$project.name )
         
  p$init.files = loadfunctions( c( "spacetime", "utility", "parallel", "bathymetry", "polygons", "substrate", "coastline" ) )
  p$libs = RLibrary( "rgdal", "maps", "mapdata", "maptools", "lattice", "parallel", "INLA", "geosphere", "geoR", "gstat", "spBayes",  
                     "sp", "raster", "colorspace" ,  "splancs", "fields",
                     "bigmemory.sri", "synchronicity", "bigmemory", "biganalytics", "bigtabulate", "bigalgebra" )
  p = spatial.parameters( type="canada.east.highres", p=p ) # highest resolution still 
  
  p = spacetime.parameters(p)  # load defaults
    
  p$substrate.bigmemory.reset = FALSE  
 
  # cluster definition
  nc = 1
  # nc = 5
  p$clusters = rep( "localhost", nc )
  # p$clusters = c( rep( "nyx", nc ), rep ("tartarus", nc), rep("kaos", nc ) )


 
  ### -----------------------------------------------------------------
  make.substrate.db = FALSE
  if (make.substrate.db) {
    substrate.db ( DS="substrate.initial.redo" ) # bring in Kostelev's data ... stored as a SpatialGridDataFrame
		substrate.db ( DS="lonlat.highres.redo" ) # in future .. additional data would be added here ...
  }
 
  ### -----------------------------------------------------------------
  spatial.covariance.redo = FALSE
  if (spatial covariance.redo) {
    p$clusters = c( rep( "nyx", 24 ), rep ("tartarus", 24), rep("kaos", 24 ) )
    # p$substrate.bigmemory.reset = TRUE  # reset needed if variables entering are changing (eg., addiing covariates with interpolation, etc)
    substrate.db( p=p, DS="covariance.spatial.redo" ) 
  }
  covSp = substrate.db( p=p, DS="covariance.spatial" ) 


  ### -----------------------------------------------------------------
  spatial.interpolation.redo = FALSE
  if (spatial.interpolation.redo) {
    # do not use all CPU's as INLA itself is partially run in parallel
    # RAM reqiurements are a function of data density and mesh density .. currently ~ 12 GB / run
    p$clusters = c( rep( "nyx", 5 ), rep ("tartarus", 5), rep("kaos", 5 ) )
    # p$substrate.bigmemory.reset = TRUE   # reset needed if variables entering are changing (eg., addiing covariates with interpolation, etc)
    substrate.db( p=p, DS="spde.redo" ) 
    # to see the raw saved versions of the the results:
    # predSp = spacetime.db( p=p, DS="predictions.redo" )  
    # statSp = spacetime.db( p=p, DS="statistics.redo" )
    # to see the assimilated data:
    # B = substrate.db( p=p, DS="substrate.spacetime.finalize" )
  }


  ### -----------------------------------------------------------------
  # as the interpolation process is so expensive, regrid based off the above run
  # if you want more, will need to add to the list and modify the selection criteria
  p$grids.new = c( "canada.east.highres", "canada.east", "SSE", "snowcrab", "SSE.mpa" )
  substrate.db( p=p, DS="complete.redo" ) 

  # test outputs/ access methods
  # plot( substrate.db( p, DS="complete", return.format="brick" )$substrate ) # raster brick
  # spplot( substrate.db( p, DS="complete", return.format="sp" ), "substrate" ) # spatial points/grid data frame
   



