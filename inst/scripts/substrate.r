
# process Substrate information using SPDE /RINLA .. no GMT dependency

  ## NOTE:: substrate size is really only relevant for SSE/snowcrab domain right now as no
  ##        other data source has been found/identified
  ##        but working at the size of canada.east.highres for compatibility with bathymetry
  ##        .. might change this in future as it is also expensive in time .. but really only done once in a while, sooo...
  ## TODO:: add data collected by snow crab survey and any others for that matter

  p = bio.substrate::substrate.parameters( DS="bio.substrate", nc=1 )

  # nc=5; p$clusters = c( rep( "nyx", nc ), rep ("tartarus", nc), rep("kaos", nc ) )


  ### -----------------------------------------------------------------

    substrate.db ( DS="substrate.initial.redo" ) # bring in Kostelev's data ... stored as a SpatialGridDataFrame
		substrate.db ( DS="lonlat.highres.redo" ) # in future .. additional data would be added here ...

    substrate.db( p=p, DS="substrate.spacetime.inputs.data.redo" )  # Warning: req ~ 15 min, 40 GB RAM (2015, Jae) data to model (with covariates if any)
    substrate.db( p=p, DS="substrate.spacetime.inputs.prediction.redo" ) # i.e, pred locations (with covariates if any )


  ### -----------------------------------------------------------------

  p$clusters = c( rep( "nyx", 24 ), rep ("tartarus", 24), rep("kaos", 24 ) )
  p = spacetime( method="spatial.covariance", p=p, overwrite=TRUE ,
    DATA=list ( input=substrate.db( p=p, DS="substrate.spacetime.inputs.data" ),
                output=substrate.db( p=p, DS="substrate.spacetime.inputs.prediction")) )

      # to see the raw saved versions of the the results:
      covSp = spacetime( p=p, DS="spatial.covariance" ) # load saved data


  ### -----------------------------------------------------------------
  # do not use all CPU's as INLA itself is partially run in parallel
  # RAM reqiurements are a function of data density and mesh density .. currently ~ 12 GB / run
  p$clusters = c( rep( "nyx", 5 ), rep ("tartarus", 5), rep("kaos", 5 ) )

  p = spacetime( method="xy.inla",
    DATA=list( input=substrate.db( p=p, DS="substrate.spacetime.inputs.data" ),
               output=substrate.db( p=p, DS="substrate.spacetime.inputs.prediction") ),
    p=p, overwrite=TRUE )

      # to see the raw saved versions of the the results:
      predSp = spacetime( p=p, DS="inla.predictions" )
      statSp = spacetime( p=p, DS="inla.statistics" )

  B = substrate.db( p=p, DS="substrate.spacetime.finalize" )

  ### -----------------------------------------------------------------
  # as the interpolation process is so expensive, regrid based off the above run
  # if you want more, will need to add to the list and modify the selection criteria
  p$grids.new = c( "canada.east.highres", "canada.east", "SSE", "snowcrab", "SSE.mpa" )
  substrate.db( p=p, DS="complete.redo" )

  # test outputs/ access methods
  # plot( substrate.db( p, DS="complete", return.format="brick" )$substrate ) # raster brick
  # spplot( substrate.db( p, DS="complete", return.format="sp" ), "substrate" ) # spatial points/grid data frame




