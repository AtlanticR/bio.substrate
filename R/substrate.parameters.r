
substrate.parameters = function(DS="bio.substrate", p=NULL, resolution="canada.east.highres", nc=1) {

  if ( is.null(p) ) p=list()
  if ( !exists("project.name", p) ) p$project.name=DS


  if (DS=="bio.substrate"){
   
    p$project.root = project.datadirectory( p$project.name )

    p$libs = bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons",
        "bio.substrate", "bio.coastline" )
    p$libs = c( p$libs, RLibrary( "rgdal", "maps", "mapdata", "maptools", "lattice", "parallel", "INLA",
                       "geosphere", "geoR", "gstat", "spBayes",
                       "sp", "raster", "colorspace" ,  "splancs", "fields", "ff", "ffbase" ) )

    p = spacetime_parameters( type=resolution, p=p ) # highest resolution still
    p = spacetime.parameters(p)  # load defaults

    # cluster definition
    p$clusters = rep( "localhost", nc )

    return(p)
  }


  if (DS=="bio.substrate.local") {
    p$spacetime_variogram_engine = "gstat"  # "geoR" seg faults frequently ..
    p$spacetime.prediction.dist.min = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    p$upsampling = c( 1.1, 1.2, 1.5, 2 )  # local block search fractions
    p$downsampling = c( 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2 ) # local block search fractions  -- need to adjust based upon data density
    p$mesh.boundary.resolution = 120 # discretization
    p$mesh.boundary.convex = -0.03  # curavature of boundary

    # p$variables = list( Y="substrate", X=c("z", "dZ", "ddZ", "Z.rangeMode" ), LOCS=c("plon", "plat") )
    p$variables = list( Y="substrate", COV=c("z", "dZ" ), LOCS=c("plon", "plat") )
    p$spacetime.link = function( X ) { log(X ) }  ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one -order of magnitude
    p$spacetime.invlink = function( X ) { exp(X) }
    p$dist.max = 100 # length scale (km) of local analysis .. for acceptance into the local analysis/model
    p$dist.min = 75 # lower than this .. subsampling occurs
    p$dist.pred = 0.95 # % of dist.max where **predictions** are retained (to remove edge effects)
    p$n.min = 30 # n.min/n.max changes with resolution: at p$pres=0.25, p$dist.max=25: the max count expected is 40000
    p$n.max = 8000 # numerical time/memory constraint -- anything larger takes too much time
    p$expected.range = 50 #+units=km km , with dependent var on log scale
    p$expected.sigma = 1e-1  # spatial standard deviation (partial sill) .. on log scale
    p$spatial.field.name = "spatial.field"  # name used in formula to index the spatal random field
    p$spacetime_engine_modelformula = formula( substrate ~ -1 + intercept
      + f( inla.group(log(z+1000) ), model="rw2")
      + f( inla.group(log(dZ+0.01)), model="rw2")
      # + f( inla.group( log(ddZ+0.01) ), model="rw2")
      # + f( inla.group( log(Z.rangeMode+0.01)), model="rw2" )
      + f( spatial.field, model=SPDE ) )
    p$spacetime.family = "gaussian"
    p$spacetime.outputs = c( "predictions.projected", "statistics" ) # "random.field", etc.
    p$statsvars = c("range", "range.sd", "spatial.error", "observation.error")
    # if not in one go, then the value must be reconstructed from the correct elements:
    p$sbbox = spacetime_db( p=p, DS="statistics.box" ) # bounding box and resoltuoin of output statistics defaults to 1 km X 1 km
    p$spacetime.stats.boundary.redo = FALSE ## estimate boundart of data to speed up stats collection? Do not need to redo if bounds have already been determined
    p$nPreds = p$nplons * p$nplats
 
    p$spacetime.posterior.extract = function(s, rnm) {
      # rnm are the rownames that will contain info about the indices ..
      # optimally the grep search should only be done once but doing so would
      # make it difficult to implement in a simple structure/manner ...
      # the overhead is minimal relative to the speed of modelling and posterior sampling
      i_intercept = grep("intercept", rnm, fixed=TRUE ) # matching the model index "intercept" above .. etc
      i_spatial.field = grep("spatial.field", rnm, fixed=TRUE )
      return(  s$latent[i_intercept,1] + s$latent[ i_spatial.field,1] )
    }

    return(p)
  }

}



