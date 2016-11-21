
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
    p = spacetime_parameters(p)  # load defaults

    # cluster definition
    p$clusters = rep( "localhost", nc )

    return(p)
  }


  if (DS=="bio.substrate.local") {
    p$spacetime_rsquared_threshold = 0.3 # lower threshold
    p$spacetime_distance_prediction = 5 # this is a half window km
    p$spacetime_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )

    p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )  # fractions of median distance scale (dist.max, dist.min)/2 to try in local block search

    p$spacetime_engine = "gam" # see model form in spacetime.r (method="xyts")
    p$spacetime_engine_modelformula = formula(
      substrate ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=40, bs="ts") +
      s(z, k=3, bs="ts") + s(dZ, k=3, bs="ts" ) + s(ddZ, k=3, bs="ts" ) )
    # p$spacetime_engine_modelformula = formula( substrate ~ -1 + intercept
    #   + f( inla.group(log(z+1000) ), model="rw2")
    #   + f( inla.group(log(dZ+0.01)), model="rw2")
    #   # + f( inla.group( log(ddZ+0.01) ), model="rw2")
    #   # + f( inla.group( log(Z.rangeMode+0.01)), model="rw2" )
    #   + f( spatial.field, model=SPDE ) )
    p$spacetime_model_distance_weighted = TRUE

    p$spacetime_covariate_modeltype="gam"
    p$spacetime_covariate_modelformula = p$spacetime_engine_modelformula

    p$variables = list( Y="t", LOCS=c("plon", "plat"), COV=c("z", "dZ", "ddZ") )
    p$spacetime_family = gaussian(link="log")
   ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one -order of magnitude
    # or to make your own
    # p$spacetime_family = function(offset=0) {
    #   structure(list(
    #     linkfun = function(mu) mu + offset,
    #     linkinv = function(eta) eta - offset,
    #     mu.eta = function(eta) NA,
    #     valideta = function(eta) TRUE,
    #     name = paste0("logexp(", offset, ")") ),
    #     class = "link-glm" )
    #

    p$dist.max = 50 # length scale (km) of local analysis .. for acceptance into the local analysis/model
    p$dist.min = 2 # lower than this .. subsampling occurs

    p$n.min = 30 # n.min/n.max changes with resolution: at p$pres=0.25, p$dist.max=25: the max count expected is 40000
    p$n.max = 8000 # numerical time/memory constraint -- anything larger takes too much time

    p$sbbox = spacetime_db( p=p, DS="statistics.box" ) # bounding box and resoltuoin of output statistics defaults to 1 km X 1 km
    p$non_convex_hull_alpha = 20  # radius in distance units (km) to use for determining boundaries
    p$theta = p$pres # FFT kernel bandwidth (SD of kernel) required for method "harmonic.1/kernel.density"

    p$spacetime.noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
    p$quantile_bounds = c(0.001, 0.999) # remove these extremes in interpolations

    return(p)
  }

}



