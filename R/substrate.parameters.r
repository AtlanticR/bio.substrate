
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

    p = spatial_parameters( type=resolution, p=p ) # highest resolution still
    p = spatial_parameters(p)  # load defaults

    # cluster definition
    p$clusters = rep( "localhost", nc )

    return(p)
  }


  if (DS=="sthm") {

    p$libs = RLibrary( c( p$libs, "sthm" ) ) # required for parallel
    
    p$sthm_rsquared_threshold = 0.3 # lower threshold
    p$sthm_distance_prediction = 5 # this is a half window km
    p$sthm_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )

    p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )  
    
    p$sthm_local_modelengine = "gam" # see model form in sthm.r (method="xyts")
    p$sthm_local_modelformula = formula(
      substrate ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=40, bs="ts") +
      s(z, k=3, bs="ts") + s(dZ, k=3, bs="ts" ) + s(ddZ, k=3, bs="ts" ) )
    # p$sthm_local_modelformula = formula( substrate ~ -1 + intercept
    #   + f( inla.group(log(z+1000) ), model="rw2")
    #   + f( inla.group(log(dZ+0.01)), model="rw2")
    #   # + f( inla.group( log(ddZ+0.01) ), model="rw2")
    #   # + f( inla.group( log(Z.rangeMode+0.01)), model="rw2" )
    #   + f( spatial.field, model=SPDE ) )
    p$sthm_local_family = gaussian(link="log")
    p$sthm_local_model_distanceweighted = TRUE

    p$sthm_global_modelengine = NULL  # this means skip global model
    p$sthm_global_modelformula = NULL
    p$sthm_global_family = NULL

    p$variables = list( Y="t", LOCS=c("plon", "plat"), COV=c("z", "dZ", "ddZ") )
   ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one -order of magnitude
    # or to make your own
    # p$sthm_local_family = function(offset=0) {
    #   structure(list(
    #     linkfun = function(mu) mu + offset,
    #     linkinv = function(eta) eta - offset,
    #     mu.eta = function(eta) NA,
    #     valideta = function(eta) TRUE,
    #     name = paste0("logexp(", offset, ")") ),
    #     class = "link-glm" )
    #
    p$n.min = 30 # n.min/n.max changes with resolution
    p$n.max = 8000 # numerical time/memory constraint -- anything larger takes too much time

    p$sthm_sbox = sthm_db( p=p, DS="statistics.box" ) # bounding box and resoltuoin of output statistics defaults to 1 km X 1 km
    p$sthm_nonconvexhull_alpha = 20  # radius in distance units (km) to use for determining boundaries
    p$sthm_theta = p$pres # FFT kernel bandwidth (SD of kernel) required for method "harmonic.1/kernel.density"
    p$sthm_noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
    p$quantile_bounds = c(0.001, 0.999) # remove these extremes in interpolations

    return(p)
  }

}



