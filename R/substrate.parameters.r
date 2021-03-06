
substrate.parameters = function(DS="bio.substrate", p=NULL, resolution="canada.east.highres") {

  if ( is.null(p) ) p=list()
  if ( !exists("project.name", p) ) p$project.name=DS


  if (DS=="bio.substrate"){

    p$project.root = project.datadirectory( p$project.name )

    p$libs = bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons",
        "bio.substrate", "bio.coastline" )
    p$libs = c( p$libs, RLibrary( "rgdal", "maps", "mapdata", "maptools", "lattice", "parallel", 
                       "geosphere", "geoR", "gstat", "spBayes", "RandomFields", 
                       "sp", "raster", "colorspace" ,  "splancs", "fields" ) )

    p = spatial_parameters( type=resolution, p=p ) # specified in param default, above ..
    p$spatial.domain.subareas = c( "canada.east.highres", "canada.east", "SSE", "snowcrab", "SSE.mpa" )
  
    return(p)
  }


  if (DS=="lbm") {

    p$libs = RLibrary( c( p$libs, "lbm" ) ) # required for parallel
    p$storage.backend="bigmemory.ram"
    if (!exists("clusters", p)) p$clusters = rep("localhost", detectCores() )
  
    p$boundary = FALSE 
    p$depth.filter = log(1) # the depth covariate is input as log(depth) so, choose andy stats locations with elevation > log(1 m) as being on land

    p$lbm_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations
    
    p$lbm_rsquared_threshold = 0.25 # lower threshold
    p$lbm_distance_prediction = 4 # this is a half window km
    p$lbm_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    p$lbm_distance_scale = 25 # km ... approx guess of 95% AC range 
    p$lbm_distance_min = p$lbm_distance_statsgrid 
    p$lbm_distance_max = 50 

    p$n.min = 100 # n.min/n.max changes with resolution
    p$n.max = 3000 # numerical time/memory constraint -- anything larger takes too much time .. anything less .. errors
    p$sampling = c( 0.5, 0.75, 0.9, 1.1, 1.25, 1.5 )  # fractions of median distance scale to try in local block search
 
    p$variables = list( Y="log.substrate.grainsize", LOCS=c("plon", "plat"), COV=c("z", "dZ", "ddZ") )
    p$varnames = c( p$variables$LOCS, p$variables$COV ) # to retain fom data sources

    # global model to permit extrapolation ... no space here only in the local model
    p$lbm_global_modelengine = "gam"  # if ==NULL, this means skip global model
    p$lbm_global_modelformula = formula(
      log.substrate.grainsize ~  s( log(z), k=3, bs="ts") + s( log(dZ), k=3, bs="ts" ) + s( log(ddZ), k=3, bs="ts" ) )
    

    p$lbm_global_family = gaussian() # Y-var already log transformed
    p$lbm_local_family = gaussian() # residuals are already log-tranformed so expect gaussian ..

    if (!exists("lbm_variogram_method", p)) p$lbm_variogram_method = "fast"
    if (!exists("lbm_local_modelengine", p)) p$lbm_local_modelengine="krige"  # currently the perferred approach 

    if ( p$lbm_local_modelengine %in% c("krige" )) { 
      
      # nothing to do ... faster than gstat

    } else if (p$lbm_local_modelengine == "gam") {
      # GAM are overly smooth .. adding more knots might be good but speed is the cost .. k=50 to 100 seems to work nicely
      ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one order of magnitude
      p$lbm_local_modelformula = formula( 
        log.substrate.grainsize ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=200, bs="ts") )  
      p$lbm_local_model_distanceweighted = TRUE  
      p$lbm_gam_optimizer ="perf"

    } else if ( p$lbm_local_modelengine == "bayesx" ) {
    
      p$lbm_local_modelformula = formula( log.substrate.grainsize ~ s(plon, bs="ps") + s(plat, bs="ps") + s(plon, plat, bs="te") )  # more detail than "gs" .. "te" is preferred
      p$lbm_local_model_bayesxmethod="MCMC"  # REML actually seems to be the same speed ... i.e., most of the time is spent in thhe prediction step ..
      p$lbm_local_model_distanceweighted = TRUE  
      p$lbm_local_family ="gaussian"

    } else if ( p$lbm_local_modelengine == "fft" ) {
    
      p$lbm_lowpass_phi = p$pres*2 # FFT based method when operating gloablly
      p$lbm_lowpass_nu = 0.5 # this is exponential covar

    } else {
    
      message( "The specified lbm_local_modelengine is not tested/supported ... you are on your own ;) ..." )

    }
    
    return(p)
  }

}



