
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

    p = spatial_parameters( type=resolution, p=p ) # highest resolution still

    return(p)
  }


  if (DS=="hivemod") {

    p$libs = RLibrary( c( p$libs, "hivemod" ) ) # required for parallel
    p$storage.backend="bigmemory.ram"
    p$boundary = TRUE 
    p$depth.filter = log(1) # depth is given as log(depth) so, choose andy stats locations with elevation > 1 m as being on land

    p$hivemod_nonconvexhull_alpha = 20  # radius in distance units (km) to use for determining boundaries
    p$hivemod_lowpass_phi = p$pres/5 # FFT based method when operating gloablly
    p$hivemod_lowpass_nu = 0.5 # this is exponential covar
    p$hivemod_noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
    p$hivemod_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations
    
    p$hivemod_rsquared_threshold = 0.05 # lower threshold
    p$hivemod_distance_prediction = 7.5 # this is a half window km
    p$hivemod_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
    p$hivemod_distance_scale = 25 # km ... approx guess of 95% AC range 
    p$hivemod_distance_min = 1
    p$hivemod_distance_max = 50 

    p$n.min = 50 # n.min/n.max changes with resolution
    p$n.max = 2000 # numerical time/memory constraint -- anything larger takes too much time
    p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )  # fractions of median distance scale to try in local block search
 
    p$variables = list( Y="grainsize", LOCS=c("plon", "plat"), COV=c("z", "dZ", "ddZ") )
 
    p$hivemod_global_modelengine = "gam"  # if ==NULL, this means skip global model
    p$hivemod_global_modelformula = formula(
      grainsize ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=100, bs="ts") +
      s(z, k=3, bs="ts") + s(dZ, k=3, bs="ts" ) + s(ddZ, k=3, bs="ts" ) )
    p$hivemod_global_family = gaussian(link="log")
    
    if (!exists("hivemod_variogram_method", p)) p$hivemod_variogram_method = "fast"
    if (!exists("hivemod_local_modelengine", p)) p$hivemod_local_modelengine="gam"  # currently the perferred approach 

    if (p$hivemod_local_modelengine == "gam") {
      # GAM are overly smooth .. adding more knots might be good but speed is the cost .. k=50 to 100 seems to work nicely
      ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one order of magnitude
      p$hivemod_local_modelformula = formula( 
        grainsize ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=100, bs="ts") )  
      p$hivemod_local_model_distanceweighted = TRUE  
      p$hivemod_local_family = gaussian()
     
    } else if ( p$hivemod_local_modelengine == "bayesx" ) {
    
      ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one order of magnitude
      p$hivemod_local_modelformula = formula( grainsize ~ s(plon, bs="ps") + s(plat, bs="ps") + s(plon, plat, bs="te") )  # more detail than "gs" .. "te" is preferred
      p$hivemod_local_model_bayesxmethod="MCMC"  # REML actually seems to be the same speed ... i.e., most of the time is spent in thhe prediction step ..
      p$hivemod_local_model_distanceweighted = TRUE  
      p$hivemod_local_family = gaussian()
      p$hivemod_local_family_bayesx ="gaussian"

    } else {
    
      message( "The specified hivemod_local_modelengine is not tested/supported ... you are on your own ;) ..." )

    }
    

      
    return(p)
  }

}



