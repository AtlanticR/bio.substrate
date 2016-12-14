
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


  if (DS=="lstfilter") {

    p$libs = RLibrary( c( p$libs, "lstfilter" ) ) # required for parallel
    p$storage.backend="bigmemory.ram"

    p$lstfilter_rsquared_threshold = 0.1 # lower threshold
    p$lstfilter_distance_prediction = 7.5 # this is a half window km
    p$lstfilter_distance_statsgrid = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )

    p$sampling = c( 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.1, 1.2, 1.5, 1.75, 2 )  
    
    p$lstfilter_local_modelengine = "gam" 
    p$lstfilter_local_modelformula = formula(
      substrate ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=100, bs="ts") )

    p$lstfilter_local_family = gaussian(link="log")
    p$lstfilter_local_model_distanceweighted = TRUE

    p$lstfilter_global_modelengine = "gam"  # if ==NULL, this means skip global model
    p$lstfilter_global_modelformula = formula(
      substrate ~ s(plon,k=3, bs="ts") + s(plat, k=3, bs="ts") + s(plon, plat, k=100, bs="ts") +
      s(z, k=3, bs="ts") + s(dZ, k=3, bs="ts" ) + s(ddZ, k=3, bs="ts" ) )
    p$lstfilter_global_family = gaussian()

    p$variables = list( Y="grainsize", LOCS=c("plon", "plat"), COV=c("z", "dZ", "ddZ") )
    p$n.min = 50 # n.min/n.max changes with resolution
    p$n.max = 1000 # numerical time/memory constraint -- anything larger takes too much time

    p$lstfilter_nonconvexhull_alpha = 20  # radius in distance units (km) to use for determining boundaries
    p$lstfilter_phi = p$pres/5  # FFT-based methods when operating globally
    p$lstfilter_nu = 0.5

    p$lstfilter_noise = 0.001  # distance units for eps noise to permit mesh gen for boundaries
    p$lstfilter_quantile_bounds = c(0.01, 0.99) # remove these extremes in interpolations

    p$boundary = TRUE 
    p$depth.filter = log(1) # depth is given as log(depth) so, choose andy stats locations with elevation > 1 m as being on land

    return(p)
  }

}



