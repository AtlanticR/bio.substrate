
substrate.parameters = function(DS, p=NULL, resolution="canada.east.highres", nc=1) {

  if ( is.null(p) ) p=list()
  if ( !exists("project.name", p) ) p$project.name=DS

  p$project.root = project.datadirectory( p$project.name )

  p$libs = bioLibrary( "bio.spacetime", "bio.utilities", "bio.bathymetry", "bio.polygons",
      "bio.substrate", "bio.coastline" )
  p$libs = c( p$libs, RLibrary( "rgdal", "maps", "mapdata", "maptools", "lattice", "parallel", "INLA",
                     "geosphere", "geoR", "gstat", "spBayes",
                     "sp", "raster", "colorspace" ,  "splancs", "fields", "bigmemory" ) )

  p = spatial.parameters( type=resolution, p=p ) # highest resolution still
  p = spacetime.parameters(p)  # load defaults
  p$substrate.bigmemory.reset = FALSE

  # cluster definition
  p$clusters = rep( "localhost", nc )

  return(p)

}



