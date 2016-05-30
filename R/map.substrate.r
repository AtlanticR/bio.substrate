
  map.substrate = function( params=NULL ) {

    substrate = substrate.db( p, DS="lonlat.grid")
          
    M.SS = subselect.xy2grid (area="cfaall", DS="file", loc="grids", fname="mask.cfa")
    # screen out the Bay of Fundy (4X) area
    M.4X = get.boxes ( area="4X", DS="file", fname=file.path(project.datadirectory("snowcrab"), "snowcrab", "R", "grids", "mask.4X.rdata" ) )
    M.SS = M.SS * M.4X
    rm(M.4X)

    Z = bathymetry.db(p=p, DS="Z.planar.grid")
    attr(Z, "dimnames") = NULL
    M.ZSS = ifelse( Z < 5 | Z > 500, NA, 1) # identify the domain of the Scotian Shelf
    M.SS = M.SS * M.ZSS
    rm(M.ZSS);gc()

    X = substrate * M.SS

    image(log(X))
    Pr( "png", file.path(project.datadirectory("substrate"), "R"), "substrate" ) 

    return (NULL)
  }


