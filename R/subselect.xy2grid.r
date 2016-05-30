
  subselect.xy2grid = function(xy=NULL, area="cfaall", DS="file", loc="grids", fname="cfa.grid" ) {
    fname = file.path(loc, paste(fname, area, "rdata", sep="."))
    if (DS=="redo") {
      inside = filter.region.polygon(xy, area)
      xy$areafilter = NA
      xy$areafilter[inside] = 1
      grid = xyz2grid ( xyz=xy, xx=lons, yy=lats)

      loc = dirname( fname )
      dir.create( path=loc, recursive=T, showWarnings=F )
      save( grid, file=file.path(loc,fname), compress=T)

    }
    if (DS=="file") load (fname)
    return(grid)
  }

