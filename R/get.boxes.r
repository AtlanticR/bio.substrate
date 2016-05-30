
get.boxes = function(area=NULL, DS="file", lons=NULL, lats=NULL, fname=NULL) {

  if (DS=="file") { load(fname) }

  if (DS=="redo") {
    ss = expand.grid(x=lons, y=lats)
    attr(ss, "out.attrs") = NULL  # remove uneeded data
    names(ss) = c("lon", "lat")
    G = xyz2grid ( xyz=ss, xx=lons, yy=lats )
    if (area== "4X") {
      # blank area 1
        k = range( which(lons <= -66 & lons >= -68 ) )
        l = range( which(lats >= 44 & lats <= 45.5 ) )
        G[ c(k[1]:k[2]), c(l[1]:l[2])] = NA
      # blank area 2
        k = range( which(lons <= -64.5 & lons >= -66 ) )
        l = range( which(lats >= 44.5 & lats <= 46 ) )
        G[ c(k[1]:k[2]), c(l[1]:l[2])] = NA
      # blank area 3
        k = range( which(lons <= -63 & lons >= -64.5 ) )
        l = range( which(lats >= 45 & lats <= 45.5 ) )
        G[ c(k[1]:k[2]), c(l[1]:l[2])] = NA
    }
    dir.create( dirname(fname) )
    save( G, file=fname, compress=T )
  }
  return( G )
}


