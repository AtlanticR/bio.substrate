
substrate.figures = function( DS=NULL, p=NULL ) {

  bioLibrary( "bio.coastline" )

  if ( DS=="predictions" ) {

    P = p$ff$P
    Ploc = p$ff$Ploc
    p$spatial.domain="canada.east"  # force isobaths to work in levelplot
    datarange = range( P[ which(is.finite(P[,2])),2] )
    dr = seq( datarange[1], datarange[2], length.out=100)
    oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="not.land", tag="predictions" )
    levelplot( log( P[oc,2] ) ~ plons + plats, Ploc[oc,], aspect="iso", main=NULL, at=dr, col.regions=rev(color.code( "seis", dr)) ,
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE),
        panel = function(x, y, subscripts, ...) {
          panel.levelplot (x, y, subscripts, aspect="iso", rez=c(1,1), ...)
          sp.lines( isobath.db( p=p, DS="isobath", depths=c(200, 400 ), crs=p$internal.crs ), col = "slateblue", cex=0.1 )
          sp.lines( coastline.db( p=p), crs=p$internal.crs, col = "steelblue", cex=0.1 )
      }
    )
    close( Ploc )
    close( P)
 }

  if ( DS=="predictions.error" ) {
    
    P = p$ff$P
    Ploc = p$ff$Ploc
    p$spatial.domain="canada.east"  # force isobaths to work in levelplot
    datarange = log( c( 2, 50 ))
    dr = seq( datarange[1], datarange[2], length.out=100)
    oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="not.land", tag="predictions" )
    levelplot( log( P[oc,3] ) ~ plons + plats, Ploc[oc,], aspect="iso", main=NULL, at=dr, col.regions=rev(color.code( "seis", dr)) ,
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE),
        panel = function(x, y, subscripts, ...) {
          panel.levelplot (x, y, subscripts, aspect="iso", rez=c(1,1), ...)
          sp.lines( isobath.db( p=p, DS="isobath", depths=c(200, 400 ), crs=p$internal.crs  ), col = "slateblue", cex=0.1 )
          sp.lines( coastline.db( p=p, crs=p$internal.crs), col = "steelblue", cex=0.1 )
      }
    )
    close( Ploc )
    close( P)
  }

  if ( DS=="statistics" ) {
    
    S = p$ff$S
    Sloc = p$ff$Sloc
    p$spatial.domain="canada.east"  # force isobaths to work in levelplot
    datarange = log( c( 5, 800 ))
    dr = seq( datarange[1], datarange[2], length.out=150)
    oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="not.land", tag="statistics" )
    levelplot( log(S[oc,1])  ~ Sloc[oc,1] + Sloc[oc,2] , aspect="iso", at=dr, col.regions=color.code( "seis", dr) ,
      contour=FALSE, labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE), cex=2,
      panel = function(x, y, subscripts, ...) {
        panel.levelplot (x, y, subscripts, aspect="iso", rez=c(5,5), ...)
        sp.lines( isobath.db( p=p, DS="isobath", depths=c( 200, 400 ), crs=p$internal.crs  ), col = "slateblue", cex=0.1 )
        sp.lines( coastline.db( p=p, crs=p$internal.crs), col = "steelblue", cex=0.1 )
      }
    )
    close( Sloc )
    close( S)
  }




}


