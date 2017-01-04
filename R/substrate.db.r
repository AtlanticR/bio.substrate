

  substrate.db = function( p=NULL, DS=NULL, varnames=NULL ) {

    if ( DS %in% c("substrate.initial", "substrate.initial.redo") ) {
      # Read in the ArcInfo ascii grid file using library maptools and output a SpatialGridDataFrame
      # data provided by Kostelev:
      # Kostylev, V.E., and Hannah, C.G., 2007, Process-driven characterization and mapping of seabed habitats,
      # in Todd, B.J.,and Greene, H.G., eds., Mapping the Seafloor for Habitat Characterization:
      # Geological Association of Canada, Special Paper 47, p. 171-184.
      # Scotian shelf gridded grain size (mm).
      # NAD83 UTM zone 20 (I think)

      rawdata.file = file.path( project.datadirectory("bio.substrate"), "data", "grainsize.txt" )
      filename = file.path( project.datadirectory("bio.substrate"), "data", "substrate.asciigrid.rdata" )

      if (DS =="substrate.initial" ) {
        load( filename )
        return ( substrate )
      }
      proj4.params = "+proj=utm +zone=20 +ellps=GRS80  +datum=NAD83 +units=m" #resolution is 500m X 500m
      substrate = sp::read.asciigrid( rawdata.file, proj4string=CRS( proj4.params ), colname="grainsize" )  ## mm
      save( substrate, file=filename, compress=T )
      return(filename)
    }

    # lon - lat converted
    if (  DS %in% c("lonlat.highres", "lonlat.highres.redo") ) {
      filename = file.path( project.datadirectory("bio.substrate"), "data", "substrate.lonlat.highres.rdata" )
      if (DS =="lonlat.highres" ) {
        load( filename)
        return( substrate)
      }
      # initial data stored in planar coords ... convert to lon/lats
      substrate = substrate.db( DS="substrate.initial" )
      substrate = as.data.frame( substrate )
      names(substrate) = c("grainsize", "plon", "plat" )
      substrate = substrate[,c("plon", "plat", "grainsize")]
      proj4.params = "+proj=utm +zone=20 +ellps=GRS80 +datum=NAD83 +units=m"  # original/raw data still in NAD83 geoid
      substrate= planar2lonlat ( substrate, proj4.params )
      substrate= substrate[ ,c("lon", "lat", "grainsize")]
      save( substrate, file=filename, compress=T   )
      return ( filename )
    }


    # ---------------------------------------


    # ---------------------------------------


    if ( DS %in% c("lonlat.interpolated", "lonlat.interpolated.redo") ) {
      # interpolation to internal grid
      # locally (internally) force the highest possible resolution to not lose data and extrapolate safely
      filename.lonlat.interp = file.path( project.datadirectory("bio.substrate"), "data",
					paste( p$spatial.domain, "substrate.lonlat.interpolated.rdata", sep=".")  )
      if (DS =="lonlat.interpolated" ) {
        load (filename.lonlat.interp )
        return( substrate)
      }
      substrate = substrate.db( p, DS="lonlat.highres" )

      p$res = "-I10s"  # 10 arc sec -- ie. all data
      p$tension = "-T1" # interpolated but minimally smoothed solutions
      rlons = range(p$lons)
      rlats = range(p$lats)
      p$region = paste("-R", rlons[1], "/", rlons[2], "/", rlats[1], "/", rlats[2], sep="")

			names.sub = colnames( substrate )
      grainsize.range = range( substrate$grainsize, na.rm=T )
      substrate = interpol.grid(xyz=substrate, params=p, getdata=T, method="tps" )
      colnames( substrate ) = names.sub  # the above function renames the 3rd var to z
      # interpolation can bring in data larger or smaller than realistic
      substrate$grainsize[ which( ( substrate$grainsize < grainsize.range[1] )) ] = grainsize.range[1]
      substrate$grainsize[ which( ( substrate$grainsize > grainsize.range[2] )) ] = grainsize.range[2]
      save( substrate, file=filename.lonlat.interp, compress=T )
      return( filename.lonlat.interp )
    }

    # ------------

    if ( DS %in% c("lonlat", "lonlat.redo", "lonlat.grid") ) {
      # interpolation to internal grid
      # locally (internally) force the highest possible resolution to not lose data
      filename.lonlat = file.path( project.datadirectory("bio.substrate"), "data",
					paste( p$spatial.domain, "substrate.lonlat.rdata", sep=".") )
      filename.lonlat.grid = file.path( project.datadirectory("bio.substrate"), "data",
					paste( p$spatial.domain, "substrate.lonlat.grid.rdata", sep=".") )

      if (DS =="lonlat.grid" ) {
        load( filename.lonlat.grid )
        return( substrate )
      }
      if (DS =="lonlat" ) {
        load( filename.lonlat )
        return( substrate )
      }

			ilons = c( p$lons, p$lons[length(p$lons)]+(p$lons[2]-p$lons[1]) )
      ilats = c( p$lats, p$lats[length(p$lats)]+(p$lats[2]-p$lats[1]) )

      substrate = substrate.db( p, DS="lonlat.interpolated" )
      substrate$lon = as.numeric(as.character(cut(substrate$lon, ilons, include.lowest=T, right=F, labels=p$lons)))
      substrate$lat = as.numeric(as.character(cut(substrate$lat, ilats, include.lowest=T, right=F, labels=p$lats)))

      gc()
      substrate = block.spatial ( xyz=substrate, function.block=block.mean )
      save( substrate, file=filename.lonlat, compress=T )

      gc()
      substrate = xyz2grid( substrate, p$lons, p$lats)
      save( substrate, file=filename.lonlat.grid, compress=T )
      return ( paste( filename.lonlat.grid, filename.lonlat, sep="\n") )
    }


    if ( DS %in% c("planar", "planar.redo", "planar.grid") ) {
      # Re-grid data to be internally consistent with the snowcrab coordinate system
      # WGS84 ellipsoid and not NAD83 ...
      filename.planar = file.path( project.datadirectory("bio.substrate"), "data", paste( p$spatial.domain, "substrate.planar.rdata", sep=".") )
      filename.planar.grid = file.path( project.datadirectory("bio.substrate"), "data", paste( p$spatial.domain, "substrate.planar.grid.rdata", sep=".") )

      if (DS =="planar.grid" ) {
        load( filename.planar.grid )
        return( substrate )
      }
      if (DS =="planar" ) {
        load( filename.planar)
        return( substrate )
      }

      substrate = substrate.db( p, DS="lonlat.interpolated" )
      substrate = lonlat2planar( substrate,  proj.type=p$internal.projection )  # utm20, WGS84 (snowcrab geoid)
      substrate = substrate[ ,c("plon", "plat", "grainsize" )]

			substrate$plon = grid.internal( substrate$plon, p$plons )
      substrate$plat = grid.internal( substrate$plat, p$plats )

      gc()
      substrate = block.spatial ( xyz=substrate, function.block=block.mean)
      save( substrate, file=filename.planar, compress=T )

      gc()
      substrate = xyz2grid(substrate, p$plons, p$plats)
      save( substrate, file=filename.planar.grid, compress=T )
      return ( paste( filename.planar.grid, filename.planar, sep="\n" ) )
    }

    #-------------------------

    if ( DS=="lbm.inputs") {

      B = bathymetry.db( p, DS="baseline", varnames=p$varnames )
      B$z = log( B$z) # ranges  are too large in some cases to use untransformed 2 orders or more (e.g. 40 to 2000 m)
      bid = lbm::array_map( "2->1", trunc(cbind(B$plon-p$plon[1], B$plat-p$plat[1])/p$pres) + 1, c(p$nplons,p$nplats) )

      S = substrate.db( p=p, DS="lonlat.highres" ) 
      S = lonlat2planar( S,  proj.type=p$internal.projection )  # utm20, WGS84 (snowcrab geoid)
      S = S[ ,c("plon", "plat", "grainsize" )]
      S$log.substrate.grainsize = log( S$grainsize )
      S$grainsize = NULL

      # merge covars into S
      sid = lbm::array_map( "2->1", trunc(cbind(S$plon-p$plon[1], S$plat-p$plat[1])/p$pres) + 1, c(p$nplons, p$nplats) )
      u = match( sid, bid )
      B_matched = B[u, ]
      S = cbind(S, B_matched )
      S = S[ is.finite( rowSums(S) ), ]
      OUT  = list( LOCS=B[, p$variables$LOCS], COV =B[,  p$variables$COV ] )         

      return(  list( input=S, output=OUT ) )

    }

    #-------------------------

    if ( DS %in% c("lbm.finalize.redo", "lbm.finalize" )) {
      #// substrate( p, DS="lbm.finalize(.redo)" return/create the
      #//   lbm interpolated method formatted and finalised for production use with predictions and statistics
      fn = file.path(  project.datadirectory("bio.substrate"), "interpolated",
        paste( "substrate", "lbm", "finalized", p$spatial.domain, "rdata", sep=".") )
      if (DS =="lbm.finalize" ) {
        B = NULL
        if ( file.exists ( fn) ) load( fn)
        return( B )
      }

      B = expand.grid( p$plons, p$plats, KEEP.OUT.ATTRS=FALSE)
      names( B ) = c("plon", "plat")
      Bmean = lbm_db( p=p, DS="lbm.prediction", ret="mean" )
      Bsd = lbm_db( p=p, DS="lbm.prediction", ret="sd" )
      B = cbind(B, Bmean, Bsd)
      rm (Bmean, Bsd); gc()
      names(B) = c( "plon", "plat", "grainsize", "grainsize.sd") 

      # remove land
      oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="land", tag="predictions" ,internal.projection=p$internal.projection)
      B$grainsize[oc] = NA
      B$grainsize.sd[oc]   = NA

      rm(preds); gc()

      # merge into statistics
      BS = lbm( p=p, DS="lbm.statistics" )
      B = cbind( B, BS )
      # names(B) = c( names(B), p$statsvars )

      save( B, file=fn, compress=TRUE)
      return(fn)

      if (0) {
        levelplot( log(grainsize) ~ plon + plat, B, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        levelplot( log(grainsize.range) ~ plon + plat, B, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        levelplot( grainsize.sd ~ plon + plat, B, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      }
    }


    # ------------


    if ( DS %in% c( "complete", "complete.redo" ) ) {
     #// substrate.db( DS="complete" .. ) returns the final form of the substrate data after
     #// regridding and selection to area of interest as specificied by girds.new=c("SSE", etc)
      Z = NULL

      if ( DS %in% c("complete") ) {

        domain = NULL
        if ( is.null(domain)) {
          if ( exists("spatial.domain", p)) {
            domain = p$spatial.domain
          } else if ( exists( "new.grids", p) )  { # over-rides p$spatial domain
            if( length( p$new.grids )== 1 ) {
              domain = p$new.grids
        } } }
        fn = file.path( project.datadirectory("bio.substrate", "interpolated"),
          paste( "substrate", "complete", domain, "rdata", sep=".") )
        if ( file.exists ( fn) ) load( fn)

        Znames = names(Z)
        if (is.null(varnames)) {
          varnames=Znames
        } else {
          varnames = intersect( Znames, varnames )
        }
        Z = Z[ , varnames]

        return( Z )
      }

      p0 = p  # the originating parameters

      Z0 = substrate.db( p=p0, DS="lbm.finalize" )
      coordinates( Z0 ) = ~ plon + plat
      crs(Z0) = crs( p0$interal.crs )
      
      grids = unique( c( p$spatial.domain, p$new.grids ))

      for (gr in grids ) {
        Z = NULL
        print(gr)
        p1 = spatial_parameters( type=gr )
        for (vn in names(Z0)) {
          Z[[vn]] = raster::projectRaster(
            from = raster::rasterize( Z0, bio.spacetime::spatial_parameters_to_raster(p0), field=vn, fun=mean),
            to   = bio.spacetime::spatial_parameters_to_raster( p1) )
        }
        Z = as( brick(Z), "SpatialPointsDataFrame" )
        Z = as.data.frame(Z)
        u = names(Z)
        names(Z)[ which( u=="x") ] = "plon"
        names(Z)[ which( u=="y") ] = "plat"
        fn = file.path( project.datadirectory("bio.substrate", "interpolated"),
          paste( "substrate", "complete", p1$spatial.domain, "rdata", sep=".") )
        save (Z, file=fn, compress=TRUE)
      }
      return ( "Completed subsets" )
    }

  }


