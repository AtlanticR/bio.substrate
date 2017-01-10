

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

      gridparams = list( dims=c(p$nplons, p$nplats), corner=c(p$plons[1], p$plats[1]), res=c(p$pres, p$pres) )

      B = bathymetry.db( p, DS="baseline", varnames=p$varnames )
      B$z = log( B$z) # ranges  are too large in some cases to use untransformed 2 orders or more (e.g. 40 to 2000 m)

      bid = lbm::array_map( "xy->1", B[,c("plon", "plat")], gridparams=gridparams )

      S = substrate.db( p=p, DS="lonlat.highres" ) 
      S = lonlat2planar( S,  proj.type=p$internal.projection )  # utm20, WGS84 (snowcrab geoid)
      S = S[ ,c("plon", "plat", "grainsize" )]
      S$log.substrate.grainsize = log( S$grainsize )
      S$grainsize = NULL

      # merge covars into S
      sid = lbm::array_map( "xy->1", S[,c("plon", "plat")], gridparams=gridparams )

      u = match( sid, bid )
      B_matched = B[u, ]
      S = cbind(S, B_matched )
      S = S[ is.finite( rowSums(S) ), ]
      OUT  = list( LOCS=B[, p$variables$LOCS], COV =B[,  p$variables$COV ] )         

      return(  list( input=S, output=OUT ) )

    }

    #-------------------------


    if ( DS %in% c( "complete", "complete.redo" ) ) {
     #// substrate.db( DS="complete" .. ) returns the final form of the substrate data after
     #// regridding and selection to area of interest as specificied by girds.new=c("SSE", etc)

      if ( DS %in% c("complete") ) {
        S = NULL
        fn = file.path( project.datadirectory("bio.substrate", "modelled"),
          paste( "substrate", "complete", p$spatial.domain, "rdata", sep=".") )
        if ( file.exists ( fn) ) load( fn)
        Snames = names(S)
        if (is.null(varnames)) {
          varnames=Snames
        } else {
          varnames = intersect( Snames, varnames )
        }
        S = S[ , varnames]
        return( S )
      }

      Smean = lbm_db( p=p, DS="lbm.prediction", ret="mean" )
      Ssd = lbm_db( p=p, DS="lbm.prediction", ret="sd" )
      S = as.data.frame( cbind( Smean, Ssd) )
      rm (Smean, Ssd); gc()
      names(S) = c( "log.substrate.grainsize", "log.substrate.grainsize.sd") 


      if (0) {
        B = bathymetry.db(p=p, DS="baseline")
        levelplot( log(S[,1]) ~ plon + plat, B, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        levelplot( log(S[,2]) ~ plon + plat, B, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      }

      if (0) {
        # remove land
        oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="land", tag="predictions" ,internal.projection=p$internal.projection)
        S$[oc,] = NA
        rm(oc); gc()
      }

      # merge into statistics
      SS = lbm_db( p=p, DS="stats.to.prediction.grid" )
      S = cbind( S, SS )

      fn = file.path( project.datadirectory("bio.substrate", "modelled"),
        paste( "substrate", "complete", p$spatial.domain, "rdata", sep=".") )
      save (S, file=fn, compress=TRUE)

      p0 = p  # the originating parameters
      S0 = S
      L0 = bathymetry.db( p=p0, DS="baseline" )
      L0i = array_map( "xy->2", L0, 
        corner=c(p0$plons[1], p0$plats[1]), res=c(p0$pres, p0$pres) )
   
      varnames = setdiff( names(S0), c("plon","plat", "lon", "lat") )  
      #using fields
      grids = setdiff( unique( p0$new.grids ), p0$spatial.domain )
      for (gr in grids ) {
        print(gr)
        p1 = spatial_parameters( type=gr ) #target projection
        L1 = bathymetry.db( p=p1, DS="baseline" )
        L1i = array_map( "xy->2", L1[, c("plon", "plat")], 
          corner=c(p1$plons[1], p1$plats[1]), res=c(p1$pres, p1$pres) )
        L1 = planar2lonlat( L1, proj.type=p1$internal.crs )
        S = L1
        L1$plon_1 = L1$plon # store original coords
        L1$plat_1 = L1$plat
        L1 = lonlat2planar( L1, proj.type=p0$internal.crs )
        p1$wght = fields::setup.image.smooth( 
          nrow=p1$nplons, ncol=p1$nplats, dx=p1$pres, dy=p1$pres,
          theta=p1$pres, xwidth=4*p1$pres, ywidth=4*p1$pres )
        for (vn in varnames) {
          S[[vn]] = spatial_warp( S0[,vn], L0, L1, p0, p1, L0i, L1i )
        }
    
        S = S[ , names(S0) ]
        fn = file.path( project.datadirectory("bio.substrate", "modelled"),
          paste( "substrate", "complete", p1$spatial.domain, "rdata", sep=".") )
        save (S, file=fn, compress=TRUE)
      }

      if(0){
        datarange = quantile( S$log.substrate.grainsize, probs=c(0.01, 0.99), na.rm=TRUE )
        dr = seq( datarange[1], datarange[2], length.out=100)
        levelplot(log.substrate.grainsize~plon+plat, S,  at=dr, col.regions=(color.code( "seis", dr)))

      }
      return ( "Completed subsets" )
    }


  }


