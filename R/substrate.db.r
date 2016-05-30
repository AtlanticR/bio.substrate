

  substrate.db = function( p=NULL, DS=NULL ) {

    if ( DS %in% c("substrate.initial", "substrate.initial.redo") ) {
      # Read in the ArcInfo ascii grid file using library maptools and output a SpatialGridDataFrame
      # data provided by Kostelev:
      # Kostylev, V.E., and Hannah, C.G., 2007, Process-driven characterization and mapping of seabed habitats,
      # in Todd, B.J.,and Greene, H.G., eds., Mapping the Seafloor for Habitat Characterization:
      # Geological Association of Canada, Special Paper 47, p. 171-184.
      # Scotian shelf gridded grain size (mm).
      # NAD83 UTM zone 20 (I think)

      rawdata.file = file.path( project.datadirectory("substrate"), "data", "grainsize.txt" )
      filename = file.path( project.datadirectory("substrate"), "data", "substrate.asciigrid.rdata" )

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
      filename = file.path( project.datadirectory("substrate"), "data", "substrate.lonlat.highres.rdata" )
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


    if (DS=="substrate.gmt.retired") {

      ### This uses GMT-based methods .. it is now deprecated

      p = list()
      p$init.files = 	loadfunctions( c("spacetime", "utility", "substrate", "bathymetry", "polygons" ) )
      p$libs = RLibrary( "maptools" , "rgdal" )

      # --------------------------------------
      # create the main database
      # some require upto 1.2 GB RAM, ~ 5 min
      # no need to run again unless the substrate data file is updated ...

      make.substrate.db = FALSE
      if (make.substrate.db) {
        substrate.db ( DS="substrate.initial.redo" ) # stored as a SpatialGridDataFrame
        substrate.db ( DS="lonlat.highres.redo" ) # simple conversion to lonlat
        for ( j in c( "SSE", "canada.east" ) ) {  # sse and snowcrab have same domains
          p = spatial.parameters( type=j )
          substrate.db ( p, DS="lonlat.interpolated.redo" )
          substrate.db ( p, DS="lonlat.redo" )
          substrate.db ( p, DS="planar.redo" )
        }
      }


      # --------------------------------------
      # substrate = substrate.db( p, DS="lonlat.highres" ) # or lonlat to refresh, planar or planar.saved
      # library(lattice)
      # levelplot( log(grainsize) ~ lon + lat, substrate, main = "ln( grainsize; mm )", aspect="iso")

      # load the imported data in a data.frame format in a snow crab-consistent coordinates framework
      p = spatial.parameters( type="SSE" )

      substrate = substrate.db( p, DS="planar" ) # or lonlat to refresh, planar or planar.saved
      i = which( substrate$plon< 990 &  substrate$plon > 220  &
                 substrate$plat< 5270 &  substrate$plat > 4675
      )
      substrate = substrate[ i, ]
      substrate$grainsize =log( substrate$grainsize )
      inside = filter.region.polygon( substrate[, c("plon", "plat") ], "cfaall", planar=T)
      datacols = c("plon", "plat", "grainsize")
      datarange = seq(-5,3, length.out=50)
      cols = color.code( "blue.black", datarange )
      outfn = "substrate.grainsize"
      annot = "ln ( Grain size; mm )"
      map( xyz=substrate[inside,datacols], cfa.regions=F, depthcontours=T, pts=NULL, annot=annot,
        fn=outfn, loc=file.path( project.datadirectory("substrate"), "R"), at=datarange , col.regions=cols )

    }


    # ---------------------------------------


    if ( DS %in% c("lonlat.interpolated", "lonlat.interpolated.redo") ) {
      # interpolation to internal grid
      # locally (internally) force the highest possible resolution to not lose data and extrapolate safely
      filename.lonlat.interp = file.path( project.datadirectory("substrate"), "data",
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
      filename.lonlat = file.path( project.datadirectory("substrate"), "data",
					paste( p$spatial.domain, "substrate.lonlat.rdata", sep=".") )
      filename.lonlat.grid = file.path( project.datadirectory("substrate"), "data",
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
      filename.planar = file.path( project.datadirectory("substrate"), "data", paste( p$spatial.domain, "substrate.planar.rdata", sep=".") )
      filename.planar.grid = file.path( project.datadirectory("substrate"), "data", paste( p$spatial.domain, "substrate.planar.grid.rdata", sep=".") )

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

    # ---------------

    if (DS %in% c("substrate.spacetime.inputs.data.redo", "substrate.spacetime.inputs.data") ) {

      datadir = project.datadirectory("substrate", "data" )
			dir.create( datadir, showWarnings=F, recursive=T )
      fn = file.path( datadir, paste( "substrate", "spacetime", p$spatial.domain, "rdata", sep=".") )

      if (DS =="substrate.spacetime.inputs.data" ) {
        load( fn)
        return( substrate )
      }

      # begin with bathymetry data and add another layer to it
      substrate = bathymetry.db( p, DS="complete", return.format = "list" )
      substrate$substrate = projectRaster(
          from=raster( substrate.db( DS="substrate.initial" ) ),
          to=spatial.parameters.to.raster( p) )
      substrate = as( brick(substrate), "SpatialGridDataFrame" )

      save (substrate, file=fn, compress=TRUE)
      return(fn)
    }


    ### ------


    if (DS %in% c("substrate.spacetime.inputs.prediction.redo", "substrate.spacetime.inputs.prediction") ) {

      datadir = project.datadirectory("substrate", "data" )
			dir.create( datadir, showWarnings=F, recursive=T )
      fn = file.path( datadir, paste( "substrate", "spacetime", p$spatial.domain, "rdata", sep=".") )

      if (DS =="substrate.spacetime.inputs.prediction" ) {
        #load( fn)
        substrate = substrate.db( p, DS="substrate.spacetime.inputs.data" )
        return( substrate )
      }

      ### prediction grids are the same as the input grid .. do nothing for now
      ### but kept separate from "*...inputs" in case thy diverge in future
      print( "This is just a placeholder for more elaborate models ..  for now, grids are the same as inputs")
      # substrate = substrate.db( p, DS="substrate.spacetime.inputs.data" )
      # save (substrate, file=fn, compress=TRUE)
      return(fn)
    }

    #-------------------------

    if ( DS %in% c("substrate.spacetime.finalize.redo", "substrate.spacetime.finalize" )) {
      #// substrate( p, DS="substrate.spacetime.finalize(.redo)" return/create the
      #//   spacetime interpolated method formatted and finalised for production use with predictions and statistics
      fn = file.path(  project.datadirectory("substrate"), "interpolated",
        paste( "substrate", "spacetime", "finalized", p$spatial.domain, "rdata", sep=".") )
      if (DS =="substrate.spacetime.finalize" ) {
        B = NULL
        if ( file.exists ( fn) ) load( fn)
        return( B )
      }

      preds = spacetime.db( p=p, DS="predictions" )
      nr = p$nplons
      nc = p$nplats

      BP = expand.grid( plon=p$plons, plat=p$plats ) # coords of full prediction area
      attr( BP, "out.attrs") = NULL
      BP$substrate.predictionMean = preds[,2]
      BP$substrate.predictionSD   = preds[,3]

      # remove land
      oc = landmask( db="worldHires", regions=c("Canada", "US"), return.value="land", tag="predictions" )
      BP$substrate.predictionMean[oc] = NA
      BP$substrate.predictionSD[oc]   = NA

      rm(preds); gc()

      # merge into statistics
      BS = spacetime.db( p=p, DS="statistics" )
      B = cbind( BP, BS )
      names(B) = c( names(BP), "substrate.rangeMode", "substrate.rangeSD", "substrate.spatialSD", "substrate.observationSD" )

      save( B, file=fn, compress=TRUE)
      return(fn)

      if (0) {
        levelplot( log(substrate.predictionMean) ~ plon + plat, B, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        levelplot( log(substrate.rangeMode) ~ plon + plat, B, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
        levelplot( substrate.rangeSD ~ plon + plat, B, aspect="iso", labels=FALSE, pretty=TRUE, xlab=NULL,ylab=NULL,scales=list(draw=FALSE) )
      }
    }


    # ------------


     if (DS=="reset.bigmemory.objects" ) {
        # note::depends upon bathymetry
        substrate.db ( p=p, DS="substrate.spacetime.inputs.data.redo" )
        substrate.db( p=p, DS="substrate.spacetime.inputs.prediction.redo" )

        # reset input data objects
        spacetime.db( p=p, DS="bigmemory.inputs.data", B=substrate.db( p=p, DS="substrate.spacetime.inputs.data" ) )
        spacetime.db( p=p, DS="bigmemory.inputs.prediction", B=substrate.db( p=p, DS="substrate.spacetime.inputs.prediction" ) ) # note this is the same as inputs

        # reset bigmemory output data objects  (e.g., if you are restarting)
        spacetime.db( p=p, DS="predictions.bigmemory.initialize" )
        spacetime.db( p=p, DS="statistics.bigmemory.initialize" )

        # define boundary polygon for data to drop land etc vary parameters until it matches data ...
        p$mesh.boundary.resolution = 120 # discretization
        p$mesh.boundary.convex = -0.03  # curavature of boundary
        spacetime.db( p, DS="boundary.redo" )
    }


    # ------------


    if ( DS %in% c( "spde", "spde.redo" ) ) {
      #// substrate.db( DS="spde" .. ) returns the spatial interpolations from inla

      p$dist.mwin = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      p$upsampling = c( 1.1, 1.2, 1.5, 2 )  # local block search fractions
      p$downsampling = c( 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25, 0.2 ) # local block search fractions  -- need to adjust based upon data density
      p$sbbox = spacetime.db( p=p, DS="statistics.box" ) # bounding box and resoltuoin of output statistics defaults to 1 km X 1 km
      # p$variables = list( Y="substrate", X=c("z", "dZ", "ddZ", "Z.rangeMode" ), LOCS=c("plon", "plat") )
      p$variables = list( Y="substrate", X=c("z", "dZ" ), LOCS=c("plon", "plat") )

      p$spacetime.link = function( X ) { log(X ) }  ## data range is from -100 to 5467 m .. 1000 shifts all to positive valued by one -order of magnitude
      p$spacetime.invlink = function( X ) { exp(X) }

      p$dist.max = 100 # length scale (km) of local analysis .. for acceptance into the local analysis/model
      p$dist.min = 75 # lower than this .. subsampling occurs
      p$dist.pred = 0.95 # % of dist.max where **predictions** are retained (to remove edge effects)
      p$n.min = 30 # n.min/n.max changes with resolution: at p$pres=0.25, p$dist.max=25: the max count expected is 40000
      p$n.max = 8000 # numerical time/memory constraint -- anything larger takes too much time

      p$expected.range = 50 #+units=km km , with dependent var on log scale
      p$expected.sigma = 1e-1  # spatial standard deviation (partial sill) .. on log scale

      p$spatial.field.name = "spatial.field"  # name used in formula to index the spatal random field

      p$modelformula = formula( substrate ~ -1 + intercept
        + f( inla.group(log(z+1000) ), model="rw2")
        + f( inla.group(log(dZ+0.01)), model="rw2")
        # + f( inla.group( log(ddZ+0.01) ), model="rw2")
        # + f( inla.group( log(Z.rangeMode+0.01)), model="rw2" )
        + f( spatial.field, model=SPDE ) )

      p$spacetime.family = "gaussian"
      p$spacetime.outputs = c( "predictions.direct", "statistics" ) # "random.field", etc.
      p$statsvars = c("range", "range.sd", "spatial.error", "observation.error")

      if (p$substrate.bigmemory.reset) substrate.db( p=p, DS="reset.bigmemory.objects" )

      # run the beast .. warning this will take a very long time! (weeks)
      sS = spacetime.db( p, DS="statistics.bigmemory.status" )
      sS$n.incomplete / ( sS$n.problematic + sS$n.incomplete + sS$n.complete)

      if (0) {
        # for checking status of outputs during parallel runs:
        bathymetry.figures( DS="statistics", p=p )
        bathymetry.figures( DS="predictions", p=p )
        bathymetry.figures( DS="predictions.error", p=p )

        S = bigmemory::attach.big.matrix(p$descriptorfile.S, path=p$tmp.datadir)  # statistical outputs
        hist(S[,1] )
        o = which( S[,1] > 600 )
        S[o,] = NA
        S[sS$problematic,] = NA
        o = which( S[,1] < 10 )
        S[o,] = NA
      }

      p = make.list( list( jj=sample( sS$incomplete ) ), Y=p ) # random order helps use all cpus
      parallel.run( spacetime.interpolate.inla.local, p=p ) # no more GMT dependency! :)
      # spacetime.interpolate.inla.local( p=p, debugrun=TRUE )  # if testing serial process

      # save to file
      spacetime.db( p=p, DS="predictions.redo" )
      spacetime.db( p=p, DS="statistics.redo" )  # this also rescales results to the full domain

      ### bring together stats and predictions and any other required computations
      substrate.db( p=p, DS="substrate.spacetime.finalize.redo" )

      # clean up bigmemory files
      spacetime.db( p=p, DS="bigmemory.cleanup" )

    }

    # ------------

    if ( DS %in% c( "covariance.spatial", "covariance.spatial.redo" ) ) {
      #// substrate.db( DS="covariance" .. ) returns the spatial covariance estimates
      Z = NULL
      rootdir = file.path( p$project.root, "spacetime" )
      fn.results =  file.path( rootdir, paste( "spatial", "covariance", p$spatial.domain, "rdata", sep=".") )

      if  (DS %in% c("covariance.spatial"))  {
        stats = NULL
        if (file.exists( fn.results) ) load( fn.results )
        return(stats)
      }

      p$variogram.engine = "gstat"  # "geoR" seg faults frequently ..
      p$dist.max = 150 # length scale (km) of local analysis .. for acceptance into the local analysis/model
      p$dist.min = 100 # length scale (km) of local analysis .. beyond which subsampling occurs
      p$dist.mwin = 5 # resolution (km) of data aggregation (i.e. generation of the ** statistics ** )
      p$n.min = 30 # n.min/n.max changes with resolution: at p$pres=0.25, p$dist.max=25: the max count expected is 40000
      p$n.max = 10000 # numerical time/memory constraint -- anything larger takes too much time
      p$upsampling = c( 1.1, 1.2, 1.5, 2 )  # local block search fractions
      p$downsampling = c( 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.25 ) # local block search fractions  -- need to adjust based upon data density
      p$sbbox = spacetime.db( p=p, DS="statistics.box" ) # bounding box and resoltuoin of output statistics defaults to 1 km X 1 km
      p$variables = list( Y="substrate", LOCS=c("plon", "plat") )
      p$spacetime.link = function( X ) { log(X)  }
      p$spacetime.invlink = function( X ) { exp(X)  }
      p$statsvars = c("varTot", "varSpatial", "varObs", "range", "phi", "kappa" )

      # set up the data and problem using bigmemory data objects
      p = spacetime.db( p=p, DS="bigmemory.filenames" )
      print( paste( "Temporary files are being created at:", p$tmp.datadir ) )

      if (p$substrate.bigmemory.reset) substrate.db( p=p, DS="reset.bigmemory.objects" )

      if (0) {
        # to reset results manually .. just a template
        # p = spacetime.db( p=p, DS="bigmemory.filenames" )
        S = bigmemory::attach.big.matrix(p$descriptorfile.S, path=p$tmp.datadir)  # statistical outputs
        hist(S[,1] )
        o = which( S[,1] > xxx )
        S[o,] = NA
        S[sS$problematic,] = NA
        o = which( S[,1] < yyy )
        S[o,] = NA
        # etc ...
      }

      sS = spacetime.db( p, DS="statistics.bigmemory.status" )
      sS$n.incomplete / ( sS$n.problematic + sS$n.incomplete + sS$n.complete)

      p = make.list( list( jj=sample( sS$incomplete ) ), Y=p ) # random order helps use all cpus
      parallel.run( spacetime.covariance.spatial, p=p ) # no more GMT dependency! :)
      # spacetime.covariance.spatial( p=p )  # if testing serial process

      print( paste( "Results are being saved to:", fn.results ) )
      stats = bigmemory::attach.big.matrix(p$descriptorfile.S, path=p$tmp.datadir)  # statistical outputs
      stats = as.data.frame( stats[] )
      save(stats, file=fn.results, compress=TRUE )

      print( paste( "Temporary files are being deleted at:", p$tmp.datadir, "tmp" ) )
      spacetime.db( p=p, DS="bigmemory.cleanup" )

      return( fn.results )
    }

    # -------------

    if ( DS %in% c( "complete", "complete.redo" ) ) {
     #// substrate.db( DS="complete" .. ) returns the final form of the substrate data after
     #// regridding and selection to area of interest as specificied by girds.new=c("SSE", etc)
      Z = NULL

      if ( DS %in% c("complete") ) {

        domain = NULL
        if ( is.null(domain)) {
          if ( exists("spatial.domain", p)) {
            domain = p$spatial.domain
          } else if ( exists( "grids.new", p) )  { # over-rides p$spatial domain
            if( length( p$grids.new )== 1 ) {
              domain = p$grids.new
        } } }
        fn = file.path( project.datadirectory("substrate", "interpolated"),
          paste( "substrate", "complete", domain, "rdata", sep=".") )
        if ( file.exists ( fn) ) load( fn)
        if ( return.format == "dataframe" ) { ## default
          Z = as( brick(Z), "SpatialPointsDataFrame" )
          Z = as.data.frame(Z)
          u = names(Z)
          names(Z)[ which( u=="x") ] = "plon"
          names(Z)[ which( u=="y") ] = "plat"
          return( Z )
        }
        if ( return.format %in% c("list") ) return( Z  )
      }

      p0 = p  # the originating parameters
      Z0 = substrate.db( p=p0, DS="substrate.spacetime.finalize" )
      coordinates( Z0 ) = ~ plon + plat
      crs(Z0) = crs( p0$interal.crs )
      above.sealevel = which( Z0$z < 0 ) # depth values < 0 are above
      if (length(above.sealevel)>0) Z0[ above.sealevel ] = NA

      Z = list()

      grids = unique( c( p$spatial.domain, p$grids.new ))

      for (gr in grids ) {
        p1 = spatial.parameters( type=gr )
        for (vn in names(Z0)) {
          Z[[vn]] = projectRaster(
            from =rasterize( Z0, spatial.parameters.to.raster(p0), field=vn, fun=mean),
            to   =spatial.parameters.to.raster( p1) )
        }
        fn = file.path( project.datadirectory("substrate", "interpolated"),
          paste( "substrate", "complete", p1$spatial.domain, "rdata", sep=".") )
        save (Z, file=fn, compress=TRUE)
        print(fn)
      }
      return ( "Completed subsets" )
    }

  }


