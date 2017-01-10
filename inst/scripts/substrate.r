
  ## NOTE:: substrate size is really only relevant for SSE/snowcrab domain right now as no
  ##        other data source has been found/identified
  ##        but working at the size of canada.east.highres for compatibility with bathymetry
  ## TODO:: add data collected by snow crab survey and any others for that matter

  p = bio.substrate::substrate.parameters()

  if ( basedata.redo ) {
    substrate.db ( DS="substrate.initial.redo" ) # bring in Kostelev's data ... stored as a SpatialGridDataFrame
		substrate.db ( DS="lonlat.highres.redo" ) # in future .. additional data would be added here
  }


  

  p = bio.substrate::substrate.parameters() # reset to defaults
  p$lbm_local_modelengine = "krige" 
  p$storage.backend="bigmemory.ram"  # filebacked metods are still too slow ..
  p = bio.substrate::substrate.parameters( p=p, DS="lbm" )
  # p$clusters = rep("localhost",  detectCores() )
   
  DATA = 'substrate.db( p=p, DS="lbm.inputs" )'
  p = lbm( p=p, tasks=c("initiate", "globalmodel" ), DATA=DATA )
  
  if (0) {
    # to summarize just the global model
    o = lbm_db( p=p, DS="global_model" )
    summary(o)  
  
# Global model results: .. dZ a lbm_db( p=p, DS="global_model" and ddZ not informative (at a global level) .. drop

# Family: gaussian 
# Link function: log 

# Formula:
# log.substrate.grainsize ~ s(plon, k = 3, bs = "ts") + s(plat, 
#     k = 3, bs = "ts") + s(plon, plat, k = 100, bs = "ts") + s(z, 
#     k = 3, bs = "ts") + s(dZ, k = 3, bs = "ts") + s(ddZ, k = 3, 
#     bs = "ts")

# Parametric coefficients:
#             Estimate Std. Error t value Pr(>|t|)    
# (Intercept)  -68.885      4.323  -15.94   <2e-16 ***
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Approximate significance of smooth terms:
#                    edf Ref.df        F  p-value    
# s(plon)      7.559e-06      2    0.000 4.54e-07 ***
# s(plat)      1.975e+00      2    8.217 0.000229 ***
# s(plon,plat) 9.687e+01     97  164.805  < 2e-16 ***
# s(z)         1.999e+00      2 1303.932  < 2e-16 ***
# s(dZ)        3.179e-04      2    0.000 0.077859 .  
# s(ddZ)       1.662e-04      2    0.000 0.352290    
# ---
# Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# R-sq.(adj) =  0.248   Deviance explained = -12.1%
# GCV = 3.0141  Scale est. = 3.0137    n = 713974

  }

  # DATA='substrate.db( p=p, DS="lbm.inputs" )'
  # p = lbm( p=p, DATA=DATA, tasks=c("initiate", "globalmodel" ) )
  p = lbm( p=p, tasks=c( "stage1" ) ) # do not need other stages  
  p = lbm( p=p, tasks=c( "save" ) )

 
  # as the interpolation process is so expensive, regrid based off the above run
  substrate.db( p=p, DS="complete.redo" )

  



