
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
  p$storage.backend="bigmemory.ram"  # filebacked metods are still too slow ..
  p$lbm_local_modelengine = "krige" 

  p = bio.substrate::substrate.parameters( p=p, DS="lbm" )
  # p$clusters = rep("localhost",  detectCores() )
   
  DATA = 'substrate.db( p=p, DS="lbm.inputs" )'
  lbm( p=p, tasks=c("initiate", "globalmodel" ), DATA=DATA ) # 30 min
  

  # DATA='substrate.db( p=p, DS="lbm.inputs" )'
  # lbm( p=p, DATA=DATA, tasks=c("initiate", "globalmodel" ) )
  lbm( p=p, tasks=c( "stage1" ) ) # do not need other stages  .. 3hrs
  lbm( p=p, tasks=c( "save" ) )

  # to view progress in terminal:
  # watch -n 120 cat /home/jae/bio.data/bio.substrate/modelled/t/canada.east/lbm_current_status

  # to view maps from an external R session:
  # lbm(p=p, tasks="debug_pred_static_map", vindex=1)
  # lbm(p=p, tasks="debug_pred_static_log_map", vindex=1)
  # lbm(p=p, tasks="debug_pred_dynamic_map", vindex=1)
  # lbm(p=p, tasks="debug_stats_map", vindex=1)

 
  # as the interpolation process is so expensive, regrid based off the above run
  substrate.db( p=p, DS="complete.redo" )

  o = substrate.db( p=p, DS="complete" )
  b = bathymetry.db(p=p, DS="baseline")
  lattice::levelplot( o$log.substrate.grainsize ~ plon +plat, data=b, aspect="iso")



# to summarize just the global model
o = lbm_db( p=p, DS="global_model" )
summary(o)  
plot(o)
  
# Global model results:

Family: gaussian 
Link function: identity 

Formula:
log.substrate.grainsize ~ s(plon, k = 3, bs = "ts") + s(plat, 
    k = 3, bs = "ts") + s(plon, plat, k = 100, bs = "ts") + s(z, 
    k = 3, bs = "ts") + s(dZ, k = 3, bs = "ts") + s(ddZ, k = 3, 
    bs = "ts")

Parametric coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept) -0.786193   0.001055  -745.5   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
                edf Ref.df        F p-value    
s(plon)       2.000      2   9627.1  <2e-16 ***
s(plat)       2.000      2   3252.9  <2e-16 ***
s(plon,plat) 97.000     97   7990.3  <2e-16 ***
s(z)          2.000      2 112281.1  <2e-16 ***
s(dZ)         1.963      2   1659.5  <2e-16 ***
s(ddZ)        1.979      2    440.7  <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.705   Deviance explained = 70.5%
GCV = 0.79419  Scale est. = 0.79407   n = 713947


