
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
  lbm( p=p, tasks=c("initiate", "globalmodel" ), DATA=DATA ) # 5 min
  

  # DATA='substrate.db( p=p, DS="lbm.inputs" )'
  # lbm( p=p, DATA=DATA, tasks=c("initiate", "globalmodel" ) )
  lbm( p=p, tasks=c( "stage1" ) ) # do not need other stages as data density is so high .. 8 hrs
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
log.substrate.grainsize ~ s(log(z), k = 3, bs = "ts") + s(log(dZ), 
    k = 3, bs = "ts") + s(log(ddZ), k = 3, bs = "ts")

Parametric coefficients:
             Estimate Std. Error t value Pr(>|t|)    
(Intercept) -0.786182   0.001686  -466.3   <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

Approximate significance of smooth terms:
              edf Ref.df        F p-value    
s(log(z))   2.000      2 114449.8  <2e-16 ***
s(log(dZ))  1.999      2   2614.3  <2e-16 ***
s(log(ddZ)) 1.999      2    746.5  <2e-16 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

R-sq.(adj) =  0.245   Deviance explained = 24.5%
GCV = 2.0298  Scale est. = 2.0298    n = 713953
---




