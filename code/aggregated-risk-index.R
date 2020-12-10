#[ status: on-going]
#[ note: ]
#@ This code is for the calculation of aggregated risk index
#@ 1. calculate geometric mean of 'hazard', 'exposure' and 'vulnerability' components 
#@ 2. classify the result [0,1] using equal interval from 'Very Low' to 'Very High'
#@ 3. calculate changes from baseline
#=====================================#
{
  #package management
  if(!require(devtools)){
    install.packages("devtools")
  }
  
  # geometric mean
  library(Compind)
  
  # data handling
  library(tidyverse)
  
  # modeling
  library(usdm)
  
  # spaital related
  library(raster)
  library(rgeos)
  library(rgdal)
  library(sp)
  library(maptools) # Vector data management (sp)
  library(maps)     # Easy access to basic map layers
  library(sf)
  
  # visualization
  library(rasterVis)
  library(gridExtra)
  library(ggplot2)
  library(corrplot)
  library(RColorBrewer)
}
#====================================#

## GLOBAL VARIABLES ##
{

  ## country polygon
  th <- readOGR('./data/shapefile/thailand.shp')
  
  ## reference raster
  ref.rs <- raster('./data/raster/landscape/historical_2000-2019/crop/_EVI_min.tif')
  ref.rs[!is.na(ref.rs)] <- 0
}

## FUNCTION ##
risk.result<- function(index, name){
  risk.index.df <- as.data.frame(index,xy=T)
  risk.index.df.noNA <- risk.index.df %>%
    drop_na()
  risk.index.df.NA <- risk.index.df[!complete.cases(risk.index.df),]
  risk.index.df.NA$geometric_mean <- NA
  risk.index.geom_estimated <- ci_geom_gen(
    risk.index.df.noNA[name],
    meth='EQUAL'
  )
  risk.index.geom_estimated.df <- risk.index.geom_estimated$ci_mean_geom_est
  risk.index.df.noNA$geometric_mean <- risk.index.geom_estimated.df
  # bind aggregation estimation to NA
  risk.index.aggregate.df <- bind_rows(risk.index.df.NA,risk.index.df.noNA)
  
  # convert to raster
  last <- ncol(risk.index.aggregate.df)
  risk.index.aggregate <- rasterFromXYZ(risk.index.aggregate.df[,c(1,2,last)])
  crs(risk.index.aggregate) <- crs(ref.rs)
  risk.index.aggregate
}
index.class <- function(raster){
  raster[raster <= 0.2] <- 10                # very low
  raster[raster > 0.2 & raster <= 0.4] <- 20 # low
  raster[raster > 0.4 & raster <= 0.6] <- 30 # moderate
  raster[raster > 0.6 & raster <= 0.8] <- 40 # high
  raster[raster > 0.8 & raster <= 1.0] <- 50 # very high
  raster
}


## READ RISK COMPONENT
{
  # hazard
  hazard.index <- stack('./output/raster/revised/hazard_baseline.tif',
                        './output/raster/revised/hazard_RCP45_SSP2_bufferzone.tif',
                        './output/raster/revised/hazard_RCP85_SSP5_bufferzone.tif',
                        './output/raster/revised/hazard_RCP45_SSP2_no_bufferzone.tif',
                        './output/raster/revised/hazard_RCP85_SSP5_no_bufferzone.tif'
                        )
  
  names(hazard.index) <- c('hazard_baseline',
                           'hazard_RCP45_SSP2_BZ',
                           'hazard_RCP85_SSP5_BZ',
                           'hazard_RCP45_SSP2_noBZ',
                           'hazard_RCP85_SSP5_noBZ'
                           )
  
  # exposure
  exposure.index <- stack('./output/raster/revised/exposure_baseline.tif',
                          './output/raster/revised/exposure_SSP2_2040.tif',
                          './output/raster/revised/exposure_SSP5_2040.tif')
  
  names(exposure.index) <- c('exposure_baseline',
                             'exposure_RCP45_SSP2',
                             'exposure_RCP85_SSP5')
  
  
  # vulnerability
  vul.index <- stack('./output/raster/revised/vulnerability_baseline.tif',
                     './output/raster/revised/vulnerability_RCP45_SSP2.tif',
                     './output/raster/revised/vulnerability_RCP85_SSP5.tif')
  
  names(vul.index) <- c('vulnerability_baseline',
                        'vulnerability_RCP45_SSP2',
                        'vulnerability_RCP85_SSP5')
}

## AGGREGATED RISK INDEX
{
  # create set of scenario with relavent components
  {
    risk.baseline <- stack(hazard.index$hazard_baseline,
                           exposure.index$exposure_baseline,
                           vul.index$vulnerability_baseline)
    levelplot(risk.baseline)
    
    risk.a1 <- stack(hazard.index$hazard_RCP45_SSP2_BZ,
                     exposure.index$exposure_RCP45_SSP2,
                     vul.index$vulnerability_RCP45_SSP2)
    levelplot(risk.a1)
    
    risk.a2 <- stack(hazard.index$hazard_RCP85_SSP5_BZ,
                     exposure.index$exposure_RCP85_SSP5,
                     vul.index$vulnerability_RCP85_SSP5)
    levelplot(risk.a2)
    
    risk.b1 <- stack(hazard.index$hazard_RCP45_SSP2_noBZ,
                     exposure.index$exposure_RCP45_SSP2,
                     vul.index$vulnerability_RCP45_SSP2)
    levelplot(risk.b1)
    
    risk.b2 <- stack(hazard.index$hazard_RCP85_SSP5_noBZ,
                     exposure.index$exposure_RCP85_SSP5,
                     vul.index$vulnerability_RCP85_SSP5)
    levelplot(risk.b2)
  }

  # apply geometric mean
  {
    risk.index.baseline <- risk.result(risk.baseline,names(risk.baseline))
    risk.index.a1 <- risk.result(risk.a1,names(risk.a1))
    risk.index.a2 <- risk.result(risk.a2,names(risk.a2))
    risk.index.b1 <- risk.result(risk.b1,names(risk.b1))
    risk.index.b2 <- risk.result(risk.b2,names(risk.b2))
    
    risk.index.combine <- stack(
      risk.index.baseline,
      risk.index.a1,
      risk.index.a2,
      risk.index.b1,
      risk.index.b2
    )
    
    names(risk.index.combine) <- c(
      'risk_baseline',
      'risk_RCP45_SSP2_BZ',
      'risk_RCP85_SSP5_BZ',
      'risk_RCP45_SSP2_noBZ',
      'risk_RCP85_SSP5_noBZ'
    )
    
    # plot
    levelplot(risk.index.combine)
    # save
    writeRaster(risk.index.combine, bylayer=T, names(risk.index.combine),format='GTiff',overwrite=T,'./output/raster/')
  }
}

## RISK CLASSIFICATION
{
  risk.class <- index.class(risk.index.combine)
  names(risk.class) <- c(
    'risk_class_baseline',
    'risk_class_RCP45_SSP2_BZ',
    'risk_class_RCP85_SSP5_BZ',
    'risk_class_RCP45_SSP2_noBZ',
    'risk_class_RCP85_SSP5_noBZ'
  )
  
  #plot
  levelplot(risk.class)
  #save
  writeRaster(risk.class, bylayer=T, names(risk.class), overwriet=T, format='GTiff','./output/raster/')
}

## CHANGE (from baseline)
{
  risk.change <- risk.index.combine[[2:5]] - risk.index.baseline
  names(risk.change) <- c(
    'risk_change_RCP45_SSP2_BZ',
    'risk_change_RCP85_SSP5_BZ',
    'risk_change_RCP45_SSP2_noBZ',
    'risk_change_RCP85_SSP5_noBZ'
  )
  
  # plot
  levelplot(risk.change)
  # save
  writeRaster(risk.change, bylayer=T, names(risk.change), overwriet=T, format='GTiff','./output/raster/')
}