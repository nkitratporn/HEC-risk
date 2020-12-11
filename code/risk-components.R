#[ status: on-going]
#[ note: ]
#@ This code is for preping and calculating risk components: hazard, exposure, and vulnerability
#@ The calculation was perofrmed for baseline and all future scenarios
#@ hazard has 5 scenarios: baseline, A1 (RCP45-SSP2-BZ), A2 (RCP85-SSP5-BZ), B1 (RCP45-SSP2-noBZ), B2 (RCP85-SSP5-noBZ)
#@ exposure has 3 scenarios: baseline, 1, 2
#@ vulnerability has 3 scenarios: baseline, 1, 2
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
  #devtools::install_github('biomodhub/biomod2')
  #  library(biomod2)   # Ensemble SDM package 
  library(ENMeval)
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
  # The palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  ## SET DIR ##
  setwd('D:/Academic/Research/0_maximus-materials/analysis/regional-land-elephant-trend/')
  
  ## TH polygon
  th <- readOGR('./data/shapefile/thailand.shp')
  
  ref.rs <- raster('./data/raster/landscape/historical_2000-2019/crop/_EVI_min.tif')
  ref.rs[!is.na(ref.rs)] <- 0
}

## FUNCTIONS
raster.transformation <- function(x, trans = "norm", smin=0, smax=255) {
  slog <- function(x) { ifelse(abs(x) <= 1, 0, sign(x)*log10(abs(x)))}
  rmin <- raster::cellStats(x, stat = "min", na.rm = TRUE)
  rmax <- raster::cellStats(x, stat = "max", na.rm = TRUE)
  rmean <- raster::cellStats(x, stat = "mean", na.rm = TRUE)
  rsd <- raster::cellStats(x, stat = "sd", na.rm = TRUE)
  
  if( trans == "slog" && rmin < 0) {
    stop(" Minimum value < 0, cannot log transform")
  }
  if( trans == "nl" && rmin < 0) {
    stop(" Minimum value < 0, cannot log transform")
  }
  if( trans == "sr" && rmin < 0) {
    stop(" Minimum value < 0, cannot log transform")
  }
  if( trans == "norm" && rmin < 0) {
    print(" Min value < 0, running row standardization instead")
    return( x / rmax )
  }
  
  if( trans == "norm") {
    message("applying normalization transformation", "\n")
    return( ( x - rmin ) / ( rmax - rmin ) )
  } else if ( trans == "rstd") {
    message("applying row-standardization transformation", "\n")
    return( x / rmax )
  } else if ( trans == "std") {
    message("applying standardization transformation", "\n")
    return( (x - rmean) / rsd )
  } else if ( trans == "stretch") {
    message("applying stretch transformation", "\n")
    return( (x - rmin) * smax / (rmax - rmin) + smin )
  } else if ( trans == "nl") {
    message("applying log transformation", "\n")
    return(  raster::calc(x, fun=log) )
  } else if ( trans == "slog") {
    message("applying singned-log10 transformation", "\n")
    return(raster::calc(x, fun=slog) )
  } else if ( trans == "sr") {
    message("applying sqare-root transformation", "\n")
    return(  raster::calc(x, fun=sqrt) )		  
  } else {
    stop("Not a valid transformation") 
  }		
}
rs.prep <- function(raster,ref,country){
  rs.resample <- resample(raster, ref, method='bilinear')
  rs.resample.mask <- raster::mask(raster::crop(raster,extent(ref)),country)
  rs.resmaple.mask
}
index.comp <- function(index, name, ref){
  index.df <- as.data.frame(index,xy=T)
  index.df.noNA <- index.df %>%
    drop_na()
  index.df.NA <- index.df[!complete.cases(index.df),]
  index.df.NA$geometric_mean <- NA
  index.geom_estimated <- ci_geom_gen(
    index.df.noNA[name],
    meth='EQUAL'
  )
  index.geom_estimated.df <- index.geom_estimated$ci_mean_geom_est
  index.df.noNA$geometric_mean <- index.geom_estimated.df
  # bind aggregation estimation to NA
  index.aggregate.df <- bind_rows(index.df.NA,index.df.noNA)
  
  # convert to raster
  last <- ncol(index.aggregate.df)
  index.aggregate <- rasterFromXYZ(index.aggregate.df[,c(1,2,last)])
  crs(index.aggregate) <- crs(ref)
  index.aggregate
}
index.class <- function(raster){
  raster[raster <= 0.2] <- 10                # very low
  raster[raster > 0.2 & raster <= 0.4] <- 20 # low
  raster[raster > 0.4 & raster <= 0.6] <- 30 # moderate
  raster[raster > 0.6 & raster <= 0.8] <- 40 # high
  raster[raster > 0.8 & raster <= 1.0] <- 50 # very high
  raster
}

#============= HAZARD ==============

# load all necessary indexes 1) probability of occurrence and 2) invert distance from known PA with elephants
# prepare dataset to ensure similar extent, resolution, projection
{
  ## 1. probability of occurrance
  # 1.1 climatic
  
  climate.result <- stack(
    stack('./AsianElephant/proj_Present/proj_Present_AsianElephant_ensemble.grd')[[3]],
    stack('./AsianElephant/proj_NearFuture_rcp45_m1/proj_NearFuture_rcp45_m1_AsianElephant_ensemble.grd')[[3]],
    stack('./AsianElephant/proj_NearFuture_rcp45_m2/proj_NearFuture_rcp45_m2_AsianElephant_ensemble.grd')[[3]],
    stack('./AsianElephant/proj_NearFuture_rcp45_m3/proj_NearFuture_rcp45_m3_AsianElephant_ensemble.grd')[[3]],
    stack('./AsianElephant/proj_NearFuture_rcp45_m4/proj_NearFuture_rcp45_m4_AsianElephant_ensemble.grd')[[3]],
    stack('./AsianElephant/proj_NearFuture_rcp45_m5/proj_NearFuture_rcp45_m5_AsianElephant_ensemble.grd')[[3]],
    stack('./AsianElephant/proj_NearFuture_rcp85_m1/proj_NearFuture_rcp85_m1_AsianElephant_ensemble.grd')[[3]],
    stack('./AsianElephant/proj_NearFuture_rcp85_m2/proj_NearFuture_rcp85_m2_AsianElephant_ensemble.grd')[[3]],
    stack('./AsianElephant/proj_NearFuture_rcp85_m3/proj_NearFuture_rcp85_m3_AsianElephant_ensemble.grd')[[3]],
    stack('./AsianElephant/proj_NearFuture_rcp85_m4/proj_NearFuture_rcp85_m4_AsianElephant_ensemble.grd')[[3]],
    stack('./AsianElephant/proj_NearFuture_rcp85_m5/proj_NearFuture_rcp85_m5_AsianElephant_ensemble.grd')[[3]]
  )
  
  names(climate.result) <- c(
    'Baseline',
    'RCP45_CESM1-BGC',
    'RCP45_MPI-ESM-MR',
    'RCP45_MIROC5',
    'RCP45_IPSL-CM5A-MR',
    'RCP45_CanESM2',
    'RCP85_CESM1-BGC',
    'RCP85_MPI-ESM-MR',
    'RCP85_MIROC5',
    'RCP85_IPSL-CM5A-MR',
    'RCP85_CanESM2'
  )
  
  climate.result <- climate.result/1000
  
  climate.result <- raster::mask(raster::crop(climate.result,extent(ref.rs)),th)
  
  climatic.prob.baseline <- climate.result$Baseline
  climatic.mean.baseline <- resample(climatic.prob.baseline, ref.rs, method='bilinear')
  climatic.mean.baseline <- raster::mask(raster::crop(climatic.mean.baseline,extent(ref.rs)),th)
  names(climatic.mean.baseline) <- c('baseline_suitability')
  
  climatic.prob.rcp45 <- climate.result[[2:6]]
  climatic.mean.rcp45 <- mean(climatic.prob.rcp45)
  climatic.mean.rcp45 <- resample(climatic.mean.rcp45, ref.rs, method='bilinear')
  climatic.mean.rcp45 <- raster::mask(raster::crop(climatic.mean.rcp45,extent(ref.rs)),th)
  names(climatic.mean.rcp45) <- c('RCP45_suitability')
  
  climatic.prob.rcp85 <- climate.result[[7:11]]
  climatic.mean.rcp85 <- mean(climatic.prob.rcp85)
  climatic.mean.rcp85 <- resample(climatic.mean.rcp85, ref.rs, method='bilinear')
  climatic.mean.rcp85 <- raster::mask(raster::crop(climatic.mean.rcp85,extent(ref.rs)),th)
  names(climatic.mean.rcp85) <- c('RCP85_suitability')
  
  # 1.2 landscape
  
  land.result <- stack(
    stack('./AsianElephant.land.noEVI/proj_Present_land_noEVI/proj_Present_land_noEVI_AsianElephant.land.noEVI_ensemble.grd')[[5]],
    stack('./AsianElephant.land.noEVI/proj_ssp2_buffer/proj_ssp2_buffer_AsianElephant.land.noEVI_ensemble.grd')[[5]],
    stack('./AsianElephant.land.noEVI/proj_ssp2_nobuffer/proj_ssp2_nobuffer_AsianElephant.land.noEVI_ensemble.grd')[[5]],
    stack('./AsianElephant.land.noEVI/proj_ssp5_buffer/proj_ssp5_buffer_AsianElephant.land.noEVI_ensemble.grd')[[5]],
    stack('./AsianElephant.land.noEVI/proj_ssp5_nobuffer/proj_ssp5_nobuffer_AsianElephant.land.noEVI_ensemble.grd')[[5]]
  )
  
  names(land.result) <- c('Baseline',
                          'RCP45_SSP2_BZ',
                          'RCP45_SSP2_noBZ',
                          'RCP85_SSP5_BZ',
                          'RCP85_SSP5_noBZ'
  )
  
  land.result <- land.result/1000
  
  landscape.prob <- raster::mask(raster::crop(land.result,extent(ref.rs)),th)
  
  ## 2. inverse distance from PAs
  invs.distance <- raster('./data/raster/hazard/elephant_dist_inv_norm_th.tif')
  invs.distance <- raster::mask(raster::crop(invs.distance,extent(ref.rs)),th)
  names(invs.distance) <- c('Dispersal_probability')
  
}

# set scenario for hazard with 5 possible scenarios
# (a) 1 baseline scenario
# (b) 2 climatic scenarios (average climatic results from 5 GCMs): RCP45=> A1/B1 and RCP85=> A2/B2
# (c) 4 landscape scenario: A1, A2, B1, B2
{
  hazard.base <- stack(climatic.mean.baseline,landscape.prob$Baseline,invs.distance)
  hazard.a1 <- stack(climatic.mean.rcp45,landscape.prob$RCP45_SSP2_BZ,invs.distance)
  hazard.a2 <- stack(climatic.mean.rcp85,landscape.prob$RCP85_SSP5_BZ,invs.distance)
  hazard.b1 <- stack(climatic.mean.rcp45,landscape.prob$RCP45_SSP2_noBZ,invs.distance)
  hazard.b2 <- stack(climatic.mean.rcp85,landscape.prob$RCP85_SSP5_noBZ,invs.distance)
}

# multicollinearity
{
  pairs(hazard.base)
  pairs(hazard.a1)
  pairs(hazard.a2)
  pairs(hazard.b1)
  pairs(hazard.b2)
}

# composite index
{
  hazard.base.i <- index.comp(hazard.base,names(hazard.base),ref.rs)
  hazard.a1.i <- index.comp(hazard.a1,names(hazard.a1),ref.rs)
  hazard.a2.i <- index.comp(hazard.a2,names(hazard.a2),ref.rs)
  hazard.b1.i <- index.comp(hazard.b1,names(hazard.b1),ref.rs)
  hazard.b2.i <- index.comp(hazard.b2,names(hazard.b2),ref.rs)
  
  hazard.index <- stack(hazard.base.i,
                        hazard.a1.i,
                        hazard.a2.i,
                        hazard.b1.i,
                        hazard.b2.i
  )
  
  names(hazard.index) <- c(
    'hazard_baseline',
    'hazard_RCP45_SSP2_BZ',
    'hazard_RCP85_SSP5_BZ',
    'hazard_RCP45_SSP2_noBZ',
    'hazard_RCP85_SSP5_noBZ'
  )
  
  #plot
  levelplot(hazard.index)
  ## save
  writeRaster(hazard.index,bylayer=T, names(hazard.index), format='GTiff', overwrite=T, './output/raster/')
}

# changes
{
  hazard.index <- stack('./output/raster/update/hazard_baseline.tif',
                        './output/raster/update/hazard_RCP45_SSP2_BZ.tif',
                        './output/raster/update/hazard_RCP85_SSP5_BZ.tif',
                        './output/raster/update/hazard_RCP45_SSP2_noBZ.tif',
                        './output/raster/update/hazard_RCP85_SSP5_noBZ.tif')
  
  
  hazard.diff <- hazard.index[[c(2:5)]] - hazard.index$hazard_baseline
  names(hazard.diff) <- c(
    'hazard_diff_RCP45_SSP2_BZ',
    'hazard_diff_RCP85_SSP5_BZ',
    'hazard_diff_RCP45_SSP2_noBZ',
    'hazard_diff_RCP85_SSP5_noBZ'
  )
  
  # plot
  levelplot(hazard.diff)
  # save
  writeRaster(hazard.diff, bylayer=T, names(hazard.diff), format='GTiff', overwrite=T,'./output/raster/')
}

## classify hazard
{
  hazard.class <- index.class(hazard.index)
  
  names(hazard.class) <- c(
    'hazard_class_baseline',
    'hazard_class_RCP45_SSP2_BZ',
    'hazard_class_RCP85_SSP5_BZ',
    'hazard_class_RCP45_SSP2_noBZ',
    'hazard_class_RCP85_SSP5_noBZ'
  )
  
  # plot
  levelplot(hazard.class,margin=F)
  # save
  writeRaster(hazard.class,bylayer=T, names(hazard.class), format='GTiff', overwrite=T,'./output/raster/')
}

#============= EXPOSURE ==============
# read population count data
{
  #rural baseline 
  popcount <- raster('./data/raster/exposure/BaseYear_1km/baseYr_rural_2000.tif')
  
  # resample to 500m resolution and mask to country boundary
  popcount_th_down <- rs.prep(popcount,ref.rs,th)
  popcount_th_down <- popcount_th_down/4
  names(popcount_th_down) <- c('rural_baseline')
  writeRaster(popcount_th_down, './data/raster/exposure/popcount_th_down.tif',format='GTiff',overwrite=T)
  
  ## SSP2
  pop.ssp2 <- raster('./data/raster/exposure/population_projection/ssp2_rural_2040.tif')
  names(pop.ssp2) <- c('rural_ssp2_2040')
  pop.ssp2 <- rs.prep(pop.ssp2,ref.rs,th)
  pop.ssp2 <- pop.ssp2/4
  
  ## SSP5
  pop.ssp5 <- raster('./data/raster/exposure/population_projection/ssp5_rural_2040.tif')
  names(pop.ssp5) <- c('rural_ssp5_2040')
  pop.ssp5 <- rs.prep(pop.ssp5,ref.rs,th)
  pop.ssp5 <- pop.ssp5/4
  
  # combine
  pop.stack <- stack(popcount_th_down,pop.ssp2,pop.ssp5)
  
  # plot
  levelplot(pop.stack)
  #save
  writeRaster(pop.stack, overwrite=T,'./data/raster/exposure/population_combine.grd')
}

# log10 transform to deal with skewness
{
  # ensure all cell is > 0
  pop.stack[pop.stack<0] <- 0
  pop.log10 <- raster.transformation(pop.stack,trans='slog')
  names(pop.log10) <- c('Ln.pop_baseline','Ln.pop_RCP45_SSP2','Ln.pop_RCP85_SSP5')
  
  # plot
  levelplot(pop.stack.log10)
  
  ## save
  writeRaster(pop.log10, overwrite=T,'./data/raster/exposure/population_combine_log10.grd')
}

## normalization of log population count => exposure index
{
  
  pop.log10.min = min(cellStats(pop.log10, "min"))
  pop.log10.max = max(cellStats(pop.log10, "max"))
  
  pop.log10.norm_minmax <- (pop.log10 - pop.log10.min) / (pop.log10.max - pop.log10.min)
  
  # plot
  levelplot(pop.stack.log10.norm_minmax)
  # save
  writeRaster(pop.log10.norm_minmax,'./data/raster/exposure/population_combine_log10_norm.grd')
  
  exposure.index <- pop.stack.log10.norm_minmax
  names(exposure.index) <- c('exposure_baseline',
                             'exposure_RCP45_SSP2_2040',
                             'exposure_RCP85_SSP5_2040')
  
  writeRaster(exposure.index, bylayer=T, names(exposure.index), overwrite=T, format='GTiff','./data/raster/exposure/')
}

# changes
{
  exposure.diff <- exposure.index[[2:3]] - exposure.index$exposure_baseline
  names(exposure.diff) <- c('exposure_diff_RCP45_SSP2_2040','exposure_diff_RCP45_SSP5_2040')
  
  # save
  writeRaster(exposure.diff, bylayer=T, names(exposure.diff), overwrite=T, format='GTiff','./output/raster/')
}

## Classification
{
  exposure.class <- index.class(exposure.index)
  names(exposure.class) <- c('exposure_class_baseline',
                             'exposure_class_RCP45_SSP2_2040',
                             'exposure_class_RCP85_SSP5_2040')
  # plot
  levelplot(exposure.class,margin=F)
  # save
  writeRaster(exposure.result, bylayer=T, names(exposure.result), overwriet=T, format='GTiff','./output/raster/')
  
} 

#============= VULNERABILITY =========
## data preparation
{
  # read in vector for socioeconomic variables
  vulnerability.vector <- readOGR('./data/shapefile/vulnerability_subindicators_joined2.shp')
  names(vulnerability.vector@data)
  
  monthly.income <- rasterize(vulnerability.vector, ref.rs, field='monthly_in', fun='last')  
  workforce.secondaryEd <- rasterize(vulnerability.vector, ref.rs, field='X.pop_secon', fun='last')
  access.internet <- rasterize(vulnerability.vector, ref.rs, field='X.household', fun='last')
  
  ## probability of drought
  probdrought <- stack('./output/raster/drought_probability.grd')
  names(probdrought)
  
  # resampling
  probdrought <- rs.prep(probdrought,ref.rs,th)
}

##min-max monthly income
{
  # invert function
  raster.invert <- function(x) {  
    if (!inherits(x, "RasterLayer")) stop("MUST BE RasterLayer OBJECT")
    rmax <- raster::cellStats(x, stat = "max", na.rm = TRUE, asSample = FALSE)
    rmin <- raster::cellStats(x, stat = "min", na.rm = TRUE, asSample = FALSE)
    return( ((x - rmax) * -1) + rmin )
  }  
  
  # inverse min-max normalization
  income.norm_minmax_invs <- raster.invert(income.norm_minmax)
  workforce.secondaryEd_invs <- raster.invert(workforce.secondaryEd/100) 
  access.internet_invs <- raster.invert(access.internet/100)
  
  # save
  writeRaster(income.norm_minmax, overwrite=T,'./data/raster/vulnerability/income_normalized_min-max_invs.tif')
  writeRaster(workforce.secondaryEd_invs, overwrite=T,'./data/raster/vulnerability/workforce_secondaryEd_invs.tif')
  writeRaster(access.internet_invs, overwrite=T,'./data/raster/vulnerability/access_internet_invs.tif')
}

# set 3 scenarios
# (a) static socioeconomic
# (b) drought probability (2 scenarios: RCP45-SSP2 and RCP85-SSP5)
{
  vul.baseline <- stack(probdrought$drought_baseline,income.norm_minmax,workforce.secondaryEd_invs,access.internet_invs)
  vul.1 <- stack(probdrought$drought_RCP45_SSP2,income.norm_minmax,workforce.secondaryEd_invs,access.internet_invs)
  vul.2 <- stack(probdrought$drought_RCP85_SSP5,income.norm_minmax,workforce.secondaryEd_invs,access.internet_invs)
}

# multicollinearity
{
  pairs(vul.baseline)
  pairs(vul.1)
  pairs(vul.2)
}

# composite
{
  vul.baseline.i <- index.comp(vul.baseline,names(vul.baseline),ref.rs)
  vul.1.i <- index.comp(vul.1,names(vul.1),ref.rs)
  vul.2.i <- index.comp(vul.2,names(vul.2),ref.rs)
  
  vul.index <- stack(vul.baseline, vul.rcp45_ssp2, vul.rcp85_ssp5)
  names(vul.index) <- c('vulnerability_baseline',
                        'vulnerability_RCP45_SSP2',
                        'vulnerability_RCP85_SSP5')
  
  # save
  writeRaster(vul.index, bylayer=T, names(vul.index),format='GTiff',overwrite=T,'./output/raster/')
}

## calculate changes
{
  vul.diff <- vul.index[[2:3]] - vul.index$vulnerability_baseline
  names(vul.diff) <- c('vulnerability_diff_RCP45_SSP2',
                       'vulnerability_diff_RCP85_SSP5')
  
  # plot
  levelplot(vul.diff)
  # save
  writeRaster(vul.diff, bylayer=T, names(vul.diff),format='GTiff',overwrite=T,'./output/raster/')
  
}
## classification
{
  vul.class <- index.class(vul.index)
  names(vul.class) <- c(
    'vul_class_baseline',
    'vul_class_RCP45_SSP2',
    'vul_class_RCP85_SSP5'
  )
  
  # save
  writeRaster(vul.class,bylayer=T, names(vul.class), format='GTiff', overwrite=T,'./output/raster/')
}

#============= VISUALIZATION =========
breaks <- seq(0, 1, by = 0.1)
index.cols <- colorRampPalette(rev(brewer.pal(5,'RdYlGn')))(length(breaks) - 1)

levelplot(hazard.index, layout=c(2,1), at=breaks, col.regions=index.cols)
levelplot(exposure.index, layout=c(2,1), at=breaks, col.regions=index.cols)
levelplot(vul.index, layout=c(2,1), at=breaks, col.regions=index.cols)