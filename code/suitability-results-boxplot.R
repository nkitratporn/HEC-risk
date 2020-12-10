#[ status: on-going]
#[ note: ]
#@ This code is for visualizing the future suitability results under climate and landscape scenario
#@ Boxplot was used to compare results from different landscape scenarios and all GCMs under climatic scenarios 
#=====================================#
{
  #package management
  if(!require(devtools)){
    install.packages("devtools")
  }
  
  # data handling
  library(tidyverse)
  
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
  
  ## TH polygon
  th <- readOGR('./data/shapefile/thailand.shp')
  
  ## reference image
  ref.rs <- raster('./data/raster/landscape/historical_2000-2019/crop/_EVI_min.tif')
  ref.rs[!is.na(ref.rs)] <- 0
}

## LANDSCAPE
{
  # stack all landscape results (5 layers: 1 baseline + 4 future scenarios)
  land.result.stack <- stack(
    stack('./AsianElephant.land.noEVI/proj_Present_land_noEVI/proj_Present_land_noEVI_AsianElephant.land.noEVI_ensemble.grd')[[5]],
    stack('./AsianElephant.land.noEVI/proj_ssp2_buffer/proj_ssp2_buffer_AsianElephant.land.noEVI_ensemble.grd')[[5]],
    stack('./AsianElephant.land.noEVI/proj_ssp2_nobuffer/proj_ssp2_nobuffer_AsianElephant.land.noEVI_ensemble.grd')[[5]],
    stack('./AsianElephant.land.noEVI/proj_ssp5_buffer/proj_ssp5_buffer_AsianElephant.land.noEVI_ensemble.grd')[[5]],
    stack('./AsianElephant.land.noEVI/proj_ssp5_nobuffer/proj_ssp5_nobuffer_AsianElephant.land.noEVI_ensemble.grd')[[5]]
  )
  
  # rename raster layer
  names(land.result.stack) <- c('Baseline',
                                'RCP45-SSP2_BZ',
                                'RCP45-SSP2_noBZ',
                                'RCP85-SSP5_BZ',
                                'RCP85-SSP5_noBZ'
  )
  
  # change range from [0,1000] to [0,1]
  land.result.stack <- land.result.stack/1000
  
  # create boxplot
  boxplot(land.result.stack,notch=TRUE,col=c("grey","#E69F00", "#56B4E9", "#009E73", "#F0E442"),
                ylab = 'suitability probability',xaxt='n')
  
  #tick <- seq_along(tmp$names)
  #axis(1, at = tick, labels = F)
  #text(tick, par("usr")[3] - 0.2, tmp$names, srt = 45, xpd = T)
}

## CLIMATIC
{
  # stack all climatic suitability results [11 raster: 1 baseline and 10 future scenarios, 5 GCMs under 2 scenarios ]
  climate.result.base <- stack(
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
  
  # rename the raster layer
  names(climate.result.base) <- c(
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
  
  # change range from [0,1000] to [0,1]
  climate.result.stack <- climate.result.base/1000
  
  # crop and mask to the country boundary
  climate.result.stack <- raster::mask(raster::crop(climate.result.stack,extent(ref.rs)),th)
  
  # create boxplot from raster
  boxplot(climate.result.stack,notch=TRUE,
          col=c("grey","#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2"),
          ylab = 'suitability probability',xaxt='n',ylim = c(0, 1))
  #invisible(lapply(1:ncol(climate.result.stack),
  #                 function(i) segments(x0 = i - 0.4,
  #                                      y0 = mean(climate.result.stack[, i]),
  #                                      x1 = i + 0.4,
  #                                      y1 = mean(climate.result.stack[, i]),
  #                                      col = "red", lwd = 2)))
  
}