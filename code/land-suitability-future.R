#[ status: on-going]
#[ note: ]
#@ This code is for projecting future suitability under landscape different landscape conditions
#@ The baseline model was previously calibrated. Hence, it can be load and use
#@ Relavent variables are read and prepared
#=====================================#
{
  #package management
  if(!require(devtools)){
    install.packages("devtools")
  }
  
  # data handling
  library(tidyverse)
  
  # modeling
  #devtools::install_github('biomodhub/biomod2')
  library(biomod2)   # Ensemble SDM package 
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
  
  ## TH polygon
  th <- readOGR('./data/shapefile/thailand.shp')
  
  ## reference image
  ref.rs <- raster('./data/raster/landscape/historical_2000-2019/crop/_EVI_min.tif')
  ref.rs[!is.na(ref.rs)] <- 0
}

## SUITABILITY PROJECTION FOR FUTURE LANDSCAPE
{
  # check the sampling in landscape model
  land.sample <- get(load('./AsianElephant.land.noEVI/.BIOMOD_DATA/land_noEVI/formated.input.data'))
  plot(land.sample)
  
  # load calibrated model
  land_out_esm <- get(load('./AsianElephant.land.noEVI/AsianElephant.land.noEVI.land_noEVIensemble.models.out'))
  land_out <- get(load('./AsianElephant.land.noEVI/AsianElephant.land.noEVI.land_noEVI.models.out'))
  
  # load variable for each scenario
  {
    # static variables
    {
      static <- stack(paste0(
        './data/raster/landscape/historical_2000-2019/lc_2016_unmixed/pre-process/',
        list.files('./data/raster/landscape/historical_2000-2019/lc_2016_unmixed/pre-process/',
                   pattern='*.tif$')))
      
      names(static) <- c(
        'climatic_suitability', 'cover_percent_1500m','cover_percent_3000m',
        'cover_percent_6000m','crop_dist','EVI_amplitude',
        'EVI_min','food_percent_1500m','food_percent_3000m',
        'food_percent_6000m','forest_dist','plantation_dist',
        'slope', 'transport_dist', 'TRI',
        'urban_dist', 'water_dist'
      )
      static <- static[[c('transport_dist','TRI','water_dist')]]
      names(static)
    }
    
    ## A1 rcp45-ssp2 buffer
    {
      ## load and prep variables
      ssp2.buffer <- stack(paste0('./data/raster/landscape/future_projection/lc_2045_variables/',
                                  list.files('./data/raster/landscape/future_projection/lc_2045_variables/',
                                             pattern='*ssp2_buffer.tif$')))
      names(ssp2.buffer)
      names(ssp2.buffer) <- c('cover_percent_1500m','cover_percent_3000m','cover_percent_6000m',
                              'crop_dist','food_percent_1500m','food_percent_3000m',
                              'food_percent_6000m', 'forest_dist','plantation_dist',
                              'urban_dist')
      ssp2.buffer <- ssp2.buffer[[c(
        'cover_percent_6000m','crop_dist','food_percent_6000m', 'forest_dist','plantation_dist',
        'urban_dist'
      )]]
      ssp2.buffer <- raster::mask(raster::crop(ssp2.buffer,extent(ref.rs)),th)
      ssp2.buffer.vars <- stack(static,ssp2.buffer)
      names(ssp2.buffer.vars)
      
      ## project by model
      land.prj.ssp2_buffer <- BIOMOD_Projection(
        modeling.output = land_out,
        new.env = ssp2.buffer.vars,
        proj.name = 'ssp2_buffer',
        selected.models = 'all',
        binary.meth = 'TSS',
        compress = 'gzip',
        clamping.mask = T,
        output.format = '.grd',
        do.stack=T
      )
      
      ## ensemble projection
      land.prjESM.ssp2_buffer <- BIOMOD_EnsembleForecasting(
        EM.output = land_out_esm,
        projection.output = land.prj.ssp2_buffer,
        binary.meth = 'TSS',
        output.format = '.grd',
        do.stack=T
      )
    }
    
    ## A2 rcp85-ssp5 buffer
    {
      ## load and prep variables
      ssp5.buffer <- stack(paste0('./data/raster/landscape/future_projection/lc_2045_variables/',
                                  list.files('./data/raster/landscape/future_projection/lc_2045_variables/',
                                             pattern='*ssp5_buffer.tif$')))
      names(ssp5.buffer) <- c('cover_percent_1500m','cover_percent_3000m','cover_percent_6000m',
                              'crop_dist','food_percent_1500m','food_percent_3000m',
                              'food_percent_6000m', 'forest_dist','plantation_dist',
                              'urban_dist')
      ssp5.buffer <- ssp5.buffer[[c(
        'cover_percent_6000m','crop_dist','food_percent_6000m', 'forest_dist','plantation_dist',
        'urban_dist'
      )]]
      ssp5.buffer <- raster::mask(raster::crop(ssp5.buffer,extent(ref.rs)),th)
      ssp5.buffer.vars <- stack(static,ssp5.buffer)
      names(ssp5.buffer.vars)
      
      ## project by model
      land.prj.ssp5_buffer <- BIOMOD_Projection(
        modeling.output = land_out,
        new.env = ssp5.buffer.vars,
        proj.name = 'ssp5_buffer',
        selected.models = 'all',
        binary.meth = 'TSS',
        compress = 'gzip',
        clamping.mask = T,
        output.format = '.grd',
        do.stack=T
      )
      
      ## ensemble projection
      land.prjESM.ssp5_buffer <- BIOMOD_EnsembleForecasting(
        EM.output = land_out_esm,
        projection.output = land.prj.ssp5_buffer,
        binary.meth = 'TSS',
        output.format = '.grd',
        do.stack=T
      )
    }
    
    ## B1 rcp45-ssp2 no buffer
    {
      ## load and prep variables
      ssp2.nobuffer <- stack(paste0('./data/raster/landscape/future_projection/lc_2045_variables/',
                                    list.files('./data/raster/landscape/future_projection/lc_2045_variables/',
                                               pattern='*ssp2_nobuffer.tif$')))
      names(ssp2.nobuffer) <- c('cover_percent_1500m','cover_percent_3000m','cover_percent_6000m',
                                'crop_dist','food_percent_1500m','food_percent_3000m',
                                'food_percent_6000m', 'forest_dist','plantation_dist',
                                'urban_dist')
      ssp2.nobuffer <- ssp2.nobuffer[[c(
        'cover_percent_6000m','crop_dist','food_percent_6000m', 'forest_dist','plantation_dist',
        'urban_dist'
      )]]
      ssp2.nobuffer <- raster::mask(raster::crop(ssp2.nobuffer,extent(ref.rs)),th)
      ssp2.nobuffer.vars <- stack(static,ssp2.nobuffer)
      names(ssp2.nobuffer.vars)
      
      ## project by model
      land.prj.ssp2_nobuffer <- BIOMOD_Projection(
        modeling.output = land_out,
        new.env = ssp2.nobuffer.vars,
        proj.name = 'ssp2_nobuffer',
        selected.models = 'all',
        binary.meth = 'TSS',
        compress = 'gzip',
        clamping.mask = T,
        output.format = '.grd',
        do.stack=T
      )
      
      ## ensemble projection
      land.prjESM.ssp2_nobuffer <- BIOMOD_EnsembleForecasting(
        EM.output = land_out_esm,
        projection.output = land.prj.ssp2_nobuffer,
        binary.meth = 'TSS',
        output.format = '.grd',
        do.stack=T
      )
    }
    
    ## B2 rcp85-ssp5 no buffer
    {
      ## load and prep variables
      ssp5.nobuffer <- stack(paste0('./data/raster/landscape/future_projection/lc_2045_variables/',
                                    list.files('./data/raster/landscape/future_projection/lc_2045_variables/',
                                               pattern='*ssp5_nobuffer.tif$')))
      names(ssp5.nobuffer) <- c('cover_percent_1500m','cover_percent_3000m','cover_percent_6000m',
                                'crop_dist','food_percent_1500m','food_percent_3000m',
                                'food_percent_6000m', 'forest_dist','plantation_dist',
                                'urban_dist')
      ssp5.nobuffer <- ssp5.nobuffer[[c(
        'cover_percent_6000m','crop_dist','food_percent_6000m', 'forest_dist','plantation_dist',
        'urban_dist'
      )]]
      ssp5.nobuffer <- raster::mask(raster::crop(ssp5.nobuffer,extent(ref.rs)),th)
      ssp5.nobuffer.vars <- stack(static,ssp5.nobuffer)
      names(ssp5.nobuffer.vars)
      
      ## project by model
      land.prj.ssp5_nobuffer <- BIOMOD_Projection(
        modeling.output = land_out,
        new.env = ssp5.nobuffer.vars,
        proj.name = 'ssp5_nobuffer',
        selected.models = 'all',
        binary.meth = 'TSS',
        compress = 'gzip',
        clamping.mask = T,
        output.format = '.grd',
        do.stack=T
      )
      
      ## ensemble projection
      land.prjESM.ssp5_nobuffer <- BIOMOD_EnsembleForecasting(
        EM.output = land_out_esm,
        projection.output = land.prj.ssp5_nobuffer,
        binary.meth = 'TSS',
        output.format = '.grd',
        do.stack=T
      )
    }
  }
}
