#[ status: on-going]
#[ note: ]
#@ This code performs comparison of PR, TASMIN, and TASMAX from ERA5 and NEX-GDDP to observed value obtained from weather stations
#@ First, the locations of weather stations were used to extract value from raster
#@ Second, combine all value under relavent cateogry (scenario, model, and dataset)
#@ Thrid, perform RMSE calculation
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
  
  ## parallel processing
  library(foreach)
  library(doParallel)
}
#=====================================#

## GLOBAL VARIABLES ##
{
  # The palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  ## SET DIR ##
  setwd('D:/Academic/Research/0_maximus-materials/analysis/regional-land-elephant-trend/')
  
  ## TH polygon
  th <- readOGR('./data/shapefile/thailand.shp')
  
  ## reference raster for cropping
  ref.rs <- raster('./data/raster/landscape/historical_2000-2019/crop/_EVI_min.tif')
  ref.rs[!is.na(ref.rs)] <- 0
  
  # first date of 2015
  first.date <- as.Date('20150101',format = '%Y%m%d')
}

## FUNCTIONS
rs.extract.df <- function(raster,points,scenario,model,dataset){
  # extract raster
  rs.extract <- extract(raster,points,cellnumbers=T,df=T)
  
  # data re-arrangement
  rs.df.m <- rs.extract %>% 
    reshape2::melt(id.vars = c('ID','cells'))
  
  if(dataset=='ERA5'){
    rs.df <- rs.df.m %>%
      mutate(date = as.Date(gsub('\\D+','\\1',variable), format = '%Y%m%d')) %>%
      dplyr::filter(!is.na(value)) %>%
      mutate(value = ifelse(value<0,0,value),
             scenario = scenario,
             model = model,
             dataset = dataset)
  }else{
    rs.df <- rs.df.m %>%
      mutate(date = first.date -1 + as.numeric(sapply(str_extract_all(variable, "\\d+"),tail,1))) %>%
      dplyr::filter(!is.na(value)) %>%
      mutate(value = ifelse(value<0,0,value),
             scenario = scenario,
             model = model,
             dataset = dataset)
  }
  rs.df
}

# weather station location
{
  weather_station <- readOGR('./data/climate_hist_2015/weather_station_location.shp')
  weather_station.df <- read_csv('./data/climate_hist_2015/weather_station_location.csv')
}

## extract value from raster [PR]
{
  ## ERA5
  {
    era5.pr <- stack(paste0('data/climate_hist_2015/era5_pr.tif'))
    era5.pr.df <- rs.extract.df(era5.pr,weather_station,'reanalysis','reanalysis','ERA5')
  }
  
  ## NEX-GDDP
  {
    # list model
    model <- c('CESM1-BGC','MPI-ESM-MR','MIROC5', 'IPSL-CM5A-MR','CanESM2')
    
    # rcp45
    rcp45.m0.pr <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m0_pr.tif'))
    rcp45.m1.pr <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m1_pr.tif'))
    rcp45.m2.pr <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m2_pr.tif'))
    rcp45.m3.pr <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m3_pr.tif'))
    rcp45.m4.pr <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m4_pr.tif'))
    
    # rcp85
    rcp85.m0.pr <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m0_pr.tif'))
    rcp85.m1.pr <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m1_pr.tif'))
    rcp85.m2.pr <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m2_pr.tif'))
    rcp85.m3.pr <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m3_pr.tif'))
    rcp85.m4.pr <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m4_pr.tif'))
    
    # extract
    rcp45.m0.pr.df <- rs.extract.df(rcp45.m0.pr,weather_station,'RCP45',model[1],'NEX-GDDP')
    rcp45.m1.pr.df <- rs.extract.df(rcp45.m1.pr,weather_station,'RCP45',model[2],'NEX-GDDP')
    rcp45.m2.pr.df <- rs.extract.df(rcp45.m2.pr,weather_station,'RCP45',model[3],'NEX-GDDP')
    rcp45.m3.pr.df <- rs.extract.df(rcp45.m3.pr,weather_station,'RCP45',model[4],'NEX-GDDP')
    rcp45.m4.pr.df <- rs.extract.df(rcp45.m4.pr,weather_station,'RCP45',model[5],'NEX-GDDP')
    
    rcp85.m0.pr.df <- rs.extract.df(rcp85.m0.pr,weather_station,'RCP85',model[1],'NEX-GDDP')
    rcp85.m1.pr.df <- rs.extract.df(rcp85.m1.pr,weather_station,'RCP85',model[2],'NEX-GDDP')
    rcp85.m2.pr.df <- rs.extract.df(rcp85.m2.pr,weather_station,'RCP85',model[3],'NEX-GDDP')
    rcp85.m3.pr.df <- rs.extract.df(rcp85.m3.pr,weather_station,'RCP85',model[4],'NEX-GDDP')
    rcp85.m4.pr.df <- rs.extract.df(rcp85.m4.pr,weather_station,'RCP85',model[5],'NEX-GDDP')
  }
  
  # combine df
  {
    pr.df.combine <- data.frame()
    pr.df.combine <- pr.df.combine %>% 
      bind_rows(era5.pr.df,
                rcp45.m0.pr.df,rcp45.m1.pr.df,rcp45.m2.pr.df,rcp45.m3.pr.df,rcp45.m4.pr.df,
                rcp85.m0.pr.df,rcp85.m1.pr.df,rcp85.m2.pr.df,rcp85.m3.pr.df,rcp85.m4.pr.df)
    pr.df.combine <- pr.df.combine %>% mutate(variable = 'pr')
    head(pr.df.combine)
  }
}

## extract value from raster [TASMIN]
{
  ## ERA5
  {
    era5.tasmin <- stack(paste0('data/climate_hist_2015/era5_tasmin.tif'))
    era5.tasmin.df <- rs.extract.df(era5.tasmin,weather_station,'reanalysis','reanalysis','ERA5')
  }
  
  ## NEX-GDDP
  {
    # list model
    model <- c('CESM1-BGC','MPI-ESM-MR','MIROC5', 'IPSL-CM5A-MR','CanESM2')
    
    # rcp45
    rcp45.m0.tasmin <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m0_tasmin.tif'))
    rcp45.m1.tasmin <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m1_tasmin.tif'))
    rcp45.m2.tasmin <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m2_tasmin.tif'))
    rcp45.m3.tasmin <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m3_tasmin.tif'))
    rcp45.m4.tasmin <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m4_tasmin.tif'))
    
    # rcp85
    rcp85.m0.tasmin <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m0_tasmin.tif'))
    rcp85.m1.tasmin <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m1_tasmin.tif'))
    rcp85.m2.tasmin <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m2_tasmin.tif'))
    rcp85.m3.tasmin <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m3_tasmin.tif'))
    rcp85.m4.tasmin <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m4_tasmin.tif'))
    
    # extract
    rcp45.m0.tasmin.df <- rs.extract.df(rcp45.m0.tasmin,weather_station,'RCP45',model[1],'NEX-GDDP')
    rcp45.m1.tasmin.df <- rs.extract.df(rcp45.m1.tasmin,weather_station,'RCP45',model[2],'NEX-GDDP')
    rcp45.m2.tasmin.df <- rs.extract.df(rcp45.m2.tasmin,weather_station,'RCP45',model[3],'NEX-GDDP')
    rcp45.m3.tasmin.df <- rs.extract.df(rcp45.m3.tasmin,weather_station,'RCP45',model[4],'NEX-GDDP')
    rcp45.m4.tasmin.df <- rs.extract.df(rcp45.m4.tasmin,weather_station,'RCP45',model[5],'NEX-GDDP')
    
    rcp85.m0.tasmin.df <- rs.extract.df(rcp85.m0.tasmin,weather_station,'RCP85',model[1],'NEX-GDDP')
    rcp85.m1.tasmin.df <- rs.extract.df(rcp85.m1.tasmin,weather_station,'RCP85',model[2],'NEX-GDDP')
    rcp85.m2.tasmin.df <- rs.extract.df(rcp85.m2.tasmin,weather_station,'RCP85',model[3],'NEX-GDDP')
    rcp85.m3.tasmin.df <- rs.extract.df(rcp85.m3.tasmin,weather_station,'RCP85',model[4],'NEX-GDDP')
    rcp85.m4.tasmin.df <- rs.extract.df(rcp85.m4.tasmin,weather_station,'RCP85',model[5],'NEX-GDDP')
  }
  
  # combine df
  {
    tasmin.df.combine <- data.frame()
    tasmin.df.combine <- tasmin.df.combine %>% 
      bind_rows(era5.tasmin.df,
                rcp45.m0.tasmin.df,rcp45.m1.tasmin.df,rcp45.m2.tasmin.df,rcp45.m3.tasmin.df,rcp45.m4.tasmin.df,
                rcp85.m0.tasmin.df,rcp85.m1.tasmin.df,rcp85.m2.tasmin.df,rcp85.m3.tasmin.df,rcp85.m4.tasmin.df)
    tasmin.df.combine <- tasmin.df.combine %>% mutate(variable = 'tasmin')
    head(tasmin.df.combine)
  }
}

## extract value from raster [TASMAX]
{
  ## ERA5
  {
    era5.tasmax <- stack(paste0('data/climate_hist_2015/era5_tasmax.tif'))
    era5.tasmax.df <- rs.extract.df(era5.tasmax,weather_station,'reanalysis','reanalysis','ERA5')
  }
  
  ## NEX-GDDP
  {
    # list model
    model <- c('CESM1-BGC','MPI-ESM-MR','MIROC5', 'IPSL-CM5A-MR','CanESM2')
    
    # rcp45
    rcp45.m0.tasmax <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m0_tasmax.tif'))
    rcp45.m1.tasmax <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m1_tasmax.tif'))
    rcp45.m2.tasmax <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m2_tasmax.tif'))
    rcp45.m3.tasmax <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m3_tasmax.tif'))
    rcp45.m4.tasmax <- stack(paste0('data/climate_hist_2015/nexgddp_rcp45_m4_tasmax.tif'))
    
    # rcp85
    rcp85.m0.tasmax <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m0_tasmax.tif'))
    rcp85.m1.tasmax <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m1_tasmax.tif'))
    rcp85.m2.tasmax <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m2_tasmax.tif'))
    rcp85.m3.tasmax <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m3_tasmax.tif'))
    rcp85.m4.tasmax <- stack(paste0('data/climate_hist_2015/nexgddp_rcp85_m4_tasmax.tif'))
    
    # extract
    rcp45.m0.tasmax.df <- rs.extract.df(rcp45.m0.tasmax,weather_station,'RCP45',model[1],'NEX-GDDP')
    rcp45.m1.tasmax.df <- rs.extract.df(rcp45.m1.tasmax,weather_station,'RCP45',model[2],'NEX-GDDP')
    rcp45.m2.tasmax.df <- rs.extract.df(rcp45.m2.tasmax,weather_station,'RCP45',model[3],'NEX-GDDP')
    rcp45.m3.tasmax.df <- rs.extract.df(rcp45.m3.tasmax,weather_station,'RCP45',model[4],'NEX-GDDP')
    rcp45.m4.tasmax.df <- rs.extract.df(rcp45.m4.tasmax,weather_station,'RCP45',model[5],'NEX-GDDP')
    
    rcp85.m0.tasmax.df <- rs.extract.df(rcp85.m0.tasmax,weather_station,'RCP85',model[1],'NEX-GDDP')
    rcp85.m1.tasmax.df <- rs.extract.df(rcp85.m1.tasmax,weather_station,'RCP85',model[2],'NEX-GDDP')
    rcp85.m2.tasmax.df <- rs.extract.df(rcp85.m2.tasmax,weather_station,'RCP85',model[3],'NEX-GDDP')
    rcp85.m3.tasmax.df <- rs.extract.df(rcp85.m3.tasmax,weather_station,'RCP85',model[4],'NEX-GDDP')
    rcp85.m4.tasmax.df <- rs.extract.df(rcp85.m4.tasmax,weather_station,'RCP85',model[5],'NEX-GDDP')
  }
  
  # combine df
  {
    tasmax.df.combine <- data.frame()
    tasmax.df.combine <- tasmax.df.combine %>% 
      bind_rows(era5.tasmax.df,
                rcp45.m0.tasmax.df,rcp45.m1.tasmax.df,rcp45.m2.tasmax.df,rcp45.m3.tasmax.df,rcp45.m4.tasmax.df,
                rcp85.m0.tasmax.df,rcp85.m1.tasmax.df,rcp85.m2.tasmax.df,rcp85.m3.tasmax.df,rcp85.m4.tasmax.df)
    tasmax.df.combine <- tasmax.df.combine %>% mutate(variable = 'tasmax')
    head(tasmax.df.combine)
  }
}

## combine and compare all dataset
{
  # combine climate dataset
  climate.dataset.df <- bind_rows(pr.df.combine,tasmin.df.combine,tasmax.df.combine)
  summary(climate.dataset.df)
  write_csv(climate.dataset.df,'./data/climate_hist_2015/climate_dataset_combine_2015.csv')
  
  # read observed dataset
  obs.df <- read_csv('./data/climate_hist_2015/observed_climatedata_th_2015.csv')
  obs.df <- obs.df %>%
    mutate(date = as.Date(date,format='%m/%d/%Y'))
  summary(obs.df)
  obs.df <- obs.df %>%
    mutate(value = ifelse(variable!='pr' & value == 0,NA,value))
  
  # combine observed and reanalysis
  climate.all.df <- bind_rows(obs.df,climate.dataset.df)
  write_csv(climate.all.df, './data/climate_hist_2015/climate_dataset_observation_combine_2015.csv')
  climate.all.df <- climate.all.df %>%
    mutate(scenario_model = paste0(scenario,'_',model))  %>%
    mutate(scenario_model = ifelse(scenario_model == 'observed_observed',
                                   '0bserved_TMD',
                                   ifelse(scenario_model == 'reanalysis_reanalysis',
                                          'ERA5_Reanalysis',scenario_model)))
  glimpse(climate.all.df)
}

## plot data for comparison using boxplot
{
  ggplot(climate.all.df %>% filter(variable=='pr' & value > 0 & ID=='41'),
         aes(x=scenario_model,y=value,fill=factor(scenario_model))) +
    geom_boxplot() +
    theme_bw() +
    labs(y = 'Daily Precipitation (mm)') +
    theme(legend.position = 'none',
          axis.title.x = element_blank(),axis.text.x = element_text(angle=90,vjust=0.5),
          text = element_text(size=14))
  
  ggplot(climate.all.df %>% filter(variable=='tasmin'),
         aes(x=scenario_model,y=value,fill=factor(scenario_model))) +
    geom_boxplot() +
    theme_bw() +
    labs(y='Min. Temperature (?C)') +
    theme(legend.position = 'none',
          axis.title.x = element_blank(),axis.text.x = element_text(angle=90,vjust=0.5),
          text = element_text(size=14))
  ggplot(climate.all.df %>% filter(variable=='tasmax'),
         aes(x=scenario_model,y=value,fill=factor(scenario_model))) +
    geom_boxplot() +
    theme_bw() +
    labs(y='Max. Temperature (?C)') +
    theme(legend.position = 'none',
          axis.title.x = element_blank(),axis.text.x = element_text(angle=90,vjust=0.5),
          text = element_text(size=14))
}

## calculate RMSE
{
  climate.all.reshape <- climate.all.df %>%
    mutate(var.name = paste0(scenario_model,'_',variable)) %>%
    dplyr::select(ID, date, variable, value, var.name) %>%
    pivot_wider(
      names_from = var.name,
      values_from = value
    )
  all_na <- function(x) any(!is.na(x))
  
  pr.all.reshape <- climate.all.reshape %>% filter(variable == 'pr') %>% select_if(all_na) %>% 
    drop_na() %>% as.data.frame()
  colnames(pr.all.reshape) <- c('ID','date','variable','Observed','ERA5', colnames(pr.all.reshape[6:15]))

  tasmin.all.reshape <- climate.all.reshape %>% filter(variable == 'tasmin') %>% select_if(all_na) %>% 
    drop_na() %>% as.data.frame()
  colnames(tasmin.all.reshape) <- c('ID','date','variable','Observed','ERA5', colnames(tasmin.all.reshape[6:15]))
  
  tasmax.all.reshape <- climate.all.reshape %>% filter(variable == 'tasmax') %>% select_if(all_na) %>% 
    drop_na() %>% as.data.frame()
  colnames(tasmax.all.reshape) <- c('ID','date','variable','Observed','ERA5', colnames(tasmax.all.reshape[6:15]))
  
  v <- 1;
  rmse <- vector(mode='numeric');
  r2 <- vector(mode='numeric');
  variable <- vector(mode='character')
  for(i in 5:ncol(pr.all.reshape)){
    rmse[v] <- sqrt(mean((pr.all.reshape[,4]-pr.all.reshape[,i])^2))
    r2[v] <- cor(pr.all.reshape[,4],pr.all.reshape[,i])^2
    variable[v] <- colnames(pr.all.reshape)[i]
    v <- v+1
  }
  pr.error <- data.frame(dataset = variable, rmse = rmse, r.square = r2) %>% mutate(variable='pr')
  
  v <- 1;
  rmse <- vector(mode='numeric');
  r2 <- vector(mode='numeric');
  variable <- vector(mode='character')
  for(i in 5:ncol(tasmax.all.reshape)){
    rmse[v] <- sqrt(mean((tasmax.all.reshape[,4]-tasmax.all.reshape[,i])^2))
    r2[v] <- cor(tasmax.all.reshape[,4],tasmax.all.reshape[,i])^2
    variable[v] <- colnames(tasmax.all.reshape)[i]
    v <- v+1
  }
  tasmax.error <- data.frame(dataset = variable, rmse = rmse, r.square = r2) %>% mutate(variable='tasmax')
  
  v <- 1;
  rmse <- vector(mode='numeric');
  r2 <- vector(mode='numeric');
  variable <- vector(mode='character')
  for(i in 5:ncol(tasmin.all.reshape)){
    rmse[v] <- sqrt(mean((tasmin.all.reshape[,4]-tasmin.all.reshape[,i])^2))
    r2[v] <- cor(tasmin.all.reshape[,4],tasmin.all.reshape[,i])^2
    variable[v] <- colnames(tasmin.all.reshape)[i]
    v <- v+1
  }
  tasmin.error <- data.frame(dataset = variable, rmse = rmse, r.square = r2) %>% mutate(variable='tasmin')
  
  rmse.combine <- bind_rows(pr.error,tasmax.error,tasmin.error)
  write_csv(rmse.combine, './data/climate_hist_2015/rmse.csv')
}
