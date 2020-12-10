#[ status: on-going]
#[ note: ]
#@ This code is for calculating number of exposed population and area with varying level of vulnerability to different hazard level
#@ The result is in csv format
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
  library(RColorBrewer)
}
#====================================#

## GLOBAL VARIABLES ##
{
  # The palette with black:
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  ## TH polygon
  th <- readOGR('./data/shapefile/thailand.shp')
}

## FUNCTIONS
quantify.df <- function(raster,scenario){
  new.df <- as.data.frame(raster)
  colnames(new.df) <- c('population','hazard.class','vulnerability.class')
  new.df <- new.df %>%
    drop_na() %>%
    mutate(scenario = scenario) %>%
    group_by(hazard.class, vulnerability.class) %>%
    summarize(population = round(sum(population),digits = 0),
              area = n()*0.5*0.5) %>%
    mutate(scenario = 'scenario')
  new.df
}

## QUANTIFY 
## 1. number of exposed population to hazard at different vulnerable level
## 2. area of exposure under varying hazard level where different vulnerable population reside
{
  # read each component
  {
    # hazard
    hazard.class <- stack('./output/raster/revised/hazard_class_baseline.tif',
                          './output/raster/revised/hazard_class_RCP45_SSP2_bufferzone.tif',
                          './output/raster/revised/hazard_class_RCP85_SSP5_bufferzone.tif',
                          './output/raster/revised/hazard_class_RCP45_SSP2_no_bufferzone.tif',
                          './output/raster/revised/hazard_class_RCP85_SSP5_no_bufferzone.tif'
    )
    
    names(hazard.class) <- c('hazard_baseline',
                             'hazard_RCP45_SSP2_BZ',
                             'hazard_RCP85_SSP5_BZ',
                             'hazard_RCP45_SSP2_noBZ',
                             'hazard_RCP85_SSP5_noBZ'
    )
    
    # vulnerability
    vul.class <- stack('./output/raster/revised/vulnerability_class_baseline.tif',
                       './output/raster/revised/vulnerability_class_RCP45_SSP2.tif',
                       './output/raster/revised/vulnerability_class_RCP85_SSP5.tif')
    
    names(vul.class) <- c('vulnerability_baseline',
                          'vulnerability_RCP45_SSP2',
                          'vulnerability_RCP85_SSP5')
    
    # number of population
    population <- stack('./data/raster/exposure/population_combine.grd')
    names(population) <- c('population.baseline','population.ssp2','population.ssp5')
  }
  
  # combine into relavent scenario
  {
    # set scenarios
    baseline <- stack(population$population.baseline,hazard.class$hazard_baseline,vul.class$vulnerability_baseline)
    a1 <- stack(population$population.ssp2,hazard.class$hazard_RCP45_SSP2_BZ,vul.class$vulnerability_RCP45_SSP2)
    a2 <- stack(population$population.ssp5,hazard.class$hazard_RCP85_SSP5_BZ,vul.class$vulnerability_RCP85_SSP5)
    b1 <- stack(population$population.ssp2,hazard.class$hazard_RCP45_SSP2_noBZ,vul.class$vulnerability_RCP45_SSP2)
    b2 <- stack(population$population.ssp5,hazard.class$hazard_RCP85_SSP5_noBZ,vul.class$vulnerability_RCP85_SSP5)
    
    # convert to df
    baseline.df <- quantify.df(baseline,'Baseline')
    a1.df <- quantify.df(a1,'A1')
    a2.df <- quantify.df(a2,'A2')
    b1.df <- quantify.df(b1,'B1')
    b2.df <- quantify.df(b2,'B2')
  }
  
  # clean df
  {
    result.df <- bind_rows(baseline.df,a1.df,a2.df,b1.df,b2.df)
    
    result.df <- result.df %>%
      mutate(hazard.class = ifelse(hazard.class==10,'Very Low',
                                   ifelse(hazard.class==20,'Low',
                                          ifelse(hazard.class==30,'Moderate',
                                                 ifelse(hazard.class==40,'High','Very High')))),
             vulnerability.class = ifelse(vulnerability.class==10,'Very Low',
                                          ifelse(vulnerability.class==20,'Low',
                                                 ifelse(vulnerability.class==30,'Moderate',
                                                        ifelse(vulnerability.class==40,'High','Very High')))))
    
    
    # save
    write_csv(result.df,'./output/exposed_population_under_hazard-vulnerability.csv') 
  }
}

## VISUALIZE
{
  ggplot(result.df, aes(x=factor(hazard.class), y=population, fill=factor(vulnerability.class))) +
    geom_bar(position='stack',stat='identity') +
    facet_wrap(~scenario)+
    theme_bw() +
    theme(axis.title.x = element_blank(), axis.text.x = element_text(angle=45,vjust=1,hjust=1),
          legend.position = 'bottom',legend.title = element_blank(),legend.direction = 'horizontal')
}