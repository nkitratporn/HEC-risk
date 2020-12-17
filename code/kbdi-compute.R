#=====================================#
## KBDI Calculation
#@ KBDI equation and symbol explained
########
#@ KBDI = KBDIy - (100  * NPR) + DF                                                              (eq.1)
#@ dQ   = (800.0-KBDIy)*(0.968*exp(0.0486*tasmax)-8.30)/1000.0/(1.0+10.88*exp(-0.0441*prAnnaul)) (eq.2)
#@ dQ   = 0  if tasmax < 50F else dQ = dQ                                                        (eq.3) 
########
#@ KBDI    : KBDI value at day t
#@ KBDIy   : KBDI value at day t-1
#@ NPR     : net rainfall at day t
#@ pr      : total 24-h rainfall at t in inch
#@ tasmax  : maximum temperature at t in Fahrenheit
#@ prAnnual: mean annual precipitation
#@ NF      : net precipitation
#@ dQ      : drought factor calculate using eq.2
########
#@ @citation: 
#@  Keetch, J.J. and Byram, G.M. (1968) A drought index for forest fire control. USDA Forest Service.
#@  Alexander, M.E., 1990. Computer calculation of the Keetch-Byram Drought Index - programmers beware. Fire Management Notes 51, 23ย–25.
#=====================================#

#========loading library============= 
## data handling
library(tidyverse)
library(lubridate)

## raster and shapefile processing ##
library(rgeos)
library(rgdal)
library(raster)
library(sp)
library(sf)

## visualization
library(rasterVis)
library(ggplot2)

## generating graph 
library(maptools)

## c++ communication (lower load with recursive function/looping)
library(Rcpp)

## parallel processing
library(foreach)
library(doParallel)

#======KBDI CPP Calc========
## create function for kbdi calculation
## solution 1 using C++ to speedup recursive function
## if(kbdi > 800.0) { kbdi = 800.0; }
cppFunction('NumericVector kbdiC(NumericVector x, NumericVector y, NumericVector z, NumericVector s){
	
	int n = x.size();
    double kbdi;
	NumericVector out(n);
    
    for (int i = 0; i < n; i++){
        if(s[i]==1) { 
            kbdi = x[i] + ((800-x[i])*y[i]) - (100*z[i]);
        }
        else if (s[i]==0){
            kbdi = kbdi + ((800-kbdi)*y[i]) - (100*z[i]); 
        }
        else {
            kbdi = NA_REAL;
        }

        if(kbdi < 0.0) { kbdi = 0.0; }
    
        out[i] = kbdi;
    }
    return out;
}')

#======HISTORICAL ========

## check number of files: pr-precipitation & tasmax-max.temperature
pr.files <- list.files('data/ERA5-mask/',pattern = '.*pr-([0-9]+).*')
tasmax.files <- list.files('data/ERA5-mask/',pattern = '.*tasmax-([0-9]+).*')

pr.files
tasmax.files

# annual average precipitation (pr)
{
  prAnnual <- data.frame()
  ## calculate sum of pr for each year
  for (i in 2:length(pr.files)) {
    start <- proc.time()[[3]]
    cat('load raster stack for ', sub(".*(\\d+{4}).*$", "\\1", pr.files[i]),'\n')
    
    pr <- stack(paste0('data/ERA5-mask/',pr.files[i]))
    first.date <- as.Date(paste0(gsub(".*(\\d+{4}).*$", "\\1",names(pr[[1]])),'0101'), format = '%Y%m%d')
    
    #========CONVERT TO DATAFRAME===========#
    
    ### 1. convert to dataframe
    pr.df <- data.frame()
    
    #Register CoreCluster
    UseCores <- detectCores() - round(detectCores()*0.15, digits = 0)  #Define how many cores you want to use [leave around 15%]
    cl <- makeCluster(UseCores)
    registerDoParallel(cl)
    
    pr.df <- foreach(j=1:nlayers(pr), .packages = 'raster',
                     .combine = 'cbind',.export = 'pr.df') %dopar% {
                       if (j > 1){
                         temp <- as.data.frame(pr[[j]],xy=T)
                         temp <- as.data.frame(temp[,-c(1,2)])
                         colnames(temp)[1] <- names(pr[[j]])
                       }else{
                         temp <- as.data.frame(pr[[j]],xy=T)
                         colnames(temp)[3] <- names(pr[[j]])
                       }
                       temp
                     }
    
    
    pr.df$location <- seq(1,nrow(pr.df))
    rm(pr)
    
    pr.df.m <- pr.df %>% 
      reshape2::melt(id.vars = c('location','x','y')) %>%
      mutate(date = first.date -1 + as.numeric(gsub("(?:[^.]+\\.){2}([^.]+).*", "\\1",variable))) %>%
      rename(pr = value) %>%
      dplyr::filter(!is.na(pr)) %>%
      mutate(pr = ifelse(pr<0,0,pr))
    
    stopCluster(cl)
    rm(pr.df)
    cat('df: COMPLETED... ')
    
    ## extract date info into day month and year
    pr.df.m <- pr.df.m %>%
      mutate(year = year(date))
    
    ## convert to nested data
    pr.df.m <- pr.df.m %>%
      nest(data= -c(location, x, y,year))
    
    ## calculate sum pr by year
    pr.df.m  <- pr.df.m  %>%
      mutate(pr_sum = map_dbl(data, ~sum(.$pr)))
    cat('annual calc: COMPLETED... ')
    
    pr.df.m <- pr.df.m %>%
      dplyr::select(-c(data))
    
    prAnnual <- prAnnual %>%
      bind_rows(pr.df.m)
    cat('prAnnual_sum: COMPLETED... ',proc.time()[[3]]-start,'\n')
    rm(pr.df.m)
  }
  glimpse(prAnnual)
  
  ## convert to nested data by location
  prAnnual.df.m <- prAnnual %>%
    nest(data=-c(location, x, y))
  
  ## calculate mean of annual pr
  prAnnual.df.m <- prAnnual.df.m %>%
    mutate(prAnnual = map_dbl(data,~mean(.$pr_sum))) %>%
    dplyr::select(-c(data))
  
  glimpse(prAnnual.df.m)
  head(prAnnual.df.m)
  
  ## save
  prAnnual.df.m %>% write_csv('data/ERA5-mask/prAnnual1990-2019.csv')
}

# daily KBDI
{
  last.kbdi <- data.frame()
  #i <- 2
  for (i in 1:length(pr.files)) { #i in 1:length(pr.files)
    start <- proc.time()[[3]]
    cat('load raster stack ', sub(".*(\\d+{4}).*$", "\\1", pr.files[i]),'\n')
    
    pr <- stack(paste0('data/ERA5-mask/',pr.files[i]))
    tasmax <- stack(paste0('data/ERA5-mask/', tasmax.files[i]))
    
    if(i>1){
      first.date <- as.Date(paste0(gsub(".*(\\d+{4}).*$", "\\1",names(pr[[1]])),'0101'), format = '%Y%m%d')
    }
    
    #========CONVERT TO DATAFRAME===========#
    
    ### 1. convert to dataframe
    pr.df <- data.frame()
    tasmax.df <- data.frame()
    ## extract regional-scale mean and stddev into df
    
    #Register CoreCluster
    UseCores <- detectCores() - round(detectCores()*0.15, digits = 0) #Define how many cores you want to use
    cl <- makeCluster(UseCores)
    registerDoParallel(cl)
    
    pr.df <- foreach(j=1:nlayers(pr), .packages = 'raster',
                     .combine = 'cbind',.export = 'pr.df') %dopar% {
                       if (j > 1){
                         temp <- as.data.frame(pr[[j]],xy=T)
                         temp <- as.data.frame(temp[,-c(1,2)])
                         colnames(temp)[1] <- names(pr[[j]])
                       }else{
                         temp <- as.data.frame(pr[[j]],xy=T)
                         colnames(temp)[3] <- names(pr[[j]])
                       }
                       temp
                     }
    
    tasmax.df <- foreach(j=1:nlayers(tasmax), .packages = 'raster',
                         .combine = 'cbind',.export = 'tasmax.df') %dopar% {
                           if (j > 1){
                             temp <- as.data.frame(tasmax[[j]],xy=T)
                             temp <- as.data.frame(temp[,-c(1,2)])
                             colnames(temp)[1] <- names(tasmax[[j]])
                           }else{
                             temp <- as.data.frame(tasmax[[j]],xy=T)
                             colnames(temp)[3] <- names(tasmax[[j]])
                           }
                           temp
                         }
    
    
    pr.df$location <- seq(1,nrow(pr.df))
    tasmax.df$location <- seq(1,nrow(tasmax.df))
    
    rm(pr,tasmax)
    
    ### 2. melt by column names and extract date from variable name, change column name
    if(i==1){
      pr.df.m <- pr.df %>% 
        reshape2::melt(id.vars = c('location','x','y')) %>%
        mutate(date = as.Date(gsub('.*X([0-9]+).*','\\1',variable), format = '%Y%m%d'),
               value = ifelse(value <0, 0, value)) %>%
        rename(pr = value)
      
      tasmax.df.m <- tasmax.df %>%
        reshape2::melt(id.vars = c('location','x','y')) %>%
        mutate(date = as.Date(gsub('.*X([0-9]+).*','\\1',variable), format = '%Y%m%d')) %>%
        rename(tasmax = value)
    }else{
      pr.df.m <- pr.df %>% 
        reshape2::melt(id.vars = c('location','x','y')) %>%
        mutate(date = first.date -1 + as.numeric(gsub("(?:[^.]+\\.){2}([^.]+).*", "\\1",variable)),
               value = ifelse(value <0, 0, value)) %>%
        rename(pr = value) %>%
        dplyr::filter(!is.na(pr))
      
      tasmax.df.m <- tasmax.df %>%
        reshape2::melt(id.vars = c('location','x','y')) %>%
        mutate(date = first.date -1 + as.numeric(gsub("(?:[^.]+\\.){2}([^.]+).*", "\\1",variable))) %>%
        rename(tasmax = value) %>%
        dplyr::filter(!is.na(tasmax))
    }
    
    stopCluster(cl)
    
    rm(pr.df,tasmax.df)
    cat('df: COMPLETED ')
    
    #========CALCULATION STEPS===========#
    # 1. calculate net precipitation (NF) dataframe
    #@ NPR is computed by subtracting 0.2 in. from the value of daily rainfall.
    #@ If there are consecutive wet days, 0.2 in. is subtracted only once on the day
    #@ when cumulative rainfall exceeds 0.2. 
    #@ A wet period ends when two rainy days are separated by one day without measurable rainfall,
    #@ thus 0.2 has to be subtracted again in the next rain period.
    
    ## calculate cumulative rainfall
    pr.df.m <- pr.df.m %>% 
      group_by(location) %>%
      arrange(date, .by_group = T) %>% 
      mutate(temp = cumsum(pr==0)) %>%
      ungroup()
    
    pr.df.m <- pr.df.m %>%
      group_by(location,temp) %>%
      mutate(cRain = cumsum(pr)) %>%
      ungroup()
    cat('cRain: COMPLETED| ')
    
    ## set counter to reduce 0.2 in. start at 1 from cRain > 0.2 and reset when 0 again
    ## create new group whenever cRain < 0.2  and set 1 to anything cRain > 0.2 and cumsum
    pr.df.m <- pr.df.m %>%
      group_by(location) %>%
      mutate(grp = cumsum(cRain < 0.2), counter = ifelse(cRain > 0.2,1,0)) %>%   
      ungroup()
    
    ## calculate cumsum under location and grp
    pr.df.m <- pr.df.m %>%
      group_by(location, grp) %>%
      mutate(counterTemp = cumsum(counter!=0)) %>% 
      ungroup() %>%
      ## deduct 0.2 from counterTemp==1 (first consecutive day with cRain > 0.2), then directly use pr afterward. otherwise nf is 0
      mutate(nf = ifelse(counterTemp==1, cRain-0.2, ifelse(counterTemp>1, pr,0))) %>%
      dplyr::select(-c(temp,grp,counter,counterTemp))
    
    cat('NF: COMPLETED ')
    
    # 2. calculate kbdi
    ## join tasmax and prAnnual to pr
    pr.df.m <- pr.df.m %>%
      full_join(dplyr::select(tasmax.df.m, c(location,date,tasmax)), by = c('location','date'))
    
    pr.df.m <- pr.df.m %>%
      left_join(dplyr::select(prAnnual.df.m, c(location,prAnnual)), by = c('location'))
    
    ## calculate PET part inside dQ formula
    #@ dQ   = 10^(-3)*(800-Q0) [0.968exp(0.0486T) - 8.3] dt / [1+10.88 exp (-0.0441 R)]
    pr.df.m <- pr.df.m %>%
      mutate(pet = ifelse(tasmax>=50,(0.968*exp(0.0486*tasmax)-8.30)/1000.0/(1.0+10.88*exp(-0.0441*prAnnual)),0))
    
    if (nrow(last.kbdi)==0) {
      pr.df.m <- pr.df.m %>%
        group_by(location) %>%
        mutate(wetspell = ifelse(cRain > 5, 1, ifelse(cRain==max(cRain), 1, 0)),
               start = cumsum(wetspell!=0)) %>%
        mutate(init = ifelse(start==1 & wetspell == 1, 1, ifelse(start==0 & wetspell==0, NA, 0))) %>%
        mutate(last_kbdi = ifelse(init==1,0,NA)) %>%
        ungroup() %>%
        dplyr::select(-c(wetspell,start))
      
    }else{
      pr.df.m <- pr.df.m %>%
        left_join(last.kbdi, by='location') %>%
        group_by(location) %>%
        mutate(init = ifelse(row_number()==1, 1,0))
    }
    cat('last:',nrow(last.kbdi),'init: COMPLETED..')
    
    ## call kbdi function
    kbdi.df.m <- pr.df.m %>%
      group_by(location) %>%
      mutate(kbdi = kbdiC(last_kbdi,pet,nf,init)) %>%
      ungroup() %>%
      dplyr::select(location,x,y,date,kbdi)
    
    cat('KBDI: COMPLETED ')
    
    last.kbdi <- kbdi.df.m %>%
      group_by(location) %>%
      slice(n()) %>%
      dplyr::select(location,kbdi) %>%
      rename(last_kbdi = kbdi) %>%
      ungroup()
    
    rm(pr.df.m)
    
    write_csv(kbdi.df.m,path = paste0('output/kbdi-df/kbdi-',sub(".*(\\d+{4}).*$", "\\1", pr.files[i]),'.csv'),
              col_names = T)
    cat('save: COMPLETED ')
    
    end <- proc.time()[[3]]
    cat(end-start,'sec \n')
  }
}

#======FUTURE ========
## set scenario and models 
scenario <- c('rcp45','rcp85')
models <- c('CESM1-BGC','MPI-ESM-MR','MIROC5', 'IPSL-CM5A-MR','CanESM2')

## list all files
pr.files <- list.files('data/NEX-mask/',pattern = '.*pr-*')
tasmax.files <- list.files('data/NEX-mask/',pattern = '.*tasmax*')

pr.files; tasmax.files;

## load last date KBDI from previous year (2019-12-31)
last.kbdi2015 <- read_csv('output/kbdi-df/kbdi-2015.csv', col_names = T, cols(
                          location = col_double(),
                          x=col_double(),
                          y=col_double(),
                          date=col_date(format = ''),
                          kbdi=col_double()
                          )) %>%
  group_by(location) %>%
  slice(n()) %>%
  dplyr::select(location,kbdi) %>%
  rename(last_kbdi = kbdi) %>%
  ungroup()

head(last.kbdi2015)

## daily KBDI
{
  for (k in seq(scenario)) {
    cat(scenario[k],': ')
    for (m in seq(models)) {
      cat(models[m],'\n')
      pr.select <- pr.files[grep(pr.files, pattern = paste0('.*',scenario[k],'.*m',m,'.*'))]
      tasmax.select <- tasmax.files[grep(tasmax.files, pattern = paste0('.*',scenario[k],'.*m',m,'.*'))]
      
      # test if number of files of pr and tasmax is equal
      if (length(pr.select)!=length(tasmax.select)) {
        cat('\n number of pr and tasmax files is not equal!')
        break
      }
      if (length(pr.select)==0 | length(tasmax.select)==0) {
        cat('\n no files!')
        break
      }
      
      for (n in 1:length(pr.select)) {# length(pr.select)
        start <- proc.time()[[3]]
        cat('load raster stacks: ', sub(".*(\\d+{4}).*$", "\\1", pr.select[n]),'\n')
        
        # load raster stacks
        pr <- stack(paste0('data/NEX-mask/', pr.select[n]))
        tasmax <- stack(paste0('data/NEX-mask/', tasmax.select[n]))
        
        # setup first date in stack
        first.date <- as.Date(paste0(gsub(".*(\\d+{4}).*$", "\\1",names(pr[[1]])),'0101'), format = '%Y%m%d')
        
        #========CONVERT TO DATAFRAME===========#
        
        ### 1. convert to dataframe
        pr.df <- data.frame()
        tasmax.df <- data.frame()
        
        #Register CoreCluster
        UseCores <- detectCores() - 10 #Define how many cores you want to use
        cl <- makeCluster(UseCores)
        registerDoParallel(cl)
        
        # convert raster to df in parallel 
        pr.df <- foreach(j=1:nlayers(pr), .packages = 'raster',
                         .combine = 'cbind',.export = 'pr.df') %dopar% {
                           if (j > 1){
                             temp <- as.data.frame(pr[[j]],xy=T)
                             temp <- as.data.frame(temp[,-c(1,2)])
                             colnames(temp)[1] <- names(pr[[j]])
                           }else{
                             temp <- as.data.frame(pr[[j]],xy=T)
                             colnames(temp)[3] <- names(pr[[j]])
                           }
                           temp
                         }
        
        tasmax.df <- foreach(j=1:nlayers(tasmax), .packages = 'raster',
                             .combine = 'cbind',.export = 'tasmax.df') %dopar% {
                               if (j > 1){
                                 temp <- as.data.frame(tasmax[[j]],xy=T)
                                 temp <- as.data.frame(temp[,-c(1,2)])
                                 colnames(temp)[1] <- names(tasmax[[j]])
                               }else{
                                 temp <- as.data.frame(tasmax[[i]],xy=T)
                                 colnames(temp)[3] <- names(tasmax[[j]])
                               }
                               temp
                             }
        
        # add location id
        pr.df$location <- seq(1,nrow(pr.df))
        tasmax.df$location <- seq(1,nrow(tasmax.df))
        
        # remove raster stack (reduce memory occupency)
        rm(pr,tasmax)
        cat('df: COMPLETED ')
        ### 2. melt by column names and extract date from variable name, change column name
        pr.df.m <- pr.df %>% 
          reshape2::melt(id.vars = c('location','x','y')) %>%
          mutate(date = first.date -1 + as.numeric(gsub("(?:[^.]+\\.){4}([^.]+).*", "\\1",variable))) %>%
          rename(pr = value) %>%
          dplyr::filter(!is.na(pr)) %>%
          mutate(pr = ifelse(pr<0,0,pr))
        
        
        tasmax.df.m <- tasmax.df %>%
          reshape2::melt(id.vars = c('location','x','y')) %>%
          mutate(date = first.date -1 + as.numeric(gsub("(?:[^.]+\\.){4}([^.]+).*", "\\1",variable))) %>%
          rename(tasmax = value) %>%
          dplyr::filter(!is.na(tasmax))
        
        stopCluster(cl)
        
        rm(pr.df,tasmax.df)
        cat('melt: COMPLETED ')
        
        #========CALCULATION STEPS===========#
        # 1. calculate net precipitation (NF) dataframe
        #@ NPR is computed by subtracting 0.2 in. from the value of daily rainfall.
        #@ If there are consecutive wet days, 0.2 in. is subtracted only once on the day
        #@ when cumulative rainfall exceeds 0.2. 
        #@ A wet period ends when two rainy days are separated by one day without measurable rainfall,
        #@ thus 0.2 has to be subtracted again in the next rain period.
        
        ## calculate cumulative rainfall
        pr.df.m <- pr.df.m %>% 
          group_by(location) %>%
          arrange(date, .by_group = T) %>% 
          mutate(temp = cumsum(pr==0)) %>%
          ungroup()
        
        pr.df.m <- pr.df.m %>%
          group_by(location,temp) %>%
          mutate(cRain = cumsum(pr)) %>%
          ungroup()
        cat('cRain: COMPLETED ')
        
        ## set counter to reduce 0.2 in. start at 1 from cRain > 0.2 and reset when 0 again
        ## create new group whenever cRain < 0.2  and set 1 to anything cRain > 0.2 and cumsum
        pr.df.m <- pr.df.m %>%
          group_by(location) %>%
          mutate(grp = cumsum(cRain < 0.2), counter = ifelse(cRain > 0.2,1,0)) %>%   
          ungroup()
        
        ## calculate cumsum under location and grp
        pr.df.m <- pr.df.m %>%
          group_by(location, grp) %>%
          mutate(counterTemp = cumsum(counter!=0)) %>% 
          ungroup() %>%
          ## deduct 0.2 from counterTemp==1 (first consecutive day with cRain > 0.2), then directly use pr afterward. otherwise nf is 0
          mutate(nf = ifelse(counterTemp==1, cRain-0.2, ifelse(counterTemp>1, pr,0))) %>%
          dplyr::select(-c(temp,grp,counter,counterTemp))
        
        cat('NF: COMPLETED ')
        
        # 2. calculate kbdi
        ## join tasmax and prAnnual to pr
        pr.df.m <- pr.df.m %>%
          full_join(dplyr::select(tasmax.df.m, c(location,date,tasmax)), by = c('location','date'))
        
        pr.df.m <- pr.df.m %>%
          left_join(dplyr::select(prAnnual.df.m, c(location,prAnnual)), by = c('location'))
        
        ## calculate PET part inside dQ formula
        #@ dQ   = 10^(-3)*(800-Q0) [0.968exp(0.0486T) - 8.3] dt / [1+10.88 exp (-0.0441 R)]
        pr.df.m <- pr.df.m %>%
          mutate(pet = ifelse(tasmax>=50,(0.968*exp(0.0486*tasmax)-8.30)/1000.0/(1.0+10.88*exp(-0.0441*prAnnual)),0))
        
        if(n==1){
          last.kbdi <- last.kbdi2015
          cat('kbdi2015... ')
        }
        
        pr.df.m <- pr.df.m %>%
          left_join(last.kbdi, by='location') %>%
          group_by(location) %>%
          mutate(init = ifelse(row_number()==1, 1,0))
        
        cat('init: COMPLETED..')
        
        ## call kbdi function
        kbdi.df.m <- pr.df.m %>%
          group_by(location) %>%
          mutate(kbdi = kbdiC(last_kbdi,pet,nf,init)) %>%
          ungroup() %>%
          dplyr::select(location,x,y,date,kbdi)
        
        cat('KBDI: COMPLETED ')
        
        last.kbdi <- kbdi.df.m %>%
          group_by(location) %>%
          slice(n()) %>%
          dplyr::select(location,kbdi) %>%
          rename(last_kbdi = kbdi) %>%
          ungroup()
        
        write_csv(kbdi.df.m,path = paste0('output/kbdi-df/kbdi-',scenario[k],'-m',m,'-',sub(".*(\\d+{4}).*$", "\\1", pr.select[n]),'.csv'),
                  col_names = T)
        cat('save: COMPLETED ')
        
        end <- proc.time()[[3]]
        cat(end-start,'sec \n')
      }
    } 
  }
}