#'---
#'title: "Red Snapper Example iSSF"
#'author: "John Fieberg"
#'date: "`r format(Sys.time(), '%d %B, %Y')`" 
#'---
#'  
#' 
#' ## Document Preamble
#+ docPreamble, warning=FALSE, message=FALSE
library(knitr) 
library(here)
library(amt)
library(lubridate)
library(dplyr)
library(readr)
library(ggplot2)
library(sf)
library(raster)
library(rgdal)

# Chunk options
opts_chunk$set(fig.width=5, fig.height=5, warning=FALSE, message=FALSE)

#' Read in Snapper data
snap<-read_csv(here("Data", "red.snapper.locations.csv"))
glimpse(snap)

#Read in data (Vairable names)

#' - trans = transmitter number
#' - datetime = date and time of the acoustic detection
#' - julian = day of the year of acoustic detection
#' - lat = latitude of fish position (degrees N)
#' - lon = longitude of fish position (degrees W)
#' - depth = depth of fish from surface (m; note that the study area is approximately 38 m deep)
#' - num = unique fish number to be used to identify individuals in this study -since one transmitter was used on two separate fish
#' - mov = movement rate of an individual fish between that and following -location (m/s)
#' - dist = distance fish moved between that and following location (m)
#' - sec = time between that and following acoustic detection 

#' Read in the seabed maps and create a raster object for each point in the map
ras <- raster(here("data/Seabed Maps/geotiff/", "BathyMean.tif"))
reefclass <- raster(here("data/Seabed Maps/geotiff/", "ChickenRock_Classification.tif"))

#' ## Coordinate reference systems


#' Useful reference links for finding the CRS 

#' - [projfinder.com] (projfinder.com)
#' - [spatialreference.org] (spatialreference.org)

#' check crs of seabed map
crs(ras)

#' Assign a coordinate reference system to the csv file of locations that we want to extract.
#' Get location points for the original data based on x and y (x = lon, y = lat) and project to
#' epsg 4326
snap <- SpatialPointsDataFrame(coords=cbind(snap$lon, snap$lat), snap,
                               proj4string = CRS("+init=epsg:4326"))

# Transform the original data crs same as the raster layer 
snap <- spTransform(snap, crs(ras))

#' ## Summarizing Movements

#tr1<-snap %>% make_track(., .x=lat, .y=lon, .t=datetime, crs=sp::CRS("+proj=longlat +datum=WGS84"), id=trans)
#tr2<-transform_coords(tr1, sp::CRS("+init=epsg:4326"))
snapdf<-as.data.frame(snap)
tr2<-snapdf %>% make_track(., .x=coords.x1, .y=coords.x2, .t=datetime, crs=sp::CRS("+init=epsg:4326"), id=num,
                           julian=julian, depth=depth)


#' Check sampling rate by fish
#' 
#' - looks like using a resampling interval equal to every 3, 6, 12, 18, 24 minutes might work...
tr2 %>% nest(data= -"id") %>%  mutate(sr = map(data, summarize_sampling_rate)) %>% dplyr::select(id, sr) %>% 
  unnest(cols=sr) %>% print(n=nrow(tr2))

#' Resample data to 1 location every 6 minutes with a tolerance of 1 minute, then calculate 
#' step length and turn angle for each observation
tr3<-tr2 %>% nest(data= -"id") %>% mutate(steps = map(data, function(x) 
  x %>% track_resample(rate = minutes(6), tolerance = minutes(1)) %>% steps_by_burst())) 


#' ## Fit SSF models to individuals and using mixed model

temp<-tr2 %>% nest(data=-id)%>%  
  mutate(data2 = map(
    data, ~ .x %>%
      track_resample(rate = minutes(6), tolerance = minutes(1)) %>%
      filter_min_n_burst() %>% 
      steps_by_burst() %>%
      random_steps() 
  ))
sl_distr(temp$data2[[1]])

#' Create random points
ssfdat <- tr2 %>% nest(data=-id)%>%  
  mutate(data2 = map(
    data, ~ .x %>%
      track_resample(rate = minutes(6), tolerance = minutes(1)) %>%
      filter_min_n_burst() %>% 
      steps_by_burst() %>% 
      random_steps() %>% 
      extract_covariates(ras, where="end") %>% 
      extract_covariates(reefclass, where="both") %>%
      rename(reefStart = ChickenRock_Classification_start) %>%
      rename(reefEnd = ChickenRock_Classification_end) %>%  
      mutate(log_sl_=log(sl_))%>%
      filter(is.infinite(log_sl_)!=TRUE) # gets rid of really small steps)
  ))


ssfdat2 <- ssfdat[,-2] %>%
    mutate(stepdist=map(data2, ~.x %>%
         sl_distr()), tadist=map(data2,~.x %>%
         ta_distr()))  
   # unnest(cols=c(data2)) %>%
   
ssfdat3<- ssfdat2 %>% dplyr::select(c(id, data2)) %>% unnest(cols=data2) %>%
         mutate(step_id=as.factor(paste(id, step_id_, sep=":"))) # create stratification variable


 
#' Fit models to individuals. Won't work for individuals that are not found in all
#' reef classes
ssffits <- ssfdat3  %>% nest(data2=-id) %>%
  mutate(mod = map(data2, function(x) (try(fit_issf(case_ ~  + sl_+log_sl_ +
                                          as.factor(reefStart):(sl_ +  log_sl_) + strata(step_id), data = x)))))

#' Pull off coefficients for models that converged and plot
coefssf<-NULL
for(i in 1:35){
 if(attr(ssffits$mod[[i]], "class")!="try-error"){
   coefssf<-rbind(coefssf, cbind(id=ssffits$id[[i]],  broom::tidy(ssffits$mod[[i]]$model)
                  ))
 } else(print(ssffits$id[[i]]))
}



ssffits <- ssfdat3  %>% filter(id!= 26) %>% filter(id!=42) %>% nest(data2=-id) %>%
  mutate(mod = map(data2, function(x){
    fit_issf(case_ ~  + sl_+log_sl_ + as.factor(reefStart):(sl_ +  log_sl_)
              + strata(step_id), data = x, model=TRUE)}) 
    )

ssffits2 <- ssffits  %>% 
  mutate(tidy=map(mod, ~broom::tidy(.x$model))) %>% unnest(tidy)

#' Merge on tentative distribution information
ssffits<-left_join(ssffits, ssfdat2[,-2])

#' Calculate speed for each class
speeds<-NULL
for(i in unique(ssffits$id)){
 tm<-ssffits %>% filter(id==i)
 hab1<-update_gamma(dist=tm$stepdist[[1]],
                    beta_sl = tm$mod[[1]]$model$coefficients["sl_"],
                    beta_log_sl = tm$mod[[1]]$model$coefficients["log_sl_"])
 hab2<-update_gamma(dist=tm$stepdist[[1]],
                    beta_sl = tm$mod[[1]]$model$coefficients["sl_"] +
                      tm$mod[[1]]$model$coefficients["sl_:as.factor(reefStart)2"], 
                    beta_log_sl = tm$mod[[1]]$model$coefficients["log_sl_"] +
                      tm$mod[[1]]$model$coefficients["log_sl_:as.factor(reefStart)2"]) 
 hab3<-update_gamma(dist=tm$stepdist[[1]],
                    beta_sl = tm$mod[[1]]$model$coefficients["sl_"] +
                      tm$mod[[1]]$model$coefficients["sl_:as.factor(reefStart)3"], 
                    beta_log_sl = tm$mod[[1]]$model$coefficients["log_sl_"] +
                      tm$mod[[1]]$model$coefficients["log_sl_:as.factor(reefStart)3"])
 hab4<-update_gamma(dist=tm$stepdist[[1]],
                    beta_sl = tm$mod[[1]]$model$coefficients["sl_"] +
                      tm$mod[[1]]$model$coefficients["sl_:as.factor(reefStart)4"], 
                    beta_log_sl = tm$mod[[1]]$model$coefficients["log_sl_"] +
                      tm$mod[[1]]$model$coefficients["log_sl_:as.factor(reefStart)4"])
 speed<-c( hab1$params$shape*hab1$params$scale,
            hab2$params$shape*hab2$params$scale,
            hab3$params$shape*hab3$params$scale,
            hab4$params$shape*hab4$params$scale)
 speeds<-rbind(speeds, cbind(id=rep(i,4), habs=1:4,speed))  
}  


#+ fig.width=16, fig.height=16
ggplot(ssffits2, aes(x=id, y=estimate))+
  geom_point(size=3.5, position=position_dodge(width=0.3))+
  xlab("")+ylab(expression(hat(beta)))+facet_wrap(~term, scales="free")

ggplot(as.data.frame(speeds), aes(x=id, y=speed))+
  geom_point(size=3.5, position=position_dodge(width=0.3))+
  xlab("")+ylab(expression(hat(speed)))+facet_wrap(~habs, scales="free")

speeds<-as.data.frame(speeds)
speeds %>% group_by(habs) %>% summarize(mean(speed, na.rm=TRUE), median(speed, na.rm=TRUE))
  
#' ## next steps
#' 
#' - do the same, but using glmmTMB so can fit model to all individuals
#' will need to update code to update gamma distribution & calculate speed
#' 
#' - can look at directional persistence for turn angles
#' - compare SSF nand ctmc 
#' 
  
#' Fit using glmmTMB
#+ messages=FALSE, warning=FALSE
library(glmmTMB) 
ssf.tmp <- glmmTMB(case_ ~ sl_+log_sl_ + as.factor(reefStart):(sl_+log_sl_)+
                      (1|step_id) +  + I(0+as.factor(reefStart):(sl_+log_sl_)|id),
                   data = ssfdat2, family=poisson(),
                    doFit=FALSE)

#ssf.tmp <- glmmTMB(case_ ~ as.factor(reefEnd) + sl_+log_sl_ +
#                     (1|step_id), family=poisson(), data = ssfdat,
#                   doFit=FALSE)
ssf.tmp$parameters$theta[1] <- log(1e3)
nvarparm<-length(ssf.tmp$parameters$theta)
ssf.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
#ssf.tmp$mapArg <- list(theta=factor(c(NA)))
snapper.ssf<- glmmTMB:::fitTMB(ssf.tmp)
summary(snapper.ssf)


 
