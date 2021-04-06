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

#' Create random points
temp<-tr2 %>% nest(data=-id)%>%  
  mutate(data2 = map(
    data, ~ .x %>%
      track_resample(rate = minutes(6), tolerance = minutes(1)) %>%
      filter_min_n_burst() %>% 
      steps_by_burst() %>%
      random_steps() 
  ))
sl_distr(temp$data2[[1]])

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


#' ## JULIANA
#' 
#' Step 1: pick 3 individuals
#' 
#' - use the above temp data to pick 3 individuals
#' - tabulate for each individual (case_ by reefStart)
#' - tabulate for each individual (case_ by reefEnd)
#' - filter (case_1 i.e., used locations), then tabulate by indivdual (reefStart by reefEdn)

#' Step 2: fit 3 SSFs
#' 
#' - Model1: fit_issf(case_ ~  + sl_+log_sl_ + 
#' as.factor(reefStart):(sl_ +  log_sl_) + strata(step_id), data = data by individual)
#' 
#' - Model 2: fit_issf(case_ ~  + sl_+log_sl_ + 
#' as.factor(reefEnd) + strata(step_id), data = data by individual)
#' 
#' - Model 3: fit_issf(case_ ~  + sl_+log_sl_ + 
#' as.factor(reefStart):(sl_ +  log_sl_) + as.factor(reefEnd) strata(step_id), data = data by individual)
#' 
#' 
#' Step 3: compare above to ctmc approach to understand differences/similarities
#' Step 4: estimate transition probabilities from the 2 models/compare simulated movements, etc
#' 
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
speeds<-data.frame()
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
 speeds<-rbind(speeds, data.frame(id=rep(i,4), habs=as.factor(c("sand", "low", "medium", "high")),speed))  
}  


#+ fig.width=16, fig.height=16
ggplot(ssffits2, aes(x=id, y=estimate))+
  geom_point(size=3.5, position=position_dodge(width=0.3))+
  xlab("")+ylab(expression(hat(beta)))+facet_wrap(~term, scales="free")

ggplot(as.data.frame(speeds), aes(x=id, y=speed))+
  geom_point(size=3.5, position=position_dodge(width=0.3))+
  xlab("")+ylab(expression(hat(speed)))+facet_wrap(~habs)

#' Calculate mean/median (across individuals) of their mean speed for each habitat
speeds<-as.data.frame(speeds)
speeds %>% group_by(habs) %>% summarize(mean(speed, na.rm=TRUE), median(speed, na.rm=TRUE))
  

#' ##' Mixed SSF
#' 
#' Repeat, but using glmmTMB
#+ messages=FALSE, warning=FALSE
library(glmmTMB) 
ssfdat3$reefStart<-as.factor(ssfdat3$reefStart)
ssfdat3$id<-as.factor(ssfdat3$id)
ssfdat3$hab2<-I(ssfdat3$reefStart==2)
ssfdat3$hab3<-I(ssfdat3$reefStart==3)
ssfdat3$hab4<-I(ssfdat3$reefStart==4)
ssfdat3$hab2_sl_<-ssfdat3$hab2*ssfdat3$sl_
ssfdat3$hab3_sl_<-ssfdat3$hab3*ssfdat3$sl_
ssfdat3$hab4_sl_<-ssfdat3$hab4*ssfdat3$sl_
ssfdat3$hab2_log_sl_<-ssfdat3$hab2*ssfdat3$log_sl_
ssfdat3$hab3_log_sl_<-ssfdat3$hab3*ssfdat3$log_sl_
ssfdat3$hab4_log_sl_<-ssfdat3$hab4*ssfdat3$log_sl_

#' Warning, this takes a *long* time to fit!
ssf.tmp <- glmmTMB(case_ ~ sl_+log_sl_ + hab2_sl_+hab3_sl_+hab4_sl_+
                     hab2_log_sl_+hab3_log_sl_+hab4_log_sl_ -1 + (1|step_id)+
                     (0+sl_+hab2_sl_+hab3_sl_+hab4_sl_|id)+
                     (0+log_sl_+hab2_log_sl_+hab3_log_sl_+hab4_log_sl_|id),
                   data = ssfdat3, family=poisson(), doFit=FALSE)

                 #  control=glmmTMBControl(optimizer = "bobyqa", optCtrl = list(maxfun = 100000)))#,control=glmmTMBControl(optimizer=optim))

ssf.tmp$parameters$theta[1] <- log(1e3)
nvarparm<-length(ssf.tmp$parameters$theta)
ssf.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
#ssf.tmp$mapArg <- list(theta=factor(c(NA)))
snapper.ssf<- glmmTMB:::fitTMB(ssf.tmp)
summary(snapper.ssf)

#' Random effects
pars<-ranef(snapper.ssf)$cond$id


#' Calculate speed for each class
speeds2<-data.frame()
for(i in unique(ssffits$id)){
  tm<-ssffits %>% filter(id==i)
  re<-pars[rownames(pars)==i,]
  hab1<-update_gamma(dist=tm$stepdist[[1]],
                     beta_sl = re["sl_"],
                     beta_log_sl =re["log_sl_"])
  hab2<-update_gamma(dist=tm$stepdist[[1]],
                     beta_sl = re["sl_"] +
                       re["hab2_sl_"], 
                     beta_log_sl = re["log_sl_"] +
                       re["hab2_log_sl_"]) 
  hab3<-update_gamma(dist=tm$stepdist[[1]],
                     beta_sl = re["sl_"] +
                       re["hab3_sl_"], 
                     beta_log_sl = re["log_sl_"] +
                       re["hab3_log_sl_"]) 
  hab4<-update_gamma(dist=tm$stepdist[[1]],
                     beta_sl = re["sl_"] +
                       re["hab4_sl_"], 
                     beta_log_sl = re["log_sl_"] +
                       re["hab4_log_sl_"]) 
  speed<-c(as.numeric(hab1$params$shape)*as.numeric(hab1$params$scale),
           as.numeric(hab2$params$shape)*as.numeric(hab2$params$scale),
           as.numeric(hab3$params$shape)*as.numeric(hab3$params$scale),
           as.numeric(hab4$params$shape)*as.numeric(hab4$params$scale))
  speeds2<-rbind(speeds2, data.frame(id=rep(i,4), habs=as.factor(c("sand", "low", "medium", "high")),speed))  
}  


#' Calculate mean/median (across individuals) of their mean speed for each habitat
speeds2<-as.data.frame(speeds2)
speeds2 %>% group_by(habs) %>% summarize(mean(speed, na.rm=TRUE), median(speed, na.rm=TRUE))

ggplot(speeds2, aes(x=id, y=speed))+
  geom_point(size=3.5, position=position_dodge(width=0.3))+
  xlab("")+ylab(expression(hat(speed)))+facet_wrap(~habs)

#' compare individual fits with glmmTMB
plot(speeds$speed, speeds2$speed)
abline(0,1)


  


 
#' TRY SCALING: still problems
#' 

#' Try scaling first
ssfdat3$ssl_<-as.numeric(scale(ssfdat3$sl_))
ssfdat3$slog_sl_<-as.numeric(scale(ssfdat3$log_sl_))
ssfdat3$hab2_ssl_<-ssfdat3$hab2*ssfdat3$ssl_
ssfdat3$hab3_ssl_<-ssfdat3$hab3*ssfdat3$ssl_
ssfdat3$hab4_ssl_<-ssfdat3$hab4*ssfdat3$ssl_
ssfdat3$hab2_slog_sl_<-ssfdat3$hab2*ssfdat3$slog_sl_
ssfdat3$hab3_slog_sl_<-ssfdat3$hab3*ssfdat3$slog_sl_
ssfdat3$hab4_slog_sl_<-ssfdat3$hab4*ssfdat3$slog_sl_

ssf.tmp <- glmmTMB(case_ ~ ssl_+slog_sl_ + hab2_ssl_+hab3_ssl_+hab4_ssl_+
                     hab2_slog_sl_+hab3_slog_sl_+hab4_slog_sl_ -1+ (1|step_id)+
                     (0+ssl_+hab2_ssl_+hab3_ssl_+hab4_ssl_|id)+
                     (0+slog_sl_+hab2_slog_sl_+hab3_slog_sl_+hab4_slog_sl_|id),
                   data = ssfdat3, family=poisson(), doFit=FALSE)
ssf.tmp$parameters$theta[1] <- log(1e6)
nvarparm<-length(ssf.tmp$parameters$theta)
ssf.tmp$mapArg <- list(theta=factor(c(NA,1:(nvarparm-1))))
#ssf.tmp$mapArg <- list(theta=factor(c(NA)))
snapper.ssf<- glmmTMB:::fitTMB(ssf.tmp)
summary(snapper.ssf)


#' ## next steps
#' 
#' - can look at directional persistence from turn angles
#' - compare SSF nand ctmc 
#' 
#' ## Save data
#' 
save(ssfdat, ssfdat2, ssfdat3, file = "data/processed_data/ssfdat.rdata")
#' 