library(amt)
library(lubridate)
library(dplyr)
library(readr)
library(ggplot2)
library(sf)
library(raster)
library(rgdal)
library(sjPlot)
library(tidyverse)
library(cowplot)
library(glmmTMB)

### Load the data and raster 
snap<-read_csv("./Data/red.snapper.locations.csv")

reefclass <- raster("./data/Seabed Maps/geotiff/ChickenRock_Classification.tif")
dst_edge <- raster("./data/Seabed Maps/geotiff/dist_edge.tif")

### Prepare the data and project the data
snap <- SpatialPointsDataFrame(coords=cbind(snap$lon, snap$lat), snap,
                               proj4string = CRS("+init=epsg:4326"))
snap <- spTransform(snap, crs(dst_edge))
snapdf<-as.data.frame(snap)

#' ## Summarizing Movements

tr2<-snapdf %>% make_track(., .x=coords.x1, .y=coords.x2, .t=datetime, crs=sp::CRS("+init=epsg:4326"), id=num,
                           julian=julian, depth=depth)

tr2 %>% nest(data= -"id") %>%  mutate(sr = map(data, summarize_sampling_rate)) %>% 
  dplyr::select(id, sr) %>% 
  unnest(cols=sr) %>% print(n=nrow(tr2))

### Prepare the data
ssfdat_6 <- tr2 %>% nest(data=-id) %>% 
  mutate(data = map(
    data, ~ .x %>% 
      track_resample(rate = minutes(6), tolerance = minutes(1)) %>%
      filter_min_n_burst() %>% 
      steps_by_burst() %>% 
      random_steps() %>% 
      extract_covariates(dst_edge, where="end") %>% 
      extract_covariates(reefclass, where="both"))) %>%
  unnest(cols=data) %>%
  rename(reefStart = ChickenRock_Classification_start) %>%
  rename(reefEnd = ChickenRock_Classification_end)  %>%
  mutate(log_sl_=log(sl_), cos_ta_= cos(ta_))

ssfdat_6r <- ssfdat_6 %>% filter(is.infinite(log_sl_)!=TRUE)
ssfdat_6rf <- subset(ssfdat_6r, abs(ssfdat_6r$ta_)<3 & ssfdat_6r$sl_<100)
ssfdat_6rf$step_id <- with(ssfdat_6rf, paste(id, step_id_, sep = "_"))
glimpse(ssfdat_6rf)

ssfdat_6rf$reefS2 <- ifelse(ssfdat_6rf$reefStart == 2, 1, 0)
ssfdat_6rf$reefS3 <- ifelse(ssfdat_6rf$reefStart == 3, 1, 0)
ssfdat_6rf$reefS4 <- ifelse(ssfdat_6rf$reefStart == 4, 1, 0)
ssfdat_6rf$reefE2 <- ifelse(ssfdat_6rf$reefEnd == 2, 1, 0)
ssfdat_6rf$reefE3 <- ifelse(ssfdat_6rf$reefEnd == 3, 1, 0)
ssfdat_6rf$reefE4 <- ifelse(ssfdat_6rf$reefEnd == 4, 1, 0)

##################################################################################
########################################### Individual model fitting in amt

ind_id <-unique(ssfdat_6rf$id)

sl_df <-data.frame(ind_id)
ta_df <-data.frame(ind_id)

for(i in 1:length(ind_id))
{
  ind_temp <-subset(ssfdat_6rf, id== ind_id[i])
  sl_fit <-fit_distr(ind_temp$sl_, "gamma", na.rm=T)
  ta_fit <-fit_distr(ind_temp$ta_, "vonmises", na.rm=T)

  sl_df[i,2:3]<-sl_fit$params
  ta_df[i,2:3]<-ta_fit$params

#sl_df

##model fitting
  ssftemp <-ssfdat_6rf %>% filter(id == ind_id[i]) %>%
    fit_issf(case_ ~   sl_ + log_sl_ + cos_ta_ + as.factor(reefStart):(sl_ +  log_sl_ + cos_ta_) + 
               dist_edge + reefE2 + reefE3 + reefE4+
               strata(step_id_))
  
  #tmp[[i]] <- summary(ssftemp)$coefficients
  #print (dim(tmp[[i]]))
  ## Update the step-length distribution
  sl_df[i,4:5] <- update_gamma(dist= sl_fit, 
                               beta_sl = ssftemp$model$coefficients["sl_"],
                               beta_log_sl = ssftemp$model$coefficients["log_sl_"])$params
  
  
  #step-length distribution for low relief reefclass
  sl_df[i,6:7] <- update_gamma(dist = sl_fit,
                               beta_sl = ssftemp$model$coefficients["sl_"] +
                                 ssftemp$model$coefficients["sl_:as.factor(reefStart)2"],
                               beta_log_sl = ssftemp$model$coefficients["log_sl_"] +
                                 ssftemp$model$coefficients["log_sl_:as.factor(reefStart)2"])$params
  
  #step-length distribution for medium relief reefclasses
  sl_df[i,8:9] <- update_gamma(dist = sl_fit,
                               beta_sl = ssftemp$model$coefficients["sl_"] +
                                 ssftemp$model$coefficients["sl_:as.factor(reefStart)3"],
                               beta_log_sl = ssftemp$model$coefficients["log_sl_"] +
                                 ssftemp$model$coefficients["log_sl_:as.factor(reefStart)3"])$params
  
  #step-length distribution for high relief reefclasses
  sl_df[i,10:11] <- update_gamma(dist = sl_fit,
                                 beta_sl = ssftemp$model$coefficients["sl_"] +
                                   ssftemp$model$coefficients["sl_:as.factor(reefStart)4"],
                                 beta_log_sl = ssftemp$model$coefficients["log_sl_"] +
                                   ssftemp$model$coefficients["log_sl_:as.factor(reefStart)4"])$params
  
  ## Turn-angle distrubtuions
  ## Update the turn angle distribution
  ta_df[i,4:5] <- update_vonmises(dist= ta_fit, 
                                  beta_cos_ta = ssftemp$model$coefficients["cos_ta_"])$params
  
  
  #step-length distribution for low relief reefclass
  ta_df[i,6:7] <- update_vonmises(dist = ta_fit,
                                  beta_cos_ta = ssftemp$model$coefficients["cos_ta_"] +
                                    ssftemp$model$coefficients["sl_:as.factor(reefStart)2"])$params
  
  #step-length distribution for medium relief reefclasses
  ta_df[i,8:9] <- update_vonmises(dist = ta_fit,
                                  beta_cos_ta = ssftemp$model$coefficients["cos_ta_"] +
                                    ssftemp$model$coefficients["cos_ta_:as.factor(reefStart)3"])$params
  
  #step-length distribution for high relief reefclasses
  ta_df[i,10:11] <- update_vonmises(dist = ta_fit,
                                    beta_cos_ta = ssftemp$model$coefficients["cos_ta_"] +
                                      ssftemp$model$coefficients["cos_ta_:as.factor(reefStart)4"])$params
  
}

tmp2<-do.call(rbind, tmp)
names(sl_df)[4:11] <-c("sand_shape", "sand_scale", "low_rel_shape", "low_rel_scale",
                       "med_rel_shape", "med_rel_scale", "high_rel_shape", "high_rel_scale") 
head(sl_df)

write.csv(sl_df, "Updated_sl.csv")

names(ta_df)[4:11] <-c("sand_kappa", "sand_mu", "low_rel_kappa", "low_rel_mu",
                       "med_rel_kappa", "med_rel_mu", "high_rel_kappa", "high_rel_mu") 
head(ta_df)

write.csv(ta_df, "Updated_ta.csv")

####################################################################################
### fitting the TMB model

ssf.tmp_6m7 <- glmmTMB(case_ ~  sl_ + log_sl_ + cos(ta_) +   reefE2 + reefE3 + reefE4   + dist_edge +
                         reefS2:sl_ + reefS3:sl_ +reefS4:sl_ + reefS2:log_sl_ +reefS3:log_sl_ +reefS4:log_sl_ +
                         reefS2:cos(ta_)+reefS3:cos(ta_)+reefS4:cos(ta_)+ (1|step_id)+
                         (0+ sl_ +(log_sl_) + cos(ta_)|id) + (0+dist_edge|id)+
                         (0+reefE2 + reefE3 + reefE4 |id), REML = TRUE,
                       family=poisson(), data =ssfdat_6rf,     doFit=FALSE)

ssf.tmp_6m7$parameters$theta[1] <- log(1e3)
nvarparm7 <- length(ssf.tmp_6m7$parameters$theta)
ssf.tmp_6m7$mapArg <- list(theta=factor(c(1:(nvarparm7))))
#ssf.tmp$mapArg <- list(theta=factor(c(NA)))
snapper.ssf_67<- glmmTMB:::fitTMB(ssf.tmp_6m7)
summary(snapper.ssf_67)

##Updating the step length distribution
#### Step length analysis plot
sl_fit <-fit_distr(ssfdat_6rf$sl_, "gamma", na.rm = T)
#fit_distr(ssfdat_6rf$sl_, "gamma", na.rm=T)
#ta_fit <-fit_distr(ind_temp$ta_, "vonmises", na.rm=T)

sl_df <- as.data.frame(sl_fit$params)
#ta_df[i,2:3]<-ta_fit$params

## Update the distribution
sl_df[2,] <- as.data.frame(update_gamma(dist = sl_fit, beta_sl = coefsum["sl_",],
                                        beta_log_sl = coefsum["log_sl_",])$params)

#step-length distribution for low relief reefclass
sl_df[3,] <- update_gamma(dist = sl_fit,
                          beta_sl = coefsum["sl_",] +  coefsum["sl_:reefS2",],
                          beta_log_sl = coefsum["log_sl_",] + coefsum["log_sl_:reefS2",])$params

#step-length distribution for medium relief reefclass
sl_df[4,] <- update_gamma(dist = sl_fit,
                          beta_sl = coefsum["sl_",] +  coefsum["sl_:reefS3",],
                          beta_log_sl = coefsum["log_sl_",] + coefsum["log_sl_:reefS3",])$params

#step-length distribution for high relief reefclass
sl_df[5,] <- update_gamma(dist = sl_fit,
                          beta_sl = coefsum["sl_",] +  coefsum["sl_:reefS4",],
                          beta_log_sl = coefsum["log_sl_",] + coefsum["log_sl_:reefS4",])$params


plot_sl <-data.frame(x = rep(NA, 200))
p <-ggplot() 
col <-c(NA,"brown","tomato", "purple", "blue")
for(i in 2:nrow(sl_df)){
  #plot_sl$x <- data.frame(x = rep(NA, 100))
  plot_sl$x <- seq(from = 0, to = 100, length.out = 200)
  
  plot_sl$y <- dgamma(x = plot_sl$x,
                      shape = sl_df[i,"shape"],
                      scale = abs(sl_df[i,"scale"]))
  
  p<-p + geom_line(data= plot_sl, aes(x = x, y = y), color= col[i], size = 1) +
    theme_bw()+labs(x="Step length", y="Density")+
    theme(axis.title  = element_text(size=16),
          axis.text =   element_text(size=10), legend.position = c(0.8,0.8))
  
}  

p <- p+ scale_color_manual(name= "Habitat type", 
                           values = c("Sand"= "brown", "Low-relief"="tomato",
                                      "Medium-relief"= "purple","High-relief"= "blue"))

p

#### Turn angle analysis plot
ta_fit <-fit_distr(ssfdat_6rf$ta_, "vonmises", na.rm = T)

ta_df <- as.data.frame(ta_fit$params)

## Update the distribution
ta_df[2,] <- update_vonmises(dist = ta_fit,
                             beta_cos_ta = coefsum["cos(ta_)",])$params

#step-length distribution for low relief reefclass
ta_df[3,] <- update_vonmises(dist = ta_fit,
                             beta_cos_ta =  (coefsum["cos(ta_)",]  +  coefsum["cos(ta_):reefS2",]))$params

#step-length distribution for medium relief reefclass
ta_df[4,] <- update_vonmises(dist = ta_fit,
                             beta_cos_ta =  (coefsum["cos(ta_)",]  +  coefsum["cos(ta_):reefS3",]))$params

#step-length distribution for high relief reefclass
ta_df[5,] <- update_vonmises(dist = ta_fit,
                             beta_cos_ta =  (coefsum["cos(ta_)",]  +  coefsum["cos(ta_):reefS4",]))$params

plot_ta <- data.frame(x = rep(NA, 400))
p1<-ggplot()

for(i in 2:nrow(ta_df)){
  plot_ta$x <- seq(from = -pi, to = pi, length.out = 400)
  
  plot_ta$y <- circular::dvonmises(x = plot_ta$x,
                                   kappa = abs(ta_df[i,"kappa"]), mu=circular(pi))
  
  p1<-p1 + geom_line(data= plot_ta, aes(x = x, y = y), color= col[i], size = 1) +
    theme_bw()+labs(x="Turn-angle", y="Density")+
    theme(axis.title  = element_text(size=16),
          axis.text =   element_text(size=10), legend.position = c(0.5,0.8))
  
}  

p1 <- p1+ scale_color_manual(name= "Habitat type", 
                             values = c("Sand"= "brown", "Low-relief"="tomato",
                                        "Medium-relief"= "purple","High-relief"= "blue"))

p1

cmb_plt <-plot_grid(p, p1, labels=c("A", "B"),ncol=2)

save_plot(cmb_plt, filename = "Habitat_model_sl_ta.jpeg", 
          base_width = 12, base_height = 7,units = "in", dpi = 300)

hab_sl_ta <-cbind(sl_df, ta_df)
hab_sl_ta$habitat <-c("tentative", "sand", "low-rel", "med-rel", "high-rel")
hab_sl_ta
write.csv(hab_sl_ta, "Distribution_parameters.csv")

