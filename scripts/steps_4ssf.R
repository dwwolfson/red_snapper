#' # Tracks and steps for for fitting SSF
#' 
#' Description: This script creates tracks from snappers location data and steps ready to fit SSF. It creates steps after choosing a sampling rate, generating random available points and extracting covariates at start and end of the step. Ids 9, 12 and 16 were selected for the analysis. This script takes inputs/chunks from John's code PrelimSSF2.R and ssf_3id.R (code that tabulates location data by individual and habitat type). Finally data is saved "snap_steps" and ready to fit SSF with the following sampling rate: rate = minutes(6), tolerance = minutes(1). snap_steps is a nested track with two columns id and a list column with the random steps (including the log_sl_ and cos_ta_) for each id.
#' 
#' Programmer: JV
#' 
#' _____________________________________________________________________________
#' ## Preamble
#' 
#' Load libraries
#' 
#+ library, warnings = "hide"
library(ezknitr)
library(knitr)
library(devtools)

library(here)
library(amt)
library(lubridate)
library(dplyr)
library(readr)
library(ggplot2)
library(sf)
library(raster)
library(rgdal)
library(purrr)

#' Clear environment and set seed
#' 
remove(list=ls())
set.seed(168846)

#' _____________________________________________________________________________
#' ## Load Data and set projections
#' 
snap<-read_csv(here("Data", "red.snapper.locations.csv"))
glimpse(snap)

#ras <- raster(here("data/Seabed Maps/geotiff/", "BathyMean.tif")) # not using it for now
reefclass <- raster(here("data/Seabed Maps/geotiff/", "ChickenRock_Classification.tif"))

snap_sp <- SpatialPointsDataFrame(coords=cbind(snap$lon, snap$lat), snap,
                                  proj4string = CRS("+init=epsg:4326"))

snap <- spTransform(snap_sp, crs(reefclass)) # use reef class projection, as it's different from the one from ras (NAD83 for bathymean.tif)
snap_crs <- projection(snap_sp) # get snap projection

#' _____________________________________________________________________________
#' ## Make tracks and steps
#' 
#' ### Make tracks and filter 3 ids
#' 
snap_tracks <- snap %>%
  as.data.frame() %>% 
  mutate(num = as.factor(num)) %>% 
  filter(num == "9" | num == "12" | num == "16") %>% 
  make_track(., .x = coords.x1, .y = coords.x2, .t = datetime, crs = sp::CRS(snap_crs), id = num)

#' ### Make steps
#' 
snap_steps <- snap_tracks %>% 
  nest(data = -id) %>%  
  mutate(data = map(data, ~ .x %>%
                      track_resample(rate = minutes(6), tolerance = minutes(1)) %>%
                      filter_min_n_burst() %>% 
                      steps_by_burst() %>% 
                      random_steps() %>% # fits a tentative gamma and von Mises distribution to sl and ta, respectively, and samples random available points from these distributions
                      extract_covariates(reefclass, where="both") %>%
                      rename(reefStart = ChickenRock_Classification_start) %>%
                      rename(reefEnd = ChickenRock_Classification_end) %>%
                      mutate(reefStart = factor(reefStart, levels = 1:4, labels = c("sand", "low_relief", "medium_relief", "high_relief"))) %>% 
                      mutate(reefEnd = factor(reefEnd, levels = 1:4, labels = c("sand", "low_relief", "medium_relief", "high_relief"))) %>% #  to factors and habitat names
                      mutate(log_sl_ = log(sl_), cos_ta_ = cos(ta_))
  ))

#' ### Get rid of small steps
#' 
snap_steps <- snap_steps %>% 
  mutate(data = map(data, ~ .x %>%
                      filter(is.infinite(log_sl_)!=TRUE))) # gets rid of really small steps

glimpse(snap_steps$data)
#' 
#' 
#' _____________________________________________________________________________
#' ## Save Data
#' 
saveRDS(snap_steps, file = "data/processed_data/steps_4ssf.RDS")
#' _____________________________________________________________________________
#' ## Footer
#' 
devtools::session_info()
#' spun with:
#' ezknitr::ezspin(file = "scripts/steps_4ssf.R", keep_md = FALSE, out_dir = "output", fig_dir = "figures")
