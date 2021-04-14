#' # Tracks and steps for for fitting SSF
#' 
#' Description: This script creates tracks from snappers location data and steps ready to fit SSF. It creates steps after choosing a sampling rate, generating random available points and extracting covariates at start and end of the step. Ids 9, 12 and 16 were selected for the analysis. This script takes inputs/chunks from John's code PrelimSSF2.R and ssf_3id.R (code that tabulates location data by individual and habitat type). Finally data is saved "snap_steps" and ready to fit SSF with different sampling rates. snap_steps is a nested track with columns id and a single list column including steps from different sampling rates and the random steps (including the log_sl_ and cos_ta_) for each id. To facilitate generating steps for different sampling rates it uses the create_steps function (steps_function.R) stored in the scripts folder.
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
library(tidyr)

#' Clear environment and set seed
#' 
remove(list=ls())
set.seed(168846)

#' Source step function
#' 
source(here("scripts", "steps_function.R"))

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
  make_track(., .x = coords.x1, .y = coords.x2, .t = datetime, crs = sp::CRS(snap_crs), id = num) %>% 
  nest(data = -id)

#' ### Summarize sampling rate for each individual
#' 
samp_rate <- snap_tracks %>% 
  mutate(samp_rate = map_df(.x = data, .f = summarize_sampling_rate)) %>% 
  pluck("samp_rate")

samp_rate

#' ### Make tracks with different sampling rates
#' 
snap_tracks_sr <- snap_tracks %>%  
  mutate(track_2mins = map(data, ~ .x %>%
                           track_resample(rate = minutes(2), tolerance = minutes(1)))) %>%
  mutate(track_5mins = map(data, ~ .x %>%
                              track_resample(rate = minutes(5), tolerance = minutes(1)))) %>% 
  mutate(track_10mins = map(data, ~ .x %>%
                              track_resample(rate = minutes(10), tolerance = minutes(1)))) %>%
  mutate(track_30mins = map(data, ~ .x %>%
                               track_resample(rate = minutes(30), tolerance = minutes(1)))) %>%
  mutate(track_60mins = map(data, ~ .x %>%
                               track_resample(rate = minutes(60), tolerance = minutes(1))))
  
#' ### Use tracks with different sampling rates to create steps
#' 
snap_steps <- create_steps(nested_track = snap_tracks_sr, column2map = snap_tracks_sr$track_2mins) %>% rename(steps_2mins = new_step) %>%
  create_steps(column2map = snap_tracks_sr$track_5mins) %>% rename(steps_5mins = new_step) %>% 
  create_steps(column2map = snap_tracks_sr$track_10mins) %>% rename(steps_10mins = new_step) %>% 
  create_steps(column2map = snap_tracks_sr$track_30mins) %>% rename(steps_30mins = new_step) %>% 
  create_steps(column2map = snap_tracks_sr$track_60mins) %>% rename(steps_60mins = new_step)
# 
snap_steps <- snap_steps %>% pivot_longer(cols = c(steps_2mins, steps_5mins, steps_10mins, steps_30mins, steps_60mins), names_to = "samp_rate", values_to = "steps") %>% mutate(id = paste(id, samp_rate, sep = "-")) %>% select(id, samp_rate, steps) # go long

glimpse(snap_steps)

#vec <- c(snap_tracks_sr$samp_rate_2m, snap_tracks_sr$samp_rate_5m)
#l <- lapply(vec, function(x) create_steps(nested_track = snap_tracks_sr, vec))
#
#test <- snap_tracks_sr %>% 
  #mutate(steps_all = pmap(l, create_steps(nested_track = snap_tracks_sr, column2map = snap_tracks_sr$samp_rate_2m))) # failed attempt to do previous code in a single line
# 
#test_1 <- snap_tracks_sr %>% map2(snap_tracks_sr$samp_rate_2m, snap_tracks_sr$samp_rate_5m, ~create_steps)
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
