#' # Update parameters after fitting reefStart SSF
#' 
#' Description: This script updates parameters of step length and turn angle distribution for each habitat class for our first model (reefStart) using purrr::map syntax. It uses as input steps generated in steps_4ssf.R and stored in process data as steps_4ssf.RDS. To update the parameters, it follows examples from here:  https://conservancy.umn.edu/handle/11299/218272  (See Appendix B)
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
#' ## Load Data
#' 
snap_steps <- readRDS(here("data", "processed_data", "steps_4ssf.RDS"))
#' 
#' _____________________________________________________________________________
#' ## Model 1: reefStart
#' 
#' Movement depends on where the individual is at the start of movement step. This is included in the model with the intraction between movement characteristics and reef habitat at the start of movement
#' 
ssf_start <- snap_steps %>%
  mutate(mod_start = map(data, function(x) (try(fit_issf(case_ ~  + sl_ + log_sl_ + cos_ta_ + (reefStart):(sl_ +  log_sl_ + cos_ta_) + strata(step_id_), data = x)))))

#glimpse(ssf_start$mod_start) # which has columns id, data, mod_start

#' _____________________________________________________________________________
#' ## Inspect model coefficients and model output structure
#' 
#' This is just inspection - no need to run

id9_mod_start <- ssf_start$mod_start[[1]] # inspect coefficients for first id. 

coefs4all <- do.call(rbind, lapply(1:length(ssf_start$mod_start), function(x) data.frame(ID = ssf_start$id[[x]], as.data.frame(ssf_start$mod_start[[x]]$model$coefficients)))) # inspect coefficients for all ids

coefs4all$coeff_name <- coefs4all %>% row.names(.) # coeff names
rownames(coefs4all) <- NULL # remove row names
#' 
#' 
#' _____________________________________________________________________________
#' ## Manually update parameters for step-length distribution and calculate speed
#' 
# Sand step-length distribution updated and stored in up_sand_sl list column. Then map the latter and calculate speed. Same for all habitats.

sand_sl <- ssf_start %>% 
  mutate(up_sand_sl = map(mod_start, ~update_gamma(
    dist = .x$sl_, # sand is the reference category
    beta_sl = .x$model$coefficients["sl_"], 
    beta_log_sl = .x$model$coefficients["log_sl_"]))) %>% 
  mutate(speed_sand = map_dbl(up_sand_sl, ~(.x$params$shape*.x$params$scale))) %>% 
  pluck("speed_sand") 

# Low_relief step-length distribution updated and stored in up_lr_sl list column

low_relief_sl <- ssf_start %>% 
  mutate(up_lr_sl = map(mod_start, ~update_gamma(
    dist = .x$sl_, 
    beta_sl = .x$model$coefficients["sl_"] +
      .x$model$coefficients["sl_:reefStartlow_relief"],
    beta_log_sl = .x$model$coefficients["log_sl_"] +
      .x$model$coefficients["log_sl_:reefStartlow_relief"]))) %>% 
  mutate(speed_low = map_dbl(up_lr_sl, ~(.x$params$shape*.x$params$scale))) %>% 
  pluck("speed_low")

# Medium_relief step-length distribution updated and stored in up_mr_sl list column

medium_relief_sl <- ssf_start %>% 
  mutate(up_mr_sl = map(mod_start, ~update_gamma(
    dist = .x$sl_, 
    beta_sl = .x$model$coefficients["sl_"] +
      .x$model$coefficients["sl_:reefStartmedium_relief"],
    beta_log_sl = .x$model$coefficients["log_sl_"] +
      .x$model$coefficients["log_sl_:reefStartmedium_relief"]))) %>% 
  mutate(speed_med = map_dbl(up_mr_sl, ~(.x$params$shape*.x$params$scale))) %>% 
  pluck("speed_med")

# High_relief step-length distribution updated and stored in up_hr_sl list column

high_relief_sl <- ssf_start %>% 
  mutate(up_hr_sl = map(mod_start, ~update_gamma(
    dist = .x$sl_, 
    beta_sl = .x$model$coefficients["sl_"] +
      .x$model$coefficients["sl_:reefStarthigh_relief"],
    beta_log_sl = .x$model$coefficients["log_sl_"] +
      .x$model$coefficients["log_sl_:reefStarthigh_relief"]))) %>% 
  mutate(speed_high = map_dbl(up_hr_sl, ~(.x$params$shape*.x$params$scale))) %>% 
  pluck("speed_high")

speeds_df <- cbind(sand_sl, low_relief_sl, medium_relief_sl, high_relief_sl) %>% as.data.frame()

speeds_df

#' _____________________________________________________________________________
#' ## Manually update parameters for turn-angle distribution

# Sand turn-angle distribution
sand_ta <- ssf_start %>% 
  mutate(up_sand_ta = map(mod_start, ~update_vonmises(
  dist = .x$ta_, 
  beta_cos_ta = .x$model$coefficients["cos_ta_"])))

# LR turn-angle distribution
low_relief_ta <- ssf_start %>% 
  mutate(up_lr_ta = map(mod_start, ~update_vonmises(
  dist = .x$ta_, 
  beta_cos_ta = .x$model$coefficients["cos_ta_"] +
    .x$model$coefficients["cos_ta_:reefStartlow_relief"])))

# MR turn-angle distribution
medium_relief_ta <- ssf_start %>% 
  mutate(up_mr_ta = map(mod_start, ~update_vonmises(
    dist = .x$ta_, 
    beta_cos_ta = .x$model$coefficients["cos_ta_"] +
      .x$model$coefficients["cos_ta_:reefStartmedium_relief"])))

# HR turn-angle distribution
high_relief_ta <- ssf_start %>% 
  mutate(up_hr_ta = map(mod_start, ~update_vonmises(
    dist = .x$ta_, 
    beta_cos_ta = .x$model$coefficients["cos_ta_"] +
      .x$model$coefficients["cos_ta_:reefStarthigh_relief"])))

#sand_ta %>% pluck("up_sand_ta") # check results
#low_relief_ta %>% pluck("up_lr_ta")
#medium_relief_ta %>% pluck("up_mr_ta")
#high_relief_ta %>% pluck("up_hr_ta")
#'
#'
#' _____________________________________________________________________________
#' ## Save Data
#' 
#' _____________________________________________________________________________
#' ## Footer
#' 
devtools::session_info()
#' spun with:
#' ezknitr::ezspin(file = "scripts/ssf_updated_params.R", keep_md = FALSE, out_dir = "output", fig_dir = "figures")
