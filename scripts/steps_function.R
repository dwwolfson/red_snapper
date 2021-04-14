# Function that creates steps using a nested_track object and a column to be mapped (e.g., tracks column with a particular sampling rate)

create_steps <- function(nested_track, column2map){
  
  #Set seed
  
  set.seed(100)
  
  #nested_track$coltouse <- nested_tracks[, column2map] # duplicated column and give name
  
  # Generate steps
  snap_steps <- nested_track %>% 
    mutate(new_step = map(column2map, ~ .x %>%
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
  
  # Get rid of small steps

  snap_steps <- snap_steps %>% 
    mutate(new_step = map(new_step, ~ .x %>%
                        filter(is.infinite(log_sl_)!=TRUE))) # gets rid of really small steps
  return(snap_steps)
}

  