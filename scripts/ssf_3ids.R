#' ## Fit SSF only for three individuals
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
library(tidyr)

remove(list=ls())
set.seed(168847)

#' 
#' Step 1: pick 3 individuals
#' 
#' - use the above temp data to pick 3 individuals
#' - tabulate for each individual (case_ by reefStart)
#' - tabulate for each individual (case_ by reefEnd)
#' - filter (case_1 i.e., used locations), then tabulate by indivdual (reefStart by reefEdn)
#' 
#' ## Load data
#' 
load(here("data/processed_data", "ssfdat.rdata")) # ssfdat from John's code (see PrelimSSF2.R).

#' ## Unnest and reformat data for tabulation
#' 
dat4tab <- ssfdat %>% 
  tidyr::unnest(data2) %>% 
  mutate(reefStart = factor(reefStart, levels = 1:4, labels = c("sand", "low_relief", "medium_relief", "high_relief"))) %>% 
  mutate(reefEnd = factor(reefEnd, levels = 1:4, labels = c("sand", "low_relief", "medium_relief", "high_relief"))) %>% #  to factors and habitat names
  as.data.frame()  %>%  
  dplyr::select(-data) %>% 
  filter(case_ == TRUE)

#' Table for reefStart
#' 
tabStart <- dat4tab %>% 
  group_by(id, case_, reefStart) %>% 
  summarise(N = n()) %>% group_by(id) %>% 
  #filter(!any(is.na(reefStart))) %>% #filter out id groups with NAs
  spread(reefStart, N) # go wide

#' Table for reefEnd
#' 
tabEnd <- dat4tab %>% 
  group_by(id, case_, reefEnd) %>% 
  summarise(N = n()) %>% group_by(id) %>% 
  #filter(!any(is.na(reefEnd))) %>% #filter out id groups with NAs
  spread(reefEnd, N)  # go wide

# ids with a high number of used locations for different habitat types can be 9, 12, 16. Id 16 moves all over the place, 12 uses a very small area and both like a lot the low and high relief habitats vs id 9 that likes the sand (see plots in Drive/figures). I wanted to capture contrast in habitat use.

#' Step 2: fit 3 SSFs
#' 
#' Use ssfdat3 as it already has the stratification variable and it's ready for fitting the models
#' 
#' Filter by the selected ids 9, 12,16
#' 
ssfdat_3ids <- ssfdat3 %>% 
  filter(id == "9" | id == "12" | id == "16")
#'
#' Model1: fit_issf(case_ ~  + sl_+log_sl_ + 
#' as.factor(reefStart):(sl_ +  log_sl_) + strata(step_id), data = data by individual)
#'
ssfm1 <- ssfdat_3ids  %>% nest(data2=-id) %>%
  mutate(mod1 = map(data2, function(x) (try(fit_issf(case_ ~  + sl_+log_sl_ +
                                                      as.factor(reefStart):(sl_ +  log_sl_) + strata(step_id), data = x)))))

#' Model 2: fit_issf(case_ ~  + sl_+log_sl_ + 
#' as.factor(reefEnd) + strata(step_id), data = data by individual)

ssfm2 <- ssfdat_3ids  %>% nest(data2=-id) %>%
  mutate(mod2 = map(data2, function(x) (try(fit_issf(case_ ~  + sl_+log_sl_ +
                                                      as.factor(reefEnd):(sl_ +  log_sl_) + strata(step_id), data = x)))))

#' Model 3: fit_issf(case_ ~  + sl_+log_sl_ + 
#' as.factor(reefStart):(sl_ +  log_sl_) + as.factor(reefEnd) strata(step_id), data = data by individual)
#' 
ssfm3 <- ssfdat_3ids  %>% nest(data2=-id) %>%
  mutate(mod3 = map(data2, function(x) (try(fit_issf(case_ ~  + sl_+log_sl_ +
                                                      as.factor(reefStart):(sl_ +  log_sl_) + as.factor(reefEnd) + strata(step_id), data = x)))))
#' 
#' Pull off coefficients for first model
coefssf1<-NULL
for(i in 1:3){
  if(attr(ssfm1$mod1[[i]], "class")!="try-error"){
    coefssf1<-rbind(coefssf1, cbind(id=ssfm1$id[[i]],  broom::tidy(ssfm1$mod1[[i]]$model)
    ))
  } else(print(ssfm1$id[[i]]))
}

#' 
#' Pull off coefficients for second model - no convergence warnings
coefssf2<-NULL
for(i in 1:3){
  if(attr(ssfm2$mod2[[i]], "class")!="try-error"){
    coefssf2<-rbind(coefssf2, cbind(id=ssfm2$id[[i]],  broom::tidy(ssfm2$mod2[[i]]$model)
    ))
  } else(print(ssfm2$id[[i]]))
}

#' Pull off coefficients for third model
coefssf3<-NULL
for(i in 1:3){
  if(attr(ssfm3$mod3[[i]], "class")!="try-error"){
    coefssf3<-rbind(coefssf3, cbind(id=ssfm3$id[[i]],  broom::tidy(ssfm3$mod3[[i]]$model)
    ))
  } else(print(ssfm3$id[[i]]))
}

#' 
#' 
#' Step 3: compare above to ctmc approach to understand differences/similarities
#' Step 4: estimate transition probabilities from the 2 models/compare simulated movements, etc