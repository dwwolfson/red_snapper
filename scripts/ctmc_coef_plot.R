#' Load libraries
#+ warning=FALSE, message=FALSE
library(shiny)
library(here)
library(dplyr)
library(broom)
library(dotwhisker)
library(ggplot2)

#' Read in the processed data from [list file here]
glm_5m <- readRDS(here("processed_data","glm_ready_3fish_5_m^2_res.Rdata"))
glm_10m <- readRDS(here("processed_data","glm_ready_3fish_10_m^2_res.Rdata"))
glm_15m <- readRDS(here("processed_data","glm_ready_3fish_15_m^2_res.Rdata"))
glm_20m <- readRDS(here("processed_data","glm_ready_3fish_20_m^2_res.Rdata"))
glm_25m <- readRDS(here("processed_data","glm_ready_3fish_25_m^2_res.Rdata"))
glm_30m <- readRDS(here("processed_data","glm_ready_3fish_30_m^2_res.Rdata"))
glm_35m <- readRDS(here("processed_data","glm_ready_3fish_35_m^2_res.Rdata"))
glm_40m <- readRDS(here("processed_data","glm_ready_3fish_40_m^2_res.Rdata"))
glm_45m <- readRDS(here("processed_data","glm_ready_3fish_45_m^2_res.Rdata"))
glm_50m <- readRDS(here("processed_data","glm_ready_3fish_50_m^2_res.Rdata"))

#' Fit models for each ID for each grid-cell size
by_snap_5m <- glm_5m %>%
  group_by(id) %>%                        # group data by transmission
  do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
  ungroup() %>% rename(model = id) %>% mutate(size = 5)

by_snap_10m <- glm_10m %>%
  group_by(id) %>%                        # group data by transmission
  do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
  ungroup() %>% rename(model = id) %>% mutate(size = 10)

by_snap_15m <- glm_15m %>%
  group_by(id) %>%                        # group data by transmission
  do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
  ungroup() %>% rename(model = id) %>% mutate(size = 15)

by_snap_20m <- glm_20m %>%
  group_by(id) %>%                        # group data by transmission
  do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
  ungroup() %>% rename(model = id) %>% mutate(size = 20)

by_snap_25m <- glm_25m %>%
  group_by(id) %>%                        # group data by transmission
  do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
  ungroup() %>% rename(model = id) %>% mutate(size = 25)

by_snap_30m <- glm_30m %>%
  group_by(id) %>%                        # group data by transmission
  do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
  ungroup() %>% rename(model = id) %>% mutate(size = 30)

by_snap_35m <- glm_35m %>%
  group_by(id) %>%                        # group data by transmission
  do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
  ungroup() %>% rename(model = id) %>% mutate(size = 35)

by_snap_40m <- glm_40m %>%
  group_by(id) %>%                        # group data by transmission
  do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
  ungroup() %>% rename(model = id) %>% mutate(size = 40)

by_snap_45m <- glm_45m %>%
  group_by(id) %>%                        # group data by transmission
  do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
  ungroup() %>% rename(model = id) %>% mutate(size = 45)

by_snap_50m <- glm_50m %>%
  group_by(id) %>%                        # group data by transmission
  do(tidy(glm(z~habitat+crw, offset=log(tau), family=poisson, data=.))) %>%         # run model on each group
  ungroup() %>% rename(model = id) %>% mutate(size = 50)

#' rbind all the grid size coefficient values
by_snap <- rbind(by_snap_5m, by_snap_10m, by_snap_15m, by_snap_20m, by_snap_25m, by_snap_30m, by_snap_35m, by_snap_40m, by_snap_45m, by_snap_50m)

by_snap$size = as.factor(by_snap$size)


#' Plot results
#+ width=800, height=550
ggplot(by_snap, aes(x=size, y=estimate, color=size))+
  geom_pointrange(aes(ymin=estimate-1.96*std.error, ymax=estimate+1.96*std.error),
                position = position_dodge2(width = 0.5, padding = 0.5))+
  geom_point()+
  theme_bw() + ylab("Coefficient Estimate") +
  geom_hline(yintercept = 0, colour = "grey60", linetype = 2) +
  ggtitle("CTMC GLM") +
  theme(plot.title = element_text(face = "bold"), 
        legend.position = "none",
        axis.text.x = element_text(angle = 60, hjust = 1)) + xlab("Snapper ID") +
        facet_wrap(model~term, scales="free_y", ncol=5)
