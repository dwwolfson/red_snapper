#' # Visualization of red snappers movement data
#' 
#' Description: Mapping red snappers movement data on chicken rock classification habitat with 10x10 resolution
#' 
#' Programmer: JV
#' 
#' _____________________________________________________________________________
#' ## Preamble
#' 
#' Load libraries
#' 
#+ library, warnings = "hide"
library(here)
library(ezknitr)
library(knitr)
library(devtools)
library(lubridate)
library(raster)
library(dplyr)
library(tidyverse)
library(sf)
library(ggplot2)

#' Clear environment and set seed
#' 
remove(list=ls())
set.seed(168846)

#' _____________________________________________________________________________
#' ## Load Data - David code for setting crs and df formatting
#' 
reefclass <- raster(here("data/Seabed Maps/geotiff/", "ChickenRock_Classification.tif")) # sand == 1, low == 2, medium == 3, high == 4

# aggregate the resolution of the raster from 1 meter by 1 meter to 10x10
ras<-aggregate(reefclass, fact=10, fun=modal)

# import red snapper data
df<-read_csv(here("data/red.snapper.locations.csv"))

# pull off coordinate reference system from land cover raster layer
crs_ras<-crs(ras)

# convert snapper points to spatial object
df<-df %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)
df<-st_transform(df, crs=crs_ras)

# pull off coordinates for later
df<-df %>% 
  mutate(x=st_coordinates(.)[,1],
         y=st_coordinates(.)[,2])

#' _____________________________________________________________________________
#' ## Plotting individual movement
#' 
#' ### Split df to summer and fall dates
#' 
#' 
ras_df <- as.data.frame(ras, xy = TRUE) #raster to data frame 10x10m resolution

df_may <- df %>% filter((datetime >= as.Date("2019-05-01 00:00:01") & datetime < as.Date("2019-08-31 23:59:59"))) %>% arrange(datetime)#summer

df_sep <- df %>% filter((datetime >= as.Date("2019-09-01 00:00:01") & datetime <= as.Date("2019-12-31 23:59:59"))) %>% arrange(datetime)#fall

#' ### Plot all
#' 
p_all <- ggplot() +
  geom_raster(data = ras_df , aes(x = x, y = y, fill = ChickenRock_Classification, alpha= 0.5)) +
  scale_fill_viridis_c() +
  #geom_tile(colour="grey20", aes(fill = ChickenRock_Classification)) + #remove geom_raster and fill argument if this is used
  coord_quickmap() +
  geom_point(data = df, aes(x=x, y=y), size = 0.05, color = "black", shape = 1) + geom_path(data = df, aes(group = num, x=x, y=y, color = as.factor(num)), alpha = 0.7, size =0.1) + theme(legend.position = "none") + facet_wrap(~num)

#' ### Plot movement data - May

p_may_aug <- ggplot() +
geom_raster(data = ras_df , aes(x = x, y = y, fill = ChickenRock_Classification, alpha= 0.5)) +
  scale_fill_viridis_c() +
  #geom_tile(colour="grey20", aes(fill = ChickenRock_Classification)) + #remove geom_raster and fill argument if this is used
  coord_quickmap() +
  geom_point(data = df_may, aes(x=x, y=y), size = 0.05, color = "black", shape = 1) + geom_path(data = df_may, aes(group = num, x=x, y=y, color = as.factor(num)), alpha = 0.7, size =0.1) + theme(legend.position = "none") + facet_wrap(~num)
  
#' ### Plot movement data - September

p_sep_dec <- ggplot() +
  geom_raster(data = ras_df , aes(x = x, y = y, fill = ChickenRock_Classification, alpha= 0.5)) +
  scale_fill_viridis_c() +
  #geom_tile(colour="grey20", aes(fill = ChickenRock_Classification)) +
  coord_quickmap() +
  geom_point(data = df_sep, aes(x=x, y=y), size = 0.05, color = "black", shape = 1) + geom_path(data = df_sep, aes(group = num, x=x, y=y, color = as.factor(num)), alpha = 0.7, size =0.1) + theme(legend.position = "none") + facet_wrap(~num)

#' _____________________________________________________________________________
#' ## Save Data
#' 
ggsave("output/figures/p_all.png", p_all, dpi = 320, width = 10, height = 8)
ggsave("output/figures/p_may_aug.png", p_may_aug, dpi = 320, width = 10, height = 8)
ggsave("output/figures/p_sep_dec.png", p_sep_dec, dpi = 320, width = 10, height = 7)

#' _____________________________________________________________________________
#' ## Footer
#' 
devtools::session_info()
#' spun with:
#' ezknitr::ezspin(file = "scripts/movement_plot.R", keep_md = FALSE, out_dir = "output", fig_dir = "figures")
