# calculate full 'distance to edge' raster
library(here)
library(sf)
library(raster)
library(tidyverse)
library(rgeos)
library(fasterize)


seabed<-st_read(here("data/Seabed Maps/Shapefile/NF_19_06_Segmentation.shp"))
crs_seabed<-st_crs(seabed)


 # plot(seabed["Class_name"])

# Merge the hardbottom categories
seabed<-seabed %>% 
  mutate(binary_hab=ifelse(Class_name%in%c("High", "Low", "Medium"), "hard-bottom", "sand"))

seabed<-seabed %>% 
  group_by(binary_hab) %>% 
  summarise(geometry=st_union(geometry))



# Convert from polygon (with area), to line (w/o area)
sea_lines<-st_cast(seabed, "MULTILINESTRING")

ras<- raster(here("data/Seabed Maps/geotiff/", "ChickenRock_Classification.tif"))


# Convert edge lines to sf polygon with a buffer and then turn to a raster
# fasterize only works with polygons, not lines
line_buf<-st_buffer(sea_lines, 1)
line_ras<-fasterize(line_buf, ras)
dist_ras<-distance(line_ras)

# save raster
writeRaster(dist_ras, here("processed_data/dist_to_edge_raster/dist_ras.tif"))

################################################
###############################################
# Merge sand with low hard-bottom
seabed<-seabed %>% 
  mutate(binary_hab=ifelse(Class_name%in%c("High", "Medium"), "med/high_hard-bottom", "low_hard-bottom/sand"))

seabed<-seabed %>% 
  group_by(binary_hab) %>% 
  summarise(geometry=st_union(geometry))



# Convert from polygon (with area), to line (w/o area)
sea_lines<-st_cast(seabed, "MULTILINESTRING")

ras<- raster(here("data/Seabed Maps/geotiff/", "ChickenRock_Classification.tif"))


# Convert edge lines to sf polygon with a buffer and then turn to a raster
# fasterize only works with polygons, not lines
line_buf<-st_buffer(sea_lines, 1)
line_ras<-fasterize(line_buf, ras)
dist_ras<-distance(line_ras)

# save raster
writeRaster(dist_ras, here("processed_data/dist_to_edge_raster/sand_and_low_vs_med_and_high/sand-low_vs_med-high.tif"))

