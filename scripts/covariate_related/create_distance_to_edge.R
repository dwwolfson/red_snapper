#' # Create distance to edge raster values
#' 
#' Description: This is my attempt to create a 'distance to edge of sand/hard bottom' predictor variable for every pixel cell in the original seabed raster
#' 
#' Programmer: DaVid
#' 
#' _____________________________________________________________________________
#' ## Preamble
#' 
#' Load libraries
#' 
#+ library, warnings = "hide"
library(here)
library(sf)
library(raster)
library(tidyverse)
library(rgeos)

seabed<-st_read(here("data/Seabed Maps/Shapefile/NF_19_06_Segmentation.shp"))

plot(seabed["Class_name"])

# Merge the hardbottom categories
seabed<-seabed %>% 
  mutate(binary_hab=ifelse(Class_name%in%c("High", "Low", "Medium"), "hard-bottom", "sand"))

seabed<-seabed %>% 
  group_by(binary_hab) %>% 
  summarise(geometry=st_union(geometry))

plot(seabed)

# Convert from polygon (with area), to line (w/o area)
sea_lines<-st_cast(seabed, "MULTILINESTRING")
sea_iso<-st_cast(sea_lines, "LINESTRING")

#bring in raster 
ras<-raster(here("data/Seabed Maps/geotiff/", "ChickenRock_Classification.tif"))

# identify which cells in the raster are on the edge of sand and hardbottom
# bob<-raster::extract(ras, sea_iso) # this takes forever, although it doesn't immediately run out of memory

extent(sea_iso)
extent(sea_lines)
extent(seabed)
# all the same

extent(ras)
# aaaalmost exactly the same
extent(ras)<-extent(seabed)
# now the same


p<-as(ras, "SpatialPoints")
sp_iso<-as(sea_iso, "Spatial")
proj4string(p)
proj4string(sp_iso)

d<-gDistance(p, sp_iso, byid=T)
# Can't execute this line
# Error: cannot allocate vector of size 15.6 Gb