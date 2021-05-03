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


#bring in snapper points
df<-read_csv(here("data/red.snapper.locations.csv"))

# convert snapper points to spatial object
df<-df %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

# transform points to crs of raster
df<-st_transform(df, crs=crs_seabed)

# double check
st_crs(df)==st_crs(sea_lines)


df$dist_edge<-NA

for(i in 1:nrow(df)){
  df$dist_edge[[i]] <- min(
    st_distance(sea_lines, df[i,]))
  print(i)
  }

# save output
st_write(df, here("processed_data/points_with_dist_edge/df.gpkg"))
write_csv(df, here("processed_data/points_with_dist_edge/df.csv"))
saveRDS(df, here("processed_data/points_with_dist_edge/df.rds"))
