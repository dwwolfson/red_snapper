# This script is to aggregate the 1m grid size habitat layer into bigger cell sizes
# Output is saved in ~/processed_data/rasters_varying_cell_sizes/


library(here)
library(raster)
library(glue)

# read in raster
ras <- raster(here("data/Seabed Maps/geotiff/", "ChickenRock_Classification.tif"))

grid_sizes<-seq(5,50, by=5)

for(i in seq_along(grid_sizes)){
  tmp_ras<-aggregate(ras, fact=grid_sizes[[i]], fun = modal)
  writeRaster(tmp_ras, 
       file=paste(here("processed_data/rasters_varying_cell_sizes"), 
                  glue("{grid_sizes[[i]]}_m^2_raster.tif"), 
                  sep="/"))
}
