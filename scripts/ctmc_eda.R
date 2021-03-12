library(here)
library(raster)
library(sf)
library(tidyverse)
library(mgcv)
library(ctmcmove)
library(rgdal)



#' Read in the seabed maps and create a raster object for each point in the map
reefclass <- raster(here("data/Seabed Maps/geotiff/", "ChickenRock_Classification.tif"))

# sand == 1
# low == 2
# medium == 3
# high == 4

# aggregate the resolution of the raster from 1 meter by 1 meter to 10x10
ras<-aggregate(reefclass, fact=10, fun=modal)

# import red snapper data
df<-read_csv(here("data/red.snapper.locations.csv"))

# trans = transmitter number
# datetime = date and time of the acoustic detection
# julian = day of the year of acoustic detection
# lat = latitude of fish position (degrees N)
# lon = longitude of fish position (degrees W)
# depth = depth of fish from surface (m; note that the study area is approximately 38 m deep)
# num = unique fish number to be used to identify individuals in this study since one transmitter was used on two separate fish
# mov = movement rate of an individual fish between that and following location (m/s)
# dist = distance fish moved between that and following location (m)
# sec = time between that and following acoustic detection 

# pull off coordinate reference system from land cover raster layer
crs_ras<-crs(ras)

# convert snapper points to spatial object
df<-df %>% st_as_sf(coords = c("lon", "lat"), crs = 4326)

# transform points to crs of raster
df<-st_transform(df, crs=crs_ras)

# double check
st_crs(df)==st_crs(ras)

# create a raster layer for the intercept
int<-ras
values(int)<-1

stack_rast<-stack(int, ras)
names(stack_rast)<-cbind(c("int", names(ras)))

######

# pull off coordinates for later
df<-df %>% 
  mutate(x=st_coordinates(.)[,1],
         y=st_coordinates(.)[,2])

# First drop points that aren't in the raster because they make ctmc choke
df$ras_value<-raster::extract(ras, df)
df<-df %>% drop_na(ras_value)
# 152 points dropped (only 0.04%)

# make sure everything is ordered correctly and don't need geom column anymore (slows things down)
df<-df %>% arrange(trans, datetime)
df<-df %>% rename(id=trans) #trans is an object in a later function so better to change

ids<-as_factor(unique(df$id))
snapper_list<-list()



for(i in seq_along(ids)){
  tmp<-df %>% filter(id==ids[[i]])
  
  ## turn time into a numeric value (in days since Jan 1, 1970)
  t<-as.numeric(strptime(tmp$datetime, format="%Y-%m-%d %H:%M:%S"))/60/60/24
  
  ## get xy values for each time point
  xy<-tmp %>% dplyr::select(x,y)
  
  snapper_list[[i]]<-data.frame(t=t, x=xy$x, y=xy$y)
}

#######
# CTMC Path #

n_snap<-length(snapper_list)
ctmc_list<-list()

for(i in 1:n_snap){
  ctmc_list[[i]]<-path2ctmc(xy = snapper_list[[i]][,-1], 
                            t = as.numeric(snapper_list[[i]][,1]),
                            method="LinearInterp",
                            rast = stack_rast
                            )}


###############
# Turn CTMC path into poisson glm data
glm_list<-list()

for(i in 1:n_snap){
  glm_list[[i]]<-ctmc2glm(ctmc_list[[i]], 
                          stack.static = stack_rast,
                          stack.grad = ras)
}
  
str(glm_list[[1]])
  
# convert from list to dataframe
snap_dat<-glm_list[[1]]  
for(i in 2:n_snap){
  snap_dat<-rbind(snap_dat, glm_list[[i]])
}  
  
length(which(snap_dat$tau==0))
# only 2,436 tau out of 4M are 0

length(which(snap_dat$tau<10^-5))
# and only 23k are less than 0.00005

# I'll remove tau under 10^-5 because that is what Hanks et al and Brennan et al use
# We probably want to think this over though

snap_dat<-snap_dat %>% filter(tau>10^-5)

# Make sure habitat is categorical
snap_dat<-snap_dat %>% rename(habitat=ChickenRock_Classification)
snap_dat$habitat<-as.factor(snap_dat$habitat)

# Fit glm for all fish
m1<-glm(z~habitat+crw,
        offset=log(tau),
        family=poisson,
        data=snap_dat)

summary(m1)

# Should we also fit a glm to each fish separately?
