library(tidyverse)
library(here)

# get file names
dat_files<-list.files(path=here("processed_data"), pattern="*.Rdata", full.names = T)

# function for import and add column with grid size
read_func<-function(x){
  out<-readRDS(x)
  grid_size<-sapply(strsplit(x, "_"), function(x) x[8]) # pull off grid size ( Not Generic approach)
  cbind(grid_size=grid_size, out)
}  

# bring in data
df<-map_dfr(dat_files, read_func) %>% 
  as_tibble()

