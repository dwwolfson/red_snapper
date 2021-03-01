# red_snapper
red snapper mvmt analysis collab

-Implement Hanks CTDS modeling approach

To do:
1) Decide on biologically appropriate scale for GPS and covariates (seafloor is 1 meter res)
    -Grid size should be big enough that the animals frequently makes multiple movements and stays inside 
    -Otherwise residence time doesn't accumulate and it's not very informative if it's always super low
    
2) Rescale raster to chosen resolution

3) Try ctmc::path2ctmc function first, if this doesn't work, the options are to use Buderman function model or Johnson crawl package to fit    continuous-time movement model 
   - crawl package looks to have better documentation and support
   
4) convert CTMC path to Poisson glm data  (ctmc::ctmc2glm)

5) fit glm