##########################
# 4a - Process migration distance for each species
#
# Using BirdLife Range maps
##########################


# set dirs -----------------------------------------------------------

dir <- 'XXXX'
run_date <- '2022-09-23'
pc_run_date <- '2022-04-08'


# load packages -----------------------------------------------------------

library(sf)


# load data ---------------------------------------------------------------

#MAPS data
tdata <- readRDS(paste0(dir, '/Results/pro-PC-', pc_run_date, 
                        '/pro-PC-data-', pc_run_date, '.rds'))$pro_data

#Birdlife range maps
BL_data <- sf::st_read(dsn = paste0(dir, 'Data/BOTW/BOTW.gdb/a00000009.gdbtable'))


# process ---------------------------------------------------------------

MAPS_usp <- unique(tdata$sci_name)
BL_usp <- unique(BL_data$sci_name)

#read in MAPS/BL key (to match species names)
MAPS_BL_key <- read.csv(paste0(dir, '/Data/BOTW/MAPS_BL_key.csv')

out <- data.frame(sci_name = rep(NA, length(usp)), dis = NA)

#loop through each species
for (k in 1:length(usp))
{
  #k <- 12
  print(paste0('Processing: ', k, ' of ', length(usp)))
  #filter by species
  if (MAPS_usp[k] %in% MAPS_BL_key$MAPS)
  {
    sp <- MAPS_BL_key$BL[which(MAPS_usp[k] %in% MAPS_BL_key$MAPS)]
  } else {
    if (MAPS_usp[k] %in% BL_usp)
    {
      sp <- MAPS_usp[k]
    } else {
      sp <- NA
    }
  }
  
  out$sci_name[k] <- MAPS_usp[k]
  
  if (!is.na(sp))
  {
    temp <- dplyr::filter(BL_data, sci_name == sp)
  
    #convert to equal area projection
    temp2 <- temp %>%
      st_transform(crs = "+proj=laea +lat_0=0 +lon_0=-70 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
  
    #filter by range
    #seasonal
    #1 = Resident
    #2 = Breeding season
    #3 = Non-breeding season
    #4 = Passage
    #5 = Uncertain
    
    #combine breeding and resident (or non-breeding and resident)
    temp_breeding <- dplyr::filter(temp2, seasonal %in% c(2,1))
    temp_non_breeding <- dplyr::filter(temp2, seasonal %in% c(3,1))
    
    if (NROW(temp_breeding) > 0 & NROW(temp_non_breeding) > 0)
    {
      if (st_geometry_type(temp_breeding) == 'MULTIPOLYGON' & 
          st_geometry_type(temp_non_breeding) == 'MULTIPOLYGON')
      {
        tb2 <- st_combine(temp_breeding)
        tnb2 <- st_combine(temp_non_breeding)
        
        #get centroid
        br_cen <- st_centroid(tb2)
        nb_cen <- st_centroid(tnb2)
        
        br_cen_4326 <- br_cen %>%
          st_transform(4326)
        nb_cen_4326 <- nb_cen %>%
          st_transform(4326)
        
        #distance between centroids (mig distance)
        dist <- as.numeric(st_distance(br_cen, nb_cen)) / 1000
        
        out$dis[k] <- round(dist, 0) #in km
      }
    }
  } 
}


# save to rds -------------------------------------------------------------

saveRDS(out, paste0(dir, 'Data/mig_distance.rds'))

