######################
# 1b - Process MODIS greenup (MCD12Q2)
#
######################


# Raw data DL instructions ------------------------------------------------

#go to: https://lpdaacsvc.cr.usgs.gov/appeears/task/area
#DL MCD12Q1
#DL MCD12Q2

#projection:
# Datum:
#   WGS84
# PROJ.4:
#   +proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs


# specify dir -------------------------------------------------------------

dir <- 'XXXX'
run_date <- '2022-03-10'


# load packages -----------------------------------------------------------

library(terra)
library(sp)
library(sf)


# Greenup -----------------------------------------------------------------

#read in MAPS stations
data <- readRDS(paste0(dir, 'Data/MAPS_data/MAPS-data.rds'))

#MAPS station lat/lon
sll <- unique(data[,c('station', 'lat', 'lng')])
pts_sp <- sp::SpatialPoints(cbind(sll$lng, sll$lat))
sf_pts <- sf::st_as_sf(pts_sp)
sf::st_crs(sf_pts) <- 4326

#reproject points add buffer around stations (10 km)
RADIUS <- 10000
sf_pts2 <- sf::st_transform(sf_pts, 
                            crs = sf::st_crs('+proj=aea +lat_0=23 +lon_0=-96 +lat_1=29.5 +lat_2=45.5 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs')) %>%
  sf::st_buffer(dist = RADIUS)

#read in forest mask
lc_for_mask <- terra::rast(paste0(dir, 'Data/environment/processed/', 
                                  run_date, '/lc_for_mask-', run_date, '.tif'))

#extract greenup for each year, filter by QA score and forest land cover type
#'MidGreenup' = mid (50%)
#forest = TRUE -> filter by forest land cover type
gr_pro_fun <- function(gr_type = 'Greenup', years = 2001:2019, forest = TRUE)
{
  out_df <- data.frame(station = rep(sll$station, length(2001:2019)),
                       year = rep(2001:2019, each = NROW(sll)),
                       lat = rep(sll$lat, length(2001:2019)),
                       lng = rep(sll$lng, length(2001:2019)),
                       gr_mn = NA,
                       n_pix = NA,
                       gr_type = gr_type)
  counter <- 1
  #process rasters
  for (i in 1:length(years))
  {
    #i <- 1
    print(paste0('Processing: ', years[i]))
    
    #QA overall
    qao <- terra::rast(paste0(dir, 'Data/environment/RAW/MCD12Q2/MCD12Q2.006_QA_Overall_0_doy', 
                              years[i], '001_aid0001.tif'))
    
    #number of veg pheno cycles
    nc <- terra::rast(paste0(dir, 'Data/environment/RAW/MCD12Q2/MCD12Q2.006_NumCycles_doy', 
                             years[i], '001_aid0001.tif'))
    
    #create mask based on QA scores and number of cycles
    #only 'best' and 'good'
    QA_mask <- qao %in% c(0, 1)
    #only cells with 1 cycle
    cy_mask <- terra::ifel(nc == 1, 1, 0)
    
    #greenup (onset, midup, or middown)
    gr <- terra::rast(paste0(dir, 'Data/environment/RAW/MCD12Q2/MCD12Q2.006_', 
                             gr_type, '_0_doy', years[i], '001_aid0001.tif'))
    
    #filter by QA score
    gr_f <- terra::mask(gr, QA_mask, maskvalue = 1,
                        inverse = TRUE)
    #filter by # cycles
    gr_f2 <- terra::mask(gr_f, cy_mask, maskvalue = 1,
                        inverse = TRUE)
    
    if (forest == TRUE)
    {
      #filter by forest mask
      gr_f3 <- terra::mask(gr_f2, lc_for_mask, maskvalue = 1,
                           inverse = TRUE)

    } else {
      gr_f3 <- gr_f2
    }
    
    #extract and take mean within polygon
    e1 <- terra::extract(gr_f3, vect(sf_pts2), fun = function(x) mean(x, na.rm = TRUE))
    
    #number of pixels
    c1 <- terra::extract(gr_f3, vect(sf_pts2), fun = function(x) sum(x > 0, na.rm = TRUE))
    
    jday <- as.numeric(julian(as.Date(e1[,2],
                                      origin = as.Date('1970-01-01')),
                              origin = as.Date(paste0(years[i], '-01-01'))))
    
    #fill df
    out_df$gr_mn[counter:(counter + NROW(sll) - 1)] <- jday
    out_df$n_pix[counter:(counter + NROW(sll) - 1)] <- c1[,2]
    
    counter <- (counter + NROW(sll))
    
    rm(gr)
    rm(QA_mask)
    rm(cy_mask)
    rm(gr_f)
    rm(gr_f2)
    rm(jday)
    rm(e1)
    rm(c1)
    gc()
  }
  return(out_df)
}


# run function ------------------------------------------------------------

tt <- proc.time()
#filter by forest cover
mid_gr_out_f <- gr_pro_fun(gr_type = 'MidGreenup', forest = TRUE)
proc.time() - tt


# save to RDS -------------------------------------------------------------

setwd(paste0(dir, 'Data/environment/processed/', run_date))
saveRDS(mid_gr_out_f, paste0('MidGreenup-', RADIUS/1000, 'km-', run_date, '-forest.rds'))
