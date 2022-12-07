######################
# 1a - Process MODIS landcover (MCD12Q1) - create forest mask for greenup
#
######################


# Raw data DL instructions ------------------------------------------------

#go to: https://lpdaacsvc.cr.usgs.gov/appeears/task/area
#select MCD12Q1 -> all layers
#select GeoTIFF and geographic projection
#DL MCD12Q1
#select MCD12Q2 -> all layers
#select GeoTIFF and geographic projection
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
library(ggplot2)
library(rgdal)


# MCD12Q1 - Landcover -----------------------------------------------------

#LCCS1 Land Cover Layer
#from user guide:
#1 - barren
#2 - perm snow and ice
#3 - water bodies
#11 - evergreen needleleaf forests
#12 - evergreen broadleaf forests
#13 - deciduous needleleaf forests
#14 - deciduous broadleaf forests
#15 - mixed broadleaf/needleleaf forests
#16 - mixed broadleaf evergreen/deciduous forests
#21 - open forests
#22 - sparse forests
#31 - dense herbaceous
#32 - sparse herbaceous
#41 - dense shrublands
#42 - shrubland/grassland mosaics
#43 - sparse shrublands
#255 - unclassified


#read in  tif - 2019 LCCS1
lc1 <- terra::rast(paste0(dir, '/Data/environment/RAW/MCD12Q1/MCD12Q1.006_LC_Prop1_doy2019001_aid0001.tif'))


# create mask of just forest land cover type ------------------------------

lc_for_mask <- lc1 %in% c(11, 12, 13, 14, 15, 16, 21, 22)


# save to RDS -------------------------------------------------------------

#create output dir if it doesn't exist
ifelse(!dir.exists(paste0(dir, 'Data/environment/processed/', run_date)),
       dir.create(paste0(dir, 'Data/environment/processed/', run_date)),
       FALSE)

#saveRDS(lc_for_mask, paste0('lc_for_mask-', run_date,'.rds'))
terra::writeRaster(lc_for_mask, paste0(dir, 'Data/environment/processed/', run_date, '/lc_for_mask-', run_date, '.tif'))

