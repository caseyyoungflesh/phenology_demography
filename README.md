# phenology_demography


[![DOI](https://zenodo.org/badge/575464073.svg)](https://zenodo.org/badge/latestdoi/575464073)


Code for Youngflesh et al. In Press *PNAS*

This repository contains code to assess the demographic consequences of phenological dynamics in North American birds.

**Associated publications:**

Youngflesh, C, GA Montgomery, JF Saracco, DAW Miller, RP Guralnick, AH Hurlbert, RB Siegel, R LaFrance, MW Tingley. Demographic consequences of phenological asynchrony for North American songbirds. In Press at *__PNAS__*


**Repository structure:**
  * `Scripts/`
      * `1-process-gr/`
        * `1a-landcover.R` - produce forest landcover mask
        * `1b-greenup.R` - calculate greenup
      * `2-pro-PC.R` - productivity as a function of pheno indices
      * `3-juv-gr.R` - breeding phenology sensitivity to green-up
      * `4-sens-juv-PC/`
        * `4a-process-range-maps.R` - calculate migration distance
        * `4-sens-juv-pc.R` - calculate breeding sens as a function of species' traits
      * `5-predict.R` - predict demographic impact of phenological change
      * `Model_files/` - Stan files
        * `2-pro-pc.stan` - Stan file for model in `2-pro-PC.R`
        * `3-juv-gr-stan` - Stan file for model in `3-juv-gr.R`
        * `4-sens-juv-pc.stan` - Stan file for model in `4-sens-juv-pc.R`
  * `Data/` (ignored)
      * `environment/`
        * `RAW/`
          * `MCD12Q1/` - MODIS-derived landcover data
          * `MCD12Q2/` - MODIS-derived phenology data
        * `processed/` - processed env data
        * `climate/` - climate model data from `climatena.ca`
      * `MAPS_data/` - MAPS data
      * `BOTW/` - BirdLife range maps
      * `bird_phylo/` - data from `birdtree.org`
  * `Results/` (ignored)
