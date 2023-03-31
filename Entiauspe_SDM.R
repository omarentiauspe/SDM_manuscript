################################################################################
# Can fieldwork driven by predictive species distribution 
# models yield new rare or relevant geographic records? 
# A case study with Neotropical snakes.
# Entiauspe-Neto, O.M. et al.
################################################################################

# Load packages #

# setwd("~/Rdir") #

library(dismo)
library(maptools)
library(maps)    
library(mapdata) 
library(dplyr)
library(CoordinateCleaner)

# Create a new directory # 
# 
# if(!dir.exists("~/SDM_Work/Species_Loc")){
#   dir.create("~/SDM_Work/Species_Loc", recursive = TRUE)
# }
# 
# setwd("~/SDM_Work/Species_Loc")

# Add individual geographic distribution files in CSV # 

geospecies <- read.csv("A_dimidiata.csv", header = TRUE)

geospecies$lon <- geospecies$lon
geospecies$lat <- geospecies$lat

# Clean coordinates for duplicates # 

capgeo <- geospecies %>%
  select(Species, lat, lon) %>% 
  filter(complete.cases(.)) %>%  
  distinct() 

str(capgeo)

# Plot a geographic distribution map # 

map(
  'worldHires',
  xlim = c(min(capgeo$lon) - 5.0, max(capgeo$lon) + 5.0),
  ylim = c(min(capgeo$lat) - 5.0, max(capgeo$lat) + 5.0),
  fill = T,
  col = "light grey"
)

box()

points(capgeo$lon,
       capgeo$lat,
       col = "orange",
       pch = 20,
       cex = 0.7)

capc <- capgeo %>% 
  dplyr::select(Species, lat, lon)

write.csv(capc, "snake_locs.csv", row.names = FALSE)

# Extract bioclimatic variables # 

require(raster)
require(maps)
require(mapdata)

# Create new directory for environmental envelopes # 

# if (!dir.exists("~/SDM_Work/Env_Data")) {
#   dir.create("~/SDM_Work/Env_Data")   
# }
# 
# setwd("~/SDM_Work/Env_Data") 

# Extract worldclim data # 

env <- getData("worldclim", var="bio", res=2.5)

# Visualize Bioclim data # 

bioclim_names <-
  c(
    "Annual_Mean_Temp",
    "Mean_Diurnal_Range",
    "Isothermality",
    "Temp_Seasonality",
    "Max_Temp_Warmest Month",
    "Min_Temp_Coldest_Month",
    "Temp_Annual_Range",
    "Mean_Temp_Wettest_Quarter",
    "Mean_Temp_Driest_Quarter",
    "Mean_Temp_Warmest_Quarter",
    "Mean_Temp_Coldest_Quarter",
    "Annual_Precip",
    "Precip_Wettest_Month",
    "Precip_Driest_Month",
    "Precip_Seasonality",
    "Precip_Wettest_Quarter",
    "Precip_Driest_Quarter",
    "Precip_Warmest_Quarter",
    "Precip_Coldest_Quarter"
  )

names(env) <- bioclim_names

plot(env)

# Cropping raster files to extent of distribution #


capg <- read.csv("snake_locs.csv")

buff <- 3   

xmin <- min(capg$lon) - 5.0
xmax <- max(capg$lon) + 5.0
ymin <- min(capg$lat) - 5.0
ymax <- max(capg$lat) + 5.0 

e <- extent(xmin, xmax, ymin, ymax)

envcrop <- crop(env, e)

plot(envcrop[[1]], main = "Annual Mean Temperature")
map(
  'worldHires',
  xlim = c(xmin, xmax),
  ylim = c(ymin, ymax),
  fill = F,
  add = T
)

points(capg$lon, capg$lat, pch = "+")


ten_div <-
  c(1, 2, 5, 6, 7, 8, 9, 10, 11)  

for (layer in ten_div) {
  envcrop[[layer]] <- envcrop[[layer]] / 5
}

# Overwriting bioclimatic variables #

layers <- c(1, 2, 3, 4, 5, 6, 7, 8,
            9, 10, 11, 12, 13, 14,
            15, 16, 17, 18, 19)

writeRaster(
  stack(envcrop[[layers]]),
  paste0("bio", layers), 
  bylayer = TRUE,
  format = 'ascii',
  overwrite = T
)

# Presence - absence map #


library(raster)
library(maps)
library(mapdata)
library(ggplot2)
library(scales)
library(rgdal)

capout <- raster("Apostolepis_dimidiata.asc")
plot(capout)
points(capg$lon, capg$lat, pch="+")

xmin <- extent(capout)[1]
xmax <- extent(capout)[2]
ymin <- extent(capout)[3]
ymax <- extent(capout)[4]


n <- c(0, 0.262, 0, 0.262, 1, 1)
rclmat2 <- matrix(n, ncol = 3, byrow = TRUE)

rclmat2

rc2 <- reclassify(capout, rclmat2)
plot(rc2)
map(
  'worldHires',
  xlim = c(xmin, xmax),
  ylim = c(ymin, ymax),
  fill = F,
  add = T
)
points(capg$lon, capg$lat, pch = "+")
