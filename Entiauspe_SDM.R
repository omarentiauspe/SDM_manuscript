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

snakgeo <- geospecies %>%
  select(Species, lat, lon) %>% 
  filter(complete.cases(.)) %>%  
  distinct() 

str(snakgeo)

# Plot a geographic distribution map # 

map(
  'worldHires',
  xlim = c(min(snakgeo$lon) - 5.0, max(snakgeo$lon) + 5.0),
  ylim = c(min(snakgeo$lat) - 5.0, max(snakgeo$lat) + 5.0),
  fill = T,
  col = "light grey"
)

box()

points(snakgeo$lon,
       snakgeo$lat,
       col = "orange",
       pch = 20,
       cex = 0.7)

sapc <- snakgeo %>% 
  dplyr::select(Species, lat, lon)

write.csv(sapc, "snake_locs.csv", row.names = FALSE)

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


sapg <- read.csv("snake_locs.csv")

buff <- 3   

xmin <- min(sapg$lon) - 5.0
xmax <- max(sapg$lon) + 5.0
ymin <- min(sapg$lat) - 5.0
ymax <- max(sapg$lat) + 5.0 

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

points(sapg$lon, sapg$lat, pch = "+")


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

snkout <- raster("Apostolepis_dimidiata.asc")
plot(snkout)
points(sapg$lon, sapg$lat, pch="+")

xmin <- extent(snkout)[1]
xmax <- extent(snkout)[2]
ymin <- extent(snkout)[3]
ymax <- extent(snkout)[4]


n <- c(0, 0.262, 0, 0.262, 1, 1)
rclmat2 <- matrix(n, ncol = 3, byrow = TRUE)

rclmat2

rc2 <- reclassify(snkout, rclmat2)
plot(rc2)
map(
  'worldHires',
  xlim = c(xmin, xmax),
  ylim = c(ymin, ymax),
  fill = F,
  add = T
)
points(sapg$lon, sapg$lat, pch = "+")


####
# Exploratory analyses with removing highly correlated variables
####
# This is optional and was NOT used in final analyses; comment these lines
# if you are not using this 
# Step 1: List and load ASC files

setwd("~/SDM_TEST/Species_Loc")

file_list <- list.files("~/SDM_TEST/Species_Loc",
                    	pattern = ".asc", full.names = TRUE)

# Stack ASC files into a raster stack and convert to data frame
bioFinal <- stack(file_list)
bioFinal <- as.data.frame(bioFinal)

# Step 2: Calculate VIF and identify highly correlated variables
vif.cor <- vifcor(bioFinal, th = 0.7)
variables_to_remove <- vif.cor@excluded  # Variables to exclude due to high correlation

# Step 3: Create a curated data frame excluding highly correlated variables
bioFinal_curated <- bioFinal[, !names(bioFinal) %in% variables_to_remove]

# Step 4: Save each curated variable as a GeoTIFF file
output_folder <- "~/SDM_TEST/Species_Loc"
if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE)

# Define a template raster for extent and resolution (optional)
template_raster <- raster(file_list[1])  # Using the first file as a template

# Loop through each curated variable
for (var in colnames(bioFinal_curated)) {
  # Create a raster layer with the same extent and resolution as the template
  layer <- raster(template_raster)
  values(layer) <- as.numeric(bioFinal_curated[[var]])
 
  # Define file path and save as GeoTIFF
  tiff_file <- paste0(output_folder, "/", var, ".tif")
  writeRaster(layer, filename = tiff_file, format = "GTiff", overwrite = TRUE)
}

# Step 5: Generate a detailed report of excluded variables and correlations
# Retrieve the correlation matrix and identify variable pairs with correlation > 0.7
cor_matrix <- vif.cor@corMatrix
correlated_pairs <- data.frame(Variable1 = character(), Variable2 = character(), Correlation = numeric())

for (i in seq_len(nrow(cor_matrix))) {
  for (j in seq(i + 1, ncol(cor_matrix))) {
	if (abs(cor_matrix[i, j]) > 0.7) {  # Use threshold as per vifcor settings
  	correlated_pairs <- rbind(correlated_pairs,
                            	data.frame(Variable1 = rownames(cor_matrix)[i],
                                       	Variable2 = colnames(cor_matrix)[j],
                                       	Correlation = cor_matrix[i, j]))
	}
  }
}

# Step 6: Write the report to a text file
report_file <- paste0(output_folder, "/correlation_report.txt")
writeLines("Highly Correlated Variable Pairs Removed:\n", con = report_file)
write.table(correlated_pairs, file = report_file, append = TRUE, row.names = FALSE, col.names = TRUE)

cat("Curated GeoTIFF variables and detailed correlation report saved successfully!")

# Export the correlation matrix as a CSV file
write.csv(cor_matrix, file = "~/SDM_TEST/Species_Loc/correlation_curated_matrix.csv", row.names = TRUE)
write.csv(variables_to_remove, file = "~/SDM_TEST/Species_Loc/correlation_raw_matrix.csv", row.names = TRUE)

###
# End of exploratory analyses
###

# Load required libraries
library(raster)
library(maps)
library(mapdata)

# Load raster
snkout <- raster("Apostolepis_dimidiata.asc")
plot(snkout)
points(sapg$lon, sapg$lat, pch="+")

# Reclassification matrix
n <- c(0, 0.262, 0, 0.262, 1, 1)
rclmat2 <- matrix(n, ncol = 3, byrow = TRUE)

# Reclassify the raster
rc2 <- reclassify(snkout, rclmat2)
plot(rc2)

# Map plotting
map('worldHires', xlim = c(extent(snkout)[1], extent(snkout)[2]), ylim = c(extent(snkout)[3], extent(snkout)[4]), fill = FALSE, add = TRUE)
points(sapg$lon, sapg$lat, pch = "+")

# Export the presence layer as an .asc file
writeRaster(rc2, "presence_layer.asc", format = "ascii", overwrite = TRUE)


####
# Alternative presence-absence raster option, with the 
# "Maximum training sensitivity plus specificity" threshold
####
# Load required libraries
library(raster)
library(maps)
library(mapdata)
library(sp)  # for spatial objects

# Load raster
snkout <- raster("Apostolepis_dimidiata.asc")
plot(snkout)
points(sapg$lon, sapg$lat, pch="+")

# Reclassification matrix
n <- c(0, 0.262, 0, 0.262, 1, 1)
rclmat2 <- matrix(n, ncol = 3, byrow = TRUE)

# Reclassify the raster
rc2 <- reclassify(snkout, rclmat2)
plot(rc2)

# Map plotting
map('worldHires', xlim = c(extent(snkout)[1], extent(snkout)[2]), ylim = c(extent(snkout)[3], extent(snkout)[4]), fill = FALSE, add = TRUE)
points(sapg$lon, sapg$lat, pch = "+")

# Convert presence points to spatial object
coordinates <- SpatialPoints(cbind(sapg$lon, sapg$lat), proj4string = crs(snkout))

# Step 1: Extract presence values
presence_values <- extract(snkout, coordinates)

# Step 2: Calculate maximum training sensitivity plus specificity threshold manually
threshold <- min(presence_values, na.rm = TRUE)

# Step 3: Reclassify raster using the threshold
binary_raster <- calc(snkout, fun = function(x) ifelse(x >= threshold, 1, 0))
plot(binary_raster)

# Step 4: Create observed and predicted classifications
# Extract binary values at presence points for comparison
observed <- ifelse(!is.na(presence_values), 1, 0)  # Assuming all presence points are observed presences
predicted <- ifelse(presence_values >= threshold, 1, 0)

# Step 5: Generate the confusion matrix
conf_matrix <- table(observed, predicted)

# Step 6: Calculate the Kappa statistic manually
# Kappa calculation based on confusion matrix
n <- sum(conf_matrix)  # total samples
p_o <- sum(diag(conf_matrix)) / n  # observed agreement
p_e <- sum(rowSums(conf_matrix) * colSums(conf_matrix)) / (n * n)  # expected agreement
kappa <- (p_o - p_e) / (1 - p_e)
print(kappa)

# Export the binary presence-absence layer as an .asc file
writeRaster(binary_raster, "binary_presence_layer.asc", format = "ascii", overwrite = TRUE)

