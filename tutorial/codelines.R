# Suggested code lines for running hydRopclim functions
# Change the input and output path and files C:/.../

library(hydRopclim)
database1 <- read.csv("C:/.../data.csv")
output <- "C:/.../output.csv"
pgridcorr(data=database1)

library(hydRopclim)
database2 <- read.csv("C:/.../data.csv")
hBS = 2479 #Station in masl
hNNR = 3012 #Grid.level in masl
hx = 2000 #The location to correct in masl
gradb <- c(-0.6, -0.6,	-0.6,	-0.6,	-0.6,	-0.6,	-0.6,	-0.6,	-0.6,	-0.6,	-0.6,	-0.6)
output <- "C:/.../output.csv"
tgridcorr(data=database2, hBS, hNNR, hx, LR=gradb)

library(reshape2)
library(ggplot2)
library(wesanderson)
library(cowplot)
library(hydRopclim)
data <- read.csv("C:/.../data.csv")
index <- data$Index
variables <- data[ ,3:ncol(data)]
rwin <- 11 #insert the sliding window in years
output1 <- "C:/.../indexes.csv"
output2 <- "C:/.../runcorr.csv"
output3 <- "C:/.../runcorr_format.csv"
p1 <- seasavg(p=index, start=9, win=6) # e.g.SONDJF
p2 <- seassum(p=variables, start=9, win=8) #e.g. SONDJFMA
indexcorrl(index.seas=p1, variable.seas=p2, rwin=11)

library(cluster)
library(sp)
library(terra)
library(hydRopclim)
database4 <- read.csv("C:/.../data.csv")
region <- vect("C:/.../region.shp")
n <- 3
output <- "C:/.../output_test.csv"
hydrocluster(file=database4, shp=region, clusters=n)

library(hydRopclim)
database5 <- read.csv("C:/.../data.csv")
output <- "C:/.../output.csv"
lat <- -5.11 # Enter basin latitude
hydrochange(data=database5, lat)

library(hydRopclim)
library(ggplot2)
database7 <- read.csv("C:/.../data.csv")
a <-2352   #Area in km2
l <- 88.3  #Main channel lenght in km
p <- 261   #Perimeter in km
output <- "C:/.../output.csv"
rindex(data=database7, a, l, p)

library(hydRopclim)
library(rgdal)
library(raster)
library(dplyr)
library(airGR)
library(parallel)
library(purrr)
library(lubridate)
library(zoo)
setwd("C:/.../")
source("https://raw.githubusercontent.com/hydrocodes/hydRopclim/main/tutorial/spatial_grad/interpolation.R")
#data upload
temp_stations <- read.csv('.../serie_tiempo_temp.csv')
prec_stations <- read.csv('.../serie_tiempo_prc.csv')
SRTM_0 <- raster('.../basin_90m.tif') #check dem_changes function below
ccas <- readOGR('.../basin.shp') #shp
grad_pr= 4*10^(-4)
grad_temp =  -6.5/1000
output <- ".../output_3basins.csv"
path_graphs <- ".../output/"


### A simple function for changing DEM m resolution to 0.1,0.25, 1, 5 or 10 km
dem_changes <- function(DEM_DATA, res_srtm, res_resultados) {
  if (res_srtm == 90) {
    if (res_resultados == 10) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(108, 108))
    } else if (res_resultados == 5) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(54, 54))
    } else if (res_resultados == 1) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(10.8, 10.8))
    } else if (res_resultados == 0.5) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(5.4, 5.4))
    } else if (res_resultados == 0.1) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(1.08, 1.08))
    }
    else {
      cat("Invalid result resolution. It must be 10km, 5km, 1km, 500m or 100m.")
    }
  } else if (res_srtm == 30) {
    if (res_resultados == 10) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(324, 324))
    } else if (res_resultados == 5) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(162, 162))
    } else if (res_resultados == 1) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(32.4, 32.4))
    } else if (res_resultados == 0.5) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(16.2, 162))
    } else if (res_resultados == 0.1) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(3.24, 3.24))
    } else {
      cat("Invalid result resolution. It must be 10km, 5km, 1km, 500m or 100m.")
    }
  } else if (res_srtm == 12.5) {
    if (res_resultados == 10) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(777.6, 777.6))
    } else if (res_resultados == 5) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(388.8, 388.8))
    } else if (res_resultados == 1) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(77.66, 77.66))
    } else if (res_resultados == 0.5) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(38.88, 38.88))
    } else if (res_resultados == 0.1) {
      DEM_DATA <- DEM_DATA %>% aggregate(c(7.766, 7.766))
    } else {
      cat("Invalid result resolution. It must be 10km, 5km, 1km, 500m or 100m.")
    }
  } else {
    cat("Invalid DEM resolution. It must be 12.5m, 30m or 90m.")
  }
  return(DEM_DATA)
}
SRTM_0 <- dem_changes(SRTM_0, 90, 10) # here need to change the dem resolution and resolution that you want
SRTM_0[SRTM_0[] < 0] <- 0
names(SRTM_0) <- 'Elevation'
ccas$NOMBRE <- 'Cuenca' #name of the basin (if you dont have "NOMBRE" or "NAME" in the .shp)

grad_spatial(SRTM_0, temp_stations, prec_stations, ccas, grad_temp, grad_pr)

library(hydRopclim)
database10 <- read.csv("C:/.../data.csv")
output <- "C:/.../output.csv"
rvm(database10, zeros = F)
