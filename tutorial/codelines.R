# Recommended code lines for running well hydRopclim functions
# Change the input and output path and files C:/.../

library(hydRopclim)
database1 <- read.csv("C:/.../data.csv",header=TRUE, check.names = F, stringsAsFactors = F)
output <- "C:/.../output.csv"
pgridcorr(data=database1)

library(hydRopclim)
database2 <- read.csv("C:/.../data.csv",header=TRUE, check.names = F, stringsAsFactors = F)
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
data <- read.csv("C:/.../data.csv",header=TRUE, check.names = F, stringsAsFactors = F)
index <- data$Index
variables <- data[ ,3:ncol(data)]
rwin <- 11 #insert the sliding window in years
output1 <- "C:/.../indexes.csv"
output2 <- "C:/.../runcorr.csv"
output3 <- "C:/.../runcorr_format.csv"
p1 <- seasavg(p=index, start=9, win=6) # e.g.SONDJF
p2 <- seassum(p=variables, start=9, win=8) #e.g. SONDJFMA
indexcorrl(index.seas=p1, variable.seas=p2, rwin=11)

library(stats)
library(cluster)
library(sp)
library(rgdal)
library(hydRopclim)
database4 <- read.csv("C:/.../data.csv", header = T, check.names = F, stringsAsFactors = F)
region <- readOGR("C:/.../region.shp")
n <- 3
output <- "C:/.../output_test.csv"
hydrocluster(file=database4, shp=region, clusters=n)

library(hydRopclim)
database5 <- read.csv("C:/.../data.csv",header=TRUE, check.names = F, stringsAsFactors = F)
output <- "C:/.../output.csv"
lat <- -5.11 # Enter basin latitude
hydrochange(data=database5, lat)

library(hydRopclim)
library(ggplot2)
database7 <- read.csv("C:/.../data.csv",header=TRUE, check.names = F, stringsAsFactors = F)
a <-2352   #Area in km2
l <- 88.3  #Main channel lenght in km
p <- 261   #Perimeter in km
output <- "C:/.../output.csv"
rindex(data=database7, a, l, p)
