#' @title pgridcorr
#' @description Correcting precipitation grid
#' @param data A monthly dataframe: %b-%Y(Date), in-situ data (Station) and grid data (Grid)
#' @return Corrected precipitation grid in function of an in-situ station, control quality, plots and an output file
#' @examples pgridcorr(data)
#' @export
pgridcorr <- function(data)
{ prec  <- data$Station
 grido <- data$Grid
 dates <- as.POSIXct(paste("01", data$Date, sep = "-"), format = "%d-%b-%Y")
 prec.matrix <- t(matrix(prec,12))
 grido.matrix <- t(matrix(grido,12))
 trans.station <- log10(colMeans(prec.matrix)+1)
 trans.grido <- log10(colMeans(grido.matrix)+1)
 f2 <- trans.grido/trans.station
 gms <- grido.matrix+1
 gridc.matrix <- matrix(NA,nrow(gms),12)
 for (i in 1:nrow(gms)) {
  gridc.matrix[i,] <- ((gms[i,])^(1/f2))-1
 }
 gridc <- matrix(t(gridc.matrix),1)
 gridcv <- as.vector(gridc)
# Time series plot
layout(matrix(c(1, 1, 2, 2), nrow = 1), widths = c(2, 1))
plot(dates,grido,col="blue",type="l", ann=F, axes=F, ylim=range(grido,prec,gridcv))
par(new=TRUE)
plot(dates,prec,col="black",type="l", ann=F, axes=F, ylim=range(grido,prec,gridcv))
par(new=TRUE)
plot(dates, gridcv, type="l",col="red", ylab="P (mm/month)", xlab='', ylim=range(grido,prec,gridcv))
legend("top",legend = c("Station", "Original grid", "Corrected Grid"),
             col = c("black","blue","red"), lty=1, cex=0.8,horiz=TRUE, xpd=TRUE, inset = -0.15, bty = "n")

plot(prec,grido,col="black", ann=F, axes=F, ylim=range(grido,prec,gridcv))
par(new=TRUE)
plot(prec,gridcv,col="red", xlab="Station (mm/month)", ylab="Grid (mm/month)", ylim=range(grido,prec,gridcv))
legend("top",legend = c("Original", "Corrected"),
             col = c("black","red"), lty=1, cex=0.8,horiz=TRUE, xpd=TRUE, inset = -0.15, bty = "n")

# Control quality: rmse(root-mean-square error %), cc(correlation coefficient), slope(from a linear regression)
 rmse <- sqrt(mean((prec - grido)^2))*100/mean(prec)
 rmse.corr <- sqrt(mean((prec - gridcv)^2))*100/mean(prec)
 cc <- cor(prec, grido)
 cc.corr <- cor(prec, gridcv)
 slope <- cc*sd(grido)/sd(prec)
 slope.corr <- cc.corr*sd(gridcv)/sd(prec)
 original <- data.frame(rmse,cc,slope)
 corrected <- data.frame(rmse.corr,cc.corr,slope.corr)
 data$Grid.corr <- gridcv
 write.csv(data,output)
 return(list(gridcv, original, corrected))
}

#' @title tgridcorr
#' @description Correcting temperature grid
#' @param data A monthly dataframe: %b-%Y(Date), in-situ data (Station) and grid data at a level pressure (Grid.level)
#' @param hBS A numeric value: station elevation in masl
#' @param hNNR A numeric value: grid level elevation in masl
#' @param hx A numeric value: Location elevation to correct in masl
#' @param LR A 12 elements vector: mean monthly lapse rate in °C/100m
#' @return Corrected temperature grid at a given elevation in function of an in-situ station, control quality, plots and an output file
#' @examples tgridcorr(data,hBS,hNNR,hx,LR)
#' @export
tgridcorr <- function(data,hBS,hNNR,hx,LR)
{ temp  <- data$Station
 grido <- data$Grid.level
 dates <- as.POSIXct(paste("01", data$Date, sep = "-"), format = "%d-%b-%Y")
 temp.matrix <- t(matrix(temp,12))
 grido.matrix <- t(matrix(grido,12))
 mean.temp <- colMeans(temp.matrix)
 mean.grido <- colMeans(grido.matrix)
 if (hBS<hNNR)
 { f1 <- (LR/100)*abs(hBS-hNNR)
 } else {
   f1 <- (-LR/100)*abs(hBS-hNNR)
 }
 if (hx>hNNR)
 { f2 <- (LR/100)*abs(hx-hNNR)
 } else {
   f2 <- (-LR/100)*abs(hx-hNNR)
 }
 f <- (mean.temp+f1+f2)/mean.grido
 gridc.matrix <- matrix(NA,nrow(grido.matrix),12)
 for (i in 1:nrow(grido.matrix)) {
   gridc.matrix[i,] <- (grido.matrix[i,])*f
 }
 gridc <- matrix(t(gridc.matrix),1)
 gridcv <- as.vector(gridc)

# Matching test
 if (hBS<hNNR)
 { fto <- (-LR/100)*abs(hBS-hNNR)
 } else {
   fto <- (LR/100)*abs(hBS-hNNR)
 }
 ft <- (mean.grido+fto)/mean.temp
 temptc.matrix <- matrix(NA,nrow(temp.matrix),12)
 for (i in 1:nrow(temp.matrix)) {
   temptc.matrix[i,] <- (temp.matrix[i,])*ft
 }
 temptc <- matrix(t(temptc.matrix),1)
 temptcv <- as.vector(temptc)

# Time series plot
layout(matrix(c(1, 1, 2, 2), nrow = 1), widths = c(1,1,2,2))
plot(temp,grido,col="black", ann=F, axes=F, ylim=range(grido,temp,temptcv))
par(new=TRUE)
plot(temp,temptcv,col="green", xlab="Station (°C)", ylab="Grid (°C)", ylim=range(grido,temp,temptcv))
legend("top",legend = c("Original", "Corrected"),
             col = c("black","green"), lty=1, cex=0.8,horiz=TRUE, xpd=TRUE, inset = -0.15, bty = "n")

plot(dates,grido,col="blue",type="l", ann=F, axes=F, ylim=range(grido,temp,gridcv))
par(new=TRUE)
plot(dates,temp,col="black",type="l", ann=F, axes=F, ylim=range(grido,temp,gridcv))
par(new=TRUE)
plot(dates, gridcv, type="l",col="red", ylab="T (°C)", xlab='', ylim=range(grido,temp,gridcv))
legend("top",legend = c("Station", "Original grid", "Target elev"),
             col = c("black","blue","red"), lty=1, cex=0.8,horiz=TRUE, xpd=TRUE, inset = -0.15, bty = "n")

# Control quality: rmse(root-mean-square error %), cc(correlation coefficient)
 rmse <- sqrt(mean((temp - grido)^2))*100/mean(temp)
 rmse.corr <- sqrt(mean((temp - temptcv)^2))*100/mean(temp)
 cc <- cor(temp, grido)
 cc.corr <- cor(temp, temptcv)
 original <- data.frame(rmse,cc)
 corrected <- data.frame(rmse.corr,cc.corr)
 data$Grid.level.corr <- gridcv
 write.csv(data,output)
 return(list(gridcv,original,corrected))
}

#' @title indexcorrl
#' @description Running correlation between a climate index versus hydrological indexes
#' @param index.seas A vector object: a yearly climate index
#' @param variable.seas A matrix object: a yearly hydrological indexes
#' @param rwin A numeric value: sliding window length in years
#' @return Running correlation, linear trend slopes, plots and 3 outputs files
#' @examples indexcorrl(index.seas,variable.seas,rwin)
#' @export
indexcorrl <- function(index.seas, variable.seas, rwin)
{ ptr <- length(index.seas)
  ptnd <- nrow(variable.seas)
  if(ptr>ptnd){
    index.seas <- index.seas[-ptr]
  } else if(ptr<ptnd){
      variable.seas <- variable.seas[-1,]
  }
 rcorr <- matrix(NA,(nrow(variable.seas)-rwin),ncol(variable.seas))
 for (i in 1:(nrow(variable.seas)-rwin)) {
  rcorr[i,] <- cor(variable.seas[i:(i+rwin-1),],index.seas[i:(i+rwin-1)])
}
# Formatting new yearly dates
data$Date <- as.POSIXct(paste("01", data$Date, sep = "-"), format = "%d-%b-%Y")
s <- substring(data$Date,1,4)
sn <- as.numeric(s[1])
sno <- round(mean(c(sn,sn+rwin)), digits=0)
fechas <- sno:(sno+nrow(rcorr)-1)

# Preparing the plot
rownames(rcorr) <- fechas
colnames(rcorr) <- colnames(variables)
rcorrplot <- rcorr[,order(ncol(rcorr):1)]
cc <- vector("numeric", ncol(rcorrplot))
slope <- vector("numeric", ncol(rcorrplot))
for (i in 1:ncol(rcorr)) {
  cc[i] <- cor(rcorr[,i], fechas)
  slope[i] <- cc[i]*sd(rcorr[,i])/sd(fechas)
}
slope
df <- data.frame(rev(slope),1:ncol(rcorr))

reord <- melt(rcorrplot)
reord <- reord[reord$value!=0,]
pal <- wes_palette("Zissou1", type = "continuous")
p1 <- ggplot(reord, aes(x=factor(Var1), y=Var2, fill = value))+
  geom_tile()+
  scale_fill_gradientn(colours = pal, name = "r") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle=90))
p2 <- ggplot(df, aes(x=X1.ncol.rcorr., y= rev.slope.)) +
  theme_light()+ geom_line()+geom_point() +
  labs(y = "Linear trend slope") +
  scale_x_discrete(breaks=c(1:ncol(rcorr))) +
  theme(axis.title.y=element_blank())+ coord_flip()
p3 <- plot_grid(p1, p2, ncol=2, rel_widths = c(5,1))
idx <- merge(index.seas, variable.seas, by = 0, sort=F, all=T)
write.csv(idx,file=output1)
write.csv(rcorr,file=output2)
write.csv(reord,file=output3)
return(list(index.seas, variable.seas, rcorr, df, p3))
}

#' @title seasavg
#' @description Calculation of seasonal average vector index for a season of n-months
#' @param p A vector object: a monthly climate index
#' @param start A numeric value: season starting month (e.g. September=9)
#' @param win A numeric value: season length in months (e.g. 3 defines September-October-November)
#' @return Seasonal average index vector
#' @examples seasavg(p,start,win)
#' @export
seasavg <- function (p,start,win)
{ nyr <- NROW(p)/12
q <- vector("numeric", nyr)
if((start+win) <= 13)
{ q[1] <- mean(p[start:(start+win-1)])
for (j in 1:(nyr-1)) {
  q[j+1] <- mean(p[(start+12*j):(start+12*j+win-1)])
}
r <- q

} else
{ q[1] <- mean(p[start:(start+win-1)])
for (j in 1:(nyr-1)) {
  q[j+1] <- mean(p[(start+12*j):(start+12*j+win-1)])
}
r <- q[1:(length(q)-1)]
}
return(r)
}

#' @title seasavg2
#' @description Calculation of seasonal average matrix indexes for a season of n-months
#' @param p A matrix object: monthly hydroclimate indexes
#' @param start A numeric value: season starting month (e.g. September=9)
#' @param win A numeric value: season length in months (e.g. 3 defines September-October-November)
#' @return Seasonal average indexes matrix
#' @examples seasavg2(p,start,win)
#' @export
seasavg2 <- function(p,start,win)
{ nyr <- nrow(p)/12
q <- matrix(NA,nyr,ncol(p))
for(i in 1:ncol(p))
{ if((start+win) <= 13)
{ for (j in 0:(nyr-1)) {
  q[j+1,i] <- mean(p[(start+12*j):(start+12*j+win-1),i])
}
  r <- q
} else
{  for (j in 0:(nyr-1)) {
  q[j+1,i] <- mean(p[(start+12*j):(start+12*j+win-1),i])
}
  r <- q[-(nrow(q)),]
}
}
return(r)
}

#' @title seassum
#' @description Calculation of seasonal sum matrix indexes for a season of n-months
#' @param p A matrix object: monthly hydroclimate indexes
#' @param start A numeric value: season starting month (e.g. September=9)
#' @param win A numeric value: season length in months (e.g. 3 defines September-October-November)
#' @return Seasonal sum indexes matrix
#' @examples seassum(p,start,win)
#' @export
seassum <- function(p,start,win)
{ nyr <- nrow(p)/12
q <- matrix(NA,nyr,ncol(p))
for(i in 1:ncol(p))
{ if((start+win) <= 13)
{ for (j in 0:(nyr-1)) {
  q[j+1,i] <- sum(p[(start+12*j):(start+12*j+win-1),i])
}
  r <- q
} else
{  for (j in 0:(nyr-1)) {
  q[j+1,i] <- sum(p[(start+12*j):(start+12*j+win-1),i])
}
  r <- q[-(nrow(q)),]
}
}
return(r)
}

#' @title zscorem
#' @description Transformation of monthly m-hydroclimatic variables into m-zscores indexes over 12-months each one
#' @param data A matrix object: monthly hydroclimatic variables
#' @return Z-scores indexes over 12-months by each hydroclimatic variable
#' @examples zscorem(data)
#' @export
zscorem <- function(data)
{values <- data[,2:ncol(data)]
nyr <- nrow(values)/12
s <- matrix(NA,12,ncol(values))
t <- matrix(NA,12,ncol(values))
for(i in 1:12){
  s[i,] <- colMeans(values[seq(i, nrow(values), 12), ])
  t[i,] <- apply(values[seq(i, nrow(values), 12), ],2,sd)
}
list.val <- list()
list.fin <- list()
for(i in 1:ncol(values))
{list.val[[i]] <- t(matrix(values[,i:i],12,))
list.fin[[i]] <- as.vector((t(list.val[[i]])-s[,i])/t[,i])
}
write.table(list.fin, col.names=colnames(values), row.names=F, file=output, sep = ",")
return(list.fin)
}

#' @title hydrocluster
#' @description K-means clustering of stations hydroclimatic time series
#' @param file A dataframe object: head station names, latitude, longitude and monthly or annual data
#' @param shp A spatial object: a polygon shapefile of study region
#' @param clusters A numeric value: number of clusters
#' @return K-means clustering spatial visualization, control quality by silhouette evaluation and output file
#' @examples hydrocluster(file,shp,clusters)
#' @export
hydrocluster <- function(file,shp,clusters)
{ data <- file[3:nrow(file),2:ncol(file)]
coords.df <- data.frame(t(file[1:2,2:ncol(file)]))
colnames(coords.df) <- file[1:2,1]
data_t<-t(data)
# K-means and silhouette evaluation
km <- kmeans(data_t,centers=clusters, iter.max=100, nstart=2)
bdd <- data.frame(km$cluster,coords.df)
dissE <- daisy(data_t)
dE2   <- dissE^2
sk2   <- silhouette(km$cl, dE2)
group <- as.factor(bdd$km.cluster)
# Reading and georeferencing stations points and shapefile region
d <- data.frame(lon=bdd$long, lat=bdd$lat)
sp::coordinates(d) <- c("lon", "lat")
sp::proj4string(d) <- sp::CRS('+proj=longlat +datum=WGS84')
shpgeo <- project(shp, crs('+proj=longlat +datum=WGS84'))

# Plot of location map and silhouette
par(mfrow=c(1,2))
sp::plot(shpgeo, col="grey", border=NA, axes=T, xlab="Longitude", ylab="Latitude")
sp::plot(d, axes=T,pch=16, col=group, add=T)
text(bdd$long,bdd$lat, row.names(bdd), cex=0.6, pos=3,col="black")
legend("right", legend=levels(group), col=levels((group)),
             pch=16, inset=0.1, cex=0.8, xpd=T, bty = "n", x.intersp=0.1, y.intersp=0.4)
plot(sk2, main="")

average.sil <- mean(sk2[,"sil_width"])
neg.sil <- sum(sk2[,"sil_width"]< 0)
p <- as.data.frame(data_t)
colnames(p) <- file[3:nrow(file),1]
q <- cbind(bdd,p)
ord.cluster <- q[order(q$km.cluster),]
write.csv(t(ord.cluster),output)
return(list(paste0("Average Silhouette: ", round(average.sil,4)), paste0("Negative Silhouettes: ", neg.sil), km$cluster))
}

#' @title hydrochange
#' @description Hydroclimatic change analysis at annual time step for a database that includes Mean Temperature
#' @param data An annual dataframe object: %Y(Date), Precipitation (P in mm), Mean temperature (Tm in ?C) and Runoff (R in mm)
#' @param lat A numerical value: watershed latitude
#' @return Estimation of potential evapotranspiration by Oudin method, simulation of actual evapotranspiration by Budyko-Zhang model, quantifying impacts of climate and human activities on runoff change and quantifying watershed sensitivity and adaptation, plots and an output file
#' @examples hydrochange(data,lat)
#' @export
hydrochange <- function(data,lat)
{ P <- data$P
Tm <- data$Tm
Rm <- data$R
AET <- P-Rm
# Annual PET by Oudin method
J <- 1:365
delta <- 0.409 * sin(0.0172 * J - 1.39)
dr <- 1 + 0.033 * cos(0.0172 * J)
latr <- lat/57.2957795
sset <- -tan(latr) * tan(delta)
omegas <- sset * 0
omegas[sset>={-1} & sset<=1] <- acos(sset[sset>={-1} & sset<=1])
omegas[sset<{-1}] <- max(omegas)
Ra <- 37.6 * dr * (omegas * sin(latr) * sin(delta) + cos(latr) * cos(delta) * sin(omegas))
Ra <- ifelse(Ra < 0, 0, Ra)
Rann <- mean(Ra)*12/2.45
for (i in 1:length(Tm)) {
  if(Tm[i]+5>0)
  { PET <- Rann*(Tm+5)*30/100
  } else {0}
}
x <- PET/P
y <- AET/P
df = data.frame(x,y)
# Budyko-Zhang coefficient (wZ)
m1 <- nls( y ~ (1+wZ*x)/(1+wZ*x+1/x), data=df, start=list(wZ=0.001), trace=T)
summary(m1)
wZvalue <- summary(m1)$coefficients[1,1]; wZvalue
rse <- 100*summary(m1)$sigma; rse
if(wZvalue<0)
{print("w<0, careful! results could not be meaningful")
}
AETsim <- (P+PET*wZvalue)/(1+wZvalue*PET/P+P/PET)
# correlation coefficient
RSS.p1 <- sum(residuals(m1)^2)
TSS <- sum(y-mean(y)^2)
ccZ <- 1-(RSS.p1/TSS); ccZ
# Variables annual changes
year <- data$Date
prec.lm = lm(P ~ year)
pet.lm = lm(PET ~ year)
q.lm = lm(Rm ~ year)
aet.lm = lm(AET ~ year)
deltaP <- summary(prec.lm)$coefficients[2,1]
deltaPET <- summary(pet.lm)$coefficients[2,1]
deltaQ <- summary(q.lm)$coefficients[2,1]
deltaAET <- summary(aet.lm)$coefficients[2,1]

# Quantifying climate and human impacts on runoff
alpha <- (1+2*mean(x)+3*wZvalue*mean(x))/(1+mean(x)+wZvalue*(mean(x))^2)^2
betha <- -(1+2*wZvalue*mean(x))/(1+mean(x)+wZvalue*(mean(x))^2)^2
deltaQclim <- alpha*deltaP + betha*deltaPET
deltaQh <- deltaQ - abs(deltaQclim)
Clim <- deltaQclim*100/deltaQ; Clim
Hum <- deltaQh*100/deltaQ; Hum
# Quantifying basin adaptation and sensitivity
Adapt <- 90+atan(deltaAET*mean(P)-mean(AET)*deltaP)/(deltaPET*mean(P)-mean(PET)*deltaP)*180/pi
Sensitv <- sqrt( ((deltaAET*mean(P)-mean(AET)*deltaP)/mean(P)^2)^2 + ((deltaPET*mean(P)-mean(PET)*deltaP)/mean(P)^2)^2 )

# Creating plots
layout(matrix(c(1,1,2,2,1,1,3,3), nrow = 4), widths = c(1, 1))
print(plot(year,P, type="l", col="blue", ann=F, axes=F, ylim=range(P,PET,Rm)))
par(new=TRUE)
print(plot(year,PET,col="red",type="l", ann=F, axes=F, ylim=range(P,PET,Rm)))
par(new=TRUE)
print(plot(year, Rm, type="l",col="black", ylab="Annual (mm)", xlab='',ylim=range(P,PET,Rm)))
print(legend("top",legend = c("P", "PET", "Runoff"),
             col = c("blue","red","black"), lty=1, cex=0.8,horiz=TRUE, xpd=TRUE, inset = -0.6, bty = "n"))
print(plot(x,y,xlab="Dryness index PET/P", ylab="Evaporative index AET/P"))
s <- seq(from=0, to=10, length=50)
print(lines(s,predict(m1,list(x=s)), col="red"))
print(wformat <- format(wZvalue,digits=3L))
print(legend("top",legend=parse(text=sprintf('w == %s',wformat)),
             col = c("red"), lty=1, cex=0.9,horiz=TRUE, xpd=TRUE, inset = -0.6, bty = "n"))
print(plot(AET,AETsim,xlab="AET (mm)", ylab="Simulated AET (mm)"))

quality <- data.frame(wZvalue,ccZ,rse)
impacts <- data.frame(Clim,Hum)
trajectories <- data.frame(Adapt,Sensitv)
data$PET <- PET
data$AET <- AET
data$AETsim <- AETsim
write.csv(data,file=output)
return(list(PET, AETsim, quality, impacts, trajectories))
}

#' @title hydrochange2
#' @description Hydroclimatic change analysis at annual time step for a database that includes Potential Evapotranspiration
#' @param data An annual dataframe object: %Y(Date), Precipitation (P in mm), Potential evapotranspiration (PET in mm) and Runoff (R in mm)
#' @param lat A numerical value: watershed latitude
#' @return Simulation of actual evapotranspiration by Budyko-Zhang model, quantifying impacts of climate and human activities on runoff change and quantifying watershed sensitivity and adaptation, plots and an output file
#' @examples hydrochange2(data,lat)
#' @export
hydrochange2 <- function(data,lat)
{ P <- data$P
PET <- data$PET
Rm <- data$R
AET <- P-Rm
x <- PET/P
y <- AET/P
df = data.frame(x,y)
# Budyko-Zhang coefficient (wZ)
m1 <- nls( y ~ (1+wZ*x)/(1+wZ*x+1/x), data=df, start=list(wZ=0.001), trace=T)
wZvalue <- summary(m1)$coefficients[1,1]
rse <- 100*summary(m1)$sigma
if(wZvalue<0)
{print("w<0, careful! results could not be meaningful")
}
AETsim <- (P+PET*wZvalue)/(1+wZvalue*PET/P+P/PET)
# correlation coefficient
RSS.p1 <- sum(residuals(m1)^2)
TSS <- sum(y-mean(y)^2)
ccZ <- 1-(RSS.p1/TSS)
# Variables annual changes
year <- data$Date
prec.lm = lm(P ~ year)
pet.lm = lm(PET ~ year)
q.lm = lm(Rm ~ year)
aet.lm = lm(AET ~ year)
deltaP <- summary(prec.lm)$coefficients[2,1]
deltaPET <- summary(pet.lm)$coefficients[2,1]
deltaQ <- summary(q.lm)$coefficients[2,1]
deltaAET <- summary(aet.lm)$coefficients[2,1]

# Quantifying climate and human impacts on runoff
alpha <- (1+2*mean(x)+3*wZvalue*mean(x))/(1+mean(x)+wZvalue*(mean(x))^2)^2
betha <- -(1+2*wZvalue*mean(x))/(1+mean(x)+wZvalue*(mean(x))^2)^2
deltaQclim <- abs(alpha*deltaP + betha*deltaPET)
deltaQh <- deltaQ - deltaQclim
Clim <- deltaQclim*100/deltaQ; Clim
Hum <- deltaQh*100/deltaQ; Hum
# Quantifying basin adaptation and sensitivity
Adapt <- 90+atan(deltaAET*mean(P)-mean(AET)*deltaP)/(deltaPET*mean(P)-mean(PET)*deltaP)*180/pi
Sensitv <- sqrt( ((deltaAET*mean(P)-mean(AET)*deltaP)/mean(P)^2)^2 + ((deltaPET*mean(P)-mean(PET)*deltaP)/mean(P)^2)^2 )

# Creating plots
layout(matrix(c(1,1,2,2,1,1,3,3), nrow = 4), widths = c(1, 1))
plot(year,P, type="l", col="blue", ann=F, axes=F, ylim=range(P,PET,Rm))
par(new=TRUE)
plot(year,PET,col="red",type="l", ann=F, axes=F, ylim=range(P,PET,Rm))
par(new=TRUE)
plot(year, Rm, type="l",col="black", ylab="Annual (mm)", xlab='',ylim=range(P,PET,Rm))
legend("top",legend = c("P", "PET", "Runoff"),
       col = c("blue","red","black"), lty=1, cex=0.8,horiz=TRUE, xpd=TRUE, inset = -0.6, bty = "n")
plot(x,y,xlab="Dryness index PET/P", ylab="Evaporative index AET/P")
s <- seq(from=0, to=10, length=50)
lines(s,predict(m1,list(x=s)), col="red")
wformat <- format(wZvalue,digits=3L)
legend("top",legend=parse(text=sprintf('w == %s',wformat)),
       col = c("red"), lty=1, cex=0.9,horiz=TRUE, xpd=TRUE, inset = -0.6, bty = "n")
plot(AET,AETsim,xlab="AET (mm)", ylab="Simulated AET (mm)")

quality <- data.frame(wZvalue,ccZ,rse)
impacts <- data.frame(Clim,Hum)
trajectories <- data.frame(Adapt,Sensitv)
data$AET <- AET
data$AETsim <- AETsim
write.csv(data,file=output)
return(list(AETsim, quality, impacts, trajectories))
}

#' @title rindex
#' @description Estimation of runoff index in an ungauged basin
#' @param data A monthly dataframe object: %b-%Y(Date), Precipitation (P) and Potential Evapotranspiration (PET) in mm
#' @param a A numeric value: Watershed area in km2
#' @param l A numeric value: Watershed main channel lenght in km
#' @param p A numeric value: Watershed perimeter in km2
#' @return Runoff index time series, plot and output file
#' @examples rindex(data,a,l,p)
#' @export
rindex <- function(data,a,l,p)
{ X1 <- (a^0.393*l^-4.107*p^4.291)/64.5
X2 <- log(0.883*a^0.369*l^-0.229*p^-0.168)
X1o <- X1/2
X2o <- 30
phi <- tanh(data$P/X1)
psi <- tanh(data$PET/X1)
phio <- tanh(data$P[1]/X1)
psio <- tanh(data$PET[1]/X1)

S1 <- c(1:nrow(data))
S1[1] <- (X1o+X1*phio)/(1+phio*X1o/X1)
for (i in 2:nrow(data)) {
  S2 <- S1*(1-psi)/(1+psi*(1-S1/X1))
  S <- S2/((1+(S2/X1)^3)^(1/3))
  S1[i] <- (S[i-1]+X1*phi[i])/(1+phi[i]*S[i-1]/X1)
  S2 <- S1*(1-psi)/(1+psi*(1-S1/X1))
  S <- S2/((1+(S2/X1)^3)^(1/3))
}
P1 <- c(1:nrow(data))
P1[1] <- data$P[1]+X1o-S1[1]
for (i in 2:nrow(data)) {
  P1[i] <- data$P[i]+S[i-1]-S1[i]
  P2 <- S2-S
  P3 <- P1+P2
}
R1 <- c(1:nrow(data))
R1[1] <- P3[1]+X2o
for (i in 2:nrow(data)) {
  R2 <- R1*X2
  Q <- R2^2/(R2+60)
  R <- R2-Q
  R1[i] <- P3[i] + R[i-1]
  R2 <- R1*X2
  Q <- R2^2/(R2+60)
  R <- R2-Q
}
data$Qsim <- Q
Qmatrix <- t(matrix(Q,12))

res <- scale(Qmatrix)
attributes(res)[3:4] <- NULL
res_vector<-as.vector(t(res))
data$Date <- as.POSIXct(paste("01", data$Date, sep = "-"), format = "%d-%b-%Y")
data$RIndex <- res_vector

# Plotting time series
theme_set(theme_bw())
plotrindex <- ggplot(data, aes(x = Date, y = RIndex, fill = RIndex >=0)) +
  geom_col(position = "identity")+
  scale_fill_manual(values = c("red", "blue"), guide = "none")+
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=12),
        plot.title=element_text(size=14))
write.csv(data,output)
return(list(res_vector, Q, plotrindex))
}

#' @title spatial_grad
#' @description Areal precipitation, temperature and potential evapotranspiration
#' @param DEM A digital elevation model in masl
#' @param temp_stations A dataframe object: head station names, long in degrees, lat in degrees, elevation in masl and monthly temperature data with dates %b-%Y
#' @param prec_stations A dataframe object: head station names, long in degrees, lat in degrees, elevation in masl and monthly precipitation data with dates %b-%Y
#' @param ccas A shapefile object containing polygons representing basins 
#' @param grad_temp A numerical value of temperature gradient (C/m) 
#' @param grad_pr A numerical value of precipitation gradiente (mm/m) 
#' @return Corrected mean areal precipitation, temperature, potential evapotranspiration time series by elevation and for each basin and plots
#' @examples spatial_grad(DEM, temp_stations, prec_stations, ccas, grad_temp, grad_pr)
#' @export
spatial_grad <- function(DEM, temp_stations, prec_stations, ccas, grad_temp, grad_pr) {
  temp_est <- t(temp_stations)
  colnames(temp_est) <- temp_est[1, ]
  temp_est <- temp_est[-1, ]
  temp_est <- temp_est[,1:3]
  temp_est <- data.frame(temp_est)
  temp_est<- data.frame(
    lat = as.numeric(temp_est$lat),
    long = as.numeric(temp_est$long),
    Z = as.numeric(temp_est$Z)
  )
  temp_time_series <- temp_stations[-(1:3),]
  fechas <- temp_time_series$station
  fechas_convertidas <- parse_date_time(fechas,orders = "my")
  fechas_final <- as.yearmon(fechas_convertidas)
  fecha <- as.Date(as.yearmon(fechas_final))
  
  prec_est <- t(prec_stations)
  colnames(prec_est) <- prec_est[1, ]
  prec_est <- prec_est[-1, ]
  prec_est <- prec_est[,1:3]
  prec_est <- data.frame(prec_est)
  prec_est<- data.frame(
    lat = as.numeric(prec_est$lat),
    long = as.numeric(prec_est$long),
    Z = as.numeric(prec_est$Z)
  )
  prec_time_series <- prec_stations[-(1:3),]
  fechas <- prec_time_series$station
  fechas_convertidas <- parse_date_time(fechas,orders = "my")
  fechas_final <- as.yearmon(fechas_convertidas)
  fecha <- as.Date(as.yearmon(fechas_final))
  prec_time_series <- prec_time_series[,-1]
  prec_time_series <- as.data.frame(lapply(prec_time_series,as.numeric))
  prec_time_series$Date <- fecha
  
  # Precipitation interpolation
  brick_pr <- valery_interpolation(DEM = SRTM_0,
                                         sta_coor = prec_est[,c('long','lat')],
                                         sta_z = prec_est$Z,
                                         data = prec_time_series[,-ncol(prec_time_series)], # sin los tiempos 
                                         grad = grad_pr,
                                         var = 'pr')
  # Temperature interpolation
  brick_tas <- valery_interpolation(DEM = SRTM_0,
                                          sta_coor = temp_est[,c('long','lat')],
                                          sta_z = temp_est$Z,
                                          data = temp_time_series[,-ncol(temp_time_series)],
                                          grad = grad_temp,
                                          var = 'temp')
  # PET-Oudin estimation 
  brick_EP <- PE_oudin_ras_mon(brick_tas = brick_tas, dates = temp_time_series$Date)

  ############### brick pr ##############
  num_capas <- nlayers(brick_pr)
  indices <- seq(1, num_capas, by = 12)
  suma_rasterbrick <- brick()

  for (i in seq_along(indices)) {
    # Calcular el índice final para cada grupo de 12 capas
    fin <- min(indices[i] + 11, num_capas)
    
    # Sumar las capas dentro del rango actual
    suma_temp <- sum(brick_pr[[indices[i]:fin]])
    suma_12_pr_rasterbrick <- addLayer(suma_rasterbrick, suma_temp)
  }
  
  ############### brick PET ##############
  num_capas <- nlayers(brick_EP)
  indices <- seq(1, num_capas, by = 12)
  suma_rasterbrick <- brick()
  for (i in seq_along(indices)) {
    # Calcular el índice final para cada grupo de 12 capas
    fin <- min(indices[i] + 11, num_capas)
    # Sumar las capas dentro del rango actual
    suma_temp <- sum(brick_EP[[indices[i]:fin]])
    suma_12_PET_rasterbrick <- addLayer(suma_rasterbrick, suma_temp)
  }
  
  promedio_P <- mean(suma_12_pr_rasterbrick, na.rm = TRUE)
  promedio_T <- mean(brick_tas, na.rm = TRUE)
  promedio_EP <- mean(suma_12_PET_rasterbrick, na.rm = TRUE)
  
  ### CORTAR CUENCAS
  brick_P_recortado <- mask(promedio_P, ccas)
  brick_T_recortado <- mask(promedio_T, ccas)
  brick_EP_recortado <- mask(promedio_EP, ccas)
  
  # Extracting time series from basins
  ccas$NAME <- ccas$NOMBRE
  ccas$NAMES <- ccas$NOMBRE
  df_pr_ccas <- raster::extract(brick_pr,ccas,fun=mean, weights=T) %>% t #saca promedio por cuenca de pr y temp
  df_pr_ccas <- df_pr_ccas %>% as.data.frame() %>% data.frame(dates = prec_time_series$Date,.)
  names(df_pr_ccas)[-1] <- ccas$NOMBRE
  rownames(df_pr_ccas) <- NULL
  
  df_tas_ccas <- raster::extract(brick_tas,ccas,fun=mean, weights=T) %>% t
  df_tas_ccas <- df_tas_ccas %>% as.data.frame() %>% data.frame(dates = temp_time_series$Date,.)
  names(df_tas_ccas)[-1] <- ccas$NOMBRE
  rownames(df_tas_ccas) <- NULL
  
  df_EP_ccas <- raster::extract(brick_EP,ccas,fun=mean, weights=T) %>% t
  df_EP_ccas <- df_EP_ccas %>% as.data.frame() %>% data.frame(dates = temp_time_series$Date,.)
  names(df_EP_ccas)[-1] <- ccas$NOMBRE
  rownames(df_EP_ccas) <- NULL
  
  # Agregar fechas a los dataframes
  df_pr_ccas$dates <- temp_time_series$Date
  df_tas_ccas$dates <- temp_time_series$Date
  df_EP_ccas$dates <- temp_time_series$Date

  cuencas_nombres <- ccas$NOMBRE

  df_series_tiempo_final <- data.frame(dates = df_pr_ccas$dates)

  # Iteracion sobre las cuencas y agregar columnas al dataframe
  for (cuenca_nombre in cuencas_nombres) {
    pr_columna_nombre <- paste("P.", cuenca_nombre, sep = "")
    tas_columna_nombre <- paste("Tm.", cuenca_nombre, sep = "")
    EP_columna_nombre <- paste("PET.", cuenca_nombre, sep = "")
    # agregando columnas al dataframe
    df_series_tiempo_final[[pr_columna_nombre]] <- df_pr_ccas[, cuenca_nombre]
    df_series_tiempo_final[[tas_columna_nombre]] <- df_tas_ccas[, cuenca_nombre]
    df_series_tiempo_final[[EP_columna_nombre]] <- df_EP_ccas[, cuenca_nombre]
  }
  df_series_tiempo_final <- df_series_tiempo_final[, colSums(!is.na(df_series_tiempo_final)) > 0]
  
  ################################## PLOTEANDO ########################
  # List of brick to plot
  bricks <- list(brick_P_recortado, brick_T_recortado, brick_EP_recortado)
  tamanio_letra <- 1.2 #font size
  # Figure dimensions
  ancho_grafico <- 6400
  alto_grafico <- 4800
  # Nombres deseados para los archivos
  nombres_archivos <- c("graph_P_mean.png", "graph_Tm_mean.png", "graph_EP_mean.png")
  # Títulos deseados para los gráficos
  titulos_graficos <- c("Annual mean Precipitation (mm)", "Annual mean Temperature (°C)", "Annual mean Potential Evapotranspiration (mm)")
  
  # Bucle para crear y guardar los gráficos
  for (i in seq_along(bricks)) {
    # Crear un nuevo plot con dimensiones ajustadas
    png(paste(path_graphs, nombres_archivos[i], sep = ""), width = ancho_grafico, height = alto_grafico, res = 800)
    plot(bricks[[i]], main = titulos_graficos[i], cex.main = tamanio_letra * 1, cex.lab = tamanio_letra)
    # Añadir las cuencas encima del plot
    plot(ccas, add = TRUE, border = "black")
    # Cerrar el gráfico actual
    dev.off()
  }
  par(mfrow = c(1, 3))
  for (i in seq_along(bricks)) {
    plot(bricks[[i]], main = titulos_graficos[i], cex.main = tamanio_letra * 1, cex.lab = tamanio_letra)
    # Superponer el shapefile en el gráfico actual
    plot(ccas, add = TRUE, border = "black")
  }
  write.csv(df_series_tiempo_final,output)
}
