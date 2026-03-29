spatial_krige <- function(data=p, shp=v)
  {x <- as.numeric(p[2,2:ncol(p)])
  y <- as.numeric(p[1,2:ncol(p)])
  z <- matrix(NA,nrow=(nrow(p)-2), ncol=(ncol(p)-1))
  geodata <- list()
  emp_variogram <- list(); FIT_VARIOGRAM <- list()
  krico <- list(); krobj <- list()
  pred.matrix <- list(); krige.raster <- list()
  r <- list(); r.basin <- list(); ts <- list()
  grid <-  expand.grid(seq(min(x),max(x),res),
                       seq(min(y),max(y),res))
  for (i in 1:(nrow(p)-2)) {
    z[i,] <- as.numeric(p[i+2,2:ncol(p)])
    geodata[[i]]  <-  as.geodata(cbind(x,y,z[i,]))
    emp_variogram[[i]]  <-  variog(geodata[[i]])
    FIT_VARIOGRAM[[i]] <-  variofit(emp_variogram[[i]])
    krico[[i]] <-  krige.control(type.krige="OK",
                                 obj.model=FIT_VARIOGRAM[[i]])
    krobj[[i]] <-  krige.conv(geodata[[i]], locations=grid, krige = krico[[i]])$predict
    pred.matrix[[i]] <- matrix(krobj[[i]], 
                               nrow = length(unique(grid$Var2)),
                               ncol = length(unique(grid$Var1)),
                               byrow =  T)
    krige.raster[[i]] <- flip(raster(pred.matrix[[i]]),"y")
    extent(krige.raster[[i]]) <- c(min(grid$Var1), max(grid$Var1), 
                                   min(grid$Var2), max(grid$Var2))
    crs(krige.raster[[i]]) <- CRS("+proj=longlat +datum=WGS84")
    r[[i]] <- rast(krige.raster[[i]])
    r.basin[[i]] <- mask(r[[i]],v)
    ts[[i]] <- extract(r.basin[[i]], v, fun=mean)
  }
  ts.krg <- unlist(lapply(ts, '[[', 2))
  par(mfrow = c(1, 2))
  plot(r.basin[[1]])
  points(x,y, pch="+")
  plot(v, add=T)
  plot(ts.krg, type="l",ylab="Regionalized")
}
