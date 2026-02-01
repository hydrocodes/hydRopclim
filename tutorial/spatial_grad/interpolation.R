valery_interpolation <- function(DEM ,sta_coor, sta_z, data, 
                                 grad, var = 'pr'){
  nDates <- nrow(data)
  w <- t( sapply(1:nrow(sta_coor),function(i){
    (1/distanceFromPoints(DEM,sta_coor[i,])^2)[]
  })  )
  fac <- matrix(rep(DEM[],length(sta_z)),nrow=length(sta_z),byrow = T) - sta_z
  if (var=='pr') fac <- exp(grad*fac)
  else if(var=='temp')   fac <- (grad*fac)
  
  raster <- brick(brick(DEM),nl = nrow(data),values=T)
  if (var=='pr'){
    for (i in 1:nrow(data)){
      ind <- which(!is.na(data[i,]))
      Prcorr <- as.numeric(data[i,ind])*fac[ind,]
      Sw <- colSums(w[ind,])
      raster[[i]] <- colSums(w[ind,]* Prcorr)/Sw
      
      cat('loading...',round(i/nDates*100,2),'%    \r')   
    }
  } else if (var=='temp'){
    for (i in 1:nrow(data)){
      ind <- which(!is.na(data[i,]))
      Tcorr <- (fac[ind,]+as.numeric(data[i,ind]))
      Sw <- colSums(w[ind,]) 
      raster[[i]] <- colSums(w[ind,]* Tcorr)/Sw
      cat('loading...',round(i/nDates*100,2),'%   ',' \r')
    }    
  }
  return(raster)
}

PE_oudin_ras_mon <- function(brick_tas,dates,cca=NULL,n_cores = 4){
  rasterEP <- brick_tas
  ### Daily dates
  dates_i <- (dates %>% as.POSIXlt())
  dates_i$mday <- rep(1,length(dates_i)) # puede obviarse si todos son 1 (dia_mes)
  dates_f <- dates_i
  dates_f$mon <- dates_f$mon +1 
  dates_f$mday <- dates_f$mday -1
  dates_d <- lapply(1:length(dates_i),
                    function(k) seq(dates_i[k],dates_f[k],by = 'days')) 
  n_dias <- sapply(dates_d,length)
  dates_d <- reduce(dates_d, c)
  
  lats <- coordinates(brick_tas)[,'y']
  dias_julianos <- as.numeric(format(dates_d,'%j'))
  
  data <- lapply(1:ncell(rasterEP),function(i) {
    cat(round(i/ncell(rasterEP),5),'   %','\r')
    temp <- as.vector(brick_tas[i])
    temp <- rep(temp,n_dias)
    lat <- lats[i]
    PE_Oudin(dias_julianos,
             temp, lat,'deg') %>%    #days$yday+1  #days11
      data.frame(dates = rep(dates_i,n_dias),.) %>% group_by(dates) %>% 
      summarize_all(sum) %>% .$.
  }) %>%
    simplify2array(., higher = (TRUE == "array"))
  for(i in 1:nlayers(rasterEP)){
    rasterEP[[i]][] <- data[i,]
  }
  if(!is.null(cca)){
    values_EP <- raster::extract(rasterEP,cca,weights=TRUE,mean) %>% as.numeric
    data.frame(dates = dates,EP = values_EP) %>% return()
  } else { return(rasterEP)}
}

