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
  fecha <- as.Date(as.yearmon(fecha_final))
  
  temp_time_series <- temp_time_series[,-1]
  temp_time_series <- as.data.frame(lapply(temp_time_series,as.numeric))
  temp_time_series$Date <- fecha
  
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
  fecha <- as.Date(as.yearmon(fecha_final))
  
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
