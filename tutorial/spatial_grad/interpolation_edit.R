valery_interpolation <- function(DEM ,sta_coor, sta_z, data , 
                                 p = 2 , Θ, var = 'pr'){
  DEMnl <- brick(brick(DEM),nl=nrow(sta_coor))
  DEMnl[] <- DEM[]
  nDates <- nrow(data)
  
  w <- brick(DEM,nl=nrow(sta_coor))
  for (i in 1:nrow(sta_coor)){
    w[[i]] <- 1/distanceFromPoints(DEM,sta_coor[i,])^p
  }
  Sw <- sum(w[[1:nrow(sta_coor)]])
  
  raster <- brick(DEM)
  if (var=='pr'){
    for (i in 1:nrow(data)){
      Prcorr <- (as.numeric(data[i,])*exp(Θ*(DEMnl-sta_z)))
      raster[[i]] <- sum(w* Prcorr)/Sw
      cat('avanzando ',i/nDates*100,'% \n')
    }
  } else if (var=='temp'){
    for (i in 1:nrow(data)){
      Tcorr <- (Θ*(DEMnl-sta_z)+as.numeric(data[i,]))
      raster[[i]] <- sum(w* Tcorr)/Sw
      cat('avanzando ',i/nDates*100,'% \n')
    }    
  }
  return(raster)
}

valery_interpolation_NA_v2 <- function(DEM ,sta_coor, sta_z, data , 
                                       p = 2 , Θ, var = 'pr'){
  # DEMnl <- brick(brick(DEM),nl=nrow(sta_coor))
  # DEMnl[] <- DEM[]
  nDates <- nrow(data)
  # ras_z <- DEM[]
  
  
  # w <- brick(DEM,nl=nrow(sta_coor))
  w <- t( sapply(1:nrow(sta_coor),function(i){
    (1/distanceFromPoints(DEM,sta_coor[i,])^p)[]
  })  )
  # for (i in 1:nrow(sta_coor)){
  #   w[[i]] <- 1/distanceFromPoints(DEM,sta_coor[i,])^p
  # }
  fac <- matrix(rep(DEM[],length(sta_z)),nrow=length(sta_z),byrow = T) - sta_z
  if (var=='pr') fac <- exp(Θ*fac)
  else if(var=='temp')   fac <- (Θ*fac)
  
  raster <- brick(brick(DEM),nl = nrow(data),values=T)
  # brick(DEM,nl = nrow(data),values=T)
  if (var=='pr'){
    for (i in 1:nrow(data)){
      ind <- which(!is.na(data[i,]))
      Prcorr <- as.numeric(data[i,ind])*fac[ind,]
      # Prcorr <- (as.numeric(data[i,ind])*exp(Θ*(DEMnl[[ind]]-sta_z[ind])))
      # Sw <- sum(w[[ind]])
      Sw <- colSums(w[ind,])
      # raster[[i]] <- sum(w[[ind]]* Prcorr)/Sw
      raster[[i]] <- colSums(w[ind,]* Prcorr)/Sw
      
      cat('avanzando ',round(i/nDates*100,2),'%    \r')   
      # round(i/nDates*100,2)
    }
    
  } else if (var=='temp'){
    for (i in 1:nrow(data)){
      ind <- which(!is.na(data[i,]))
      # Tcorr <- (Θ*(DEMnl[[ind]]-sta_z[ind])+as.numeric(data[i,ind]))
      Tcorr <- (fac[ind,]+as.numeric(data[i,ind]))
      # Sw <- sum(w[[ind]])
      Sw <- colSums(w[ind,])
      # raster[[i]] <- sum(w[[ind]]* Tcorr)/Sw
      raster[[i]] <- colSums(w[ind,]* Tcorr)/Sw
      cat('avanzando ',round(i/nDates*100,2),'%   ',' \r') #,'i=',i
      # round(i/nDates*100,2)
      # format(round(i/nDates*100, 2), nsmall = 2)
    }    
  }
  return(raster)
}


PE_oudin_ras_mon_v2 <- function(brick_tas,dates,cca = NULL,n_cores = 4){
  require(dplyr)
  require(airGR)
  require(purrr)
  #require(parallel)
  rasterEP <- brick_tas #brick(brick_tas)
  # names(rasterEP) <- names(rasterT)
  
  ### ---- OBTENER FECHAS DIARIAS  ---- ###
  dates_i <- (dates %>% as.POSIXlt())
  dates_i$mday <- rep(1,length(dates_i))         # puede obviarse si todos son 1 (dia_mes)
  dates_f <- dates_i
  dates_f$mon <- dates_f$mon +1 
  dates_f$mday <- dates_f$mday -1
  # dates_i <- as.POSIXct(dates_i)
  # dates_f <- as.POSIXct(dates_f)
  
  dates_d <- lapply(1:length(dates_i),
                    function(k) seq(dates_i[k], dates_f[k], by = 'days')) 
  n_dias <- sapply(dates_d, length)
  dates_d <- reduce(dates_d, c)
  
  lats <- coordinates(brick_tas)[,'y']
  dias_julianos <- as.numeric(format(dates_d,'%j'))
  
data <- lapply(1:ncell(rasterEP), function(i) {
    cat(round(i/ncell(rasterEP)),'%','\r')
    temp <- as.vector(brick_tas[i])
    temp <- rep(temp, n_dias)
    # lat <- coordinates(brick_tas)[,'y'][i]
    lat <- lats[i]
    # rasterEP[i] <- 
    PE_Oudin(dias_julianos,
             temp, lat,'deg') %>%    #days$yday+1  #days11
      data.frame(dates = rep(dates_i, n_dias), .) %>% group_by(dates) %>% 
      summarize_all(sum) %>% .$. #%>% matrix(nrow = 1)
}) %>% 
    simplify2array(., higher = (TRUE == "array"))
  
  for(i in 1:nlayers(rasterEP)){
    rasterEP[[i]][] <- data[i,]
  }
  
  if(!is.null(cca)){
    values_EP <- raster::extract(rasterEP,cca,weights=TRUE,mean) %>% as.numeric
    # df_EP <- 
    data.frame(dates = dates,EP = values_EP) %>% return()
  } else { return(rasterEP)}
}

                    
spatial_grad <- function(DEM, temp_stations, prec_stations, ccas,gradiente_temp,gradiente_pp) {
  # preprocesamiento - datos 
  ## TEMPERATURA 
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
  fecha <- temp_time_series$station
  fecha <- as.Date(paste0("01-", fecha), format = "%d-%b-%Y")
  temp_time_series <- temp_time_series[,-1]
  #temp_time_series <- na.omit(temp_time_series)
  temp_time_series <- as.data.frame(lapply(temp_time_series,as.numeric))

  temp_time_series$Date <- fecha
  
  
  ## PRECIPITACION 
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
  fecha <- prec_time_series$station
  fecha <- as.Date(paste0("01-", fecha), format = "%d-%b-%Y")
  prec_time_series <- prec_time_series[,-1]
  #prec_time_series <- na.omit(prec_time_series)
  prec_time_series <- as.data.frame(lapply(prec_time_series,as.numeric))

  prec_time_series$Date <- fecha
  
  
  

  
  # Interpolación de la precipitación
  brick_pr <- valery_interpolation_NA_v2(DEM = SRTM_0,
                                         sta_coor = prec_est[,c('long','lat')],
                                         sta_z = prec_est$Z,
                                         data = prec_time_series[,-ncol(prec_time_series)], # sin los tiempos 
                                         p = 2,
                                         Θ = 4*10^(-4),
                                         var = 'pr')
  
  
  # Interpolación de la temperatura
  brick_tas <- valery_interpolation_NA_v2(DEM = SRTM_0,
                                          sta_coor = temp_est[,c('long','lat')],
                                          sta_z = temp_est$Z,
                                          data = temp_time_series[,-ncol(temp_time_series)],
                                          p = 2,
                                          Θ = -6.5/1000,
                                          var = 'temp')
  
  # Cálculo de la evapotranspiración potencial (EP)
  brick_EP <- PE_oudin_ras_mon_v2(brick_tas = brick_tas, dates = temp_time_series$Date)
  
  #exportando archivos raster

  brick_tas_cropped <- mask(crop(brick_tas, ccas), ccas)
  brick_pr_cropped <- mask(crop(brick_pr, ccas), ccas)
  brick_EP_cropped <- mask(crop(brick_EP, ccas), ccas)
  
  
  # Promedio ANUAL para temperatura y evapotranspiración (variables promedio)
  mean_tas <- mean(brick_tas_cropped)
  
  # Para precipitación, calculamos el promedio mensual primero y luego la suma para acumulado anual
  fechas <- temp_time_series$Date
  
  promedio_mensual <- function(brick, fechas) {
    months <- format(fechas, "%m")
    stack_mensual <- stack()
    for (m in 1:12) {
      capas_mes <- which(months == sprintf("%02d", m))
      if (length(capas_mes) > 0) {
        promedio_mes <- mean(brick[[capas_mes]])
        stack_mensual <- addLayer(stack_mensual, promedio_mes)
      }
    }
    return(stack_mensual)
  }
  
  tas_mensual <- promedio_mensual(brick_tas_cropped, fechas)
  pr_mensual  <- promedio_mensual(brick_pr_cropped, fechas)
  etp_mensual <- promedio_mensual(brick_EP_cropped, fechas)
  
  # Nombres de los meses
  month_names <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
  names(tas_mensual) <- month_names
  names(pr_mensual)  <- month_names
  names(etp_mensual) <- month_names
  
  # Precipitación anual acumulada como suma de promedios mensuales
  pr_anual_acumulada <- calc(pr_mensual, sum)
  names(pr_anual_acumulada) <- "Precipitacion_Anual_Acumulada"
  
  etp_anual_acumulada <- calc(etp_mensual, sum)
  names(etp_anual_acumulada) <- "Evapotranspiracion_Anual_Acumulada"
  
  
  # ============================
  # EXPORTAR NETCDF
  # ============================
  
  # Promedios ANUALES (tas, etp)
  writeRaster(mean_tas, file.path(path_graphs, "tas_promedio_anual.nc"),
              format = "CDF", varname = "tas", zname = "time", overwrite = TRUE)
  
  writeRaster(etp_anual_acumulada, file.path(path_graphs, "etp_promedio_anual.nc"),
              format = "CDF", varname = "etp", zname = "time", overwrite = TRUE)
  
  # Precipitación anual acumulada (suma de promedios mensuales)
  writeRaster(pr_anual_acumulada, file.path(path_graphs, "pr_promedio_anual.nc"),
              format = "CDF", varname = "pr", zname = "time", overwrite = TRUE)
  
  # Promedios MENSUALES
  writeRaster(tas_mensual, file.path(path_graphs, "tas_promedio_mensual.nc"),
              format = "CDF", varname = "tas", zname = "month", overwrite = TRUE)
  
  writeRaster(pr_mensual, file.path(path_graphs, "pr_promedio_mensual.nc"),
              format = "CDF", varname = "pr", zname = "month", overwrite = TRUE)
  
  writeRaster(etp_mensual, file.path(path_graphs, "etp_promedio_mensual.nc"),
              format = "CDF", varname = "etp", zname = "month", overwrite = TRUE)

  
  
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
  
  # Ver el resultado
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
  
  
  ##########
  
  promedio_P <- mean(suma_12_pr_rasterbrick, na.rm = TRUE)
  promedio_T <- mean(brick_tas, na.rm = TRUE)
  promedio_EP <- mean(suma_12_PET_rasterbrick, na.rm = TRUE)
  
  ### CORTAR CUENCAS
  brick_P_recortado <- mask(promedio_P, ccas)
  brick_T_recortado <- mask(promedio_T, ccas)
  brick_EP_recortado <- mask(promedio_EP, ccas)
  
  #SACANDO PARA SERIES DE TIEMPO
  # Extracción de datos por cuenca
  ccas$NAME <- ccas$NOMBRE
  ccas$NAMES <- ccas$NOMBRE
  
  df_pr_ccas <- raster::extract(brick_pr,ccas,fun=mean, weights=T) %>% t #saca promedio por cuenca de pr y temp
  df_pr_ccas <- df_pr_ccas %>% as.data.frame() %>% data.frame(dates =fecha,.)
  names(df_pr_ccas)[-1] <- ccas$NOMBRE
  rownames(df_pr_ccas) <- NULL
  
  df_tas_ccas <- raster::extract(brick_tas,ccas,fun=mean, weights=T) %>% t
  df_tas_ccas <- df_tas_ccas %>% as.data.frame() %>% data.frame(dates = fecha,.)
  names(df_tas_ccas)[-1] <- ccas$NOMBRE
  rownames(df_tas_ccas) <- NULL
  
  df_EP_ccas <- raster::extract(brick_EP,ccas,fun=mean, weights=T) %>% t
  df_EP_ccas <- df_EP_ccas %>% as.data.frame() %>% data.frame(dates = fecha,.)
  names(df_EP_ccas)[-1] <- ccas$NOMBRE
  rownames(df_EP_ccas) <- NULL
  
  # Agregar fechas a los dataframes
  df_pr_ccas$dates <- temp_time_series$Date
  df_tas_ccas$dates <- temp_time_series$Date
  
  df_EP_ccas$dates <- temp_time_series$Date

###########################
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
  
  # Lista de bricks a plotear
  bricks <- list(brick_P_recortado, brick_T_recortado, brick_EP_recortado)
  
  # Tamaño de letra deseado
  tamanio_letra <- 1.2
  
  # Dimensiones del gráfico 
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
  # Bucle 
  for (i in seq_along(bricks)) {
    plot(bricks[[i]], main = titulos_graficos[i], cex.main = tamanio_letra * 1, cex.lab = tamanio_letra)
    # Superponer el shapefile en el gráfico actual
    plot(ccas, add = TRUE, border = "black")
  }#################################################################
  library(openxlsx)
  write.csv(df_series_tiempo_final,output)

}
