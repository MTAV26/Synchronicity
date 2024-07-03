rm(list = ls())
graphics.off()
gc()

library(ncdf4)
library(fields)
library(maptools)
library(RColorBrewer)
library(viridis)
library(RobustLinearReg)
library(wql)
library(raster)
library(sf)
library(dplyr)
library(ggplot2)
library(tidyr)
library(cowplot)
library(lubridate)
data(wrld_simpl)

# Directorios
dir_fwi = '/home/miguel/FWI_T_EU/'
dir_data = '/home/miguel/Dropbox/synchronicity/data'
dir_out = '/home/miguel/Dropbox/synchronicity/csvs/'

# Configuraciones
seasons <- c("DJF", "MAM", "JJA", "SON")
num_1seasons <- c("DJF", "MAM", "JJA", "SON")
region <- "NEU" #CEU or MED
perturbado = c("T0.0_P0.6", "T0.0_P0.8", "T0.0_P1.0", "T0.0_P1.2", "T0.0_P1.4",
               "T0.0_P1.6", "T1.0_P0.6", "T1.0_P0.8", "T1.0_P1.0", "T1.0_P1.2",
               "T1.0_P1.4", "T1.0_P1.6", "T2.0_P0.6", "T2.0_P0.8", "T2.0_P1.0",
               "T2.0_P1.2", "T2.0_P1.4", "T2.0_P1.6", "T3.0_P0.6", "T3.0_P0.8",
               "T3.0_P1.0", "T3.0_P1.2", "T3.0_P1.4", "T3.0_P1.6", "T4.0_P0.6",
               "T4.0_P0.8", "T4.0_P1.0", "T4.0_P1.2", "T4.0_P1.4", "T4.0_P1.6",
               "T5.0_P0.6", "T5.0_P0.8", "T5.0_P1.0", "T5.0_P1.2", "T5.0_P1.4",
               "T5.0_P1.6")

num_datasets <- length(perturbado)
num_seasons <- length(seasons)

# Cargar datos
fname <- file.path(dir_data, 'ERA5_1981-2022-025-daily.nc')
obs.nc <- nc_open(fname)
obs.nc$dim$lon$vals -> lon
obs.nc$dim$lat$vals -> lat
lat=rev(lat)
ptn <- expand.grid(lon, lat)
FWI <- ncvar_get(obs.nc,"fwinx")
ni = dim(FWI)[1]
nj = dim(FWI)[2]

year_fwi <- 1981:2001
fechas <- seq(as.Date("1981-01-01"), as.Date("2001-12-31"), by = "day")
FWI_obs <- FWI[, nj:1, 1:length(fechas)]

# Máscara de desierto
shapefile_path <- file.path(dir_data, 'biomas_europa.shp')
r <- raster(t(FWI_obs[,,1]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +datum=WGS84"))
r <- flip(r, direction='y')
biomas <- st_read(shapefile_path)
biomas <- st_transform(biomas, crs = crs(r))
biomas_raster <- rasterize(biomas, raster(r), field = "BIOME", fun='first', background=NA)
mask_13_matrix <- as.matrix(biomas_raster)
mask_13_matrix[mask_13_matrix != 13] <- NA
mask_13_matrix[mask_13_matrix == 13] <- 1
mask_13_matrix <- t(mask_13_matrix)[, length(lat):1]

FWI_obs <- array(apply(FWI_obs, 3, function(slice) {
  slice[mask_13_matrix == 1] <- NA
  return(slice)
}), dim = dim(FWI_obs))

# Cargar área de celdas
fname <- file.path(dir_data, 'area_km2.nc')
obs.nc <- nc_open(fname)
gridarea <- ncvar_get(obs.nc, "cell_area")
gridarea <- gridarea[, dim(gridarea)[2]:1]
nc_close(obs.nc)

if (length(fechas) != dim(FWI_obs)[3]) {
  stop("La longitud de la secuencia de fechas no coincide con la dimensión temporal de la matriz FWI")
}

# Definir estaciones
get_season <- function(date) {
  month <- month(date)
  if (month %in% c(12, 1, 2)) {
    return("Winter")
  } else if (month %in% c(3, 4, 5)) {
    return("Spring")
  } else if (month %in% c(6, 7, 8)) {
    return("Summer")
  } else {
    return("Autumn")
  }
}

 # Filtrar las fechas y datos hasta el 31 de diciembre de 2001
end_date_filter <- as.Date("2001-12-31")
fechas_filtered <- fechas[fechas <= end_date_filter]
FWI_filtered <- FWI_obs[, , fechas <= end_date_filter]

# Excluir solo enero y febrero de 1981
exclude_jan_feb_1981 <- which(year(fechas_filtered) == 1981 & month(fechas_filtered) %in% c(1, 2))

# Reemplazar valores en lugar de eliminarlos
if (length(exclude_jan_feb_1981) > 0) {
  fechas_filtered[exclude_jan_feb_1981] <- NA
  FWI_filtered[,,exclude_jan_feb_1981] <- NA
}

# Reemplazar valores en diciembre de 2001
exclude_dec_2001 <- which(year(fechas_filtered) == 2001 & month(fechas_filtered) == 12)

if (length(exclude_dec_2001) > 0) {
  fechas_filtered[exclude_dec_2001] <- NA
  FWI_filtered[,,exclude_dec_2001] <- NA
}

# Asignar estaciones a las fechas filtradas
seasons_filtered <- factor(sapply(fechas_filtered, get_season), levels = c("Winter", "Spring", "Summer", "Autumn"))

# Crear un array para almacenar los percentiles 95 por estación
percentiles_95 <- array(NA, dim = c(dim(FWI_filtered)[1], dim(FWI_filtered)[2],
                                    length(levels(seasons_filtered))), dimnames = list(NULL, NULL, levels(seasons_filtered)))

# Crear una lista para almacenar las fechas utilizadas para cada estación
fechas_por_estacion <- list()

# Calcular percentiles 95 por estación y recopilar fechas
for (i in 1:dim(FWI_filtered)[1]) {
  for (j in 1:dim(FWI_filtered)[2]) {
    for (k in 1:length(levels(seasons_filtered))) {
      season_data <- FWI_filtered[i, j, seasons_filtered == levels(seasons_filtered)[k]]
      season_dates <- fechas_filtered[seasons_filtered == levels(seasons_filtered)[k]]

      # Guardar las fechas utilizadas para cada estación en la lista
      if (length(season_dates) > 0) {
        if (is.null(fechas_por_estacion[[levels(seasons_filtered)[k]]])) {
          fechas_por_estacion[[levels(seasons_filtered)[k]]] <- season_dates
        } else {
          fechas_por_estacion[[levels(seasons_filtered)[k]]] <- unique(c(fechas_por_estacion[[levels(seasons_filtered)[k]]], season_dates))
        }
      }

      # Calcular el percentil 95
      if (length(season_data) > 0) {
        percentiles_95[i, j, k] <- quantile(season_data, probs = 0.95, na.rm = TRUE)
      }
    }
  }
}

# Mostrar el listado de las fechas por estación
for (season in names(fechas_por_estacion)) {
  cat("Estación:", season, "\n")
  cat("Fechas:\n")
  print(fechas_por_estacion[[season]])
  cat("\n")
}

# Crear matriz binaria de FWI por estaciones
seasons_all <- factor(sapply(fechas, get_season), levels = c("Winter", "Spring", "Summer", "Autumn"))
FWI_binary <- array(NA, dim = c(dim(FWI_filtered)[1], dim(FWI_filtered)[2], length(fechas_filtered)))

for (i in 1:dim(FWI_filtered)[1]) {
  for (j in 1:dim(FWI_filtered)[2]) {
    for (k in 1:length(levels(seasons_filtered))) {
      season <- levels(seasons_filtered)[k]
      percentil_95_value <- percentiles_95[i, j, k]
      indices <- which(seasons_filtered == season)
      # if(length(indices) > 0 && all(indices <= dim(FWI_filtered)[3])) {
        FWI_binary[i, j, indices] <- ifelse(FWI_filtered[i, j, indices] > percentil_95_value, 1, 0)
      # }
    }
  }
}

#===============================================================================
# Función para crear máscara de las regiones
create_mask <- function(lon, lat, polygon) {
  lon_grid <- rep(lon, each=length(lat))
  lat_grid <- rep(lat, times=length(lon))
  mask <- matrix(point.in.polygon(lon_grid, lat_grid, polygon$lon, polygon$lat), nrow=length(lat), ncol=length(lon))
  return(mask)
}

sid <- array(NA, dim = dim(FWI_binary))
for (i in 1:dim(FWI_binary)[3]) {
  sid[, , i] <- FWI_binary[, , i] * gridarea
}

# Definir región
polygon_reg <- switch(region,
                      "MED" = list(lon=c(-10, -10, 40, 40), lat=c(30, 45, 45, 30)),
                      "CEU" = list(lon=c(-10, 40, 40, -10), lat=c(45, 45, 61.3, 48)),
                      "NEU" = list(lon=c(-10, 40, 40, -10), lat=c(48, 61.3, 75, 75))
)
mask_reg <- create_mask(lon, lat, polygon_reg)
mask_reg <- ifelse(mask_reg > 0, 1, NA)
mask_reg <- t(mask_reg)

FWI_obs <- array(NA, dim = dim(sid))
for (m in 1:dim(sid)[3]) {
  FWI_obs[,,m] <- sid[,,m] * mask_reg
}

sid_mask<-apply(FWI_obs, c(1,2), mean, na.rm=T)
sid_mask[!is.na(sid_mask)]=1

# Mostrar la máscara

sid2_obs <- apply(FWI_obs, 3, sum, na.rm=TRUE)
data <- data.frame(date = fechas, value = sid2_obs)
data_new2 <- data %>%
  mutate(year = as.numeric(format(date, "%Y")),
         month = format(date, "%m"))

# Matrices para resultados
data_matrix <- matrix(NA, nrow = num_datasets, ncol = num_seasons,
                      dimnames = list(perturbado, levels(seasons_filtered)))
data_matrix_pval <- matrix(NA, nrow = num_datasets, ncol = num_seasons,
                           dimnames = list(perturbado, levels(seasons_filtered)))

dif_data_matrix <- matrix(NA, nrow = num_datasets, ncol = num_seasons,
                          dimnames = list(perturbado, levels(seasons_filtered)))
dif_data_matrix_pval <- matrix(NA, nrow = num_datasets, ncol = num_seasons,
                               dimnames = list(perturbado, levels(seasons_filtered)))

fut_matrix <- matrix(NA, nrow = num_datasets, ncol = num_seasons, 
                     dimnames = list(perturbado,levels(seasons_filtered)))
porc_fut_matrix <- matrix(NA, nrow = num_datasets, ncol = num_seasons,
                          dimnames = list(perturbado, levels(seasons_filtered)))

sid2_mask<-sid_mask* gridarea
total_reg=sum(sid2_mask , na.rm=T)

# Procesar datos perturbados
#===============================================================================
for (d in 1:num_datasets) {
  perturb <- perturbado[d]

  # fname <- file.path(dir_fwi, paste0('FWI_T5.0_P0.8_EU_1981_2010_int.nc'))
  fname <- file.path(dir_fwi, paste0('FWI_', perturb, '_EU_1981_2010_int.nc'))
  obs.nc <- nc_open(fname)
  obs.nc$dim$lon$vals -> lon
  obs.nc$dim$lat$vals -> lat
  lat=rev(lat)

  FWI_TP <- ncvar_get(obs.nc,"fwinx")
  ni = dim(FWI_TP)[1]
  nj = dim(FWI_TP)[2]
  FWI_TP=FWI_TP[, nj:1, 1:length(fechas)]


  FWI_TP <- array(apply(FWI_TP, 3, function(slice) {
    slice[mask_13_matrix == 1] <- NA
    return(slice)
  }), dim = dim(FWI_TP))


  if (length(fechas) != dim(FWI_TP)[3]) {
    stop("La longitud de la secuencia de fechas no coincide con la dimensión temporal de la matriz FWI")
  }


  # Filtrar las fechas y datos hasta el 31 de diciembre de 2001
  end_date_filter <- as.Date("2001-12-31")
  fechas_filtered <- fechas[fechas <= end_date_filter]
  FWI_filtered <- FWI_TP[, , fechas <= end_date_filter]

  # Excluir solo enero y febrero de 1981
  exclude_jan_feb_1981 <- which(year(fechas_filtered) == 1981 & month(fechas_filtered) %in% c(1, 2))

  # Reemplazar valores en lugar de eliminarlos
  if (length(exclude_jan_feb_1981) > 0) {
    fechas_filtered[exclude_jan_feb_1981] <- NA
    FWI_filtered[,,exclude_jan_feb_1981] <- NA
  }

  # Reemplazar valores en diciembre de 2001
  exclude_dec_2001 <- which(year(fechas_filtered) == 2001 & month(fechas_filtered) == 12)

  if (length(exclude_dec_2001) > 0) {
    fechas_filtered[exclude_dec_2001] <- NA
    FWI_filtered[,,exclude_dec_2001] <- NA
  }

  # Crear matriz binaria de FWI por estaciones
  seasons_all <- factor(sapply(fechas, get_season), levels = c("Winter", "Spring", "Summer", "Autumn"))
  FWI_binary <- array(NA, dim = dim(FWI_filtered))
  for (i in 1:dim(FWI_filtered)[1]) {
    for (j in 1:dim(FWI_filtered)[2]) {
      for (k in 1:length(levels(seasons_all))) {
        season <- levels(seasons_all)[k]
        percentil_95_value <- percentiles_95[i, j, k]
        indices <- which(seasons_all == season)
        FWI_binary[i, j, indices] <- ifelse(FWI_filtered[i, j, indices] > percentil_95_value, 1, 0)
      }
    }
  }
 
  sid2 <- array(NA, dim = dim(FWI_binary))

  for (g in 1:dim(FWI_binary)[3]) {
    sid2[, , g] <- FWI_binary[, , g] * gridarea
  }

  sid2_reg <- array(NA, dim = dim(sid2))
  for (h in 1:dim(sid2)[3]) {
    sid2_reg[, , h] <- sid2[, , h] * mask_reg
  }

  sid2_tp <- apply(sid2_reg, 3, sum, na.rm=TRUE)
  data2 <- data.frame(date = fechas, value = sid2_tp)
  data_new1 <- data2 %>%
    mutate(year = as.numeric(format(date, "%Y")),
           month = format(date, "%m"))

  for (m in 1:length(num_1seasons)) {
    season1 <- num_1seasons[m]
      if (season1 == 'DJF') {
       
        aux1 <- data_new1 %>%
          filter((month == "12" & year >= 1981 & year <= 2000) | 
                   (month %in% c("01", "02") & year >= 1982 & year <= 2001)) %>%
          mutate(season_year = if_else(month == "12", year + 1, year))
        
        aux2 <- data_new2 %>%
          filter((month == "12" & year >= 1981 & year <= 2000) | 
                   (month %in% c("01", "02") & year >= 1982 & year <= 2001)) %>%
          mutate(season_year = if_else(month == "12", year + 1, year))
        
        data_aggr1 <- aux1 %>%
          group_by(season_year) %>%
          summarize(value = mean(value, na.rm = TRUE))
        
        data_aggr2 <- aux2 %>%
          group_by(season_year) %>%
          summarize(value = mean(value, na.rm = TRUE))
        
        diff_data_cha <- data_aggr1$value - data_aggr2$value
        average_future_fwi_tp <- mean(data_aggr1$value, na.rm = TRUE)
        average_past_fwi_obs <- mean(data_aggr2$value, na.rm = TRUE)

        porct_past<-(100*(average_past_fwi_obs)/total_reg)
        porct_fut<-(100*(average_future_fwi_tp)/total_reg)

        perc_change <- round(100 * (average_future_fwi_tp - average_past_fwi_obs) / average_past_fwi_obs)
        t_test_result <- t.test(data_aggr1$value[2:21], data_aggr2$value[2:21])
        pval_chan <- ifelse(t_test_result$p.value < 0.05, "*", "")

        dif_change <- round(average_future_fwi_tp - average_past_fwi_obs, 2)
        dif_chan = pval_chan

      } else {
        months <- switch(season1,
                         "JJA" = c("06", "07", "08"),
                         "MAM" = c("03", "04", "05"),
                         "SON" = c("09", "10", "11")
        )
        
        df_sub1 <- subset(data_new1, month %in% months)
        df_sub2 <- subset(data_new2, month %in% months)

        data_aggr1 <- aggregate(value ~ year, df_sub1, FUN = mean)
        data_aggr2 <- aggregate(value ~ year, df_sub2, FUN = mean)
        diff_data_cha <- data_aggr1$value - data_aggr2$value

        average_future_fwi_tp <- mean(data_aggr1$value, na.rm = TRUE)
        average_past_fwi_obs <- mean(data_aggr2$value, na.rm = TRUE)

        porct_past<-(100*(average_past_fwi_obs)/total_reg)
        porct_fut<-(100*(average_future_fwi_tp)/total_reg)

        perc_change <- round(100 * (average_future_fwi_tp - average_past_fwi_obs) / average_past_fwi_obs)
        t_test_result <- t.test(data_aggr1$value, data_aggr2$value)
        pval_chan <- ifelse(t_test_result$p.value < 0.05, "*", "")

        dif_change <- round(average_future_fwi_tp - average_past_fwi_obs, 2)
        dif_chan = pval_chan
      }

    print(paste0(season1, ": porcent_change:  ", perturb, ": ", perc_change, "% ", pval_chan))
    print(paste0("diferencia_change:  ", perturb, ": ", dif_change, dif_chan))

    data_matrix[d, m] <- perc_change
    data_matrix_pval[d, m] <- pval_chan
    data_current_season <- data_matrix[, m]
    data_pval_current_season <- data_matrix_pval[, m]

    data_df <- data.frame(perturbado, Change = data_current_season, pval = data_pval_current_season)
    data_df <- data_df %>% separate(perturbado, into = c("Temperature", "Precipitation"), sep = "_")
    data_df$Temperature <- factor(data_df$Temperature, levels = unique(gsub("_.*", "", perturbado)))
    data_df$Precipitation <- factor(data_df$Precipitation, levels = unique(gsub(".*_", "", perturbado)))
    output_file <- file.path(dir_out, paste0(season1, "_percentile95_", region, "_resultados_perc_change.csv"))
    write.csv(data_df, file = output_file, row.names = TRUE)

    dif_data_matrix[d, m] <- dif_change
    dif_data_current_season <- dif_data_matrix[, m]

    dif_data_matrix_pval[d, m] <- dif_chan
    dif_data_pval_current_season <- dif_data_matrix_pval[, m]

    fut_matrix[d, m] <- average_future_fwi_tp
    dif_fut_matrix <- fut_matrix[, m]

    porc_fut_matrix[d, m] <- porct_fut
    dif_porc_fut_matrix <- porc_fut_matrix[, m]

    dif_data_df <- data.frame(perturbado, Change = dif_data_current_season,
                              pval = dif_data_pval_current_season,
                              past = average_past_fwi_obs,
                              future = dif_fut_matrix,
                              total_region = total_reg,
                              porcent_past = porct_past,
                              porcent_fut = dif_porc_fut_matrix)
    dif_data_df <- dif_data_df %>% separate(perturbado, into = c("Temperature", "Precipitation"), sep = "_")
    dif_data_df$Temperature <- factor(dif_data_df$Temperature, levels = unique(gsub("_.*", "", perturbado)))
    dif_data_df$Precipitation <- factor(dif_data_df$Precipitation, levels = unique(gsub(".*_", "", perturbado)))
    dif_output_file <- file.path(dir_out, paste0("dif_", season1, "_percentile95_", region, "_resultados_dif_change.csv"))
    write.csv(dif_data_df, file = dif_output_file, row.names = TRUE)

  }
}
