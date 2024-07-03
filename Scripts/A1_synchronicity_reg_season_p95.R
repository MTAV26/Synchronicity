
# Limpiar el entorno de trabajo
rm(list = ls())
graphics.off()
gc()

# Cargar librerías necesarias
library(ncdf4)
library(fields)
library(viridis)
library(raster)
library(dplyr)
library(ggplot2)
library(lubridate)
library(sf)

# Rutas de datos
dir_data <- '/home/miguel/Dropbox/synchronicity/data/'
dir_out <- '/home/miguel/Dropbox/synchronicity/Figures/'

#======================= Cargar datos ==========================================
fname <- file.path(dir_data, 'ERA5_1981-2022-025-daily.nc')
obs.nc <- nc_open(fname)
lon <- obs.nc$dim$lon$vals
lat <- rev(obs.nc$dim$lat$vals)
FWI <- ncvar_get(obs.nc, "fwinx")
FWI <- FWI[, dim(FWI)[2]:1, ]
nc_close(obs.nc)

################################################################################
# MASCARA DESIERTO
################################################################################
shapefile_path <- file.path(dir_data, 'biomas_europa.shp')
FWI_R <- FWI[, ncol(FWI):1,]

r <- raster(t(FWI_R[,,1]), xmn=min(lon), xmx=max(lon), ymn=min(lat), ymx=max(lat), crs=CRS("+proj=longlat +datum=WGS84"))
r <- flip(r, direction='y')# Invertir la orientación de latitudes
biomas <- st_read(shapefile_path)# Leer la capa de biomas
biomas <- st_transform(biomas, crs = crs(r))
biomas_raster <- raster(r)
biomas_raster <- rasterize(biomas, biomas_raster, field = "BIOME", fun='first', background=NA)
mask_13_matrix <- as.matrix(biomas_raster)
mask_13_matrix[mask_13_matrix != 13] <- NA
mask_13_matrix[mask_13_matrix == 13] <- 1
mask_13_matrix <- t(mask_13_matrix)
mask_13_matrix <- mask_13_matrix[, length(lat):1] 
for (i in 1:dim(FWI)[3]) {
  FWI[,,i][mask_13_matrix == 1] <- NA
}

# Cargar el archivo de áreas de celdas
fname <- file.path(dir_data, 'area_km2.nc')
obs.nc <- nc_open(fname)
gridarea <- ncvar_get(obs.nc, "cell_area")
gridarea <- gridarea[, dim(gridarea)[2]:1]
nc_close(obs.nc)


fechas <- seq(as.Date("1981-01-01"), as.Date("2022-12-31"), by = "day")
if (length(fechas) != dim(FWI)[3]) {
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
FWI_filtered <- FWI[, , fechas <= end_date_filter]

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


fechas_por_estacion <- list()
for (i in 1:dim(FWI_filtered)[1]) {
  for (j in 1:dim(FWI_filtered)[2]) {
    for (k in 1:length(levels(seasons_filtered))) {
      season_data <- FWI_filtered[i, j, seasons_filtered == levels(seasons_filtered)[k]]
      season_dates <- fechas_filtered[seasons_filtered == levels(seasons_filtered)[k]]
      if (length(season_dates) > 0) {
        if (is.null(fechas_por_estacion[[levels(seasons_filtered)[k]]])) {
          fechas_por_estacion[[levels(seasons_filtered)[k]]] <- season_dates
        } else {
          fechas_por_estacion[[levels(seasons_filtered)[k]]] <- unique(c(fechas_por_estacion[[levels(seasons_filtered)[k]]], season_dates))
        }
      }
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
FWI_binary <- array(NA, dim = dim(FWI))
for (i in 1:dim(FWI)[1]) {
  for (j in 1:dim(FWI)[2]) {
    for (k in 1:length(levels(seasons_all))) {
      season <- levels(seasons_all)[k]
      percentil_95_value <- percentiles_95[i, j, k]
      indices <- which(seasons_all == season)
      FWI_binary[i, j, indices] <- ifelse(FWI[i, j, indices] > percentil_95_value, 1, 0)
    }
  }
}

# Multiplicar los valores binarios por el área de la celda
sid4 <- array(NA, dim = dim(FWI_binary))
for (i in 1:dim(FWI_binary)[3]) {
  sid4[, , i] <- FWI_binary[, , i] * gridarea
}

# Función para crear máscaras
create_mask <- function(lon, lat, polygon) {
  lon_grid <- rep(lon, each = length(lat))
  lat_grid <- rep(lat, times = length(lon))
  mask <- matrix(point.in.polygon(lon_grid, lat_grid, polygon$lon, polygon$lat), nrow = length(lat), ncol = length(lon))
  return(mask)
}

# Definir regiones y estaciones
regions <- c(
  "MED",
  "CEU",
  "NEU"
)
seasons <- c(
  "DJF" ,
  "MAM",
  "JJA",
  "SON"
)
# library(dplyr)
for (region in regions) {
  polygon_reg <- switch(region,
                        "MED" = list(lon = c(-10, -10, 40, 40), lat = c(30, 45, 45, 30)),
                        "CEU" = list(lon = c(-10, 40, 40, -10), lat = c(45, 45, 61.3, 48)),
                        "NEU" = list(lon = c(-10, 40, 40, -10), lat = c(48, 61.3, 75, 75))
  )
  
  if( region == "NEU"){
    ip1="a"
    ip2="b"
    ip3="c"
    ip4="d"
  } else if (region == "CEU"){
    ip1="e"
    ip2="f"
    ip3="g"
    ip4="h"
  } else if (region == "MED"){
    ip1="i"
    ip2="j"
    ip3="k"
    ip4="l"
  }

  # Crear máscaras
  mask_reg <- create_mask(lon, lat, polygon_reg)
  mask_reg <- ifelse(mask_reg > 0, 1, NA)
  mask_reg <- t(mask_reg)
  
  sid=sid4
  for (i in 1:dim(sid)[3]) {
    sid[, , i] <- sid[, , i] * mask_reg
  }
  
  sid_mask<-apply(sid, c(1,2), mean, na.rm=T)
  sid_mask[!is.na(sid_mask)]=1
  sid2_mask<-sid_mask* gridarea
  total_reg=sum(sid2_mask , na.rm=T)

  sid2 <- apply(sid, 3, sum, na.rm = TRUE)
  # Agregar y preparar datos
  data <- data.frame(date = fechas, value = sid2)
  data_new1 <- data %>%
    mutate(year = as.numeric(format(date, "%Y")),
           month = format(date, "%m"))
  
  # Eliminar diciembre de 2022
  data_new1 <- data_new1 %>%
    filter(!(year == 2022 & month == "12"))
  
  for (season in seasons) {
    title <- switch(season,
                    'DJF' = paste0(ip1,") Change of Synchronicity (DJF; ", region, "): "),
                    'MAM' = paste0(ip2,") Change of Synchronicity (MAM; ", region, "): "),
                    'JJA' = paste0(ip3,") Change of Synchronicity (JJA; ", region, "): "),
                    'SON' = paste0(ip4,") Change of Synchronicity (SON; ", region, "): "),
    )
    
    df_sub <- switch(season,
                     'MAM' = filter(data_new1, month %in% c("03", "04", "05")),
                     'JJA' = filter(data_new1, month %in% c("06", "07", "08")),
                     'SON' = filter(data_new1, month %in% c("09", "10", "11")),
                     'DJF' = {
                       # Filtrar diciembre de cada año y enero, febrero del siguiente año
                       data_new1 %>%
                         filter((month == "12" & year >= 1981 & year <= 2021) | 
                                  (month %in% c("01", "02") & year >= 1982 & year <= 2022)) %>%
                         mutate(season_year = if_else(month == "12", year + 1, year))
                     })
    
    if(season =="DJF"){
      data_aggr1 <- aggregate(value ~ season_year, df_sub, FUN = mean)
    } else {
      data_aggr1 <- aggregate(value ~ year, df_sub, FUN = mean)
    }
    
    average_1981_2001 <- mean(data_aggr1$value[if (season == 'DJF') 1:20 else 1:21])
    average_2002_2022 <- mean(data_aggr1$value[if (season == 'DJF') 21:41 else 22:42])
    
    #===========================================================================
    #para porcentaje de superficie
    porct_past<-(100*(average_1981_2001)/total_reg)
    porct_fut<-(100*(average_2002_2022)/total_reg)
    porct_change<-porct_fut-porct_past
    print(paste0(" p95.", region," ", season," = past: ", round(porct_past, 1),
                 "; present: ", round(porct_fut,1), " Diff: ", round(porct_change,1), sep=""))
    #===========================================================================
    
    sd_1981_2001_mas <- average_1981_2001 + sd(data_aggr1$value[if (season == 'DJF') 1:20 else 1:21])
    sd_2002_2022_mas <- average_2002_2022 + sd(data_aggr1$value[if (season == 'DJF') 21:41 else 22:42])
    sd_1981_2001_menos <- average_1981_2001 - sd(data_aggr1$value[if (season == 'DJF') 1:20 else 1:21])
    sd_2002_2022_menos <- average_2002_2022 - sd(data_aggr1$value[if (season == 'DJF') 21:41 else 22:42])
    
    perc_change <- round(100 * (average_2002_2022 - average_1981_2001) / average_1981_2001)
    t_test_result <- t.test(data_aggr1$value[if (season == 'DJF') 1:20 else 1:21], data_aggr1$value[if (season == 'DJF') 21:41 else 22:42])
    pval_chan <- ifelse(t_test_result$p.value < 0.05, "*", "")
    
    prim <- average_1981_2001
    seg <- average_2002_2022
    FWI_reg <- data_aggr1$value
    
    if(season =="DJF"){
      year_fwi <- as.numeric(data_aggr1$season_year)
    } else {
      year_fwi <-as.numeric(data_aggr1$year)
    }
    

    # Create a data frame for ggplot
    plot_data <- data.frame(year = year_fwi, FWI_reg = FWI_reg)
    period_data1 <- data.frame(x = c(min(year_fwi), year_fwi[if (season == 'DJF') 20 else 21]), y = c(prim, prim))
    period_data2 <- data.frame(x = c(year_fwi[if (season == 'DJF') 20 else 21], max(year_fwi)), y = c(seg, seg))
    
    period_data1_sd_mas <- data.frame(x = c(min(year_fwi), year_fwi[if (season == 'DJF') 20 else 21]), y = c(sd_1981_2001_mas, sd_1981_2001_mas))
    period_data2_sd_mas <- data.frame(x = c(year_fwi[if (season == 'DJF') 20 else 21], max(year_fwi)), y = c(sd_2002_2022_mas, sd_2002_2022_mas))
    
    period_data1_sd_menos <- data.frame(x = c(min(year_fwi), year_fwi[if (season == 'DJF') 20 else 21]), y = c(sd_1981_2001_menos, sd_1981_2001_menos))
    period_data2_sd_menos <- data.frame(x = c(year_fwi[if (season == 'DJF') 20 else 21], max(year_fwi)), y = c(sd_2002_2022_menos, sd_2002_2022_menos))
    
    average_1981_2001_label <- paste0(if (season == "DJF") "1982-2001" else "1981-2001", " Average (", round(average_1981_2001, 0), " km²)")
    average_2002_2022_label <- paste0("2002-2022 Average (", round(average_2002_2022, 0), " km²)")
    porcentaje_cambio <- paste0(perc_change, "%", pval_chan)
    
    
    lim1 <- "percentil_95"
    
    file_name <- paste0("Figure2/sincro_", season, "_", region, "_", lim1, ".pdf")
    file_path <- file.path(dir_out, file_name)
    
    round_to_nearest <- function(x, base) {
      base * round(x / base)
    }

    max_value <- max(plot_data$FWI_reg)
    min_value <- min(plot_data$FWI_reg)
    

    ylim_lower <- min_value - 100
    ylim_upper <- max_value + 100
    tamaño_titulo <- 2

    # Crear el gráfico y guardarlo en un archivo PDF
    pdf(file = file_path, width = 8.5, height = 6)
    plot(plot_data$year, plot_data$FWI_reg, col = "black", pch = 19, cex = 1.5,
         xlab = "Years", ylab = paste0("Surface [km²] with Synchronicity for FWI > p95",sep = ""),
         main = paste0(title, porcentaje_cambio), ylim = c(ylim_lower, ylim_upper), axes = FALSE,
         cex.main = tamaño_titulo)
    grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
    abline(v = year_fwi[if (season == 'DJF') 20 else 21], col = "black", lty = 2, lwd = 2)
    lines(period_data1$x, period_data1$y, col = "blue4", lwd = 2)
    lines(period_data2$x, period_data2$y, col = "darkred", lwd = 2)
    lines(period_data1_sd_mas$x, period_data1_sd_mas$y, col = "blue4", lwd = 2, lty = 2)
    lines(period_data1_sd_menos$x, period_data1_sd_menos$y, col = "blue4", lwd = 2, lty = 2)
    lines(period_data2_sd_mas$x, period_data2_sd_mas$y, col = "darkred", lwd = 2, lty = 2)
    lines(period_data2_sd_menos$x, period_data2_sd_menos$y, col = "darkred", lwd = 2, lty = 2)
    axis(1) # Eje x
    if(region == "NEU"){
      axis(2, at = seq(ylim_lower, ylim_upper, by = 50000),
           labels = format(seq(ylim_lower, ylim_upper, by = 50000), big.mark = ",", scientific = FALSE))
    } else{
      axis(2, at = seq(ylim_lower, ylim_upper, by = 50000),
           labels = format(seq(ylim_lower, ylim_upper, by = 50000), big.mark = ",", scientific = FALSE))
    }
    box()
    legend("topleft", legend = c("Observed Data", average_1981_2001_label, average_2002_2022_label),
           col = c("black", "blue4", "darkred"), pch = c(19, NA, NA), lty = c(NA, 1, 1), lwd = 1, bty = "n", cex = 1.15)
    dev.off()
  }
}
