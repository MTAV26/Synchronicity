
# Limpiar entorno y cargar librerías
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
library(ncdf4)
library(sp)
library(maptools) # loads sp library too
library(RColorBrewer) # creates nice color schemes
library(classInt) # finds class intervals for continuous variables
library(fields)

library(maps)
library(pracma)
library(verification)
library(psych)
library(sf)
library(oce)
library(raster)
library(ggplot2)
library(terra)
library(graticule)
library(rgdal)
library(rworldmap)
library(rnaturalearth)
data(wrld_simpl)

# Directorios
dir_fwi = '/home/miguel/FWI_T_EU/'
dir_data = '/home/miguel/Dropbox/synchronicity/data'
dir_out = '/home/miguel/Dropbox/synchronicity/Figures/'

# Configuraciones
seasons <- c("DJF", "MAM", "JJA", "SON")
num_1seasons <- c("DJF", "MAM", "JJA", "SON")
region <- "MED"


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
rm(FWI)
gc()
image(FWI_obs[,,9])

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

image(FWI_obs[,,9])
FWI_obs <- array(apply(FWI_obs, 3, function(slice) {
  slice[mask_13_matrix == 1] <- NA
  return(slice)
}), dim = dim(FWI_obs))

image(FWI_obs[,,40])
rm(biomas, biomas_raster, r, ptn)
gc()

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
seasons_filtered <- factor(sapply(fechas_filtered, get_season),
                           levels = c("Winter", "Spring", "Summer", "Autumn"))

# Crear un array para almacenar los percentiles 95 por estación
percentiles_95 <- array(NA, dim = c(dim(FWI_filtered)[1], dim(FWI_filtered)[2],
                                    length(levels(seasons_filtered))), 
                        dimnames = list(NULL, NULL, levels(seasons_filtered)))

fechas_por_estacion <- list()
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


#===============================================================================
#PLOT
data("coastlineWorld")
crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m"
worldmap <- rworldmap::getMap(resolution = "coarse") # Load countries 
raster::crs(worldmap)
worldmap <- sp::spTransform(worldmap,crs) # Project to Robinson

lati <- c(-90, -45, 0, 45, 90)
long <- c(-180, -135,-90, -45, 0, 45, 90, 135, 180)
labs <- graticule::graticule_labels(lons = long, lats = lati, xline = -180,
                                    yline = 90, proj = crs) # labels for the graticules
lines <- graticule::graticule(lons = long, lats = lati, proj = crs) # graticules
shapename <- read_sf('/home/miguel/Dropbox/capas/ne_10m_ocean.shp')
ocean <- sf::st_transform(shapename , crs = crs) 

library(fields)
library(viridis)
library(maps)
library(oce)

# Crear la paleta inferno
paleta_inferno <- viridis(256, option = "inferno")
zlim_comun <- c(0.05, 130.18)
my_breaks <- seq(zlim_comun[1], zlim_comun[2], length.out = 257)
my_col <- paleta_inferno

# Dibujar los mapas

pdf(paste(dir_out,"/FWI_95_",season,".pdf",sep = ""),width=8, height=7.5)
par(mfrow = c(2, 2), mar = c(4, 4, 2, 2), oma = c(0, 0, 3, 0))  # Ajustar márgenes y espacio

for (isc in 1:dim(percentiles_95)[3]) {
  
  if (isc == 1){
    tit1 = "a) Percentile 95 (DJF) of the FWI"
  } else if (isc == 2){
    tit1 = "b) Percentile 95 (MAM) of the FWI"
  } else if (isc == 3){
    tit1 = "c) Percentile 95 (JJA) of the FWI"
  } else if (isc == 4){
    tit1 = "d) Percentile 95 (SON) of the FWI"
  }
  
  maxFWI=130.18
  data2plot=percentiles_95[,,isc]

  zlim <- c(0, maxFWI)
  my_breaks <- seq(zlim[1],zlim[2], length.out = 21)
  my_col <-rev(inferno(length(my_breaks)-1))
  def_breaks = seq(0,zlim[2],length.out=length(my_breaks))
  
  # Plot del mapa
  mapPlot(coastlineWorld, col="white",
          projection="+proj=robin",
          longitudelim=c(-0, 29), # Ajusta este rango según necesites
          latitudelim=c(31.5, 69.5), # Ajusta este rango según necesites
          drawBox = TRUE, main = tit1, cex.main = 1,
          # line = 1, adj = 0.5,
          axes = TRUE, grid = FALSE)
  
  mapImage(lon, lat, data2plot, col = my_col, breaks = my_breaks)
  mapImage(lon, lat, mask_13_matrix, col = "gray")
  plot(ocean, col = "#E0FFFF", add =TRUE)
  mapGrid(dlongitude = 15, dlatitude = 15,  lty = 1,lwd = 1, col="black")
  
}

dev.off()
# 
# pdf(paste(dir_out,"/paleta_p95.pdf",sep = ""),width=8, height=4)
# image.plot(
#   # zlim=c2,
#   zlim=c(0, 1),
#   legend.only=TRUE,
#   smallplot=c(.1,.9, .1,.2),
#   col=my_col,
#   breaks = my_breaks,
#   horizontal = TRUE)
# 
# dev.off()

