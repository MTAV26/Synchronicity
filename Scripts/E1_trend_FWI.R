rm(list = ls())
graphics.off()
gc()

library(maps)
library(ncdf4)
library(fields)
library(maptools)
library(RColorBrewer)
library(viridis)
library(RobustLinearReg)
library(wql) 
library(sf)
library(dplyr)
library(raster)
library(oce)
library(terra)
library(sp)
library(rgdal)
library(rworldmap)
library(rnaturalearth)
library(ocedata)
data("coastlineWorld")
library(graticule) # coordenadas en el mapa
data(wrld_simpl)

dir_data = '/home/miguel/Dropbox/synchronicity/data/'
dir_out = '/home/miguel/Dropbox/synchronicity/Figures/'

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
ocean <- sf::st_transform(shapename , crs = crs) # 
shapefile_path <- file.path(dir_data, 'biomas_europa.shp')

seasons= c('DJF','MAM', 'JJA', 'SON')
for (ss in 1:length(seasons)) {
  season = seasons[ss]
  
      fname <- file.path(dir_data, 'ERA5_1981-2022-025-season.nc')
      obs.nc <- nc_open(fname)
      obs.nc$dim$lon$vals -> lon
      obs.nc$dim$lat$vals -> lat
      lat=rev(lat)
      FWI <- ncvar_get(obs.nc,"fwinx") 
      FWI=FWI[,length(lat):1,]
      ptn <- expand.grid(lon, lat)
      
      ni = dim(FWI)[1]
      nj = dim(FWI)[2]
      
      #============================================================================#
      # MASCARA DESIERTO
      #============================================================================#
      FWI_R <- FWI[, ncol(FWI):1,]
      r <- raster(t(FWI_R[,,1]), 
                  xmn=min(lon), xmx=max(lon), 
                  ymn=min(lat), ymx=max(lat),
                  crs=CRS("+proj=longlat +datum=WGS84"))
      
      r <- flip(r, direction='y')# Invertir la orientación de latitudes
      biomas <- st_read(shapefile_path)# Leer la capa de biomas
      biomas <- st_transform(biomas, crs = crs(r))
      biomas_raster <- raster(r)
      biomas_raster <- rasterize(biomas, biomas_raster, field = "BIOME",
                                 fun='first', background=NA)
      # Convertir el raster biomas_raster == 13 a matriz 1 y sino NA
      mask_13_matrix <- as.matrix(biomas_raster)
      mask_13_matrix[mask_13_matrix != 13] <- NA
      mask_13_matrix[mask_13_matrix == 13] <- 1
      mask_13_matrix <- t(mask_13_matrix)
      mask_13_matrix <- mask_13_matrix[, length(lat):1] 
      # Aplicar la máscara al raster FWI
      for (i in 1:dim(FWI)[3]) {
        FWI[,,i][mask_13_matrix == 1] <- NA
      }
      #======================== season ============================================# 
      
      if (season == 'JJA') {
        FWI=FWI[,,seq(3,dim(FWI)[3],4)] 
        year_fwi=1981:2022
      } else if (season == 'MAM') {
        FWI=FWI[,,seq(2,dim(FWI)[3],4)]
        year_fwi=1981:2022
      } else if (season == 'SON') {
        FWI=FWI[,,seq(4,dim(FWI)[3],4)]
        year_fwi=1981:2022
      } else if (season == 'DJF') {
        FWI=FWI[,,seq(5,dim(FWI)[3]-1,4)]
        year_fwi=1982:2022
      }
      
      #========= TREND FWI MAP =============#
      
      changes <- matrix(data = NA,nrow = ni,ncol = nj)
      pval   <- matrix(data = NA,nrow = ni,ncol = nj)
      
      for (i in 1:ni) {
        for (j in 1:nj) {
          if (!is.na(sum(FWI[i,j,]))) {
            dum = mannKen(FWI[i,j,])
            changes[i,j] = round(100*(dum$sen.slope*length(year_fwi))/mean(FWI[i,j,]))
            pval[i,j] = dum$p.value
            rm(dum)
          }
        }
      }
      
      pmaxFWI=100
      data2plot=changes
      data2plot[data2plot>pmaxFWI]=pmaxFWI
      data2plot[data2plot< -pmaxFWI]=-pmaxFWI
      
      zlim <- c(-pmaxFWI, pmaxFWI)
      my_breaks <- seq(zlim[1],zlim[2], length.out = 11)
      my_col <- rev(colorRampPalette(brewer.pal(11, "Spectral"))(10))
      def_breaks = seq(0,zlim[2],length.out=length(my_breaks))
      
      pval1=pval
      pval1[pval>0.05]=NA
      pval1[!is.na(pval1)]=1
      # image(pval1)
      indices <- which(pval1 == 1, arr.ind = TRUE)
      coords <- SpatialPoints(cbind(lon[indices[,1]],
                                    lat[indices[,2]]),
                              proj4string=CRS("+proj=longlat +datum=WGS84"))
      coords_robin <- spTransform(coords, CRS("+proj=robin"))
      
      if (season == "DJF"){
        tit1 = paste0("a) Seasonal trend (DJF) of the FWI", sep="")
      } else if (season == "MAM"){
        tit1 = paste0("b) Seasonal trend (MAM) of the FWI", sep="")
      } else if (season == "JJA"){
        tit1 = paste0("c) Seasonal trend (JJA) of the FWI", sep="")
      } else if (season == "SON"){
        tit1 = paste0("d) Seasonal trend (SON) of the FWI", sep="")
      }
      
      pdf(paste(dir_out,"trend_fwi_",season,".pdf",sep = ""),width=7.6,
          height=7.6)  
      mapPlot(coastlineWorld, col="white",
              projection="+proj=robin",
              longitudelim=c(-0, 29), # Ajusta este rango según necesites
              latitudelim=c(31.5, 69.5), # Ajusta este rango según necesites
              drawBox = TRUE, 
              main = tit1,
              cex.main = 2,
              # line = 1, adj = 0.5,
              axes = TRUE, grid = FALSE)
      
      mapImage(lon, lat, data2plot, col = my_col, breaks = my_breaks)
      points(coords_robin, pch = 20, cex = 0.2, col = "black") # Ajusta el color y el tamaño de los puntos según tus necesidades
      mapImage(lon, lat, mask_13_matrix, col = "gray")
      plot(ocean, col = "#E0FFFF", add =TRUE)
      mapGrid(dlongitude = 15, dlatitude = 15,  lty = 1,lwd = 1, col="black")
      dev.off()
      
      # if(season == "SON"){
      #   pdf(paste(dir_out,"Figure3/paleta.pdf",sep = ""),width=8, height=4)
      # 
      #   image.plot(legend.only = TRUE, zlim = zlim, col = my_col, breaks = my_breaks,
      #              horizontal = TRUE, legend.width = 1.5, legend.shrink = 0.75,
      #              legend.args = list(text = "TREND", side = 1, line = 2, cex = 1.5)
      #   )
      #   dev.off()
      # }

}
