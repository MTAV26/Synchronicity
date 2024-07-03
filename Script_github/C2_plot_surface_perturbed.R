rm(list = ls())
graphics.off()
gc()

library(ggplot2)
library(scales)
library(cowplot)

breaks<-c(-15, -10, -5, -2.5, -1, 0, 1, 2.5, 5, 10, 15)
num_intervals <- length(breaks)+1  # No sumar nada porque "NA" se maneja en otro lugar
colors <- c('#3e2401','#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3',
            '#c7eae5','#80cdc1','#35978f','#01665e','#003c30', '#001f14')
custom_palette <- rev(colors)

# Función para truncar los valores de Change a un máximo de 300 y manejar NA
truncar_cambio <- function(data_df) {
  data_df$dif_pocet_past_fut <- as.numeric(data_df$dif_pocet_past_fut)
  data_df$dif_pocet_past_fut[is.infinite(data_df$dif_pocet_past_fut)] <- NA
  data_df$dif_pocet_past_fut[data_df$dif_pocet_past_fut == 0] <- NA
  return(data_df)
}

crear_grafico <- function(data_df, title) {
  data_df <- truncar_cambio(data_df)
  labels <- c("< -15%", "-15 to -10%","-10 to -5%", "-5 to -2.5%", "-2.5 to 1%", "1 to 0%",
              "0 to 1%"," 1 to 2.5%", "2.5 to 5%", "5 to 10%", "10 to 15%", "> 15%")
  
  # Crear una columna adicional para definir el color del texto basado en el valor de pval
  data_df$text_color <- ifelse(data_df$dif_pocet_past_fut > 10,  "white", "black")
  
  ggplot(data_df, aes(x = Temperature, y = Precipitation, fill = cut(dif_pocet_past_fut, breaks = c(-Inf, breaks, Inf), include.lowest = TRUE))) +
    geom_tile(color = "black", size = 0.1) +
    scale_fill_manual(values = custom_palette, drop = FALSE, na.value = "white", name = "", labels = c(labels)) +
    geom_text(aes(label = pval, color = text_color), size = 4, vjust = 0.5) +  # Ajustar el color del texto según la nueva columna
    labs(title = title,
         x = "Temperature change (°C)",
         y = "Precipitation change (%)") +
    scale_x_discrete(labels = c("1", "2", "3", "4", "5", "6")) +
    scale_y_discrete(labels = c("-40", "-20", "0", "20", "40", "60")) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),  # Tamaño del título ajustado
      axis.title.x = element_text(size = 11),
      axis.title.y = element_text(size = 11),
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      legend.key.width = unit(3, "cm"),
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.5),
      axis.ticks.length = unit(0.2, "cm")
    ) +
    scale_color_identity()  # Utilizar la identidad de color para aplicar los colores directamente
}

# Función para crear el gráfico con leyenda
crear_grafico_con_leyenda <- function(data_df, title) {
  data_df <- truncar_cambio(data_df)
  labels <- c("< -15%", "-15 to -10%","-10 to -5%", "-5 to -2.5%", "-2.5 to 1%", "1 to 0%",
              "0 to 1%"," 1 to 2.5%", "2.5 to 5%", "5 to 10%", "10 to 15%", "> 15%")
  
  ggplot(data_df, aes(x = Temperature, y = Precipitation, fill = cut(dif_pocet_past_fut, breaks = c(-Inf, breaks, Inf), include.lowest = TRUE))) +
    geom_tile(color = "black", size = 0.1) +
    scale_fill_manual(values = custom_palette, drop = FALSE, na.value = "white", name = "", labels = c(labels),guide = guide_legend(ncol = 6)) +
    geom_text(aes(label = pval), color = "black", size = 1, vjust = 0.5) +
    labs(title = title,
         x = "Temperature change (°C)",
         y = "Precipitation change (%)") +
    scale_x_discrete(labels = c("1", "2", "3", "4", "5", "6")) +
    scale_y_discrete(labels = c("-40", "-20", "0", "20", "40", "60")) +
    theme_minimal(base_size = 15) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 15),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom",
      legend.key.width = unit(1.5, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black"),
      axis.line.x.top = element_blank(),
      axis.line.y.right = element_blank(),
      axis.ticks.length = unit(0.2, "cm")
    )
}

# Función personalizada para extraer la leyenda
get_legend <- function(myggplot) {
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}

# Cargar los archivos CSV y convertir Change a numérico
load_and_prepare_data <- function(file_path) {
  data_df <- read.csv(file_path)
  data_df$Change <- as.numeric(data_df$Change)
  return(data_df)
}

syn <- "percentile95"

if(syn == "percentile95"){
  synp = "p95"
}else{
  synp = syn
}

dir_data <- "/home/miguel/Dropbox/synchronicity/csvs/"
dir_out <-"/home/miguel/Dropbox/synchronicity/csvs/"
# Cargar los archivos CSV
djf_df_neu <- read.csv(paste0(dir_data,"dif_DJF_",syn,"_NEU_resultados_dif_change.csv", sep=""))
djf_df_neu$dif_pocet_past_fut<-(djf_df_neu$porcent_fut - djf_df_neu$porcent_past)
mam_df_neu <- read.csv(paste0(dir_data,"dif_MAM_",syn,"_NEU_resultados_dif_change.csv", sep=""))
mam_df_neu$dif_pocet_past_fut<-(mam_df_neu$porcent_fut - mam_df_neu$porcent_past)
jja_df_neu <- read.csv(paste0(dir_data,"dif_JJA_",syn,"_NEU_resultados_dif_change.csv", sep=""))
jja_df_neu$dif_pocet_past_fut<-(jja_df_neu$porcent_fut - jja_df_neu$porcent_past)
son_df_neu <- read.csv(paste0(dir_data,"dif_SON_",syn,"_NEU_resultados_dif_change.csv", sep=""))
son_df_neu$dif_pocet_past_fut<-(son_df_neu$porcent_fut - son_df_neu$porcent_past)

djf_df_ceu <- read.csv(paste0(dir_data,"dif_DJF_",syn,"_CEU_resultados_dif_change.csv", sep=""))
djf_df_ceu$dif_pocet_past_fut<-(djf_df_ceu$porcent_fut - djf_df_ceu$porcent_past)
mam_df_ceu <- read.csv(paste0(dir_data,"dif_MAM_",syn,"_CEU_resultados_dif_change.csv", sep=""))
mam_df_ceu$dif_pocet_past_fut<-(mam_df_ceu$porcent_fut - mam_df_ceu$porcent_past)
jja_df_ceu <- read.csv(paste0(dir_data,"dif_JJA_",syn,"_CEU_resultados_dif_change.csv", sep=""))
jja_df_ceu$dif_pocet_past_fut<-(jja_df_ceu$porcent_fut - jja_df_ceu$porcent_past)
son_df_ceu <- read.csv(paste0(dir_data,"dif_SON_",syn,"_CEU_resultados_dif_change.csv", sep=""))
son_df_ceu$dif_pocet_past_fut<-(son_df_ceu$porcent_fut - son_df_ceu$porcent_past)

djf_df_med <- read.csv(paste0(dir_data,"dif_DJF_",syn,"_MED_resultados_dif_change.csv", sep=""))
djf_df_med$dif_pocet_past_fut<-(djf_df_med$porcent_fut - djf_df_med$porcent_past)
mam_df_med <- read.csv(paste0(dir_data,"dif_MAM_",syn,"_MED_resultados_dif_change.csv", sep=""))
mam_df_med$dif_pocet_past_fut<-(mam_df_med$porcent_fut - mam_df_med$porcent_past)
jja_df_med <- read.csv(paste0(dir_data,"dif_JJA_",syn,"_MED_resultados_dif_change.csv", sep=""))
jja_df_med$dif_pocet_past_fut<-(jja_df_med$porcent_fut - jja_df_med$porcent_past)
son_df_med <- read.csv(paste0(dir_data,"dif_SON_",syn,"_MED_resultados_dif_change.csv", sep=""))
son_df_med$dif_pocet_past_fut<-(son_df_med$porcent_fut - son_df_med$porcent_past)

change_vector <- c(
  djf_df_neu$dif_pocet_past_fut,
  mam_df_neu$dif_pocet_past_fut,
  jja_df_neu$dif_pocet_past_fut,
  son_df_neu$dif_pocet_past_fut,
  djf_df_ceu$dif_pocet_past_fut,
  mam_df_ceu$dif_pocet_past_fut,
  jja_df_ceu$dif_pocet_past_fut,
  son_df_ceu$dif_pocet_past_fut,
  djf_df_med$dif_pocet_past_fut,
  mam_df_med$dif_pocet_past_fut,
  jja_df_med$dif_pocet_past_fut,
  son_df_med$dif_pocet_past_fut
)

percentiles <- c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.90, 0.95)
percentil_values <- quantile(change_vector, percentiles, na.rm = TRUE)

# Crear los gráficos
djf_plot_neu <- crear_grafico(djf_df_neu, paste0("Change DJF (NEU; FWI > ",synp,")", sep=""))
mam_plot_neu <- crear_grafico(mam_df_neu, paste0("Change MAM (NEU; FWI > ",synp,")", sep=""))
jja_plot_neu <- crear_grafico(jja_df_neu, paste0("Change JJA (NEU; FWI > ",synp,")", sep=""))
son_plot_neu <- crear_grafico(son_df_neu, paste0("Change SON (NEU; FWI > ",synp,")", sep=""))

djf_plot_ceu <- crear_grafico(djf_df_ceu, paste0("Change DJF (CEU; FWI > ",synp,")", sep=""))
mam_plot_ceu <- crear_grafico(mam_df_ceu, paste0("Change MAM (CEU; FWI > ",synp,")", sep=""))
jja_plot_ceu <- crear_grafico(jja_df_ceu, paste0("Change JJA (CEU; FWI > ",synp,")", sep=""))
son_plot_ceu <- crear_grafico(son_df_ceu, paste0("Change SON (CEU; FWI > ",synp,")", sep=""))

djf_plot_med <- crear_grafico(djf_df_med, paste0("Change DJF (MED; FWI > ",synp,")", sep=""))
mam_plot_med <- crear_grafico(mam_df_med, paste0("Change MAM (MED; FWI > ",synp,")", sep=""))
jja_plot_med <- crear_grafico(jja_df_med, paste0("Change JJA (MED; FWI > ",synp,")", sep=""))
son_plot_med <- crear_grafico(son_df_med, paste0("Change SON (MED; FWI > ",synp,")", sep=""))

leyenda <- get_legend(crear_grafico_con_leyenda(djf_df_med, ""))
grafico_final <- plot_grid(djf_plot_neu, mam_plot_neu, jja_plot_neu, son_plot_neu,
                           djf_plot_ceu, mam_plot_ceu, jja_plot_ceu, son_plot_ceu,
                           djf_plot_med, mam_plot_med, jja_plot_med, son_plot_med,
                           ncol = 4, align = "hv")
grafico_final_con_leyenda <- plot_grid(grafico_final, leyenda, ncol = 1, rel_heights = c(1, 0.1))

print(grafico_final_con_leyenda)
ggsave(paste0(dir_out, "sup_diff_resultados_",syn,"_change.pdf",sep=""), grafico_final_con_leyenda, width = 15, height = 9)

