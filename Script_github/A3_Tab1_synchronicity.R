rm(list = ls())
graphics.off()
gc()

library(ggplot2)
library(reshape2)
library(knitr)
library(kableExtra)
library(rmarkdown)

valores_p95 <- c(
  "ns", "ns", "ns", "ns",  # NEU: DJF, MAM, JJA, SON
  "ns", 73, 110, 90,      # CEU: DJF, MAM, JJA, SON
  "ns", 44, 66, 32        # MED: DJF, MAM, JJA, SON
)
matriz_p95 <- matrix(valores_p95, nrow = 3, ncol = 4, byrow = TRUE,
                     dimnames = list(c("NEU", "CEU", "MED"), c("DJF", "MAM", "JJA", "SON")))

# Matriz para Figure1_50.pdf
valores_50 <- c(
  "ns",  "ns",  "ns",  "ns",  # NEU: DJF, MAM, JJA, SON
  "ns", 314, 271, 389,   # CEU: DJF, MAM, JJA, SON
    47, 23, 19,  "ns"         # MED: DJF, MAM, JJA, SON
)
matriz_50 <- matrix(valores_50, nrow = 3, ncol = 4, byrow = TRUE,
                    dimnames = list(c("NEU", "CEU", "MED"), c("DJF", "MAM", "JJA", "SON")))

# Matriz para Figure1_38.pdf
valores_38 <- c(
  "ns",  "ns", "ns", "ns", # NEU: DJF, MAM, JJA, SON
  "ns", 224, 160, 260,    # CEU: DJF, MAM, JJA, SON
    33, 20, 16,  "ns"         # MED: DJF, MAM, JJA, SON
)
matriz_38 <- matrix(valores_38, nrow = 3, ncol = 4, byrow = TRUE,
                    dimnames = list(c("NEU", "CEU", "MED"), c("DJF", "MAM", "JJA", "SON")))

# Crear una matriz combinada con las estaciones en las columnas y las regiones en las filas
matriz_combinada <- cbind(
  matriz_50[, "DJF"], matriz_38[, "DJF"], matriz_p95[, "DJF"],
  matriz_50[, "MAM"], matriz_38[, "MAM"], matriz_p95[, "MAM"],
  matriz_50[, "JJA"], matriz_38[, "JJA"], matriz_p95[, "JJA"],
  matriz_50[, "SON"], matriz_38[, "SON"], matriz_p95[, "SON"]
)

# Nombrar las filas y columnas de la matriz combinada
rownames(matriz_combinada) <- c("NEU", "CEU", "MED")
colnames(matriz_combinada) <- c("50", "38", "p95", "50", "38", "p95", "50", "38", "p95", "50", "38", "p95")

# Convertir la matriz en un data frame
df <- as.data.frame(matriz_combinada)

# dev.off()
# Crear la tabla bonita usando kableExtra
tabla <- kable(df, format = "html", caption = "", align = "c") %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"), 
                full_width = FALSE, font_size = 14) %>%
  add_header_above(c("." = 1, "DJF" = 3, "MAM" = 3, "JJA" = 3, "SON" = 3), 
                   background = "black", color = "white") %>%
  row_spec(0, bold = TRUE, background = "grey25", color = "white") %>%
  column_spec(1, bold = TRUE, background = "grey25", color = "white") %>%
  column_spec(2:13, width = "5em", background = "#F5F5F5") %>%
  kable_styling(font = "Times New Roman", position = "center") %>%
  add_header_above(c(" " = 1, " " = 12), background = "grey15", color = "white")

dir_out <- '/home/miguel/Dropbox/synchronicity/Figures/'
file_path <- file.path(dir_out, tabla)
pdf(file = file_path, width = 8.5, height = 6)

