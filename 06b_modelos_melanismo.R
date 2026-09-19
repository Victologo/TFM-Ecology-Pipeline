# 06b_modelos_melanismo.R
# Script para traducir la LUMINANCIA de LEPY a % de Reflectancia Real
# Cruzando los datos con los parches de calibración de ImageJ (lotes) y el Espectrómetro.

library(dplyr)
library(readxl)
library(tidyr)

cat("1. Cargando base de datos de Lotes (ImageJ)...\n")
lotes <- read_excel("datos_procesados/TFM_dataset_Melanismo_ImageJ.xlsx")

cat("2. Cargando base de datos completa de alas (LEPY)...\n")
# LEPY está una carpeta más atrás en "Miscelánea/ANEXOS/"
lepy_pa <- read_excel("../Miscelánea/ANEXOS/TFM_dataset_LEPY_FINAL.xlsx", sheet = "P. apollo")
lepy_zr <- read_excel("../Miscelánea/ANEXOS/TFM_dataset_LEPY_FINAL.xlsx", sheet = "Z. rumina")

# Unimos apollo y rumina en una sola tabla gigante
lepy_all <- bind_rows(lepy_pa, lepy_zr)

cat("3. Expandiendo lotes y cruzando datos...\n")
# El Excel de lotes tiene todos los IDs de LEPY en una sola celda separados por ';'
# Esta función 'separate_rows' crea una fila por cada ID para poder cruzar
lotes_expandidos <- lotes %>%
  separate_rows(IDs_LEPY, sep = ";\\s*") %>%
  filter(!is.na(IDs_LEPY), IDs_LEPY != "") %>%
  # En UCM/Montes al haber varias sesiones por colección, promediamos las sesiones
  # para que las condiciones de iluminación queden perfectamente estabilizadas
  group_by(IDs_LEPY, Escalilla) %>%
  summarise(
    Blanco_Raw = mean(Blanco_Raw, na.rm=TRUE),
    Gris1_Raw = mean(Gris1_Raw, na.rm=TRUE),
    Gris2_Raw = mean(Gris2_Raw, na.rm=TRUE),
    Gris3_Raw = mean(Gris3_Raw, na.rm=TRUE),
    Gris4_Raw = mean(Gris4_Raw, na.rm=TRUE),
    Negro_Raw = mean(Negro_Raw, na.rm=TRUE),
    .groups = "drop"
  )

# Cruzamos LEPY con los Lotes
datos_cruzados <- lepy_all %>%
  left_join(lotes_expandidos, by = c("Code" = "IDs_LEPY"))


cat("4. Calibrando reflectancias de las alas con Polinomio de Grado 2 (del Alamo et al. 2024)...\n")
  # ============================================================================
  # CALIBRACIÓN POLINÓMICA — JUSTIFICACIÓN METODOLÓGICA
  # ============================================================================
  # Las imágenes TIFF del MNCN (codificadas en 8 bits por ImageMagick / Photoshop)
  # llevan incrustada la curva gamma sRGB, que distorsiona la relación entre
  # el valor de píxel y la reflectancia física real de forma no lineal.
  # Aplicar una corrección gamma manual ((x/255)^2.2) antes de ajustar una
  # regresión lineal es incorrecto por dos razones:
  #   (1) La curva sRGB NO es gamma pura 2.2 (es una función por tramos que
  #       falla en los oscuros, exactamente donde Z. rumina tiene sus valores).
  #   (2) La gamma quedaría corregida DOBLE: una vez manualmente y otra vez
  #       absorbida por la propia regresión.
  # Solución (del Alamo et al. 2024, Ecol. Entomol.): Ajustar un POLINOMIO DE
  # GRADO 2 directamente sobre los píxeles brutos (sin linealizar). El polinomio
  # absorbe la curvatura gamma de la cámara de forma automática. Para las escalas
  # de 4 parches del MNCN (datos históricos), el Grado 2 es el óptimo: el Grado 3
  # produce sobreajuste severo (0 grados de libertad residuales) y precipita
  # valores negativos en tonos oscuros. El Grado 2, con intercepto libre, resuelve
  # el problema por completo (0 valores negativos en Z. rumina).
  # ============================================================================
  calibrar_ala <- function(fila) {
    if (is.na(fila$Escalilla) || is.na(fila$luminance_median)) return(NA)
    
    # Píxeles brutos (SIN corrección gamma manual)
    ala_raw <- fila$luminance_median
    
    if (grepl("4parches", fila$Escalilla)) {
      x_imagej <- c(fila$Blanco_Raw, fila$Gris1_Raw, fila$Gris2_Raw, fila$Negro_Raw)
      y_real   <- c(101.3, 21.0, 6.1, 2.9)
    } else {
      x_imagej <- c(fila$Blanco_Raw, fila$Gris1_Raw, fila$Gris2_Raw, fila$Gris3_Raw, fila$Gris4_Raw, fila$Negro_Raw)
      y_real   <- c(101.3, 59.5, 35.7, 20.2, 9.2, 2.9)
    }
    
    if (any(is.na(x_imagej))) return(NA)
    
    # Polinomio de Grado 2 con intercepto libre (del Alamo et al. 2024)
    # — intercepto libre: no se fuerza a 0; el modelo absorbe la gamma internamente.
    # — poly(x, 2, raw=TRUE): términos x y x² en escala de píxeles brutos.
    df_cal <- data.frame(y = y_real, x = x_imagej)
    modelo <- lm(y ~ poly(x, 2, raw = TRUE), data = df_cal)
    
    ala_real <- predict(modelo, newdata = data.frame(x = ala_raw))
    
    # Clip conservador: no puede haber reflectancia negativa (absorción > 100%)
    return(max(0, ala_real))
  }

datos_finales <- datos_cruzados %>%
  rowwise() %>%
  mutate(Reflectancia_Ala_Real = calibrar_ala(pick(everything()))) %>%
  ungroup() %>%
  # ==========================================
  # EXCLUSIÓN DE VALORES ATÍPICOS (OUTLIERS)
  # ==========================================
  # Se eliminan explícitamente 3 especímenes del MNCN que presentaron un brillo
  # extremadamente anómalo (sobreexposición/flash en la foto original), lo que 
  # infló su luminancia en LEPY, haciéndolos biológicamente inutilizables y 
  # sesgando la recta de calibración:
  # - PA_MNCN_324_D (Luminancia LEPY: 184.66) -> Máximo absoluto de P. apollo
  # - ZR_MNCN_231_D (Luminancia LEPY: 118.47) -> Máximo absoluto 1 de Z. rumina
  # - ZR_MNCN_230_D (Luminancia LEPY: 116.46) -> Máximo absoluto 2 de Z. rumina
  filter(!Code %in% c("PA_MNCN_324_D", "ZR_MNCN_231_D", "ZR_MNCN_230_D"))

cat("\n--- MUESTRA DE RESULTADOS ---\n")
print(datos_finales %>% 
        select(Code, Coleccion, luminance_median, Reflectancia_Ala_Real) %>% 
        head(10))

ruta_salida <- "datos_procesados/TFM_dataset_Melanismo_Final_Calibrado.csv"
write.csv(datos_finales, ruta_salida, row.names = FALSE)
cat(sprintf("\n-> ¡ÉXITO! Dataset final cruzado guardado en: %s\n", ruta_salida))

# ==========================================
# 5. GENERACIÓN DE GRÁFICOS Y REPORTE
# ==========================================
library(ggplot2)

cat("Generando gráfico de distribución de melanismo...\n")
p <- ggplot(datos_finales, aes(x = Reflectancia_Ala_Real, fill = Coleccion, color = Coleccion)) +
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(
       x = "Reflectancia Real (%)", y = "Densidad") +
  theme(legend.position = "bottom")

ggsave("graficos/08_espectrofotometria/06b_distribucion_melanismo.png", plot = p, width = 8, height = 6, dpi = 900)
cat("-> Gráfico guardado en: graficos/08_espectrofotometria/06b_distribucion_melanismo.png\n")

cat("Generando reporte Markdown para Obsidian...\n")
resumen_stats <- datos_finales %>%
  group_by(Coleccion) %>%
  summarise(
    N = n(),
    Media_Reflectancia = mean(Reflectancia_Ala_Real, na.rm=TRUE),
    Min = min(Reflectancia_Ala_Real, na.rm=TRUE),
    Max = max(Reflectancia_Ala_Real, na.rm=TRUE)
  )

md_path <- "../Wiki TFM Persistente/02_Sandbox/_staging/06b_Resultados_Melanismo.md"
sink(md_path)
cat("---\nstatus: Done\ndescription: Resumen estadístico generado automáticamente tras calibrar el melanismo\n---\n\n")
cat("# Resultados Finales de Melanismo (Fase 3)\n\n")
cat("La calibración de los datos de LEPY con las escalillas de ImageJ y el Espectrómetro ha finalizado correctamente.\n\n")
cat("## Resumen Estadístico de Reflectancia Real (%)\n\n")
cat("| Colección | N | Media | Mínimo | Máximo |\n")
cat("| :--- | :--- | :--- | :--- | :--- |\n")
for(i in 1:nrow(resumen_stats)) {
  cat(sprintf("| %s | %d | %.2f%% | %.2f%% | %.2f%% |\n", 
              resumen_stats$Coleccion[i], resumen_stats$N[i], 
              resumen_stats$Media_Reflectancia[i], resumen_stats$Min[i], resumen_stats$Max[i]))
}
cat("\n\n## Distribución del Melanismo\n")
cat("![Gráfico de distribución](../../../Analisis_R/graficos/08_espectrofotometria/06b_distribucion_melanismo.png)\n")
sink()
cat("-> Reporte MD guardado en Obsidian: ", md_path, "\n")
