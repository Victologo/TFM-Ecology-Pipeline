# ==============================================================================
# Script 01: Preparación, Unificación de Datos y Visualización para el TFM
# ==============================================================================
# Este script realiza la preparación de datos de captura, morfometría (LEPY) y
# clima (ClimateDT), calcula las ventanas ecológicas y genera las primeras
# visualizaciones de tendencias alares para Parnassius apollo y Zerynthia rumina.
#
# Utiliza rutas relativas, funcionando idénticamente en PC sobremesa y portátil.

# 1. Cargar librerías necesarias
library(readxl)
library(tidyverse)

# 2. Definición de rutas relativas
ruta_combinado <- "../Miscelánea/ANEXOS/TFM_dataset_maestro_con_lepy.xlsx"

# 3. Funciones de curación y normalización
# Limpieza de nombres para resolver discrepancias de encoding de ClimateDT
clean_names <- function(x) {
  gsub("[^a-zA-Z0-9_ -]", "", x)
}

# Conversión y limpieza de la altitud de captura a formato numérico
parse_altitude <- function(x) {
  x_char <- as.character(x)
  # Quitar puntos que actúan como separador de miles (ej. 1.500 m -> 1500)
  x_clean <- gsub("\\.", "", x_char)
  # Quedarse con el primer valor en rangos (ej. 600-667 -> 600)
  x_clean <- str_split(x_clean, "[-/]") %>% map_chr(1)
  # Extraer solo dígitos
  x_digits <- gsub("[^0-9]", "", x_clean)
  as.numeric(x_digits)
}

# ==============================================================================
# SECCIÓN A: PROCESAMIENTO DE ZERYNTHIA RUMINA
# ==============================================================================
print("Cargando y procesando Zerynthia rumina...")

# A.1. Cargar desde dataset combinado (ya tiene MAESTRO + LEPY unidos, col B = code_lepy correcto)
maestro_zr <- read_excel(ruta_combinado, sheet = "Zerynthia_rumina") %>%
  mutate(
    localidad_clean = clean_names(localidad),
    altitud_num = parse_altitude(altitud_m)
  )

# A.2. Calcular forewing_length a partir de las columnas LEPY ya presentes en el combinado
# Solo usar filas con lepy_status == "valid" y cara dorsal (code_lepy termina en _D)
dataset_final_zr <- maestro_zr %>%
  filter(lepy_status == "valid") %>%
  mutate(
    forewing_length = case_when(
      is.na(poi_dist_inner_outer_l) | poi_dist_inner_outer_l == 0 ~ poi_dist_inner_outer_r,
      is.na(poi_dist_inner_outer_r) | poi_dist_inner_outer_r == 0 ~ poi_dist_inner_outer_l,
      abs(poi_dist_inner_outer_l - poi_dist_inner_outer_r) <= 2.5 ~ (poi_dist_inner_outer_l + poi_dist_inner_outer_r) / 2,
      TRUE ~ pmax(poi_dist_inner_outer_l, poi_dist_inner_outer_r)
    ),
    wingspan_estimated = contour_width_calibrated
  )

# ==============================================================================
# SECCIÓN B: PROCESAMIENTO DE PARNASSIUS APOLLO
# ==============================================================================
print("Cargando y procesando Parnassius apollo...")

# B.1. Cargar desde dataset combinado
maestro_pa <- read_excel(ruta_combinado, sheet = "Parnassius_apollo") %>%
  mutate(
    localidad_clean = clean_names(localidad_2),
    altitud_num = parse_altitude(altitud_m)
  )

# B.2. Calcular forewing_length a partir de columnas LEPY ya presentes
dataset_final_pa <- maestro_pa %>%
  filter(lepy_status == "valid") %>%
  filter(id_specimen != "PA_MNCN_161" & code_lepy != "PA_MNCN_161_D") %>%
  mutate(
    forewing_length = case_when(
      is.na(poi_dist_inner_outer_l) | poi_dist_inner_outer_l == 0 ~ poi_dist_inner_outer_r,
      is.na(poi_dist_inner_outer_r) | poi_dist_inner_outer_r == 0 ~ poi_dist_inner_outer_l,
      abs(poi_dist_inner_outer_l - poi_dist_inner_outer_r) <= 2.5 ~ (poi_dist_inner_outer_l + poi_dist_inner_outer_r) / 2,
      TRUE ~ pmax(poi_dist_inner_outer_l, poi_dist_inner_outer_r)
    ),
    wingspan_estimated = contour_width_calibrated
  )


# ==============================================================================
# SECCIÓN C: EXPORTACIÓN DE DATASETS CONSOLIDADOS
# ==============================================================================
if (!dir.exists("datos_procesados")) {
  dir.create("datos_procesados")
}
write_csv(dataset_final_pa, "datos_procesados/dataset_final_apollo.csv")
write_csv(dataset_final_zr, "datos_procesados/dataset_final_rumina.csv")

# Identificación de Outliers Estadísticos basados en IQR por Colección
detectar_outliers_iqr <- function(df, especie_name) {
  df %>%
    filter(!is.na(forewing_length) & !is.na(fuente)) %>%
    group_by(fuente) %>%
    mutate(
      q25 = quantile(forewing_length, 0.25, na.rm = TRUE),
      q75 = quantile(forewing_length, 0.75, na.rm = TRUE),
      iqr = q75 - q25,
      lower_limit = q25 - 1.5 * iqr,
      upper_limit = q75 + 1.5 * iqr,
      es_outlier = forewing_length < lower_limit | forewing_length > upper_limit
    ) %>%
    filter(es_outlier) %>%
    select(id_specimen, code_lepy, fuente, forewing_length, lower_limit, upper_limit) %>%
    mutate(Especie = especie_name) %>%
    ungroup()
}

outliers_pa <- detectar_outliers_iqr(dataset_final_pa, "Parnassius apollo")
outliers_zr <- detectar_outliers_iqr(dataset_final_zr, "Zerynthia rumina")
outliers_total <- bind_rows(outliers_pa, outliers_zr)

# Escribir a CSV para que el usuario pueda consultarlo directamente
write_csv(outliers_total, "datos_procesados/outliers_detectados.csv")

# Mostrar reporte en la consola de R
print("=== REPORTANDO OUTLIERS ESTADÍSTICOS DETECTADOS (IQR) ===")
if (nrow(outliers_total) > 0) {
  print(as.data.frame(outliers_total))
  print(paste("Se han guardado", nrow(outliers_total), "outliers en: datos_procesados/outliers_detectados.csv"))
} else {
  print("No se ha detectado ningún outlier en el dataset.")
}

print("¡Datos consolidados y exportados a la carpeta datos_procesados/!")

# --- REPORTE ESTADÍSTICO COMPLETO DE LOS DATASETS FINALES ---
cat("\n============================================================\n")
cat(" RESUMEN DE DATASETS FINALES (Script 01)\n")
cat("============================================================\n")

cat("\n--- PARNASSIUS APOLLO ---\n")
cat("N total en maestro:", nrow(maestro_pa), "\n")
cat("N con forewing_length (LEPY válido):", sum(!is.na(dataset_final_pa$forewing_length)), "\n")
cat("N sin forewing_length (sin match LEPY):", sum(is.na(dataset_final_pa$forewing_length)), "\n")
print(table(dataset_final_pa$fuente, useNA = "always"))
cat("\nEstadísticos de forewing_length (mm):\n")
print(summary(dataset_final_pa$forewing_length))
cat("SD:", round(sd(dataset_final_pa$forewing_length, na.rm = TRUE), 3), "\n")
cat("\nDistribución por sexo:\n")
print(table(dataset_final_pa$sexo, useNA = "always"))
cat("\nRango de años:\n")
print(range(dataset_final_pa$año, na.rm = TRUE))
cat("\nRango de altitud (m):\n")
print(range(dataset_final_pa$altitud_num, na.rm = TRUE))

cat("\n--- ZERYNTHIA RUMINA ---\n")
cat("N total en maestro:", nrow(maestro_zr), "\n")
cat("N con forewing_length (LEPY válido):", sum(!is.na(dataset_final_zr$forewing_length)), "\n")
cat("N sin forewing_length (sin match LEPY):", sum(is.na(dataset_final_zr$forewing_length)), "\n")
print(table(dataset_final_zr$fuente, useNA = "always"))
cat("\nEstadísticos de forewing_length (mm):\n")
print(summary(dataset_final_zr$forewing_length))
cat("SD:", round(sd(dataset_final_zr$forewing_length, na.rm = TRUE), 3), "\n")
cat("\nDistribución por sexo:\n")
print(table(dataset_final_zr$sexo, useNA = "always"))
cat("\nRango de años:\n")
print(range(dataset_final_zr$año, na.rm = TRUE))
cat("\nRango de altitud (m):\n")
print(range(dataset_final_zr$altitud_num, na.rm = TRUE))
cat("============================================================\n\n")


