# ==============================================================================
# Script 08: Modelos Estadísticos de Desfase Térmico y Altitud para Melanismo
# ==============================================================================
# Este script analiza la variación de la Reflectancia Alar Calibrada
# (Melanismo) en función del clima histórico (ClimateDT: Año 0 vs. Año -1),
# el gradiente altitudinal y el sexo, contrastando la hipótesis del "eco térmico"
# y la regla del melanismo térmico (Bogert) en P. apollo y Z. rumina.
#
# Genera modelos lineales mixtos (LMM), diagnósticos DHARMa, matrices de
# correlación, paneles gráficos estacionales y forest plots de coeficientes.
# ==============================================================================

print("=== INICIANDO MÓDULO 08: MODELOS DE MELANISMO Y DESFASE TÉRMICO ===")

# 1. Carga de Paquetes
suppressPackageStartupMessages({
  library(tidyverse)
  library(lme4)
  library(lmerTest)
  library(MuMIn)
  library(corrplot)
  library(cowplot)
  library(DHARMa)
})

# 2. Configuración de Directorios y Estilo Gráfico
if (!dir.exists("graficos/08_melanismo_clima")) {
  dir.create("graficos/08_melanismo_clima", recursive = TRUE)
}
if (!dir.exists("graficos/08_melanismo_clima/diagnosticos_dharma")) {
  dir.create("graficos/08_melanismo_clima/diagnosticos_dharma", recursive = TRUE)
}
if (!dir.exists("datos_procesados")) {
  dir.create("datos_procesados", recursive = TRUE)
}

tema_tfm <- theme_minimal(base_size = 13) +
  theme(
    plot.title = element_blank(),
    plot.subtitle = element_blank(),
    axis.title = element_text(face = "bold", size = 11),
    legend.position = "bottom",
    panel.background = element_rect(fill = "#FAFAFA", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

colores_sexo <- c("Macho" = "#2C7BB6", "Hembra" = "#D7191C", "Indeterminado" = "gray60")
colores_especie <- c("Parnassius apollo" = "#4A505A", "Zerynthia rumina" = "#E6A100")
opacidad_puntos <- 0.65

clean_names <- function(x) { gsub("[^a-zA-Z0-9_ -]", "", x) }

# ==============================================================================
# SECCIÓN 1: CARGA Y CRUCE DE DATOS
# ==============================================================================
print("-> Cargando datasets biológicos, fenotípicos y climáticos...")

# Especímenes atípicos por desgaste extremo a excluir
outliers_desgaste <- c(
  "ZR_MON_052_D", "ZR_MON_001_D", "ZR_MON_008_D", 
  "ZR_UCM_019_D", "ZR_UCM_020_D", "ZR_MNCN_319_D",
  "ZR_UCM_067_D", "ZR_MON_051_D", "ZR_UCM_100_D", "ZR_UCM_087_D",
  "PA_MNCN_161_D"
)

# 1.1 Dataset Maestro de Melanismo Calibrado (Fase 3)
df_melanismo <- read_csv("datos_procesados/TFM_dataset_Melanismo_Final_Calibrado.csv", show_col_types = FALSE) %>%
  filter(!Code %in% outliers_desgaste)

# 1.2 Datasets Biológicos Base (Script 01)
df_base_pa <- read_csv("datos_procesados/dataset_final_apollo.csv", show_col_types = FALSE)
df_base_zr <- read_csv("datos_procesados/dataset_final_rumina.csv", show_col_types = FALSE)

# 1.3 Series Climáticas Históricas (ClimateDT)
clima_pa <- read_csv("../Miscelánea/ANEXOS/ClimateDT_output_apollo.csv", show_col_types = FALSE) %>% 
  mutate(ID_clean = clean_names(ID))
clima_zr <- read_csv("../Miscelánea/ANEXOS/ClimateDT_output_rumina.csv", show_col_types = FALSE) %>% 
  mutate(ID_clean = clean_names(ID))

# Extensión de Pedriza para Z. rumina
clima_zr_pedriza <- clima_zr %>% 
  filter(ID == "La_Pedriza") %>% 
  mutate(ID_clean = "PedrizadeManzanares", ID = "Pedriza_de_Manzanares")
clima_zr_extended <- bind_rows(clima_zr, clima_zr_pedriza) %>% 
  distinct(ID_clean, Year, .keep_all = TRUE)
clima_pa_extended <- clima_pa %>% 
  distinct(ID_clean, Year, .keep_all = TRUE)

# ==============================================================================
# SECCIÓN 2: EXTRACCIÓN DE VENTANAS CLIMÁTICAS ESTACIONALES
# ==============================================================================
print("-> Extrayendo ventanas climáticas estacionales (Año 0 vs. Año -1)...")

# 2.1 Zerynthia rumina
extraer_clima_zr <- function(row) {
  loc <- clean_names(row$localidad_clean)
  yr <- as.numeric(row$año)
  if (is.na(loc) | is.na(yr)) {
    return(data.frame(T_AMJ1=NA, P_AMJ1=NA, T_JAS1=NA, P_JAS1=NA, T_OND1=NA, P_OND1=NA, T_EFM0=NA, P_EFM0=NA, T_AMJ0=NA, P_AMJ0=NA))
  }
  c_prev <- clima_zr_extended %>% filter(ID_clean == loc, Year == (yr - 1))
  c_curr <- clima_zr_extended %>% filter(ID_clean == loc, Year == yr)
  
  if (nrow(c_prev) == 0 | nrow(c_curr) == 0) {
    return(data.frame(T_AMJ1=NA, P_AMJ1=NA, T_JAS1=NA, P_JAS1=NA, T_OND1=NA, P_OND1=NA, T_EFM0=NA, P_EFM0=NA, T_AMJ0=NA, P_AMJ0=NA))
  }
  
  T_AMJ1 <- mean(as.numeric(c_prev[1, paste0(c(rep("tmn",3), rep("tmx",3)), c("04","05","06","04","05","06"))]), na.rm = TRUE)
  P_AMJ1 <- sum(as.numeric(c_prev[1, paste0("prc", c("04","05","06"))]), na.rm = TRUE)
  
  T_JAS1 <- mean(as.numeric(c_prev[1, paste0(c(rep("tmn",3), rep("tmx",3)), c("07","08","09","07","08","09"))]), na.rm = TRUE)
  P_JAS1 <- sum(as.numeric(c_prev[1, paste0("prc", c("07","08","09"))]), na.rm = TRUE)
  
  T_OND1 <- mean(as.numeric(c_prev[1, paste0(c(rep("tmn",3), rep("tmx",3)), c("10","11","12","10","11","12"))]), na.rm = TRUE)
  P_OND1 <- sum(as.numeric(c_prev[1, paste0("prc", c("10","11","12"))]), na.rm = TRUE)
  
  T_EFM0 <- mean(as.numeric(c_curr[1, paste0(c(rep("tmn",3), rep("tmx",3)), c("01","02","03","01","02","03"))]), na.rm = TRUE)
  P_EFM0 <- sum(as.numeric(c_curr[1, paste0("prc", c("01","02","03"))]), na.rm = TRUE)
  
  T_AMJ0 <- mean(as.numeric(c_curr[1, paste0(c(rep("tmn",3), rep("tmx",3)), c("04","05","06","04","05","06"))]), na.rm = TRUE)
  P_AMJ0 <- sum(as.numeric(c_curr[1, paste0("prc", c("04","05","06"))]), na.rm = TRUE)
  
  return(data.frame(T_AMJ1, P_AMJ1, T_JAS1, P_JAS1, T_OND1, P_OND1, T_EFM0, P_EFM0, T_AMJ0, P_AMJ0))
}

clima_zr_lista <- lapply(1:nrow(df_base_zr), function(i) extraer_clima_zr(df_base_zr[i, ]))
df_base_zr <- bind_cols(df_base_zr, bind_rows(clima_zr_lista))

# 2.2 Parnassius apollo
extraer_clima_pa <- function(row) {
  loc <- clean_names(row$localidad_clean)
  yr <- as.numeric(row$año)
  if (is.na(loc) | is.na(yr)) {
    return(data.frame(T_JJA1=NA, P_JJA1=NA, T_SON1=NA, P_SON1=NA, T_DEF0=NA, P_DEF0=NA, T_MAM0=NA, P_MAM0=NA, T_JJA0=NA, P_JJA0=NA))
  }
  c_prev <- clima_pa_extended %>% filter(ID_clean == loc, Year == (yr - 1))
  c_curr <- clima_pa_extended %>% filter(ID_clean == loc, Year == yr)
  
  if (nrow(c_prev) == 0 | nrow(c_curr) == 0) {
    return(data.frame(T_JJA1=NA, P_JJA1=NA, T_SON1=NA, P_SON1=NA, T_DEF0=NA, P_DEF0=NA, T_MAM0=NA, P_MAM0=NA, T_JJA0=NA, P_JJA0=NA))
  }
  
  T_JJA1 <- mean(as.numeric(c_prev[1, paste0(c(rep("tmn",3), rep("tmx",3)), c("06","07","08","06","07","08"))]), na.rm = TRUE)
  P_JJA1 <- sum(as.numeric(c_prev[1, paste0("prc", c("06","07","08"))]), na.rm = TRUE)
  
  T_SON1 <- mean(as.numeric(c_prev[1, paste0(c(rep("tmn",3), rep("tmx",3)), c("09","10","11","09","10","11"))]), na.rm = TRUE)
  P_SON1 <- sum(as.numeric(c_prev[1, paste0("prc", c("09","10","11"))]), na.rm = TRUE)
  
  t_def0_vals <- c(as.numeric(c_prev[1, c("tmn12", "tmx12")]), as.numeric(c_curr[1, c("tmn01", "tmx01", "tmn02", "tmx02")]))
  p_def0_vals <- c(as.numeric(c_prev[1, "prc12"]), as.numeric(c_curr[1, c("prc01", "prc02")]))
  T_DEF0 <- mean(t_def0_vals, na.rm = TRUE)
  P_DEF0 <- sum(p_def0_vals, na.rm = TRUE)
  
  T_MAM0 <- mean(as.numeric(c_curr[1, paste0(c(rep("tmn",3), rep("tmx",3)), c("03","04","05","03","04","05"))]), na.rm = TRUE)
  P_MAM0 <- sum(as.numeric(c_curr[1, paste0("prc", c("03","04","05"))]), na.rm = TRUE)
  
  T_JJA0 <- mean(as.numeric(c_curr[1, paste0(c(rep("tmn",3), rep("tmx",3)), c("06","07","08","06","07","08"))]), na.rm = TRUE)
  P_JJA0 <- sum(as.numeric(c_curr[1, paste0("prc", c("06","07","08"))]), na.rm = TRUE)
  
  return(data.frame(T_JJA1, P_JJA1, T_SON1, P_SON1, T_DEF0, P_DEF0, T_MAM0, P_MAM0, T_JJA0, P_JJA0))
}

clima_pa_lista <- lapply(1:nrow(df_base_pa), function(i) extraer_clima_pa(df_base_pa[i, ]))
df_base_pa <- bind_cols(df_base_pa, bind_rows(clima_pa_lista))

# ==============================================================================
# SECCIÓN 3: UNIFICACIÓN CON REFLECTANCIA CALIBRADA
# ==============================================================================
print("-> Unificando reflectancia alar calibrada con metadatos y clima...")

# Cruce para Zerynthia rumina
df_zr_melanismo <- df_melanismo %>% 
  filter(grepl("^ZR", Code)) %>% 
  select(Code, Reflectancia_Ala_Real, Escalilla, Blanco_Raw, Gris1_Raw, Negro_Raw) %>% 
  inner_join(df_base_zr, by = c("Code" = "code_lepy")) %>% 
  mutate(
    sexo_factor = case_when(
      sexo %in% c("Macho", "Hembra") ~ sexo,
      TRUE ~ "Indeterminado"
    )
  )

# Cruce para Parnassius apollo
df_pa_melanismo <- df_melanismo %>% 
  filter(grepl("^PA", Code)) %>% 
  select(Code, Reflectancia_Ala_Real, Escalilla, Blanco_Raw, Gris1_Raw, Negro_Raw) %>% 
  inner_join(df_base_pa, by = c("Code" = "code_lepy")) %>% 
  mutate(
    sexo_factor = case_when(
      sexo %in% c("Macho", "Hembra") ~ sexo,
      TRUE ~ "Indeterminado"
    )
  )

# Guardar datasets consolidados de melanismo + clima
write_csv(df_zr_melanismo, "datos_procesados/dataset_melanismo_clima_rumina.csv")
write_csv(df_pa_melanismo, "datos_procesados/dataset_melanismo_clima_apollo.csv")

  # ==============================================================================
  # 1.5. VALIDACIÓN VISUAL: SESGO DE COLECCIONES (REFLECTANCIA)
  # ==============================================================================
  df_val_col <- bind_rows(
    df_zr_melanismo %>% mutate(Especie = "Zerynthia rumina",  fecha_captura = as.character(fecha_captura)),
    df_pa_melanismo %>% mutate(Especie = "Parnassius apollo", fecha_captura = as.character(fecha_captura))
  )
  
  p_val <- ggplot(df_val_col, aes(x = fuente, y = Reflectancia_Ala_Real, fill = fuente)) +
    geom_boxplot(alpha = 0.7, outlier.color = "red", outlier.size = 2) +
    facet_wrap(~ Especie, scales = "free_y") +
    tema_tfm +
    labs(x = "Colección Fotográfica", y = "Reflectancia Alar Calibrada (%)") +
    theme(legend.position = "none")
    
  ggsave("graficos/08_melanismo_clima/08_00_boxplot_colecciones_melanismo.png", p_val, width = 8, height = 5, dpi = 300)
  
cat("Muestras unificadas con Reflectancia Calibrada:\n")
cat(" - Z. rumina: Total =", nrow(df_zr_melanismo), 
    "| Sexados =", sum(df_zr_melanismo$sexo_factor != "Indeterminado"), 
    "| Indeterminados =", sum(df_zr_melanismo$sexo_factor == "Indeterminado"), "\n")
cat(" - P. apollo: Total =", nrow(df_pa_melanismo), 
    "| Sexados =", sum(df_pa_melanismo$sexo_factor != "Indeterminado"), 
    "| Indeterminados =", sum(df_pa_melanismo$sexo_factor == "Indeterminado"), "\n")

# ==============================================================================
# SECCIÓN 4: MODELIZACIÓN ESTADÍSTICA (LMM)
# ==============================================================================
print("-> Ajustando Modelos Lineales Mixtos (LMM) de Melanismo...")
ctrl <- lmerControl(check.conv.singular = "ignore")

# --- 4.1 Zerynthia rumina: Modelos Globales (Todos los especímenes) y Sexados ---
df_zr_mod_t <- df_zr_melanismo %>% 
  drop_na(Reflectancia_Ala_Real, altitud_num, T_AMJ1, T_JAS1, T_OND1, T_EFM0, T_AMJ0, localidad_clean, año)
df_zr_mod_p <- df_zr_melanismo %>% 
  drop_na(Reflectancia_Ala_Real, altitud_num, P_AMJ1, P_JAS1, P_OND1, P_EFM0, P_AMJ0, localidad_clean, año)

# Modelo Temperatura Total (Todos los individuos)
m_zr_temp_all <- lmer(Reflectancia_Ala_Real ~ fuente + altitud_num + año + T_AMJ1 + T_JAS1 + T_OND1 + T_EFM0 + T_AMJ0 + (1|localidad_clean), 
                      data = df_zr_mod_t, control = ctrl)

# Modelo Temperatura Sexados (Solo Machos y Hembras para evaluar Dimorfismo)
df_zr_mod_t_sex <- df_zr_mod_t %>% filter(sexo %in% c("Macho", "Hembra"))
m_zr_temp_sex <- lmer(Reflectancia_Ala_Real ~ sexo + fuente + altitud_num + año + T_AMJ1 + T_JAS1 + T_OND1 + T_EFM0 + T_AMJ0 + (1|localidad_clean), 
                      data = df_zr_mod_t_sex, control = ctrl)

# Modelo Precipitación Total
m_zr_prec_all <- lmer(Reflectancia_Ala_Real ~ fuente + altitud_num + año + P_AMJ1 + P_JAS1 + P_OND1 + P_EFM0 + P_AMJ0 + (1|localidad_clean), 
                      data = df_zr_mod_p, control = ctrl)

# Modelo Precipitación Sexados
df_zr_mod_p_sex <- df_zr_mod_p %>% filter(sexo %in% c("Macho", "Hembra"))
m_zr_prec_sex <- lmer(Reflectancia_Ala_Real ~ sexo + fuente + altitud_num + año + P_AMJ1 + P_JAS1 + P_OND1 + P_EFM0 + P_AMJ0 + (1|localidad_clean), 
                      data = df_zr_mod_p_sex, control = ctrl)


# --- 4.2 Parnassius apollo: Modelos de Temperatura y Precipitación ---
df_pa_mod_t <- df_pa_melanismo %>% 
  drop_na(Reflectancia_Ala_Real, altitud_num, T_JJA1, T_SON1, T_DEF0, T_MAM0, T_JJA0, localidad_clean, año)
df_pa_mod_p <- df_pa_melanismo %>% 
  drop_na(Reflectancia_Ala_Real, altitud_num, P_JJA1, P_SON1, P_DEF0, P_MAM0, P_JJA0, localidad_clean, año)

df_pa_mod_t_sex <- df_pa_mod_t %>% filter(sexo %in% c("Macho", "Hembra"))
m_pa_temp_sex <- lmer(Reflectancia_Ala_Real ~ sexo + fuente + altitud_num + año + T_JJA1 + T_SON1 + T_DEF0 + T_MAM0 + T_JJA0 + (1|localidad_clean), 
                      data = df_pa_mod_t_sex, control = ctrl)

df_pa_mod_p_sex <- df_pa_mod_p %>% filter(sexo %in% c("Macho", "Hembra"))
m_pa_prec_sex <- lmer(Reflectancia_Ala_Real ~ sexo + fuente + altitud_num + año + P_JJA1 + P_SON1 + P_DEF0 + P_MAM0 + P_JJA0 + (1|localidad_clean), 
                      data = df_pa_mod_p_sex, control = ctrl)


# Guardar resúmenes numéricos completos en Markdown
sink("graficos/08_melanismo_clima/08_01_LMM_Resultados_Melanismo.md")
cat("=========================================================================\n")
cat("  1. ZERYNTHIA RUMINA — LMM TEMPERATURA (MUESTRA COMPLETA, N =", nrow(df_zr_mod_t), ")\n")
cat("=========================================================================\n")
print(summary(m_zr_temp_all))
cat("\n\n=========================================================================\n")
cat("  2. ZERYNTHIA RUMINA — LMM TEMPERATURA CON SEXO (N =", nrow(df_zr_mod_t_sex), ")\n")
cat("=========================================================================\n")
print(summary(m_zr_temp_sex))
cat("\n\n=========================================================================\n")
cat("  3. ZERYNTHIA RUMINA — LMM PRECIPITACIÓN CON SEXO (N =", nrow(df_zr_mod_p_sex), ")\n")
cat("=========================================================================\n")
print(summary(m_zr_prec_sex))
cat("\n\n=========================================================================\n")
cat("  4. PARNASSIUS APOLLO — LMM TEMPERATURA CON SEXO (N =", nrow(df_pa_mod_t_sex), ")\n")
cat("=========================================================================\n")
print(summary(m_pa_temp_sex))
cat("\n\n=========================================================================\n")
cat("  5. PARNASSIUS APOLLO — LMM PRECIPITACIÓN CON SEXO (N =", nrow(df_pa_mod_p_sex), ")\n")
cat("=========================================================================\n")
print(summary(m_pa_prec_sex))
sink()

print("-> Resúmenes de modelos LMM guardados en graficos/08_melanismo_clima/08_01_LMM_Resultados_Melanismo.md")

# ==============================================================================
# SECCIÓN 5: DIAGNÓSTICO DE SUPUESTOS CON DHARMa
# ==============================================================================
print("-> Ejecutando diagnósticos DHARMa...")

modelos_dharma <- list(
  "ZR_Melanismo_Temp_Global" = m_zr_temp_all,
  "ZR_Melanismo_Temp_Sexados" = m_zr_temp_sex,
  "ZR_Melanismo_Prec_Sexados" = m_zr_prec_sex,
  "PA_Melanismo_Temp_Sexados" = m_pa_temp_sex,
  "PA_Melanismo_Prec_Sexados" = m_pa_prec_sex
)

sink("graficos/08_melanismo_clima/diagnosticos_dharma/DHARMa_Melanismo_Tests.md")
cat("# Diagnóstico de Supuestos LMM Melanismo — DHARMa\n\n")

for (nombre in names(modelos_dharma)) {
  cat("##", nombre, "\n")
  sim <- simulateResiduals(fittedModel = modelos_dharma[[nombre]], n = 1000, plot = FALSE)
  # Exportar gráfico de residuos (Desactivado por petición)
  # png(paste0("graficos/08_melanismo_clima/diagnosticos_dharma/08_00_dharma_", nombre, ".png"),
  #     width = 5400, height = 2700, res = 900)
  # plot(sim, main = paste("DHARMa Melanismo —", nombre))
  # dev.off()
  
  t_unif <- testUniformity(sim, plot = FALSE)
  t_disp <- testDispersion(sim, plot = FALSE)
  t_out  <- testOutliers(sim, plot = FALSE)
  
  cat("- **Uniformidad de residuos** (KS-test): D =", round(t_unif$statistic, 4),
      "| p =", round(t_unif$p.value, 4),
      ifelse(t_unif$p.value >= 0.05, "✅ OK", "⚠️ REVISAR"), "\n")
  cat("- **Dispersión** (ratio observado/esperado):", round(t_disp$statistic, 4),
      "| p =", round(t_disp$p.value, 4),
      ifelse(t_disp$p.value >= 0.05, "✅ OK", "⚠️ REVISAR"), "\n")
  cat("- **Outliers**: p =", round(t_out$p.value, 4),
      ifelse(t_out$p.value >= 0.05, "✅ OK", "⚠️ REVISAR"), "\n\n")
}
sink()

# ==============================================================================
# SECCIÓN 6: VISUALIZACIONES Y FOREST PLOTS
# ==============================================================================
print("-> Generando visualizaciones y paneles estacionales...")

# 6.1 Gráfico Altitudinal (Regla de Bogert / Melanismo Térmico)
p_alt_pa <- ggplot(df_pa_melanismo %>% filter(!is.na(altitud_num) & sexo_factor %in% c("Macho", "Hembra")), 
                   aes(x = altitud_num, y = Reflectancia_Ala_Real, color = sexo_factor)) +
  geom_point(alpha = opacidad_puntos, size = 2.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.2, linewidth = 1) +
  scale_color_manual(values = colores_sexo, name = "Sexo:") +
  labs(title = expression(bold("A") ~ " — " ~ bolditalic("Parnassius apollo")), 
       x = "Altitud de captura (m s.n.m.)",
       y = "Reflectancia alar real (%)") +
  tema_tfm +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0, margin = margin(b = 6)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.4),
    panel.background = element_rect(fill = "#FAFAFA", color = "gray80", linewidth = 0.5),
    legend.position = "none"
  )

p_alt_zr <- ggplot(df_zr_melanismo %>% filter(!is.na(altitud_num) & sexo_factor %in% c("Macho", "Hembra")), 
                   aes(x = altitud_num, y = Reflectancia_Ala_Real, color = sexo_factor)) +
  geom_point(alpha = opacidad_puntos, size = 2.2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.2, linewidth = 1) +
  scale_color_manual(values = colores_sexo, name = "Sexo:") +
  labs(title = expression(bold("B") ~ " — " ~ bolditalic("Zerynthia rumina")), 
       x = "Altitud de captura (m s.n.m.)",
       y = "Reflectancia alar real (%)") +
  tema_tfm +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0, margin = margin(b = 6)),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.4),
    panel.background = element_rect(fill = "#FAFAFA", color = "gray80", linewidth = 0.5),
    legend.position = "none"
  )

# Leyenda unificada inferior para el panel A y B
leg_altitud <- cowplot::get_legend(
  ggplot(df_pa_melanismo %>% filter(sexo_factor %in% c("Hembra", "Macho")), 
         aes(x = altitud_num, y = Reflectancia_Ala_Real, color = factor(sexo_factor, levels = c("Hembra", "Macho")))) +
    geom_point(size = 3) +
    scale_color_manual(values = colores_sexo, name = "Sexo:") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 11),
      legend.margin = margin(t = 4, b = 2)
    )
)

grid_altitud <- cowplot::plot_grid(p_alt_pa, p_alt_zr, ncol = 2, align = "h", rel_widths = c(1, 1))
panel_altitud <- cowplot::plot_grid(grid_altitud, leg_altitud, ncol = 1, rel_heights = c(1, 0.09))

# 6.1b Gráficos Individuales con título descriptivo propio
p_alt_pa_indiv <- p_alt_pa +
  labs(title = expression(bolditalic("Parnassius apollo")), subtitle = "Melanismo vs. Altitud de captura") +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 10, color = "gray30"))

p_alt_zr_indiv <- p_alt_zr +
  labs(title = expression(bolditalic("Zerynthia rumina")), subtitle = "Melanismo vs. Altitud de captura") +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 10, color = "gray30"))

# 6.1c Gráfico Tendencia Histórica (Melanismo vs Año)
p_hist_pa <- ggplot(df_pa_melanismo %>% filter(!is.na(año) & sexo_factor %in% c("Macho", "Hembra")), 
                   aes(x = año, y = Reflectancia_Ala_Real, color = sexo_factor)) +
  geom_point(alpha = opacidad_puntos, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.2) +
  scale_color_manual(values = colores_sexo, name = "Sexo:") +
  labs(title = expression(bolditalic("Parnassius apollo")), subtitle = "Tendencia histórica de reflectancia alar", x = "Año de captura",
       y = "Reflectancia alar real (%)") +
  tema_tfm +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.4),
    panel.background = element_rect(fill = "#FAFAFA", color = "gray80", linewidth = 0.5)
  )

p_hist_zr <- ggplot(df_zr_melanismo %>% filter(!is.na(año) & sexo_factor %in% c("Macho", "Hembra")), 
                   aes(x = año, y = Reflectancia_Ala_Real, color = sexo_factor)) +
  geom_point(alpha = opacidad_puntos, size = 2) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.2) +
  scale_color_manual(values = colores_sexo, name = "Sexo:") +
  labs(title = expression(bolditalic("Zerynthia rumina")), subtitle = "Tendencia histórica de reflectancia alar", x = "Año de captura",
       y = "Reflectancia alar real (%)") +
  tema_tfm +
  theme(
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 10, color = "gray30"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(color = "gray90", linewidth = 0.4),
    panel.background = element_rect(fill = "#FAFAFA", color = "gray80", linewidth = 0.5)
  )

# Guardar paneles y gráficos individuales
ggsave("graficos/08_melanismo_clima/08_02_Melanismo_vs_Altitud.png", plot = panel_altitud, width = 11, height = 5.2, dpi = 900)
ggsave("graficos/08_melanismo_clima/08_02_pa_altitud.png", plot = p_alt_pa_indiv, width = 7, height = 5, dpi = 900)
ggsave("graficos/08_melanismo_clima/08_02_zr_altitud.png", plot = p_alt_zr_indiv, width = 7, height = 5, dpi = 900)

ggsave("graficos/08_melanismo_clima/08_03_pa_tendencia_historica.png", plot = p_hist_pa, width = 7, height = 5, dpi = 900)
ggsave("graficos/08_melanismo_clima/08_03_zr_tendencia_historica.png", plot = p_hist_zr, width = 7, height = 5, dpi = 900)

# 6.2 Función para Gráficos de Dispersión Estacionales con Soporte de Indeterminados
crear_grafico_melanismo_estacional <- function(df, especie, est_name, var_col, tipo_var, color_esp, mostrar_indet = TRUE) {
  df_sub <- df %>% filter(!is.na(.data[[var_col]]) & !is.na(Reflectancia_Ala_Real))
  if (!mostrar_indet) {
    df_sub <- df_sub %>% filter(sexo_factor != "Indeterminado")
  }
  if (nrow(df_sub) < 5) return(NULL)
  
  eje_x_title <- if (tipo_var == "T") paste("Tª Media", est_name, "(°C)") else paste("Precipitación Acumulada", est_name, "(mm)")
  titulo <- if (tipo_var == "T") paste(especie, ": Reflectancia vs. Tª", est_name) else paste(especie, ": Reflectancia vs. Prec.", est_name)
  
  fit_global <- lm(as.formula(paste("Reflectancia_Ala_Real ~", var_col)), data = df_sub)
  r2_glob <- summary(fit_global)$r.squared
  p_glob <- summary(fit_global)$coefficients[2, 4]
  
  sub_text <- paste0("Global (N=", nrow(df_sub), "): R² = ", sprintf("%.3f", r2_glob), 
                     " (p=", ifelse(p_glob < 0.001, "<0.001", sprintf("%.3f", p_glob)), ")")
  
  g <- ggplot(df_sub, aes(x = .data[[var_col]], y = Reflectancia_Ala_Real)) +
    geom_point(aes(color = sexo_factor), alpha = opacidad_puntos, size = 1.9) +
    geom_smooth(method = "lm", formula = y ~ x, color = "gray20", se = TRUE, alpha = 0.2) +
    scale_color_manual(values = colores_sexo, name = "Sexo:") +
    labs(title = titulo, subtitle = sub_text, x = eje_x_title, y = "Reflectancia Real (%)") +
    tema_tfm
  
  return(g)
}

# 6.3 Paneles Estacionales de Z. rumina (Sexados + No sexados en todos los gráficos)
estaciones_zr <- list(
  c("Primavera(Año-1)", "T_AMJ1", "P_AMJ1"),
  c("Verano(Año-1)", "T_JAS1", "P_JAS1"),
  c("Otoño(Año-1)", "T_OND1", "P_OND1"),
  c("Invierno(Año 0)", "T_EFM0", "P_EFM0"),
  c("Primavera(Año 0)", "T_AMJ0", "P_AMJ0")
)

plots_zr_mel <- list()
for (est in estaciones_zr) {
  p_t <- crear_grafico_melanismo_estacional(df_zr_melanismo, "Z. rumina", est[1], est[2], "T", "#E6A100", mostrar_indet = FALSE)
  p_p <- crear_grafico_melanismo_estacional(df_zr_melanismo, "Z. rumina", est[1], est[3], "P", "#E6A100", mostrar_indet = FALSE)
  plots_zr_mel[[length(plots_zr_mel) + 1]] <- p_t
  plots_zr_mel[[length(plots_zr_mel) + 1]] <- p_p
}

panel_zr_mel_grid <- plot_grid(plotlist = plots_zr_mel, ncol = 2, labels = "AUTO")
title_zr_mel <- NULL
panel_zr_mel_final <- plot_grid(title_zr_mel, panel_zr_mel_grid, ncol = 1, rel_heights = c(0.03, 1))
ggsave("graficos/08_melanismo_clima/08_03_Zerynthia_rumina_Melanismo_PANEL.png", plot = panel_zr_mel_final, width = 14, height = 22, dpi = 900)

# 6.4 Paneles Estacionales de P. apollo
estaciones_pa <- list(
  c("Verano(Año-1)", "T_JJA1", "P_JJA1"),
  c("Otoño(Año-1)", "T_SON1", "P_SON1"),
  c("Invierno(Año 0)", "T_DEF0", "P_DEF0"),
  c("Primavera(Año 0)", "T_MAM0", "P_MAM0"),
  c("Verano(Año 0)", "T_JJA0", "P_JJA0")
)

plots_pa_mel <- list()
for (est in estaciones_pa) {
  p_t <- crear_grafico_melanismo_estacional(df_pa_melanismo, "P. apollo", est[1], est[2], "T", "#4A505A", mostrar_indet = FALSE)
  p_p <- crear_grafico_melanismo_estacional(df_pa_melanismo, "P. apollo", est[1], est[3], "P", "#4A505A", mostrar_indet = FALSE)
  plots_pa_mel[[length(plots_pa_mel) + 1]] <- p_t
  plots_pa_mel[[length(plots_pa_mel) + 1]] <- p_p
}

panel_pa_mel_grid <- plot_grid(plotlist = plots_pa_mel, ncol = 2, labels = "AUTO")
title_pa_mel <- NULL
panel_pa_mel_final <- plot_grid(title_pa_mel, panel_pa_mel_grid, ncol = 1, rel_heights = c(0.03, 1))
ggsave("graficos/08_melanismo_clima/08_04_Parnassius_apollo_Melanismo_PANEL.png", plot = panel_pa_mel_final, width = 14, height = 22, dpi = 900)

# 6.5 Forest Plot de Coeficientes LMM (Melanismo vs Clima)
print("-> Generando Forest Plots de Coeficientes de LMM para Melanismo...")

nombres_cronologicos <- c(
  "T_AMJ1" = "1. Primavera (Año-1)", "P_AMJ1" = "1. Primavera (Año-1)",
  "T_JAS1" = "2. Verano (Año-1)",    "P_JAS1" = "2. Verano (Año-1)",
  "T_OND1" = "3. Otoño (Año-1)",     "P_OND1" = "3. Otoño (Año-1)",
  "T_EFM0" = "4. Invierno (Año 0)",   "P_EFM0" = "4. Invierno (Año 0)",
  "T_AMJ0" = "5. Primavera (Año 0)",  "P_AMJ0" = "5. Primavera (Año 0)",
  
  "T_JJA1" = "1. Verano (Año-1)",    "P_JJA1" = "1. Verano (Año-1)",
  "T_SON1" = "2. Otoño (Año-1)",     "P_SON1" = "2. Otoño (Año-1)",
  "T_DEF0" = "3. Invierno (Año 0)",   "P_DEF0" = "3. Invierno (Año 0)",
  "T_MAM0" = "4. Primavera (Año 0)",  "P_MAM0" = "4. Primavera (Año 0)",
  "T_JJA0" = "5. Verano (Año 0)",    "P_JJA0" = "5. Verano (Año 0)"
)

extraer_coefs_mel <- function(model, especie, tipo) {
  s <- summary(model)
  coefs <- as.data.frame(s$coefficients)
  coefs$Variable <- rownames(coefs)
  coefs <- coefs %>% filter(!Variable %in% c("(Intercept)", "sexoMacho", "año", "altitud_num"))
  coefs$Especie <- especie
  coefs$Tipo <- tipo
  colnames(coefs)[1:2] <- c("Estimate", "StdError")
  colnames(coefs)[ncol(coefs)-3] <- "p_value"
  coefs$CI_low <- coefs$Estimate - 1.96 * coefs$StdError
  coefs$CI_high <- coefs$Estimate + 1.96 * coefs$StdError
  coefs$Estacion_Clean <- nombres_cronologicos[coefs$Variable]
  return(coefs)
}

coefs_mel_all <- bind_rows(
  extraer_coefs_mel(m_zr_temp_sex, "Zerynthia rumina", "Temperatura (°C)"),
  extraer_coefs_mel(m_zr_prec_sex, "Zerynthia rumina", "Precipitación (mm)"),
  extraer_coefs_mel(m_pa_temp_sex, "Parnassius apollo", "Temperatura (°C)"),
  extraer_coefs_mel(m_pa_prec_sex, "Parnassius apollo", "Precipitación (mm)")
)

write_csv(coefs_mel_all, "datos_procesados/res_melanismo_lmm_coefs.csv")

orden_niveles <- c(
  "5. Verano (Año 0)", "5. Primavera (Año 0)", "4. Primavera (Año 0)", "4. Invierno (Año 0)", 
  "3. Invierno (Año 0)", "3. Otoño (Año-1)", "2. Otoño (Año-1)", "2. Verano (Año-1)", 
  "1. Verano (Año-1)", "1. Primavera (Año-1)"
)

coefs_mel_all$Estacion_Clean <- factor(coefs_mel_all$Estacion_Clean, levels = unique(orden_niveles))

g_forest_mel <- ggplot(coefs_mel_all, aes(x = Estimate, y = Estacion_Clean, color = Especie)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_pointrange(aes(xmin = CI_low, xmax = CI_high), size = 0.7, linewidth = 0.9) +
  facet_wrap(~ Especie + Tipo, scales = "free", ncol = 2) +
  scale_color_manual(values = c("Zerynthia rumina" = "#E6A100", "Parnassius apollo" = "#4A505A")) +
  labs(x = "Efecto Estimado sobre Reflectancia Real (%)",
       y = "Secuencia Cronológica Estacional") +
  tema_tfm +
  theme(legend.position = "none")

ggsave("graficos/08_melanismo_clima/08_05_ForestPlot_Melanismo_LMM.png", plot = g_forest_mel, width = 11, height = 7, dpi = 900)

print("=== MÓDULO 08 EJECUTADO CON ÉXITO ===")


