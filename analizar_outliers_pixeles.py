import pandas as pd
import os
import sys
from pathlib import Path
import traceback
import matplotlib.pyplot as plt

lepy_dir = r"C:\Users\Victor\Desktop\Victor_TFM_2026_07-07\LEPY-main"
sys.path.append(lepy_dir)
os.chdir(lepy_dir)

from lepy.image import Image
from lepy import utils
from lepy.outputs.writer.plotter import Plotter
from lepy.outputs.writer import OutputWriter

base_dir = r"C:\Users\Victor\Desktop\Victor_TFM_2026_07-07"
excel_path = r"C:\Users\Victor\Mi unidad (victoraa2002@gmail.com)\Bóveda Victólogo\03 Academia\Máster ZOO UCM\TFM\Miscelánea\ANEXOS\TFM_dataset_LEPY_FINAL.xlsx"

out_dir = os.path.join(base_dir, "output", "outliers_pixeles")
os.makedirs(out_dir, exist_ok=True)

print("Leyendo Excel...")
df_apollo = pd.read_excel(excel_path, sheet_name="P. apollo")
df_rumina = pd.read_excel(excel_path, sheet_name="Z. rumina")

# Filtrar outliers
outliers_apollo = df_apollo[
    (df_apollo["contour_width_calibrated"] > 85) |
    (df_apollo["contour_width_calibrated"] < 67)
]
outliers_rumina = df_rumina[
    (df_rumina["contour_width_calibrated"] > 50) |
    (df_rumina["contour_width_calibrated"] < 35)
]

targets = []
for _, row in outliers_apollo.iterrows():
    targets.append((row['Code'], "apollo", str(row['Coleccion'])))
for _, row in outliers_rumina.iterrows():
    targets.append((row['Code'], "rumina", str(row['Coleccion'])))

def get_config_file(especie, coleccion):
    col = coleccion.upper()
    if "MNCN" in col:
        return f"config_{especie}_MNCN.yml"
    elif "MONTES" in col:
        return f"config_{especie}_Montes_esc2.yml"
    else:
        return f"config_{especie}.yml"

def find_image(code):
    for root, dirs, files in os.walk(base_dir):
        # Saltar carpetas de output y LEPY-main
        dirs[:] = [d for d in dirs if d not in ('output', 'LEPY-main', 'herramientas_y_temporales')]
        for f in files:
            if code in f and f.lower().endswith(('.tif', '.tiff')):
                return Path(root), f
    return None, None

# Crear writer usando la API correcta de OutputWriter
writer = OutputWriter(out_dir, config="config_apollo.yml")
plotter = Plotter(out_dir, plot_interm=False, save_contours=False)

print(f"Total de outliers a analizar: {len(targets)}")

for code, especie, coleccion in targets:
    print(f"\n--- Procesando {code} ---")

    img_folder, img_name = find_image(code)
    if not img_name:
        print(f"  AVISO: imagen no encontrada para {code}")
        continue

    conf_file = get_config_file(especie, coleccion)
    conf_path = os.path.join(lepy_dir, conf_file)
    if not os.path.exists(conf_path):
        conf_path = os.path.join(lepy_dir, f"config_{especie}.yml")

    try:
        config = utils.read_config(conf_path)

        # Ajustar ruta de template al directorio de LEPY
        if hasattr(config, 'calibration') and hasattr(config.calibration, 'template_path'):
            tp = config.calibration.template_path
            if tp and not os.path.isabs(tp):
                config.calibration.template_path = os.path.join(lepy_dir, tp)

        config.segmentation.method = "flatbug"

        img = Image(root=img_folder, key=code, rgb_fname=img_name)
        img.read()
        img.segment(config.segmentation)

        # Calibración real (en mm) para el CSV
        calib = img.calibrate(config.calibration)
        real_scale = None
        if calib and calib.scale and calib.scale > 0:
            real_scale = calib.scale
        else:
            real_scale = float(config.calibration.get("fallback_scale", 33.0))

        # POIs en mm para el CSV
        pois_mm = img.pois(config.points_of_interest, scale=real_scale)
        img.color_stats()

        # Guardar CSV (en mm, compatible con el Excel maestro)
        writer(img_name, img.stats)

        # Plot con escala 1.0 (píxeles) solo para la imagen visual
        import copy
        calib_px = copy.copy(calib) if calib else None
        pois_px = img.pois(config.points_of_interest, scale=1.0)
        if calib_px is not None:
            calib_px.px_per_square = 1.0 * getattr(calib_px, 'size_per_square', 1.0)
        plotter.plot(img, pois=pois_px, calib_result=calib_px, show_scalebar=False)

        print(f"  OK: CSV en mm + plot en px guardados para {code}")
        plt.close('all')

    except Exception as e:
        print(f"  Error procesando {code}: {e}")
        traceback.print_exc()

print(f"\n=== COMPLETADO ===")
print(f"CSV de datos (en mm):   {out_dir}\\stats.csv")
print(f"Plots en píxeles:       {out_dir}\\visualisations\\")
