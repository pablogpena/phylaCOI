#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recorre múltiples carpetas con archivos FASTA y ejecuta el análisis
automatizado definido en run_pipeline() desde procesar_fasta_unico.py.
"""

from pathlib import Path
from procesar_fasta_unico import run_pipeline  # Asegúrate de que el nombre del archivo es correcto

# Ruta al directorio raíz donde están las 60 carpetas
ROOT_DIR = Path("/mnt/e/Gonzalez-Miguens preprint 2025/000000_version2_maria/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS")  # <--- ¡MODIFICA ESTA RUTA!

# Nombres esperados de los archivos dentro de cada carpeta
NAMES_FILE = "nombres_secuencias.txt"
METADATA_FILE = "metadata_eKOI_ver2.csv"

# Recorremos cada carpeta
for folder in ROOT_DIR.iterdir():
    if folder.is_dir():
        fasta_files = list(folder.glob("*.fasta"))
        if not fasta_files:
            print(f"[SALTO] No se encontró archivo FASTA en {folder.name}")
            continue

        fasta_path = fasta_files[0]
        names_path = folder / NAMES_FILE
        meta_path = folder / METADATA_FILE
        out_dir = folder / "output"

        if names_path.exists() and meta_path.exists():
            print(f"[EJECUTANDO] Procesando carpeta: {folder.name}")
            try:
                run_pipeline(fasta_path, names_path, meta_path, out_dir)
            except Exception as e:
                print(f"[ERROR] Falló en {folder.name}: {e}")
        else:
            print(f"[OMITIDO] Archivos necesarios no encontrados en {folder.name}")
