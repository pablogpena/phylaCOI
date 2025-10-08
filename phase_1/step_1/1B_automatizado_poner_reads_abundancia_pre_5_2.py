import pandas as pd
import os
import re

# Ruta del directorio raíz con las carpetas de los filos
base_dir = r"E:\Gonzalez-Miguens preprint 2025\000000_version2_maria\FASE1\3_magia_nombre_carpeta\eKOI_metabarcoding_database_FILOS"

# Iterar sobre cada carpeta de filo
for phylum_folder in os.listdir(base_dir):
    phylum_path = os.path.join(base_dir, phylum_folder)
    
    if not os.path.isdir(phylum_path):
        continue  # Saltar si no es una carpeta

    print(f"\nProcesando filo: {phylum_folder}")

    output_dir = os.path.join(phylum_path, "output")
    nombres_secuencias_path = os.path.join(output_dir, "nombres_secuencias.txt")
    abundances_unique_path = os.path.join(output_dir, "abundances_unique.csv")

    # Verificar que existan ambos archivos
    if not os.path.isfile(nombres_secuencias_path) or not os.path.isfile(abundances_unique_path):
        print(f"  ❌ Archivos faltantes en: {output_dir}")
        continue

    # Cargar nombres_secuencias.txt
    secuencia_abundancias = {}
    with open(nombres_secuencias_path, "r", encoding="utf-8") as f:
        for line in f:
            match = re.match(r"(\S+)\s+merged_sample=\{(.*)\}", line)
            if match:
                secuencia, localidades = match.groups()
                localidades_dict = re.findall(r"'([^']+)'\s*:\s*(\d+)", localidades)
                secuencia_abundancias[secuencia] = {loc.lower(): int(reads) for loc, reads in localidades_dict}

    # Cargar abundances_unique.csv
    abundances_unique = pd.read_csv(abundances_unique_path)

    # Añadir columna "Abundancia"
    abundancias = []
    for _, row in abundances_unique.iterrows():
        seq_id = row["ID"]
        locality = str(row["Localities"]).lower()
        abundancia = secuencia_abundancias.get(seq_id, {}).get(locality)
        abundancias.append(abundancia)

    abundances_unique["Abundancia"] = abundancias

    # Guardar archivo actualizado
    output_file = os.path.join(output_dir, "abundances_unique_actualizado.csv")
    abundances_unique.to_csv(output_file, index=False)
    print(f"  ✅ Archivo guardado: {output_file}")

