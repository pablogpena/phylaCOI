# -*- coding: utf-8 -*-
import os
import subprocess
import pandas as pd

# Definir el nombre de la base de datos de referencia
base_de_datos = "eKOI_ver1.fasta"

# Definir el porcentaje mínimo de similitud aceptable
porcentaje_minimo = 84

# Obtener la lista de archivos FASTA en el directorio actual, excluyendo la base de datos de referencia
archivos_fasta = [archivo for archivo in os.listdir() if archivo.endswith('.fasta') and archivo != base_de_datos]

# Iterar sobre cada archivo FASTA
for archivo_fasta in archivos_fasta:
    # Crear una carpeta con el nombre del archivo FASTA
    nombre_carpeta = os.path.splitext(archivo_fasta)[0]
    os.makedirs(nombre_carpeta, exist_ok=True)

    # Nombre del archivo sin extensión
    nombre_archivo = os.path.splitext(archivo_fasta)[0]

    # Ejecutar vsearch para comparar el archivo FASTA con la base de datos
    comando = f"vsearch --usearch_global {archivo_fasta} --db {base_de_datos} --id 0.84 --blast6out {nombre_carpeta}/{nombre_archivo}_resultados.txt --quiet"
    subprocess.run(comando, shell=True)

    # Lista para almacenar los resultados de la búsqueda
    resultados = []

    # Leer los resultados del archivo de salida
    with open(f"{nombre_carpeta}/{nombre_archivo}_resultados.txt", 'r') as f:
        for linea in f:
            campos = linea.strip().split('\t')
            nombre_secuencia = campos[0]
            porcentaje_similitud = float(campos[2])
            secuencia_asignada = campos[1]
            resultados.append((nombre_secuencia, porcentaje_similitud, secuencia_asignada))

    # Crear un DataFrame con los resultados
    df = pd.DataFrame(resultados, columns=['Nombre de la Secuencia', 'Porcentaje de Similitud', 'Secuencia Asignada'])

    # Guardar los resultados en un archivo de Excel
    nombre_excel = f"{nombre_carpeta}/{nombre_archivo}_resultados.xlsx"
    df.to_excel(nombre_excel, index=False)

    # Filtrar las secuencias con un porcentaje de asignación mayor o igual al 84%
    df_filtrado = df[df['Porcentaje de Similitud'] >= porcentaje_minimo]

    # Separar la columna taxonómica por ";"
    df_filtrado['Secuencia Asignada'] = df_filtrado['Secuencia Asignada'].str.split(';')

    # Agrupar las secuencias por la quinta columna formada y guardarlas en archivos FASTA separados
    for group_name, group_df in df_filtrado.groupby(df_filtrado['Secuencia Asignada'].str[4]):
        # Crear un nuevo archivo FASTA para el grupo taxonómico
        fasta_filename = f"{nombre_carpeta}/{group_name}_{nombre_archivo}.fasta"
        with open(fasta_filename, 'w') as fasta_file:
            # Obtener los nombres de las secuencias del grupo
            secuencias_grupo = group_df['Nombre de la Secuencia'].tolist()
            # Escribir las secuencias correspondientes al grupo en el nuevo archivo FASTA
            with open(archivo_fasta, 'r') as original_fasta:
                for linea in original_fasta:
                    if linea.startswith('>'):
                        nombre_secuencia = linea.strip()[1:]
                        if nombre_secuencia in secuencias_grupo:
                            fasta_file.write(linea)
                            fasta_file.write(next(original_fasta))  # Escribir la secuencia

    # Guardar los resultados filtrados en un nuevo archivo de Excel
    nombre_excel_filtrado = f"{nombre_carpeta}/{nombre_archivo}_resultados_filtrados.xlsx"
    df_filtrado.to_excel(nombre_excel_filtrado, index=False)
