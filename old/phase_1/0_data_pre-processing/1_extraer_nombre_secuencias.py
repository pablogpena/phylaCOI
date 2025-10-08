# -*- coding: utf-8 -*-
# Abre el archivo fasta en modo lectura
with open("eKOI_metabarcoding_database_ver2.fasta", "r") as fasta_file:
    # Abre un nuevo archivo de texto en modo escritura
    with open("nombres_secuencias.txt", "w") as txt_file:
        # Lee el archivo fasta línea por línea
        for line in fasta_file:
            # Si la línea comienza con ">" (indicando el nombre de la secuencia)
            if line.startswith(">"):
                # Escribe el nombre de la secuencia en el archivo de texto
                txt_file.write(line.strip()[1:] + "\n")
