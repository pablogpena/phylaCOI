from Bio import SeqIO

# Lee el archivo fasta
input_file = "eKOI_metabarcoding_database_ver2.fasta"
output_file = "eKOI_metabarcoding_database_clean.fasta"

with open(output_file, "w") as output_handle:
    for record in SeqIO.parse(input_file, "fasta"):
        # Divide el identificador de la secuencia por el espacio y toma solo la primera parte
        identifier = record.id.split(" ")[0]
        # Escribe la nueva secuencia en el nuevo archivo fasta
        output_handle.write(">" + identifier + "\n")
        output_handle.write(str(record.seq) + "\n")

print("Se ha creado el nuevo archivo fasta:", output_file)
