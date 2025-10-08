#!/usr/bin/env Rscript
# Automate haplotype assignment across multiple phylum folders

suppressPackageStartupMessages({
  library(Biostrings)  # readDNAStringSet
  library(ape)         # as.DNAbin
  library(pegas)       # haplotype
})

# -----------------------------------------------------------------------------#
# Paths
# -----------------------------------------------------------------------------#
base_dir <- "E:/Gonzalez-Miguens preprint 2025/000000_version2_maria/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS"

# Obtener subcarpetas de filos
phyla_dirs <- list.dirs(base_dir, recursive = FALSE)

# Función principal
process_phylum <- function(phylum_path) {
  phylum_name <- basename(phylum_path)
  message("Procesando: ", phylum_name)
  
  # Rutas a los archivos necesarios
  abundances_file    <- file.path(phylum_path, "output", "abundances_unique.csv")
  aligned_fasta_file <- file.path(phylum_path, "output", "aligned_sequences_mafft.fasta")
  otus_map_file      <- file.path(phylum_path, "output", "otus", "otus_mapping.txt")
  otus_keep_file     <- file.path(phylum_path, "output", "otus", "informative_OTU.txt")
  output_file        <- file.path(phylum_path, "output", "abundances_filtered_haplotypes.csv")
  
  # Verificación de existencia de archivos
  required_files <- c(abundances_file, aligned_fasta_file, otus_map_file, otus_keep_file)
  if (!all(file.exists(required_files))) {
    warning("Faltan archivos en: ", phylum_name)
    return(NULL)
  }
  
  # -----------------------------------------------------------------------------#
  # Load data
  # -----------------------------------------------------------------------------#
  abundances <- read.csv(abundances_file, stringsAsFactors = FALSE)
  otus_map   <- read.table(otus_map_file,  col.names = c("OTU", "UniqueID"))
  otus_keep  <- scan(otus_keep_file, what = character())
  
  # Añadir OTU a la tabla de abundancia
  abundances$OTU <- otus_map$OTU[match(abundances$UniqueID, otus_map$UniqueID)]
  
  # Filtrar OTUs informativas
  abundances_filt <- abundances[abundances$OTU %in% otus_keep, ]
  if (nrow(abundances_filt) == 0) {
    message("No hay datos que procesar en: ", phylum_name)
    return(NULL)
  }

  # -----------------------------------------------------------------------------#
  # Read alignment and convert to DNAbin
  # -----------------------------------------------------------------------------#
  alignment <- readDNAStringSet(aligned_fasta_file)
  alignment_dnabin <- as.DNAbin(alignment)
  
  # -----------------------------------------------------------------------------#
  # Assign haplotypes within each OTU
  # -----------------------------------------------------------------------------#
  # Inicializar columna de haplotipos
  abundances_filt$Haplotype <- NA
  
  # Asignar haplotipos por OTU
  for (otu in unique(abundances_filt$OTU)) {
    ids <- abundances_filt$UniqueID[abundances_filt$OTU == otu]
    aln_otu <- alignment_dnabin[names(alignment_dnabin) %in% ids]
    
    if (length(aln_otu) > 1) {
      haps <- haplotype(aln_otu)
      idx  <- attr(haps, "index")
      seq_nm <- names(aln_otu)
      
      for (h in seq_along(idx)) {
        abundances_filt$Haplotype[abundances_filt$UniqueID %in% seq_nm[idx[[h]]]] <- h
      }
    }
  }
  
  
  # -----------------------------------------------------------------------------#
  # Export result
  # -----------------------------------------------------------------------------#
  write.csv(abundances_filt, output_file, row.names = FALSE)
  message("Resultado guardado en: ", output_file)
}

# Procesar cada filo
for (phylum_path in phyla_dirs) {
  tryCatch({
    process_phylum(phylum_path)
  }, error = function(e) {
    message("Error en ", phylum_path, ": ", e$message)
  })
}
