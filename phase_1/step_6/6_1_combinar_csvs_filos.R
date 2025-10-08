# Carpeta raíz donde están las carpetas de los filos
root_folder <- "C:/TEMP/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS"

list.files(folder, recursive = TRUE, pattern = "div_abun_conn_combined.csv")




# Listamos las carpetas de cada filo dentro de root_folder (no recursivo)
phyla_folders <- list.dirs(root_folder, recursive = FALSE)

# Lista para guardar los data.frames de cada filo
all_data <- list()

for (folder in phyla_folders) {
  # Ruta al archivo div_abun_conn_combined.csv dentro de la subcarpeta "output"
  file_path <- file.path(folder, "output", "div_abun_conn_combined.csv")
  
  if (file.exists(file_path)) {
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Añadir columna con el nombre del filo (nombre de la carpeta)
    df$Filo <- basename(folder)
    
    all_data[[length(all_data) + 1]] <- df
  } else {
    message("Archivo no encontrado en: ", file_path)
  }
}

# Combinar todos los data.frames
combined_df <- bind_rows(all_data)

# Guardar el archivo combinado en la carpeta raíz
write.csv(combined_df, file = file.path(root_folder, "div_abun_conn_combined_master.csv"), row.names = FALSE)
