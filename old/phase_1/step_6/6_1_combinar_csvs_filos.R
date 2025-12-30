# +
# ======================================================================
# Combine per-phylum div_abun_conn_combined.csv files into a master table
# ======================================================================

# Root directory containing one subfolder per phylum
root_folder <- "C:/TEMP/FASE1/3_magia_nombre_carpeta/eKOI_metabarcoding_database_FILOS"

# (Optional) List all matching files recursively, for sanity checking
# This line does not affect the workflow
list.files(root_folder, recursive = TRUE, pattern = "div_abun_conn_combined.csv")

# List first-level subdirectories (one per phylum; non-recursive)
phyla_folders <- list.dirs(root_folder, recursive = FALSE)

# Initialize an empty list to store data frames from each phylum
all_data <- list()

# Loop over each phylum folder
for (folder in phyla_folders) {
  
  # Path to the combined per-phylum file inside the "output" subdirectory
  file_path <- file.path(folder, "output", "div_abun_conn_combined.csv")
  
  # Check whether the file exists
  if (file.exists(file_path)) {
    
    # Read the per-phylum combined dataset
    df <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Add a column identifying the phylum (folder name)
    df$Phylum <- basename(folder)
    
    # Store the data frame in the list
    all_data[[length(all_data) + 1]] <- df
    
  } else {
    # Inform the user if the expected file is missing
    message("File not found: ", file_path)
  }
}

# Combine all per-phylum data frames into a single master table
combined_df <- bind_rows(all_data)

# Write the master dataset to the root directory
write.csv(
  combined_df,
  file = file.path(root_folder, "div_abun_conn_combined_master.csv"),
  row.names = FALSE
)
