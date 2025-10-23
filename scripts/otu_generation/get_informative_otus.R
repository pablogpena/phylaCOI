# +
# #!/usr/bin/env Rscript
# -

suppressPackageStartupMessages({
  library(dplyr)
})

# +
# Process the input directory flag 
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2 || !(args[1] %in% c("-i", "--input"))) {
  stop("Usage: Rscript run_distance_pipeline.R -i <phylums_directory>")
    base_dir <- normalizePath(args[2], mustWork = TRUE)
message("Base phylum directory: ", base_dir)
}
# -

#Load the funcions 
source("/workspace/src/otu_utils/informative_otus_utils.R")

# +
#Get all the phyllums in the working directory
phyla_dirs <- list.dirs(base_dir, recursive = FALSE)

if (length(phyla_dirs) == 0) {
  stop("No phylum folders were found in the specified base directory: ", base_dir)
}

# +
#Process each phylums otus 
for (phylum_path in phyla_dirs) {
  tryCatch({
    process_phylum(phylum_path)
  }, error = function(e) {
    message("Error processing ", phylum_path, ": ", e$message)
  })
}

message("Pipeline completed successfully for all phyla.")
