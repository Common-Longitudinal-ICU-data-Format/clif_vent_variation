# Source renv activation first to ensure renv is properly loaded
source("renv/activate.R")

# Get all library paths after renv has set up the environment
all_lib_paths <- .libPaths()

# Find the one that contains "ventilation_variation"
matching_path <- all_lib_paths[grepl("ventilation_variation", all_lib_paths)]

# If a matching path is found, set it as the only library path
if (length(matching_path) > 0) {
  .libPaths(matching_path)
}

# Print the final library path (for debugging)
print(.libPaths())
