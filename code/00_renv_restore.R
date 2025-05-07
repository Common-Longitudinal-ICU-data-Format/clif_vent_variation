# 00_renv_restore.R
# Ensure renv is installed and load the project environment

################################################ OLD BELOW
# if (!requireNamespace("renv", quietly = TRUE)) {
#   install.packages("renv")
# }
# 
# 
# ##### TROUBLESHOOTING START #####
# ##### Try running below and then RERUN this file 

# renv::deactivate()
# unlink("renv/library", recursive = TRUE, force = TRUE)  # Deletes the renv library
# renv::activate()

# ##### TROUBLESHOOTING  END  #####
# 
# 
# # Restore the project's package environment
# renv::restore(prompt = FALSE, rebuild = TRUE)
# renv::activate()
################################################ OLD ABOVE



###############################################################
###############################################################

# 3/5/2025 Trying to fix things after renv broke again
# Things I THINK worked
#           1) Updated to the newest RENV
#           2) deleted RENV folder (KEEP THE renv.lock)
#                  --then manually install with renv::install(c"") (see below)
#           3) My R was diff from renv.lock file and it MIGHT have been causing issues
#                  --ultimately i don't think it was as big of an issue
#           4) the setting below might have helped at some point but ultimately might not have mattered
#                  --commenting out for now
#             ~~~~~~~~~~~~~~~~~~~ could run the code which is great, but still having issues with status not syncing
#           5) might be an issue with 2 diff renv places.  Added to Rprofile to state we need to use renv in THIS project
#           6) might be an issue with which packages need sourcing.  see code below to help
#           
#           
#           
#           
# Ensure renv is installed
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}



# # Disable cache to force full package installations
# renv::settings$use.cache(FALSE)
# 
# # Ensure we install everything from binary
# options(pkgType = "binary")
# options(renv.config.ppm.enabled = FALSE)
# options(renv.config.install.staged = FALSE)
# 
# # Set renv to strict snapshot mode (ensures exact package versions)
# renv::settings$snapshot.type("explicit")
# 
# # Disable Posit Package Manager to prevent automatic transformations
# renv::settings$ppm.enabled(FALSE)


# Capture renv status to detect issues
# status_info <- capture.output(renv::status())
# 
# # If renv is inconsistent, perform cleanup & rebuild
# if (any(grepl("inconsistent|missing", status_info, ignore.case = TRUE))) {
#   message("\n Issues detected in renv environment! Cleaning and rebuilding...\n")
#   
#   # Step 1: Remove stale locks and clean the renv library
#   # renv::clean()
#   
#   # Step 2: Remove any cached versions of packages
#   renv::purge("", prompt = FALSE)
#   
#   # Step 3: Fully restore and force rebuild all packages
#   renv::restore(prompt = FALSE, rebuild = TRUE)
#   
#   message("\n renv environment successfully rebuilt!\n")
# } else {
#   message("\n renv environment is already consistent. Skipping restore.\n")
# }
# 
# # Final activation to ensure renv is fully loaded
# # renv::activate()
# 
# # Troubleshooting instructions (if renv keeps failing)
# message("\n If renv still fails, try running the following manually:\n")
# message("renv::deactivate()")
# message('unlink("renv/library", recursive = TRUE, force = TRUE)')
# message("renv::activate()")
# 
# 
###############################################################################
#   How to find Packages that need SOURCE

installed_pkgs <- installed.packages()

# Convert R version to numeric correctly
r_version <- as.numeric(paste(R.version$major, R.version$minor, sep = "."))

# Extract 'Built' column safely, handling missing values
built_version <- installed_pkgs[, "Built"]
built_version <- suppressWarnings(as.numeric(sub("-.*", "", built_version)))  # Extract major.minor from 'Built'

# Handle missing values in Built column
built_version[is.na(built_version)] <- 0  # Replace NA with 0 to force mismatch

# Identify packages built under a different R version
needs_source_version <- rownames(installed_pkgs)[built_version != r_version]

# Get packages that require compilation (C/C++/Fortran)
needs_compilation <- rownames(installed_pkgs)[installed_pkgs[, "NeedsCompilation"] == "yes"]

# Exclude base R packages (which should not be rebuilt)
base_packages <- rownames(installed_pkgs)[installed_pkgs[, "Priority"] %in% c("base", "recommended")]

# Remove corrupted packages (those missing `package.rds` in Meta/)
corrupted_pkgs <- rownames(installed_pkgs)[!file.exists(file.path(installed_pkgs[, "LibPath"], installed_pkgs[, "Package"], "Meta/package.rds"))]

# Combine all checks and remove NAs
needs_source <- setdiff(unique(c(needs_source_version, needs_compilation, corrupted_pkgs)), base_packages)
needs_source <- needs_source[!is.na(needs_source)]  # Remove any NA values

# Print the list of packages that should be installed from source
print(needs_source)




##### 
#     only run below if you are having issues when you KNOW you have package installed 
#     but its not loading in the right way for some reason

# renv::remove(needs_source)
# renv::install(c("glue", "magrittr", "rlang", "vctrs", "knitr"))  ## Key dependencies
# renv::rebuild(needs_source, type = "source", prompt = FALSE)

###############################################################################
###############################################################################
###############################################################################
# #   If you need a full reboot, this should work as long as lock.file is there!
# 
# # Completely remove the renv cache
# renv::deactivate()
# unlink("Y:/DataStageData/Nick/CLIF/CLIF_ventilation_variation_2024/renv", recursive = TRUE, force = TRUE)
# 
# # Manually clean staging & remove broken libraries
# unlink("Y:/DataStageData/Nick/CLIF/CLIF_ventilation_variation_2024/renv/staging", recursive = TRUE, force = TRUE)
# unlink("Y:/DataStageData/Nick/CLIF/CLIF_ventilation_variation_2024/renv/library", recursive = TRUE, force = TRUE)
# 
# # Clear Renv cache
# Sys.setenv(RENV_PATHS_CACHE = "")
# 
# renv::activate()
# renv::clean()
# 
# ## These are essential, get these installed first
# renv::install(c("glue", "magrittr", "rlang", "vctrs", "knitr"))
# 
# renv::restore()

###############################################################################
###############################################################################
###############################################################################

# # Initialize renv for the project:
# renv::init(bare = TRUE, settings = list(use.cache = FALSE))
# 
# # Install required packages:
# renv::install("BiocManager")
# BiocManager::install("IRanges")

# # renv::install(c("tidyverse", "ggthemes", "systemfonts",  "styler", "readxl", "writexl", "DBI", "dbplyr", "knitr", "pandoc", "janitor", "data.table", "duckdb" ,"powerjoin", "collapse", "tidyfast", "datapasta", "fst", "dtplyr", "bit64", "zoo", "tzdb", "fuzzyjoin", "arrow", "hrbrthemes", "here", "table1", "rvest", "tidymodels", "pscl", "survminer", "gt", "gtsummary", "broom.helpers", "broom.mixed", "lme4", "lmerTest", "merTools", "finalfit"))
# 
#   
# renv::install(c("tidyverse", "systemfonts", "readxl", "DBI", "dbplyr", "knitr", "janitor", "data.table", "powerjoin", "collapse", "tidyfast", "datapasta", "fst", "dtplyr", "bit64", "tzdb", "arrow", "hrbrthemes", "here", "table1", "rvest", "tidymodels", "pscl", "survival", "survminer", "gt", "gtsummary", "broom.helpers", "finalfit"))
# 
# # # Save the project's package state:
# renv::snapshot()
################################################################################

library(lmerTest)
library(tzdb)
library(tidyverse)
library(systemfonts)
library(readxl)
library(DBI)
library(dbplyr)
library(knitr)
library(janitor)
library(data.table)
library(powerjoin)
library(collapse)
library(tidyfast)
library(datapasta)
library(fst)
library(dtplyr)
library(bit64)
library(arrow)
library(hrbrthemes)
library(here)
library(table1)
library(rvest)
library(gt)
library(gtsummary)
library(broom.helpers)
library(broom.mixed)
library(tidymodels)
library(pscl)
library(survival)
library(survminer)
# library(merTools)
library(lme4)
library(finalfit)
library(furrr)
library(rms)
library(ggeffects)
library(rspiro)


options(arrow.unsafe_metadata = TRUE)
select <- dplyr::select

setDTthreads(0)

# Plan for parallel execution
try(plan(multisession), silent = TRUE)



## create output directories 
create_directories <- function() {
  # Define the directories
  dirs <- c("output", "output/intermediate", "output/final",  "output/intermediate/clean_db")
  
  # Loop through and create directories if they don't exist
  for (dir in dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      message(paste("Directory created:", dir))
    } else {
      message(paste("Directory already exists:", dir))
    }
  }
}

# Call the function to create the directories
create_directories()



# Load the configuration utility
source("utils/config.R")

# Access configuration parameters
site_name <- config$site_name
tables_path <- config$tables_path
file_type <- config$file_type
time_zone <- config$time_zone


valid_time_zones <- c('US/Alaska', 'US/Aleutian', 'US/Arizona', 'US/Central', 'US/East-Indiana', 'US/Eastern', 'US/Hawaii', 'US/Indiana-Starke', 'US/Michigan', 'US/Mountain', 'US/Pacific', 'US/Samoa')

if (is.null(time_zone) || !(time_zone %in% valid_time_zones)) {
  stop(paste0("Please update the config script to specify a valid time_zone. ",
              "Valid options are: ", paste(valid_time_zones, collapse = ", "), "."))
} else {
  message(paste("Using time zone:", time_zone))
}




# Print the configuration parameters
print(paste("Site Name:", site_name))
print(paste("Tables Path:", tables_path))
print(paste("File Type:", file_type))
print(paste("Time Zone:", time_zone))



# List of all table names from the CLIF 2.0 ERD
tables <- c(
  "patient",
  "hospitalization",
  "vitals",
  "labs",
  "medication_admin_continuous",
  "adt",
  "patient_assessments",
  "respiratory_support",
  "position",
  "dialysis",
  "intake_output",
  "ecmo_mcs",
  "procedures",
  "admission_diagnosis",
  "provider",
  "sensitivity",
  "medication_orders",
  "medication_admin_intermittent",
  "therapy_details",
  "microbiology_culture",
  "sensitivity",
  "microbiology_nonculture"
)

# Tables that should be set to TRUE for this project
true_tables <- c(
  "hospitalization",
  "vitals",
  "labs",
  # "medication_admin_continuous",
  "adt",
  # "patient_assessments",
  "respiratory_support",
  # "position",
  # "dialysis",
  # "intake_output",
  # "ecmo_mcs",
  # "procedures",
  # "admission_diagnosis",
  # "provider",
  # "sensitivity",
  # "medication_orders",
  # "medication_admin_intermittent",
  # "therapy_details",
  # "microbiology_culture",
  # "sensitivity",
  # "microbiology_nonculture",
  "patient"
)

# Create a named vector and set the boolean values
table_flags <- setNames(tables %in% true_tables, tables)


# List all CLIF files in the directory
clif_table_filenames <- list.files(path = tables_path, 
                                   pattern = paste0("^clif_.*\\.", file_type, "$"), 
                                   full.names = TRUE)

# Extract the base names of the files (without extension)
clif_table_basenames <- basename(clif_table_filenames) %>%
  str_remove(paste0("\\.", file_type, "$"))

# Create a lookup table for required files based on table_flags
required_files <- paste0("clif_", names(table_flags)[table_flags])

# Check if all required files are present
missing_tables <- setdiff(required_files, clif_table_basenames)
if (length(missing_tables) > 0) {
  stop(paste("Missing required tables:", 
             paste(missing_tables, collapse = ", ")))
}

# Filter only the filenames that are required
required_filenames <- clif_table_filenames[clif_table_basenames %in% required_files] 


# Quickly look at the data like we would in stata, default is 100
ni_peek <- function(x, n=100){
  view(head(x, n))
}

# check missing variables
ni_count_missing <- function(df, group_vars, x_var) {
  df |> 
    group_by(pick({{ group_vars }})) |> 
    summarise(
      n_miss = sum(is.na({{ x_var }})),
      .groups = "drop"
    )
}
# flights |>  count_missing(c(year, month, day), dep_time)

# look at the variables
ni_check_variables <- function(df, n=500) {
  check <- tibble(
    col_name = names(df), 
    col_type = map_chr(df, vctrs::vec_ptype_full),
    n_miss = map_int(df, \(x) sum(is.na(x)))
  )
  print(check, n=n)
}

ni_count_prop <- function(df, var, sort = TRUE) {
  df |>
    count({{ var }}, sort = sort) |>
    mutate(prop = n / sum(n)
    )
}

# Identify and create a table of duplicate rows based on clif_hospitalizations_joined_id
ni_duplicate_finder <- function(df, group_vars = c(clif_hospitalizations_joined_id), n=1){
  df |> 
    dplyr::group_by(pick({{ group_vars }})) |> 
    filter(n() > {{ n }}) |> 
    ungroup()
}

ni_tic <- function() {
  .GlobalEnv$time_start_temp <- proc.time()
  .GlobalEnv$time_start_sys <- Sys.time()
  return(Sys.time())
}

ni_toc <- function() {
  time_end_sys <- Sys.time()
  .GlobalEnv$time_diff.time <- round(time_end_sys - time_start_sys,2)
  mylist <- list(time_diff.time, cat("Finished in",timetaken(time_start_temp),"\n"
  ))
  return(print(mylist[[1]]))
}

fio2warning <- function() {
  warning("fio2 variable needed to be fixed, it was multiplied by 100!!!")
}

labwarning <- function() {
  warning("lab values are character and needed to be fixed, they were forced to numeric!!!")
}

vitalwarning <- function() {
  warning("vital values are character and needed to be fixed, they were forced to numeric!!!")
}



read_data <- function(file_path) {
  file_path_temp <- paste0(tables_path,"/", file_path, ".", file_type)
  local_tz = time_zone
  print(file_path_temp)
  print(local_tz)
  
  if (grepl("\\.csv$", file_path_temp)) {
    return(read.csv(file_path_temp) |> 
             mutate(across(ends_with("_id"), as.character)) |> 
             mutate(across(where(~ inherits(.x, c("POSIXct", "POSIXlt"))),
                           ~ with_tz(.x, local_tz))) |> 
             ## make sure DATES are DATES
             mutate(across(ends_with("_date"), ~ as.Date(.x))) |> 
             mutate(across(starts_with("date_"), ~ as.Date(.x))) 
    )
    
  } else if (grepl("\\.parquet$", file_path_temp)) {
    return(arrow::open_dataset(file_path_temp)  |> 
             mutate(across(ends_with("_id"), as.character)) |> 
             mutate(across(where(~ inherits(.x, c("POSIXct", "POSIXlt"))),
                           ~ with_tz(.x, local_tz))) |>
             ## make sure DATES are DATES
             mutate(across(ends_with("_date"), ~ as.Date(.x))) |> 
             mutate(across(starts_with("date_"), ~ as.Date(.x))) |> 
             collect()
    )
    
  } else if (grepl("\\.fst$", file_path_temp)) {
    return(fst::read.fst(file_path_temp) |> 
             mutate(across(ends_with("_id"), as.character)) |> 
             mutate(across(where(~ inherits(.x, c("POSIXct", "POSIXlt"))),
                           ~ with_tz(.x, local_tz))) |> 
             ## make sure DATES are DATES
             mutate(across(ends_with("_date"), ~ as.Date(.x))) |> 
             mutate(across(starts_with("date_"), ~ as.Date(.x))) 
    )
  } else {
    stop("Unsupported file format")
  }
}



# for quick opening and filtering...
#     ni_open_dataset_clif(adt) |> filter(patient_id == "asdf") |> collect() |> View()
ni_open_dataset_clif <- function(file){
  file_quoted <- deparse(substitute(file))
  open_dataset(paste0(tables_path, "/",  "clif_", file_quoted, ".parquet"), thrift_string_size_limit = 1000000000) |> 
    mutate(across(where(is.character), as.character)) |> 
    mutate(across(where(is.character), str_to_lower)) |> 
    mutate(across(where(~ inherits(.x, c("POSIXct", "POSIXlt"))),
                  ~ with_tz(.x, local_tz))) |> 
    ## make sure DATES are DATES
    mutate(across(ends_with("_date"), ~ as.Date(.x))) |> 
    mutate(across(starts_with("date_"), ~ as.Date(.x))) 
}



### Make sure to save in UTC
write_parquet_tz <- function(df, file_path) {
  # Identify datetime columns
  datetime_cols <- names(df)[sapply(df, function(col) inherits(col, c("POSIXct", "POSIXlt")))]
  
  if (length(datetime_cols) > 0) {
    df <- df %>%
      mutate(across(all_of(datetime_cols), ~ with_tz(., "UTC")))  # Attach UTC correctly
    
    # df <- as_arrow_table(df) 
    
    message("Stored timestamps correctly in UTC for columns: ", 
            paste(datetime_cols, collapse = ", "))
  } else {
    # df <- as_arrow_table(df)
    message("No datetime columns required conversion.")
  }
  
  # Write to Parquet
  arrow::write_parquet(df, file_path)
  
  message("Saved ", file_path, " successfully.")
}



## Date range
start_date <- "2018-01-01"
end_date <- "2024-09-30"




# Function to source a .qmd file
source_qmd <- function(file, local = FALSE, ...){
  options(knitr.duplicate.label = 'allow')
  
  tempR <- tempfile(tmpdir = ".", fileext = ".R")
  on.exit(unlink(tempR))
  knitr::purl(file, output=tempR)
  
  envir <- globalenv()
  source(tempR, local = envir, ...)
}


message("##########################################
You have completed setup!!! Woohooo!!! \n Now proceeding on to 01_cohorting script \n 
##########################################")


# Example usage: 
# Assuming "my_data_analysis.qmd" is the .qmd file you want to source
# source_qmd(here("code", "01_cohorting.qmd"), echo = TRUE)

# 
# 
# message("##########################################
# You have completed Cohorting!!! Woohooo!!! \n Now proceeding on to 02_statistical_analysis script \n 
# ##########################################")


# source_qmd(here("code", "02_statistical_analysis.qmd"), echo = TRUE)



output_dir_hosp = "output/intermediate/clean_db"
df_name_hosp <- "hospital_order"
date_file_hosp <- file.path(output_dir_hosp, "most_recent_save_resp.Rdata")
load(date_file_hosp)  # This will load 'most_recent_save_resp' variable
input_file_hosp <- file.path(output_dir_hosp, paste0(df_name_hosp, "_.", most_recent_save_resp, ".parquet"))
df_hosp <- arrow::read_parquet(input_file_hosp)
df_hosp |> View()
write.csv(df_hosp, paste0("output/final/hospital_order_table__",site_name,"_", Sys.Date(),".csv"), row.names = FALSE, fileEncoding = "UTF-8")



  