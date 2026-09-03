#!/usr/bin/env Rscript

# Merge two sample-level CSV files by their shared sample_id column.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2 || length(args) > 3) {
  stop(
    "Usage: Rscript scripts/merge_sample_metadata.R \\\n  <species_calls.csv> <metadata.csv> [output.csv]",
    call. = FALSE
  )
}

species_calls_path <- args[[1]]
metadata_path <- args[[2]]
output_path <- if (length(args) == 3) {
  args[[3]]
} else {
  "merged_sample_data.csv"
}

input_paths <- c(species_calls = species_calls_path, metadata = metadata_path)
missing_files <- input_paths[!file.exists(input_paths)]
if (length(missing_files) > 0) {
  stop(
    "Input file(s) do not exist: ",
    paste(names(missing_files), missing_files, sep = " = ", collapse = "; "),
    call. = FALSE
  )
}

species_calls <- read.csv(species_calls_path, stringsAsFactors = FALSE, check.names = FALSE)
metadata <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)

if (!"sample_id" %in% names(species_calls) || !"sample_id" %in% names(metadata)) {
  stop("Both input CSVs must contain a column named 'sample_id'.", call. = FALSE)
}

for (input_name in names(input_paths)) {
  input_data <- if (input_name == "species_calls") species_calls else metadata
  if (anyNA(input_data$sample_id) || any(!nzchar(trimws(input_data$sample_id)))) {
    stop(input_name, " contains missing or empty sample_id values.", call. = FALSE)
  }
  if (anyDuplicated(input_data$sample_id)) {
    duplicate_ids <- unique(input_data$sample_id[duplicated(input_data$sample_id)])
    stop(
      input_name,
      " contains duplicate sample_id values: ",
      paste(head(duplicate_ids, 5), collapse = ", "),
      call. = FALSE
    )
  }
}

# Keep every sample from either input; unmatched fields are written as NA.
merged_data <- merge(
  species_calls,
  metadata,
  by = "sample_id",
  all = TRUE,
  sort = FALSE
)

output_directory <- dirname(output_path)
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}
write.csv(merged_data, output_path, row.names = FALSE, na = "")

message(
  "Wrote ", nrow(merged_data), " merged rows and ", ncol(merged_data),
  " columns to: ", normalizePath(output_path, mustWork = FALSE)
)
