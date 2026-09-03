#!/usr/bin/env Rscript

# Generate one map of sampling locations for each aim_species value.

required_packages <- "maps"
missing_packages <- required_packages[!vapply(
  required_packages,
  requireNamespace,
  logical(1),
  quietly = TRUE
)]
if (length(missing_packages) > 0) {
  stop(
    "Install the required package(s) before running this script: ",
    paste(missing_packages, collapse = ", "),
    "\nFor example: install.packages(c(\"",
    paste(missing_packages, collapse = "\", \""),
    "\"))",
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
script_path <- commandArgs()[grepl("^--file=", commandArgs())]
script_dir <- if (length(script_path) == 1) {
  dirname(normalizePath(sub("^--file=", "", script_path)))
} else {
  getwd()
}
project_root <- normalizePath(file.path(script_dir, ".."), mustWork = TRUE)

metadata_path <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(project_root, "data", "vo_agam_release", "v3", "metadata", "merged_sample_data.csv")
}
output_directory <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(project_root, "sampling_location_maps")
}

if (!file.exists(metadata_path)) {
  stop("Metadata file does not exist: ", metadata_path, call. = FALSE)
}

samples <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
coordinate_columns <- c("latitude", "longitude")
required_columns <- c(coordinate_columns, "aim_species", "year", "month")
if (!all(required_columns %in% names(samples))) {
  stop(
    "The input must contain columns named ",
    paste(required_columns, collapse = ", "),
    call. = FALSE
  )
}

samples$latitude <- suppressWarnings(as.numeric(samples$latitude))
samples$longitude <- suppressWarnings(as.numeric(samples$longitude))
samples$year <- suppressWarnings(as.numeric(samples$year))
samples$month <- suppressWarnings(as.numeric(samples$month))
valid_coordinates <- with(
  samples,
  is.finite(latitude) & is.finite(longitude) &
    latitude >= -90 & latitude <= 90 &
    longitude >= -180 & longitude <= 180
)
if (!all(valid_coordinates)) {
  warning(sum(!valid_coordinates), " rows have invalid coordinates and will be omitted.")
  samples <- samples[valid_coordinates, , drop = FALSE]
}
if (nrow(samples) == 0) {
  stop("No valid sampling coordinates were found.", call. = FALSE)
}
samples$aim_species <- trimws(samples$aim_species)
if (any(!nzchar(samples$aim_species) | is.na(samples$aim_species))) {
  warning("Rows without an aim_species value will be omitted.")
  samples <- samples[!is.na(samples$aim_species) & nzchar(samples$aim_species), , drop = FALSE]
}
if (nrow(samples) == 0) {
  stop("No samples with valid coordinates and aim_species values were found.", call. = FALSE)
}

plot_species <- c("arabiensis", "coluzzii", "gambiae")
ignored_samples <- !samples$aim_species %in% plot_species
if (any(ignored_samples)) {
  warning(
    sum(ignored_samples),
    " intermediate or unsupported samples will be omitted from the plots."
  )
  samples <- samples[!ignored_samples, , drop = FALSE]
}
if (nrow(samples) == 0) {
  stop("No samples remain after excluding intermediate species.", call. = FALSE)
}

valid_time <- is.finite(samples$year) & samples$year > 0 &
  is.finite(samples$month) & samples$month >= 1 & samples$month <= 12
if (any(!valid_time)) {
  warning(sum(!valid_time), " samples have unknown or invalid year/month values and will be omitted.")
  samples <- samples[valid_time, , drop = FALSE]
}
if (nrow(samples) == 0) {
  stop("No samples remain after excluding unknown or invalid year/month values.", call. = FALSE)
}

if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

species_values <- plot_species[plot_species %in% unique(samples$aim_species)]
for (species in species_values) {
  species_samples <- samples[samples$aim_species == species, , drop = FALSE]

  longitude_range <- range(species_samples$longitude)
  latitude_range <- range(species_samples$latitude)
  longitude_padding <- max(diff(longitude_range) * 0.05, 2)
  latitude_padding <- max(diff(latitude_range) * 0.05, 2)
  x_limits <- c(longitude_range[1] - longitude_padding, longitude_range[2] + longitude_padding)
  y_limits <- c(latitude_range[1] - latitude_padding, latitude_range[2] + latitude_padding)
  safe_species_name <- gsub("[^A-Za-z0-9]+", "_", tolower(species))
  output_path <- file.path(output_directory, paste0("sampling_locations_", safe_species_name, ".png"))
  years <- sort(unique(species_samples$year))
  panel_columns <- min(4, max(1, ceiling(sqrt(length(years)))))
  panel_rows <- ceiling(length(years) / panel_columns)
  month_colors <- grDevices::hcl.colors(12, palette = "Dynamic")

  png(output_path, width = 600 * panel_columns, height = 520 * panel_rows, res = 150)
  par(
    mfrow = c(panel_rows, panel_columns),
    mar = c(2.5, 2.5, 2.5, 1),
    oma = c(3, 0, 3, 0)
  )
  for (year in years) {
    year_samples <- species_samples[species_samples$year == year, , drop = FALSE]
    # Keep separate month-colored points when a site was sampled in multiple months.
    site_key <- paste(year_samples$latitude, year_samples$longitude, year_samples$month, sep = "_")
    site_counts <- as.data.frame(table(site_key), stringsAsFactors = FALSE)
    site_coordinates <- unique(year_samples[c("latitude", "longitude", "month")])
    site_coordinates$site_key <- paste(site_coordinates$latitude, site_coordinates$longitude, site_coordinates$month, sep = "_")
    sites <- merge(site_coordinates, site_counts, by = "site_key", sort = FALSE)
    names(sites)[names(sites) == "Freq"] <- "sample_count"

    maps::map(
      "world",
      xlim = x_limits,
      ylim = y_limits,
      fill = TRUE,
      col = "white",
      border = "grey55",
      bg = "aliceblue"
    )
    points(
      sites$longitude,
      sites$latitude,
      pch = 21,
      bg = month_colors[sites$month],
      col = "white",
      lwd = 0.7,
      cex = 0.8 + 1.8 * sqrt(sites$sample_count / max(sites$sample_count))
    )
    title(main = as.character(year))
    box()
  }
  mtext(paste("Anopheles gambiae sampling locations:", species), outer = TRUE, line = 1)
  legend(
    "bottom",
    legend = month.abb,
    col = month_colors,
    pch = 21,
    pt.bg = month_colors,
    pt.cex = 1,
    ncol = 6,
    inset = 0.01,
    bty = "n",
    xpd = NA,
    title = "Sampling month"
  )
  dev.off()

  message(
    "Wrote ", length(years), " year panels to: ",
    normalizePath(output_path, mustWork = FALSE)
  )
}
