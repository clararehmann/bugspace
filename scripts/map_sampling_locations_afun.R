#!/usr/bin/env Rscript

# Generate year-panel sampling maps for A. funestus metadata.

if (!requireNamespace("maps", quietly = TRUE)) {
  stop(
    "Install the required package before running this script: maps\n",
    "For example: install.packages(\"maps\")",
    call. = FALSE
  )
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2) {
  stop(
    "Usage: Rscript scripts/map_sampling_locations_afun.R \\\n  [metadata.csv] [output.png]",
    call. = FALSE
  )
}

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
  file.path(project_root, "data", "vo_afun_release", "v1.0", "metadata", "all_samples_metadata.csv")
}
output_path <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(project_root, "sampling_locations_afun.png")
}

if (!file.exists(metadata_path)) {
  stop("Metadata file does not exist: ", metadata_path, call. = FALSE)
}
samples <- read.csv(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("year", "month", "latitude", "longitude")
if (!all(required_columns %in% names(samples))) {
  stop(
    "The input must contain columns named ",
    paste(required_columns, collapse = ", "),
    call. = FALSE
  )
}

for (column in required_columns) {
  samples[[column]] <- suppressWarnings(as.numeric(samples[[column]]))
}
valid_rows <- with(
  samples,
  is.finite(latitude) & latitude >= -90 & latitude <= 90 &
    is.finite(longitude) & longitude >= -180 & longitude <= 180 &
    is.finite(year) & year > 0 &
    is.finite(month) & month >= 1 & month <= 12
)
if (any(!valid_rows)) {
  warning(sum(!valid_rows), " rows have invalid coordinates, year, or month and will be omitted.")
  samples <- samples[valid_rows, , drop = FALSE]
}
if (nrow(samples) == 0) {
  stop("No valid sampling records were found.", call. = FALSE)
}

longitude_range <- range(samples$longitude)
latitude_range <- range(samples$latitude)
longitude_padding <- max(diff(longitude_range) * 0.05, 2)
latitude_padding <- max(diff(latitude_range) * 0.05, 2)
x_limits <- c(longitude_range[1] - longitude_padding, longitude_range[2] + longitude_padding)
y_limits <- c(latitude_range[1] - latitude_padding, latitude_range[2] + latitude_padding)
years <- sort(unique(samples$year))
panel_columns <- min(4, max(1, ceiling(sqrt(length(years)))))
panel_rows <- ceiling(length(years) / panel_columns)
month_colors <- grDevices::hcl.colors(12, palette = "Dynamic")

output_directory <- dirname(output_path)
if (!dir.exists(output_directory)) {
  dir.create(output_directory, recursive = TRUE)
}

png(output_path, width = 600 * panel_columns, height = 520 * panel_rows, res = 150)
par(mfrow = c(panel_rows, panel_columns), mar = c(2.5, 2.5, 2.5, 1), oma = c(3, 0, 3, 0))
for (year in years) {
  year_samples <- samples[samples$year == year, , drop = FALSE]
  site_key <- paste(year_samples$latitude, year_samples$longitude, year_samples$month, sep = "_")
  site_counts <- as.data.frame(table(site_key), stringsAsFactors = FALSE)
  site_coordinates <- unique(year_samples[c("latitude", "longitude", "month")])
  site_coordinates$site_key <- paste(site_coordinates$latitude, site_coordinates$longitude, site_coordinates$month, sep = "_")
  sites <- merge(site_coordinates, site_counts, by = "site_key", sort = FALSE)
  names(sites)[names(sites) == "Freq"] <- "sample_count"

  maps::map("world", xlim = x_limits, ylim = y_limits, fill = TRUE, col = "white", border = "grey55", bg = "aliceblue")
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
mtext("Anopheles funestus sampling locations", outer = TRUE, line = 1)
legend("bottom", legend = month.abb, col = month_colors, pch = 21, pt.bg = month_colors, pt.cex = 1, ncol = 6, inset = 0.01, bty = "n", xpd = NA, title = "Sampling month")
dev.off()

message(
  "Wrote ", length(years), " year panels to: ",
  normalizePath(output_path, mustWork = FALSE)
)
