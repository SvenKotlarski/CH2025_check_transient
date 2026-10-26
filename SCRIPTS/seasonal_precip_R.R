#!/usr/bin/env Rscript

# Swiss-mean seasonal precipitation (mm/day) for all RCP8.5 model chains, read
# from the pre-computed seasonal means (1981-2099). One txt file per model
# chain with columns year, DJF, MAM, JJA, SON.
# DJF is labelled by the year of its December (DJF 2050 = Dec 2050 - Feb 2051).
# The incomplete first (Jan-Feb 1981) and last (Dec 2099) winter are dropped.

suppressPackageStartupMessages(library(ncdf4))

BASEDIR <- "/net/stratus/c2sm-data/MCH_CH2025/indicators/output/qm/gridded/EUR11/pr/RCP85"
RESDIR <- "/highres_nobackup/svenk/CH2025_check_transient/RESULTS"
SEASONS <- c("DJF", "MAM", "JJA", "SON")

swiss_mean_seasonal <- function(path) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  dims <- vapply(nc$var[["pr"]]$dim, function(d) d$name, character(1))
  tpos <- match("time", dims)
  pr <- ncvar_get(nc, "pr", collapse_degen = FALSE)
  pr <- aperm(pr, c(tpos, setdiff(seq_along(dims), tpos)))
  dim(pr) <- c(dim(pr)[1], prod(dim(pr)[-1]))
  # cells outside Switzerland are missing; equal-area 1 km grid
  means <- rowMeans(pr, na.rm = TRUE)

  origin <- as.Date(sub("days since ([0-9-]+).*", "\\1",
                        ncatt_get(nc, "time", "units")$value))
  bnds <- ncvar_get(nc, "time_bnds")
  length_days <- bnds[2, ] - bnds[1, ] + 1
  dates <- origin + floor(ncvar_get(nc, "time"))
  month <- as.integer(format(dates, "%m"))
  year <- as.integer(format(dates, "%Y"))

  # seasons are centred in Jan (DJF), Apr (MAM), Jul (JJA), Oct (SON);
  # the truncated first/last DJF steps are centred in Jan/Dec but are short
  season <- SEASONS[match(month, c(1, 4, 7, 10), nomatch = NA)]
  season[month == 12] <- "DJF"
  complete <- length_days >= 85
  year[season == "DJF" & month == 1] <- year[season == "DJF" & month == 1] - 1L
  data.frame(year = year, season = season, value = means)[complete, ]
}

files <- Sys.glob(file.path(BASEDIR, "*", "*", "seasonal", "pr_*_seasonal_mean.nc"))
if (length(files) == 0L) stop("No seasonal files found in ", BASEDIR, call. = FALSE)
dir.create(RESDIR, showWarnings = FALSE, recursive = TRUE)

for (f in files) {
  chain <- tolower(paste(basename(dirname(dirname(dirname(f)))),
                         basename(dirname(dirname(f))), sep = "-"))
  res <- swiss_mean_seasonal(f)
  years <- sort(unique(res$year))
  out <- data.frame(year = years)
  for (s in SEASONS) {
    sub <- res[res$season == s, ]
    out[[s]] <- round(sub$value[match(years, sub$year)], 3)
  }
  outfile <- file.path(RESDIR, sprintf("seasonal_pr_CH_%s_rcp85.txt", chain))
  write.table(out, outfile, quote = FALSE, row.names = FALSE, sep = " ", na = "NA")
  message("Wrote ", outfile)
}
