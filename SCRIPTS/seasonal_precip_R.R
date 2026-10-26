#!/usr/bin/env Rscript

# Swiss-mean seasonal precipitation (mm/day) for all RCP8.5 model chains, read
# from the pre-computed seasonal means (1981-2099). One txt file per model
# chain with columns year, DJF, MAM, JJA, SON.
# DJF is labelled by the year of its December (DJF 2050 = Dec 2050 - Feb 2051).
# The incomplete first (Jan-Feb 1981) and last (Dec 2099) winter are dropped.
#
# The observational grid (MeteoSwiss daily precipitation on a 0.02 deg lon/lat
# grid) is processed the same way: cos(lat)-weighted Swiss mean of the daily
# values, averaged over complete seasons 1981-2020, written to
# seasonal_pr_CH_obs_1981-2020.txt.

suppressPackageStartupMessages(library(ncdf4))

BASEDIR <- "/net/stratus/c2sm-data/MCH_CH2025/indicators/output/qm/gridded/EUR11/pr/RCP85"
RESDIR <- "/highres_nobackup/svenk/CH2025_check_transient/RESULTS"
SEASONS <- c("DJF", "MAM", "JJA", "SON")
OBSFILE <- "/net/co2/c2sm-data/ch202X/Obs_Data/MCH/processed/pr_ch02.lonlat_1971-2021.nc"
OBS_YEARS <- 1981:2020

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

# ---- observations -----------------------------------------------------------
obs_swiss_mean_daily <- function(path, chunk = 1000L) {
  nc <- nc_open(path)
  on.exit(nc_close(nc))
  lat <- ncvar_get(nc, "lat")
  nlon <- nc$dim$lon$len
  nlat <- nc$dim$lat$len
  ntime <- nc$dim$time$len
  w <- rep(cos(lat * pi / 180), each = nlon) # lon varies fastest
  origin <- as.Date(sub("days since ([0-9-]+).*", "\\1",
                        ncatt_get(nc, "time", "units")$value))
  dates <- origin + floor(ncvar_get(nc, "time"))

  means <- numeric(ntime)
  for (start in seq(1L, ntime, by = chunk)) {
    n <- min(chunk, ntime - start + 1L)
    pr <- ncvar_get(nc, "pr", start = c(1, 1, start), count = c(nlon, nlat, n),
                    collapse_degen = FALSE)
    dim(pr) <- c(nlon * nlat, n)
    ok <- !is.na(pr)
    pr[!ok] <- 0
    means[start:(start + n - 1L)] <- colSums(pr * w) / colSums(ok * w)
  }
  data.frame(date = dates, value = means)
}

obs_seasonal <- function(daily, years) {
  month <- as.integer(format(daily$date, "%m"))
  year <- as.integer(format(daily$date, "%Y"))
  season <- SEASONS[(month %% 12) %/% 3 + 1L]
  year[month %in% 1:2] <- year[month %in% 1:2] - 1L # DJF labelled by December
  d <- data.frame(year = year, season = season, value = daily$value)
  out <- data.frame(year = years)
  for (s in SEASONS) {
    sub <- d[d$season == s & !is.na(d$value), ]
    agg <- tapply(sub$value, sub$year, function(v) if (length(v) >= 85) mean(v) else NA)
    out[[s]] <- round(as.numeric(agg[as.character(years)]), 3)
  }
  out
}

obs <- obs_seasonal(obs_swiss_mean_daily(OBSFILE), OBS_YEARS)
obsfile <- file.path(RESDIR, sprintf("seasonal_pr_CH_obs_%d-%d.txt",
                                     min(OBS_YEARS), max(OBS_YEARS)))
write.table(obs, obsfile, quote = FALSE, row.names = FALSE, sep = " ", na = "NA")
message("Wrote ", obsfile)
