#!/usr/bin/env Rscript

# 2x2 panel plot (DJF, MAM / JJA, SON) of Swiss-mean seasonal precipitation
# as percentage difference w.r.t. each series' own 1991-2020 mean, for all RCP8.5 model chains, from the txt files written by
# seasonal_precip_R.R. The multi-chain mean is added as a bold black line and
# the observed Swiss mean (1981-2020) as a red line.

RESDIR <- "/highres_nobackup/svenk/CH2025_check_transient/RESULTS"
PLOTDIR <- "/highres_nobackup/svenk/CH2025_check_transient/PLOTS"
OUTFILE <- file.path(PLOTDIR, "seasonal_pr_CH_rcp85_pctdiff_1991-2020.png")
OBSFILE <- file.path(RESDIR, "seasonal_pr_CH_obs_1981-2020.txt")
SEASONS <- c("DJF", "MAM", "JJA", "SON") # row-wise: UL, UR, LL, LR

files <- Sys.glob(file.path(RESDIR, "seasonal_pr_CH_*_rcp85.txt"))
if (length(files) == 0L) stop("No seasonal result files in ", RESDIR, call. = FALSE)
chains <- toupper(sub("^seasonal_pr_CH_(.*)_rcp85\\.txt$", "\\1", basename(files)))
data <- lapply(files, read.table, header = TRUE)
names(data) <- chains

obs <- read.table(OBSFILE, header = TRUE)

years <- sort(unique(unlist(lapply(data, `[[`, "year"))))
REF <- 1991:2020
pct <- function(y, v) 100 * (v / mean(v[y %in% REF], na.rm = TRUE) - 1)
cols <- grDevices::hcl.colors(length(chains), "Dark 3")

dir.create(PLOTDIR, showWarnings = FALSE, recursive = TRUE)
png(OUTFILE, width = 2200, height = 1700, res = 200)
layout(rbind(matrix(1:4, 2, 2, byrow = TRUE), 5), heights = c(1, 1, 0.45))
par(mar = c(4, 4.5, 2.5, 1), oma = c(0, 0, 3, 0))

for (s in SEASONS) {
  mat <- sapply(data, function(d) pct(d$year, d[[s]])[match(years, d$year)])
  plot(range(years), range(mat, na.rm = TRUE, finite = TRUE), type = "n", xlab = "Year",
       ylab = "Precipitation change vs 1991-2020 (%)", main = s)
  grid()
  abline(h = 0, lty = 2)
  for (i in seq_along(chains)) lines(years, mat[, i], col = cols[i], lwd = 1)
  lines(years, rowMeans(mat, na.rm = TRUE), col = "black", lwd = 3)
  lines(obs$year, pct(obs$year, obs[[s]]), col = "red", lwd = 2)
}

par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", legend = c(chains, "MEAN", "OBS (1981-2020)"),
       col = c(cols, "black", "red"), lwd = c(rep(1.5, length(chains)), 3, 2), ncol = 5, cex = 0.7, bty = "n")
mtext("Swiss mean precipitation, % diff. to 1991-2020, CH2025 transient (non-BTS-shifted)",
      outer = TRUE, side = 3, line = 0.5, cex = 1.2, font = 2)
invisible(dev.off())
message("Wrote ", OUTFILE)
