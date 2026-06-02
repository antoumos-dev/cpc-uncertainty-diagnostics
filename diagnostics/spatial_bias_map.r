.libPaths("/store_new/mch/msclim/share/CATs/cats/lib-R4.4.0/")
library(geocors)
library(raster)
library(dplyr)
.libPaths("/store_new/mch/msclim/sideris/R/lib/")
library(rgdal)
library(nowprecip)
library(nowtrack)
library(ced)
#library(mchradIO)
library(log4r)
library(RNetCDF)
library(epitools)
library(png)
library(RColorBrewer)
library(gstat)
library(rhdf5)
library(animation)
library(ncdf4)
library(lubridate)
library(dplyr)
library(data.table)
options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx4096m"))
library(XLConnect)


project_root    <- "/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats"
utils_file      <- file.path(project_root, "R", "utils.r")
plot_utils_file <- file.path(project_root, "R", "plot_utils.r")
source(utils_file)
source(plot_utils_file)
setwd(project_root)

threshold <- 5
years     <- sprintf("%02d", 16:25)

##### Load and combine all years
flatten_sp <- function(cv_data) {
  do.call(rbind, lapply(names(cv_data), function(t) {
    sp <- cv_data[[t]]
    if (!isS4(sp) || !inherits(sp, "Spatial")) return(NULL)
    df           <- sp@data
    xy           <- coordinates(sp)
    df$x         <- xy[, 1]; df$y <- xy[, 2]
    df$timestamp <- t
    df
  }))
}

all_on  <- vector("list", length(years))
all_off <- vector("list", length(years))

for (i in seq_along(years)) {
  YY       <- years[i]
  yr_label <- paste0("20", YY)
  f_on     <- file.path(project_root, paste0("cross_val_active_",   yr_label, ".rds"))
  f_off    <- file.path(project_root, paste0("cross_val_inactive_", yr_label, ".rds"))
  if (!file.exists(f_on) || !file.exists(f_off)) { message("Skipping ", yr_label); next }
  message("Loading ", yr_label, " ...")
  all_on[[i]]  <- flatten_sp(readRDS(f_on))
  all_off[[i]] <- flatten_sp(readRDS(f_off))
}

cv_on  <- do.call(rbind, Filter(Negate(is.null), all_on))
cv_off <- do.call(rbind, Filter(Negate(is.null), all_off))
yr_range <- paste0("20", years[1], "-20", years[length(years)])

##### Build the three subsets (ON and OFF) ---------------------------------

# 1. All wet hours
on_all  <- cv_on  %>% filter(raingauge > threshold)
off_all <- cv_off %>% filter(raingauge > threshold)

# 2. Active timestamps — timestamps where any correction fired
active_ts <- unique(cv_on$timestamp[cv_on$radar != cv_on$radar.orig])
on_active  <- cv_on  %>% filter(timestamp %in% active_ts, raingauge > threshold)
off_active <- cv_off %>% filter(timestamp %in% active_ts, raingauge > threshold)

# 3. Corrected gauge-timestamp pairs — only rows where the correction fired at that gauge
on_corrected  <- cv_on %>% filter(radar != radar.orig, raingauge > threshold)
corrected_key <- paste(on_corrected$nat_abbr, on_corrected$timestamp, sep = "_")
off_key       <- paste(cv_off$nat_abbr,       cv_off$timestamp,       sep = "_")
off_corrected <- cv_off %>% filter(off_key %in% corrected_key, raingauge > threshold)

subsets <- list(
  all_wet       = list(on = on_all,       off = off_all,       min_n_per_year = 10,
                       label = "All wet hours"),
  active_ts     = list(on = on_active,    off = off_active,    min_n_per_year = 5,
                       label = "Active timestamps"),
  corrected_pts = list(on = on_corrected, off = off_corrected, min_n_per_year =  2,
                       label = "Corrected gauge-timestamp pairs")
)

xlim_ch <- c(485, 835)
ylim_ch <- c(75,  295)
ch      <- load.map.elements(xlim = xlim_ch, ylim = ylim_ch)

##### Produce one PNG per subset -------------------------------------------

for (sname in names(subsets)) {
  s      <- subsets[[sname]]
  sb_on  <- station_bias(s$on,  s$min_n_per_year)
  sb_off <- station_bias(s$off, s$min_n_per_year)

  # Difference: ON bias − OFF bias at stations present in both
  sb_diff <- inner_join(
    sb_on  %>% select(nat_abbr, x, y, bias, rel_bias, mu_rg, n) %>%
               rename(bias_on = bias, rb_on = rel_bias, mu_rg_on = mu_rg, n_on = n),
    sb_off %>% select(nat_abbr,       bias, rel_bias, mu_rg, n) %>%
               rename(bias_off = bias, rb_off = rel_bias, mu_rg_off = mu_rg, n_off = n),
    by = "nat_abbr"
  ) %>%
    mutate(bias     = bias_on  - bias_off,
           rel_bias = rb_on    - rb_off,
           n        = pmin(n_on, n_off))

  message(sname, " — ON: ", nrow(sb_on), "  OFF: ", nrow(sb_off),
          "  matched diff: ", nrow(sb_diff), " stations")

  # ── Absolute bias PNG ──────────────────────────────────────────────────────
  cs_bias <- make_colorscale(c(sb_on$bias, sb_off$bias))
  cs_diff <- make_colorscale(sb_diff$bias)

  thr_txt <- gsub("\\.", "p", sprintf("%.2f", threshold))
  out_png <- file.path(project_root, paste0("spatial_bias_", sname, "_", yr_range, "_thr_", thr_txt, ".png"))
  png(out_png, width = 4800, height = 1400, res = 150, pointsize = 20)
  par(oma = c(0, 0, 5, 0), cex = 1.4, cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.3)
  layout(matrix(c(1, 2, 3, 4, 5), 1, 5), widths = c(10, 10, 1.2, 10, 1.2))
  draw_map_panel(sb_on,   cs_bias, paste("ON |",       s$label), ch, xlim_ch, ylim_ch)
  draw_map_panel(sb_off,  cs_bias, paste("OFF |",      s$label), ch, xlim_ch, ylim_ch)
  draw_colorbar(cs_bias)
  draw_map_panel(sb_diff, cs_diff, paste("ON − OFF |", s$label), ch, xlim_ch, ylim_ch)
  draw_colorbar(cs_diff)
  mtext(paste0("Average annual bias (kriging − raingauge) | ",
               s$label, " | raingauge > ", threshold, " mm/h | ", yr_range),
        outer = TRUE, line = 1.5, cex = 1.5)
  dev.off()
  message("Saved: ", out_png)

  # ── Relative bias PNG ──────────────────────────────────────────────────────
  cs_rb_bias <- make_colorscale(c(sb_on$rel_bias, sb_off$rel_bias))
  cs_rb_diff <- make_colorscale(sb_diff$rel_bias)

  sb_on_rb   <- sb_on;   sb_on_rb$bias   <- sb_on$rel_bias
  sb_off_rb  <- sb_off;  sb_off_rb$bias  <- sb_off$rel_bias
  sb_diff_rb <- sb_diff; sb_diff_rb$bias <- sb_diff$rel_bias

  out_png_rb <- file.path(project_root, paste0("spatial_relbias_", sname, "_", yr_range, "_thr_", thr_txt, ".png"))
  png(out_png_rb, width = 4800, height = 1400, res = 150, pointsize = 20)
  par(oma = c(0, 0, 5, 0), cex = 1.4, cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.3)
  layout(matrix(c(1, 2, 3, 4, 5), 1, 5), widths = c(10, 10, 1.2, 10, 1.2))
  draw_map_panel(sb_on_rb,   cs_rb_bias, paste("ON |",       s$label), ch, xlim_ch, ylim_ch)
  draw_map_panel(sb_off_rb,  cs_rb_bias, paste("OFF |",      s$label), ch, xlim_ch, ylim_ch)
  draw_colorbar(cs_rb_bias, unit = "–")
  draw_map_panel(sb_diff_rb, cs_rb_diff, paste("ON − OFF |", s$label), ch, xlim_ch, ylim_ch)
  draw_colorbar(cs_rb_diff, unit = "–")
  mtext(paste0("Average annual relative bias (bias / mean raingauge) | ",
               s$label, " | raingauge > ", threshold, " mm/h | ", yr_range),
        outer = TRUE, line = 1.5, cex = 1.5)
  dev.off()
  message("Saved: ", out_png_rb)
}
