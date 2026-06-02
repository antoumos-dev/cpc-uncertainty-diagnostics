.libPaths("/store_new/mch/msclim/share/CATs/cats/lib-R4.4.0/")
#library(mchdwh)
library(geocors)
library(fields)
library(robustbase)
library(datefuns)
#library(lattice)
.libPaths("/store_new/mch/msclim/sideris/R/lib/")
library(rgdal)
library(raster)
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

##### Helper functions
rmse <- function(x, y) sqrt(mean((x - y)^2, na.rm = TRUE))
bias <- function(x, y) mean(x - y,           na.rm = TRUE)
mae  <- function(x, y) mean(abs(x - y),      na.rm = TRUE)

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

metrics_df <- function(on_kr, on_rg, off_kr, off_rg) {
  data.frame(
    Metric = c("Bias", "MAE", "MSE", "RMSE", "Corr", "N"),
    ON  = c(bias(on_kr, on_rg), mae(on_kr, on_rg),
            mean((on_kr - on_rg)^2, na.rm = TRUE), rmse(on_kr, on_rg),
            cor(on_kr, on_rg, use = "complete.obs"), length(on_kr)),
    OFF = c(bias(off_kr, off_rg), mae(off_kr, off_rg),
            mean((off_kr - off_rg)^2, na.rm = TRUE), rmse(off_kr, off_rg),
            cor(off_kr, off_rg, use = "complete.obs"), length(off_kr))
  )
}

ru_df <- function(on_df, off_df) {
  on_idr  <- compute_idr_points(on_df$kriging,  on_df$variance)
  off_idr <- compute_idr_points(off_df$kriging, off_df$variance)
  on_rv   <- on_idr[ on_df$kriging  > 0.1] / on_df$kriging[ on_df$kriging  > 0.1]
  off_rv  <- off_idr[off_df$kriging > 0.1] / off_df$kriging[off_df$kriging > 0.1]
  df <- as.data.frame(rbind(as.numeric(summary(on_rv))[1:6], as.numeric(summary(off_rv))[1:6]))
  colnames(df) <- c("Min", "Q1", "Median", "Mean", "Q3", "Max")
  cbind(Case = c("ON", "OFF"), df)
}

# Bin kriging values and compute mean variance and IQR per bin — reveals how
# uncertainty tracks intensity and whether ON/OFF diverge at high mu
mu_bins <- c(0, 0.5, 1, 2, 5, 10, 20, Inf)
bin_labels <- c("0-0.5", "0.5-1", "1-2", "2-5", "5-10", "10-20", ">20")

varmu_df <- function(on_df, off_df) {
  bin_one <- function(df, case) {
    df <- df[df$kriging > 0 & !is.na(df$kriging) & !is.na(df$variance), ]
    b  <- cut(df$kriging, breaks = mu_bins, labels = bin_labels, right = FALSE)
    ag <- aggregate(cbind(variance, kriging) ~ b, data = df, FUN = mean,  na.rm = TRUE)
    n  <- aggregate(kriging ~ b,               data = df, FUN = length)
    ag$N    <- n$kriging
    ag$ru   <- ag$variance / ag$kriging
    ag$Case <- case
    ag
  }
  rbind(bin_one(on_df, "ON"), bin_one(off_df, "OFF"))
}

##### Year loop and data processing
years     <- sprintf("%02d", 16:25)
threshold <- 5.0 
sheet_names <- c("All_ts", "Active_ts", "Active_points", "Ticino_all", "Ticino_active")

metrics_by_year <- setNames(lapply(sheet_names, function(s) list()), sheet_names)
ru_by_year      <- setNames(lapply(sheet_names, function(s) list()), sheet_names)
varmu_by_year   <- setNames(lapply(sheet_names, function(s) list()), sheet_names)

for (YY in years) {
  yr_label <- paste0("20", YY)
  f_on  <- file.path(project_root, paste0("cross_val_active_",   yr_label, ".rds"))
  f_off <- file.path(project_root, paste0("cross_val_inactive_", yr_label, ".rds"))

  if (!file.exists(f_on) || !file.exists(f_off)) {
    message("Skipping ", yr_label, ": files not found")
    next
  }
  message("Processing ", yr_label, " ...")

  cv_all_on  <- flatten_sp(readRDS(f_on))
  cv_all_off <- flatten_sp(readRDS(f_off))

  common_ts  <- intersect(unique(cv_all_on$timestamp), unique(cv_all_off$timestamp))
  on_common  <- cv_all_on[ cv_all_on$timestamp  %in% common_ts, ]
  off_common <- cv_all_off[cv_all_off$timestamp %in% common_ts, ]

  wet_on  <- on_common[ on_common$raingauge  > threshold, ]
  wet_off <- off_common[off_common$raingauge > threshold, ]

  ts_active  <- unique(cv_all_on$timestamp[cv_all_on$radar != cv_all_on$radar.orig])
  on_active  <- on_common[ on_common$timestamp  %in% ts_active & on_common$raingauge  > threshold, ]
  off_active <- off_common[off_common$timestamp %in% ts_active & off_common$raingauge > threshold, ]

  on_corrected  <- on_common[on_common$radar != on_common$radar.orig & on_common$raingauge > threshold, ]
  corrected_key <- paste(on_common$nat_abbr[on_common$radar != on_common$radar.orig],
                         on_common$timestamp[on_common$radar != on_common$radar.orig], sep = "_")
  off_key       <- paste(off_common$nat_abbr, off_common$timestamp, sep = "_")
  off_corrected <- off_common[off_key %in% corrected_key & off_common$raingauge > threshold, ]

  in_ticino <- function(df) df$x >= 670 & df$x <= 750 & df$y >= 70 & df$y <= 150
  on_ticino  <- on_common[ in_ticino(on_common)  & on_common$raingauge  > threshold, ]
  off_ticino <- off_common[in_ticino(off_common) & off_common$raingauge > threshold, ]

  ts_active_ticino <- unique(
    on_common$timestamp[on_common$radar != on_common$radar.orig & in_ticino(on_common)]
  )
  on_ticino_active  <- on_ticino[ on_ticino$timestamp  %in% ts_active_ticino, ]
  off_ticino_active <- off_ticino[off_ticino$timestamp %in% ts_active_ticino, ]

  metrics_by_year[["All_ts"]][[YY]]        <- metrics_df(wet_on$kriging,            wet_on$raingauge,            wet_off$kriging,            wet_off$raingauge)
  metrics_by_year[["Active_ts"]][[YY]]     <- metrics_df(on_active$kriging,         on_active$raingauge,         off_active$kriging,         off_active$raingauge)
  metrics_by_year[["Active_points"]][[YY]] <- metrics_df(on_corrected$kriging,      on_corrected$raingauge,      off_corrected$kriging,      off_corrected$raingauge)
  metrics_by_year[["Ticino_all"]][[YY]]    <- metrics_df(on_ticino$kriging,         on_ticino$raingauge,         off_ticino$kriging,         off_ticino$raingauge)
  metrics_by_year[["Ticino_active"]][[YY]] <- metrics_df(on_ticino_active$kriging,  on_ticino_active$raingauge,  off_ticino_active$kriging,  off_ticino_active$raingauge)

  ru_by_year[["All_ts"]][[YY]]        <- ru_df(wet_on,            wet_off)
  ru_by_year[["Active_ts"]][[YY]]     <- ru_df(on_active,         off_active)
  ru_by_year[["Active_points"]][[YY]] <- ru_df(on_corrected,      off_corrected)
  ru_by_year[["Ticino_all"]][[YY]]    <- ru_df(on_ticino,         off_ticino)
  ru_by_year[["Ticino_active"]][[YY]] <- ru_df(on_ticino_active,  off_ticino_active)

  varmu_by_year[["All_ts"]][[YY]]        <- varmu_df(wet_on,            wet_off)
  varmu_by_year[["Active_ts"]][[YY]]     <- varmu_df(on_active,         off_active)
  varmu_by_year[["Active_points"]][[YY]] <- varmu_df(on_corrected,      off_corrected)
  varmu_by_year[["Ticino_all"]][[YY]]    <- varmu_df(on_ticino,         off_ticino)
  varmu_by_year[["Ticino_active"]][[YY]] <- varmu_df(on_ticino_active,  off_ticino_active)
}

##### Average metrics across years (N is summed; all other metrics are averaged)
mean_metrics <- function(mlist) {
  dfs <- Filter(Negate(is.null), mlist) # Remove NULL entries (for years that were skipped)
  if (length(dfs) == 0) return(NULL)
  metric_names <- dfs[[1]]$Metric
  is_N         <- metric_names == "N"
  on_mat       <- sapply(dfs, function(d) d$ON)
  off_mat      <- sapply(dfs, function(d) d$OFF)
  on_agg       <- rowMeans(on_mat,  na.rm = TRUE)
  off_agg      <- rowMeans(off_mat, na.rm = TRUE)
  on_agg[is_N]  <- rowSums(on_mat[ is_N, , drop = FALSE], na.rm = TRUE)
  off_agg[is_N] <- rowSums(off_mat[is_N, , drop = FALSE], na.rm = TRUE)
  data.frame(Metric = metric_names, ON = on_agg, OFF = off_agg)
}

mean_ru <- function(rulist) {
  dfs <- Filter(Negate(is.null), rulist)
  if (length(dfs) == 0) return(NULL)
  avg_case <- function(case) {
    rows <- do.call(rbind, lapply(dfs, function(d) d[d$Case == case, -1]))
    as.data.frame(lapply(rows, mean, na.rm = TRUE))
  }
  cbind(Case = c("ON", "OFF"), rbind(avg_case("ON"), avg_case("OFF")))
}

mean_varmu <- function(vmlist) {
  dfs <- Filter(Negate(is.null), vmlist)
  if (length(dfs) == 0) return(NULL)
  combined <- do.call(rbind, dfs)
  # Average variance, kriging, ru per bin and case across years (weight by N)
  result <- do.call(rbind, lapply(c("ON", "OFF"), function(case) {
    sub <- combined[combined$Case == case, ]
    ag  <- aggregate(cbind(variance, kriging, ru, N) ~ b, data = sub, FUN = mean, na.rm = TRUE)
    ag$Case <- case
    ag
  }))
  result
}

##### Write to Excel
yr_first  <- paste0("20", years[1])
yr_last   <- paste0("20", years[length(years)])
yr_range  <- paste0(yr_first, "-", yr_last)

thr_label <- paste0("raingauge > ", threshold, " mm")
sheet_labels <- list(
  All_ts        = paste("All timestamps —",       thr_label, "(mean", yr_range, ")"),
  Active_ts     = paste("Active timestamps —",    thr_label, "(mean", yr_range, ")"),
  Active_points = paste("Corrected points only —",thr_label, "(mean", yr_range, ")"),
  Ticino_all    = paste("Ticino — all timestamps (mean", yr_range, ")"),
  Ticino_active = paste("Ticino — active timestamps (mean", yr_range, ")")
)

write_sheet <- function(wb, sheet, df_metrics, df_ru, df_varmu, label) {
  createSheet(wb, name = sheet)
  writeWorksheet(wb, data.frame(Section = label),                                    sheet = sheet, startRow =  1, startCol = 1, header = FALSE)
  writeWorksheet(wb, df_metrics,                                                      sheet = sheet, startRow =  3, startCol = 1)
  writeWorksheet(wb, data.frame(Section = "Relative uncertainty (variance/kriging)"), sheet = sheet, startRow = 12, startCol = 1, header = FALSE)
  writeWorksheet(wb, df_ru,                                                           sheet = sheet, startRow = 13, startCol = 1)
  if (!is.null(df_varmu)) {
    writeWorksheet(wb, data.frame(Section = "Variance vs mu (binned kriging)"),       sheet = sheet, startRow = 17, startCol = 1, header = FALSE)
    writeWorksheet(wb, df_varmu,                                                       sheet = sheet, startRow = 18, startCol = 1)
  }
}

out_xlsx <- file.path(project_root, paste0("cv_scores_mean_", yr_range, "_", threshold, ".xlsx"))
if (file.exists(out_xlsx)) file.remove(out_xlsx)
wb <- loadWorkbook(out_xlsx, create = TRUE)

for (sh in sheet_names) {
  write_sheet(wb, sh,
              mean_metrics(metrics_by_year[[sh]]),
              mean_ru(ru_by_year[[sh]]),
              mean_varmu(varmu_by_year[[sh]]),
              sheet_labels[[sh]])
}

saveWorkbook(wb)
message("Done. Output: ", out_xlsx)


########
######## Plotting variance vs mu (binned kriging) to visualize how uncertainty tracks intensity and whether ON/OFF diverge at high mu 
########

##### variance ~ mu plots

plot_varmu <- function(vm, title) {
  if (is.null(vm)) return(invisible(NULL))
  on_vm  <- vm[vm$Case == "ON",  ]
  off_vm <- vm[vm$Case == "OFF", ]
  bins_present <- levels(factor(vm$b, levels = bin_labels))
  x <- seq_along(bins_present)

  ylim_var <- range(c(on_vm$variance, off_vm$variance), na.rm = TRUE)
  ylim_ru  <- range(c(on_vm$ru,       off_vm$ru),       na.rm = TRUE)

  par(mfrow = c(1, 2), mar = c(5, 4, 3, 1))

  # Panel 1: mean variance vs mean kriging (mu)
  plot(on_vm$kriging,  on_vm$variance,  type = "b", pch = 16, col = "steelblue",
       xlab = "Mean kriging (mm)", ylab = "Mean variance (mm)",
       main = paste(title, "\nVariance vs mu"), ylim = ylim_var)
  lines(off_vm$kriging, off_vm$variance, type = "b", pch = 17, col = "tomato")
  legend("topleft", legend = c("ON", "OFF"), col = c("steelblue", "tomato"),
         pch = c(16, 17), lty = 1, bty = "n")

  # Panel 2: relative uncertainty (variance/kriging) vs bin
  on_x  <- match(on_vm$b,  bins_present)
  off_x <- match(off_vm$b, bins_present)
  plot(on_x,  on_vm$ru,  type = "b", pch = 16, col = "steelblue",
       xaxt = "n", xlab = "Kriging bin (mm)", ylab = "Mean variance / kriging",
       main = paste(title, "\nRelative uncertainty vs mu"), ylim = ylim_ru)
  lines(off_x, off_vm$ru, type = "b", pch = 17, col = "tomato")
  axis(1, at = x, labels = bins_present, cex.axis = 0.8)
  legend("topleft", legend = c("ON", "OFF"), col = c("steelblue", "tomato"),
         pch = c(16, 17), lty = 1, bty = "n")
}

out_pdf <- file.path(project_root, paste0("varmu_", yr_range, "_thr", threshold, ".pdf"))
pdf(out_pdf, width = 10, height = 5)
for (sh in sheet_names) {
  vm <- mean_varmu(varmu_by_year[[sh]])
  plot_varmu(vm, sh)
}
dev.off()
message("Variance-mu plots: ", out_pdf)
