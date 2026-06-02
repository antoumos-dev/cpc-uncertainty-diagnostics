.libPaths("/store_new/mch/msclim/share/CATs/cats/lib-R4.4.0/")
library(geocors)
library(fields)
library(robustbase)
library(datefuns)
.libPaths("/store_new/mch/msclim/sideris/R/lib/")
library(rgdal)
library(raster)
library(nowprecip)
library(nowtrack)
library(ced)
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

# ── Config ────────────────────────────────────────────────────────────────────
YEARS     <- sprintf("%02d", 23:24)   # "23" … "24"
threshold <- 0.5
period    <- "2023-2024"

project_root <- "/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats"
source(file.path(project_root, "R", "utils.r"))
source(file.path(project_root, "R", "plot_utils.r"))
setwd(project_root)

path_on <- "/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/"

parse_date_from_filename <- function(fname) {
  base <- sub("\\.rda$", "", basename(fname))
  as.POSIXct(strptime(paste0("20", substr(base,4,5), " ",
                              substr(base,6,8), " ", substr(base,9,10)),
                      format = "%Y %j %H", tz = "UTC"))
}

# as.numeric() on quantile strips the "25%" name so column headers stay Q1/Q3
fmt6 <- function(x) round(c(Mean   = mean(x, na.rm=TRUE),
                             Median = median(x, na.rm=TRUE),
                             Q1     = as.numeric(quantile(x, .25, na.rm=TRUE)),
                             Q3     = as.numeric(quantile(x, .75, na.rm=TRUE)),
                             SD     = sd(x, na.rm=TRUE),
                             N      = sum(!is.na(x))), 4)

# ── Per-year accumulators ─────────────────────────────────────────────────────
activations_by_year <- list()   # for Sheet 1
corr_rows_by_year   <- list()   # for Sheet 2
compare_long_all    <- list()   # for Sheet 3
result_full_all     <- list()   # for Sheets 4 & 5

for (YR in YEARS) {
  yr_label <- paste0("20", YR)
  message("Loading year: ", yr_label)

  all_data     <- readRDS(paste0("all_data_",          yr_label, ".rds"))
  all_data_off <- readRDS(paste0("all_data_conv_off_", yr_label, ".rds"))
  result       <- readRDS(paste0("result_conv_on_",    yr_label, ".rds"))

  rda_files <- list.files(path_on, pattern = paste0("CPC", YR, ".*\\.rda$"), full.names = FALSE)
  dates <- as.POSIXct(sapply(rda_files, parse_date_from_filename),
                      origin = "1970-01-01", tz = "UTC")

  # ── Sheet 1 data ───────────────────────────────────────────────────────────
  coeff_matrix <- do.call(cbind, lapply(all_data, function(df) {
    if (!is.data.frame(df)) return(NULL)
    df$coef.var
  }))
  at <- data.frame(
    year     = as.integer(yr_label),
    date     = dates,
    count    = as.numeric(colSums(coeff_matrix > threshold, na.rm = TRUE)),
    fraction = as.numeric(colMeans(coeff_matrix > threshold, na.rm = TRUE))
  )
  at$month  <- as.integer(format(dates, "%m"))
  at$season <- cut(at$month,
                   breaks = c(0, 2, 5, 8, 11, 12),
                   labels = c("DJF","MAM","JJA","SON","DJF2"), right = TRUE)
  levels(at$season)[levels(at$season) == "DJF2"] <- "DJF"
  activations_by_year[[YR]] <- at

  # ── Sheet 2 data ───────────────────────────────────────────────────────────
  cr <- do.call(rbind, lapply(seq_along(all_data), function(i) {
    df     <- all_data[[i]]
    if (!is.data.frame(df)) return(NULL)
    active <- !is.na(df$coef.var) & df$coef.var > threshold
    if (sum(active) == 0) return(NULL)
    corr     <- df$radar[active] - df$radar.orig[active]
    orig     <- df$radar.orig[active]
    rel_corr <- ifelse(orig > 0.1, corr / orig, NA_real_)
    data.frame(n_active      = sum(active),
               mean_corr     = mean(corr,          na.rm = TRUE),
               abs_mean_corr = mean(abs(corr),     na.rm = TRUE),
               min_corr      = min(corr,           na.rm = TRUE),
               max_corr      = max(corr,           na.rm = TRUE),
               mean_rel_corr = mean(rel_corr,      na.rm = TRUE),
               abs_rel_corr  = mean(abs(rel_corr), na.rm = TRUE))
  }))
  corr_rows_by_year[[YR]] <- cr

  # ── Sheet 3 data ───────────────────────────────────────────────────────────
  dates_on  <- names(all_data)
  dates_off <- names(all_data_off)
  common_on  <- which(dates_on  %in% dates_off)
  common_off <- which(dates_off %in% dates_on)

  cl <- rbindlist(lapply(seq_along(common_on), function(i) {
    df_on  <- all_data[[common_on[i]]]
    df_off <- all_data_off[[common_off[i]]]
    if (!is.data.frame(df_on) || !is.data.frame(df_off)) return(NULL)
    data.table(station_id = seq_len(nrow(df_on)),
               time       = dates_on[common_on[i]],
               radar_on   = df_on$radar,
               radar_off  = df_off$radar,
               cv_on      = df_on$coef.var,
               cv_off     = df_off$coef.var)
  }), fill = TRUE)
  compare_long_all[[YR]] <- cl

  # ── Sheets 4/5 data ────────────────────────────────────────────────────────
  coefvar_long <- rbindlist(lapply(seq_along(all_data), function(i) {
    df <- all_data[[i]]
    if (!is.data.frame(df)) return(NULL)
    data.table(station_id = seq_len(nrow(df)),
               time       = as.POSIXct(names(all_data)[i], format = "%Y-%m-%d %H:%M", tz = "UTC"),
               coef.var   = df$coef.var)
  }))

  setDT(result)
  rf <- merge(result, coefvar_long, by = c("station_id", "time"))
  rf[, active  := coef.var > threshold & !is.na(coef.var)]
  rf[, rel_unc := iqr / mu]
  result_full_all[[YR]] <- rf
}

# ── Pool all years ────────────────────────────────────────────────────────────
corr_rows_pool    <- rbindlist(corr_rows_by_year,   fill = TRUE)
compare_long_pool <- rbindlist(compare_long_all,    fill = TRUE)
result_full_pool  <- rbindlist(result_full_all,     fill = TRUE)

compare_long_pool[, delta_radar := radar_on - radar_off]
compare_long_pool[, delta_cv    := cv_on    - cv_off   ]

# ═══════════════════════════════════════════════════════════════════════════════
# SHEET 1 — Seasonal activation summary
# ═══════════════════════════════════════════════════════════════════════════════

seasonal_per_year <- do.call(rbind, lapply(YEARS, function(YR) {
  yr_label <- paste0("20", YR)
  at <- activations_by_year[[YR]]
  do.call(rbind, lapply(c("DJF","MAM","JJA","SON"), function(s) {
    sub <- at[at$season == s, ]
    data.frame(
      Year          = as.integer(yr_label),
      Season        = s,
      N_Timesteps   = nrow(sub),
      Mean_Count    = round(mean(sub$count), 3),
      Mean_Fraction = round(mean(sub$fraction, na.rm = TRUE), 4),
      Max_Count     = if (nrow(sub) > 0) max(sub$count) else NA_integer_,
      N_Active_Ts   = sum(sub$count > 0),
      Pct_Active_Ts = round(mean(sub$count > 0) * 100, 1)
    )
  }))
}))

setDT(seasonal_per_year)
seasonal_avg <- seasonal_per_year[, .(
  N_Years            = .N,
  Mean_N_Timesteps   = round(mean(N_Timesteps,   na.rm = TRUE), 1),
  Mean_Mean_Count    = round(mean(Mean_Count,    na.rm = TRUE), 3),
  Mean_Mean_Fraction = round(mean(Mean_Fraction, na.rm = TRUE), 4),
  Mean_Max_Count     = round(mean(Max_Count,     na.rm = TRUE), 1),
  Mean_N_Active_Ts   = round(mean(N_Active_Ts,   na.rm = TRUE), 1),
  Mean_Pct_Active_Ts = round(mean(Pct_Active_Ts, na.rm = TRUE), 1)
), by = Season][order(match(Season, c("DJF","MAM","JJA","SON")))]

# ═══════════════════════════════════════════════════════════════════════════════
# SHEET 2 — Correction statistics (pooled)
# ═══════════════════════════════════════════════════════════════════════════════

correction_stats <- data.frame(
  Metric = c(
    "Total timesteps with activations",
    "Overall mean correction (mm/h)",
    "Overall mean abs correction (mm/h)",
    "Overall min correction (mm/h)",
    "Overall max correction (mm/h)",
    "Weighted mean correction (mm/h)",
    "Weighted mean abs correction (mm/h)",
    "Weighted mean relative correction",
    "Weighted mean abs relative correction"
  ),
  Value = round(c(
    nrow(corr_rows_pool),
    mean(corr_rows_pool$mean_corr),
    mean(corr_rows_pool$abs_mean_corr),
    min(corr_rows_pool$min_corr),
    max(corr_rows_pool$max_corr),
    weighted.mean(corr_rows_pool$mean_corr,     corr_rows_pool$n_active),
    weighted.mean(corr_rows_pool$abs_mean_corr, corr_rows_pool$n_active),
    weighted.mean(corr_rows_pool$mean_rel_corr, corr_rows_pool$n_active, na.rm = TRUE),
    weighted.mean(corr_rows_pool$abs_rel_corr,  corr_rows_pool$n_active, na.rm = TRUE)
  ), 4),
  Notes = c(
    paste0(nrow(corr_rows_pool), " timesteps pooled across ", length(YEARS), " years"),
    "negative = radar decreased", "magnitude regardless of sign",
    "single grid cell, single timestep (any year)", "single grid cell, single timestep (any year)",
    "weighted by n_active per timestep", "weighted magnitude",
    "relative to radar.orig, cells with orig > 0.1", "weighted magnitude"
  )
)

# ═══════════════════════════════════════════════════════════════════════════════
# SHEET 3 — Direct ON vs OFF: delta radar and delta CoV (pooled)
# ═══════════════════════════════════════════════════════════════════════════════

wet          <- compare_long_pool[!is.na(radar_on) & !is.na(radar_off) & (radar_on > 0 | radar_off > 0)]
active_cells <- wet[cv_on > threshold]

on_vs_off <- rbind(
  data.frame(Subset="All wet cells",          Variable="delta_radar (ON-OFF, mm/h)", t(fmt6(wet$delta_radar))),
  data.frame(Subset="All wet cells",          Variable="delta_CoV   (ON-OFF)",       t(fmt6(wet$delta_cv))),
  data.frame(Subset="Active cells (CoV>0.5)", Variable="delta_radar (ON-OFF, mm/h)", t(fmt6(active_cells$delta_radar))),
  data.frame(Subset="Active cells (CoV>0.5)", Variable="delta_CoV   (ON-OFF)",       t(fmt6(active_cells$delta_cv)))
)

# ═══════════════════════════════════════════════════════════════════════════════
# SHEET 4 — μ, IQR, IQR/μ: active vs inactive (pooled)
# ═══════════════════════════════════════════════════════════════════════════════

mu_iqr_summary <- rbind(
  data.frame(Variable="mu (mm/h)",  Group="Active",   t(fmt6(result_full_pool$mu[ result_full_pool$active]))),
  data.frame(Variable="mu (mm/h)",  Group="Inactive", t(fmt6(result_full_pool$mu[!result_full_pool$active]))),
  data.frame(Variable="IQR (mm/h)", Group="Active",   t(fmt6(result_full_pool$iqr[ result_full_pool$active]))),
  data.frame(Variable="IQR (mm/h)", Group="Inactive", t(fmt6(result_full_pool$iqr[!result_full_pool$active]))),
  data.frame(Variable="IQR/mu",     Group="Active",   t(fmt6(result_full_pool$rel_unc[ result_full_pool$active]))),
  data.frame(Variable="IQR/mu",     Group="Inactive", t(fmt6(result_full_pool$rel_unc[!result_full_pool$active])))
)

# ═══════════════════════════════════════════════════════════════════════════════
# SHEET 5 — Binned intensity: IQR/mu by decile at top 10% of μ (pooled)
# ═══════════════════════════════════════════════════════════════════════════════

mu_p90  <- quantile(result_full_pool$mu, 0.90, na.rm = TRUE)
intense <- result_full_pool[!is.na(mu) & mu >= mu_p90]
intense[, mu_bin := cut(mu,
                        breaks = quantile(mu, probs = seq(0, 1, 0.1), na.rm = TRUE),
                        include.lowest = TRUE, labels = FALSE)]

bin_comparison <- intense[, .(
  Bin_mu_range           = paste0("[", round(min(mu, na.rm=TRUE),2), ", ",
                                       round(max(mu, na.rm=TRUE),2), "]"),
  Median_mu              = round(median(mu,  na.rm = TRUE), 3),
  N_Active               = sum(active == TRUE,  na.rm = TRUE),
  N_Inactive             = sum(active == FALSE, na.rm = TRUE),
  Median_relunc_active   = round(median(rel_unc[active == TRUE],  na.rm = TRUE), 4),
  Median_relunc_inactive = round(median(rel_unc[active == FALSE], na.rm = TRUE), 4)
), by = mu_bin][order(mu_bin)]

bin_comparison[, Delta_relunc    := round(Median_relunc_active - Median_relunc_inactive, 4)]
bin_comparison[, Pct_diff_relunc := round((Median_relunc_active - Median_relunc_inactive) /
                                           Median_relunc_inactive * 100, 2)]

# ═══════════════════════════════════════════════════════════════════════════════
# WRITE EXCEL
# ═══════════════════════════════════════════════════════════════════════════════

out_xlsx <- file.path(project_root, paste0("conv_control_stats_", period, ".xlsx"))
wb <- loadWorkbook(out_xlsx, create = TRUE)

write_section <- function(wb, sheet, header, df) {
  createSheet(wb, name = sheet)
  writeWorksheet(wb, data.frame(V1 = header), sheet = sheet,
                 startRow = 1, startCol = 1, header = FALSE)
  writeWorksheet(wb, df, sheet = sheet, startRow = 3, startCol = 1)
}

write_section(wb, "1_Activation_Seasonal",
  paste("Seasonal activation | multi-year average | CoV >", threshold, "|", period),
  as.data.frame(seasonal_avg))

write_section(wb, "1b_Activation_PerYear",
  paste("Seasonal activation | per-year breakdown | CoV >", threshold),
  as.data.frame(seasonal_per_year))

write_section(wb, "2_Correction_Stats",
  paste("Correction statistics | radar − radar.orig at active cells | pooled", period),
  correction_stats)

write_section(wb, "3_ON_vs_OFF",
  paste("Direct ON vs OFF | delta radar and delta CoV | pooled", period),
  on_vs_off)

write_section(wb, "4_mu_IQR_summary",
  paste("mu, IQR, IQR/mu: active vs inactive | pooled", period),
  mu_iqr_summary)

write_section(wb, "5_Intensity_bins",
  paste("IQR/mu by intensity decile (top 10% of mu) | active vs inactive | pooled", period),
  as.data.frame(bin_comparison))

saveWorkbook(wb)
message("Saved: ", out_xlsx)
