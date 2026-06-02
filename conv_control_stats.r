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
YEAR      <- "24"
threshold <- 0.5
yr_label  <- paste0("20", YEAR)

project_root    <- "/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats"
utils_file      <- file.path(project_root, "R", "utils.r")
plot_utils_file <- file.path(project_root, "R", "plot_utils.r")
source(utils_file)
source(plot_utils_file)
setwd(project_root)

# ── Load pre-processed data (produced by conv_control.r) ──────────────────────
all_data     <- readRDS(paste0("all_data_", yr_label, ".rds"))          # CC ON
all_data_off <- readRDS(paste0("all_data_conv_off_", yr_label, ".rds")) # CC OFF
result       <- readRDS(paste0("result_conv_on_", yr_label, ".rds"))                 # nearest-grid μ/IQR, CC ON
result_off   <- readRDS(paste0("result_conv_off_", yr_label, ".rds"))                # nearest-grid μ/IQR, CC OFF

path_on  <- "/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/"
rda_files <- list.files(path_on, pattern = paste0("CPC", YEAR, ".*\\.rda$"), full.names = FALSE)

parse_date_from_filename <- function(fname) {
  base <- sub("\\.rda$", "", basename(fname))
  as.POSIXct(strptime(paste0("20", substr(base,4,5), " ",
                              substr(base,6,8), " ", substr(base,9,10)),
                      format = "%Y %j %H", tz = "UTC"))
}
dates <- as.POSIXct(sapply(rda_files, parse_date_from_filename),
                    origin = "1970-01-01", tz = "UTC")

# ═══════════════════════════════════════════════════════════════════════════════
# SHEET 1 — Seasonal activation summary
# ═══════════════════════════════════════════════════════════════════════════════

coeff_matrix <- do.call(cbind, lapply(all_data, function(df) df$coef.var))

activations_ts <- data.frame(
  date     = dates,
  count    = as.numeric(colSums(coeff_matrix > threshold, na.rm = TRUE)),
  fraction = as.numeric(colMeans(coeff_matrix > threshold, na.rm = TRUE))
)
activations_ts$month  <- as.integer(format(dates, "%m"))
activations_ts$season <- cut(activations_ts$month,
                             breaks = c(0, 2, 5, 8, 11, 12),
                             labels = c("DJF","MAM","JJA","SON","DJF2"), right = TRUE)
levels(activations_ts$season)[levels(activations_ts$season) == "DJF2"] <- "DJF"

seasonal_summary <- do.call(rbind, lapply(c("DJF","MAM","JJA","SON"), function(s) {
  sub <- activations_ts[activations_ts$season == s, ]
  data.frame(
    Season        = s,
    N_Timesteps   = nrow(sub),
    Mean_Count    = round(mean(sub$count),                3),
    Mean_Fraction = round(mean(sub$fraction, na.rm = TRUE), 4),
    Max_Count     = max(sub$count),
    N_Active_Ts   = sum(sub$count > 0),
    Pct_Active_Ts = round(mean(sub$count > 0) * 100,    1)
  )
}))

# ═══════════════════════════════════════════════════════════════════════════════
# SHEET 2 — Correction statistics (radar − radar.orig at active cells, CC ON)
# ═══════════════════════════════════════════════════════════════════════════════

corr_rows <- do.call(rbind, lapply(seq_along(all_data), function(i) {
  df     <- all_data[[i]]
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
    nrow(corr_rows),
    mean(corr_rows$mean_corr),
    mean(corr_rows$abs_mean_corr),
    min(corr_rows$min_corr),
    max(corr_rows$max_corr),
    weighted.mean(corr_rows$mean_corr,     corr_rows$n_active),
    weighted.mean(corr_rows$abs_mean_corr, corr_rows$n_active),
    weighted.mean(corr_rows$mean_rel_corr, corr_rows$n_active, na.rm = TRUE),
    weighted.mean(corr_rows$abs_rel_corr,  corr_rows$n_active, na.rm = TRUE)
  ), 4),
  Notes = c(
    paste0(nrow(corr_rows), " / ", length(all_data),
           " (", round(nrow(corr_rows)/length(all_data)*100, 1), "% of year)"),
    "negative = radar decreased", "magnitude regardless of sign",
    "single grid cell, single timestep", "single grid cell, single timestep",
    "weighted by n_active per timestep", "weighted magnitude",
    "relative to radar.orig, cells with orig > 0.1", "weighted magnitude"
  )
)

# ═══════════════════════════════════════════════════════════════════════════════
# SHEET 3 — Direct ON vs OFF: delta radar and delta CoV
# ═══════════════════════════════════════════════════════════════════════════════

dates_on  <- names(all_data)
dates_off <- names(all_data_off)
common_on  <- which(dates_on  %in% dates_off)
common_off <- which(dates_off %in% dates_on)

compare_long <- rbindlist(lapply(seq_along(common_on), function(i) {
  df_on  <- all_data[[common_on[i]]]
  df_off <- all_data_off[[common_off[i]]]
  data.table(station_id = seq_len(nrow(df_on)),
             time       = dates_on[common_on[i]],
             radar_on   = df_on$radar,
             radar_off  = df_off$radar,
             cv_on      = df_on$coef.var,
             cv_off     = df_off$coef.var)
}), fill = TRUE)

compare_long[, delta_radar := radar_on - radar_off]
compare_long[, delta_cv    := cv_on    - cv_off   ]

wet          <- compare_long[!is.na(radar_on) & !is.na(radar_off) & (radar_on > 0 | radar_off > 0)]
active_cells <- wet[cv_on > threshold]

fmt6 <- function(x) round(c(Mean=mean(x,na.rm=T), Median=median(x,na.rm=T),
                             Q1=quantile(x,.25,na.rm=T), Q3=quantile(x,.75,na.rm=T),
                             SD=sd(x,na.rm=T), N=sum(!is.na(x))), 4)

on_vs_off <- rbind(
  data.frame(Subset="All wet cells",          Variable="delta_radar (ON-OFF, mm/h)", t(fmt6(wet$delta_radar))),
  data.frame(Subset="All wet cells",          Variable="delta_CoV   (ON-OFF)",       t(fmt6(wet$delta_cv))),
  data.frame(Subset="Active cells (CoV>0.5)", Variable="delta_radar (ON-OFF, mm/h)", t(fmt6(active_cells$delta_radar))),
  data.frame(Subset="Active cells (CoV>0.5)", Variable="delta_CoV   (ON-OFF)",       t(fmt6(active_cells$delta_cv)))
)

# ═══════════════════════════════════════════════════════════════════════════════
# SHEET 4 — μ, IQR, IQR/μ: active vs inactive
# ═══════════════════════════════════════════════════════════════════════════════

coefvar_long <- rbindlist(lapply(seq_along(all_data), function(i) {
  df <- all_data[[i]]
  data.table(station_id = seq_len(nrow(df)), time = dates[i], coef.var = df$coef.var)
}))

result_full        <- merge(result, coefvar_long, by = c("station_id", "time"))
result_full$active <- result_full$coef.var > threshold & !is.na(result_full$coef.var)
result_full[, rel_unc := iqr / mu]

mu_iqr_summary <- rbind(
  data.frame(Variable="mu (mm/h)",  Group="Active",   t(fmt6(result_full$mu[ result_full$active]))),
  data.frame(Variable="mu (mm/h)",  Group="Inactive", t(fmt6(result_full$mu[!result_full$active]))),
  data.frame(Variable="IQR (mm/h)", Group="Active",   t(fmt6(result_full$iqr[ result_full$active]))),
  data.frame(Variable="IQR (mm/h)", Group="Inactive", t(fmt6(result_full$iqr[!result_full$active]))),
  data.frame(Variable="IQR/mu",     Group="Active",   t(fmt6(result_full$rel_unc[ result_full$active]))),
  data.frame(Variable="IQR/mu",     Group="Inactive", t(fmt6(result_full$rel_unc[!result_full$active])))
)

wx <- wilcox.test(result_full$rel_unc[ result_full$active],
                  result_full$rel_unc[!result_full$active], conf.int = TRUE)
wilcox_result <- data.frame(
  Test        = "Wilcoxon rank-sum: IQR/mu active vs inactive",
  W_statistic = round(wx$statistic, 1),
  p_value     = formatC(wx$p.value, format = "e", digits = 3),
  Estimate    = round(wx$estimate,    4),
  CI_low      = round(wx$conf.int[1], 4),
  CI_high     = round(wx$conf.int[2], 4)
)

# ═══════════════════════════════════════════════════════════════════════════════
# SHEET 5 — Binned intensity: IQR/mu by decile at top 10% of μ
# ═══════════════════════════════════════════════════════════════════════════════

mu_p90    <- quantile(result_full$mu, 0.90, na.rm = TRUE)
intense   <- result_full[!is.na(mu) & mu >= mu_p90]
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
  Median_relunc_inactive = round(median(rel_unc[active == FALSE], na.rm = TRUE), 4),
  Wilcox_p               = round(tryCatch({
    na  <- sum(active == TRUE  & !is.na(rel_unc))
    ni  <- sum(active == FALSE & !is.na(rel_unc))
    if (na < 5 || ni < 5) NA_real_   # too few obs for reliable test
    else wilcox.test(rel_unc[active == TRUE], rel_unc[active == FALSE])$p.value
  }, error = function(e) NA_real_), 4)
), by = mu_bin][order(mu_bin)]

bin_comparison[, Delta_relunc   := round(Median_relunc_active - Median_relunc_inactive, 4)]
bin_comparison[, Pct_diff_relunc := round((Median_relunc_active - Median_relunc_inactive) /
                                           Median_relunc_inactive * 100, 2)]

# ═══════════════════════════════════════════════════════════════════════════════
# SHEET 6 — Log-log regression: IQR ~ mu elasticity, active vs inactive
# ═══════════════════════════════════════════════════════════════════════════════

station_comparison <- result_full[!is.na(mu), .(
  x                 = x_station[1],
  y                 = y_station[1],
  n_active          = sum(active == TRUE,  na.rm = TRUE),
  n_inactive        = sum(active == FALSE, na.rm = TRUE),
  mean_mu_active    = mean(mu[active == TRUE],    na.rm = TRUE),
  mean_mu_inactive  = mean(mu[active == FALSE],   na.rm = TRUE),
  mean_iqr_active   = mean(iqr[active == TRUE],   na.rm = TRUE),
  mean_iqr_inactive = mean(iqr[active == FALSE],  na.rm = TRUE)
), by = station_id]

reg_data <- station_comparison[!is.na(mean_mu_active) & !is.na(mean_iqr_active) & n_active >= 10]
reg_data[, log_mu_active    := log(mean_mu_active)]
reg_data[, log_iqr_active   := log(mean_iqr_active)]
reg_data[, log_mu_inactive  := log(mean_mu_inactive)]
reg_data[, log_iqr_inactive := log(mean_iqr_inactive)]

reg_combined <- data.table(
  log_mu  = c(reg_data$log_mu_active,  reg_data$log_mu_inactive),
  log_iqr = c(reg_data$log_iqr_active, reg_data$log_iqr_inactive),
  active  = c(rep(1L, nrow(reg_data)), rep(0L, nrow(reg_data)))
)

fit_combined <- lm(log_iqr ~ log_mu * active, data = reg_combined)
cf <- coef(fit_combined)
su <- summary(fit_combined)$coefficients

elasticity_inactive <- cf["log_mu"]
elasticity_active   <- cf["log_mu"] + cf["log_mu:active"]

loglog_results <- data.frame(
  Parameter = c(
    "Intercept (inactive)", "Elasticity ε (inactive)",
    "Intercept shift (active − inactive)",
    "Elasticity ε (active)",
    "Interaction p-value (log_mu:active)",
    "R-squared"
  ),
  Estimate = c(
    round(cf["(Intercept)"],  4),
    round(elasticity_inactive, 4),
    round(cf["active"],        4),
    round(elasticity_active,   4),
    formatC(su["log_mu:active", "Pr(>|t|)"], format = "e", digits = 3),
    round(summary(fit_combined)$r.squared, 4)
  ),
  Interpretation = c(
    "Baseline IQR at log(mu)=0, inactive periods",
    "IQR scales as mu^ε during inactive periods (ε<1: IQR grows slower than mu)",
    "Negative = IQR systematically lower during active periods (log space)",
    "IQR elasticity during convective control activation",
    "Significance of active vs inactive slope difference",
    "Variance in log(IQR) explained by log(mu) and active flag"
  )
)

spearman_r <- cor(station_comparison$mean_iqr_active   / station_comparison$mean_mu_active,
                  station_comparison$mean_iqr_inactive  / station_comparison$mean_mu_inactive,
                  use = "complete.obs", method = "spearman")

spearman_df <- data.frame(
  Metric = "Spearman r: per-station IQR/mu (active vs inactive)",
  Value  = round(spearman_r, 4),
  Interpretation = "Spatial consistency of relative uncertainty between regimes"
)

# ═══════════════════════════════════════════════════════════════════════════════
# SHEET 7 — Cross-regime: IQR/mu at active pairs, ON vs OFF
# For gauge-timestamp pairs where CC fired (active in ON), compare IQR/mu
# between the ON run and the same pairs in the OFF run.
# ═══════════════════════════════════════════════════════════════════════════════

setDT(result_off)
result_off_sub <- result_off[, .(station_id, time, mu_off = mu, iqr_off = iqr)]

cross_regime <- merge(
  result_full[active == TRUE,
              .(station_id, time, mu_on = mu, iqr_on = iqr, rel_unc_on = rel_unc)],
  result_off_sub,
  by = c("station_id", "time")
)
cross_regime[, rel_unc_off  := iqr_off / mu_off]
cross_regime[, delta_relunc := rel_unc_on - rel_unc_off]

# Complete pairs only — so ON, OFF, and delta rows all reflect the same N
complete_pairs <- cross_regime[!is.na(rel_unc_on) & !is.na(rel_unc_off)]

cross_summary <- rbind(
  data.frame(Variable = "IQR/mu  — ON",        t(fmt6(complete_pairs$rel_unc_on))),
  data.frame(Variable = "IQR/mu  — OFF",       t(fmt6(complete_pairs$rel_unc_off))),
  data.frame(Variable = "IQR/mu  ON − OFF",    t(fmt6(complete_pairs$delta_relunc))),
  data.frame(Variable = "mu ON   (mm/h)",      t(fmt6(complete_pairs$mu_on))),
  data.frame(Variable = "mu OFF  (mm/h)",      t(fmt6(complete_pairs$mu_off))),
  data.frame(Variable = "IQR ON  (mm/h)",      t(fmt6(complete_pairs$iqr_on))),
  data.frame(Variable = "IQR OFF (mm/h)",      t(fmt6(complete_pairs$iqr_off)))
)

wx7 <- wilcox.test(complete_pairs$rel_unc_on, complete_pairs$rel_unc_off,
                   paired = TRUE, conf.int = TRUE)
wilcox_cross <- data.frame(
  Test        = "Wilcoxon signed-rank: IQR/mu ON vs OFF at active pairs (paired)",
  W_statistic = round(wx7$statistic, 1),
  p_value     = formatC(wx7$p.value, format = "e", digits = 3),
  Estimate    = round(wx7$estimate,    4),
  CI_low      = round(wx7$conf.int[1], 4),
  CI_high     = round(wx7$conf.int[2], 4),
  N_pairs     = nrow(complete_pairs)
)

# ═══════════════════════════════════════════════════════════════════════════════
# WRITE EXCEL
# ═══════════════════════════════════════════════════════════════════════════════

out_xlsx <- file.path(project_root, paste0("conv_control_stats_", yr_label, ".xlsx"))
wb <- loadWorkbook(out_xlsx, create = TRUE)

write_section <- function(wb, sheet, header, df) {
  createSheet(wb, name = sheet)
  writeWorksheet(wb, data.frame(V1 = header), sheet = sheet,
                 startRow = 1, startCol = 1, header = FALSE)
  writeWorksheet(wb, df, sheet = sheet, startRow = 3, startCol = 1)
}

write_section(wb, "1_Activation_Seasonal",
  paste("Seasonal activation summary | CoV threshold >", threshold, "|", yr_label),
  seasonal_summary)

write_section(wb, "2_Correction_Stats",
  paste("Correction statistics (radar − radar.orig at active cells) |", yr_label),
  correction_stats)

write_section(wb, "3_ON_vs_OFF",
  paste("Direct ON vs OFF | delta radar and delta CoV | all wet cells and active cells |", yr_label),
  on_vs_off)

write_section(wb, "4_mu_IQR_summary",
  paste("mu, IQR, IQR/mu: active vs inactive timesteps |", yr_label),
  mu_iqr_summary)

write_section(wb, "4b_Wilcoxon",
  "Wilcoxon rank-sum test: IQR/mu active vs inactive",
  wilcox_result)

write_section(wb, "5_Intensity_bins",
  paste("IQR/mu by intensity decile (top 10% of mu) | active vs inactive |", yr_label),
  as.data.frame(bin_comparison))

# write_section(wb, "6_LogLog_elasticity",
#   paste("Log-log regression: log(IQR) ~ log(mu) * active |", yr_label),
#   loglog_results)

# write_section(wb, "6b_Spearman",
#   "Spearman correlation: per-station IQR/mu between active and inactive regimes",
#   spearman_df)

write_section(wb, "7_CrossRegime_active",
  paste("Cross-regime: IQR/mu at active (CC-fired) pairs — ON vs OFF |", yr_label),
  cross_summary)

write_section(wb, "7b_CrossRegime_Wilcoxon",
  "Wilcoxon signed-rank (paired): IQR/mu ON vs OFF at active pairs",
  wilcox_cross)

saveWorkbook(wb)
message("Saved: ", out_xlsx)
