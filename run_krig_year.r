#!/usr/bin/env Rscript

### load libs
.libPaths("/store_new/mch/msclim/share/CATs/cats/lib-R4.4.0/")
library(raster)
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
library(gstat)                                                      #but I had to install sf first
library(rhdf5)
library(animation)
library(ncdf4)                     
library(lubridate)
library(dplyr)

#year <- 2023
# -----------------------------
# CONFIG
# -----------------------------

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Provide a year, e.g. Rscript run_year_thresholds.R 2021")

year <- as.character(args[1])

# "all"    = run  existing plots
# "relunc" = compute/plot only relative uncertainty (IQR/μ)
mode <- if (length(args) >= 2) tolower(args[2]) else "all"
if (!mode %in% c("all", "relunc")) stop("mode must be 'all' or 'relunc'")

# minimum μ to avoid exploding ratios
mu_min <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
if (!is.finite(mu_min) || mu_min < 0) stop("mu_min must be a non-negative number")


project_root <- "/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats"
data_root    <- "/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project"
lib_root     <- "/store_new/mch/msclim/antoumos/R/lib"

thresholds <- c(0.1, 0.5, 1, 2)
xlim <- c(480, 840) # Swiss crop window
ylim <- c(60, 300)

base_out_dir <- file.path(project_root, "out_plots")
utils_file   <- file.path(project_root, "R", "utils.r")

rda_file <- file.path(
  data_root,
  sprintf("precip_transformed_results_conv_off_new_%s.rda", year)
)

rda_file_base <- file.path(
  data_root,
  sprintf("precip_transformed_results_new_%s.rda", year)
)

message("Using input file: ", rda_file)

out_dir <- file.path(base_out_dir, paste0("year_", year))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


.libPaths(lib_root)
plot_utils_file <- file.path(project_root, "R", "plot_utils.r")
source(utils_file)      # helpers function
source(plot_utils_file) # plotting functions

suppressPackageStartupMessages({
  library(matrixStats)
  library(abind)
})

#-------------------------------
# Compute stats for one threshold
# -----------------------------

# compute_year_threshold_stats_simple <- function(rda_file, min_threshold, max_threshold = Inf, mu_min = 0.1, rel_uncert_method = c("A", "B_median","B_IQR")) {

#   rel_uncert_method <- match.arg(rel_uncert_method)
#   load(rda_file)  # expects: kriging_crop_list, variance_crop_list, timestamps

#   dates <- as.POSIXct(timestamps, origin = "1970-01-01", tz = "UTC")
#   season <- get_season(dates)

#   # Keep a reference list WITH x/y for plotting
#   variance_ref_list <- variance_crop_list

#   # Strip heavy attrs for compute (keeps dim only, per your utility)
#   variance_raw_list <- lapply(variance_crop_list, strip_attrs_keep_dim)

#   kriging_thr_list <- lapply(kriging_crop_list, function(m) {
#     m[m < min_threshold] <- NA_real_
#     strip_attrs_keep_dim(m)
#   })

#   # Annual fields
#   annual_mean_mu   <- aggregate_mean(kriging_thr_list)$mean
#   annual_accum_mu  <- aggregate_sum(kriging_thr_list)$sum  # hourly sums = mm
#   annual_wet_hours <- aggregate_wet_hours(kriging_thr_list, min_threshold, max_threshold)

#   iqr_list <- compute_iqr_list(kriging_thr_list, variance_raw_list)
#   annual_mean_iqr <- aggregate_mean(iqr_list)$mean
  
#   ### Relative uncertainty ##

#   # Option A (ratio of annual means)
#   annual_rel_A <- annual_mean_iqr / annual_mean_mu
#   annual_rel_A[!is.finite(annual_rel_A)] <- NA_real_
#   annual_rel_A[annual_mean_mu < mu_min] <- NA_real_
  
#   # Option B (median of hourly ratios)
#   annual_rel_Bmed <- aggregate_rel_uncert_B(
#     mu_list   = kriging_thr_list,
#     iqr_list  = iqr_list,
#     mu_min    = mu_min,
#     method    = "median",
#     min_hours = 50
#   )

#   annual_rel_B_iqr <- aggregate_rel_uncert_iqr(
#   mu_list   = kriging_thr_list,
#   iqr_list  = iqr_list,
#   mu_min    = mu_min,
#   min_hours = 10
#   )
# # Seasonal fields
#   seasons  <- levels(season)
#   seasonal <- setNames(vector("list", length(seasons)), seasons)

#   for (s in seasons) {
#     idx <- which(season == s)

#     s_mu_list  <- kriging_thr_list[idx]
#     s_iqr_list <- iqr_list[idx]

#     s_mean_mu   <- aggregate_mean(s_mu_list)$mean
#     s_accum_mu  <- aggregate_sum(s_mu_list)$sum
#     s_wet_hours <- aggregate_wet_hours(s_mu_list, min_threshold, max_threshold)
#     s_mean_iqr  <- aggregate_mean(s_iqr_list)$mean
    
#  #    Option A
#   s_rel_A <- s_mean_iqr / s_mean_mu
#   s_rel_A[!is.finite(s_rel_A)] <- NA_real_
#   s_rel_A[s_mean_mu < mu_min] <- NA_real_

# # Option B median
#   s_rel_Bmed <- aggregate_rel_uncert_B(
#    mu_list   = s_mu_list,
#     iqr_list  = s_iqr_list,
#     mu_min    = mu_min,
#     method    = "median",
#     min_hours = 20
#   )
# # Option B IQR of hourly relative uncertainty
#   s_rel_B_iqr <- aggregate_rel_uncert_iqr(
#     mu_list   = s_mu_list,
#     iqr_list  = s_iqr_list,
#     mu_min    = mu_min,
#     min_hours = 10
#   )
    
#     seasonal[[s]] <- list(
#       mean_mu    = s_mean_mu,
#       accum_mu   = s_accum_mu,
#       wet_hours  = s_wet_hours,
#       mean_iqr   = s_mean_iqr,
#       rel_uncert_A = s_rel_A,
#       rel_uncert_B_median = s_rel_Bmed,
#       iqr_relative_uncertainty = s_rel_B_iqr
#     )
#   }

#   list(
#     min_threshold = min_threshold,
#     max_threshold = max_threshold,
#     rel_uncert_method = rel_uncert_method,
#     annual = list(
#       mean_mu    = annual_mean_mu,
#       mean_iqr   = annual_mean_iqr,
#       rel_uncert_A = annual_rel_A,
#       rel_uncert_B_median  = annual_rel_Bmed,
#       iqr_relative_uncertainty = annual_rel_B_iqr,
#       accum_mu   = annual_accum_mu,
#       wet_hours  = annual_wet_hours
#     ),
#     seasonal = seasonal,
#     variance_ref_list = variance_ref_list
#   )
# }



# -----------------------------
# Plotting (annual + seasonal)
# -----------------------------
plot_stats_simple <- function(res, year_label, out_dir, xlim, ylim,
                              mode = c("all", "relunc"),
                              cap_quant = 0.99, palette_end = 0.99,
                              cap_quant_relunc = 0.98) {

  mode <- match.arg(mode)
  thr_txt <- gsub("\\.", "p", sprintf("%.2f", res$min_threshold))
 
  # helper to avoid repeating plot calls
  plot_one <- function(Z, fname, ttl, capq = cap_quant, pal_end = palette_end) {
    if (is.null(Z) || length(dim(Z)) != 2) return(invisible(FALSE))
    plot_cropped_field(
      Z = Z,
      variance_list = res$variance_ref_list,
      xlim = xlim, ylim = ylim,
      title = ttl,
      cap_quant = capq,
      palette_end = pal_end,
      output_file = file.path(out_dir, fname)
    )
    invisible(TRUE)
  }

  plot_relunc_variants <- function(obj, season_name = NULL) {
  # obj is either res$annual or res$seasonal[[s]]
  variants <- c(
    #A         = "rel_uncert_A",
    #B_median  = "rel_uncert_B_median",
    B_IQR     = "iqr_relative_uncertainty"
  )

  for (lab in names(variants)) {
    key <- variants[[lab]]
    if (is.null(obj[[key]])) {
      warning("plot_relunc_variants: '", key, "' is NULL, skipping plot")
      next
    }

    s_tag <- if (!is.null(season_name)) paste0("SEASON_", season_name, "_") else "ANNUAL_"
    ttl_season <- if (!is.null(season_name)) paste0(" ", season_name) else " annual"
     
    ttl <- if (lab == "B_IQR") {
      sprintf("%s%s IQR of relative uncertainty, thr=%.2f",
              year_label, ttl_season, res$min_threshold)
    }
    else {
      sprintf("%s%s relative uncertainty (%s) (IQR/\u03bc), thr=%.2f",
              year_label, ttl_season, lab, res$min_threshold)
    }
  
  fname <- sprintf("%sRELUNC_%s_%s_thr_%s.png", s_tag, lab, year_label, thr_txt)

    plot_one(
      obj[[key]],
      fname,
      ttl,
      capq = cap_quant_relunc
    )
  }
}
#Annual
  plot_relunc_variants(res$annual, season_name = NULL)

  # Seasonal
  for (s in names(res$seasonal)) {
    plot_relunc_variants(res$seasonal[[s]], season_name = s)
  }

  # If user asked ONLY relunc, stop here
  if (mode == "relunc") return(invisible(TRUE))

  # -----------------------------
  # Existing plots: only when mode == "all"
  # -----------------------------

  # Annual
  plot_one(res$annual$mean_mu,
           sprintf("ANNUAL_MU_%s_thr_%s.png", year_label, thr_txt),
           sprintf("%s annual kriging mean (μ), thr=%.2f", year_label, res$min_threshold))

  plot_one(res$annual$mean_iqr,
           sprintf("ANNUAL_IQR_%s_thr_%s.png", year_label, thr_txt),
           sprintf("%s annual IQR(90-10), thr=%.2f", year_label, res$min_threshold))

  plot_one(res$annual$accum_mu,
           sprintf("ANNUAL_ACCUM_%s_thr_%s.png", year_label, thr_txt),
           sprintf("%s annual accumulation (mm), thr=%.2f", year_label, res$min_threshold))

  plot_one(res$annual$wet_hours,
           sprintf("ANNUAL_WET_HOURS_%s_thr_%s.png", year_label, thr_txt),
           sprintf("%s annual wet hours (>%.2f)", year_label, res$min_threshold))

  # Seasonal
  for (s in names(res$seasonal)) {

    plot_one(res$seasonal[[s]]$mean_mu,
             sprintf("SEASON_%s_MU_%s_thr_%s.png", s, year_label, thr_txt),
             sprintf("%s %s kriging mean (μ), thr=%.2f", year_label, s, res$min_threshold))

    plot_one(res$seasonal[[s]]$mean_iqr,
             sprintf("SEASON_%s_IQR_%s_thr_%s.png", s, year_label, thr_txt),
             sprintf("%s %s IQR(90-10), thr=%.2f", year_label, s, res$min_threshold))

    plot_one(res$seasonal[[s]]$accum_mu,
             sprintf("SEASON_%s_ACCUM_%s_thr_%s.png", s, year_label, thr_txt),
             sprintf("%s %s accumulation (mm), thr=%.2f", year_label, s, res$min_threshold))
    plot_one(res$seasonal[[s]]$wet_hours,
            sprintf("SEASON_%s_WET_HOURS_%s_thr_%s.png", s, year_label, thr_txt),
            sprintf("%s %s wet hours (>%.2f)", year_label, s, res$min_threshold) )
  }

  invisible(TRUE)
}

# -----------------------------
# Difference plots (conv_off minus base)
# -----------------------------
plot_diff_stats <- function(res_off, res_base, year_label, out_dir, xlim, ylim,
                             cap_quant = 0.99) {
  thr_txt  <- gsub("\\.", "p", sprintf("%.2f", res_off$min_threshold))
  vref     <- res_off$variance_ref_list
  diff_pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(150)

  plot_diff <- function(Z_off, Z_base, fname, ttl) {
    if (is.null(Z_off) || is.null(Z_base)) return(invisible(FALSE))
    dZ   <- Z_off - Z_base
    vmax <- quantile(abs(dZ), cap_quant, na.rm = TRUE)
    vmax <- max(vmax, 1e-9)
    plot_cropped_field(
      Z             = dZ,
      variance_list = vref,
      xlim = xlim, ylim = ylim,
      vmin = -vmax, vmax = vmax,
      pal  = diff_pal,
      title = ttl,
      output_file = file.path(out_dir, fname)
    )
    invisible(TRUE)
  }

  fields <- list(
    list(key = "mean_mu",             label = "mean μ",                 unit = "mm/h"),
    list(key = "mean_iqr",            label = "mean IQR(90-10)",             unit = "mm/h"),
    list(key = "accum_mu",            label = "accumulation",                unit = "mm"),
    list(key = "wet_hours",           label = "wet hours",                   unit = "h"),
    list(key = "rel_uncert_B_median", label = "rel. uncertainty (B_median)", unit = "")
  )

  for (f in fields) {
    unit_str <- if (nchar(f$unit) > 0) paste0(" (", f$unit, ")") else ""
    plot_diff(
      res_off$annual[[f$key]], res_base$annual[[f$key]],
      fname = sprintf("DIFF_ANNUAL_%s_%s_thr_%s.png", toupper(f$key), year_label, thr_txt),
      ttl   = sprintf("%s annual DIFF conv_off-base %s%s, thr=%.2f",
                      year_label, f$label, unit_str, res_off$min_threshold)
    )
  }

  for (s in names(res_off$seasonal)) {
    if (!s %in% names(res_base$seasonal)) next
    for (f in fields) {
      unit_str <- if (nchar(f$unit) > 0) paste0(" (", f$unit, ")") else ""
      plot_diff(
        res_off$seasonal[[s]][[f$key]], res_base$seasonal[[s]][[f$key]],
        fname = sprintf("DIFF_SEASON_%s_%s_%s_thr_%s.png", s, toupper(f$key), year_label, thr_txt),
        ttl   = sprintf("%s %s DIFF conv_off-base %s%s, thr=%.2f",
                        year_label, s, f$label, unit_str, res_off$min_threshold)
      )
    }
  }

  invisible(TRUE)
}

# Interactive run config
# -----------------------------

# Plot / compute options
mode   <- "all"     # "all" or "relunc"
mu_min <- 0.1       # mask for tiny mu in relunc

thresholds <- c(0.1)  # can be c(0.1, 0.5, 1.0) etc


# # -----------------------------
# # RUN
# # -----------------------------
stopifnot(file.exists(rda_file))

out_dir <- ("/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats/")

message("Running year=", year, " file=", rda_file)
message("Mode=", mode, " | mu_min=", mu_min)
message("Plots -> ", normalizePath(out_dir, mustWork = FALSE))

for (thr in thresholds) {
  message("Threshold: ", thr)

  res_off <- compute_year_threshold_stats_simple(
    rda_file = rda_file,
    min_threshold = thr,
    mu_min = mu_min,
    rel_uncert_method = "B_IQR"
  )

  #plot_stats_simple(
  #  res = res_off,
  #  year_label = year,
  #  out_dir = out_dir,
#xlim = xlim,
  #  ylim = ylim,
  #  mode = mode
  #)

  if (file.exists(rda_file_base)) {
    message("Computing diff with base file: ", rda_file_base)

    res_base <- compute_year_threshold_stats_simple(
      rda_file = rda_file_base,
      min_threshold = thr,
      mu_min = mu_min,
      rel_uncert_method = "B_IQR"
    )

    plot_diff_stats(
      res_off    = res_off,
      res_base   = res_base,
      year_label = year,
      out_dir    = out_dir,
      xlim       = xlim,
      ylim       = ylim
    )
  } else {
    message("Base rda not found, skipping diff plots: ", rda_file_base)
  }
}

message("Done.")





