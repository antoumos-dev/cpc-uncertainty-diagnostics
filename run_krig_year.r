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

# Optional 2nd argument: mode
# "all"    = run your existing plots
# "relunc" = compute/plot only relative uncertainty (IQR/μ)
mode <- if (length(args) >= 2) tolower(args[2]) else "all"
if (!mode %in% c("all", "relunc")) stop("mode must be 'all' or 'relunc'")

# Optional 3rd argument: minimum μ to avoid exploding ratios (default sensible-ish)
mu_min <- if (length(args) >= 3) as.numeric(args[3]) else 0.05
if (!is.finite(mu_min) || mu_min < 0) stop("mu_min must be a non-negative number")

# Input file
rda_file <- sprintf(
  "/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/precip_transformed_results_new_%s.rda",
  year
)

thresholds <- c(0.1, 0.5, 1, 2)

# Plot domain (your Swiss crop window)
xlim <- c(480, 840)
ylim <- c(60, 300)

base_out_dir <- "/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats/out_plots"
# Output folder for plots (separate folder per run)
out_dir <- file.path(base_out_dir,
                     paste0("year_", year))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

.libPaths("/store_new/mch/msclim/antoumos/R/lib")

suppressPackageStartupMessages({
  library(matrixStats)
  library(abind)
})

# Load helper functions
source("/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats/R/utils.r")
#-------------------------------
# Compute stats for one threshold
# -----------------------------

compute_year_threshold_stats_simple <- function(rda_file, min_threshold, max_threshold , mu_min = 0.1, rel_uncert_method = c("A", "B_median","B_IQR")) {

  rel_uncert_method <- match.arg(rel_uncert_method) ## we have 2 methods to compute rel uncertainty

  load(rda_file)  # expects: kriging_crop_list, variance_crop_list, timestamps

  dates <- as.POSIXct(timestamps, origin = "1970-01-01", tz = "UTC")
  season <- get_season(dates)

  # Strip heavy attrs for compute (keeps dim only, per your utility)
  variance_raw_list <- lapply(variance_crop_list, strip_attrs_keep_dim)

  kriging_thr_list <- lapply(kriging_crop_list, function(m) {
    m[m < min_threshold] <- NA_real_
    strip_attrs_keep_dim(m)
  })

  # Annual fields
  annual_mean_mu   <- aggregate_mean(kriging_thr_list)$mean
  annual_accum_mu  <- aggregate_sum(kriging_thr_list)$sum  
  annual_wet_hours <- aggregate_wet_hours(kriging_thr_list, min_threshold, max_threshold)

  iqr_list <- compute_iqr_list(kriging_thr_list, variance_raw_list)
  annual_mean_iqr <- aggregate_mean(iqr_list)$mean
  
  ### Relative uncertainty ##

  # Option A (ratio of annual means)
  annual_rel_A <- annual_mean_iqr / annual_mean_mu
  annual_rel_A[!is.finite(annual_rel_A)] <- NA_real_
  annual_rel_A[annual_mean_mu < mu_min] <- NA_real_
  
  # Option B (median of hourly ratios)
  annual_rel_Bmed <- aggregate_rel_uncert_B(
    mu_list   = kriging_thr_list,
    iqr_list  = iqr_list,
    mu_min    = mu_min,
    method    = "median",
    min_hours = 50
  )

  annual_rel_B_iqr <- aggregate_rel_uncert_iqr(
  mu_list   = kriging_thr_list,
  iqr_list  = iqr_list,
  mu_min    = mu_min,
  min_hours = 10
  )
# Seasonal fields
  seasons  <- levels(season)
  seasonal <- setNames(vector("list", length(seasons)), seasons)

  for (s in seasons) {
    idx <- which(season == s)

    s_mu_list  <- kriging_thr_list[idx]
    s_iqr_list <- iqr_list[idx]

    s_mean_mu   <- aggregate_mean(s_mu_list)$mean
    s_accum_mu  <- aggregate_sum(s_mu_list)$sum
    s_wet_hours <- aggregate_wet_hours(s_mu_list, min_threshold, max_threshold)
    s_mean_iqr  <- aggregate_mean(s_iqr_list)$mean
    
 #    Option A
  s_rel_A <- s_mean_iqr / s_mean_mu
  s_rel_A[!is.finite(s_rel_A)] <- NA_real_
  s_rel_A[s_mean_mu < mu_min] <- NA_real_

# Option B median
  s_rel_Bmed <- aggregate_rel_uncert_B(
   mu_list   = s_mu_list,
    iqr_list  = s_iqr_list,
    mu_min    = mu_min,
    method    = "median",
    min_hours = 20
  )
    
    seasonal[[s]] <- list(
      mean_mu    = s_mean_mu,
      accum_mu   = s_accum_mu,
      wet_hours  = s_wet_hours,
      mean_iqr   = s_mean_iqr,
      rel_uncert_A = s_rel_A,
      rel_uncert_B_median = s_rel_Bmed
    )
  }

  list(
    min_threshold = min_threshold,
    max_threshold = max_threshold,
    rel_uncert_method = rel_uncert_method,
    annual = list(
      mean_mu    = annual_mean_mu,
      mean_iqr   = annual_mean_iqr,
      rel_uncert_A = annual_rel_A,
      rel_uncert_B_median  = annual_rel_Bmed,
      iqr_relative_uncertainty = annual_rel_B_iqr,
      accum_mu   = annual_accum_mu,
      wet_hours  = annual_wet_hours
    ),
    seasonal = seasonal,
    variance_ref_list = variance_crop_list
  )
}


# -----------------------------
# Plotting (annual + seasonal)
# -----------------------------
plot_stats_simple <- function(res, year_label, out_dir, xlim, ylim,
                              mode = c("all", "relunc"),
                              cap_quant = 0.99, palette_end = 0.99,
                              cap_quant_relunc = 0.98) {

  mode <- match.arg(mode)
 # thr_txt <- gsub("\\.", "p", sprintf("%.2f", res$threshold))
min_thr <- res$min_threshold
max_thr <- res$max_threshold

min_ok <- length(min_thr) == 1 && is.finite(min_thr)
max_ok <- length(max_thr) == 1 && is.finite(max_thr)

thr_label <- if (min_ok && max_ok) {
  sprintf("%.2f-%.2f", min_thr, max_thr)
} else if (min_ok) {
  sprintf(">%.2f", min_thr)
} else if (max_ok) {
  sprintf("<=%.2f", max_thr)
} else {
  "all"
}

thr_txt <- if (min_ok && max_ok) {
  sprintf("%s_to_%s",
          gsub("\\.", "p", sprintf("%.2f", min_thr)),
          gsub("\\.", "p", sprintf("%.2f", max_thr)))
} else if (min_ok) {
  sprintf("gt_%s", gsub("\\.", "p", sprintf("%.2f", min_thr)))
} else if (max_ok) {
  sprintf("lte_%s", gsub("\\.", "p", sprintf("%.2f", max_thr)))
} else {
  "all"
}

  # helper to avoid repeating plot calls
  plot_one <- function(Z, fname, ttl, capq = cap_quant, pal_end = palette_end) {
    if (is.null(Z) || length(dim(Z)) != 2) return(invisible(FALSE))
    plot_cropped_field(
      Z = Z,
      variance_list = res$variance_ref_list,
      xlim = xlim, ylim = ylim,
      main_title = ttl,
      cap_quant = capq,
      palette_end = pal_end,
      output_file = file.path(out_dir, fname)
    )
    invisible(TRUE)
  }

  plot_relunc_variants <- function(obj, season_name = NULL) {
  # obj is either res$annual or res$seasonal[[s]]
  variants <- c(
    A         = "rel_uncert_A",
    B_median  = "rel_uncert_B_median",
    B_IQR     = "iqr_relative_uncertainty"
    # you can add B_weighted etc later
  )

  for (lab in names(variants)) {
    key <- variants[[lab]]
    if (is.null(obj[[key]])) next

    s_tag <- if (!is.null(season_name)) paste0("SEASON_", season_name, "_") else "ANNUAL_"
    ttl_season <- if (!is.null(season_name)) paste0(" ", season_name) else " annual"
     
    ttl <- if (lab == "B_IQR") {
      sprintf("%s%s IQR of relative uncertainty, thr=%.2f",
              year_label, ttl_season, res$threshold)
    }
    else {
      sprintf("%s%s relative uncertainty (%s) (IQR/\u03bc), thr=%.2f",
              year_label, ttl_season, lab, res$threshold)
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
  # Annual
  title_txt <- sprintf("%s annual wet hours (%s)", year_label, thr_label)

    plot_one(
      res$annual$wet_hours,
      sprintf("ANNUAL_WET_HOURS_%s_thr_%s.png", year_label, thr_txt),
      title_txt
    )

  # Seasonal
  for (s in names(res$seasonal)) {

    season_title_txt <- sprintf("%s %s wet hours (%s)", year_label, s, thr_label)

    plot_one(
    res$seasonal[[s]]$wet_hours,
    sprintf("SEASON_%s_WET_HOURS_%s_thr_%s.png", s, year_label, thr_txt),
    season_title_txt
  )
  invisible(TRUE)
}}

print("here")

res <- compute_year_threshold_stats_simple(
  rda_file = rda_file,
  min_threshold = 0.1,
  max_threshold = 0.5,
  mu_min = mu_min
)

plot_stats_simple(
    res = res,
    year_label = year,
    out_dir = out_dir,
    xlim = xlim,
    ylim = ylim,
    mode = mode
  )

message("Done.")

# Interactive run config
# -----------------------------

# Plot / compute options
# mode   <- "relunc"     # "all" or "relunc"
# mu_min <- 0.1       # mask for tiny mu in relunc

# thresholds <- c(1.0)  # can be c(0.1, 0.5, 1.0) etc


# # -----------------------------
# # RUN
# # -----------------------------
# stopifnot(file.exists(rda_file))

# message("Running year=", year, " file=", rda_file)
# message("Mode=", mode, " | mu_min=", mu_min)
# message("Plots -> ", normalizePath(out_dir, mustWork = FALSE))

# for (thr in thresholds) {
#   message("Threshold: ", thr)

#   res <- compute_year_threshold_stats_simple(
#     rda_file = rda_file,
#     threshold = thr,
#     mu_min = mu_min,
#     rel_uncert_method = "B_median"   # or "A"
#   )

#   plot_stats_simple(
#     res = res,
#     year_label = year,
#     out_dir = out_dir,
#     xlim = xlim,
#     ylim = ylim,
#     mode = mode
#   )
# }


   # plot_one(res$seasonal[[s]]$mean_mu,
    #          sprintf("SEASON_%s_MU_%s_thr_%s.png", s, year_label, thr_txt),
    #          sprintf("%s %s kriging mean (μ), thr=%.2f", year_label, s, res$threshold))

    # plot_one(res$seasonal[[s]]$mean_iqr,
    #          sprintf("SEASON_%s_IQR_%s_thr_%s.png", s, year_label, thr_txt),
    #          sprintf("%s %s IQR(90-10), thr=%.2f", year_label, s, res$threshold))

    # plot_one(res$seasonal[[s]]$accum_mu,
    #          sprintf("SEASON_%s_ACCUM_%s_thr_%s.png", s, year_label, thr_txt),
    #          sprintf("%s %s accumulation (mm), thr=%.2f", year_label, s, res$threshold))



  # plot_one(res$annual$wet_hours,
  #          sprintf("ANNUAL_WET_HOURS_%s_thr_%s.png", year_label, thr_txt),
  #          sprintf("%s annual wet hours (>%.2f)", year_label, res$threshold))

 # Annual
  # plot_one(res$annual$mean_mu,
  #          sprintf("ANNUAL_MU_%s_thr_%s.png", year_label, thr_txt),
  #          sprintf("%s annual kriging mean (μ), thr=%.2f", year_label, res$threshold))

  # plot_one(res$annual$mean_iqr,
  #          sprintf("ANNUAL_IQR_%s_thr_%s.png", year_label, thr_txt),
  #          sprintf("%s annual IQR(90-10), thr=%.2f", year_label, res$threshold))

  # plot_one(res$annual$accum_mu,
  #          sprintf("ANNUAL_ACCUM_%s_thr_%s.png", year_label, thr_txt),
  #          sprintf("%s annual accumulation (mm), thr=%.2f", year_label, res$threshold))

  # Annual
  # plot_relunc_variants(res$annual, season_name = NULL)

  # # Seasonal
  # for (s in names(res$seasonal)) {
  #   plot_relunc_variants(res$seasonal[[s]], season_name = s)
  # }

  # # If user asked ONLY relunc, stop here
  # if (mode == "relunc") return(invisible(TRUE))

  # Existing plots: only when mode == "all"
