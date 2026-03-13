
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

##### Config ###
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
stop("Usage: Rscript run_interannual.R <start_year> <end_year> [threshold_bands] [out_dir]")
}

start_year <- as.integer(args[1])
end_year   <- as.integer(args[2])

if (!is.finite(start_year) || !is.finite(end_year))
  stop("Years must be integers")

years <- seq(start_year, end_year)

# thresholds bands (comma separated like 0.01-0.1,0.1-0.5)
threshold_bands <- if (length(args) >= 3) {
  band_strings <- strsplit(args[3], ",")[[1]]

  lapply(band_strings, function(x) {
    vals <- as.numeric(strsplit(x, "-")[[1]])
    if (length(vals) != 2 || any(!is.finite(vals))) {
      stop("Invalid threshold band: ", x,
           ". Use format like 0.01-0.1,0.1-0.5")
    }
    if (vals[1] >= vals[2]) {
      stop("Threshold band must satisfy min < max: ", x)
    }
    vals
  })
} else {
  list(c(0.1, Inf))
}

# Optional
band_labels <- vapply(
  threshold_bands,
  function(x) paste0(x[1], "-", x[2]),
  character(1)
)

# Optional output directory
out_dir <- if (length(args) >= 4) {
  args[4]
} else {
  sprintf("/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats/out_plots/interannual_bands_%d_%d/",
          start_year, end_year)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Years: ", paste(years, collapse = ", "))
message("Threshold bands: ", paste(band_labels, collapse = ", "))
message("Output dir: ", normalizePath(out_dir, mustWork = FALSE))

# Load helper functions
source("/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats/R/utils.r")

compute_interannual_stats <- function(
  years,
  threshold_bands,
  rda_pattern = "/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/precip_transformed_results_new_%s.rda"
) {
  # Name each band like "0.01_0.1"
  band_names <- vapply(
    threshold_bands,
    function(x) paste0(x[1], "_", x[2]),
    character(1)
  )

  out <- vector("list", length(threshold_bands))
  names(out) <- band_names

  for (i in seq_along(threshold_bands)) {
    band <- threshold_bands[[i]]
    min_threshold <- band[1]
    max_threshold <- band[2]

    message("=== Threshold band ", min_threshold, " to ", max_threshold, " ===")

    variance_ref_list <- NULL

    # Initialize after first successful year
    acc_annual <- NULL
    acc_season <- NULL

    for (yr in years) {
      rda_file <- sprintf(rda_pattern, yr)

      if (!file.exists(rda_file)) {
        warning("Missing: ", rda_file)
        next
      }

      message(" Year ", yr)

      res <- compute_year_threshold_stats_simple(
        rda_file = rda_file,
        min_threshold = min_threshold,
        max_threshold = max_threshold
      )

      # Keep plotting coords from first valid year
      if (is.null(variance_ref_list)) {
        variance_ref_list <- res$variance_ref_list
      }

      # Initialize accumulators from first valid year
      if (is.null(acc_annual)) {
        acc_annual <- list(
          mean_mu   = acc_init(res$annual$mean_mu),
          mean_iqr  = acc_init(res$annual$mean_iqr),
          accum_mu  = acc_init(res$annual$accum_mu),
          wet_hours = acc_init(res$annual$wet_hours)
        )

        acc_season <- setNames(
          vector("list", length(names(res$seasonal))),
          names(res$seasonal)
        )

        for (s in names(res$seasonal)) {
          acc_season[[s]] <- list(
            mean_mu   = acc_init(res$seasonal[[s]]$mean_mu),
            mean_iqr  = acc_init(res$seasonal[[s]]$mean_iqr),
            accum_mu  = acc_init(res$seasonal[[s]]$accum_mu),
            wet_hours = acc_init(res$seasonal[[s]]$wet_hours)
          )
        }
      }

      # Annual add
      acc_annual$mean_mu   <- acc_add(acc_annual$mean_mu,   res$annual$mean_mu)
      acc_annual$mean_iqr  <- acc_add(acc_annual$mean_iqr,  res$annual$mean_iqr)
      acc_annual$accum_mu  <- acc_add(acc_annual$accum_mu,  res$annual$accum_mu)
      acc_annual$wet_hours <- acc_add(acc_annual$wet_hours, res$annual$wet_hours)

      # Seasonal add
      for (s in names(res$seasonal)) {
        acc_season[[s]]$mean_mu   <- acc_add(acc_season[[s]]$mean_mu,   res$seasonal[[s]]$mean_mu)
        acc_season[[s]]$mean_iqr  <- acc_add(acc_season[[s]]$mean_iqr,  res$seasonal[[s]]$mean_iqr)
        acc_season[[s]]$accum_mu  <- acc_add(acc_season[[s]]$accum_mu,  res$seasonal[[s]]$accum_mu)
        acc_season[[s]]$wet_hours <- acc_add(acc_season[[s]]$wet_hours, res$seasonal[[s]]$wet_hours)
      }
    }

    # Finalize means
    if (is.null(acc_annual)) {
      warning("No valid years processed for threshold band ", min_threshold, " to ", max_threshold)
      next
    }

    annual_mean <- list(
      mean_mu   = acc_mean(acc_annual$mean_mu),
      mean_iqr  = acc_mean(acc_annual$mean_iqr),
      accum_mu  = acc_mean(acc_annual$accum_mu),
      wet_hours = acc_mean(acc_annual$wet_hours)
    )

    seasonal_mean <- lapply(acc_season, function(a) {
      list(
        mean_mu   = acc_mean(a$mean_mu),
        mean_iqr  = acc_mean(a$mean_iqr),
        accum_mu  = acc_mean(a$accum_mu),
        wet_hours = acc_mean(a$wet_hours)
      )
    })

    out[[i]] <- list(
      min_threshold = min_threshold,
      max_threshold = max_threshold,
      years = years,
      interannual_mean = annual_mean,
      interseasonal_mean = seasonal_mean,
      variance_ref_list = variance_ref_list
    )
  }

  out
}


message("Start computation")

res_all <- compute_interannual_stats(
  years = years,
  threshold_bands = threshold_bands
)

plot_interannual_products <- function(
  res_thr, out_dir,
  xlim = c(480, 840), ylim = c(60, 300),
  cap_quant = 0.99, palette_end = 0.99
) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  min_thr <- res_thr$min_threshold
  max_thr <- res_thr$max_threshold

  band_txt <- paste0(
    gsub("\\.", "p", sprintf("%.2f", min_thr)),
    "_to_",
    gsub("\\.", "p", sprintf("%.2f", max_thr))
  )

  band_label <- sprintf("(%.2f, %.2f]", min_thr, max_thr)

  # annual IQR
  plot_cropped_field(
    res_thr$interannual_mean$mean_iqr,
    res_thr$variance_ref_list,
    xlim, ylim,
    main_title = sprintf(
      "10-year mean IQR (years %s-%s), band %s",
      min(res_thr$years), max(res_thr$years), band_label
    ),
    cap_quant = cap_quant,
    palette_end = palette_end,
    output_file = file.path(out_dir, sprintf("10-year_IQR_band_%s.png", band_txt))
  )

  # annual wet hours
  plot_cropped_field(
    res_thr$interannual_mean$wet_hours,
    res_thr$variance_ref_list,
    xlim, ylim,
    main_title = sprintf(
      "10-year mean wet hours (years %s-%s), band %s",
      min(res_thr$years), max(res_thr$years), band_label
    ),
    cap_quant = cap_quant,
    palette_end = palette_end,
    output_file = file.path(out_dir, sprintf("10-year_WETHOURS_band_%s.png", band_txt))
  )

  # seasonal
  for (s in names(res_thr$interseasonal_mean)) {
    ss <- res_thr$interseasonal_mean[[s]]

    plot_cropped_field(
      ss$mean_iqr,
      res_thr$variance_ref_list,
      xlim, ylim,
      main_title = sprintf("%s mean IQR across 2016-2025, band %s", s, band_label),
      cap_quant = cap_quant,
      palette_end = palette_end,
      output_file = file.path(out_dir, sprintf("10_year_%s_IQR_band_%s.png", s, band_txt))
    )

    plot_cropped_field(
      ss$wet_hours,
      res_thr$variance_ref_list,
      xlim, ylim,
      main_title = sprintf("%s mean wet hours across 2016-2025, band %s", s, band_label),
      cap_quant = cap_quant,
      palette_end = palette_end,
      output_file = file.path(out_dir, sprintf("10_year_%s_WETHOURS_band_%s.png", s, band_txt))
    )
  }
}

message("Start Plotting")

for (band_name in names(res_all)) {

  plot_interannual_products(
    res_all[[band_name]],
    out_dir = file.path(out_dir, paste0("band_", band_name))
  )

}

message("Done.")


#years <- 2016:2025
#thresholds <- c(0.1)

#all_thr <- compute_interannual_stats(years, thresholds)

#out_dir <- file.path("out_plots", sprintf("interannual_%d_%d", min(years), max(years)))
#dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)