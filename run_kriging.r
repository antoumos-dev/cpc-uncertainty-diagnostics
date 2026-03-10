.libPaths("/store_new/mch/msclim/share/CATs/cats/lib-R4.4.0/")
library(raster)
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
library(gstat)                                                      #but I had to install sf first
library(rhdf5)
library(animation)
library(ncdf4)                     
library(lubridate)
library(dplyr)

setwd("/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/")

############################################### modular template ####################################################

## install.packages("matrixStats")
.libPaths("/store_new/mch/msclim/antoumos/R/lib")

library(matrixStats)

source("out_stats/R/utils.r") ### load helpers function

##### accumulation helper ###

aggregate_sum <- function(field_list) {
  ny <- nrow(field_list[[1]])
  nx <- ncol(field_list[[1]])

  s <- matrix(0, ny, nx)
  n <- matrix(0L, ny, nx)

  for (x in field_list) {
    ok <- is.finite(x)
    s[ok] <- s[ok] + x[ok]
    n[ok] <- n[ok] + 1L
  }

  s[n == 0] <- NA_real_
  list(sum = s, n_valid = n)
}

compute_year_threshold_stats <- function(rda_file, threshold) {

  load(rda_file)

  dates <- as.POSIXct(timestamps, origin = "1970-01-01", tz = "UTC")
  season <- get_season(dates)

  variance_raw_list <- lapply(variance_crop_list, strip_attrs_keep_dim)

  kriging_thr_list <- lapply(kriging_crop_list, function(m) {
    m[m < threshold] <- NA_real_
    strip_attrs_keep_dim(m)
  })

  # ---- Annual
  annual_mean_mu  <- aggregate_mean(kriging_thr_list)$mean
  annual_wet_hours <- aggregate_wet_hours(kriging_thr_list, threshold)

  iqr_list <- compute_iqr_list(kriging_thr_list, variance_raw_list)
  annual_mean_iqr <- aggregate_mean(iqr_list)$mean
  
  annual_skew_mu_log  <- temporal_skewness_blocked(kriging_thr_list, min_n = 100, eps = 1e-6)
  annual_skew_iqr_log <- temporal_skewness_blocked(iqr_list, min_n = 100, eps = 1e-6)
   

  annual_accum_mu <- aggregate_sum(kriging_thr_list)$sum

  # ---- Seasonal
  seasons <- levels(season)
  seasonal <- setNames(vector("list", length(seasons)), seasons)

  for (s in seasons) {
    idx <- which(season == s)

    seasonal[[s]] <- list(
      mean_mu   = aggregate_mean(kriging_thr_list[idx])$mean,
      accum_mu  = aggregate_sum(kriging_thr_list[idx])$sum ,
      wet_hours = aggregate_wet_hours(kriging_thr_list[idx], threshold),
      mean_iqr  = aggregate_mean(iqr_list[idx])$mean
    )
  }

  list(
    threshold = threshold,
    annual = list(
      mean_mu = annual_mean_mu,
      accum_mu = annual_accum_mu,
      wet_hours = annual_wet_hours,
      mean_iqr = annual_mean_iqr,
      skew_mu_log = annual_skew_mu_log,
      skew_iqr_log = annual_skew_iqr_log
    ),
    seasonal = seasonal,
    variance_raw_list = variance_raw_list  # needed for plotting
  )
}


plot_seasonal_fields <- function(stats, year_label, out_dir,
                                 thr_txt,
                                 xlim = c(480,840),
                                 ylim = c(60,300),
                                 cap_quant = 0.99,
                                 palette_end = 0.99) {

  if (is.null(stats$seasonal) || length(stats$seasonal) == 0) return(invisible(FALSE))

  for (s in names(stats$seasonal)) {

    # Seasonal kriging mean
    plot_cropped_field(
      Z = stats$seasonal[[s]]$mean_mu,
      variance_list = stats$variance_raw_list,
      xlim = xlim, ylim = ylim,
      title = sprintf("%s %s kriging mean, thr=%.2f", year_label, s, stats$threshold),
      cap_quant = cap_quant,
      palette_end = palette_end,
      output_file = file.path(out_dir, sprintf("MU_mean_%s_%s_thr_%s.png", year_label, s, thr_txt))
    )

    # Seasonal IQR
    plot_cropped_field(
      Z = stats$seasonal[[s]]$mean_iqr,
      variance_list = stats$variance_raw_list,
      xlim = xlim, ylim = ylim,
      title = sprintf("%s %s IQR(90-10), thr=%.2f", year_label, s, stats$threshold),
      cap_quant = cap_quant,
      palette_end = palette_end,
      output_file = file.path(out_dir, sprintf("IQR_mean_%s_%s_thr_%s.png", year_label, s, thr_txt))
    )
  }

  # Seasonal wet hours 


  # Accumulations


  invisible(TRUE)
}




export_year_threshold_results <- function(stats, year_label,
                                          out_dir = "out_stats",
                                          cap_quant = 0.99,
                                          palette_end = 0.99,
                                          xlim = c(480,840),
                                          ylim = c(60,300)) {

  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  thr_txt <- gsub("\\.", "p", sprintf("%.2f", stats$threshold))  # safer filenames like 0p10
  
 plot_cropped_field(
  Z = stats$annual$accum_mu,
  variance_list = stats$variance_raw_list,
  xlim = xlim,
  ylim = ylim,
  title = sprintf("%s annual accumulation (mm), thr=%.2f",
                  year_label, stats$threshold),
  output_file = file.path(out_dir,
                          sprintf("ANNUAL_ACCUM_%s_thr_%s.png",
                                  year_label, thr_txt))
)

}
 

###### Interannual means #####

compute_interannual_means <- function(years, threshold,
                                      file_pattern = "../data_new_project/precip_transformed_results_new_%d.rda") {

  # We'll accumulate sum + count per pixel
  mu_sum <- NULL
  mu_n   <- NULL
  iqr_sum <- NULL
  iqr_n   <- NULL

  variance_ref <- NULL  # keep one ref list for plotting coords if needed

  for (yr in years) {

    rda_file <- sprintf(file_pattern, yr)
    if (!file.exists(rda_file)) {
      warning("Missing file: ", rda_file)
      next
    }

    message("Year ", yr, " (thr=", threshold, ")")
    stats <- compute_year_threshold_stats(rda_file, threshold = threshold)

    mu  <- stats$annual$mean_mu
    iqr <- stats$annual$mean_iqr

    # initialize accumulators on first successful year
    if (is.null(mu_sum)) {
      mu_sum  <- matrix(0, nrow(mu), ncol(mu))
      mu_n    <- matrix(0L, nrow(mu), ncol(mu))
      iqr_sum <- matrix(0, nrow(iqr), ncol(iqr))
      iqr_n   <- matrix(0L, nrow(iqr), ncol(iqr))
      variance_ref <- stats$variance_raw_list
    }

    ok_mu <- is.finite(mu)
    mu_sum[ok_mu] <- mu_sum[ok_mu] + mu[ok_mu]
    mu_n[ok_mu]   <- mu_n[ok_mu] + 1L

    ok_iqr <- is.finite(iqr)
    iqr_sum[ok_iqr] <- iqr_sum[ok_iqr] + iqr[ok_iqr]
    iqr_n[ok_iqr]   <- iqr_n[ok_iqr] + 1L
  }

  if (is.null(mu_sum)) stop("No valid years processed. Check file paths.")

  interannual_mu  <- mu_sum  / pmax(mu_n, 1L)
  interannual_iqr <- iqr_sum / pmax(iqr_n, 1L)

  interannual_mu[mu_n == 0]   <- NA_real_
  interannual_iqr[iqr_n == 0] <- NA_real_

  list(
    threshold = threshold,
    years = years,
    interannual = list(
      mean_mu  = interannual_mu,
      mean_iqr = interannual_iqr,
      n_mu     = mu_n,
      n_iqr    = iqr_n
    ),
    variance_ref_list = variance_ref
  )
}

### call the function for interannual ###

years <- 2016:2025
res_01 <- compute_interannual_means(years, threshold = 0.1)

plot_cropped_field(
  Z = res_01$interannual$mean_mu,
  variance_list = res_01$variance_ref_list,
  xlim = c(480, 840),
  ylim = c(60, 300),
  title = sprintf("Interannual mean kriging mean (2016–2025), thr=%.2f", res_01$threshold),
  output_file = "out_stats/INTERANNUAL_MU_thr_0p10.png"
)

plot_cropped_field(
  Z = res_01$interannual$mean_iqr,
  variance_list = res_01$variance_ref_list,
  xlim = c(480, 840),
  ylim = c(60, 300),
  title = sprintf("Interannual mean IQR(90-10) (2016–2025), thr=%.2f", res_01$threshold),
  output_file = "out_stats/INTERANNUAL_IQR_thr_0p10.png"


)

# --------------------------
# top-level driver
# --------------------------

years <- c("2021")
rda_files <- c(
#  "precip_transformed_results_new_2021.rda",
  "data/precip_transformed_results_new_2021.rda"
)

thresholds <- c(0.1)

run_pipeline <- function(years, rda_files, thresholds) {
  for (i in seq_along(years)) {
    for (thr in thresholds) {
      message("Processing: ", years[i], " threshold=", thr)

      res <- compute_year_threshold_stats(rda_files[i], threshold = thr)

      export_year_threshold_results(res,
                                    year_label = years[i],
                                    out_dir = "out_stats")
    }
  }
}

run_pipeline(years, rda_files, thresholds)




 # Save numeric results
  # saveRDS(
  #   stats,
  #   file = file.path(out_dir,
  #                    sprintf("stats_%s_thr_%s.rds",
  #                            year_label, stats$threshold))
  # )

  # # Annual plots
  # plot_cropped_field(
  #   Z = stats$annual$mean_iqr,
  #   variance_list = stats$variance_raw_list,
  #   xlim = xlim, ylim = ylim,
  #   title = sprintf("%s mean IQR thr=%.2f",
  #                   year_label, stats$threshold),
  #   cap_quant = 0.99,
  #   palette_end = 0.99,
  #   output_file = file.path(out_dir,
  #                           sprintf("iqr_%s_thr_%s.png",
  #                                   year_label, stats$threshold))
  # )

    # Mean kriging mean
  # plot_cropped_field(
  #   Z = stats$annual$mean_mu,
  #   variance_list = stats$variance_raw_list,
  #   xlim = xlim, ylim = ylim,
  #   title = sprintf("%s annual kriging mean, thr=%.2f", year_label, stats$threshold),
  #   cap_quant = cap_quant,
  #   palette_end = palette_end,
  #   output_file = file.path(out_dir, sprintf("MU_mean_%s_thr_%s.png", year_label, thr_txt))
  # )

  # Wet hours
  # plot_cropped_field(
  #   Z = stats$annual$wet_hours,
  #   variance_list = stats$variance_raw_list,
  #   xlim = xlim, ylim = ylim,
  #   title = sprintf("%s annual wet hours (> %.2f)", year_label, stats$threshold),
  #   cap_quant = cap_quant,
  #   palette_end = palette_end,
  #   output_file = file.path(out_dir, sprintf("WET_hours_%s_thr_%s.png", year_label, thr_txt))
  # )

# Annual skewness of log(mu)
# plot_cropped_field(
#   Z = stats$annual$skew_mu_log,
#   variance_list = stats$variance_raw_list,
#   xlim = xlim, ylim = ylim,
#   title = sprintf("%s annual skewness log1p(mu), thr=%.2f",
#                   year_label, stats$threshold),
#   cap_quant = cap_quant,
#   palette_end = palette_end,
#   output_file = file.path(out_dir,
#                           sprintf("SKEW_logMU_%s_thr_%s.png",
#                                   year_label, thr_txt))
# )

# Annual skewness of log(IQR)
# plot_cropped_field(
#   Z = stats$annual$skew_iqr_log,
#   variance_list = stats$variance_raw_list,
#   xlim = xlim, ylim = ylim,
#   title = sprintf("%s annual skewness log1p(IQR), thr=%.2f",
#                   year_label, stats$threshold),
#   cap_quant = cap_quant,
#   palette_end = palette_end,
#   output_file = file.path(out_dir,
#                           sprintf("SKEW_logIQR_%s_thr_%s.png",
#                                   year_label, thr_txt))
# )

# plot_seasonal_fields(stats, year_label, out_dir, thr_txt,
                      # xlim = xlim, ylim = ylim,
                      #  cap_quant = cap_quant, palette_end = palette_end)


# plot_cropped_field(
#   Z = res_01$interannual$mean_mu,
#   variance_list = res_01$variance_ref_list,
#   xlim = c(480, 840),
#   ylim = c(60, 300),
#   title = sprintf("kriging mean (2016–2025), thr=%.2f", res_01$threshold),
#   output_file = "out_stats/INTERANNUAL_MU_thr_0p10.png"
# )

# plot_cropped_field(
#   Z = res_01$interannual$mean_iqr,
#   variance_list = res_01$variance_ref_list,
#   xlim = c(480, 840),
#   ylim = c(60, 300),
#   title = sprintf("Mean IQR(90-10) (2016–2025), thr=%.2f", res_01$threshold),
#   output_file = "out_stats/INTERANNUAL_IQR_thr_0p10.png"
# )


