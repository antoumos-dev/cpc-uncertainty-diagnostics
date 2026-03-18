
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
  stop("Usage: Rscript run_interannual.R <start_year> <end_year> [thresholds] [out_dir]")
}

start_year <- as.integer(args[1])
end_year   <- as.integer(args[2])

if (!is.finite(start_year) || !is.finite(end_year))
  stop("Years must be integers")

years <- seq(start_year, end_year)

# Optional thresholds (comma separated like 0.1,1,5)
thresholds <- if (length(args) >= 3) {
  as.numeric(strsplit(args[3], ",")[[1]])
} else {
  c(0.1)
}

if (any(!is.finite(thresholds)))
  stop("Invalid thresholds")

# Optional output directory
out_dir <- if (length(args) >= 4) {
  args[4]
} else {
  sprintf("/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats/out_plots/interannual_%d_%d/",
          start_year, end_year)
}

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Years: ", paste(years, collapse = ", "))
message("Thresholds: ", paste(thresholds, collapse = ", "))
message("Output dir: ", normalizePath(out_dir, mustWork = FALSE))


# Load helper functions
source("/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats/R/utils.r")

res_all <- compute_interannual_stats(
  years = years,
  thresholds = thresholds
)

#for (thr in names(res_all)) {
# plot_interannual_products(res_all[[thr]], out_dir)
#plot_interannual_sd_products(res_all[[thr]], out_dir)
#}

out_dir <- ("/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats/out_plots/interannual_2016_2025_2/")



#### function for interactive ###
compute_interannual_stats <- function(years, thresholds,rda_pattern = "/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/precip_transformed_results_new_%s.rda") {

  out <- vector("list", length(thresholds))
  names(out) <- as.character(thresholds)

  for (thr in thresholds) {
    message("=== Threshold ", thr, " ===")
    variance_ref_list <- NULL

    # We will initialize these after first successful year (when we know matrix dims)
    acc_annual <- NULL
    acc_annual_sd <- NULL
    acc_season <- NULL
    acc_season_sd <- NULL

    for (yr in years) {
      rda_file <- sprintf(rda_pattern, yr)
      if (!file.exists(rda_file)) {
        warning("Missing: ", rda_file)
        next
      }

      message(" Year ", yr)
      res <- compute_year_threshold_stats_simple(rda_file, threshold = thr,rel_uncert_method = "A")
      #years_done <- c(years_done, yr)

      # keep plotting coords from first valid year
      if (is.null(variance_ref_list)) variance_ref_list <- res$variance_ref_list

      # initialize accumulators from first year's matrix shapes
      if (is.null(acc_annual)) {
        acc_annual <- list(
          mean_mu   = acc_init(res$annual$mean_mu),
          mean_iqr  = acc_init(res$annual$mean_iqr),
          accum_mu  = acc_init(res$annual$accum_mu),
          wet_hours = acc_init(res$annual$wet_hours)
        )
        #SD accumulators (Welford)
        acc_annual_sd <- list(
        mean_mu   = acc_sd_init(res$annual$mean_mu),
        mean_iqr  = acc_sd_init(res$annual$mean_iqr),
        accum_mu  = acc_sd_init(res$annual$accum_mu),
        wet_hours = acc_sd_init(res$annual$wet_hours)
        )
        acc_season <- setNames(vector("list", length(names(res$seasonal))), names(res$seasonal))
        acc_season_sd <- setNames(vector("list", length(names(res$seasonal))), names(res$seasonal))
        for (s in names(res$seasonal)) {
          acc_season[[s]] <- list(
          mean_mu   = acc_init(res$seasonal[[s]]$mean_mu),
          mean_iqr  = acc_init(res$seasonal[[s]]$mean_iqr),
          accum_mu  = acc_init(res$seasonal[[s]]$accum_mu),
          wet_hours = acc_init(res$seasonal[[s]]$wet_hours)
          )
          acc_season_sd[[s]] <- list(
          mean_mu   = acc_sd_init(res$seasonal[[s]]$mean_mu),
          mean_iqr  = acc_sd_init(res$seasonal[[s]]$mean_iqr),
          accum_mu  = acc_sd_init(res$seasonal[[s]]$accum_mu),
          wet_hours = acc_sd_init(res$seasonal[[s]]$wet_hours)
        )
      }
    }

      # annual add
      acc_annual$mean_mu   <- acc_add(acc_annual$mean_mu,   res$annual$mean_mu)
      acc_annual$mean_iqr  <- acc_add(acc_annual$mean_iqr,  res$annual$mean_iqr)
      acc_annual$accum_mu  <- acc_add(acc_annual$accum_mu,  res$annual$accum_mu)
      acc_annual$wet_hours <- acc_add(acc_annual$wet_hours, res$annual$wet_hours)

      # seasonal add
      for (s in names(res$seasonal)) {
        acc_season[[s]]$mean_mu   <- acc_add(acc_season[[s]]$mean_mu,   res$seasonal[[s]]$mean_mu)
        acc_season[[s]]$mean_iqr  <- acc_add(acc_season[[s]]$mean_iqr,  res$seasonal[[s]]$mean_iqr)
        acc_season[[s]]$accum_mu  <- acc_add(acc_season[[s]]$accum_mu,  res$seasonal[[s]]$accum_mu)
        acc_season[[s]]$wet_hours <- acc_add(acc_season[[s]]$wet_hours, res$seasonal[[s]]$wet_hours)
    }
    
      # SD add
      # annual
      acc_annual_sd$mean_mu   <- acc_sd_add(acc_annual_sd$mean_mu,   res$annual$mean_mu)
      acc_annual_sd$mean_iqr  <- acc_sd_add(acc_annual_sd$mean_iqr,  res$annual$mean_iqr)
      acc_annual_sd$accum_mu  <- acc_sd_add(acc_annual_sd$accum_mu,  res$annual$accum_mu)
      acc_annual_sd$wet_hours <- acc_sd_add(acc_annual_sd$wet_hours, res$annual$wet_hours)

      # seasonal
      for (s in names(res$seasonal)) {
      acc_season_sd[[s]]$mean_mu   <- acc_sd_add(acc_season_sd[[s]]$mean_mu,   res$seasonal[[s]]$mean_mu)
      acc_season_sd[[s]]$mean_iqr  <- acc_sd_add(acc_season_sd[[s]]$mean_iqr,  res$seasonal[[s]]$mean_iqr)
      acc_season_sd[[s]]$accum_mu  <- acc_sd_add(acc_season_sd[[s]]$accum_mu,  res$seasonal[[s]]$accum_mu)
      acc_season_sd[[s]]$wet_hours <- acc_sd_add(acc_season_sd[[s]]$wet_hours, res$seasonal[[s]]$wet_hours)
    }
  }
    # finalize means
    if (is.null(acc_annual)) {
      warning("No valid years processed for threshold ", thr)
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
    # finalize SD (new)
    # annual_sd <- list(
    #   mean_mu   = acc_sd_finalize(acc_annual_sd$mean_mu),
    #   mean_iqr  = acc_sd_finalize(acc_annual_sd$mean_iqr),
    #   accum_mu  = acc_sd_finalize(acc_annual_sd$accum_mu),
    #   wet_hours = acc_sd_finalize(acc_annual_sd$wet_hours)
    # )

    # seasonal_sd <- lapply(acc_season_sd, function(a) {
    #   list(
    #     mean_mu   = acc_sd_finalize(a$mean_mu),
    #     mean_iqr  = acc_sd_finalize(a$mean_iqr),
    #     accum_mu  = acc_sd_finalize(a$accum_mu),
    #     wet_hours = acc_sd_finalize(a$wet_hours)
    #   )
    # })

      out[[as.character(thr)]] <- list(
      threshold = thr,
      years = years,
      interannual_mean = annual_mean,
      interseasonal_mean = seasonal_mean,
      interannual_sd = annual_sd,
      interseasonal_sd = seasonal_sd,
      variance_ref_list = variance_ref_list
    )
}
  out
}

get_scale_for_variable <- function(interannual_field, interseasonal_list, field_name,
                                   probs = c(0.02, 0.9)) {
  fields <- c(
    list(interannual_field),
    lapply(interseasonal_list, function(ss) ss[[field_name]])
  )
  get_common_color_scale(fields, probs = probs)
}


plot_interannual_products <- function(res_thr, out_dir,
                                      xlim = c(480,840), ylim = c(60,300),
                                      cap_quant = 0.99, palette_end = 0.99) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  thr_txt <- gsub("\\.", "p", sprintf("%.2f", res_thr$threshold))


relunc_scale <- get_scale_for_variable(
  res_thr$interannual$rel_uncert_Bmed,
  res_thr$interseasonal,
  "rel_uncert_Bmed"
)

mu_scale <- get_scale_for_variable(
  res_thr$interannual$mean_mu,
  res_thr$interseasonal,
  "mean_mu",
  probs = c(0.03,0.95)

)

iqr_scale <- get_scale_for_variable(
  res_thr$interannual$mean_iqr,
  res_thr$interseasonal,
  "mean_iqr",
  probs = c(0.03,0.95)
)

accum_scale <- get_scale_for_variable(
  res_thr$interannual$accum_mu,
  res_thr$interseasonal,
  "accum_mu",
  probs = c(0.05,0.97)
)

wet_scale <- get_scale_for_variable(
  res_thr$interannual$wet_hours,
  res_thr$interseasonal,
  "wet_hours",
  probs = c(0.05,0.97)
)

  # annual
 plot_cropped_field(
  res_thr$interannual$mean_mu,
  res_thr$variance_ref_list,
  xlim, ylim,
  vmin = mu_scale$vmin,
  vmax = mu_scale$vmax,
  title = sprintf("Interannual mean μ (years %s-%s), thr=%.2f",
                  min(res_thr$years), max(res_thr$years), res_thr$threshold),
  cap_quant = cap_quant,
  palette_end = palette_end,
  output_file = file.path(out_dir, sprintf("INTERANNUAL_MU_thr_%s.png", thr_txt))
)

plot_cropped_field(
  res_thr$interannual$mean_iqr,
  res_thr$variance_ref_list,
  xlim, ylim,
  vmin = iqr_scale$vmin,
  vmax = iqr_scale$vmax,
  title = sprintf("Interannual mean IQR(90-10), thr=%.2f", res_thr$threshold),
  cap_quant = cap_quant,
  palette_end = palette_end,
  output_file = file.path(out_dir, sprintf("INTERANNUAL_IQR_thr_%s.png", thr_txt))
)

plot_cropped_field(
  res_thr$interannual$accum_mu,
  res_thr$variance_ref_list,
  xlim, ylim,
  vmin = accum_scale$vmin,
  vmax = accum_scale$vmax,
  title = sprintf("Interannual mean accumulation (mm), thr=%.2f", res_thr$threshold),
  cap_quant = cap_quant,
  palette_end = palette_end,
  output_file = file.path(out_dir, sprintf("INTERANNUAL_ACCUM_thr_%s.png", thr_txt))
)

plot_cropped_field(
  res_thr$interannual$wet_hours,
  res_thr$variance_ref_list,
  xlim, ylim,
  vmin = wet_scale$vmin,
  vmax = wet_scale$vmax,
  title = sprintf("Interannual mean wet hours, thr=%.2f", res_thr$threshold),
  cap_quant = cap_quant,
  palette_end = palette_end,
  output_file = file.path(out_dir, sprintf("INTERANNUAL_WETHOURS_thr_%s.png", thr_txt))
)

  plot_cropped_field(
    res_thr$interannual$rel_uncert_Bmed,
    res_thr$variance_ref_list,
    xlim, ylim,
    vmin = relunc_scale$vmin,
    vmax = relunc_scale$vmax,
    title = sprintf("10 year mean of annual median(IQR/μ), thr=%.2f", res_thr$threshold),
    cap_quant = cap_quant,
    palette_end = palette_end,
    output_file = file.path(out_dir, sprintf("INTERANNUAL_RELUNC_thr_%s.png", thr_txt))
  )

#### different scales for the seasons
  accum_seasonal_scale <- get_scale_for_variable(
  res_thr$interannual$accum_mu,
  res_thr$interseasonal,
  "accum_mu",
  probs = c(0.02,0.8)
  )

  wet_seasonal_scale <- get_scale_for_variable(
  res_thr$interannual$wet_hours,
  res_thr$interseasonal,
  "accum_mu",
  probs = c(0.02,0.8)
  )

  # seasonal
  for (s in names(res_thr$interseasonal)) {
    ss <- res_thr$interseasonal[[s]]

plot_cropped_field(
    ss$mean_mu,
    res_thr$variance_ref_list,
    xlim, ylim,
    vmin = mu_scale$vmin,
    vmax = mu_scale$vmax,
    title = sprintf("%s mean μ across years, thr=%.2f", s, res_thr$threshold),
    cap_quant = cap_quant,
    palette_end = palette_end,
    output_file = file.path(out_dir, sprintf("INTERSEASON_%s_MU_thr_%s.png", s, thr_txt))
  )

  plot_cropped_field(
    ss$mean_iqr,
    res_thr$variance_ref_list,
    xlim, ylim,
    vmin = iqr_scale$vmin,
    vmax = iqr_scale$vmax,
    title = sprintf("%s mean IQR across years, thr=%.2f", s, res_thr$threshold),
    cap_quant = cap_quant,
    palette_end = palette_end,
    output_file = file.path(out_dir, sprintf("INTERSEASON_%s_IQR_thr_%s.png", s, thr_txt))
  )

  plot_cropped_field(
    ss$accum_mu,
    res_thr$variance_ref_list,
    xlim, ylim,
    vmin = accum_seasonal_scale$vmin,
    vmax = accum_seasonal_scale$vmax,
    title = sprintf("%s mean accumulation across years, thr=%.2f", s, res_thr$threshold),
    cap_quant = cap_quant,
    palette_end = palette_end,
    output_file = file.path(out_dir, sprintf("INTERSEASON_%s_ACCUM_thr_%s.png", s, thr_txt))
  )

  plot_cropped_field(
    ss$wet_hours,
    res_thr$variance_ref_list,
    xlim, ylim,
    vmin = wet_seasonal_scale$vmin,
    vmax = wet_seasonal_scale$vmax,
    title = sprintf("%s mean wet hours across years, thr=%.2f", s, res_thr$threshold),
    cap_quant = cap_quant,
    palette_end = palette_end,
    output_file = file.path(out_dir, sprintf("INTERSEASON_%s_WETHOURS_thr_%s.png", s, thr_txt))
  )

   plot_cropped_field(
  ss$rel_uncert_Bmed,
  res_thr$variance_ref_list,
  xlim, ylim,
  vmin = relunc_scale$vmin,
  vmax = relunc_scale$vmax,
  title = sprintf("%s 10 year mean of annual median(IQR/μ), thr=%.2f", s, res_thr$threshold),
  cap_quant = cap_quant,
  palette_end = palette_end,
  output_file = file.path(out_dir, sprintf("INTERSEASON_%s_RELUNC_thr_%s.png", s, thr_txt))
  )
    
  }
}





plot_interannual_sd_products <- function(res_thr, out_dir,xlim = c(480,840), ylim = c(60,300),cap_quant = 0.99, palette_end = 0.99) {

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  thr_txt <- gsub("\\.", "p", sprintf("%.2f", res_thr$threshold))

  yr_min <- min(res_thr$years)
  yr_max <- max(res_thr$years)

  # -----------------------------
  # ANNUAL SD
  # -----------------------------

  plot_cropped_field(res_thr$interannual_sd$mean_mu$sd,
                     res_thr$variance_ref_list, xlim, ylim,
                     title = sprintf("Interannual SD μ (%s–%s), thr=%.2f", yr_min, yr_max, res_thr$threshold),
                     cap_quant = cap_quant, palette_end = palette_end,
                     output_file = file.path(out_dir, sprintf("INTERANNUAL_SD_MU_thr_%s.png", thr_txt)))

  plot_cropped_field(res_thr$interannual_sd$mean_iqr$sd,
                     res_thr$variance_ref_list, xlim, ylim,
                     title = sprintf("Interannual SD IQR(90-10), thr=%.2f",res_thr$threshold),
                     cap_quant = cap_quant, palette_end = palette_end,
                     output_file = file.path(out_dir,sprintf("INTERANNUAL_SD_IQR_thr_%s.png", thr_txt)))

  plot_cropped_field(res_thr$interannual_sd$accum_mu$sd,
                     res_thr$variance_ref_list, xlim, ylim,
                     title = sprintf("Interannual SD accumulation (mm), thr=%.2f",res_thr$threshold),
                     cap_quant = cap_quant, palette_end = palette_end,
                     output_file = file.path(out_dir,sprintf("INTERANNUAL_SD_ACCUM_thr_%s.png", thr_txt)))

  plot_cropped_field(res_thr$interannual_sd$wet_hours$sd,
                     res_thr$variance_ref_list, xlim, ylim,
                     title = sprintf("Interannual SD wet hours, thr=%.2f",res_thr$threshold),
                     cap_quant = cap_quant, palette_end = palette_end,
                     output_file = file.path(out_dir,sprintf("INTERANNUAL_SD_WETHOURS_thr_%s.png", thr_txt)))

  # -----------------------------
  # SEASONAL SD
  # -----------------------------

  for (s in names(res_thr$interseasonal_sd)) {

    ss_sd <- res_thr$interseasonal_sd[[s]]

    plot_cropped_field(ss_sd$mean_mu$sd,
                       res_thr$variance_ref_list, xlim, ylim,
                       title = sprintf("%s interannual SD μ, thr=%.2f",
                                       s, res_thr$threshold),
                       cap_quant = cap_quant, palette_end = palette_end,
                       output_file = file.path(out_dir,
                         sprintf("INTERSEASON_%s_SD_MU_thr_%s.png", s, thr_txt)))

    plot_cropped_field(ss_sd$mean_iqr$sd,
                       res_thr$variance_ref_list, xlim, ylim,
                       title = sprintf("%s interannual SD IQR, thr=%.2f",
                                       s, res_thr$threshold),
                       cap_quant = cap_quant, palette_end = palette_end,
                       output_file = file.path(out_dir,
                         sprintf("INTERSEASON_%s_SD_IQR_thr_%s.png", s, thr_txt)))

    plot_cropped_field(ss_sd$accum_mu$sd,
                       res_thr$variance_ref_list, xlim, ylim,
                       title = sprintf("%s interannual SD accumulation, thr=%.2f",
                                       s, res_thr$threshold),
                       cap_quant = cap_quant, palette_end = palette_end,
                       output_file = file.path(out_dir,
                         sprintf("INTERSEASON_%s_SD_ACCUM_thr_%s.png", s, thr_txt)))

    plot_cropped_field(ss_sd$wet_hours$sd,
                       res_thr$variance_ref_list, xlim, ylim,
                       title = sprintf("%s interannual SD wet hours, thr=%.2f",
                                       s, res_thr$threshold),
                       cap_quant = cap_quant, palette_end = palette_end,
                       output_file = file.path(out_dir,
                         sprintf("INTERSEASON_%s_SD_WETHOURS_thr_%s.png", s, thr_txt)))
  }

  invisible(TRUE)
}

#years <- 2016:2025
#thresholds <- c(0.1)

#all_thr <- compute_interannual_stats(years, thresholds)

#out_dir <- file.path("out_plots", sprintf("interannual_%d_%d", min(years), max(years)))
#dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Start plotting")

#### could move plotting function into plotting utils later

for (thr in names(res_all)) {

  plot_interannual_products(
    res_all[[thr]],
    out_dir = file.path(out_dir, paste0("thr_", gsub("\\.", "p", thr)))
  )

}


message("Done.")