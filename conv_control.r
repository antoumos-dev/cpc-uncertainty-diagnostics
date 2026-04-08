library(raster)
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
library(gstat)                                                      #but I had to install sf first
library(rhdf5)
library(animation)
library(ncdf4)                     
library(lubridate)
library(dplyr)
library(data.table)



setwd("/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/")

project_root <- "/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats"
utils_file   <- file.path(project_root, "R", "utils.r")
source(utils_file) # helpers function


rda_files <- list.files(pattern = "CPC23.*\\.rda$") 

#load(paste("precip_transformed_results_new_2023.rda"))

#dates <- as.POSIXct(timestamps, origin = "1970-01-01", tz = "UTC") ### have the dates


# Format: CPC YY DDD HH FFF .rda
#          4   2  3   2  3
parse_date_from_filename <- function(fname) {
  base <- sub("\\.rda$", "", basename(fname))  # strip extension
  yy  <- substr(base, 4, 5)   # characters 4-5: year
  doy <- substr(base, 6, 8)   # characters 6-8: day of year
  hh  <- substr(base, 9, 10)  # characters 9-10: hour
  
  as.POSIXct(
    strptime(paste0("20", yy, " ", doy, " ", hh), 
             format = "%Y %j %H", tz = "UTC")
  )
}

dates <- as.POSIXct(sapply(rda_files, parse_date_from_filename), 
                    origin = "1970-01-01", tz = "UTC")


# ── 1. Build the 8760-element list ────────────────────────────────────────────
all_data <- vector("list", length(rda_files))   # length should be 8760

for (i in seq_along(rda_files)) {
  load(rda_files[i])                            # loads 'output'
  all_data[[i]] <- output[[8]][[2]]@data
  if (i %% 100 == 0) cat(sprintf("Progress: %d / %d\n", i, length(rda_files)))
}
names(all_data) <- format(dates, "%Y-%m-%d %H:%M") ### add dates as a character string


### save all-data 

saveRDS(all_data, file = "../new_project/out_stats/all_data_2023.rds")

setwd("/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats")

all_data <- readRDS("all_data_2023.rds")

# Exceedance frequency of 0.5 per grid cell for convective control activation ──────────────────────────────
threshold <- 0.5

coeff_matrix <- do.call(cbind, lapply(all_data, function(df) df$coef.var))

activations_ts <- data.frame(
  date    = dates,
  count   = as.numeric(colSums(coeff_matrix > threshold, na.rm = TRUE)),
  fraction = as.numeric(colMeans(coeff_matrix > threshold, na.rm = TRUE))
)

activations_ts$month  <- as.integer(format(dates, "%m"))

activations_ts$season <- cut(activations_ts$month,
                             breaks = c(0, 2, 5, 8, 11, 12),
                             labels = c("DJF", "MAM", "JJA", "SON", "DJF2"),
                             right  = TRUE)
# Merge DJF2 (December) back into DJF
levels(activations_ts$season)[levels(activations_ts$season) == "DJF2"] <- "DJF"


seasonal_summary <- do.call(rbind, lapply(c("DJF","MAM","JJA","SON"), function(s) {
  sub_ts <- activations_ts[activations_ts$season == s, ]
  data.frame(
    season         = s,
    n_timesteps    = nrow(sub_ts),
    mean_count     = round(mean(sub_ts$count),    3),
    mean_fraction  = round(mean(sub_ts$fraction,na.rm=T), 4),
    max_count      = max(sub_ts$count),
    n_active_ts    = sum(sub_ts$count > 0),         # timesteps with any activation
    pct_active_ts  = round(mean(sub_ts$count > 0) * 100, 1)  # % of timesteps active
  )
}))
print(seasonal_summary)


##### Quantify the difference due to convective control ###


correction_summary <- do.call(rbind, lapply(seq_along(all_data), function(i) {
  
  df <- all_data[[i]]
  
  # Find grid cells where coef.var exceeds threshold
  active <- !is.na(df$coef.var) & df$coef.var > threshold
  
  # Only proceed if there are any active cells
  if (sum(active) == 0) return(NULL)
  
  # Compute correction only for active cells
  correction <- df$radar[active] - df$radar.orig[active]
  
  data.frame(
    date      = names(all_data)[i],
    n_active  = sum(active),
    min_corr  = round(min(correction,  na.rm = TRUE), 4),
    max_corr  = round(max(correction,  na.rm = TRUE), 4),
    mean_corr = round(mean(correction, na.rm = TRUE), 4),
    abs_mean  = round(mean(abs(correction), na.rm = TRUE), 4)  # magnitude regardless of sign
  )
}))

### Print results ###

cat("── Overall correction statistics ──\n")

cat(sprintf("Total timesteps with activations: %d / %d\n",
            nrow(correction_summary), length(all_data)))

cat(sprintf("Overall mean correction: %.4f\n", mean(correction_summary$mean_corr)))
cat(sprintf("Overall mean abs correction: %.4f\n", mean(correction_summary$abs_mean)))
cat(sprintf("Overall min correction:      %.4f\n", min(correction_summary$min_corr)))
cat(sprintf("Overall max correction:      %.4f\n", max(correction_summary$max_corr)))


true_overall_mean <- weighted.mean(correction_summary$mean_corr,
                                   correction_summary$n_active)
cat(sprintf("True weighted mean correction: %.4f\n", true_overall_mean))


true_overall_abs_mean <- weighted.mean(correction_summary$abs_mean,
                                   correction_summary$n_active)
cat(sprintf("True weighted mean correction: %.4f\n", true_overall_abs_mean))


######## Connection with Kriging variance #####

#### Expected value + Variance timeseries ###

nearest_grid <- function(
  rda_file,
  station_coords,
  mu_min= 0.1
) {

  load(rda_file)

  if (!exists("kriging_crop_list")) stop("kriging_crop_list not found in rda_file")
  if (!exists("variance_crop_list")) stop("variance_crop_list not found in rda_file")

  # ── 2. Get grid coordinate vectors from reference matrix ─────────────────────

  ref <- kriging_crop_list[[1]]
    
  get_xy_from_matrix <- function(mat) {
    x <- attr(mat, "x")
    y <- attr(mat, "y")

    # fallback options if attrs are missing
    if (is.null(x) && !is.null(colnames(mat))) {
      suppressWarnings(x <- as.numeric(colnames(mat)))
    }
    if (is.null(y) && !is.null(rownames(mat))) {
      suppressWarnings(y <- as.numeric(rownames(mat)))
    }

    if (is.null(x) || is.null(y)) {
      stop("Could not find x/y coordinates in matrix attributes or dimnames.")
    }

    list(x = x, y = y)
  }

  xy <- get_xy_from_matrix(ref)
  x_vec <- xy$x
  y_vec <- xy$y

  # ── 3. Time vector ─────────────────────────────────────────────────────────────
  if (exists("timestamps")) {
    time_vec <- as.POSIXct(timestamps, origin = "1970-01-01", tz = "UTC")
  } else {
    time_vec <- seq_along(kriging_crop_list)
  }

  # ── 4. Loop over all stations ──────────────────────────────────────────────────
  
  n_stations <- nrow(station_coords)
  results    <- vector("list", n_stations)
  
  for (i in seq_len(n_stations)) {
    
    target_x <- station_coords[i, 1]
    target_y <- station_coords[i, 2]
    
    # ── 5. Find nearest grid point ───────────────────────────────────────────────
    ix <- which.min(abs(x_vec - target_x))
    iy <- which.min(abs(y_vec - target_y))
    
    nearest_x <- x_vec[ix]
    nearest_y <- y_vec[iy]
    dist      <- sqrt((nearest_x - target_x)^2 + (nearest_y - target_y)^2)
    
    # ── 6. Extract mu and variance time series at this grid point ────────────────
    mu_raw  <- sapply(kriging_crop_list,  function(m) m[ix, iy])  # note: row=y col=x
    var_raw <- sapply(variance_crop_list, function(v) v[ix, iy])
    
    # ── 7. Apply mu threshold ────────────────────────────────────────────────────
    mu_thr <- mu_raw
    mu_thr[mu_thr < mu_min] <- NA_real_
    
    # ── 8. Compute IQR ───────────────────────────────────────────────────────────
    iqr_vals <- compute_iqr_list(mu_thr, var_raw)
    
    # ── 9. Store result for this station ─────────────────────────────────────────
    results[[i]] <- data.table(
      station_id = i,
      x_station  = target_x,
      y_station  = target_y,
      x_grid     = nearest_x,
      y_grid     = nearest_y,
      dist       = dist,
      time       = time_vec,
      mu         = mu_thr,
      iqr        = iqr_vals
    )
    
    if (i %% 50 == 0) cat(sprintf("Progress: %d / %d stations\n", i, n_stations))
  }

  cat("Combining results...\n")
  t1     <- Sys.time()
  out_dt <- rbindlist(results)
  t2     <- Sys.time()
  cat(sprintf("Done. Rows: %d  Cols: %d  Time: %.1f sec\n", 
              nrow(out_dt), ncol(out_dt), as.numeric(t2 - t1)))
  return(out_dt)

}

# ── Call the function ──────────────────────────────────────────────────────────
load("/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/CPC2300100002.rda")
station_coords <- output[[8]][[1]]@coords

result <- nearest_grid(
  rda_file       = "/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/precip_transformed_results_new_2023.rda",
  station_coords = station_coords,
  mu_min         = 0.1
)




# ── The connection is through: station_id + timestamp ─────────────────────────

# You have:
# activations_ts  → date, count, fraction        (per timestep, all stations aggregated)
# result          → station_id, time, mu, iqr    (per station per timestep)
# all_data        → coef.var per station per timestep

# ── 1. Add coef.var to the result dataframe ────────────────────────────────────
# Build a long format coef.var dataframe
coefvar_long <- do.call(rbind, lapply(seq_along(all_data), function(i) {
  df <- all_data[[i]]
  data.frame(
    station_id = seq_len(nrow(df)),
    time       = dates[i],
    coef.var   = df$coef.var
  )
}))

# ── 2. Merge with result (mu, iqr) ────────────────────────────────────────────
result_full <- merge(result, coefvar_long, by = c("station_id", "time"))

# ── 3. Flag active vs non-active ──────────────────────────────────────────────
result_full$active <- result_full$coef.var > threshold & 
                      !is.na(result_full$coef.var)

# ── 4. Compare mu and iqr between active and non-active ───────────────────────
cat("── mu when convective control ACTIVE ──\n")
print(summary(result_full$mu[result_full$active == TRUE]))

cat("── mu when convective control NOT active ──\n")
print(summary(result_full$mu[result_full$active == FALSE]))

cat("── iqr when convective control ACTIVE ──\n")
print(summary(result_full$iqr[result_full$active == TRUE]))

cat("── iqr when convective control NOT active ──\n")
print(summary(result_full$iqr[result_full$active == FALSE]))

# ── 5. Per station: how does mu differ when active vs not ─────────────────────
station_comparison <- do.call(rbind, lapply(1:275, function(s) {
  sub    <- result_full[result_full$station_id == s, ]
  active <- sub[sub$active == TRUE,  ]
  noact  <- sub[sub$active == FALSE, ]
  data.frame(
    station_id   = s,
    x            = sub$x_station[1],
    y            = sub$y_station[1],
    n_active     = nrow(active),
    mean_mu_active   = mean(active$mu,  na.rm = TRUE),
    mean_mu_inactive = mean(noact$mu,   na.rm = TRUE),
    mean_iqr_active  = mean(active$iqr, na.rm = TRUE),
    mean_iqr_inactive= mean(noact$iqr,  na.rm = TRUE)
  )
}))

print(summary(station_comparison))
