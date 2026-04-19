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



setwd("/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/conv_control_off/")

project_root <- "/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats"
utils_file   <- file.path(project_root, "R", "utils.r")
plot_utils_file <- file.path(project_root,"R","plot_utils.r")
source(utils_file) # helpers function
source(plot_utils_file) # plot helpers function


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

saveRDS(all_data, file = "../new_project/out_stats/all_data_conv_off_2023.rds")

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

## nearest_grid() and get_xy_from_matrix() are defined in R/utils.r

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

# Add coef.var to the result dataframe ────────────────────────────────────
# Build a long format coef.var dataframe
coefvar_long <- rbindlist(lapply(seq_along(all_data), function(i) {
  df <- all_data[[i]]
  data.table(
    station_id = seq_len(nrow(df)),
    time       = dates[i],
    coef.var   = df$coef.var
  )
}))

#  Merge with result (mu, iqr) ────────────────────────────────────────────
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
station_comparison <- result_full[!is.na(mu), .(
  x                 = x_station[1],
  y                 = y_station[1],
  n_active          = sum(active == TRUE,  na.rm = TRUE),
  n_inactive        = sum(active == FALSE, na.rm = TRUE),
  mean_mu_active    = mean(mu[active == TRUE],   na.rm = TRUE),
  mean_mu_inactive  = mean(mu[active == FALSE],  na.rm = TRUE),
  mean_iqr_active   = mean(iqr[active == TRUE],  na.rm = TRUE),
  mean_iqr_inactive = mean(iqr[active == FALSE], na.rm = TRUE)
), by = station_id]

print(summary(station_comparison))

################################  spatial representation   ###################################################

# ── Prepare data: one row per station with total activations ──────────────────

station_activations <- station_comparison[, .(
  station_id,
  x,
  y,
  n_active,
  n_inactive,
  activation_rate = n_active/ (n_active + n_inactive) ,# fraction of wet hours active

 # ── intensity difference ──────────────────────────────────────────────────
  delta_mu         = mean_mu_active  - mean_mu_inactive,   # absolute difference
  delta_iqr        = mean_iqr_active - mean_iqr_inactive,
  ratio_mu         = mean_mu_active  / mean_mu_inactive,   # relative difference
  ratio_iqr        = mean_iqr_active / mean_iqr_inactive,

  # ── uncertainty relative to intensity ─────────────────────────────────────
  cv_active        = mean_iqr_active  / mean_mu_active,    # IQR/mu when active
  cv_inactive      = mean_iqr_inactive / mean_mu_inactive  # IQR/mu when inactive

)]

#### higher intensity correlated with higher uncertainty? ####
r <- cor(station_activations$ratio_mu,
         station_activations$ratio_iqr,
         use    = "complete.obs",
         method = "spearman")
cat(sprintf("Spearman r = %.3f\n", r))

## Significance testing - IDR/mu ##

wilcox.test(station_activations$cv_active, station_activations$cv_inactive, conf.int = TRUE )


# ── Plot ───────────────────────────────────────────────────────────────────────
png(filename = "topography_activations_2023.png", 
    width = 650, height = 650, units = "px", pointsize = 18)

plot(1, type = "n", xlim = c(320, 900), ylim = c(-140, 450),
    xlab = "Swiss easting (km)", ylab = "Swiss northing (km)",
     main = "Convective Control Activations 2023", asp = 1)

# Add topography
# NOTE: dem.fixed must be loaded beforehand (e.g. from plot_utils.r or a separate load())
rasterImage(dem.fixed, 255, -165, 965, 480)

# Add Switzerland border
# NOTE: Switzerland (SpatialPolygons) must be loaded beforehand
plot(Switzerland, add = TRUE, border = "white", lwd = 2.2)

# ── Color scale ────────────────────────────────────────────────────────────────
n_colors   <- 100
colormap   <- heat.colors(n_colors)
val_range  <- range(station_activations$n_active, na.rm = TRUE)
color_idx  <- cut(station_activations$n_active, 
                  breaks = seq(val_range[1], val_range[2], length.out = n_colors + 1),
                  labels = FALSE, include.lowest = TRUE)
color_idx[is.na(color_idx)] <- 1

# Add station points
points(station_activations$x, station_activations$y,
       pch = 21, bg = colormap[color_idx],
       cex = 0.8)

# ── Legend ─────────────────────────────────────────────────────────────────────
legend_vals   <- pretty(val_range, n = 5)
legend_colors <- colormap[cut(legend_vals, 
                              breaks = seq(val_range[1], val_range[2], 
                                           length.out = n_colors + 1),
                              labels = FALSE, include.lowest = TRUE)]

x_left   <- 310;  x_right  <- 320
y_bottom <- 100;  y_top    <- 300
legend_height <- y_top - y_bottom
step_height   <- legend_height / length(legend_colors)

for (i in seq_along(legend_colors)) {
  rect(x_left,  y_bottom + (i - 1) * step_height,
       x_right, y_bottom + i       * step_height,
       col = legend_colors[i], border = NA)
}

label_pos <- seq(y_bottom, y_top, length.out = length(legend_vals))
text(x = x_right + 5, y = label_pos, 
     labels = round(legend_vals, 2), pos = 4, cex = 0.6)
text(x = x_left + 5,  y = y_top + 15, 
     labels = "N activations", pos = 3, cex = 0.7)

dev.off()




########## log-log regression ############
############# Does intensity increases faster than uncertainty during convection ########

# ── 1. Prepare data ────────────────────────────────────────────────────────────
# Use filtered stations (n_active >= 10) to remove unreliable estimates
reg_data <- station_comparison[ 
                                 !is.na(station_comparison$mean_mu_active) & 
                                 !is.na(station_comparison$mean_iqr_active)]

# Log transform
reg_data[, log_mu  := log(station_comparison$mean_mu_active)]
reg_data[, log_iqr := log(station_comparison$mean_iqr_active)]

fit_inactive <- lm(log(station_comparison$mean_iqr_inactive) ~ log(station_comparison$mean_mu_inactive), 
                   data = reg_data)


# ── 2. Log-log regression ──────────────────────────────────────────────────────
fit <- lm(log_iqr ~ log_mu, data = reg_data)
summary(fit)


# Log-log regression: slope = elasticity
fit <- lm(log(mean_iqr_active) ~ log(mean_mu_active), 
          data = station_activations_filtered)
summary(fit)


# Stack active and inactive together with a flag
reg_combined <- data.table(
  log_mu  = c(log(reg_data$mean_mu_active), 
              log(station_comparison$mean_mu_inactive)),
  log_iqr = c(log(reg_data$mean_iqr_active),  
              log(station_comparison$mean_iqr_inactive)),
  active  = c(rep(1, nrow(reg_data)), 
              rep(0, nrow(reg_data)))
)

# Does the IQR/mu relationship differ between active and inactive?
fit_combined <- lm(log_iqr ~ log_mu * active, data = reg_combined)
summary(fit_combined)


######### Plot ##########

# ── Data for regression lines ──────────────────────────────────────────────────
mu_range <- seq(min(reg_combined$log_mu), 
                max(reg_combined$log_mu), 
                length.out = 100)

# Predicted lines from fit_combined
pred_inactive <- coef(fit_combined)["(Intercept)"] + 
                 coef(fit_combined)["log_mu"] * mu_range

pred_active   <- (coef(fit_combined)["(Intercept)"] + coef(fit_combined)["active"]) + 
                 (coef(fit_combined)["log_mu"]       + coef(fit_combined)["log_mu:active"]) * mu_range

# ── Plot ───────────────────────────────────────────────────────────────────────
png("loglog_active_vs_inactive.png", width = 650, height = 600, 
    units = "px", pointsize = 14)

plot(NULL, 
     xlim = range(reg_combined$log_mu),
     ylim = range(reg_combined$log_iqr),
     xlab = "log(mu)  —  log precipitation intensity",
     ylab = "log(IQR)  —  log uncertainty",
     main = "IQR vs mu scaling: active vs inactive")

# ── Points ─────────────────────────────────────────────────────────────────────
points(reg_combined$log_mu[reg_combined$active == 0],
       reg_combined$log_iqr[reg_combined$active == 0],
       pch = 21, bg = "steelblue", col = "steelblue", 
       cex = 0.7, alpha = 0.5)

points(reg_combined$log_mu[reg_combined$active == 1],
       reg_combined$log_iqr[reg_combined$active == 1],
       pch = 21, bg = "tomato", col = "tomato", 
       cex = 0.7)

# ── Regression lines ───────────────────────────────────────────────────────────
lines(mu_range, pred_inactive, col = "steelblue", lwd = 2.5)
lines(mu_range, pred_active,   col = "tomato",    lwd = 2.5)

# ── Reference line slope = 1 ──────────────────────────────────────────────────
abline(a = coef(fit_combined)["(Intercept)"], b = 1, 
       col = "gray50", lwd = 1.5, lty = 2)

# ── Intercept shift annotation ────────────────────────────────────────────────
# Draw vertical arrow showing the -0.253 shift at a mid-range mu value
x_ann  <- median(mu_range)
y_top  <- coef(fit_combined)["(Intercept)"] + 
          coef(fit_combined)["log_mu"] * x_ann
y_bot  <- (coef(fit_combined)["(Intercept)"] + coef(fit_combined)["active"]) + 
          (coef(fit_combined)["log_mu"] + coef(fit_combined)["log_mu:active"]) * x_ann

arrows(x_ann + 0.05, y_top, x_ann + 0.05, y_bot,
       length = 0.08, code = 3, col = "black", lwd = 1.5)
text(x_ann + 0.12, (y_top + y_bot) / 2, 
     labels = sprintf("Δ = %.3f", coef(fit_combined)["active"]),
     cex = 0.8, adj = 0)

# ── Legend ─────────────────────────────────────────────────────────────────────
legend("topleft",
       legend = c(
         sprintf("inactive  ε=%.3f", coef(fit_combined)["log_mu"]),
         sprintf("active    ε=%.3f", coef(fit_combined)["log_mu"] + 
                                     coef(fit_combined)["log_mu:active"]),
         "reference ε=1"
       ),
       col = c("steelblue", "tomato", "gray50"),
       lwd = c(2.5, 2.5, 1.5),
       lty = c(1, 1, 2),
       pch = c(21, 21, NA),
       pt.bg = c("steelblue", "tomato", NA),
       cex  = 0.85)

dev.off()


