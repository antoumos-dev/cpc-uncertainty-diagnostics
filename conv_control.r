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


setwd("/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/")


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
    mean_fraction  = round(mean(sub_ts$fraction), 4),
    max_count      = max(sub_ts$count),
    n_active_ts    = sum(sub_ts$count > 0),         # timesteps with any activation
    pct_active_ts  = round(mean(sub_ts$count > 0) * 100, 1)  # % of timesteps active
  )
}))
print(seasonal_summary)


##### Quantify the difference due to convective control ###













