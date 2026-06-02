
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


# Load helper functions
source("/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats/R/utils.r")

### load data ####

year <- 2023

rda_file <- sprintf(
  "/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/precip_transformed_results_new_%s.rda",
  year
)
load(rda_file)  # expects: kriging_crop_list, variance_crop_list, timestamps


mu_list  <- lapply(kriging_crop_list, strip_attrs_keep_dim)
var_list <- lapply(variance_crop_list, strip_attrs_keep_dim)

iqr_list <- compute_iqr_list(mu_list, var_list)


## Define bins ##

bins <- c(0.1, 1, 5, 10, Inf)
bin_labels <- c("0.1-1", "1-5", "5-10", ">10")

## pool pixel hours ##

mu_all  <- unlist(lapply(mu_list, as.vector), use.names = FALSE)
iqr_all <- unlist(lapply(iqr_list, as.vector), use.names = FALSE)

ok <- is.finite(mu_all) & is.finite(iqr_all) & mu_all >= 0.1
mu_all  <- mu_all[ok]
iqr_all <- iqr_all[ok]

#### Assign the bins ###

bin_id <- cut(mu_all,
              breaks = c(0.1, 1, 5, 10, Inf),
              labels = bin_labels,
              right = FALSE,
              include.lowest = TRUE)


result <- data.frame(
  bin = bin_labels,
  wet_pixel_hours = as.integer(table(bin_id)),
  mean_iqr = tapply(iqr_all, bin_id, mean, na.rm = TRUE),
  median_iqr = tapply(iqr_all, bin_id, median, na.rm = TRUE)
)

print(result)

#### Relative uncertainty per bin ####

rel_all <- iqr_all / mu_all
rel_all[!is.finite(rel_all)] <- NA_real_

result$mean_rel <- tapply(rel_all, bin_id, mean, na.rm = TRUE)
result$median_rel <- tapply(rel_all, bin_id, median, na.rm = TRUE)






