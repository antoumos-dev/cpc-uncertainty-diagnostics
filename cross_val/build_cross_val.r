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
library(gstat)                                                      
library(rhdf5)
library(animation)
library(ncdf4)                     
library(lubridate)
library(dplyr)
library(data.table)



args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: Rscript build_cross_val.r <YY>  (e.g. 23)")
YEAR <- sprintf("%02d", as.integer(args[1]))

path_on      <- "/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/"
path_off     <- "/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/conv_control_off/"
project_root <- "/store_new/mch/msclim/antoumos/R/develop/CPC/new_project/out_stats"

parse_date_from_filename <- function(fname) {
  base <- sub("\\.rda$", "", basename(fname))
  yy   <- substr(base, 4, 5)
  doy  <- substr(base, 6, 8)
  hh   <- substr(base, 9, 10)
  as.POSIXct(strptime(paste0("20", yy, " ", doy, " ", hh),
                      format = "%Y %j %H", tz = "UTC"))
}

cat(sprintf("===== Processing year 20%s =====\n", YEAR))

rda_files     <- list.files(path_on,  pattern = paste0("CPC", YEAR, ".*\\.rda$"), full.names = TRUE)
rda_files_off <- list.files(path_off, pattern = paste0("CPC", YEAR, ".*\\.rda$"), full.names = TRUE)

if (length(rda_files)     == 0) stop(sprintf("No ON  files found for 20%s", YEAR))
if (length(rda_files_off) == 0) stop(sprintf("No OFF files found for 20%s", YEAR))

dates     <- as.POSIXct(sapply(rda_files,     parse_date_from_filename), origin = "1970-01-01", tz = "UTC")
dates_off <- as.POSIXct(sapply(rda_files_off, parse_date_from_filename), origin = "1970-01-01", tz = "UTC")

# # conv control ON
# cross_val_active <- vector("list", length(rda_files))
# for (i in seq_along(rda_files)) {
#   load(rda_files[i])
#   cross_val_active[[i]] <- output[[6]]
#   if (i %% 100 == 0) cat(sprintf("ON  Progress: %d / %d\n", i, length(rda_files)))
# }
# names(cross_val_active) <- format(dates, "%Y-%m-%d %H:%M")
# saveRDS(cross_val_active, file = file.path(project_root, paste0("cross_val_active_20", YEAR, ".rds")))

# conv control OFF
cross_val_inactive <- vector("list", length(rda_files_off))
for (i in seq_along(rda_files_off)) {
  f <- rda_files_off[i]
  if (file.size(f) == 0) {
    cat(sprintf("Skipping zero-byte file: %s\n", basename(f)))
    next
  }
  ok <- withCallingHandlers(
    tryCatch({ load(f); TRUE }, error = function(e) {
      cat(sprintf("ERROR loading %s: %s\n", basename(f), conditionMessage(e)))
      FALSE
    }),
    warning = function(w) {
      cat(sprintf("WARN  loading %s: %s\n", basename(f), conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  if (!ok) next
  cross_val_inactive[[i]] <- output[[6]]
  if (i %% 100 == 0) cat(sprintf("OFF Progress: %d / %d\n", i, length(rda_files_off)))
}
names(cross_val_inactive) <- format(dates_off, "%Y-%m-%d %H:%M")
saveRDS(cross_val_inactive, file = file.path(project_root, paste0("cross_val_inactive_20", YEAR, ".rds")))
cat(sprintf("Saved cross_val_inactive_20%s\n", YEAR))
