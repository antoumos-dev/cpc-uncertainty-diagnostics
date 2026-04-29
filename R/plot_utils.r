
dim.x=710
dim.y=640

swiss.corners=c(485,75,835,295)
radar.corners <- c(255,-160,965,480)
h=75
INCA.corners=c(swiss.corners[1]-h,swiss.corners[2]-h,swiss.corners[3]+h,swiss.corners[4]+h)
x.INCA.1 = INCA.corners[1] - radar.corners[1]                                           
x.INCA.2 = x.INCA.1 + (INCA.corners[3]-INCA.corners[1])                                               
y.INCA.1 = INCA.corners[2] - radar.corners[2]                                           
y.INCA.2 = y.INCA.1 + (INCA.corners[4]-INCA.corners[2])                                               
n = (x.INCA.2-x.INCA.1+1) * (y.INCA.2-y.INCA.1+1)

INCA.corners <- c(410,0,910,370)

# Set dimensions
dim.x <- 710
dim.y <- 640
input.shape <- "/users/antoumos/SHAPE/"
input.shape.file <- "CHE_adm0"
input.dem.file <- "/users/antoumos/ccs4.png"

# Load the Switzerland shapefile
Switzerland <- readOGR(input.shape, input.shape.file, verbose = FALSE)

t <- geocors.trafo(x = Switzerland@polygons[[1]]@Polygons[[1]]@coords[,1],
                  y = Switzerland@polygons[[1]]@Polygons[[1]]@coords[,2],
                  from.type = "lonlat", to.type = "swisscors")

# Transform coordinates
t[[1]] <- t[[1]] / 1000
t[[2]] <- t[[2]] / 1000
tt <- matrix(c(t[[1]], t[[2]]), ncol = 2)
Switzerland@polygons[[1]]@Polygons[[1]]@coords <- tt
Switzerland@bbox[[1]] <-  255
Switzerland@bbox[[2]] <- -165
Switzerland@bbox[[3]] <-  255 + dim.x
Switzerland@bbox[[4]] <- -160 + dim.y
proj4string(Switzerland) <- ""

# Load DEM data
dem.m <- readPNG(input.dem.file)

# Adjust dimensions if necessary
dem.m <- t(dem.m)[,dim(dem.m)[[1]]:1]		
dem.m <- assignCoordsToPrecip(dem.m)


		
color.ramp=c( "#FFFFFF00"
				 ,"#640064","#AF00AF","#DC00DC","#3232C8","#0064FF"
				 ,"#009696","#00C832","#64FF00","#96FF00","#C8FF00"
				 ,"#FFFF00","#FFC800","#FFA000","#FF7D00","#E11900","#000000" )
				 				
color.levels=c( 0.00, 0.16, 0.25, 0.40, 0.63								
				   ,1.00, 1.60, 2.50, 4.00, 6.30								
				   ,10.0, 16.0, 25.0, 40.0, 63.0								
				   ,100., 160., 250.)	



## rotate ##
dem.fixed <- t(dem.m[nrow(dem.m):1, ])
dem.fixed <- dem.m[nrow(dem.m):1, ncol(dem.m):1]

### do this 3 times?
dem.fixed <- t(dem.fixed[nrow(dem.fixed):1, ])
dem.fixed <- dem.fixed[nrow(dem.fixed):1, ncol(dem.fixed):1]
			

# ── background-only function ───────────────────────────────────────────────
topography_background <- function() {
  plot.image(
    x = dem.m,
    time.stamp = NULL,
    color.ramp = color.ramp,
    color.levels = color.levels,
    dem.m = dem.m,
    Switzerland = Switzerland,
    scale.title = "Topography",
    mlayout = TRUE,
    plot.dem = TRUE,
    plot.color.scale = TRUE,
    border.col = "black"
  )
}

# --------------------------
# Swiss border loader
# --------------------------
load.map.elements <- function(
  xlim = c(0, 710),
  ylim = c(0, 640),
  input.shape = "/users/antoumos/SHAPE/",
  input.shape.file = "CHE_adm0"
) {
  library(sp)
  library(rgdal)
  library(raster)

  Switzerland <- readOGR(dsn = input.shape, layer = input.shape.file, verbose = FALSE)

  t <- geocors.trafo(
    x = Switzerland@polygons[[1]]@Polygons[[1]]@coords[, 1],
    y = Switzerland@polygons[[1]]@Polygons[[1]]@coords[, 2],
    from.type = "lonlat",
    to.type   = "swisscors"
  )

  Switzerland@polygons[[1]]@Polygons[[1]]@coords <- cbind(t[[1]] / 1000, t[[2]] / 1000)

  Switzerland <- raster::crop(
    Switzerland,
    raster::extent(xlim[1], xlim[2], ylim[1], ylim[2])
  )

  return(Switzerland)
}

# --------------------------
# Spatial field plotter
# --------------------------
plot_cropped_field <- function(
  Z,
  variance_list,
  xlim, ylim,
  vmin = NULL, vmax = NULL,
  cap_quant = 0.99,
  palette_end = 0.99,
  n_colors = 150,
  pal = NULL,
  title = "Cropped Field",
  xlab = "Swiss easting (km)",
  ylab = "Swiss northing (km)",
  output_file = NULL
) {
  xs_full <- attr(variance_list[[1]], "x")
  ys_full <- attr(variance_list[[1]], "y")

  x_idx <- which(xs_full >= xlim[1] & xs_full <= xlim[2])
  y_idx <- which(ys_full >= ylim[1] & ys_full <= ylim[2])

  full_rows <- length(xs_full)
  full_cols <- length(ys_full)

  if (nrow(Z) == full_rows && ncol(Z) == full_cols) {
    Z_crop  <- Z[x_idx, y_idx]
    xs_crop <- xs_full[x_idx]
    ys_crop <- ys_full[y_idx]
  } else {
    Z_crop  <- Z
    xs_crop <- seq(xlim[1], xlim[2], length.out = nrow(Z))
    ys_crop <- seq(ylim[1], ylim[2], length.out = ncol(Z))
  }

  if (is.null(vmin)) vmin <- min(Z_crop, na.rm = TRUE)
  if (is.null(vmax)) vmax <- quantile(Z_crop, cap_quant, na.rm = TRUE)

  Zplot <- pmin(pmax(Z_crop, vmin), vmax)

  cat(sprintf("%s | crop range: [%.3f, %.3f] | scale: [%.3f, %.3f]\n",
              title, min(Z_crop, na.rm = TRUE), max(Z_crop, na.rm = TRUE), vmin, vmax))

  if (is.null(pal)) pal <- viridisLite::viridis(n_colors, end = palette_end)
  swiss_border <- load.map.elements(xlim = xlim, ylim = ylim)

  if (!is.null(output_file)) {
    png(output_file, width = 1600, height = 1400, res = 150)
    on.exit(dev.off(), add = TRUE)
  }

  x_buffer  <- 0.09 * diff(range(xs_crop))
  xlim_plot <- c(min(xs_crop), max(xs_crop) + x_buffer)
  ylim_plot <- range(ys_crop)

  plot(NULL, xlim = xlim_plot, ylim = ylim_plot, asp = 1,
       xlab = xlab, ylab = ylab, main = title, xaxs = "i", yaxs = "i")

  image(xs_crop, ys_crop, Zplot, add = TRUE, useRaster = TRUE,
        col = pal, zlim = c(vmin, vmax))

  lx0 <- max(xs_crop) + 5
  lx1 <- max(xs_crop) + 15
  ly  <- seq(min(ys_crop), max(ys_crop), length.out = n_colors + 1)

  for (i in seq_len(n_colors))
    rect(lx0, ly[i], lx1, ly[i + 1], col = pal[i], border = NA)

  rect(lx0, min(ys_crop), lx1, max(ys_crop), border = "grey30", lwd = 0.8)

  lab_vals <- pretty(c(vmin, vmax), n = 5)
  lab_pos  <- min(ys_crop) + (lab_vals - vmin) / (vmax - vmin) * diff(range(ys_crop))
  text(lx1 + 3, lab_pos, labels = formatC(lab_vals, digits = 3, format = "fg"),
       adj = 0, cex = 0.9)

  plot(swiss_border, add = TRUE, border = "white", lwd = 2)

  invisible(list(Zcrop = Z_crop, xs = xs_crop, ys = ys_crop, vmin = vmin, vmax = vmax))
}

# --------------------------
# Shared color scale helper
# --------------------------
get_common_color_scale <- function(fields, probs = NULL) {
  vals <- unlist(lapply(fields, c), use.names = FALSE)
  vals <- vals[is.finite(vals)]

  if (!length(vals))
    return(list(vmin = NA_real_, vmax = NA_real_, n = 0))

  if (is.null(probs)) {
    vmin <- min(vals)
    vmax <- max(vals)
  } else {
    qs   <- quantile(vals, probs = probs, na.rm = TRUE)
    vmin <- qs[1]
    vmax <- qs[2]
  }

  list(vmin = vmin, vmax = vmax, n = length(vals))
}

# --------------------------
# Interannual SD plotter
# --------------------------
plot_interannual_sd_products <- function(res_thr, out_dir,
                                         xlim = c(480, 840), ylim = c(60, 300),
                                         cap_quant = 0.99, palette_end = 0.99) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  thr_txt <- gsub("\\.", "p", sprintf("%.2f", res_thr$threshold))
  yr_min  <- min(res_thr$years)
  yr_max  <- max(res_thr$years)

  plot_cropped_field(res_thr$interannual_sd$mean_mu$sd,
                     res_thr$variance_ref_list, xlim, ylim,
                     title = sprintf("Interannual SD μ (%s–%s), thr=%.2f", yr_min, yr_max, res_thr$threshold),
                     cap_quant = cap_quant, palette_end = palette_end,
                     output_file = file.path(out_dir, sprintf("INTERANNUAL_SD_MU_thr_%s.png", thr_txt)))

  plot_cropped_field(res_thr$interannual_sd$mean_iqr$sd,
                     res_thr$variance_ref_list, xlim, ylim,
                     title = sprintf("Interannual SD IQR(90-10), thr=%.2f", res_thr$threshold),
                     cap_quant = cap_quant, palette_end = palette_end,
                     output_file = file.path(out_dir, sprintf("INTERANNUAL_SD_IQR_thr_%s.png", thr_txt)))

  plot_cropped_field(res_thr$interannual_sd$accum_mu$sd,
                     res_thr$variance_ref_list, xlim, ylim,
                     title = sprintf("Interannual SD accumulation (mm), thr=%.2f", res_thr$threshold),
                     cap_quant = cap_quant, palette_end = palette_end,
                     output_file = file.path(out_dir, sprintf("INTERANNUAL_SD_ACCUM_thr_%s.png", thr_txt)))

  plot_cropped_field(res_thr$interannual_sd$wet_hours$sd,
                     res_thr$variance_ref_list, xlim, ylim,
                     title = sprintf("Interannual SD wet hours, thr=%.2f", res_thr$threshold),
                     cap_quant = cap_quant, palette_end = palette_end,
                     output_file = file.path(out_dir, sprintf("INTERANNUAL_SD_WETHOURS_thr_%s.png", thr_txt)))

  for (s in names(res_thr$interseasonal_sd)) {
    ss_sd <- res_thr$interseasonal_sd[[s]]

    plot_cropped_field(ss_sd$mean_mu$sd,  res_thr$variance_ref_list, xlim, ylim,
                       title = sprintf("%s interannual SD μ, thr=%.2f", s, res_thr$threshold),
                       cap_quant = cap_quant, palette_end = palette_end,
                       output_file = file.path(out_dir, sprintf("INTERSEASON_%s_SD_MU_thr_%s.png", s, thr_txt)))

    plot_cropped_field(ss_sd$mean_iqr$sd, res_thr$variance_ref_list, xlim, ylim,
                       title = sprintf("%s interannual SD IQR, thr=%.2f", s, res_thr$threshold),
                       cap_quant = cap_quant, palette_end = palette_end,
                       output_file = file.path(out_dir, sprintf("INTERSEASON_%s_SD_IQR_thr_%s.png", s, thr_txt)))

    plot_cropped_field(ss_sd$accum_mu$sd, res_thr$variance_ref_list, xlim, ylim,
                       title = sprintf("%s interannual SD accumulation, thr=%.2f", s, res_thr$threshold),
                       cap_quant = cap_quant, palette_end = palette_end,
                       output_file = file.path(out_dir, sprintf("INTERSEASON_%s_SD_ACCUM_thr_%s.png", s, thr_txt)))

    plot_cropped_field(ss_sd$wet_hours$sd, res_thr$variance_ref_list, xlim, ylim,
                       title = sprintf("%s interannual SD wet hours, thr=%.2f", s, res_thr$threshold),
                       cap_quant = cap_quant, palette_end = palette_end,
                       output_file = file.path(out_dir, sprintf("INTERSEASON_%s_SD_WETHOURS_thr_%s.png", s, thr_txt)))
  }

  invisible(TRUE)
}