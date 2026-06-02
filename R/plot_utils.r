
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
# Shared diff plotter — used by both plot_diff_stats() and plot_mean_diff_products()
# annual_off / annual_on : named lists with keys mean_mu, mean_iqr, accum_mu, wet_hours, rel_uncert_Bmed
# seasonal_off / seasonal_on : named list of seasons, each a named list of the same fields
# --------------------------
plot_diff_fields <- function(annual_off, annual_on,
                              seasonal_off, seasonal_on,
                              variance_ref_list,
                              label, threshold,
                              out_dir,
                              prefix      = "DIFF",
                              xlim        = c(480, 840),
                              ylim        = c(60, 300),
                              cap_quant   = 0.99,
                              include_mad = FALSE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  thr_txt  <- gsub("\\.", "p", sprintf("%.2f", threshold))
  diff_pal <- colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(150)
  mad_pal  <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(100)

  plot_diff <- function(Z_off, Z_on, fname, ttl) {
    if (is.null(Z_off) || is.null(Z_on)) return(invisible(FALSE))
    dZ   <- Z_off - Z_on
    vmax <- max(quantile(abs(dZ), cap_quant, na.rm = TRUE), 1e-9)
    plot_cropped_field(
      Z             = dZ,
      variance_list = variance_ref_list,
      xlim = xlim, ylim = ylim,
      vmin = -vmax, vmax = vmax,
      pal  = diff_pal,
      title = ttl,
      output_file = file.path(out_dir, fname)
    )
    invisible(TRUE)
  }

  plot_rel_diff <- function(Z_off, Z_on, fname, ttl) {
    if (is.null(Z_off) || is.null(Z_on)) return(invisible(FALSE))
    dZ_rel <- ifelse(abs(Z_on) < 1e-6, NA, (Z_off - Z_on) / Z_on * 100)
    vmax <- max(quantile(abs(dZ_rel), cap_quant, na.rm = TRUE), 1e-9)
    plot_cropped_field(
      Z             = dZ_rel,
      variance_list = variance_ref_list,
      xlim = xlim, ylim = ylim,
      vmin = -vmax, vmax = vmax,
      pal  = diff_pal,
      title = ttl,
      output_file = file.path(out_dir, fname)
    )
    invisible(TRUE)
  }

  plot_mad <- function(Z_off, Z_on, fname, ttl) {
    if (is.null(Z_off) || is.null(Z_on)) return(invisible(FALSE))
    plot_cropped_field(
      Z             = abs(Z_off - Z_on),
      variance_list = variance_ref_list,
      xlim = xlim, ylim = ylim,
      pal  = mad_pal,
      cap_quant = cap_quant,
      title = ttl,
      output_file = file.path(out_dir, fname)
    )
    invisible(TRUE)
  }

  fields <- list(
    list(key = "mean_mu",         label = "mean μ",          unit = "mm/h"),
    list(key = "mean_iqr",        label = "mean IQR(90-10)", unit = "mm/h"),
    list(key = "accum_mu",        label = "accumulation",    unit = "mm"),
    list(key = "wet_hours",       label = "wet hours",       unit = "h"),
    list(key = "rel_uncert_Bmed", label = "rel. uncertainty",unit = "")
  )

  for (f in fields) {
    unit_str <- if (nchar(f$unit) > 0) paste0(" (", f$unit, ")") else ""
    plot_diff(
      annual_off[[f$key]], annual_on[[f$key]],
      fname = sprintf("%s_ANNUAL_%s_thr_%s.png", prefix, toupper(f$key), thr_txt),
      ttl   = sprintf("%s annual Δ%s%s (off − on), thr=%.2f",
                      label, f$label, unit_str, threshold)
    )
    if (f$key != "rel_uncert_Bmed")
      plot_rel_diff(
        annual_off[[f$key]], annual_on[[f$key]],
        fname = sprintf("%s_ANNUAL_REL_%s_thr_%s.png", prefix, toupper(f$key), thr_txt),
        ttl   = sprintf("%s annual %s (off − on) / on [%%], thr=%.2f",
                        label, f$label, threshold)
      )
    if (include_mad)
      plot_mad(
        annual_off[[f$key]], annual_on[[f$key]],
        fname = sprintf("%s_ANNUAL_MAD_%s_thr_%s.png", prefix, toupper(f$key), thr_txt),
        ttl   = sprintf("%s annual MAD |off-base| %s%s, thr=%.2f",
                        label, f$label, unit_str, threshold)
      )
  }

  for (s in names(seasonal_off)) {
    if (!s %in% names(seasonal_on)) next
    for (f in fields) {
      unit_str <- if (nchar(f$unit) > 0) paste0(" (", f$unit, ")") else ""
      plot_diff(
        seasonal_off[[s]][[f$key]], seasonal_on[[s]][[f$key]],
        fname = sprintf("%s_SEASON_%s_%s_thr_%s.png", prefix, s, toupper(f$key), thr_txt),
        ttl   = sprintf("%s %s %s%s (off − on), thr=%.2f",
                        label, s, f$label, unit_str, threshold)
      )
      if (f$key != "rel_uncert_Bmed")
        plot_rel_diff(
          seasonal_off[[s]][[f$key]], seasonal_on[[s]][[f$key]],
          fname = sprintf("%s_SEASON_%s_REL_%s_thr_%s.png", prefix, s, toupper(f$key), thr_txt),
          ttl   = sprintf("%s %s %s (off − on) / on [%%], thr=%.2f",
                          label, s, f$label, threshold)
        )
      if (include_mad)
        plot_mad(
          seasonal_off[[s]][[f$key]], seasonal_on[[s]][[f$key]],
          fname = sprintf("%s_SEASON_%s_MAD_%s_thr_%s.png", prefix, s, toupper(f$key), thr_txt),
          ttl   = sprintf("%s %s MAD |off-base| %s%s, thr=%.2f",
                          label, s, f$label, unit_str, threshold)
        )
    }
  }

  invisible(TRUE)
}

# --------------------------
# MAD accumulation plotter
# --------------------------
plot_mad_accum_products <- function(res_mad_thr, out_dir,
                                    xlim = c(480, 840), ylim = c(60, 300),
                                    cap_quant = 0.99) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  thr_txt <- gsub("\\.", "p", sprintf("%.2f", res_mad_thr$threshold))
  period  <- paste0(min(res_mad_thr$years), "-", max(res_mad_thr$years))
  mad_pal <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(100)

  plot_one <- function(Z, fname, ttl) {
    if (is.null(Z)) return(invisible(FALSE))
    plot_cropped_field(Z = Z, variance_list = res_mad_thr$variance_ref_list,
                       xlim = xlim, ylim = ylim, pal = mad_pal, cap_quant = cap_quant,
                       title = ttl, output_file = file.path(out_dir, fname))
    invisible(TRUE)
  }

  plot_one(res_mad_thr$annual_mad_accum,
           sprintf("MAD_ANNUAL_ACCUM_thr_%s.png", thr_txt),
           sprintf("%s annual MAD |off-base| accumulation (mm), thr=%.2f", period, res_mad_thr$threshold))

  for (s in names(res_mad_thr$seasonal_mad_accum)) {
    plot_one(res_mad_thr$seasonal_mad_accum[[s]],
             sprintf("MAD_SEASON_%s_ACCUM_thr_%s.png", s, thr_txt),
             sprintf("%s %s MAD |off-base| accumulation (mm), thr=%.2f", period, s, res_mad_thr$threshold))
  }

  invisible(TRUE)
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

# ── Spatial bias helpers (used by spatial_bias_map.r) ─────────────────────────

station_bias <- function(cv_df, min_n_per_year = 30) {
  cv_df %>%
    mutate(year = substr(timestamp, 1, 4)) %>%
    group_by(nat_abbr, x, y, year) %>%
    summarise(bias_yr  = mean(kriging - raingauge, na.rm = TRUE),
              mu_rg_yr = mean(raingauge,           na.rm = TRUE),
              n_yr     = n(),
              .groups  = "drop") %>%
    filter(n_yr >= min_n_per_year) %>%
    group_by(nat_abbr, year) %>%
    summarise(x        = mean(x),
              y        = mean(y),
              bias_yr  = mean(bias_yr),
              mu_rg_yr = mean(mu_rg_yr),
              n_yr     = sum(n_yr),
              .groups  = "drop") %>%
    group_by(nat_abbr) %>%
    summarise(x       = mean(x),
              y       = mean(y),
              bias    = mean(bias_yr),
              mu_rg   = mean(mu_rg_yr),
              n       = sum(n_yr),
              n_years = n(),
              .groups = "drop") %>%
    mutate(rel_bias = bias / mu_rg) %>%
    filter(n_years >= 3)
}

make_colorscale <- function(vals, n_col = 101) {
  lim  <- max(abs(vals), na.rm = TRUE)
  pal  <- colorRampPalette(c("#2166AC", "#F7F7F7", "#D6604D"))(n_col)
  brks <- seq(-lim, lim, length.out = n_col + 1)
  list(lim = lim, pal = pal, brks = brks)
}

to_col <- function(vals, cs) {
  idx <- findInterval(vals, cs$brks, rightmost.closed = TRUE)
  cs$pal[pmax(1L, pmin(idx, length(cs$pal)))]
}

draw_colorbar <- function(cs, unit = "mm/h") {
  par(mar = c(4, 0.5, 3, 4))
  n_col <- length(cs$pal)
  cb_y  <- seq(-cs$lim, cs$lim, length.out = n_col + 1)
  plot(NULL, xlim = c(0, 1), ylim = c(-cs$lim, cs$lim),
       xaxt = "n", xlab = "", ylab = "", main = unit, cex.main = 1.2)
  for (i in seq_len(n_col))
    rect(0, cb_y[i], 1, cb_y[i + 1], col = cs$pal[i], border = NA)
  rect(0, -cs$lim, 1, cs$lim, border = "grey30", lwd = 0.8)
  lab_vals <- pretty(c(-cs$lim, cs$lim), n = 6)
  axis(4, at = lab_vals, labels = formatC(lab_vals, digits = 2, format = "f"),
       las = 1, cex.axis = 1.1)
}

draw_map_panel <- function(sb, cs, title, ch, xlim_ch, ylim_ch) {
  par(mar = c(5, 5, 4, 1))
  plot(NULL, xlim = xlim_ch, ylim = ylim_ch, asp = 1,
       xlab = "Easting (km LV03)", ylab = "Northing (km LV03)", main = title,
       cex.main = 1.6, cex.lab = 1.4, cex.axis = 1.3)
  points(sb$x, sb$y,
         pch = 21, bg = to_col(sb$bias, cs), col = "grey30",
         cex = pmin(1.0 + sqrt(sb$n / 500), 2.5), lwd = 0.3)
  plot(ch, add = TRUE, border = "black", lwd = 1.5)
}