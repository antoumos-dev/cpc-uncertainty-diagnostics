input.shape = "/users/antoumos/SHAPE/"
input.shape.file = "CHE_adm0"

load.map.elements <- function(
  xlim = c(0, 710),
  ylim = c(0, 640),
  input.shape = "/users/antoumos/SHAPE/",
  input.shape.file = "CHE_adm0"
) {
  # --- Required packages ---
  library(sp)
  library(rgdal)
  library(raster)

  # --- Load shapefile (Swiss boundary) ---
  shapefile_path <- file.path(input.shape, paste0(input.shape.file, ".shp"))
  Switzerland <- readOGR(dsn = input.shape, layer = input.shape.file, verbose = FALSE)

  # --- Transform to Swiss coordinates (LV95, km) ---
  t <- geocors.trafo(
    x = Switzerland@polygons[[1]]@Polygons[[1]]@coords[, 1],
    y = Switzerland@polygons[[1]]@Polygons[[1]]@coords[, 2],
    from.type = "lonlat",
    to.type   = "swisscors"
  )

  Switzerland@polygons[[1]]@Polygons[[1]]@coords <- cbind(t[[1]] / 1000, t[[2]] / 1000)

  # --- Crop to selected domain ---
  Switzerland <- raster::crop(
    Switzerland,
    raster::extent(xlim[1], xlim[2], ylim[1], ylim[2])
  )

  # --- Return shapefile only ---
  return(Switzerland)
}


plot_cropped_field <- function(
  Z,
  variance_list,       # only used for coordinates
  xlim, ylim,
  vmin = NULL, vmax = NULL,
  cap_quant = 0.99,
  palette_end = 0.99,
  n_colors = 150,
  main_title = "Cropped Field",
  xlab = "Swiss easting (km)",
  ylab = "Swiss northing (km)",
  output_file = NULL
) {
  xs_full <- attr(variance_list[[1]], "x")
  ys_full <- attr(variance_list[[1]], "y")

  # ---------------------------------------------------------
  # Indices for the crop window in the FULL grid
  # ---------------------------------------------------------
  x_idx <- which(xs_full >= xlim[1] & xs_full <= xlim[2])
  y_idx <- which(ys_full >= ylim[1] & ys_full <= ylim[2])

  # ---------------------------------------------------------
  # CASE A: Z is full 710×640 → crop
  # CASE B: Z is already cropped → use as-is
  # ---------------------------------------------------------
  full_rows <- length(xs_full)
  full_cols <- length(ys_full)

  if (nrow(Z) == full_rows && ncol(Z) == full_cols) {
    # FULL DOMAIN → crop (rows=y, cols=x)
    Z_crop  <- Z[x_idx, y_idx]
    xs_crop <- xs_full[x_idx]
    ys_crop <- ys_full[y_idx]
  } else {
    # ALREADY CROPPED
    Z_crop  <- Z
    xs_crop <- seq(xlim[1], xlim[2], length.out = nrow(Z))
    ys_crop <- seq(ylim[1], ylim[2], length.out = ncol(Z))
  }

  # ---------------------------------------------------------
  # Value limits
  # ---------------------------------------------------------
  if (is.null(vmin)) vmin <- min(Z_crop, na.rm = TRUE)
  if (is.null(vmax)) vmax <- quantile(Z_crop, cap_quant, na.rm = TRUE)

  Zplot <- pmin(pmax(Z_crop, vmin), vmax)

  # ---------------------------------------------------------
  # Palette
  # ---------------------------------------------------------
  pal <- viridisLite::viridis(n_colors, end = palette_end)

  # ---------------------------------------------------------
  # Border (optional)
  # ---------------------------------------------------------
  swiss_border <- load.map.elements(xlim = xlim, ylim = ylim)

  # ---------------------------------------------------------
  # Output file
  # ---------------------------------------------------------
  if (!is.null(output_file)) {
    png(output_file, width = 1600, height = 1400, res = 150)
    on.exit(dev.off(), add = TRUE)
  }

  # ---------------------------------------------------------
  # Plot window
  # ---------------------------------------------------------
  x_buffer <- 0.09 * diff(range(xs_crop))
  xlim_plot <- c(min(xs_crop), max(xs_crop) + x_buffer)
  ylim_plot <- range(ys_crop)

  plot(NULL,
       xlim = xlim_plot,
       ylim = ylim_plot,
       asp  = 1,
       xlab = xlab,
       ylab = ylab,
       main = main_title,
       xaxs = "i",
       yaxs = "i")

  # ---------------------------------------------------------
  # Main image
  # ---------------------------------------------------------
  image(xs_crop, ys_crop, Zplot,
        add      = TRUE,
        useRaster = TRUE,
        col      = pal,
        zlim     = c(vmin, vmax))
  # ---------------------------------------------------------
  # Legend
  # ---------------------------------------------------------
  lx0 <- max(xs_crop) + 5
  lx1 <- max(xs_crop) + 15
  ly  <- seq(min(ys_crop), max(ys_crop), length.out = n_colors + 1)

  for (i in seq_len(n_colors)) {
    rect(lx0, ly[i], lx1, ly[i+1], col = pal[i], border = NA)
  }

  rect(lx0, min(ys_crop), lx1, max(ys_crop),
       border = "grey30", lwd = 0.8)

  lab_vals <- pretty(c(vmin, vmax), n = 5)
  lab_pos <- min(ys_crop) +
    (lab_vals - vmin) / (vmax - vmin) * diff(range(ys_crop))

  text(lx1 + 3, lab_pos,
       labels = formatC(lab_vals, digits = 3, format = "fg"),
       adj = 0, cex = 0.9)

  plot(swiss_border, add = TRUE, border = "white", lwd = 2)
  
  # # Draw rectangle
  # rect(
  # xleft   = 680,
  # ybottom = 100,
  # xright  = 730,
  # ytop    = 140,
  # border  = "black",
  # lwd     = 2,
  # lty     = 2,
  # col     = NA
  # )

  invisible(list(
    Zcrop = Z_crop,
    xs = xs_crop,
    ys = ys_crop,
    vmin = vmin,
    vmax = vmax
  ))
}


# --------------------------
# helpers
# --------------------------

#### function to produce a 3D array
stack_mats_to_array <- function(mat_list) {
  stopifnot(length(mat_list) > 0)

  # Basic dimension check
  nr <- nrow(mat_list[[1]])
  nc <- ncol(mat_list[[1]])

  ok <- vapply(mat_list,
               function(m) is.matrix(m) && nrow(m) == nr && ncol(m) == nc,
               logical(1))

  if (!all(ok)) {
    stop("All elements must be matrices with identical dimensions.")
  }

  # Stack along 3rd dimension (time)
  abind::abind(mat_list, along = 3)
}

skew_log_fun <- function(x, min_n = 100, eps = 1e-6) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < min_n) return(NA_real_)

  xlog <- log(x + eps)

  s <- sd(xlog)
  if (!is.finite(s) || s == 0) return(NA_real_)

  m <- mean(xlog)
  mean(((xlog - m) / s)^3)
}


strip_attrs_keep_dim <- function(m) {
  out <- m
  attributes(out) <- list(dim = dim(m))
  out
}

get_season <- function(dates_utc) {
  mon <- as.integer(format(dates_utc, "%m"))
  s <- character(length(mon))
  s[mon %in% c(12, 1, 2)]  <- "DJF"
  s[mon %in% c(3, 4, 5)]   <- "MAM"
  s[mon %in% c(6, 7, 8)]   <- "JJA"
  s[mon %in% c(9, 10, 11)] <- "SON"
  factor(s, levels = c("DJF","MAM","JJA","SON"))
}

aggregate_mean <- function(field_list) {
  ny <- nrow(field_list[[1]])
  nx <- ncol(field_list[[1]])

  s <- matrix(0, ny, nx)
  n <- matrix(0L, ny, nx)

  for (x in field_list) {
    ok <- is.finite(x)
    s[ok] <- s[ok] + x[ok]
    n[ok] <- n[ok] + 1L
  }

  out <- s / pmax(n, 1L)
  out[n == 0] <- NA_real_
  list(mean = out, n_valid = n)
}

aggregate_wet_hours <- function(field_list, min_threshold, max_threshold = Inf) {
  # counts timesteps where min_threshold < field <= max_threshold (per pixel)

  ny <- nrow(field_list[[1]])
  nx <- ncol(field_list[[1]])

  wh <- matrix(0L, ny, nx)

  for (x in field_list) {
    ok <- is.finite(x)
    cond <- ok & x > min_threshold & x <= max_threshold
    wh[cond] <- wh[cond] + 1L
  }

  wh
}

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

### mean across years accumulators ###

acc_init <- function(template_mat) {
  list(
    sum = matrix(0, nrow(template_mat), ncol(template_mat)),
    n   = matrix(0L, nrow(template_mat), ncol(template_mat))
  )
}

acc_add <- function(acc, mat) {
  ok <- is.finite(mat)
  acc$sum[ok] <- acc$sum[ok] + mat[ok]
  acc$n[ok]   <- acc$n[ok] + 1L
  acc
}

acc_mean <- function(acc) {
  out <- acc$sum / pmax(acc$n, 1L)
  out[acc$n == 0] <- NA_real_
  out
}

acc_sd_init <- function(mat) {
  list(
    n    = matrix(0L, nrow(mat), ncol(mat)),
    mean = matrix(0,  nrow(mat), ncol(mat)),
    M2   = matrix(0,  nrow(mat), ncol(mat))
  )
}

#### for interannual/seasonal standard deviation ###

acc_sd_add <- function(acc, mat) {
  ok <- is.finite(mat)
  if (!any(ok)) return(acc)

  n_old <- acc$n[ok]
  n_new <- n_old + 1L

  delta  <- mat[ok] - acc$mean[ok]
  mean_n <- acc$mean[ok] + delta / n_new
  delta2 <- mat[ok] - mean_n

  acc$mean[ok] <- mean_n
  acc$M2[ok]   <- acc$M2[ok] + delta * delta2
  acc$n[ok]    <- n_new

  acc
}

acc_sd_finalize <- function(acc) {
  sd <- matrix(NA_real_, nrow(acc$mean), ncol(acc$mean))
  ok <- acc$n > 1L
  sd[ok] <- sqrt(acc$M2[ok] / (acc$n[ok] - 1L))
  list(mean = acc$mean, sd = sd, n = acc$n)
}

# Temporal median over time per pixel.
# This is expensive; we do it in blocks of pixels to avoid nt*(ny*nx) RAM blowups.
temporal_median_blocked <- function(field_list, block_size = 50000) {
  ny <- nrow(field_list[[1]])
  nx <- ncol(field_list[[1]])
  nt <- length(field_list)
  npx <- ny * nx

  # pre-flatten matrices into vectors (cheap) and store list of vectors
  vec_list <- lapply(field_list, function(m) as.vector(m))

  out <- rep(NA_real_, npx)
  idx_all <- seq_len(npx)

  for (start in seq(1, npx, by = block_size)) {
    end <- min(start + block_size - 1, npx)
    idx <- idx_all[start:end]

    # build block matrix: rows = pixels, cols = time
    # yes, still big, but controlled by block_size
    blk <- do.call(cbind, lapply(vec_list, `[`, idx))

    out[idx] <- matrixStats::rowMedians(blk, na.rm = TRUE)
  }

  matrix(out, ny, nx)
}

temporal_skewness_blocked <- function(field_list, min_n = 100, eps = 1e-6, block_size = 50000) {
  ny <- nrow(field_list[[1]])
  nx <- ncol(field_list[[1]])
  npx <- ny * nx

  vec_list <- lapply(field_list, function(m) as.vector(m))
  out <- rep(NA_real_, npx)
  idx_all <- seq_len(npx)

  for (start in seq(1, npx, by = block_size)) {
    end <- min(start + block_size - 1, npx)
    idx <- idx_all[start:end]

    blk <- do.call(cbind, lapply(vec_list, `[`, idx))  # pixels x time

    out[idx] <- apply(blk, 1, skew_log_fun, min_n = min_n, eps = eps)
  }

  matrix(out, ny, nx)
}

aggregate_rel_uncert_B <- function(mu_list, iqr_list, mu_min = 0.1,
                                  method = c("median","trim","weighted"),
                                  trim = 0.1, min_hours = 0) {
  method <- match.arg(method)

  # per-hour ratio
  rel_list <- Map(function(mu, iqr) {
    r <- iqr / mu
    r[!is.finite(r)] <- NA_real_
    r[mu < mu_min] <- NA_real_
    r
  }, mu_list, iqr_list)

  # stack to 3D array [row, col, time]
  rel_arr <- stack_mats_to_array(rel_list)  # you may already have something like this
  mu_arr  <- stack_mats_to_array(mu_list)

  n_valid <- apply(is.finite(rel_arr), c(1,2), sum)

  out <- switch(method,
    median = apply(rel_arr, c(1,2), median, na.rm = TRUE),
    trim   = apply(rel_arr, c(1,2), function(v) mean(v, trim = trim, na.rm = TRUE)),
    weighted = {
      num <- apply(rel_arr * mu_arr, c(1,2), sum, na.rm = TRUE)
      den <- apply(mu_arr, c(1,2), sum, na.rm = TRUE)
      x <- num / den
      x[!is.finite(x)] <- NA_real_
      x
    }
  )
  
  out[n_valid < min_hours] <- NA_real_
  out
}

aggregate_rel_uncert_iqr <- function(mu_list, iqr_list, mu_min = 0.1, min_hours = 0) {

  rel_list <- Map(function(mu, iqr) {
    r <- iqr / mu
    r[!is.finite(r)] <- NA_real_
    r[mu <= mu_min]  <- NA_real_
    r
  }, mu_list, iqr_list)

  rel_arr <- stack_mats_to_array(rel_list)
  n_valid <- apply(is.finite(rel_arr), c(1,2), sum)

iqr_rel <- apply(rel_arr, c(1,2), function(v) {
  quantile(v, 0.9, na.rm = TRUE) -
  quantile(v, 0.1, na.rm = TRUE)
})
  iqr_rel[n_valid < min_hours] <- NA_real_

  iqr_rel
}


# --------------------------
# IQR in transformed space (your method, wrapped)
# --------------------------
compute_iqr_list <- function(kriging_raw_list, variance_raw_list) {
  z10 <- qnorm(0.10)
  z90 <- qnorm(0.90)

  kriging_mu_sqrt_list <- mapply(
    FUN = function(mu_orig, var_sqrt) {
      mu_sqrt <- sqrt(pmax(mu_orig - var_sqrt, 0))
      mu_sqrt[is.na(mu_orig)] <- NA_real_
      mu_sqrt
    },
    mu_orig  = kriging_raw_list,
    var_sqrt = variance_raw_list,
    SIMPLIFY = FALSE
  )

  q10_sqrt_list <- mapply(
    FUN = function(mu_sqrt, var_sqrt) {
      sigma <- sqrt(var_sqrt)
      q10 <- mu_sqrt + z10 * sigma
      pmax(q10, 0)
    },
    mu_sqrt  = kriging_mu_sqrt_list,
    var_sqrt = variance_raw_list,
    SIMPLIFY = FALSE
  )

  q90_sqrt_list <- mapply(
    FUN = function(mu_sqrt, var_sqrt) {
      sigma <- sqrt(var_sqrt)
      q90 <- mu_sqrt + z90 * sigma
      pmax(q90, 0)
    },
    mu_sqrt  = kriging_mu_sqrt_list,
    var_sqrt = variance_raw_list,
    SIMPLIFY = FALSE
  )

  q10_list <- lapply(q10_sqrt_list, function(x) x^2)
  q90_list <- lapply(q90_sqrt_list, function(x) x^2)

  mapply(function(q90, q10) q90 - q10,
         q90 = q90_list, q10 = q10_list,
         SIMPLIFY = FALSE)
}

compute_year_threshold_stats_simple <- function(rda_file, min_threshold, max_threshold) {

  load(rda_file)  # expects: kriging_crop_list, variance_crop_list, timestamps

  dates  <- as.POSIXct(timestamps, origin = "1970-01-01", tz = "UTC")
  season <- get_season(dates)

  # Strip heavy attrs for compute
  variance_raw_list <- lapply(variance_crop_list, strip_attrs_keep_dim)

  # Keep only values in the band: min_threshold < mu <= max_threshold
  kriging_thr_list <- lapply(kriging_crop_list, function(m) {
    m[!(m > min_threshold & m <= max_threshold)] <- NA_real_
    strip_attrs_keep_dim(m)
  })

  # Annual fields
  annual_mean_mu   <- aggregate_mean(kriging_thr_list)$mean
  annual_accum_mu  <- aggregate_sum(kriging_thr_list)$sum
  annual_wet_hours <- aggregate_wet_hours(kriging_crop_list, min_threshold, max_threshold)

  iqr_list <- compute_iqr_list(kriging_thr_list, variance_raw_list)
  annual_mean_iqr <- aggregate_mean(iqr_list)$mean

  # Seasonal fields
  seasons  <- levels(season)
  seasonal <- setNames(vector("list", length(seasons)), seasons)

  for (s in seasons) {
    idx <- which(season == s)

    s_raw_mu_list <- kriging_crop_list[idx]
    s_mu_list     <- kriging_thr_list[idx]
    s_iqr_list    <- iqr_list[idx]

    s_mean_mu   <- aggregate_mean(s_mu_list)$mean
    s_accum_mu  <- aggregate_sum(s_mu_list)$sum
    s_wet_hours <- aggregate_wet_hours(s_raw_mu_list, min_threshold, max_threshold)
    s_mean_iqr  <- aggregate_mean(s_iqr_list)$mean

    seasonal[[s]] <- list(
      mean_mu   = s_mean_mu,
      accum_mu  = s_accum_mu,
      wet_hours = s_wet_hours,
      mean_iqr  = s_mean_iqr
    )
  }

  list(
    min_threshold = min_threshold,
    max_threshold = max_threshold,
    annual = list(
      mean_mu   = annual_mean_mu,
      mean_iqr  = annual_mean_iqr,
      accum_mu  = annual_accum_mu,
      wet_hours = annual_wet_hours
    ),
    seasonal = seasonal,
    variance_ref_list = variance_crop_list
  )
}

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

      message("  Year ", yr)
      res <- compute_year_threshold_stats_simple(rda_file, threshold = thr)

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
    annual_sd <- list(
      mean_mu   = acc_sd_finalize(acc_annual_sd$mean_mu),
      mean_iqr  = acc_sd_finalize(acc_annual_sd$mean_iqr),
      accum_mu  = acc_sd_finalize(acc_annual_sd$accum_mu),
      wet_hours = acc_sd_finalize(acc_annual_sd$wet_hours)
    )

    seasonal_sd <- lapply(acc_season_sd, function(a) {
      list(
        mean_mu   = acc_sd_finalize(a$mean_mu),
        mean_iqr  = acc_sd_finalize(a$mean_iqr),
        accum_mu  = acc_sd_finalize(a$accum_mu),
        wet_hours = acc_sd_finalize(a$wet_hours)
      )
    })

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