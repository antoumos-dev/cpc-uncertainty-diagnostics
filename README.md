# CombiPrecip Geostatistical Uncertainty

Analysis of the spatiotemporal structure and algorithmic sensitivity of kriging-based expected value and uncertainty in the operational CombiPrecip (CPC) radar-raingauge QPE over Switzerland, based on a decade-long record (2016–2025).

## Project structure

```
out_stats/
├── R/
│   ├── utils.r                        # shared helper functions
│   └── plot_utils.r                   # shared plotting functions
├── kriging/
│   ├── run_krig_year.r                # kriging stats and plots for a single year
│   ├── run_krig_year.sh               # SLURM wrapper for run_krig_year.r
│   ├── run_krig_multiyear.r           # multi-year aggregation
│   └── run_krig_multiyear.sh          # SLURM wrapper for multi-year run
├── cross_val/
│   ├── build_cross_val.r              # build cross-validation datasets
│   └── cross_val_data.r               # cross-validation analysis
├── conv_control/
│   ├── conv_control.r                 # conv-on vs conv-off comparison
│   ├── conv_control.sh                # SLURM wrapper for conv_control.r
│   ├── conv_control_stats.r           # convection control statistics
│   └── conv_control_stats_multiyear.r # multi-year convection control stats
├── diagnostics/
│   ├── intensity_bins.r               # intensity-bin frequency analysis
│   └── spatial_bias_map.r             # spatial bias mapping
├── out_plots/                         # generated figures (per year and interannual)
├── logs/                              # SLURM job logs
└── data/                              # intermediate result files
```

## Data inputs

| Path | Content |
|------|---------|
| `/store_new/mch/msclim/antoumos/R/develop/CPC/data_new_project/` | Conv-on `.rda` files (`CPC<YY>*.rda`) |
| `.../data_new_project/conv_control_off/` | Conv-off `.rda` files |
| `precip_transformed_results_new_<year>.rda` | Kriging input (conv-on) |
| `precip_transformed_results_conv_off_new_<year>.rda` | Kriging input (conv-off) |

## Scripts

### `run_krig_year.r`
Computes kriging-based precipitation statistics for a single year and generates plots.

```bash
Rscript run_krig_year.r <year> [mode] [mu_min]
```

| Argument | Default | Description |
|----------|---------|-------------|
| `year` | required | Four-digit year (e.g. `2023`) |
| `mode` | `all` | `all` = full run; `relunc` = relative uncertainty (IQR/μ) only |
| `mu_min` | `0.05` | Minimum μ threshold to avoid exploding ratios |

Precipitation thresholds: `0.1, 0.5, 1, 2` mm. Swiss domain crop: x ∈ [480, 840], y ∈ [60, 300].
Output written to `out_plots/year_<year>/`.

### `run_krig_year.sh`
SLURM job array (3 tasks) running `run_krig_year.r` for years 2016–2018 in `relunc` mode with `mu_min=0.1`.

```bash
sbatch run_krig_year.sh
```

### `run_krig_multiyear.r` / `run_krig_multiyear.sh`
Aggregates results across multiple years. Submit with:

```bash
sbatch run_krig_multiyear.sh
```

### `conv_control.r` / `cross_val_data.r`
Compare conv-on vs conv-off experiments and perform cross-validation. Set `YEAR` at the top of each script before running interactively.

### `Intesity_bins.r`
Computes frequency distributions across precipitation intensity bins.

## Output plots

Figures are written to `out_plots/`:
- `year_<year>/` — per-year kriging maps and uncertainty plots
- `interannual_*/` — multi-year summary plots

## Environment

Scripts require the MCH CATs R environment. Load via:

```bash
source /users/antoumos/.local/bin/activate-uenv
uenv start --view=climana climana/24.10:rc1
module load r gdal geos hdf5 cats proj sqlite udunits
```

R library paths: `/store_new/mch/msclim/share/CATs/cats/lib-R4.4.0/` and `/store_new/mch/msclim/sideris/R/lib/`.
