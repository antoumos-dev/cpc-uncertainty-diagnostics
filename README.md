# CombiPrecip Geostatistical Uncertainty

Quantification of geostatistical uncertainty in the operational CombiPrecip (CPC) precipitation estimation algorithm over Switzerland.
Uncertainty arises from the kriging process used to merge radar and rain-gauge observations, and is analysed across precipitation intensity bins, years, and under convection-control-on vs. convection-control-off experimental conditions.
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
