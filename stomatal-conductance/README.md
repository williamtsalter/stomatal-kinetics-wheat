# Stomatal conductance — dynamic *g*s, *A* and *W*i kinetics

Bayesian non-linear modelling of the *dynamic* responses of stomatal conductance (*g*s),
photosynthetic assimilation (*A*) and intrinsic water-use efficiency (*W*i) to step changes in
irradiance for **the wheat heat-and-drought stomatal-kinetics experiment**. See the
[repository root README](../README.md) for study context, the model equations, and the Stan
toolchain requirement.

There are **two parallel sub-streams**, each with its own `data/` and `src/`:

- **`porometer-li600/`** — high-throughput LI-600 porometer data (all genotypes; *g*s only).
- **`irga-li6800/`** — LI-6800 IRGA validation data (cv. Sunchaser; *g*s, *A* and derived *W*i).

Both fit the McAusland et al. (2016) model with **`brms`/Stan**, separately for the **light-on
(increase, `_i`)** and **light-off (decrease, `_d`)** portions of the protocol:

```
gs(t) = g0 + (g1 - g0) * exp( -exp( ((lambda - t) / k) + 1 ) )
```

Parameters: `g0` (start asymptote; light-on = *g*smin), `g1` (end asymptote; light-on = *g*smax),
`k` (time constant; `ki`/`kd`), `lambda` (lag; `λi`/`λd`). Fitting uses 6 chains × 10,000
iterations (6,000 warmup), `adapt_delta = 0.99`; fixed-effect structures are chosen by **LOOIC**;
effects are summarised as posterior medians with 95 % HPD intervals. The maximum rate of change
is derived as **Slmax = (g1 − g0)/(e·k)**.

> The reference implementation supplied by Dr Silvère Vialet-Chabrand (WUR) is under
> `reference-material/gs_fit/` (`fit_per_curve.r`, `fit_per_species.r`). Those are minimal working
> examples on the published McAusland (2016) dataset and use a stretched-exponential
> parameterisation; the project's own models below use the McAusland (2016) lagged form.

---

## Porometer (LI-600) sub-stream

### Contents & run order

Set the working directory to `porometer-li600/src/`; scripts use relative paths and render with
`knitr::spin()`.

| Step | Script | Loads | Saves |
|---|---|---|---|
| 1. Import | `porometer.import.v1.r` | `data/Porometer Data - 20250701.xlsx` | `data/porometer.v1.Rdata` |
| 2a. Fit light-on | `porometer.anal_i.v1.r` | `porometer.v1.Rdata` | `porometer.anal_i.v1.Rdata` |
| 2b. Fit light-off | `porometer.anal_d.v1.r` | `porometer.v1.Rdata` | `porometer.anal_d.v1.Rdata` |
| 3a. Parameters (on) | `porometer.stomatalCondParamEstimation_i.v1.r` | `porometer.anal_i.v1.Rdata` | `stomatalCondParamEstimation_porometer_i.v1.Rdata` |
| 3b. Parameters (off) | `porometer.stomatalCondParamEstimation_d.v1.r` | `porometer.anal_d.v1.Rdata` | `stomatalCondParamEstimation_porometer_d.v1.Rdata` |
| 4. Plot parameters | `porometer.stomatalCondParamPlotting.v2.r` | both param files | `stomatalCondParamPlotting_porometer.v2.Rdata` |
| 5. *g*s vs. anatomy | `porometer.stomatalCondvsAnatomicalTraits.v2.r` | both param files + `../../../stomatal-anatomy/data/anat.diffs.full.rda` | `stomatalCondvsAnatTraits.v1.Rdata` |

Step 5 requires the anatomy export `anat.diffs.full.rda` — **run the stomatal-anatomy stream
first** (its step 4). The parameter-estimation scripts extract EMMs and posterior draws per
parameter (via `emmeans`/`tidybayes`) and compute `Slmax`.

---

## LI-6800 (IRGA) sub-stream

### Contents & run order

Set the working directory to `irga-li6800/src/`.

| Step | Script | Loads | Saves |
|---|---|---|---|
| 1. Import | `li6800.import.v1.r` | 12 × `L01_*.xlsx` + `PlantPositionsLookupTable.csv` | `data/li6800.v1.Rdata` |
| 2a. *g*s light-on | `li6800.stomatalCondModel_i.v2.r` | `li6800.v1.Rdata` | `stomatalCondModel_i.v2.Rdata` |
| 2b. *g*s light-off | `li6800.stomatalCondModel_d.v2.r` | `li6800.v1.Rdata` | `stomatalCondModel_d.v2.RData` |
| 2c. *A* light-on | `li6800.stomatalCondModel_A_i.v1.r` | `li6800.v1.Rdata` | `stomatalCondModel_A_i.v1.Rdata` |
| 3a. Parameters (*g*s on) | `li6800.stomatalCondParamEstimation_i.v1.r` | `stomatalCondModel_i.v2.Rdata` | `stomatalCondParamEstimation_i.v1.Rdata` |
| 3b. Parameters (*g*s off) | `li6800.stomatalCondParamEstimation_d.v1.r` | `stomatalCondModel_d.v2.Rdata` | `stomatalCondParamEstimation_d.v1.Rdata` |
| 3c. Parameters (*A*) | `li6800.stomatalCondParamEstimation_A_i.v1.r` | `stomatalCondModel_A_i.v1.Rdata` | `stomatalCondParamEstimation_A_i.v1.Rdata` |
| 4. Intrinsic WUE | `li6800.iWUE_i.v1.r` | `stomatalCondParamEstimation_A_i.v1.Rdata` + `..._i.v1.Rdata` | `iWUE.v1.Rdata` |
| 5. Plot parameters | `li6800.stomatalCondParamPlotting.v1.r` | the three param files | `stomatalCondParamPlotting.v1.Rdata` |

**Model notes.** The IRGA *g*s models (`_i`, `_d`) use the lagged McAusland form; the *A* model
(`_A_i`) uses the same form **without the lag** (induction shows little lag). Fixed effects are in
`Room_Temperature × Watering` (no genotype term — a single validation cultivar). Derived from the
chosen models: `Slmax`, `T95%A`, `A95`, `gslimitA`, and `Wi = A/gs` (giving Wimax and Wi95).

> **Filename-case caveat (v2 models).** `stomatalCondModel_d.v2.r` saves
> `stomatalCondModel_d.v2.**RData**` while the downstream parameter script loads
> `...v2.**Rdata**`. On case-sensitive filesystems (Linux) this mismatch will fail to load —
> normalise the extension when re-running.

---

## Data provenance

- **LI-600:** `Porometer Data - 20250701.xlsx` (sheet *All Data Compiled*); key variable `gsw`
  (referred to as *g*s), light state in the `Light` (On/Off) column; a `time.elapsed`/`time.adj`
  variable is derived per operator with time reset to 0 at each light step.
- **LI-6800:** twelve workbooks `L01_<Temp>_<Water>_Re<n>_<date>.xlsx` (logged every 30 s) plus
  `PlantPositionsLookupTable.csv`.
- `.Rdata` files are **intermediate/fitted objects** and are regenerable; fitted brms objects are
  large.

## Requirements

R with `brms` (**+ a working Stan toolchain**), `dae`, `tidyverse` (`dplyr`, `tidyr`, `ggplot2`),
`emmeans`, `tidybayes`, `pracma`, `patchwork`, `GGally`, `lattice`, `knitr`, `readxl`. `asreml`
appears only in commented exploratory code and is not required to reproduce the chosen models.
Record versions with `sessionInfo()` when re-running. Expect multi-hour run times for the full
6-chain × 10,000-iteration fits.
