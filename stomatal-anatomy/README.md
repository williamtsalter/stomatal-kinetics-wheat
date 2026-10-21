# Stomatal anatomy — linear mixed model analysis

Linear mixed-model (LMM) analysis of the stomatal anatomy traits for the **wheat
heat-and-drought stomatal-kinetics experiment**, using **ASReml-R** and **asremlPlus**. See the
[repository root README](../README.md)
for study context, software notes and the ASReml-R licence requirement.

## Traits analysed

`StomatalCount`, `StomatalDensity` (SD), `StomatalArea` (SA), `StomatalWidth` (guard-cell width,
GCW), `StomatalLength` (guard-cell length, GCL), `PoreArea`, and `gsmax` (anatomically-derived
maximum conductance). Count, density and *g*smax are **pot-level** traits; area, width and length
are measured on **individual stomata** within a pot.

## Directory contents

```
stomatal-anatomy/
├── data/
│   ├── Stomatal Anatomy Data.xlsx        <- raw input (2 sheets: pot-level + stoma-level)
│   ├── indv.design.csv                   <- bundled design factors (see "Design data" below)
│   ├── NewGenotypeLabels.csv             <- genotype labels assigned in manuscript
│   ├── anat.dat.v1.Rdata / v2            <- imported, tidied data (v2 = outliers removed)
│   ├── anat.stoma.dat.* / anat.pot.dat.* <- stoma- and pot-level frames
│   ├── anat.Analyses.v1.RData / v2       <- fitted models (list, one per trait)
│   ├── anat.Preds.v1.Rdata / v2          <- predictions for plotting (v2 = new genotype labels)
│   └── anat.diffs.full.rda               <- full-model predictions exported for correlations
└── src/
    ├── anat.import.v1.r / v2.r    <- 1. import & tidy
    ├── load.indv.design.R                <-    rebuilds indv.dat from indv.design.csv
    ├── anat.anal.v1.r / v2.r      <- 2. fit models (parallel over traits)
    ├── anat.anal.v1.stub.r        <-    shared setup + showIndvModelResults()
    ├── anal.stub.funcs.r                 <-    reusable model/diagnostic functions (do not edit*)
    ├── anat.anout.v1.r / v2.r     <- 3. render results (Wald tables, variances, outliers)
    ├── anat.pred.v1.r / v2.r      <- 4. predictions + publication plots
    ├── utilities.R                       <-    plotting helpers (grid layout etc.)
    └── figures/                          <-    generated figures (regenerable)
```

\* `anal.stub.funcs.r` is a shared library of Chris Brien's model-fitting, variance-model
selection and diagnostic functions ("CB" in the header). Per its own note, its defaults and
functions should not be changed without CB's approval.

## Run order

Scripts use **relative paths**, so set the working directory to `src/` before running. Render any
script to its report with `knitr::spin("<script>.r")`.

1. **`anat.import.v2.r`** — reads `data/Stomatal Anatomy Data.xlsx`, rebuilds the design
   factors from `data/indv.design.csv` (via `source("load.indv.design.R")`), joins them on, and
   writes `data/anat.dat.v2.Rdata`.
2. **`anat.anal.v2.r`** — loads `anat.dat.v2.Rdata`, sources `anat.anal.v1.stub.r`
   and `anal.stub.funcs.r`, fits one LMM per trait **in parallel** (`foreach`/`doParallel`), and
   writes `data/anat.Analyses.v2.RData`.
3. **`anat.anout.v2.r`** — loads the fitted models and renders results: per-trait
   summaries, Wald F-test tables, chosen terms, residual variances, and extremal-residual
   (outlier) diagnostics.
4. **`anat.pred.v2.r`** — loads `anat.Analyses.v2.RData`, produces EMM predictions and the
   publication figure (`figures/AnatomicalTraitPredictions.*`), and exports
   `data/anat.diffs.full.rda` **which the conductance-vs-anatomy correlation script consumes.**

The `v1` scripts are the earlier pass (before outlier removal) and are retained for provenance;
**v2 underpins the manuscript**. `anat.pred.v1.r` deliberately loads the **v2** analyses.

## Model

For each trait,

```
y = mean + core-fixed (Room_Temperature nested-factorial with Watering, Genotype)
        + spatial-fixed + spatial-random + residual
```

Residual variances are allowed to differ across `Room_Temperature × Watering`; four variance
models are compared by **AIC** and the most parsimonious kept. Terms are tested with **Wald
F-tests** (α = 0.05), from the `Watering:Genotype` interaction down to main effects. Residual-vs-
fitted and normal-probability plots are inspected. The chosen model yields **EMMs** with LSD(5%)
and Least Significant Intervals for error bars. Full design detail is in the manuscript's
Supporting Information Methods S1.

## Design data (`indv.design.csv`)

The import needs experimental **design factors** that are not in the anatomy spreadsheet:
`Cell.ID` (used in the outlier-removal block and for ordering), the spatial `Room`/`Lane`/
`Position`/`Block`/`SubUnit`/`Cell` terms, `Line`/genotype, `Replicate`, `Rm_Temperature`,
`Watering`, and `WU.Total`. These originate from the earlier water-use analysis.

Originally the import loaded them from an external object, `../../data/indv.v2.Rdata`, which is
**not distributed** with the project. To make the stream self-contained, that 160-row design
table was extracted and bundled as **`data/indv.design.csv`**. The helper
**`src/load.indv.design.R`** reads the CSV and reconstructs `indv.dat`, setting every column's
type and **factor level order explicitly** so the object is identical to the original — the
fitted models are therefore unchanged. The import scripts call `source("load.indv.design.R")` in
place of the old `load()`; a commented line shows how to revert to the original object if you
ever have it.

> Provenance note: `indv.design.csv` was derived from the `indv.dat` object stored inside the
> supplied `anat.dat.v2.Rdata`; a value-level check confirmed an exact round-trip (0 mismatches
> across all 160 rows × 26 columns).

## Requirements

R with `asreml` (**licensed — see root README**), `asremlPlus`, `dae`, `tidyverse`, `knitr`,
`foreach`, `doParallel`, `parallel`, `emmeans`, `patchwork`, `qqplotr`, `gridExtra`. Record exact
versions with `sessionInfo()` when re-running.
