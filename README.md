# Stomatal responsiveness and photosynthetic induction in wheat under fluctuating light

Analysis code accompanying the manuscript:

> **Heat accelerates stomatal responsiveness and photosynthetic induction in wheat under
> fluctuating light.**
> Chaplin E., Forknall C., Brien C., Merchant A., Trethowan R., Salter W.
> *(submitted to* New Phytologist*).*

This repository contains the R code used to import, model, summarise and visualise the physiological and anatomical data from the **wheat heat-and-drought stomatal-kinetics experiment** outlined in the above study. The analyses quantify how heat and water limitation reshape the *dynamic* responses of stomatal conductance (*g*s) and photosynthesis (*A*) to step changes in irradiance in a diverse panel of heat-tolerant wheat, and test whether stomatal anatomy explains the observed kinetic variation.

---

## Contents

1. [Study overview](#1-study-overview)
2. [What the analysis does](#2-what-the-analysis-does) — detailed description
3. [Repository layout](#3-repository-layout)
4. [Software requirements](#4-software-requirements)
5. [Reproduction](#5-reproduction)
6. [Data provenance and availability](#6-data-provenance-and-availability)
7. [Attribution, contact and citation](#7-attribution-contact-and-citation)

---

## 1. Study overview

Ten heat-tolerant wheat genotypes were grown in a controlled-environment experiment at The Plant
Accelerator (Australian Plant Phenomics Network, Adelaide University) under a factorial
combination of:

- **Temperature** — ambient (22/15 °C) vs. elevated (28/21 °C), imposed in separate
  DroughtSpotter growth chambers; and
- **Water** — well-watered (20 % gravimetric water content) vs. water-limited (12 % GWC).

Because there is a single room per temperature regime, Temperature is confounded with Room and is
treated throughout as **`Room_Temperature`**. For the LI-600 porometer data, Watering is
confounded with Operator/Device and is treated as **`Device_Watering`** (D1 = Device 1,
D2 = Device 2). These confoundings are structural to the design and must be respected when
interpreting any treatment comparison.

Three data streams were collected:

| Stream | Instrument | What it measures | Modelling approach |
|---|---|---|---|
| Dynamic *g*s (high-throughput) | LI-600 porometer | *g*s response to a 0→1000→0 µmol m⁻² s⁻¹ PPFD step protocol | Bayesian non-linear (brms/Stan) |
| Dynamic *g*s, *A*, *W*i (validation) | LI-6800 IRGA | *g*s and *A* over the same light protocol (cv. Sunchaser) | Bayesian non-linear (brms/Stan) |
| Stomatal anatomy | Dino-Lite microscopy + deep learning | SD, SA, GCW, GCL, and derived *g*smax | Linear mixed models (ASReml-R) |

---

## 2. What the analysis does

This section describes, in detail, what the code in this repository actually does — the
processing steps, the statistical models, the quantities estimated, and the reproducibility
engineering applied when preparing the repository. The work was carried out by APPN biometry staff
(Clayton Forknall and Chris Brien) using a reproducible, script-per-step workflow.

### 2.1 Design and measurement context

The experiment ran in two phases: a glasshouse establishment phase and a controlled-environment
(DroughtSpotter) phase in which the temperature and watering treatments were imposed. The
experimental designs for both phases were generated and randomised with the `dae` package, and
water use was derived from the DroughtSpotter load-cell logs with `growthPheno`. Gas-exchange
measurements followed a three-stage light protocol (1 h dark at 0 µmol m⁻² s⁻¹ → 1 h at 1000 µmol m⁻² s⁻¹ → 1 h
dark at 0 µmol m⁻² s⁻¹), captured at high throughput on the LI-600 porometer and, for validation on the check
cultivar Sunchaser, in detail on the LI-6800 IRGA. Stomatal anatomy was measured in situ from
microscopy images analysed with a deep-learning pipeline (Chaplin et al., 2025; separate
repository).

### 2.2 Stomatal anatomy — linear mixed models (`stomatal-anatomy/`)

Raw microscopy-derived traits are read from the anatomy spreadsheet and joined to the
experimental design. Five traits are analysed: stomatal density (SD), stomatal
area (SA), guard-cell width (GCW), guard-cell length (GCL), and the anatomically-
derived maximum conductance *g*smax. Density and *g*smax are pot-level traits; area, width
and length are measured on individual stomata within each pot.

Each trait is analysed with a linear mixed model of the form

```
y = mean + core-fixed + spatial-fixed + spatial-random + residual
```

fitted with **ASReml-R** and **asremlPlus**. The **core fixed** model is a nested-factorial
treatment structure in `Room_Temperature`, `Watering` and `Genotype` (reflecting the
Temperature/Room confounding); **spatial** fixed and random terms come from the design. Residual
variances are allowed to differ across `Room_Temperature × Watering` combinations, and four
variance models (full heterogeneity, temperature-only, watering-only, homogeneous) are compared
by **AIC** to choose the most parsimonious. Model assumptions are checked with residual-vs-fitted
and normal-probability plots and extremal-residual (outlier) diagnostics; a small number of
individual-stomata outliers identified this way are set to missing in the `v2` re-analysis.
Treatment effects are tested with **Wald F-tests** (α = 0.05), from the highest-order interaction
downward, and a chosen model of the significant terms yields **estimated marginal means (EMMs)**
with LSD(5 %) and Least Significant Intervals for plotting. Full-model predictions for every
trait are exported for the conductance-vs-anatomy correlation analysis.

### 2.3 Dynamic *g*s — Bayesian non-linear model (`stomatal-conductance/`)

The temporal response of *g*s to a single step change in PPFD is modelled with the sigmoidal
function of McAusland et al. (2016), fitted **separately for the light-on (increase, `_i`) and
light-off (decrease, `_d`) portions**, with time reset to 0 at each step:

```
gs(t) = g0 + (g1 - g0) * exp( -exp( ((lambda - t) / k) + 1 ) )
```

with parameters

- **g0** — starting asymptote (light-on: minimum *g*s, reported as *g*smin);
- **g1** — final asymptote (light-on: maximum *g*s, reported as *g*smax);
- **k** — time constant (opening `ki`, closing `kd`);
- **lambda (λ)** — initial lag in the response (`λi`, `λd`).

Models are fitted with **`brms`** (Stan back-end). The maximal model estimates population-level
(fixed) treatment effects for every parameter — genotype (`Line`) × water (`Device_Watering` for
LI-600 / `Watering` for LI-6800) × `Room_Temperature` — with subject-level (random) effects on
`g0` and `g1` for the design and spatial terms. Priors are informed by McAusland et al. (2016)
and Woning et al. (2026). Sampling uses **6 chains × 10,000 iterations (6,000 warmup)** with
`adapt_delta = 0.99`, and convergence/behaviour is checked with the usual Stan diagnostics and
posterior-predictive checks. Competing fixed-effects structures are compared by **leave-one-out
cross-validation (LOOIC)**, working from the highest-order interaction down, to arrive at a
parsimonious chosen model for each parameter. Effects are summarised as **posterior medians with
95 % Highest Posterior Density (HPD) intervals**, and the maximum rate of change is derived as
**Slmax = (g1 − g0) / (e · k)** (after Woning et al., 2026).

Two instruments are analysed in parallel: the **LI-600 porometer** (all genotypes; high
throughput; `porometer-li600/`) and the **LI-6800 IRGA** (validation on cv. Sunchaser;
`irga-li6800/`).

### 2.4 Photosynthetic induction (*A*) and intrinsic WUE (*W*i)

For the LI-6800 data, *A* during the light-on portion is modelled with the same functional form
**without a lag term** (induction shows little lag relative to *g*s):

```
A(t) = A0 + (Amax - A0) * exp( -exp( (-t / k) + 1 ) )
```

From the chosen *A* model, **T95%A** (time to reach 0.95·Amax) is obtained by solving for *t*,
**A95** by back-substitution, and **gslimitA** (the limitation of *A* imposed by *g*s) as in
McAusland et al. (2016). **Intrinsic water-use efficiency** is derived per time point as
**Wi = A / g**s from the posterior predictions, giving Wimax and Wi95.

### 2.5 Linking kinetics to anatomy

Finally, the full-model anatomy EMMs are plotted against the estimated *g*s parameters and the
strength of each relationship is quantified with **Pearson correlations**, computed within
environment (`Room_Temperature × Device_Watering`). This is the only place the two streams meet:
the porometer correlation script reads the anatomy export `anat.diffs.full.rda`.

### 2.6 Reproducibility engineering applied in this repository

Beyond the original analysis, the following was done to make the repository self-contained and
publication-ready:

- **Removed an external dependency.** The anatomy import originally loaded design factors
  (`Cell.ID`, spatial terms, genotype, replicate, `WU.Total`, …) from an external object,
  `indv.v2.Rdata`, produced by the water-use analysis and not distributed with the project. A copy
  of that `indv.dat` object was recovered from inside the supplied `anat.dat.v2.Rdata`, exported
  to `stomatal-anatomy/data/indv.design.csv`, and the import rewritten to rebuild it via
  `stomatal-anatomy/src/load.indv.design.R`. Column types and **factor level orders** are set
  explicitly, so the reconstructed object — and therefore every downstream model — is identical
  to the original (verified by a value-level round-trip: 0 mismatches across 160 × 26 cells).
- **Cleaned the directory structure** into two clearly named analysis streams (see below) and
  fixed the single cross-stream path reference accordingly.
- **Separated third-party material.** Reference code from Dr Silvère Vialet-Chabrand is retained
  under `reference-material/gs_fit/` with attribution; copyrighted manuscript PDFs, the published
  Excel macro, duplicate raw-data copies and correspondence are **not** redistributed and are
  instead cited in `reference-material/REFERENCES.md`.

---

## 3. Repository layout

Two top-level analysis streams, each with its own `data/` and `src/` directories and its own
README:

```
.
├── README.md                         <- this file
├── .gitignore
├── figures/                          <- protocol / parameter schematic (manuscript Figure 1)
│
├── stomatal-anatomy/                 <- stomatal anatomy LMM analysis (ASReml-R)
│   ├── README.md
│   ├── data/                         <- anatomy .xlsx, bundled indv.design.csv, intermediate .Rdata
│   └── src/                          <- import / anal / anout / pred scripts (v1, v2)
│
└── stomatal-conductance/             <- dynamic gs / A / Wi analysis (brms)
    ├── README.md
    ├── porometer-li600/              <- LI-600 porometer (all genotypes; gs)
    │   ├── data/
    │   └── src/
    ├── irga-li6800/                  <- LI-6800 IRGA (cv. Sunchaser; gs, A, Wi)
    │   ├── data/
    │   └── src/
    └── reference-material/           <- Vialet-Chabrand gs_fit reference code + REFERENCES.md
```

### Script naming and reproducible output

Every script follows the convention **`<domain>.<step>.v<N>.r`** and is written in the
`knitr::spin` "roxygen comment" style (`#'` for prose, `#+` for chunk labels). Running a script
through `knitr::spin()` produces a matched set of outputs with the same stem — `.md`/`.html` (the
rendered report), `.pdf` (print version) and `.log` (console log) — with plots under `figures/`.
**Only the `.r` files are source;** the rendered files and `figures/` are regenerable.

Where two versions exist (`v1`, `v2`), **v2 is the later revision** and underpins the manuscript
(for anatomy, v2 removes identified outliers) unless the stream README notes otherwise.

Files, scripts and all material on this repository are named according to the response in question: those referring to stomatal opening response are denoted with **_i**, those for the closure response with **_d**, for the _A_ response with **_A** and for the W<sub>i</sub> response with **.iWUE**

---

## 4. Software requirements

Analyses were run under **R** (R Core Team, 2025). Install only what you need for the stream you
intend to run.

**Common / general** — `knitr`, `tidyverse` (incl. `dplyr`, `tidyr`, `ggplot2`, `stringr`,
`readxl`), `patchwork`, `gridExtra`, `qqplotr`, **`dae`** (Brien) for design utilities, and
**`growthPheno`** (Brien) for the water-use traits.

**Anatomy stream** — **`asreml`** (ASReml-R; Butler et al., 2023), **`asremlPlus`** (Brien),
`foreach`, `doParallel`, `parallel`, `emmeans`.

**Conductance stream** — **`brms`** (Bürkner, 2017) + a working **Stan** toolchain (RStan or
`cmdstanr` + a C++ compiler), `emmeans`, `tidybayes`, `pracma`, `GGally`, `lattice`.

> ⚠️ **ASReml-R is commercial software.** It is sold and supported by VSN International and
> requires a licence code validated online at load time (academic licences start from ~£375).
> The Bayesian conductance analyses use the open-source `brms`/Stan stack and have no such
> restriction.

> ⚠️ **Stan toolchain & run time.** `brms` compiles models via Stan, so a C++ toolchain must be
> installed and configured. Fitted brms objects are large and the full models (6 chains × 10,000
> iterations) can take hours; expect substantial memory use.

Record exact package versions with `sessionInfo()` (or use `renv`) when re-running, as several
packages (`asremlPlus`, `dae`, `growthPheno`, `brms`) are actively developed.

---

## 5. Reproduction

See each stream README for the exact run order. In brief:

1. Install R and the packages for the stream(s) you need (Section 4), including an ASReml-R
   licence and/or a Stan toolchain as applicable.
2. Set the working directory to the relevant `src/` folder (scripts use **relative paths** such
   as `../data/...`).
3. Run scripts **in the documented order** — each `import` produces the `.Rdata` consumed by the
   next step. Render with `knitr::spin("<script>.r")` to reproduce the reports and figures.

Both streams are self-contained; no external input files are required (the anatomy design data is
bundled as `indv.design.csv`). The only cross-stream link is that the porometer
conductance-vs-anatomy script consumes the anatomy export `anat.diffs.full.rda`, so run the
anatomy stream first if you need that step.

---

## 6. Data provenance and availability

Raw inputs supplied with the analysis (under each stream's `data/`):

- **LI-600 porometer** — `Porometer Data - 20250701.xlsx` (sheet *All Data Compiled*).
- **LI-6800 IRGA** — twelve workbooks `L01_<Temp>_<Water>_Re<n>_<date>.xlsx` plus
  `PlantPositionsLookupTable.csv`.
- **Stomatal anatomy** — `Stomatal Anatomy Data.xlsx` (pot-level and stoma-level sheets) and the
  bundled design table `indv.design.csv`.

The `.Rdata`/`.rda` files under each `data/` folder are **intermediate objects** (imported data,
fitted models, parameter exports, predictions). They are provided for convenience but are all
regenerable from the raw inputs and scripts; fitted brms objects are large.

Anatomical traits were extracted from microscopy images using the deep-learning pipeline of
Chaplin et al. (2025), maintained separately at
<https://github.com/williamtsalter/FieldDinoMicroscopy>. The manuscript also refers to a Zenodo
deposit for the raw data and LI-6800 instrument settings — **add the DOI here once minted**.

Third-party manuscripts, the published Excel macro, correspondence and duplicate raw-data copies
that were originally supplied are kept out of the published repository (see
`stomatal-conductance/reference-material/REFERENCES.md`); locally they are retained in a
git-ignored `_local/` folder.

---

## 7. Attribution, contact and citation

**Statistical analysis and code:** Clayton Forknall and Chris Brien (The Plant Accelerator,
Australian Plant Phenomics Network / Adelaide University).
**Data collection and experiment:** Edward Chaplin (The University of Sydney).
**Reference *g*s curve-fitting code:** adapted from material provided by Dr Silvère
Vialet-Chabrand (Wageningen University & Research), implementing the McAusland et al. (2016)
model.

**Please credit:** *Clayton Forknall, Edward Chaplin and William Salter.*

**Contact / repository maintainer:** William Salter — william.salter@sydney.edu.au
(corresponding author for the manuscript: Edward Chaplin — edward.chaplin@sydney.edu.au).

**How to cite this repository:** cite the manuscript above (update with volume/DOI on
acceptance). Key methods/software to cite alongside it:

- McAusland L, Vialet-Chabrand S, Davey P, Baker NR, Brendel O, Lawson T. 2016. Effects of
  kinetics of light-induced stomatal responses on photosynthesis and water-use efficiency.
  *New Phytologist* 211, 1209–1220.
- Bürkner P-C. 2017. brms: An R package for Bayesian multilevel models using Stan. *Journal of
  Statistical Software* 80(1), 1–28.
- Butler DG, Cullis BR, Gilmour AR, Gogel BJ, Thompson R. 2023. *ASReml-R Reference Manual
  Version 4.2.* VSN International.
- Brien C. `asremlPlus`, `dae`, `growthPheno` R packages.
- R Core Team. 2025. *R: A language and environment for statistical computing.*
- Woning N, Al-Salman Y, Kaiser E, et al. 2026. *New Phytologist* 249, 2338–2354; Franks PJ,
  Beerling DJ. 2009. *PNAS* 106, 10343–10347 — derived-trait definitions.

**Licence:** _to be selected._ A common choice is a permissive code licence (e.g. MIT) for the
scripts together with a data licence (e.g. CC-BY 4.0) for any data deposited here; add a
`LICENSE` file before making the repository public.
