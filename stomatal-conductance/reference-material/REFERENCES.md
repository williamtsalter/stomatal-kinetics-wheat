# Reference material

This folder holds third-party reference code that the stomatal-conductance analysis was based
on, plus citations for the source material used during method development. Copyrighted
manuscripts and duplicate data files that were originally supplied with the project are **not**
redistributed here; they are cited below instead.

## Included: `gs_fit/` (reference implementation)

`gs_fit/` contains the original minimal working examples provided by **Dr Silvère
Vialet-Chabrand** (Wageningen University & Research) for fitting the McAusland et al. (2016)
stomatal-conductance response model with `brms`/Stan:

- `fit_per_curve.r` — fits each response curve independently.
- `fit_per_species.r` — fits `tau`/`theta` per group (fixed) and `g0`/`g1` per curve (random).
- `McAusland_compiled.csv` — the example dataset (from McAusland et al., 2016) used to
  illustrate the expected input layout.

These scripts are provided **as received**, for provenance and attribution. They are *not* part
of the reproducible pipeline (the project's own models are under
`../porometer-li600/src/` and `../irga-li6800/src/`) and contain a hard-coded
`setwd("C:\\Users\\viale002\\...")` that must be edited before they will run locally. Please
credit Dr Vialet-Chabrand if you reuse this code.

## Cited (not redistributed here)

Key sources that informed the measurement protocol, the response model and the derived traits:

- **McAusland L, Vialet-Chabrand S, Davey P, Baker NR, Brendel O, Lawson T. 2016.** Effects of
  kinetics of light-induced stomatal responses on photosynthesis and water-use efficiency.
  *New Phytologist* 211, 1209–1220. — *source of the gs/A response model and derived traits.*
- **Vialet-Chabrand S, Matthews JSA, Simkin AJ, Raines CA, Lawson T. 2017.** Importance of
  fluctuations in light on plant photosynthetic acclimation. *Plant Physiology* 173, 2163–2179.
  — *the review accompanying the exponential/sigmoidal Excel macro (`PP2017-00125R1`).*
- **Vialet-Chabrand S, Dreyer E, Brendel O. 2013.** Performance of a new dynamic model for
  predicting diurnal time courses of stomatal conductance at the leaf level. *Plant, Cell &
  Environment* 36, 1529–1546.
- **Woning N, Al-Salman Y, Kaiser E, et al. 2026.** Revisiting the relationship between stomatal
  size and speed across species – a meta-analysis. *New Phytologist* 249, 2338–2354. —
  *informed priors and the Slmax derived trait.*
- **Faralli M, Cockram J, Ober E, Wall S, Galle A, Van Rie J, Raines C, Lawson T. 2019.**
  Genotypic, developmental and environmental effects on the rapidity of gs in wheat. *Frontiers
  in Plant Science.*
- **Acevedo-Siaca LG, McAusland L. 2025.** A guide to understanding and measuring photosynthetic
  induction: considerations and recommendations. *New Phytologist.*
- **Franks PJ, Beerling DJ. 2009.** Maximum leaf conductance driven by CO₂ effects on stomatal
  size and density over geologic time. *PNAS* 106, 10343–10347. — *anatomical g*smax.
- **Franks PJ, Farquhar GD. 2007.** The mechanical diversity of stomata and its significance in
  gas-exchange control. *Plant Physiology* 143, 78–87.

Software citations (`brms`, `ASReml-R`, `asremlPlus`, `dae`, `growthPheno`, R) are listed in the
[repository root README](../../README.md).
