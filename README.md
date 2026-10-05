# Stable cores and marginal churn

**Reliability of a single-year exceedance indicator across India's river monitoring network, 2016-2024**

India's Central Pollution Control Board (CPCB) lists polluted river stretches when a station's reported annual
maximum of biochemical oxygen demand (BOD) exceeds a criterion. This repository holds the data, code and results of a
study that tests how reliable that single-year exceedance status is as an indicator. It uses nine years of reported
annual extremes for 1,891 harmonized stations of the National Water Quality Monitoring Programme (NWMP).

The analysis covers:

- persistence of adverse status against three permutation nulls;
- reliability under analytical error;
- sensitivity to sampling variability and sampling effort, using monthly records for 2020-2021;
- an out-of-sample test of a multi-year, error-aware alternative;
- associations with population density, rainfall and river size (associations, not effects).

## Main findings

| Finding | Value |
|---|---|
| Probability that an adverse station stays adverse next year, P(1\|1) | FC 0.88, BOD 0.81, DO 0.69, pH 0.41 |
| Excess of P(1\|1) over each station's own long-run share (null N2) | 0.05 (FC, BOD, DO), 0.08 (pH) |
| BOD status changes: observed vs expected from analytical error alone | 1,074 vs 902 |
| BOD status changes 2020-2021: observed vs expected from sampling variability | 14 vs 10.9 (95% range 9-13) |
| Three-class indicator (2016-2020 classes): adverse share in 2021-2024 | 81% (reliably adverse) vs 2% (reliably compliant) |
| Annual turnover of the BOD list: single-year vs multi-year class | 0.36 vs 0.18 |
| States where FC compliance cannot be assessed (FC never reported above 2,500 MPN/100 mL) | 11 |

In short, single-year exceedance identifies a stable core of adverse stations reliably but churns at the margin;
listing should rest on multi-year, error-aware status.

## Repository layout

```
.
├── README.md                this file
├── run_all.sh               runs the whole pipeline
├── requirements.txt         Python packages
├── install_r_packages.R     R packages
├── CITATION.cff             how to cite
├── LICENSE                  code licence (MIT)
├── LICENSE-DATA.md          data licence (CC BY 4.0) and third-party terms
├── data/                    inputs (see data/README.md)
│   ├── raw/                 workbook extracted from the nine CPCB reports (never edited)
│   ├── processed/           harmonized station-year table, station dictionary, row-to-station map
│   ├── context/             rainfall, population, basins, river network and upstream population
│   ├── monthly/             links to the Esri India monthly records (records not included)
│   └── geo/                 state boundaries for the map (Natural Earth)
├── analysis/                analysis datasets built by the pipeline (status coding, context table)
├── code/
│   ├── r20/                 analysis pipeline, steps 00-15 (Python and R)
│   ├── local/               river-network step run on a local computer (HydroSHEDS, WorldPop)
│   └── tools/               download helper for the monthly records
├── results/                 every result table (CSV) behind the paper's numbers
├── figures/                 figures written by the pipeline (PDF and PNG)
└── manuscript/              LaTeX source of the manuscript and supplementary material
```

## Reproducing the analysis

Requirements: Python 3.11 or newer and R 4.3 or newer. The pipeline was last run with Python 3.13 (pandas 3.0,
numpy 2.5, scipy 1.18, geopandas 1.2, matplotlib 3.11) and R 4.3.3 (lme4 1.1-35, blme 1.0-5, survival 3.5-8,
sandwich 3.1-0, lmtest 0.9-40).

```sh
git clone <repository-url>
cd stable-cores-marginal-churn

python3 -m venv .venv
. .venv/bin/activate                      # Windows: .venv\Scripts\activate
pip install -r requirements.txt
Rscript install_r_packages.R

# optional: monthly records for steps 04, 05b and part of 07
python code/tools/download_esri_monthly.py

sh run_all.sh                             # Windows: run in Git Bash or WSL
```

`run_all.sh` rewrites `analysis/`, `results/`, `figures/`, `manuscript/generated/` and
`manuscript/supplementary.tex`, and writes a log per step to `logs/`. Random steps use fixed seeds, so a rerun
reproduces the deposited results. Without the monthly records the steps that need them are skipped and their
deposited result files are kept. On a laptop the full run takes roughly half an hour; most of it is the bootstrap
in step 01 and the R context models in step 09.

Each script can also be run on its own from the repository root, for example `python code/r20/03_noise.py`. Set the
environment variable `WQ_ROOT` to run the scripts against a copy of the data elsewhere.

## Pipeline

| Step | Script | What it does | Main outputs |
|---|---|---|---|
| 00 | `00_build.py` | Status coding from the harmonized table (bathing criteria), panels | `analysis/analysis_r20.csv` |
| 02 | `02_fc_censoring_heaping.py` | FC censoring (ceiling states, values of exactly 1,600), heaping at the criteria, microbial vs organic typology | `F1-F3`, `H1-H2` |
| 01 | `01_core_reproduce.py` | Persistence, transitions, null models N1, N1s, N2, overlap, station classes, clustered intervals | `T3_core`, `T3_by_panel` |
| 03 | `03_noise.py` | Analytical-error benchmark: expected vs observed status changes | `E1`, `E2` |
| 04 | `04_monthly.py` | Monthly records: count reconciliation, sampling-variability nulls, COVID-19 sensitivity, subsampling | `M1-M8` |
| 05 | `05_regime.py`, `05b_fixed_effort.py` | Persistence by sampling frequency, the 2020-2021 pair, fixed sampling effort | `R1`, `R2`, `R5` |
| 06 | `06_models.R` | Mixed models: ICC, dynamic model with initial conditions, entry cohort and frequency | `R3`, `R4` |
| 07 | `07_indicator.py` | Three-class, multi-year, error-aware indicator; out-of-sample test; list turnover; compliance rules | `I0-I4` |
| 08 | `08_context_prep.py` | Station-year table for the context models; main spatial sample | `analysis/context_r20.csv`, `analysis/stations_r20.csv` |
| 09 | `09_context.R`, `09b_loso_glm.R` | Population density, rainfall, basins, river size and upstream population (GLMM, GLM, IPW, leave-one-state-out, conditional logit) | `C0-C2` |
| 10 | `10_robustness.py` | All core metrics under 19 analysis variants | `S_robustness_full` |
| 11 | `11_figure_data.py` | Data behind figures not produced elsewhere | `G1-G4` |
| 12-15 | `12_figures_a.py`, `12_figures_b.py`, `13_tables.py`, `14_context_outputs.py`, `15_supplement.py` | Figures, LaTeX tables and the supplementary material, generated from the result files | `figures/`, `manuscript/` |

`core.py` holds the shared functions (station-by-year matrices, transitions, nulls, overlap, classes, cluster
bootstrap) and `figstyle.py` the figure style.

## Definitions in brief

- **Adverse status** (Primary Water Quality Criteria for Bathing Water): DO minimum below 5 mg/L, BOD maximum above
  3 mg/L, pH outside 6.5-8.5, FC maximum above 2,500 MPN/100 mL. Values at a limit are compliant.
- **FC censoring**: FC status is missing in the 11 states whose reported FC never exceeds the criterion and wherever the
  reported maximum is exactly 1,600 MPN/100 mL.
- **Panels**: A (1,891 stations, 11,111 station-years), A_long (stations observed in at least two years; 1,724 stations,
  10,944 station-years), C (all nine years; 467 stations).
- **Nulls**: N1 reshuffles statuses across India by year, N1s within state and year, N2 permutes each station's
  statuses across its own years.
- **Three-class indicator**: each station-year gets its distance from the criterion in analytical-error units; a station
  is reliably adverse when at least 60% of its years are confidently adverse, reliably compliant when at least 60% are
  confidently compliant and none confidently adverse, and borderline otherwise.

## Data sources

| Source | Used for | In this repository |
|---|---|---|
| CPCB annual river water-quality reports, 2016-2024 | Annual minima and maxima per station | Extracted workbook and harmonized tables |
| Esri India Living Atlas, "Time Series River Water Quality 2020_2025" | Monthly records, 2020-2021 | Links only; download with `code/tools/download_esri_monthly.py` |
| IMD gridded daily rainfall, 0.25 degree | Station-year rainfall | Station-level summaries |
| WorldPop constrained population, 100 m (release 2025A) | Population density within 1 and 5 km; upstream population | Station-level summaries |
| HydroSHEDS: HydroRIVERS v1.0, HydroBASINS v1c | Basins, river reach, discharge, upstream catchment | Station-level attributes |
| Natural Earth admin-1, 1:10m | Map in Fig. 2 | State polygons used by the figure |

Map lines in the figures delineate study areas and do not necessarily depict accepted national boundaries.

## Manuscript

`manuscript/` holds the LaTeX source (Elsevier CAS double-column class, numbered references). Compile with
`pdflatex manuscript; bibtex manuscript; pdflatex manuscript; pdflatex manuscript`. The tables in `manuscript.tex`
are copies of the generated fragments in `manuscript/generated/`; see `manuscript/README.md`.

## Citation

Please cite the article once published and the data deposit; `CITATION.cff` holds the details. The Zenodo DOI will be
added here when it is reserved.

## Licence

Code: MIT (see `LICENSE`). Data and results produced in this project: CC BY 4.0 (see `LICENSE-DATA.md`, which also
lists the terms of the third-party sources).

## Author

Ranjit Singh, Leuphana University Lüneburg. Questions and corrections are welcome as GitHub issues.
