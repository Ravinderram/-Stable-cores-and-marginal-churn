#!/bin/sh
# Full analysis pipeline (revision 20). Run from anywhere:  sh run_all.sh
# Python steps take a few minutes; the R context models (step 09) take longest.
# Steps that need the Esri monthly records (04, 05b and the last part of 07) are skipped when
# data/monthly/esri_timeseries_2020_2025.csv is missing; the deposited result files are then kept.
set -e
cd "$(dirname "$0")"
export PYTHONWARNINGS=ignore
P=code/r20
mkdir -p logs results figures analysis manuscript/generated

echo "== core: status coding, persistence, nulls"
python3 $P/00_build.py                  > logs/00_build.log
python3 $P/02_fc_censoring_heaping.py   > logs/02_fc_censoring_heaping.log
python3 $P/01_core_reproduce.py         > logs/01_core_reproduce.log
python3 $P/03_noise.py                  > logs/03_noise.log

if [ -f data/monthly/esri_timeseries_2020_2025.csv ]; then
  echo "== monthly records: sampling variability and fixed effort"
  python3 $P/04_monthly.py              > logs/04_monthly.log
  python3 $P/05b_fixed_effort.py        > logs/05b_fixed_effort.log
else
  echo "== monthly records not found: steps 04 and 05b skipped (see data/README.md)"
fi

echo "== monitoring regime, mixed models, status indicator, robustness"
python3 $P/05_regime.py                 > logs/05_regime.log
Rscript $P/06_models.R                  > logs/06_models.log 2>&1
python3 $P/07_indicator.py              > logs/07_indicator.log
python3 $P/10_robustness.py             > logs/10_robustness.log
python3 $P/11_figure_data.py            > logs/11_figure_data.log

echo "== context models (R; two indicators at a time)"
python3 $P/08_context_prep.py           > logs/08_context_prep.log
printf "ADV_BOD\nADV_DO\nADV_PH\nADV_FC_cens\nADV_MULTI_GE2_cens\n" | \
  xargs -P 2 -I{} sh -c "Rscript $P/09_context.R {} > logs/09_context_{}.log 2>&1"
python3 -c "
import pandas as pd, glob
pd.concat([pd.read_csv(f) for f in sorted(glob.glob('results/C1_part_*.csv'))]).to_csv('results/C1_context_models.csv', index=False)"
Rscript $P/09b_loso_glm.R               > logs/09b_loso_glm.log 2>&1

echo "== figures and tables"
python3 $P/12_figures_a.py              > logs/12_figures_a.log
python3 $P/12_figures_b.py              > logs/12_figures_b.log
python3 $P/14_context_outputs.py        > logs/14_context_outputs.log
python3 $P/13_tables.py                 > logs/13_tables.log
python3 $P/15_supplement.py             > logs/15_supplement.log
echo "done: results/, figures/, manuscript/generated/ and manuscript/supplementary.tex updated"
