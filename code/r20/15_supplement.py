"""Rev 20: Supplementary material (supplementary.tex), tables generated from result files."""
import pandas as pd, numpy as np
import os
B = os.environ.get('WQ_ROOT', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))); R = f'{B}/results/'; L = f'{B}/manuscript/'   # supplementary.tex sits next to figs/
def esc(s): return str(s).replace('&', '\\&').replace('%', '\\%').replace('_', ' ').replace('>=', '$\\ge$').replace('<', '$<$').replace('>', '$>$')
def tab(df, caption, label, fmt=None, size='\\scriptsize', colspec=None):
    fmt = fmt or {}
    cols = list(df.columns); n = len(cols)
    first = 4.2 if df[cols[0]].astype(str).str.len().max() > 18 else 2.6
    w = (25.6 - first - n * 0.45) / max(n - 1, 1)
    spec = f'>{{\\raggedright\\arraybackslash}}p{{{first}cm}}' + ''.join(f'>{{\\raggedright\\arraybackslash}}p{{{w:.2f}cm}}' for _ in range(n - 1))
    lines = [f'\\begin{{longtable}}{{{spec}}}', f'\\caption{{{caption}}}\\label{{{label}}}\\\\', '\\toprule',
             ' & '.join(esc(c) for c in cols) + ' \\\\', '\\midrule', '\\endfirsthead', '\\toprule', ' & '.join(esc(c) for c in cols) + ' \\\\', '\\midrule', '\\endhead']
    for _, r in df.iterrows():
        cells = []
        for c in cols:
            v = r[c]
            if isinstance(v, (float, np.floating)):
                cells.append('' if pd.isna(v) else (fmt.get(c, '{:.3f}').format(v)))
            else: cells.append(esc(v))
        lines.append(' & '.join(cells) + ' \\\\')
    lines += ['\\bottomrule', '\\end{longtable}']
    return '{' + size + '\n' + '\n'.join(lines) + '\n}\n'
parts = [r'''\documentclass[a4paper,11pt]{article}
\usepackage[a4paper,landscape,margin=1.6cm]{geometry}\usepackage{array}\usepackage{booktabs,longtable,graphicx}\usepackage[T1]{fontenc}
\renewcommand{\thetable}{S\arabic{table}}\renewcommand{\thefigure}{S\arabic{figure}}
\title{Supplementary material\\[4pt]\large Stable cores and marginal churn: reliability of a single-year exceedance indicator across India's river monitoring network, 2016-2024}
\date{}
\begin{document}\maketitle
\noindent All tables are generated from the deposited result files by \texttt{r20/15\_supplement.py}. FC is censored (missing) in 11 states whose reported FC never exceeds 2,500 MPN/100 mL and at exactly 1,600 MPN/100 mL, unless stated. Panel A\_long unless stated.
''']
P = pd.read_csv(R + 'T3_by_panel.csv')
parts.append(tab(P.rename(columns={'min_years': 'min. years'}), 'Observed \\textit{P}(1|1) minus the station-propensity null N2 by minimum number of observed years (Panels A\\_long, B4, B5, B6, C).', 'tbl:s1'))
S = pd.read_csv(R + 'S_robustness_full.csv')
S['multi_contrast'] = S.P11_with_other - S.P11_alone
parts.append(tab(S[['variant', 'indicator', 'stations', 'prevalence', 'P11', 'P01', 'N2', 'obs_minus_N2', 'jaccard1', 'persistent', 'multi_contrast']],
                 'Full robustness grid: every core metric under every variant. multi\\_contrast: \\textit{P}(1|1) when another primary indicator is adverse minus \\textit{P}(1|1) when the indicator is adverse alone.', 'tbl:s2',
                 fmt={'stations': '{:.0f}'}, colspec='p{4.2cm}lrrrrrrrrr'))
H = pd.read_csv(R + 'H1_heaping_state.csv').rename(columns={'Unnamed: 0': 'state'})
H = H[H.station_years >= 30]
parts.append(tab(H, 'Value heaping at the criteria by state: share of reported values exactly at each value (neighbour columns: mean share at 2.9 and 3.1 mg/L, and at 4.9 and 5.1 mg/L). States with at least 30 station-years; ALL = Panel A.', 'tbl:s3',
                 fmt={'station_years': '{:.0f}'}, colspec='p{3.2cm}rrrrrrrrr'))
HY = pd.read_csv(R + 'H2_heaping_year.csv')
parts.append(tab(HY, 'Value heaping at the criteria by year (Panel A).', 'tbl:s3b', fmt={'Year': '{:.0f}'}))
R1 = pd.read_csv(R + 'R1_persistence_by_frequency.csv'); R3 = pd.read_csv(R + 'R3_models_regime.csv'); R5 = pd.read_csv(R + 'R5_fixed_effort.csv')
parts.append(tab(R1, 'Persistence by sampling frequency reported in the CPCB 2023 table (not reported: mostly unlocated stations).', 'tbl:s4', fmt={'stations': '{:.0f}', 'states': '{:.0f}'}))
R3s = R3[['indicator', 'n_sy', 'stations', 'ICC_M1', 'ICC_M1_freq', 'M2_rows', 'rows_after_gap_dropped', 'lagOR_M2', 'lagOR_M2_lo', 'lagOR_M2_hi', 'lagOR_M2c', 'cohort2021_OR']]
parts.append(tab(R3s, 'Mixed models with monitoring-regime covariates. M1: year + station intercept; M1\\_freq: + 2023 sampling frequency; M2: lag + initial condition; M2c: + entry cohort + frequency. FC models fitted with blme.', 'tbl:s4b',
                 fmt={'n_sy': '{:.0f}', 'stations': '{:.0f}', 'M2_rows': '{:.0f}', 'rows_after_gap_dropped': '{:.0f}'}, colspec='p{2.2cm}rrrrrrrrrrr'))
parts.append(tab(R5, 'Between-station structure under standardized sampling effort (monthly subset, 2020-2021): all samples versus one sample per quarter (50 draws).', 'tbl:s4c', fmt={'stations': '{:.0f}'}))
F1 = pd.read_csv(R + 'F1_ceiling_states.csv'); F2 = pd.read_csv(R + 'F2_fc_censored_core.csv'); F3 = pd.read_csv(R + 'F3_typology.csv')
parts.append(tab(F1, 'Reported FC maxima by state; ceiling states are those whose reported FC maximum never exceeded 2,500 MPN/100 mL.', 'tbl:s5', fmt={'FC_station_years': '{:.0f}', 'FC_max_reported': '{:,.0f}'}))
parts.append(tab(F2[['rule','stations','station_years','states','prevalence','P11','CI_station','CI_state','P01','N2','obs_minus_N2','share_beyond_state','jaccard1','persistent','never']], 'FC persistence under alternative codings of non-assessable values.', 'tbl:s5b', fmt={'stations': '{:.0f}', 'station_years': '{:.0f}', 'states': '{:.0f}', 'classified': '{:.0f}'}, colspec='p{2.6cm}rrrrlllrrrrrrrrrr', size='\\tiny'))
parts.append(tab(F3, 'Microbial versus organic typology (50\\% cut-off) with the ceiling states coded compliant (rev 19) and with FC censored.', 'tbl:s5c', fmt={'stations': '{:.0f}'}, size='\\tiny'))
E1 = pd.read_csv(R + 'E1_mpn_precision.csv'); E2 = pd.read_csv(R + 'E2_error_benchmark.csv')
parts.append(tab(E1, 'Standard deviation of log$_{10}$ MPN: expected-information value from Haldane\'s likelihood for a 10, 1, 0.1 mL series (median over the central density range, computed) and the approximation $0.58/\\sqrt{n}$ for a tenfold series.', 'tbl:s6'))
parts.append(tab(E2, 'Analytical-error benchmark under three error models. expected\\_changes\\_error\\_only: status changes expected if each pair of consecutive years kept its level; share\\_indistinguishable\\_rev19: the measure used in rev 19 (share of observed changes whose two values differ by less than $1.96\\sqrt{2}$ SD), kept for continuity.', 'tbl:s6b',
                 fmt={'pairs': '{:.0f}', 'observed_changes': '{:.0f}', 'lo95': '{:.0f}', 'hi95': '{:.0f}', 'expected_changes_error_only': '{:.1f}'}, size='\\tiny'))
M1 = pd.read_csv(R + 'M1_reconciliation.csv').fillna('')
parts.append(tab(M1, 'Reconciliation of the monthly-record counts.', 'tbl:s7', fmt={'value': '{:.0f}'}, colspec='p{12cm}rl'))
M3 = pd.read_csv(R + 'M3_agreement.csv'); M4 = pd.read_csv(R + 'M4_within_station.csv'); M5 = pd.read_csv(R + 'M5_subsampling.csv'); M6 = pd.read_csv(R + 'M6_null_flips.csv')
parts.append(tab(M3, 'Agreement of monthly-derived annual extremes with the CPCB tables (station-years in both sources).', 'tbl:s7b', fmt={'station_years': '{:.0f}', 'stations': '{:.0f}'}))
parts.append(tab(M4, 'Within-station comparison of adverse status in the year with more versus fewer samples, 2020-2021.', 'tbl:s7c', fmt={'stations_both_years': '{:.0f}', 'stations_unequal_counts': '{:.0f}', 'discordant_more_only': '{:.0f}', 'discordant_fewer_only': '{:.0f}'}))
parts.append(tab(M5, 'Status from $k$ random samples (or one per quarter) against status from all samples, station-years with at least ten samples.', 'tbl:s7d', fmt={'station_years': '{:.0f}'}))
parts.append(tab(M6, 'Sampling-variability null: status changes 2020 to 2021 against changes expected if nothing changed, for the main set and sensitivity variants. pooled: pooled-split null; swap: month-matched swap null; p\\_upper: one-sided Monte Carlo p-value.', 'tbl:s7e',
                 fmt={'stations': '{:.0f}', 'adverse_2020': '{:.0f}', 'observed_changes': '{:.0f}', 'pooled_expected': '{:.1f}', 'swap_expected': '{:.1f}', 'pooled_lo': '{:.0f}', 'pooled_hi': '{:.0f}', 'swap_lo': '{:.0f}', 'swap_hi': '{:.0f}'},
                 colspec='lp{3.4cm}rrrrrrrrrrrr', size='\\tiny'))
I4 = pd.read_csv(R + 'I4_rules_monthly.csv')
parts.append(tab(I4, 'Compliance rules in the monthly subset, 2020-2021.', 'tbl:s8', fmt={'stations': '{:.0f}', 'station_years': '{:.0f}', 'classifications_differing_from_single_sample': '{:.0f}', 'changes_2020_2021': '{:.0f}', 'expected_changes_no_change_null': '{:.1f}', 'null_lo': '{:.0f}', 'null_hi': '{:.0f}'},
                 colspec='lp{3.2cm}rrrrrrrr', size='\\tiny'))
I1 = pd.read_csv(R + 'I1_tci_out_of_sample.csv')
parts.append(tab(I1.drop(columns=['class_2021_2024_missing']), 'Three-class multi-year indicator: classes from 2016-2020 and outcomes in 2021-2024.', 'tbl:s8b', fmt={'stations': '{:.0f}', 'station_years_2021_2024': '{:.0f}', 'class_2021_2024_missing': '{:.0f}'}, size='\\tiny'))
try:
    C = pd.read_csv(R + 'SI_context_all_models.csv')
    C = C[['indicator', 'model', 'term', 'OR', 'lo', 'hi', 'n', 'stations']]
    parts.append(tab(C, 'All context models (main spatial sample): within-state models, population-averaged models, selection weights, leave-one-state-out, state-specific slopes, persistent class, rainfall (exposure per SD = 29\\% wetter year), variance shares (OR column holds the share) and upstream population and discharge. LOSO rows: conditional model refitted with the fast nAGQ = 0 approximation, comparable with each other but not with the Laplace estimate; Table S21 gives leave-one-state-out estimates for the population-averaged GLM used in the main text.', 'tbl:s9',
                     fmt={'n': '{:.0f}', 'stations': '{:.0f}'}, colspec='lp{6cm}p{3cm}rrrrr', size='\\tiny'))
except FileNotFoundError:
    pass
try:
    LG = pd.read_csv(R + 'C2_loso_glm.csv')
    parts.append(tab(LG, 'Leave-one-state-out estimates of the population-averaged within-state GLM (state and year fixed effects, station-clustered errors); none = all states.', 'tbl:s10loso', fmt={'stations': '{:.0f}'}))
except FileNotFoundError:
    pass
parts.append(r'''\begin{figure}[h]\centering\includegraphics[width=\textwidth]{figs/FigS1_station_classes.pdf}
\caption{Station trajectory classes among stations with at least four observed years, with the range of the persistent share across 12 rules (recurrence 50, 60, 67, 75\%; runs of 2, 3, 4 years).}\end{figure}
\begin{figure}[h]\centering\includegraphics[width=\textwidth]{figs/FigS2_FC_ceiling_by_state.pdf}
\caption{Reported FC maxima by state (rev 18 figure, unchanged data).}\end{figure}
\end{document}''')
open(L + 'supplementary.tex', 'w').write('\n'.join(parts)); print('written')
