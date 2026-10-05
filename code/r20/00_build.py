"""Rev 20, step 0: rebuild the analysis dataset from the enriched station-year table (local kit output, out.zip).
The enriched table carries the harmonized station IDs, exclusions, parsed values with status flags, coordinates and
context variables produced in phases 1-8, G and H; original values are untouched. This script only recodes.
Adverse status (bathing criteria, MoEF 2000): DO min < 5 (BDL adverse); BOD max > 3; FC max > 2,500; pH min < 6.5 or
max > 8.5. Values at a limit are compliant; BDL for BOD/FC is compliant; qc_excluded and missing are missing.
Output: analysis/analysis_r20.csv"""
import pandas as pd, numpy as np
import os
B = os.environ.get('WQ_ROOT', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
d = pd.read_csv(f'{B}/data/processed/integrated_station_year_enriched.csv', low_memory=False)
d['in_panel_A'] = d.exclusion_reason.isna() | (d.exclusion_reason == 'single-year station')
d['in_panel_A_long'] = d.use_in_longitudinal == True
fix = {'ANDHRA': 'ANDHRA PRADESH', 'RD_BANGALORE': 'KARNATAKA', '-': np.nan}
d['state_clean'] = d.State.replace(fix)

def adv(val, st, fn, bdl):
    out = pd.Series(np.nan, index=val.index)
    num = st == 'numeric'
    out[num] = fn(val[num]).astype(float)
    out[st == 'BDL'] = bdl
    return out
d['ADV_DO'] = adv(d.DO_min, d.DO_min_status, lambda x: x < 5, 1.0)
d['ADV_BOD'] = adv(d.BOD_max, d.BOD_max_status, lambda x: x > 3, 0.0)
d['ADV_FC'] = adv(d.FC_max, d.FC_max_status, lambda x: x > 2500, 0.0)
lo = adv(d.PH_min, d.PH_min_status, lambda x: x < 6.5, 0.0); hi = adv(d.PH_max, d.PH_max_status, lambda x: x > 8.5, 0.0)
d['ADV_PH'] = np.where((lo == 1) | (hi == 1), 1.0, np.where(lo.notna() & hi.notna(), 0.0, np.nan))
P = ['ADV_FC', 'ADV_BOD', 'ADV_DO', 'ADV_PH']
allobs = d[P].notna().all(1)
d['ADV_MULTI_GE2'] = np.where(allobs, (d[P].sum(1) >= 2).astype(float), np.nan)
# station-level state: modal state over the station's rows (INTER-STATE kept as its own label)
st = d[d.in_panel_A].groupby('station_uid').state_clean.agg(lambda s: s.mode().iat[0] if s.notna().any() else 'UNKNOWN')
d['station_state'] = d.station_uid.map(st).fillna('UNKNOWN')
d.to_csv(f'{B}/analysis/analysis_r20.csv', index=False)
A = d[d.in_panel_A]
print('Panel A', len(A), A.station_uid.nunique(), '| A_long', d.in_panel_A_long.sum(), d[d.in_panel_A_long].station_uid.nunique())
print('prevalence Panel A:', {p: round(A[p].mean(), 4) for p in P + ['ADV_MULTI_GE2']})
print('n obs Panel A:', {p: int(A[p].notna().sum()) for p in P})
