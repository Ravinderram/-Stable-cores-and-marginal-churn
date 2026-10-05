"""Rev 20, step 4 (Major Comment 4): monthly-record benchmark rebuilt so that every count reconciles.
Changes against rev 19 (script 37):
 (a) the benchmark sample is restricted to station-years present in BOTH the CPCB annual table (Panel A_long) and the
     linked monthly records; rev 19 used every linked station with monthly records in both years, including years in
     which the station was absent from the CPCB table, which is why 143 stations exceeded the 106 possible;
 (b) samples are counted as distinct records after removing exact duplicates; sampled months (max 12) are reported
     separately; 149 station-months carry two records with different values (no day is given, so they are kept as two
     samples in the main analysis and collapsed to their mean in a sensitivity run);
 (c) the 2020 COVID-19 lockdown: (i) months June-December only in both years, (ii) a month-matched swap null that
     exchanges a station's 2020 and 2021 sample of the same calendar month, which keeps seasonality.
Links are the one-to-one links of rev 19 (data/monthly/esri_monthly_links.csv). Values never made a link.
Output: results/M1_reconciliation.csv, M2_counts.csv, M3_agreement.csv, M4_within_station.csv, M5_subsampling.csv,
        M6_null_flips.csv, M7_persistence_4samples.csv, M8_fixed_effort_status.csv"""
import pandas as pd, numpy as np, sys
from core import B
rng = np.random.default_rng(20261005)
O = f'{B}/results/'
E = pd.read_csv(f'{B}/data/monthly/esri_timeseries_2020_2025.csv', low_memory=False)
E = E[E.agency == 'CPCB'].copy()
E['year'] = E.date_.str[:4].astype(int); E['month'] = E.date_.str[5:7].astype(int)
E['key'] = E.agency + '|' + E.station + '|' + E.lat.round(6).astype(str) + '|' + E.long.round(6).astype(str)
n_raw = len(E)
E = E.drop_duplicates(subset=[c for c in E.columns if c != 'objectid'])
n_exactdup = n_raw - len(E)
LK = pd.read_csv(f'{B}/data/monthly/esri_monthly_links.csv')
LK = LK[LK.status == 'linked']
E = E.merge(LK[['key', 'station_uid']], on='key')
A = pd.read_csv(f'{B}/analysis/analysis_r20.csv', low_memory=False)
A = A[A.in_panel_A_long & A.Year.isin([2020, 2021])]
TAB = A.set_index(['station_uid', 'Year'])

IND = {'BOD': ('bod_mgl', lambda v: float(np.max(v) > 3), 'BOD_max', 'ADV_BOD'),
       'FC': ('fe_col_mpn', lambda v: float(np.max(v) > 2500), 'FC_max', 'ADV_FC_cens'),   # FC censored (step 2)
       'PH': ('ph_mgl', lambda v: float((np.min(v) < 6.5) | (np.max(v) > 8.5)), 'PH_max', 'ADV_PH')}

def build_vals(E, collapse_months=False):
    vals = {}
    for (u, y), g in E.groupby(['station_uid', 'year']):
        for ind, (col, *_r) in IND.items():
            x = g[[col, 'month']].dropna()
            if collapse_months: x = x.groupby('month', as_index=False)[col].mean()
            if len(x): vals[(u, y, ind)] = (x[col].values.astype(float), x.month.values)
    return vals
vals = build_vals(E)
in_table = lambda u, y, ind: (u, y) in TAB.index and pd.notna(TAB.loc[(u, y), IND[ind][3]])

# ---------------- M1 reconciliation of counts
rows = []
sy_m = E.groupby(['station_uid', 'year']).size().reset_index()
sy_t = sy_m[[ (u, y) in TAB.index for u, y in zip(sy_m.station_uid, sy_m.year)]]
both_t = sy_t.groupby('station_uid').year.nunique()
rows += [dict(item='monthly records (CPCB agency, 2020-2021), after removing exact duplicates', value=len(E) + 0, note=f'{n_exactdup} exact duplicates removed'),
         dict(item='linked monitoring keys (one-to-one)', value=LK.key.nunique()),
         dict(item='linked station-years with any monthly record', value=len(sy_m)),
         dict(item='linked stations with any monthly record', value=sy_m.station_uid.nunique()),
         dict(item='linked station-years also in the CPCB table (Panel A_long)', value=len(sy_t)),
         dict(item='stations among them', value=sy_t.station_uid.nunique()),
         dict(item='stations in the CPCB table in both 2020 and 2021 with monthly records in both', value=int((both_t == 2).sum())),
         dict(item='check: station-years - stations', value=len(sy_t) - sy_t.station_uid.nunique())]
for ind in IND:
    allm = sorted({k[0] for k in vals if k[2] == ind and (k[0], 2020, ind) in vals and (k[0], 2021, ind) in vals
                   and len(vals[(k[0], 2020, ind)][0]) >= 4 and len(vals[(k[0], 2021, ind)][0]) >= 4})
    tab = [u for u in allm if in_table(u, 2020, ind) and in_table(u, 2021, ind)]
    rows += [dict(item=f'{ind}: stations with >= 4 monthly samples in both years (rev 19 benchmark set)', value=len(allm)),
             dict(item=f'{ind}: of these, in the CPCB table with {ind} status in both years (rev 20 benchmark set)', value=len(tab))]
M1 = pd.DataFrame(rows); M1.to_csv(O + 'M1_reconciliation.csv', index=False); print(M1.to_string(index=False))

# table-restricted value dictionary
tv = {k: v for k, v in vals.items() if in_table(k[0], k[1], k[2])}

# ---------------- M2 sample counts (BOD), table-linked station-years
rows = []
for (u, y), g in E.groupby(['station_uid', 'year']):
    if (u, y) not in TAB.index: continue
    rows.append(dict(station_uid=u, year=y, records=len(g), months=g.month.nunique(),
                     n_BOD=int(g.bod_mgl.notna().sum()), n_FC=int(g.fe_col_mpn.notna().sum()), n_PH=int(g.ph_mgl.notna().sum()),
                     months_with_two_records=int((g.groupby('month').size() > 1).sum())))
M2 = pd.DataFrame(rows); M2.to_csv(O + 'M2_counts.csv', index=False)
nb = M2.n_BOD[M2.n_BOD > 0]
print('BOD samples per station-year: n', len(nb), 'median', nb.median(), 'IQR', nb.quantile([.25, .75]).tolist(), 'range', nb.min(), nb.max(),
      '=12:', round((nb == 12).mean(), 3), '<10:', round((nb < 10).mean(), 3))
print('months per station-year: max', M2.months.max(), '; station-years with any month holding two records:', int((M2.months_with_two_records > 0).sum()),
      '; records > 12:', int((M2.records > 12).sum()))

# ---------------- M3 agreement of monthly extremes with the table
rows = []
for ind, (col, f, tcol, adv) in IND.items():
    keys = [k for k in tv if k[2] == ind]
    t = np.array([TAB.loc[(k[0], k[1]), tcol] for k in keys], dtype=float)
    m = np.array([tv[k][0].max() for k in keys])
    tol = 0.01 * np.clip(np.abs(t), 1, None) if ind == 'FC' else 0.051
    same = np.abs(t - m) <= tol
    if ind == 'PH':
        tmin = np.array([TAB.loc[(k[0], k[1]), 'PH_min'] for k in keys], dtype=float)
        same &= np.abs(tmin - np.array([tv[k][0].min() for k in keys])) <= 0.051
    st_t = np.array([TAB.loc[(k[0], k[1]), adv] for k in keys]); st_m = np.array([f(tv[k][0]) for k in keys])
    ok = ~np.isnan(t)
    rows.append(dict(indicator=ind, station_years=int(ok.sum()), stations=len({k[0] for k, o in zip(keys, ok) if o}),
                     value_identical=round(same[ok].mean(), 3), status_identical=round((st_t == st_m)[ok].mean(), 4)))
M3 = pd.DataFrame(rows); M3.to_csv(O + 'M3_agreement.csv', index=False); print(M3.to_string(index=False))

# ---------------- M4 within-station: year with more samples vs fewer
rows = []
for ind, (col, f, tcol, adv) in IND.items():
    st = sorted({k[0] for k in tv if k[2] == ind and (k[0], 2020, ind) in tv and (k[0], 2021, ind) in tv})
    hi = []; lo = []
    for u in st:
        n0, n1 = len(tv[(u, 2020, ind)][0]), len(tv[(u, 2021, ind)][0])
        if n0 == n1: continue
        ym, yl = (2020, 2021) if n0 > n1 else (2021, 2020)
        hi.append(TAB.loc[(u, ym), adv]); lo.append(TAB.loc[(u, yl), adv])
    hi, lo = np.array(hi), np.array(lo)
    rows.append(dict(indicator=ind, stations_both_years=len(st), stations_unequal_counts=len(hi),
                     adverse_more_samples=round(hi.mean(), 3), adverse_fewer_samples=round(lo.mean(), 3),
                     discordant_more_only=int(((hi == 1) & (lo == 0)).sum()), discordant_fewer_only=int(((hi == 0) & (lo == 1)).sum())))
M4 = pd.DataFrame(rows); M4.to_csv(O + 'M4_within_station.csv', index=False); print(M4.to_string(index=False))

# ---------------- M5 subsampling (station-years with >= 10 samples, table-linked)
rows = []; R = 500
for ind, (col, f, *_r) in IND.items():
    keys = [k for k in tv if k[2] == ind and len(tv[k][0]) >= 10]
    full = np.array([f(tv[k][0]) for k in keys])
    for kk in [1, 2, 3, 4, 6, 8, 'quarterly']:
        P = np.zeros(len(keys))
        for i, k in enumerate(keys):
            v, m = tv[k]; q = (m - 1) // 3; s = 0.0
            for _ in range(R):
                idx = [rng.choice(np.where(q == qq)[0]) for qq in np.unique(q)] if kk == 'quarterly' else rng.choice(len(v), kk, replace=False)
                s += f(v[idx])
            P[i] = s / R
        rows.append(dict(indicator=ind, samples_used=kk, station_years=len(keys), adverse_share_all=round(full.mean(), 3),
                         adverse_share_k=round(P.mean(), 3), missed_given_adverse=round(np.mean(1 - P[full == 1]), 3) if (full == 1).any() else np.nan))
M5 = pd.DataFrame(rows); M5.to_csv(O + 'M5_subsampling.csv', index=False); print(M5.to_string(index=False))

# ---------------- M6 null flips: pooled split (rev 19), months Jun-Dec only, month-matched swap; M7 four-sample persistence
def stay_entry(s0, s1):
    s0, s1 = np.asarray(s0), np.asarray(s1); a = s0 == 1; b = s0 == 0
    return (s1[a] == 1).mean() if a.any() else np.nan, (s1[b] == 1).mean() if b.any() else np.nan, int(a.sum()), int((s0 != s1).sum())

def benchmark(V, ind, label, months=None, R=1000):
    f = IND[ind][1]
    def sel(u, y):
        v, m = V[(u, y, ind)]
        if months is not None: k = np.isin(m, months); v, m = v[k], m[k]
        return v, m
    st = sorted({k[0] for k in V if k[2] == ind and (k[0], 2020, ind) in V and (k[0], 2021, ind) in V})
    st = [u for u in st if len(sel(u, 2020)[0]) >= 4 and len(sel(u, 2021)[0]) >= 4]
    s0 = np.array([f(sel(u, 2020)[0]) for u in st]); s1 = np.array([f(sel(u, 2021)[0]) for u in st])
    obs = stay_entry(s0, s1)
    out = dict(indicator=ind, variant=label, stations=len(st), adverse_2020=obs[2], observed_changes=obs[3], observed_P11=round(obs[0], 3))
    # pooled split null
    F = np.zeros((R, len(st)))
    for i, u in enumerate(st):
        v0, v1 = sel(u, 2020)[0], sel(u, 2021)[0]; pool = np.concatenate([v0, v1]); n0 = len(v0)
        for r in range(R):
            p = rng.permutation(pool); F[r, i] = f(p[:n0]) != f(p[n0:])
    tot = F.sum(1)
    out.update(pooled_expected=round(tot.mean(), 1), pooled_lo=float(np.percentile(tot, 2.5)), pooled_hi=float(np.percentile(tot, 97.5)),
               pooled_p_upper=round((1 + (tot >= obs[3]).sum()) / (R + 1), 4))
    # month-matched swap null: for each calendar month sampled in both years, swap the two years' values with p = 0.5
    F2 = np.zeros((R, len(st)))
    for i, u in enumerate(st):
        (v0, m0), (v1, m1) = sel(u, 2020), sel(u, 2021)
        common = np.intersect1d(m0, m1)
        for r in range(R):
            a, b = v0.copy(), v1.copy()
            for mm in common:
                if rng.random() < 0.5:
                    i0 = np.where(m0 == mm)[0]; i1 = np.where(m1 == mm)[0]
                    if len(i0) == 1 and len(i1) == 1: a[i0[0]], b[i1[0]] = v1[i1[0]], v0[i0[0]]
            F2[r, i] = f(a) != f(b)
    tot2 = F2.sum(1)
    out.update(swap_expected=round(tot2.mean(), 1), swap_lo=float(np.percentile(tot2, 2.5)), swap_hi=float(np.percentile(tot2, 97.5)),
               swap_p_upper=round((1 + (tot2 >= obs[3]).sum()) / (R + 1), 4))
    return out, st, s0, s1

rows6, rows7 = [], []
vals_c = {k: v for k, v in build_vals(E, collapse_months=True).items() if in_table(k[0], k[1], k[2])}
for ind in IND:
    for V, lab, months in [(tv, 'main: table-linked, all months', None), (tv, 'COVID-19: June-December only, both years', list(range(6, 13))),
                           (vals_c, 'two records in a month collapsed to their mean', None),
                           ({k: v for k, v in vals.items()}, 'rev 19 set (not table-restricted), for reference', None)]:
        out, st, s0, s1 = benchmark(V, ind, lab, months); rows6.append(out); print(out)
        if lab.startswith('main'):
            f = IND[ind][1]; sims = []
            for _ in range(1000):
                a = [f(rng.choice(tv[(u, 2020, ind)][0], 4, replace=False)) for u in st]
                b = [f(rng.choice(tv[(u, 2021, ind)][0], 4, replace=False)) for u in st]
                sims.append(stay_entry(a, b))
            sims = np.array(sims, dtype=float)
            rows7.append(dict(indicator=ind, stations=len(st), P11_all=round(stay_entry(s0, s1)[0], 3), P11_4=round(np.nanmean(sims[:, 0]), 3),
                              P11_4_lo=round(np.nanpercentile(sims[:, 0], 2.5), 3), P11_4_hi=round(np.nanpercentile(sims[:, 0], 97.5), 3),
                              changes_all=stay_entry(s0, s1)[3], changes_4=round(np.nanmean(sims[:, 3]), 1)))
M6 = pd.DataFrame(rows6); M6.to_csv(O + 'M6_null_flips.csv', index=False); print(M6.to_string(index=False))
M7 = pd.DataFrame(rows7); M7.to_csv(O + 'M7_persistence_4samples.csv', index=False); print(M7.to_string(index=False))
# export fixed-effort status draws for the ICC comparison in R (MC3b): status from one sample per quarter, 50 draws
rows = []
for ind, (col, f, *_r) in IND.items():
    st = sorted({k[0] for k in tv if k[2] == ind and (k[0], 2020, ind) in tv and (k[0], 2021, ind) in tv})
    for u in st:
        for y in (2020, 2021):
            v, m = tv[(u, y, ind)]; q = (m - 1) // 3
            r = dict(indicator=ind, station_uid=u, year=y, n=len(v), quarters=len(np.unique(q)), status_all=f(v))
            for dd in range(50):
                idx = [rng.choice(np.where(q == qq)[0]) for qq in np.unique(q)]
                r[f'q{dd}'] = f(v[idx])
            rows.append(r)
pd.DataFrame(rows).to_csv(O + 'M8_fixed_effort_status.csv', index=False)
