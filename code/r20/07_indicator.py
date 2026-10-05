"""Rev 20, step 7 (Major Comments 1 and 2): a multi-year, error-aware status indicator, tested out of sample, and a
comparison of compliance rules in the monthly subset.
Year level. Each station-year gets a signed distance from the criterion in error units, z (positive = adverse side):
  BOD z = (ln x - ln 3)/0.143 (CV 15.4%); FC z = (log10 x - log10 2500)/0.26 (FC censored as in step 2);
  pH z = max((6.5 - min), (max - 8.5))/0.07; DO z = (5 - min)/0.30 (assumed field-level error; DO is secondary here).
  A year is confidently adverse if z > 1.645, confidently compliant if z < -1.645, otherwise within the error band.
  BDL counts as confidently compliant (BOD, FC) or confidently adverse (DO).
Station level, over a window with >= 3 observed years (three-class indicator, 'TCI'):
  reliably adverse   confidently adverse in >= 60% of observed years;
  reliably compliant confidently compliant in >= 60% of observed years and never confidently adverse;
  borderline         all other stations.
Out-of-sample test: classes from 2016-2020, outcomes 2021-2024 (single-year status as published by the current rule).
Listing churn: a rolling three-year TCI list (reliably adverse) against the single-year list (adverse this year).
Compliance rules (monthly subset 2020-2021, stations in the CPCB table both years): single-sample rule (annual maximum
beyond the criterion, current practice), raw-score rule (> 10% of samples beyond), binomial rule (list when
P[Bin(n, 0.10) >= k] <= 0.05, Smith et al. 2001), median rule (> 50% of samples beyond).
Output: results/I1_tci_out_of_sample.csv, I2_tci_classes.csv, I3_listing_churn.csv, I4_rules_monthly.csv"""
import sys, os
from core import *
from scipy.stats import binom
rng = np.random.default_rng(3)
d = load('A_long')
num = lambda c: d[c].where(d[c + '_status'] == 'numeric')
Z = pd.DataFrame(index=d.index)
Z['BOD'] = np.where(d.BOD_max_status == 'BDL', -9, (np.log(num('BOD_max').clip(lower=1e-6)) - np.log(3)) / 0.143)
fc = num('FC_max').clip(lower=1e-6)
Z['FC'] = np.where(d.FC_max_status == 'BDL', -9, (np.log10(fc) - np.log10(2500)) / 0.26)
Z.loc[d.ADV_FC_cens.isna(), 'FC'] = np.nan
Z['PH'] = np.fmax((6.5 - num('PH_min')) / 0.07, (num('PH_max') - 8.5) / 0.07)
Z.loc[d.ADV_PH.isna(), 'PH'] = np.nan
Z['DO'] = np.where(d.DO_min_status == 'BDL', 9, (5 - num('DO_min')) / 0.30)
for c in Z: d[f'z_{c}'] = Z[c]
STATUS = {'BOD': 'ADV_BOD', 'FC': 'ADV_FC_cens', 'PH': 'ADV_PH', 'DO': 'ADV_DO'}

def tci(zrow):
    v = zrow[~np.isnan(zrow)]
    if len(v) < 3: return None
    a = (v > 1.645).mean(); c = (v < -1.645).mean()
    if a >= 0.6: return 'reliably adverse'
    if c >= 0.6 and (v > 1.645).sum() == 0: return 'reliably compliant'
    return 'borderline'

def persistent_rule(srow):
    v = srow[~np.isnan(srow)]
    if len(v) < 3: return None
    return 'adverse majority' if v.mean() >= 0.6 else ('never adverse' if v.sum() == 0 else 'other')

rows1, rows2 = [], []
for ind, adv in STATUS.items():
    Zw = d.pivot_table(index='station_uid', columns='Year', values=f'z_{ind}', aggfunc='first').reindex(columns=YEARS)
    Sw = d.pivot_table(index='station_uid', columns='Year', values=adv, aggfunc='first').reindex(index=Zw.index, columns=YEARS)
    early, late = [2016, 2017, 2018, 2019, 2020], [2021, 2022, 2023, 2024]
    cls_e = Zw[early].apply(lambda r: tci(r.values.astype(float)), axis=1)
    cls_l = Zw[late].apply(lambda r: tci(r.values.astype(float)), axis=1)
    last = Sw[2020]
    out_late = Sw[late]
    keep = cls_e.notna() & (out_late.notna().sum(1) >= 2)
    for cname in ['reliably adverse', 'borderline', 'reliably compliant']:
        k = keep & (cls_e == cname)
        y = out_late[k].values.astype(float)
        rows1.append(dict(indicator=ind, class_2016_2020=cname, stations=int(k.sum()), share_of_classified=round(k.sum() / keep.sum(), 3),
                          station_years_2021_2024=int((~np.isnan(y)).sum()), adverse_share_2021_2024=round(np.nanmean(y), 3),
                          stations_always_adverse_later=round(np.mean([np.nanmin(r) == 1 for r in y]), 3),
                          stations_never_adverse_later=round(np.mean([np.nanmax(r) == 0 for r in y]), 3),
                          same_class_2021_2024=round((cls_l[k] == cname).mean(), 3) if cname else np.nan,
                          class_2021_2024_missing=int(cls_l[k].isna().sum())))
    # single-year rule as a predictor: status in 2020 -> 2021-2024 single-year status
    k = keep & last.notna()
    y = out_late[k].values.astype(float); p = np.repeat(last[k].values[:, None], 4, 1)
    ok = ~np.isnan(y)
    acc_single = (y[ok] == p[ok]).mean()
    # multi-year majority rule (60% of 2016-2020 years adverse), decisive for all stations
    maj = (Sw[early].mean(1) >= 0.6).astype(float)[k].values
    acc_maj = (y[ok] == np.repeat(maj[:, None], 4, 1)[ok]).mean()
    # TCI: reliably classes predicted, borderline predicted by majority but reported separately
    ck = cls_e[k].values
    rel = np.isin(ck, ['reliably adverse', 'reliably compliant'])
    pt = np.repeat((ck == 'reliably adverse').astype(float)[:, None], 4, 1)
    okr = ok & rel[:, None]
    rows2.append(dict(indicator=ind, stations=int(k.sum()), station_years=int(ok.sum()),
                      accuracy_single_year_2020=round(acc_single, 3), accuracy_majority_2016_2020=round(acc_maj, 3),
                      tci_reliable_station_share=round(rel.mean(), 3), accuracy_tci_reliable=round((y[okr] == pt[okr]).mean(), 3),
                      accuracy_single_year_on_borderline=round((y[ok & ~rel[:, None]] == p[ok & ~rel[:, None]]).mean(), 3),
                      share_late_changes_at_borderline=round(np.nansum(np.abs(np.diff(y, axis=1))[~rel]) / max(np.nansum(np.abs(np.diff(y, axis=1))), 1), 3)))
    print(rows2[-1])
I1 = pd.DataFrame(rows1); I1.to_csv(f'{B}/results/I1_tci_out_of_sample.csv', index=False); print(I1.to_string(index=False))
I2 = pd.DataFrame(rows2); I2.to_csv(f'{B}/results/I2_tci_prediction.csv', index=False); print(I2.to_string(index=False))

# listing churn: single-year list vs rolling three-year TCI list (BOD and FC), stations observed in both compared windows
rows = []
for ind, adv in STATUS.items():
    Zw = d.pivot_table(index='station_uid', columns='Year', values=f'z_{ind}', aggfunc='first').reindex(columns=YEARS)
    Sw = d.pivot_table(index='station_uid', columns='Year', values=adv, aggfunc='first').reindex(index=Zw.index, columns=YEARS)
    ent_s = ex_s = size_s = 0; ent_t = ex_t = size_t = 0
    for t in range(2018, 2024):
        a, b = Sw[t], Sw[t + 1]; ok = a.notna() & b.notna()
        ent_s += int(((a == 0) & (b == 1) & ok).sum()); ex_s += int(((a == 1) & (b == 0) & ok).sum()); size_s += int(((a == 1) & ok).sum())
        w0 = Zw[[t - 2, t - 1, t]].apply(lambda r: tci(r.values.astype(float)), axis=1)
        w1 = Zw[[t - 1, t, t + 1]].apply(lambda r: tci(r.values.astype(float)), axis=1)
        ok2 = w0.notna() & w1.notna()
        A0 = (w0 == 'reliably adverse') & ok2; A1 = (w1 == 'reliably adverse') & ok2
        ent_t += int((~A0 & A1).sum()); ex_t += int((A0 & ~A1).sum()); size_t += int(A0.sum())
    rows.append(dict(indicator=ind, single_year_list_size=size_s, single_year_entries=ent_s, single_year_exits=ex_s,
                     single_year_turnover=round((ent_s + ex_s) / size_s, 3),
                     tci_list_size=size_t, tci_entries=ent_t, tci_exits=ex_t, tci_turnover=round((ent_t + ex_t) / max(size_t, 1), 3)))
I3 = pd.DataFrame(rows); I3.to_csv(f'{B}/results/I3_listing_churn.csv', index=False); print(I3.to_string(index=False))

d[['station_uid', 'Year'] + [f'z_{c}' for c in Z]].to_csv(f'{B}/results/I0_station_year_z.csv', index=False)

# compliance rules in the monthly subset (Major Comment 2)
MONTHLY = f'{B}/data/monthly/esri_timeseries_2020_2025.csv'
if not os.path.exists(MONTHLY):
    print('Monthly records not found (see data/README.md); results/I4_rules_monthly.csv is left as deposited.')
    sys.exit(0)
E = pd.read_csv(MONTHLY, low_memory=False)
E = E[E.agency == 'CPCB'].copy(); E = E.drop_duplicates(subset=[c for c in E.columns if c != 'objectid'])
E['year'] = E.date_.str[:4].astype(int); E['month'] = E.date_.str[5:7].astype(int)
E['key'] = E.agency + '|' + E.station + '|' + E.lat.round(6).astype(str) + '|' + E.long.round(6).astype(str)
LK = pd.read_csv(f'{B}/data/monthly/esri_monthly_links.csv'); LK = LK[LK.status == 'linked']
E = E.merge(LK[['key', 'station_uid']], on='key')
T = d[d.Year.isin([2020, 2021])].set_index(['station_uid', 'Year'])
EXC = {'BOD': ('bod_mgl', lambda v: v > 3, 'ADV_BOD'), 'FC': ('fe_col_mpn', lambda v: v > 2500, 'ADV_FC_cens'),
       'PH': ('ph_mgl', lambda v: (v < 6.5) | (v > 8.5), 'ADV_PH')}
RULES = {'single sample (current)': lambda k, n: k >= 1, 'raw score > 10%': lambda k, n: k / n > 0.10,
         'binomial, 10% exceedance, alpha 0.05': lambda k, n: binom.sf(k - 1, n, 0.10) <= 0.05, 'median (> 50%)': lambda k, n: k / n > 0.5}
rows = []
for ind, (col, exc, adv) in EXC.items():
    recs = {}
    for (u, y), g in E.groupby(['station_uid', 'year']):
        v = g[col].dropna().values
        if len(v) >= 4 and (u, y) in T.index and pd.notna(T.loc[(u, y), adv]): recs[(u, y)] = v
    st = sorted({u for u, y in recs if (u, 2020) in recs and (u, 2021) in recs})
    for rname, rule in RULES.items():
        s = {(u, y): float(rule(exc(recs[(u, y)]).sum(), len(recs[(u, y)]))) for u in st for y in (2020, 2021)}
        s0 = np.array([s[(u, 2020)] for u in st]); s1 = np.array([s[(u, 2021)] for u in st])
        # sampling null for this rule: pooled split of the two years' samples
        tot = []
        for _ in range(500):
            c = 0
            for u in st:
                pool = np.concatenate([recs[(u, 2020)], recs[(u, 2021)]]); n0 = len(recs[(u, 2020)]); p = rng.permutation(pool)
                a = rule(exc(p[:n0]).sum(), n0); b = rule(exc(p[n0:]).sum(), len(pool) - n0); c += a != b
            tot.append(c)
        tot = np.array(tot)
        single = np.array([float(exc(recs[(u, y)]).sum() >= 1) for u in st for y in (2020, 2021)])
        this = np.array([s[(u, y)] for u in st for y in (2020, 2021)])
        rows.append(dict(indicator=ind, rule=rname, stations=len(st), station_years=2 * len(st), adverse_share=round(this.mean(), 3),
                         classifications_differing_from_single_sample=int((this != single).sum()),
                         changes_2020_2021=int((s0 != s1).sum()), expected_changes_no_change_null=round(tot.mean(), 1),
                         null_lo=float(np.percentile(tot, 2.5)), null_hi=float(np.percentile(tot, 97.5))))
        print(rows[-1])
I4 = pd.DataFrame(rows); I4.to_csv(f'{B}/results/I4_rules_monthly.csv', index=False); print(I4.to_string(index=False))
