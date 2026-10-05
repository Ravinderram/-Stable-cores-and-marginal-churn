"""Rev 20, step 11: data behind the figures that are not produced by other steps.
G1 Jaccard overlap by lag (observed, N1, N2; 200 permutations), G2 probability of a status change by the value in
year t (25 quantile bins), G3 station classes and the persistent share across 12 rules, G4 coverage by year/first year.
Output: results/G1_overlap_lag.csv, G2_threshold_proximity.csv, G3_classes.csv, G4_coverage.csv"""
import sys
from core import *
rng = np.random.default_rng(23)
d = load('A_long'); dA = load('A')
IND = ['ADV_FC_cens', 'ADV_BOD', 'ADV_DO', 'ADV_PH', 'ADV_MULTI_GE2_cens']
NAME = {'ADV_FC_cens': 'FC', 'ADV_BOD': 'BOD', 'ADV_DO': 'DO', 'ADV_PH': 'pH', 'ADV_MULTI_GE2_cens': '>=2 of 4'}
rows = []
for ind in IND:
    M = wide(d, ind).values.astype(float)
    n1 = [null_N1(M, rng) for _ in range(200)]; n2 = [null_N2(M, rng) for _ in range(200)]
    for lag in range(1, 9):
        rows.append(dict(indicator=NAME[ind], lag=lag, observed=jaccard_lag(M, lag)[0], N1=np.mean([jaccard_lag(X, lag)[0] for X in n1]),
                         N2=np.mean([jaccard_lag(X, lag)[0] for X in n2]), N2_lo=np.percentile([jaccard_lag(X, lag)[0] for X in n2], 2.5),
                         N2_hi=np.percentile([jaccard_lag(X, lag)[0] for X in n2], 97.5)))
pd.DataFrame(rows).to_csv(f'{B}/results/G1_overlap_lag.csv', index=False)
# G2 threshold proximity
num = lambda c: d[c].where(d[c + '_status'] == 'numeric')
rows = []
for nm, col, adv, crit in [('BOD', 'BOD_max', 'ADV_BOD', 3.0), ('DO', 'DO_min', 'ADV_DO', 5.0), ('FC', 'FC_max', 'ADV_FC_cens', 2500)]:
    V = d.assign(v=num(col)).pivot_table(index='station_uid', columns='Year', values='v', aggfunc='first').reindex(columns=YEARS)
    S = d.pivot_table(index='station_uid', columns='Year', values=adv, aggfunc='first').reindex(index=V.index, columns=YEARS)
    x = V.values[:, :-1].ravel(); a = S.values[:, :-1].ravel(); b = S.values[:, 1:].ravel()
    ok = ~np.isnan(x) & ~np.isnan(a) & ~np.isnan(b); x, ch = x[ok], (a[ok] != b[ok]).astype(float)
    qs = np.unique(np.quantile(x, np.linspace(0, 1, 26)))
    bins = np.clip(np.searchsorted(qs, x, side='right') - 1, 0, len(qs) - 2)
    for k in range(len(qs) - 1):
        m = bins == k
        if m.sum() < 10: continue
        p = ch[m].mean(); n = m.sum(); se = np.sqrt(p * (1 - p) / n)
        rows.append(dict(indicator=nm, criterion=crit, bin_lo=qs[k], bin_hi=qs[k + 1], x_median=np.median(x[m]), n=int(n), p_change=p,
                         lo=max(0, p - 1.96 * se), hi=min(1, p + 1.96 * se)))
pd.DataFrame(rows).to_csv(f'{B}/results/G2_threshold_proximity.csv', index=False)
# G3 classes + rule range
rows = []
for ind in IND[:4]:
    M = wide(d, ind).values.astype(float); c, n = class_shares(M)
    rng_p = [class_shares(M, rec=r, run=k)[0].get('persistent', 0) for r in (0.5, 0.6, 0.67, 0.75) for k in (2, 3, 4)]
    rows.append(dict(indicator=NAME[ind], classified=n, **{k: c.get(k, 0) for k in ['persistent', 'intermittent', 'diminishing', 'emerging', 'never']},
                     persistent_min=min(rng_p), persistent_max=max(rng_p)))
G3 = pd.DataFrame(rows); G3.to_csv(f'{B}/results/G3_classes.csv', index=False); print(G3.round(3).to_string(index=False))
# G4 coverage
first = dA.groupby('station_uid').Year.min()
cov = dA.assign(first=dA.station_uid.map(first)).groupby(['Year', 'first']).size().unstack(fill_value=0)
cov.to_csv(f'{B}/results/G4_coverage.csv')
nobs = dA.groupby('station_uid').Year.nunique().value_counts().sort_index(); nobs.to_csv(f'{B}/results/G4b_years_per_station.csv')
print(cov.sum(1).to_dict(), nobs.to_dict())
