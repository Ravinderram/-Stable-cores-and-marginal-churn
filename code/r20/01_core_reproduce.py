"""Rev 20, step 1: reproduce rev 19 Table 3 from the rebuilt dataset and add (Minor 4) observed - N2 by panel and
(Minor 5) station-, state- and basin-clustered intervals. Output: results/T3_core.csv, T3_by_panel.csv"""
import sys
from core import *
import os; os.makedirs(f'{B}/results', exist_ok=True)
rng = np.random.default_rng(20261005)
d = load(None)
A = d[d.in_panel_A_long]
st = A.groupby('station_uid').station_state.first()
bas = A.groupby('station_uid').hybas_lev06.first()
rows = []
for ind in ['ADV_FC_cens', 'ADV_BOD', 'ADV_DO', 'ADV_PH', 'ADV_MULTI_GE2_cens', 'ADV_FC']:
    W = wide(A, ind); M = W.values.astype(float); S = st.reindex(W.index).values
    Bs = bas.reindex(W.index); Bs = np.where(Bs.isna(), 'NB_' + pd.Series(W.index).values, Bs.astype(str)).astype(str)
    p11, p01 = stay_entry(M)
    R = 1000
    n1 = np.array([stay_entry(null_N1(M, rng))[0] for _ in range(R)])
    n2 = np.array([stay_entry(null_N2(M, rng))[0] for _ in range(R)])
    n1s = np.array([stay_entry(null_N1s(M, S, rng))[0] for _ in range(300)])
    ci_st = cluster_boot(M, np.array(W.index), lambda X: stay_entry(X)[0], R=1000, rng=rng)
    ci_state = cluster_boot(M, S, lambda X: stay_entry(X)[0], R=1000, rng=rng)
    ci_basin = cluster_boot(M, Bs, lambda X: stay_entry(X)[0], R=1000, rng=rng)
    j1 = jaccard_lag(M, 1); cs, ncl = class_shares(M)
    rows.append(dict(indicator=LABEL.get(ind, ind), stations=M.shape[0], station_years=int((~np.isnan(M)).sum()),
                     prevalence=round(np.nanmean(M), 3), P11=round(p11, 3), P01=round(p01, 3),
                     CI_station=f'{ci_st[0]:.3f}-{ci_st[1]:.3f}', CI_state=f'{ci_state[0]:.3f}-{ci_state[1]:.3f}',
                     CI_subbasin=f'{ci_basin[0]:.3f}-{ci_basin[1]:.3f}',
                     N1=round(n1.mean(), 3), N1s=round(n1s.mean(), 3), N2=round(n2.mean(), 3), N2_hi=round(np.percentile(n2, 97.5), 3),
                     obs_minus_N2=round(p11 - n2.mean(), 3), share_beyond_state=round((p11 - n1s.mean()) / (p11 - n1.mean()), 3),
                     jaccard1_mean=round(j1[0], 3), jaccard1_pooled=round(j1[1], 3),
                     persistent=round(cs.get('persistent', 0), 3), classified=ncl,
                     never=round(cs.get('never', 0), 3), emerging=round(cs.get('emerging', 0), 3), diminishing=round(cs.get('diminishing', 0), 3)))
    print(rows[-1])
T = pd.DataFrame(rows); T.to_csv(f'{B}/results/T3_core.csv', index=False)
print(T.to_string(index=False))

# Minor 4: observed - N2 by panel (stations with >= k observed years)
rows = []
for ind in ['ADV_FC_cens', 'ADV_BOD', 'ADV_DO', 'ADV_PH', 'ADV_MULTI_GE2_cens', 'ADV_FC']:
    Wall = wide(A, ind)
    for k in [2, 4, 5, 6, 9]:
        W = Wall[Wall.notna().sum(1) >= k]; M = W.values.astype(float)
        p11 = stay_entry(M)[0]
        n2 = np.array([stay_entry(null_N2(M, rng))[0] for _ in range(500)])
        rows.append(dict(indicator=LABEL.get(ind, ind), min_years=k, stations=M.shape[0], P11=round(p11, 3), N2=round(n2.mean(), 3),
                         obs_minus_N2=round(p11 - n2.mean(), 3), N2_95=f'{np.percentile(n2, 2.5):.3f}-{np.percentile(n2, 97.5):.3f}'))
TP = pd.DataFrame(rows); TP.to_csv(f'{B}/results/T3_by_panel.csv', index=False)
print(TP.to_string(index=False))
