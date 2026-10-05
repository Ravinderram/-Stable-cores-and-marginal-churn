"""Rev 20, step 5 (Major Comment 3a, Minor 12): persistence by reported sampling frequency (2023 snapshot) and the
2020-2021 year pair in the whole network versus the monthly subset.
reported_frequency comes from the CPCB 2023 table republished by Esri India (monthly / quarterly / yearly); it is a
station-level attribute observed for located stations only and refers to 2023.
Output: results/R1_persistence_by_frequency.csv, R2_pair_2020_2021.csv"""
import sys
from core import *
rng = np.random.default_rng(5)
d = load('A_long')
fr = d.groupby('station_uid').reported_frequency.agg(lambda s: s.dropna().mode().iat[0] if s.notna().any() else 'not reported')
st = d.groupby('station_uid').station_state.first()
rows = []
for ind in ['ADV_FC_cens', 'ADV_FC', 'ADV_BOD', 'ADV_DO', 'ADV_PH']:
    W = wide(d, ind)
    for f in ['monthly', 'quarterly', 'yearly', 'not reported']:
        Wf = W[fr.reindex(W.index).values == f]
        if len(Wf) < 10: continue
        M = Wf.values.astype(float); p11, p01 = stay_entry(M)
        n2 = np.mean([stay_entry(null_N2(M, rng))[0] for _ in range(300)])
        ci = cluster_boot(M, np.array(Wf.index), lambda X: stay_entry(X)[0], R=500, rng=rng)
        rows.append(dict(indicator=ind.replace('ADV_', ''), frequency_2023=f, stations=len(Wf), states=st.reindex(Wf.index).nunique(),
                         prevalence=round(np.nanmean(M), 3), P11=round(p11, 3), CI=f'{ci[0]:.2f}-{ci[1]:.2f}', P01=round(p01, 3),
                         N2=round(n2, 3), obs_minus_N2=round(p11 - n2, 3)))
R1 = pd.DataFrame(rows); R1.to_csv(f'{B}/results/R1_persistence_by_frequency.csv', index=False); print(R1.to_string(index=False))

# Minor 12: 2020-2021 pair, network vs monthly subset states
lk = pd.read_csv(f'{B}/results/M2_counts.csv').station_uid.unique()
lk_states = st.reindex(lk).value_counts()
rows = []
for ind in ['ADV_BOD', 'ADV_FC', 'ADV_PH']:
    W = wide(d, ind, 1)[[2020, 2021]].dropna(); a, b = W[2020].values, W[2021].values
    def pr(mask): 
        aa, bb = a[mask], b[mask]; return round((bb[aa == 1] == 1).mean(), 3), int((aa == 1).sum())
    for lab, mask in [('network, 2020-2021', np.ones(len(W), bool)),
                      ('network, all pairs 2016-2024', None),
                      ('states of the monthly subset, 2020-2021', st.reindex(W.index).isin(lk_states.index).values),
                      ('linked stations, 2020-2021 (table status)', W.index.isin(lk))]:
        if mask is None:
            p = stay_entry(wide(d, ind).values.astype(float))[0]; rows.append(dict(indicator=ind, set=lab, P11=round(p, 3))); continue
        p, n = pr(mask); rows.append(dict(indicator=ind, set=lab, P11=p, adverse_in_2020=n))
R2 = pd.DataFrame(rows); R2.to_csv(f'{B}/results/R2_pair_2020_2021.csv', index=False); print(R2.to_string(index=False))
print(lk_states.head(10).to_dict())
# P11 by year pair for BOD
W = wide(d, 'ADV_BOD').values.astype(float)
print('BOD P11 by year pair:', [round(np.nansum((W[:, t] == 1) & (W[:, t + 1] == 1)) / np.nansum((W[:, t] == 1) & ~np.isnan(W[:, t + 1])), 3) for t in range(8)])
