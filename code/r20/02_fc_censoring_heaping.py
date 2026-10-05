"""Rev 20, step 2: (Major Comment 3c) FC compliance that cannot be assessed is treated as censored, not compliant;
(Minor 6) value heaping at the criteria.
Censoring rule (primary in rev 20): FC status is set to missing (i) in every station-year of a state whose reported FC
maximum never exceeded 2,500 MPN/100 mL in 2016-2024 ('ceiling states'), and (ii) in any station-year whose FC maximum
is exactly 1,600 MPN/100 mL, the upper limit of an undiluted five-tube series, which can only be read as '>= 1,600'.
Sensitivity: also censor exactly 2,400 and 1,100 (three-tube / tray caps) and only rule (i).
Output: results/H1_heaping_state.csv, H2_heaping_year.csv, F1_ceiling_states.csv, F2_fc_censored_core.csv,
        F3_typology.csv; analysis/analysis_r20.csv gains ADV_FC_cens, ADV_MULTI_GE2_cens, ADV_ORG"""
import sys
from core import *
rng = np.random.default_rng(7)
d = pd.read_csv(f'{B}/analysis/analysis_r20.csv', low_memory=False)
A = d[d.in_panel_A]
num = lambda c: A[c].where(A[c + '_status'] == 'numeric')
# ---------- heaping
H = pd.DataFrame({'state': A.station_state, 'Year': A.Year,
                  'BOD_eq_3': (num('BOD_max') == 3.0).where(num('BOD_max').notna()),
                  'DO_eq_5': (num('DO_min') == 5.0).where(num('DO_min').notna()),
                  'FC_eq_1600': (num('FC_max') == 1600).where(num('FC_max').notna()),
                  'FC_eq_2400': (num('FC_max') == 2400).where(num('FC_max').notna()),
                  'FC_eq_2500': (num('FC_max') == 2500).where(num('FC_max').notna()),
                  'FC_eq_1100': (num('FC_max') == 1100).where(num('FC_max').notna())})
# neighbouring-value reference: share at 2.9 and 3.1 for BOD, 4.9 and 5.1 for DO
H['BOD_eq_2.9_or_3.1'] = num('BOD_max').isin([2.9, 3.1]).where(num('BOD_max').notna()) / 2
H['DO_eq_4.9_or_5.1'] = num('DO_min').isin([4.9, 5.1]).where(num('DO_min').notna()) / 2
cols = [c for c in H.columns if c not in ('state', 'Year')]
hs = H.groupby('state')[cols].mean().round(4); hs['station_years'] = H.groupby('state').size()
hs.loc['ALL'] = list(H[cols].mean().round(4)) + [len(H)]
hs.to_csv(f'{B}/results/H1_heaping_state.csv'); print(hs.to_string())
hy = H.groupby('Year')[cols].mean().round(4); hy.to_csv(f'{B}/results/H2_heaping_year.csv'); print(hy.to_string())

# ---------- ceiling states
fc = num('FC_max')
cs = pd.DataFrame({'state': A.station_state, 'fc': fc}).groupby('state').fc.agg(['count', 'max', lambda s: (s == 1600).mean()])
cs.columns = ['FC_station_years', 'FC_max_reported', 'share_at_1600']
cs['ceiling_state'] = (cs.FC_station_years > 0) & (cs.FC_max_reported <= 2500)
cs.loc['UNKNOWN', 'ceiling_state'] = False if 'UNKNOWN' in cs.index else None  # station without a resolvable state
cs.to_csv(f'{B}/results/F1_ceiling_states.csv'); print(cs.sort_values('FC_max_reported').to_string())
ceil = set(cs.index[cs.ceiling_state == True])
print('ceiling states:', len(ceil), sorted(ceil), 'station-years', int(A.station_state.isin(ceil).sum()))

def censor(df, rule):
    v = df.ADV_FC.copy(); f = df.FC_max.where(df.FC_max_status == 'numeric')
    if rule in ('main', 'states only', 'main + 2400/1100'): v[df.station_state.isin(ceil)] = np.nan
    if rule in ('main', 'main + 2400/1100'): v[f == 1600] = np.nan
    if rule == 'main + 2400/1100': v[f.isin([2400, 1100])] = np.nan
    return v
d['ADV_FC_cens'] = censor(d, 'main')
P4 = ['ADV_FC_cens', 'ADV_BOD', 'ADV_DO', 'ADV_PH']
allobs = d[P4].notna().all(axis=1)
d['ADV_MULTI_GE2_cens'] = np.where(allobs, (d[P4].sum(axis=1) >= 2).astype(float), np.nan)
d['ADV_ORG'] = np.where(d[['ADV_BOD', 'ADV_DO']].notna().any(axis=1), d[['ADV_BOD', 'ADV_DO']].max(axis=1), np.nan)
d.to_csv(f'{B}/analysis/analysis_r20.csv', index=False)

# ---------- FC core metrics under each rule
AL = d[d.in_panel_A_long]; st = AL.groupby('station_uid').station_state.first()
rows = []
for rule in ['none (rev 19)', 'main', 'states only', 'main + 2400/1100']:
    AL = AL.assign(FCx=AL.ADV_FC if rule.startswith('none') else censor(AL, rule))
    W = wide(AL, 'FCx'); M = W.values.astype(float); S = st.reindex(W.index).values
    p11, p01 = stay_entry(M)
    n1 = np.mean([stay_entry(null_N1(M, rng))[0] for _ in range(300)])
    n1s = np.mean([stay_entry(null_N1s(M, S, rng))[0] for _ in range(200)])
    n2 = np.mean([stay_entry(null_N2(M, rng))[0] for _ in range(300)])
    ci = cluster_boot(M, np.array(W.index), lambda X: stay_entry(X)[0], R=500, rng=rng)
    cis = cluster_boot(M, S, lambda X: stay_entry(X)[0], R=500, rng=rng)
    c, n = class_shares(M)
    rows.append(dict(rule=rule, stations=M.shape[0], station_years=int((~np.isnan(M)).sum()), states=len(set(S)),
                     prevalence=round(np.nanmean(M), 3), P11=round(p11, 3), CI_station=f'{ci[0]:.3f}-{ci[1]:.3f}', CI_state=f'{cis[0]:.3f}-{cis[1]:.3f}',
                     P01=round(p01, 3), N1=round(n1, 3), N1s=round(n1s, 3), N2=round(n2, 3), obs_minus_N2=round(p11 - n2, 3),
                     share_beyond_state=round((p11 - n1s) / (p11 - n1), 3), jaccard1=round(jaccard_lag(M, 1)[0], 3),
                     classified=n, persistent=round(c.get('persistent', 0), 3), never=round(c.get('never', 0), 3)))
    print(rows[-1])
pd.DataFrame(rows).to_csv(f'{B}/results/F2_fc_censored_core.csv', index=False)

# ---------- typology (microbial vs organic), uncensored vs censored
rows = []
for lab, col in [('rev 19 (ceiling coded compliant)', 'ADV_FC'), ('rev 20 (censored)', 'ADV_FC_cens')]:
    Wf = wide(AL, col, 4); Wo = wide(AL, 'ADV_ORG', 4)
    idx = Wf.index.intersection(Wo.index); Mf = Wf.loc[idx].values.astype(float); Mo = Wo.loc[idx].values.astype(float)
    rf = np.nanmean(Mf, 1); ro = np.nanmean(Mo, 1)
    pf = np.array([classify(r) == 'persistent' for r in Mf]); po = np.array([classify(r) == 'persistent' for r in Mo])
    from scipy.stats import pearsonr, spearmanr
    rows.append(dict(version=lab, stations=len(idx), jaccard_persistent=round((pf & po).sum() / max((pf | po).sum(), 1), 3),
                     corr_recurrence=round(pearsonr(rf, ro).statistic, 3), spearman=round(spearmanr(rf, ro).statistic, 3),
                     persistentFC_also_org=round((pf & po).sum() / pf.sum(), 3), persistentORG_also_FC=round((pf & po).sum() / po.sum(), 3),
                     multi_domain=round(np.mean((rf >= .5) & (ro >= .5)), 3), predominantly_microbial=round(np.mean((rf >= .5) & (ro < .5)), 3),
                     predominantly_organic=round(np.mean((rf < .5) & (ro >= .5)), 3), neither=round(np.mean((rf < .5) & (ro < .5)), 3)))
    print(rows[-1])
pd.DataFrame(rows).to_csv(f'{B}/results/F3_typology.csv', index=False)
