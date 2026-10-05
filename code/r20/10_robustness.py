"""Rev 20, step 10 (Major Comment 7, Fig. 11): full robustness grid; every core metric for every variant.
Variants (base = Panel A_long, regulatory thresholds, FC censored):
 thresholds: DO < 4; pH 6-9; FC > 500; values at the limit adverse; empirical P80; empirical P90 (DO, BOD, FC);
 FC coded as in rev 19 (ceiling states and 1,600 compliant);
 minimum observed years 4, 5, 6, 9; harmonization: low-confidence and ambiguous stations excluded, all-high-confidence
 stations only, raw CPCB codes; years dropped: 2016, 2017, 2020, 2016-2017; one-year gaps bridged.
Metrics: prevalence, P(1|1), P(1|0), N2 mean P(1|1) (200 permutations), observed - N2, Jaccard next year,
persistent share (60% + run of 3, >= 4 observed years), and the multi-parameter contrast (P11 with another primary
indicator adverse minus P11 alone).
Output: results/S_robustness_full.csv"""
import sys
from core import *
rng = np.random.default_rng(19)
d0 = load(None)
num = lambda D, c: D[c].where(D[c + '_status'] == 'numeric')
IND = ['ADV_FC_cens', 'ADV_BOD', 'ADV_DO', 'ADV_PH']

def recode(D, kind):
    D = D.copy()
    if kind == 'DO < 4': D['ADV_DO'] = np.where(D.DO_min_status == 'BDL', 1, np.where(num(D, 'DO_min').notna(), (num(D, 'DO_min') < 4).astype(float), np.nan))
    if kind == 'pH 6-9':
        lo, hi = num(D, 'PH_min'), num(D, 'PH_max')
        D['ADV_PH'] = np.where((lo < 6) | (hi > 9), 1.0, np.where(lo.notna() & hi.notna(), 0.0, np.nan))
    if kind == 'FC > 500':
        f = num(D, 'FC_max'); v = np.where(D.FC_max_status == 'BDL', 0, np.where(f.notna(), (f > 500).astype(float), np.nan))
        D['ADV_FC_cens'] = np.where(D.ADV_FC_cens.isna(), np.nan, v)
    if kind == 'values at the limit adverse':
        D['ADV_BOD'] = np.where(D.BOD_max_status == 'BDL', 0, np.where(num(D, 'BOD_max').notna(), (num(D, 'BOD_max') >= 3).astype(float), np.nan))
        D['ADV_DO'] = np.where(D.DO_min_status == 'BDL', 1, np.where(num(D, 'DO_min').notna(), (num(D, 'DO_min') <= 5).astype(float), np.nan))
        f = num(D, 'FC_max'); v = np.where(D.FC_max_status == 'BDL', 0, np.where(f.notna(), (f >= 2500).astype(float), np.nan))
        D['ADV_FC_cens'] = np.where(D.ADV_FC_cens.isna(), np.nan, v)
        lo, hi = num(D, 'PH_min'), num(D, 'PH_max')
        D['ADV_PH'] = np.where((lo <= 6.5) | (hi >= 8.5), 1.0, np.where(lo.notna() & hi.notna(), 0.0, np.nan))
    if kind in ('empirical P80', 'empirical P90'):
        q = 0.8 if kind.endswith('80') else 0.9; A = D[D.in_panel_A]
        for col, adv, low in [('BOD_max', 'ADV_BOD', False), ('DO_min', 'ADV_DO', True), ('FC_max', 'ADV_FC_cens', False)]:
            v = num(D, col); ref = num(A, col)
            t = ref.quantile(1 - q) if low else ref.quantile(q)
            s = (v < t) if low else (v > t)
            new = np.where(v.notna(), s.astype(float), np.where(D[col + '_status'] == 'BDL', 1.0 if low else 0.0, np.nan))
            D[adv] = np.where(D[adv].isna() & (adv == 'ADV_FC_cens'), np.nan, new) if adv == 'ADV_FC_cens' else new
    if kind == 'FC coded as in rev 19': D['ADV_FC_cens'] = D.ADV_FC
    return D

def panel(kind):
    D = d0[d0.in_panel_A_long].copy()
    if kind.startswith('min years'):
        k = int(kind.split()[-1]); n = D.groupby('station_uid').Year.nunique(); D = D[D.station_uid.map(n) >= k]
    if kind == 'harmonization: low/ambiguous excluded': D = D[~D.harmonization_confidence.isin(['low', 'ambiguous'])]
    if kind == 'harmonization: all-high stations only':
        ok = D.groupby('station_uid').harmonization_confidence.agg(lambda s: (s == 'high').all()); D = D[D.station_uid.map(ok)]
    if kind == 'raw CPCB codes':
        D = d0[d0.in_panel_A & d0['Station Code'].notna()].copy()
        D['station_uid'] = 'C' + D['Station Code'].astype(int).astype(str)
        D = D.sort_values('row_id').drop_duplicates(['station_uid', 'Year'])
        n = D.groupby('station_uid').Year.nunique(); D = D[D.station_uid.map(n) >= 2]
    if kind.startswith('drop '):
        ys = [int(y) for y in kind.replace('drop ', '').split('-')]; ys = list(range(ys[0], ys[-1] + 1))
        D = D[~D.Year.isin(ys)]
    return D

def bridge(M):
    """carry transitions across single missing years: compress each row's observed values with gap <= 2 years"""
    X = M.copy()
    for i in range(X.shape[0]):
        for t in range(1, X.shape[1] - 1):
            if np.isnan(X[i, t]) and not np.isnan(X[i, t - 1]) and not np.isnan(X[i, t + 1]):
                pass
    return X

def bridged_counts(M):
    n10 = n11 = n01 = n00 = 0
    for r in M:
        idx = np.where(~np.isnan(r))[0]
        for a, b in zip(idx[:-1], idx[1:]):
            if b - a <= 2:
                x, y = r[a], r[b]
                n11 += x == 1 and y == 1; n10 += x == 1 and y == 0; n01 += x == 0 and y == 1; n00 += x == 0 and y == 0
    return n11 / max(n11 + n10, 1), n01 / max(n01 + n00, 1)

def metrics(D, ind, bridge_gaps=False):
    W = wide(D, ind); M = W.values.astype(float)
    if bridge_gaps:
        p11, p01 = bridged_counts(M)
        n2 = np.mean([bridged_counts(null_N2(M, rng))[0] for _ in range(100)])
    else:
        p11, p01 = stay_entry(M); n2 = np.mean([stay_entry(null_N2(M, rng))[0] for _ in range(200)])
    c, ncl = class_shares(M)
    # multi-parameter contrast: P11 when another primary indicator adverse in year t vs not
    others = [o for o in IND if o != ind]
    Wo = sum(wide(D, o, 1).reindex(index=W.index, columns=YEARS).fillna(0).values for o in others) >= 1
    a, b = M[:, :-1], M[:, 1:]; ok = (a == 1) & ~np.isnan(b)
    withm = ok & Wo[:, :-1]; alone = ok & ~Wo[:, :-1]
    return dict(stations=M.shape[0], station_years=int((~np.isnan(M)).sum()), prevalence=round(np.nanmean(M), 3), P11=round(p11, 3), P01=round(p01, 3),
                N2=round(n2, 3), obs_minus_N2=round(p11 - n2, 3), jaccard1=round(jaccard_lag(M, 1)[0], 3), persistent=round(c.get('persistent', 0), 3),
                P11_with_other=round((b[withm] == 1).mean(), 3) if withm.any() else np.nan, P11_alone=round((b[alone] == 1).mean(), 3) if alone.any() else np.nan)

VARIANTS = [('base', None, None), ('DO < 4', 'DO < 4', None), ('pH 6-9', 'pH 6-9', None), ('FC > 500', 'FC > 500', None),
            ('values at the limit adverse', 'values at the limit adverse', None), ('empirical P80', 'empirical P80', None),
            ('empirical P90', 'empirical P90', None), ('FC coded as in rev 19', 'FC coded as in rev 19', None),
            ('min years 4', None, 'min years 4'), ('min years 5', None, 'min years 5'), ('min years 6', None, 'min years 6'), ('min years 9', None, 'min years 9'),
            ('harmonization: low/ambiguous excluded', None, 'harmonization: low/ambiguous excluded'),
            ('harmonization: all-high stations only', None, 'harmonization: all-high stations only'), ('raw CPCB codes', None, 'raw CPCB codes'),
            ('drop 2016', None, 'drop 2016'), ('drop 2017', None, 'drop 2017'), ('drop 2020', None, 'drop 2020'), ('drop 2016-2017', None, 'drop 2016-2017'),
            ('one-year gaps bridged', None, 'bridge')]
rows = []
for name, rc, pn in VARIANTS:
    D = panel(pn) if pn and pn != 'bridge' else d0[d0.in_panel_A_long].copy()
    if rc: D = recode(D, rc)
    for ind in IND:
        r = dict(variant=name, indicator=LABEL.get(ind, 'FC')); r.update(metrics(D, ind, bridge_gaps=(pn == 'bridge'))); rows.append(r)
    print(name, [(x['indicator'], x['P11'], x['obs_minus_N2']) for x in rows[-4:]], flush=True)
S = pd.DataFrame(rows); S.to_csv(f'{B}/results/S_robustness_full.csv', index=False)
P = S.pivot(index='variant', columns='indicator', values='P11').reindex([v[0] for v in VARIANTS]); print(P.to_string())
