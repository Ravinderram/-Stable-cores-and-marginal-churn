"""Rev 20 shared functions: wide station x year matrices, transitions, nulls, overlap, classes, cluster bootstrap."""
import pandas as pd, numpy as np
import os
B = os.environ.get('WQ_ROOT', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..')))
YEARS = list(range(2016, 2025))
PRIMARY = ['ADV_FC', 'ADV_BOD', 'ADV_DO', 'ADV_PH']
LABEL = {'ADV_FC_cens': 'FC', 'ADV_MULTI_GE2_cens': '>=2 of 4', 'ADV_FC': 'FC (rev 19 coding)', 'ADV_BOD': 'BOD', 'ADV_DO': 'DO', 'ADV_PH': 'pH', 'ADV_MULTI_GE2': '>=2 of 4'}

def load(panel='A_long'):
    d = pd.read_csv(f'{B}/analysis/analysis_r20.csv', low_memory=False)
    return d[d[f'in_panel_{panel}']].copy() if panel else d

def wide(D, col, min_obs=2):
    W = D.pivot_table(index='station_uid', columns='Year', values=col, aggfunc='first').reindex(columns=YEARS)
    W = W[W.notna().sum(1) >= min_obs]
    return W

def counts(M):
    """per-station transition counts between consecutive calendar years (both observed)"""
    a, b = M[:, :-1], M[:, 1:]
    ok = ~np.isnan(a) & ~np.isnan(b)
    n11 = ((a == 1) & (b == 1) & ok).sum(1); n10 = ((a == 1) & (b == 0) & ok).sum(1)
    n01 = ((a == 0) & (b == 1) & ok).sum(1); n00 = ((a == 0) & (b == 0) & ok).sum(1)
    return n10, n11, n01, n00

def stay_entry(M):
    n10, n11, n01, n00 = counts(M)
    return n11.sum() / max((n10 + n11).sum(), 1), n01.sum() / max((n01 + n00).sum(), 1)

def null_N1(M, rng):
    X = M.copy()
    for j in range(X.shape[1]):
        o = np.where(~np.isnan(X[:, j]))[0]; X[o, j] = rng.permutation(X[o, j])
    return X

def null_N1s(M, groups, rng):
    X = M.copy()
    for g in np.unique(groups):
        r = np.where(groups == g)[0]
        for j in range(X.shape[1]):
            o = r[~np.isnan(X[r, j])]
            if len(o) > 1: X[o, j] = rng.permutation(X[o, j])
    return X

def null_N2(M, rng):
    X = M.copy()
    for i in range(X.shape[0]):
        o = np.where(~np.isnan(X[i]))[0]; X[i, o] = rng.permutation(X[i, o])
    return X

def jaccard_lag(M, lag):
    inter = union = 0; vals = []
    for t in range(M.shape[1] - lag):
        a, b = M[:, t], M[:, t + lag]; ok = ~np.isnan(a) & ~np.isnan(b)
        i = ((a == 1) & (b == 1) & ok).sum(); u = (((a == 1) | (b == 1)) & ok).sum()
        if u: vals.append(i / u)
        inter += i; union += u
    return (np.mean(vals) if vals else np.nan), (inter / union if union else np.nan)

def longest_run(row):
    best = cur = 0
    for v in row:
        if v == 1: cur += 1; best = max(best, cur)
        elif v == 0: cur = 0
        # missing years neither extend nor break a run of observed years
    return best

def classify(row, rec=0.6, run=3, min_obs=4):
    v = row[~np.isnan(row)]
    if len(v) < min_obs: return None
    if v.sum() == 0: return 'never'
    if v.mean() >= rec and longest_run(row) >= run: return 'persistent'
    h = len(v) // 2
    first, second = v[:h].mean(), v[len(v) - h:].mean()
    if second - first >= 0.5: return 'emerging'
    if first - second >= 0.5: return 'diminishing'
    return 'intermittent'

def class_shares(M, **kw):
    c = pd.Series([classify(r, **kw) for r in M]).dropna()
    return c.value_counts(normalize=True), len(c)

def cluster_boot(M, clusters, fn, R=1000, rng=None):
    """bootstrap over clusters (rows sharing a cluster label are resampled together)"""
    rng = rng or np.random.default_rng(1)
    labs = np.unique(clusters); idx = {l: np.where(clusters == l)[0] for l in labs}
    out = []
    for _ in range(R):
        pick = rng.choice(labs, len(labs), replace=True)
        rows = np.concatenate([idx[l] for l in pick])
        out.append(fn(M[rows]))
    out = np.array(out, dtype=float)
    return np.nanpercentile(out, [2.5, 97.5], axis=0)
