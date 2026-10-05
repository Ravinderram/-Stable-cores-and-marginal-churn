"""Rev 20, step 5b (Major Comment 3b): does between-station structure survive standardized sampling effort?
Monthly subset (stations in the CPCB table in 2020 and 2021 with monthly records both years). Status from all samples
versus from one sample per quarter (50 random draws). With two years per station the latent GLMM ICC is unstable, so
the between-station share is measured by the one-way ANOVA ICC for binary data (k = 2 years) and by P(1|1) and P(1|0).
Output: results/R5_fixed_effort.csv"""
import pandas as pd, numpy as np
import os
O = os.path.join(os.environ.get('WQ_ROOT', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))), 'results') + '/'
Q = pd.read_csv(O + 'M8_fixed_effort_status.csv')
def anova_icc(w):
    k = 2; m = w.mean(1); g = w.values.mean()
    msb = k * ((m - g) ** 2).sum() / (len(w) - 1); msw = ((w.sub(m, axis=0)) ** 2).values.sum() / (len(w) * (k - 1))
    return (msb - msw) / (msb + (k - 1) * msw) if (msb + msw) > 0 else np.nan
def stats(w):
    a, b = w[2020].values, w[2021].values
    return anova_icc(w), (b[a == 1] == 1).mean() if (a == 1).any() else np.nan, (b[a == 0] == 1).mean(), w.values.mean()
rows = []
for ind, g in Q.groupby('indicator'):
    W = g.pivot(index='station_uid', columns='year', values='status_all').dropna()
    full = stats(W)
    dr = np.array([stats(g.pivot(index='station_uid', columns='year', values=f'q{k}').reindex(W.index)) for k in range(50)])
    rows.append(dict(indicator=ind, stations=len(W), prevalence_all=round(full[3], 3), prevalence_quarterly=round(np.nanmean(dr[:, 3]), 3),
                     ICC_all=round(full[0], 3), ICC_quarterly_median=round(np.nanmedian(dr[:, 0]), 3),
                     ICC_quarterly_range=f'{np.nanpercentile(dr[:, 0], 5):.2f}-{np.nanpercentile(dr[:, 0], 95):.2f}',
                     P11_all=round(full[1], 3), P11_quarterly=round(np.nanmedian(dr[:, 1]), 3), P01_all=round(full[2], 3), P01_quarterly=round(np.nanmedian(dr[:, 2]), 3)))
R5 = pd.DataFrame(rows); R5.to_csv(O + 'R5_fixed_effort.csv', index=False); print(R5.to_string(index=False))
