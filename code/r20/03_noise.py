"""Rev 20, step 3 (Major Comment 4 v, Minor 11): analytical-error benchmark as expected versus observed counts.
H0 for a pair of consecutive observed years: the true value is the same in both years (estimated by the pair mean on
the error scale) and each reported value adds independent analytical error. Under H0 the pair changes status with
probability 2p(1-p), p = P(measured value beyond the criterion). The sum over all pairs is the number of status changes
analytical error alone would produce if nothing changed; its 95% range comes from simulating the Poisson-binomial sum.
Also: expected share of station-years whose status a repeat measurement of the same samples would change
(an indicator-quality metric), and the rev 19 'indistinguishable' share for continuity (SI only).
Error models: BOD multiplicative CV 10 / 15.4 / 25% (Standard Methods 5210 B GGA limits); pH SD 0.05 / 0.07 / 0.10;
FC SD of log10 MPN: 0.23 (Haldane 1939 expected information, 5 tubes x 10/1/0.1 mL, central range; computed here),
0.26 and 0.33 (Cochran 1950 approximation 0.58/sqrt(n) for tenfold series, n = 5 and 3 tubes). The tube configuration
of each laboratory is unknown, so the range is reported with the main value. DO is not benchmarked: the annual DO
minimum is dominated by diel and field variability, which laboratory Winkler precision does not describe.
Output: results/E1_mpn_precision.csv, E2_error_benchmark.csv"""
import sys
from core import *
from scipy.stats import norm
rng = np.random.default_rng(11)
d = load('A_long')

def se_log10_mpn(lam, vols=(10, 1, 0.1), tubes=5):
    m = np.array(vols); T = np.exp(-lam * m)
    info = (tubes * m ** 2 * T / (1 - T)).sum()
    return 1 / (2.303 * lam * np.sqrt(info))
rows = []
for tubes in (5, 3):
    lam = np.exp(np.linspace(np.log(0.05), np.log(2.0), 200))
    se = np.array([se_log10_mpn(l, tubes=tubes) for l in lam])
    rows.append(dict(tubes=tubes, haldane_expected_info_median=round(np.median(se), 3), haldane_min=round(se.min(), 3), haldane_max=round(se.max(), 3),
                     cochran_0p58_over_sqrt_n=round(0.58 / np.sqrt(tubes), 3)))
E1 = pd.DataFrame(rows); E1.to_csv(f'{B}/results/E1_mpn_precision.csv', index=False); print(E1.to_string(index=False))

SC = {'BOD': ('log', 3.0, [np.log(1.10), np.log(1.154), np.log(1.25)], 'ADV_BOD', ['BOD_max']),
      'FC': ('log10', 2500, [0.23, 0.26, 0.33], 'ADV_FC', ['FC_max']),
      'FC (censored)': ('log10', 2500, [0.23, 0.26, 0.33], 'ADV_FC_cens', ['FC_max']),
      'pH': ('lin', None, [0.05, 0.07, 0.10], 'ADV_PH', ['PH_min', 'PH_max'])}
tr = {'log': np.log, 'lin': lambda x: x, 'log10': np.log10}
rows = []
for ind, (sc, T, sds, adv, cols) in SC.items():
    S = wide(d, adv, 1).values.astype(float); idx = wide(d, adv, 1).index
    V = {c: d.pivot_table(index='station_uid', columns='Year', values=c, aggfunc='first').reindex(index=idx, columns=YEARS).values.astype(float) for c in cols}
    a, b = S[:, :-1], S[:, 1:]; ok = ~np.isnan(a) & ~np.isnan(b); obs = int((ok & (a != b)).sum()); npairs = int(ok.sum())
    for lab, sd in zip(['low', 'main', 'high'], sds):
        if ind == 'pH':
            lo, hi = V['PH_min'], V['PH_max']
            p_y = 1 - (1 - norm.cdf((6.5 - lo) / sd)) * (1 - norm.cdf((hi - 8.5) / sd))       # P(measured adverse | reported)
            mlo = (lo[:, :-1] + lo[:, 1:]) / 2; mhi = (hi[:, :-1] + hi[:, 1:]) / 2
            p_pair = 1 - (1 - norm.cdf((6.5 - mlo) / sd)) * (1 - norm.cdf((mhi - 8.5) / sd))
        else:
            X = tr[sc](np.clip(V[cols[0]], 1e-9, None)); t = tr[sc](T)
            p_y = norm.cdf((X - t) / sd)
            p_pair = norm.cdf(((X[:, :-1] + X[:, 1:]) / 2 - t) / sd)
        pc = np.where(ok, 2 * p_pair * (1 - p_pair), 0.0)
        pc = np.nan_to_num(pc)
        sims = (rng.random((2000, pc.size)) < pc.ravel()).sum(1)
        remeas = np.nansum(np.where(~np.isnan(S), np.where(S == 1, 1 - p_y, p_y), np.nan)) / np.sum(~np.isnan(S))
        # rev 19 share indistinguishable (kept for the SI)
        if ind == 'pH':
            near_lo = np.minimum(np.abs(lo[:, :-1] - 6.5), np.abs(lo[:, 1:] - 6.5)) < np.minimum(np.abs(hi[:, :-1] - 8.5), np.abs(hi[:, 1:] - 8.5))
            diff = np.where(near_lo, np.abs(lo[:, 1:] - lo[:, :-1]), np.abs(hi[:, 1:] - hi[:, :-1])) / sd
        else:
            diff = np.abs(X[:, 1:] - X[:, :-1]) / sd
        flip = ok & (a != b)
        rows.append(dict(indicator=ind, scenario=lab, error_sd=round(float(sd), 3), pairs=npairs, observed_changes=obs,
                         expected_changes_error_only=round(pc.sum(), 1), lo95=int(np.percentile(sims, 2.5)), hi95=int(np.percentile(sims, 97.5)),
                         expected_over_observed=round(pc.sum() / obs, 3),
                         repeat_measurement_changes_station_years=round(float(remeas), 4),
                         share_indistinguishable_rev19=round((flip & (diff < 1.96 * np.sqrt(2))).sum() / obs, 3)))
        print(rows[-1])
E2 = pd.DataFrame(rows); E2.to_csv(f'{B}/results/E2_error_benchmark.csv', index=False); print(E2.to_string(index=False))
