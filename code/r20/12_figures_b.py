"""Rev 20 figures, part B: Fig 6 noise and sampling, Fig 7 multi-year indicator, Fig 8 microbial vs organic,
Fig 10 robustness."""
import sys
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import Patch
from figstyle import *
from core import load, wide, classify, B
R = f'{B}/results/'
pct = plt.FuncFormatter(lambda v, _: f'{v:.0%}')

def fig6():
    E = pd.read_csv(R + 'E2_error_benchmark.csv'); M6 = pd.read_csv(R + 'M6_null_flips.csv'); M5 = pd.read_csv(R + 'M5_subsampling.csv')
    fig, axs = plt.subplots(1, 3, figsize=(W2, 2.55), gridspec_kw=dict(width_ratios=[1.15, 1, 1], wspace=0.42))
    ax = axs[0]; grid(ax, 'x')
    for i, (k, lab) in enumerate([('BOD', 'BOD'), ('FC (censored)', 'FC'), ('pH', 'pH')]):
        g = E[E.indicator == k].set_index('scenario'); y = 2 - i; n = g.pairs.iat[0]
        obs = 100 * g.observed_changes.iat[0] / n
        lo, mn, hi = (100 * g.loc[s, 'expected_changes_error_only'] / n for s in ('low', 'main', 'high'))
        ax.plot([lo, hi], [y - 0.12] * 2, color=C[lab], lw=1.2); ax.scatter(mn, y - 0.12, marker=MK[lab], s=28, facecolor='white', edgecolor=C[lab], lw=1.2, zorder=3)
        ax.scatter(obs, y + 0.12, marker=MK[lab], s=30, color=C[lab], edgecolor='white', lw=0.5, zorder=3)
        ax.text(max(hi, obs) + 0.6, y, f'{int(g.observed_changes.iat[0]):,} observed\n{mn * n / 100:,.0f} expected', fontsize=5.4, color=INK2, va='center')
    ax.set_yticks([2, 1, 0]); ax.set_yticklabels(['BOD', 'FC', 'pH']); ax.set_xlim(0, 32); ax.set_ylim(-0.5, 2.5)
    ax.set_xlabel('Status changes per 100 year pairs')
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=INK2, ms=4, label='observed'),
                       Line2D([], [], marker='o', ls='-', mfc='white', mec=INK2, color=INK2, ms=4, lw=1, label='expected from analytical\nerror alone (range: error models)')],
              loc='upper center', bbox_to_anchor=(0.5, -0.25), fontsize=5.4, ncol=2)
    panel(ax, 'a', 'Analytical error, all stations')
    ax = axs[1]; grid(ax, 'y')
    rows = [('BOD', 'main: table-linked, all months', 'BOD\nall'), ('BOD', 'COVID-19: June-December only, both years', 'BOD\nJun-Dec'),
            ('PH', 'main: table-linked, all months', 'pH\nall'), ('PH', 'COVID-19: June-December only, both years', 'pH\nJun-Dec')]
    for x, (k, v, lab) in enumerate(rows):
        r = M6[(M6.indicator == k) & (M6.variant == v)].iloc[0]; col = C['pH' if k == 'PH' else k]
        ax.bar(x - 0.19, r.observed_changes, width=0.36, color=col, edgecolor='white')
        ax.bar(x + 0.19, r.pooled_expected, width=0.36, color='white', edgecolor=col, lw=1.0, hatch='//////')
        ax.errorbar(x + 0.19, r.pooled_expected, yerr=[[r.pooled_expected - r.pooled_lo], [r.pooled_hi - r.pooled_expected]], fmt='none', ecolor=INK2, elinewidth=0.8, capsize=2)
    ns = [int(M6[(M6.indicator == k) & (M6.variant == v)].stations.iat[0]) for k, v, _ in rows]
    ax.set_xticks(range(len(rows))); ax.set_xticklabels([f'{r[2]}\nn={n}' for r, n in zip(rows, ns)], fontsize=5.6); ax.set_ylim(0, 31); ax.set_ylabel('Status changes, 2020 to 2021')
    ax.legend(handles=[Patch(color=INK2, label='observed'), Patch(facecolor='white', edgecolor=INK2, hatch='//////', label='expected if nothing changed (95%)')],
              loc='upper left', fontsize=5.2)
    panel(ax, 'b', 'Sampled months, 2020-2021')
    ax = axs[2]; grid(ax, 'y')
    for k, lab in [('BOD', 'BOD'), ('FC', 'FC'), ('PH', 'pH')]:
        g = M5[(M5.indicator == k) & (M5.samples_used != 'quarterly')].copy(); g['k'] = g.samples_used.astype(int)
        xs = list(g.k) + [12]; ys = list(g.adverse_share_k / g.adverse_share_all) + [1.0]
        ax.plot(xs, ys, color=C[lab], marker=MK[lab], ms=3, lw=1.2, mec='white', mew=0.3)
        ax.text(12.3, ys[-2] if lab != 'FC' else ys[-2] - 0.03, '', fontsize=6)
        ax.annotate(lab, (xs[2], ys[2]), xytext=(-4, {'BOD': 6, 'FC': -12, 'pH': 6}[lab]), textcoords='offset points', ha='right', color=C[lab], fontsize=6.3, fontweight='bold')
    ax.axvline(4, color=AXIS, lw=0.8, ls='--'); ax.text(4.25, 0.24, 'baseline\nprotocol', fontsize=5.4, color=MUTED)
    ax.set_xticks([1, 2, 3, 4, 6, 8, 12]); ax.set_xticklabels(['1', '2', '3', '4', '6', '8', 'all']); ax.set_ylim(0.15, 1.05)
    ax.yaxis.set_major_formatter(pct); ax.set_xlabel('Samples per year used'); ax.set_ylabel('Adverse share vs all samples')
    panel(ax, 'c', 'Samples per year')
    save(fig, 'Fig06_noise_sampling')

def fig7():
    I1 = pd.read_csv(R + 'I1_tci_out_of_sample.csv'); I3 = pd.read_csv(R + 'I3_listing_churn.csv'); I4 = pd.read_csv(R + 'I4_rules_monthly.csv')
    fig, axs = plt.subplots(1, 3, figsize=(W2, 2.5), gridspec_kw=dict(width_ratios=[1.1, 1, 1.05], wspace=0.45))
    ax = axs[0]; grid(ax, 'y'); cls = ['reliably compliant', 'borderline', 'reliably adverse']
    for j, k in enumerate(['BOD', 'FC', 'DO', 'PH']):
        g = I1[I1.indicator == k].set_index('class_2016_2020').reindex(cls); lab = 'pH' if k == 'PH' else k
        xs = np.arange(3) + (j - 1.5) * 0.13
        ax.plot(xs, g.adverse_share_2021_2024, color=C[lab], lw=0.9, alpha=0.6)
        ax.scatter(xs, g.adverse_share_2021_2024, marker=MK[lab], s=22, color=C[lab], edgecolor='white', lw=0.4, zorder=3, label=lab)
    ax.set_xticks(range(3)); ax.set_xticklabels(['reliably\ncompliant', 'borderline', 'reliably\nadverse'], fontsize=6)
    ax.yaxis.set_major_formatter(pct); ax.set_ylim(0, 1); ax.set_ylabel('Adverse share of 2021-2024 station-years')
    ax.set_xlabel('Class from 2016-2020'); ax.legend(loc='upper left', fontsize=5.6, ncol=2, columnspacing=0.6, handletextpad=0.2)
    panel(ax, 'a', 'Out-of-sample test')
    ax = axs[1]; grid(ax, 'y')
    for x, k in enumerate(['BOD', 'FC', 'DO', 'PH']):
        r = I3[I3.indicator == k].iloc[0]; lab = 'pH' if k == 'PH' else k
        ax.bar(x - 0.19, r.single_year_turnover, width=0.36, color='white', edgecolor=C[lab], hatch='//////', lw=1)
        ax.bar(x + 0.19, r.tci_turnover, width=0.36, color=C[lab], edgecolor='white')
        ax.text(x - 0.19, r.single_year_turnover + 0.02, f'{r.single_year_turnover:.2f}', ha='center', fontsize=5.3, color=INK2)
        ax.text(x + 0.19, r.tci_turnover + 0.02, f'{r.tci_turnover:.2f}', ha='center', fontsize=5.3, color=INK2)
    ax.set_xticks(range(4)); ax.set_xticklabels(['BOD', 'FC', 'DO', 'pH']); ax.set_ylabel('Entries + exits per listed station per year')
    ax.legend(handles=[Patch(facecolor='white', edgecolor=INK2, hatch='//////', label='single-year list'), Patch(color=INK2, label='multi-year: reliably adverse')],
              loc='upper left', fontsize=5.4)
    panel(ax, 'b', 'List turnover, 2018-2024')
    ax = axs[2]; grid(ax, 'y'); g = I4[I4.indicator == 'BOD']
    names = {'single sample (current)': 'single\nsample', 'raw score > 10%': '> 10% of\nsamples', 'binomial, 10% exceedance, alpha 0.05': 'binomial\n(10%)', 'median (> 50%)': '> 50% of\nsamples'}
    for x, (_, r) in enumerate(g.iterrows()):
        ax.bar(x - 0.19, r.changes_2020_2021, width=0.36, color=C['BOD'], edgecolor='white')
        ax.bar(x + 0.19, r.expected_changes_no_change_null, width=0.36, color='white', edgecolor=C['BOD'], hatch='//////', lw=1)
        ax.errorbar(x + 0.19, r.expected_changes_no_change_null, yerr=[[r.expected_changes_no_change_null - r.null_lo], [r.null_hi - r.expected_changes_no_change_null]],
                    fmt='none', ecolor=INK2, elinewidth=0.8, capsize=2)
    ax.set_xticks(range(len(g))); ax.set_xticklabels([f'{names[r.rule]}\n{r.adverse_share:.0%} adv.' for _, r in g.iterrows()], fontsize=5.6); ax.set_ylim(0, 17)
    ax.set_ylabel('BOD status changes, 2020 to 2021')
    ax.legend(handles=[Patch(color=C['BOD'], label='observed'), Patch(facecolor='white', edgecolor=C['BOD'], hatch='//////', label='expected if nothing changed')],
              loc='upper right', fontsize=5.3)
    panel(ax, 'c', 'Compliance rules, 73 stations')
    save(fig, 'Fig07_multiyear_indicator')

def fig8():
    d = load('A_long'); F = pd.read_csv(R + 'F3_typology.csv')
    Wf = wide(d, 'ADV_FC_cens', 4); Wo = wide(d, 'ADV_ORG', 4); idx = Wf.index.intersection(Wo.index)
    rf = Wf.loc[idx].mean(1); ro = Wo.loc[idx].mean(1)
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.7), gridspec_kw=dict(width_ratios=[1, 1.15], wspace=0.35))
    ax = axs[0]; rng = np.random.default_rng(1)
    ax.scatter(ro + rng.uniform(-0.015, 0.015, len(ro)), rf + rng.uniform(-0.015, 0.015, len(rf)), s=5, color='#3987e5', alpha=0.45, lw=0)
    ax.axhline(0.5, color=AXIS, lw=0.7, ls='--'); ax.axvline(0.5, color=AXIS, lw=0.7, ls='--')
    for (x, y, t) in [(0.75, 0.97, 'multi-domain'), (0.18, 0.97, 'mostly microbial'), (0.75, 0.03, 'mostly organic'), (0.18, 0.03, 'neither')]:
        ax.text(x, y, t, ha='center', va='center', fontsize=5.8, color=INK2)
    ax.set_xlabel('Share of years with DO < 5 or BOD > 3'); ax.set_ylabel('Share of years with FC > 2,500')
    ax.set_xlim(-0.03, 1.03); ax.set_ylim(-0.03, 1.03); panel(ax, 'a', f'Station recurrence (n = {len(idx)})')
    ax = axs[1]; cats = ['multi_domain', 'predominantly_microbial', 'predominantly_organic', 'neither']
    labs = ['multi-domain', 'mostly microbial', 'mostly organic', 'neither']; cols = ['#0d366b', '#eb6834', '#3987e5', '#e4e3dc']
    for i, (_, r) in enumerate(F.iterrows()):
        left = 0; y = 1 - i
        for c, l, col in zip(cats, labs, cols):
            ax.barh(y, r[c], left=left, color=col, height=0.55, edgecolor='white', lw=1)
            if r[c] > 0.07: ax.text(left + r[c] / 2, y, f'{r[c]:.0%}', ha='center', va='center', fontsize=6, color='white' if col in ('#0d366b', '#eb6834', '#3987e5') else INK)
            left += r[c]
    ax.set_yticks([1, 0]); ax.set_yticklabels(['FC ceiling coded\ncompliant (rev 19)', 'FC censored\n(rev 20, primary)'], fontsize=6)
    ax.xaxis.set_major_formatter(pct); ax.set_xlim(0, 1)
    ax.legend(handles=[Patch(color=c, label=l) for c, l in zip(cols, labs)], ncol=4, loc='upper center', bbox_to_anchor=(0.45, -0.14), fontsize=5.6, columnspacing=0.8)
    for i, (_, r) in enumerate(F.iterrows()):
        ax.text(1.0, 1 - i + 0.36, f'n = {int(r.stations)}, Spearman {r.spearman:.2f}', ha='right', fontsize=5.6, color=MUTED)
    panel(ax, 'b', 'Station typology (50% cut-off)')
    save(fig, 'Fig08_microbial_vs_organic')

def fig10():
    S = pd.read_csv(R + 'S_robustness_full.csv'); groups = {
        'thresholds': ['DO < 4', 'pH 6-9', 'FC > 500', 'values at the limit adverse', 'empirical P80', 'empirical P90', 'FC coded as in rev 19'],
        'panel': ['min years 4', 'min years 5', 'min years 6', 'min years 9'],
        'identity': ['harmonization: low/ambiguous excluded', 'harmonization: all-high stations only', 'raw CPCB codes'],
        'years': ['drop 2016', 'drop 2017', 'drop 2020', 'drop 2016-2017', 'one-year gaps bridged']}
    gcol = {'thresholds': '#0b0b0b', 'panel': '#52514e', 'identity': '#8a8983', 'years': '#256abf'}
    gmk = {'thresholds': 'D', 'panel': 'o', 'identity': '^', 'years': 's'}
    fig, axs = plt.subplots(1, 2, figsize=(W2, 2.4), sharey=True, gridspec_kw=dict(wspace=0.16))
    order = ['FC', 'BOD', 'DO', 'pH']
    for ax, m, xl in [(axs[0], 'P11', 'P(1|1)'), (axs[1], 'obs_minus_N2', 'Observed P(1|1) minus N2 (ordering of adverse years)')]:
        grid(ax, 'x')
        for i, k in enumerate(order):
            y = 3 - i; base = S[(S.variant == 'base') & (S.indicator == k)][m].iat[0]
            for gi, (g, vs) in enumerate(groups.items()):
                v = S[(S.indicator == k) & S.variant.isin(vs)][m]
                ax.scatter(v, np.full(len(v), y + (gi - 1.5) * 0.11), marker=gmk[g], s=14, facecolor='white', edgecolor=gcol[g], lw=0.9, zorder=3,
                           label=g if i == 0 else None)
            ax.scatter(base, y, marker='|', s=200, color='#eb6834', lw=2.0, zorder=4, label='base' if i == 0 else None)
        ax.set_xlabel(xl)
    axs[0].set_yticks(range(4)); axs[0].set_yticklabels(order[::-1]); axs[0].set_xlim(0, 1); axs[1].set_xlim(0, 0.17)
    axs[0].legend(loc='upper center', bbox_to_anchor=(1.0, -0.2), ncol=5, fontsize=5.8)
    panel(axs[0], 'a', 'Persistence across 19 variants'); panel(axs[1], 'b', 'Beyond station propensity')
    save(fig, 'Fig10_robustness')

if __name__ == '__main__':
    for f in sys.argv[1:] or ['fig6', 'fig7', 'fig8', 'fig10']:
        globals()[f](); print('done', f)
