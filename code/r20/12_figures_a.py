"""Rev 20 figures, part A: Fig 1 workflow, Fig 2 network and map, Fig 3 persistence vs nulls, Fig 4 overlap by lag,
Fig 5 threshold proximity, Fig S1 station classes."""
import sys
import numpy as np, pandas as pd, geopandas as gpd, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch
from figstyle import *
from core import load, wide, classify, B
R = f'{B}/results/'

# ---------------- Fig 1 workflow
def fig1():
    fig, ax = plt.subplots(figsize=(W2, 2.55)); ax.set_xlim(0, 100); ax.set_ylim(0, 40); ax.axis('off')
    boxes = [('CPCB annual\nreports', '9 reports\n2016-2024\n11,157 rows', 'tables extracted to one\nworkbook, never edited', '#f0efec'),
             ('Station\nharmonization', '1,891 stations\nfrom 1,840 codes', '100 code series switch\nlocation; 68 pairs\nreviewed by hand', '#cde2fb'),
             ('Value QC and\nstatus', 'bathing criteria\nDO, BOD, pH, FC', '155 values excluded;\nFC censored where it\ncannot be assessed', '#f0efec'),
             ('Persistence\nand nulls', 'P(1|1), overlap\nN1, N1s, N2\nmixed models', 'Panel A_long:\n1,724 stations', '#cde2fb'),
             ('Indicator\nreliability', 'analytical error\nsampling effort\nmonitoring regime', 'monthly records\n2020-2021,\n255 stations', '#f0efec'),
             ('Multi-year\nindicator', 'three classes\nout-of-sample\n2016-20 -> 2021-24', 'listing turnover;\ncompliance rules', '#cde2fb')]
    w, gap, x0 = 14.2, 2.6, 0.6
    for i, (h, body, foot, col) in enumerate(boxes):
        x = x0 + i * (w + gap)
        ax.add_patch(FancyBboxPatch((x, 13), w, 22, boxstyle='round,pad=0.25,rounding_size=1.2', fc=col, ec=AXIS, lw=0.9))
        ax.text(x + w / 2, 31.5, h, ha='center', va='top', fontsize=7.2, fontweight='bold', color=INK, linespacing=1.05)
        ax.text(x + w / 2, 21.8, body, ha='center', va='center', fontsize=6.2, color=INK2, linespacing=1.15)
        ax.text(x + w / 2, 11.2, foot, ha='center', va='top', fontsize=5.7, color=MUTED, linespacing=1.15)
        if i < len(boxes) - 1:
            ax.annotate('', xy=(x + w + gap - 0.4, 24), xytext=(x + w + 0.4, 24), arrowprops=dict(arrowstyle='-|>', color=MUTED, lw=1))
    ax.text(0.6, 38.8, 'Context for 1,187 located stations: IMD rainfall, WorldPop density, HydroBASINS, HydroRIVERS discharge; associations only',
            fontsize=6.2, color=INK2, style='italic')
    save(fig, 'Fig01_workflow')

# ---------------- Fig 2 network coverage and map
def fig2():
    cov = pd.read_csv(R + 'G4_coverage.csv', index_col=0); nobs = pd.read_csv(R + 'G4b_years_per_station.csv', index_col=0).iloc[:, 0]
    d = load('A_long'); S = pd.read_csv(f'{B}/analysis/stations_r20.csv')
    cls = {}
    for ind in ['ADV_BOD', 'ADV_FC_cens']:
        W = wide(d, ind)
        c = pd.Series([classify(r) for r in W.values.astype(float)], index=W.index)
        cls[ind] = c.map(lambda x: None if x is None else ('persistent' if x == 'persistent' else ('never adverse' if x == 'never' else 'changing')))
    states = gpd.read_file(f'{B}/data/geo/india_states_ne10m.gpkg'); nb = gpd.read_file(f'{B}/data/geo/neighbours_ne10m.gpkg')
    fig = plt.figure(figsize=(W2, 6.2))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.85], width_ratios=[1, 1, 0.78], hspace=0.38, wspace=0.32)
    ax = fig.add_subplot(gs[0, 0:2]); grid(ax, 'y')
    years = cov.index.values; bottom = np.zeros(len(years)); cohorts = [c for c in cov.columns]
    ramp = ['#0d366b', '#184f95', '#256abf', '#3987e5', '#5598e7', '#6da7ec', '#86b6ef', '#9ec5f4', '#b7d3f6']
    for k, c in enumerate(cohorts):
        v = cov[c].values; ax.bar(years, v, bottom=bottom, color=ramp[k], width=0.72, edgecolor='white', lw=0.6, label=c); bottom += v
    for x, t in zip(years, bottom): ax.text(x, t + 25, f'{int(t):,}', ha='center', fontsize=5.8, color=INK2)
    ax.set_ylabel('Stations reported'); ax.set_xticks(years); ax.set_ylim(0, 1780)
    ax.legend(title='year first reported', ncol=9, loc='upper center', fontsize=5.6, title_fontsize=5.8, handlelength=1, columnspacing=0.9,
              bbox_to_anchor=(0.5, -0.13))
    panel(ax, 'a', 'Stations reported per year, by year first reported')
    ax = fig.add_subplot(gs[0, 2]); grid(ax, 'y')
    ax.bar(nobs.index, nobs.values, color='#3987e5', width=0.7, edgecolor='white')
    ax.set_xticks(range(1, 10)); ax.set_xlabel('Years observed'); ax.set_ylabel('Stations'); panel(ax, 'b', 'Years per station')
    ceil = set(pd.read_csv(R + 'F1_ceiling_states.csv').query('ceiling_state == True').state)
    name_fix = {'JAMMU & KASHMIR': 'Jammu and Kashmir'}
    ceil_ne = {name_fix.get(s, s.title().replace(' And ', ' and ')) for s in ceil}
    colors = {'persistent': '#0d366b', 'changing': '#5598e7', 'never adverse': '#cde2fb'}
    for j, (ind, title) in enumerate([('ADV_BOD', 'BOD classes'), ('ADV_FC_cens', 'FC classes')]):
        ax = fig.add_subplot(gs[1, j])
        nb.plot(ax=ax, color='#f4f3ef', edgecolor='#d6d5ce', lw=0.3)
        st = states.copy(); st['ceil'] = st.name.isin(ceil_ne)
        st.plot(ax=ax, color=np.where(st.ceil & (ind == 'ADV_FC_cens'), '#e4e3dc', '#fcfcfb'), edgecolor='#a9a8a0', lw=0.35)
        if ind == 'ADV_FC_cens':
            st[st.ceil].plot(ax=ax, facecolor='none', edgecolor='#a9a8a0', lw=0.35, hatch='////')
        P = S[S.spatial_main.astype(str) == 'True'].set_index('station_uid').join(cls[ind].rename('cls')).dropna(subset=['cls'])
        for k in ['never adverse', 'changing', 'persistent']:
            q = P[P.cls == k]; ax.scatter(q.lon, q.lat, s=5.5, color=colors[k], edgecolor='white', lw=0.25, zorder=3,
                                          label=f'{k} ({len(q)})')
        ax.set_xlim(67.5, 97.8); ax.set_ylim(7.5, 37.3); ax.set_aspect('equal'); ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values(): s.set_visible(False)
        h, l = ax.get_legend_handles_labels()
        if ind == 'ADV_FC_cens':
            h.append(Patch(facecolor='#e4e3dc', edgecolor='#a9a8a0', hatch='////', lw=0.35)); l.append('not assessable')
        ax.legend(h, l, loc='upper center', bbox_to_anchor=(0.5, 0.0), ncol=2, fontsize=5.6, markerscale=1.8, handletextpad=0.3, columnspacing=0.5, labelspacing=0.3)
        ax.set_xlabel(''); ax.set_ylabel('')
        panel(ax, 'cd'[j], title)
    # located vs unlocated by state
    ax = fig.add_subplot(gs[1, 2])
    A = load('A'); s1 = A.groupby('station_uid').station_state.first()
    loc = S.set_index('station_uid').spatial_main.astype(str).eq('True').reindex(s1.index).fillna(False)
    t = pd.DataFrame({'state': s1, 'loc': loc}).groupby('state').loc.agg(['sum', 'size']); t = t[t['size'] >= 15].sort_values('size')
    y = np.arange(len(t))
    ax.barh(y, t['size'], color='#e4e3dc', height=0.72, label='not in spatial sample')
    ax.barh(y, t['sum'], color='#3987e5', height=0.72, label='main spatial sample')
    ax.set_yticks(y); ax.set_yticklabels([s.title().replace('&', 'and').replace('Pradesh', 'Pr.').replace('Jammu and Kashmir', 'J and K') for s in t.index], fontsize=5.3)
    ax.set_xlabel('Stations'); ax.legend(loc='lower right', fontsize=5.4); grid(ax, 'x'); panel(ax, 'e', 'Stations by state (>= 15)')
    fig.text(0.01, 0.005, 'Map lines delineate study areas and do not necessarily depict accepted national boundaries. Classes: stations with >= 4 observed years; '
             'changing = intermittent, emerging or diminishing.', fontsize=5.4, color=MUTED)
    save(fig, 'Fig02_network_map')

# ---------------- Fig 3 persistence vs nulls
def fig3():
    T = pd.read_csv(R + 'T3_core.csv').set_index('indicator')
    order = ['FC', 'BOD', 'DO', 'pH', '>=2 of 4']
    fig, ax = plt.subplots(figsize=(W2, 2.9)); grid(ax, 'x')
    for i, k in enumerate(order):
        r = T.loc[k]; y = len(order) - 1 - i; col = C[k]
        lo, hi = [float(v) for v in r.CI_station.split('-')]
        slo, shi = [float(v) for v in r.CI_state.split('-')]
        ax.plot([r.N1, r.P11], [y, y], color=GRID, lw=3, zorder=1, solid_capstyle='round')
        ax.plot([slo, shi], [y - 0.22, y - 0.22], color=col, lw=0.8, alpha=0.55, zorder=2)
        ax.scatter(r.N1, y, marker='|', s=80, color=MUTED, lw=1.5, zorder=3)
        ax.scatter(r.N1s, y, marker='|', s=80, color=INK, lw=1.5, zorder=3)
        ax.scatter(r.N2, y, marker='o', s=30, facecolor='white', edgecolor=col, lw=1.3, zorder=3)
        ax.errorbar(r.P11, y, xerr=[[r.P11 - lo], [hi - r.P11]], fmt=MK[k], ms=5.5, color=col, mec='white', mew=0.7, elinewidth=1.1, capsize=0, zorder=4)
        ax.text(min(shi, 0.99) + 0.012, y - 0.05, f'{r.P11:.2f}', fontsize=6.3, color=INK2, va='center')
    ax.set_yticks(range(len(order))); ax.set_yticklabels([LAB[k] for k in order][::-1]); ax.set_xlim(0, 1.04); ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel('Probability of being adverse next year, given adverse this year, P(1|1)')
    ax.legend(handles=[Line2D([], [], marker='o', ls='', color=INK2, mec='white', ms=5, label='observed (95% CI, station clusters)'),
                       Line2D([], [], color=INK2, lw=0.8, alpha=0.55, label='95% CI, state clusters'),
                       Line2D([], [], marker='o', ls='', mfc='white', mec=INK2, ms=5, label='N2: station propensity only'),
                       Line2D([], [], marker='|', ls='', color=MUTED, ms=8, mew=1.5, label='N1: reshuffled across India'),
                       Line2D([], [], marker='|', ls='', color=INK, ms=8, mew=1.5, label='N1s: reshuffled within state and year')],
              loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=3, fontsize=6)
    save(fig, 'Fig03_persistence_vs_nulls')

# ---------------- Fig 4 overlap by lag
def fig4():
    G = pd.read_csv(R + 'G1_overlap_lag.csv')
    fig, axs = plt.subplots(1, 4, figsize=(W2, 2.0), sharey=True)
    for ax, k in zip(axs, ['FC', 'BOD', 'DO', 'pH']):
        g = G[G.indicator == k]; grid(ax, 'y')
        ax.fill_between(g.lag, g.N2_lo, g.N2_hi, color=C[k], alpha=0.15, lw=0)
        ax.plot(g.lag, g.N2, color=C[k], lw=1, ls='--', label='N2 propensity')
        ax.plot(g.lag, g.N1, color=MUTED, lw=1, ls=':', label='N1 reshuffled')
        ax.plot(g.lag, g.observed, color=C[k], lw=1.6, marker=MK[k], ms=3.2, mec='white', mew=0.4, label='observed')
        ax.set_title(k, loc='left', fontsize=7.5, fontweight='bold', color=C[k]); ax.set_xticks([1, 2, 4, 6, 8]); ax.set_xlabel('Years apart')
        ax.set_ylim(0, 0.9)
    axs[0].set_ylabel('Jaccard overlap of adverse sets')
    axs[-1].legend(handles=[Line2D([], [], color=INK2, lw=1.6, marker='o', ms=3, label='observed'), Line2D([], [], color=INK2, lw=1, ls='--', label='N2 propensity (95% band)'),
                            Line2D([], [], color=MUTED, lw=1, ls=':', label='N1 reshuffled')], loc='upper right', fontsize=5.6)
    save(fig, 'Fig04_overlap_by_lag')

# ---------------- Fig 5 threshold proximity
def fig5():
    G = pd.read_csv(R + 'G2_threshold_proximity.csv')
    fig, axs = plt.subplots(1, 3, figsize=(W2, 2.2), sharey=True)
    for ax, k in zip(axs, ['BOD', 'DO', 'FC']):
        g = G[G.indicator == k]; grid(ax, 'y')
        ax.fill_between(g.x_median, g.lo, g.hi, color=C[k], alpha=0.18, lw=0, step=None)
        ax.plot(g.x_median, g.p_change, color=C[k], lw=1.4, marker=MK[k], ms=2.8, mec='white', mew=0.3)
        ax.axvline(g.criterion.iat[0], color=INK2, lw=0.8, ls='--')
        if k == 'BOD':
            ax.set_xscale('log'); ax.set_xticks([0.5, 1, 3, 10, 30]); ax.set_xticklabels(['0.5', '1', '3', '10', '30']); ax.minorticks_off()
        if k == 'FC':
            ax.set_xscale('log'); ax.set_xticks([10, 100, 1000, 1e4, 1e5]); ax.set_xticklabels(['10', '100', '1,000', '10,000', '100,000']); ax.minorticks_off()
        ax.set_xlabel({'BOD': 'BOD maximum in year t (mg/L)', 'DO': 'DO minimum in year t (mg/L)', 'FC': 'FC maximum in year t (MPN/100 mL)'}[k])
        ax.set_title(k, loc='left', fontsize=7.5, fontweight='bold', color=C[k])
        if k == 'DO': ax.set_xlim(0, 9)
    axs[0].set_ylabel('P(status changes next year)')
    save(fig, 'Fig05_threshold_proximity')

# ---------------- Fig S1 station classes
def figS1():
    G = pd.read_csv(R + 'G3_classes.csv').set_index('indicator')
    order = ['BOD', 'FC', 'DO', 'pH']; cols = {'persistent': '#0d366b', 'intermittent': '#3987e5', 'diminishing': '#86b6ef', 'emerging': '#eb6834', 'never': '#e4e3dc'}
    fig, ax = plt.subplots(figsize=(W2, 2.3))
    for i, k in enumerate(order):
        y = len(order) - 1 - i; left = 0
        for c in cols:
            v = G.loc[k, c]; ax.barh(y, v, left=left, color=cols[c], height=0.62, edgecolor='white', lw=1)
            if v > 0.06: ax.text(left + v / 2, y, f'{v:.0%}', ha='center', va='center', fontsize=6, color='white' if c in ('persistent', 'intermittent') else INK)
            left += v
        ax.text(1.02, y, f'persistent {G.loc[k, "persistent_min"]:.0%}-{G.loc[k, "persistent_max"]:.0%}\nacross 12 rules (n = {int(G.loc[k, "classified"])})',
                va='center', fontsize=5.7, color=MUTED)
    ax.set_yticks(range(len(order))); ax.set_yticklabels([LAB[k] for k in order][::-1]); ax.set_xlim(0, 1); ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v:.0%}'))
    ax.legend(handles=[Patch(color=v, label=k) for k, v in cols.items()], ncol=5, loc='upper center', bbox_to_anchor=(0.5, -0.18))
    save(fig, 'FigS1_station_classes')

if __name__ == '__main__':
    for f in sys.argv[1:] or ['fig1', 'fig2', 'fig3', 'fig4', 'fig5', 'figS1']:
        globals()[f](); print('done', f)
