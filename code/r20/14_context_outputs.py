"""Rev 20: Fig 9 (context), Table 6 and Supplementary context tables from results/C1_context_models.csv."""
import sys
import re
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from figstyle import *
from core import load, wide, classify, B
R = f'{B}/results/'; L = f'{B}/manuscript/generated/'
C = pd.read_csv(R + 'C1_context_models.csv')
IND = {'ADV_BOD': 'BOD', 'ADV_DO': 'DO', 'ADV_FC_cens': 'FC', 'ADV_PH': 'pH', 'ADV_MULTI_GE2_cens': '>=2 of 4'}
def g(ind, model, term=None):
    x = C[(C.indicator == ind) & (C.model == model)]
    if term: x = x[x.term == term]
    return x.iloc[0] if len(x) else None
M_W = 'C1 within-state, GLMM (rev 19 specification)'; M_F = 'C1 + sampling frequency 2023'
M_PA = 'C2 GLM state+year FE (population-averaged), unweighted'; M_IPW = 'C2 GLM state+year FE, inverse-probability-of-location weights'
M_PAF = 'C2 GLM state+year FE + sampling frequency 2023'; M_CL = 'C3 persistent class, conditional logit by state'
fmt = lambda r: f'{r.OR:.2f} ({r.lo:.2f}-{r.hi:.2f})' if r is not None else '--'

def fig9():
    d = load('A_long'); S = pd.read_csv(f'{B}/analysis/stations_r20.csv').set_index('station_uid')
    S = S[S.spatial_main.astype(str) == 'True']
    fig, axs = plt.subplots(1, 3, figsize=(W2, 2.7), gridspec_kw=dict(width_ratios=[0.9, 1.25, 1.05], wspace=0.55))
    ax = axs[0]; grid(ax, 'y')
    q = pd.qcut(S.lpop5, 5, labels=False)
    for k, lab in [('ADV_BOD', 'BOD'), ('ADV_FC_cens', 'FC'), ('ADV_DO', 'DO'), ('ADV_PH', 'pH')]:
        W = wide(d, k, 4); W = W[W.index.isin(S.index)]
        p = pd.Series([classify(r) == 'persistent' for r in W.values.astype(float)], index=W.index)
        y = p.groupby(q.reindex(W.index)).mean()
        ax.plot(y.index + 1, y.values, color=C_[lab], marker=MK[lab], ms=3.5, lw=1.3, mec='white', mew=0.4)
        ax.annotate(lab, (5, y.values[-1]), xytext=(4, 0), textcoords='offset points', color=C_[lab], fontsize=6.3, fontweight='bold', va='center')
    ax.set_xticks(range(1, 6)); ax.set_xlabel('Density fifth (5 km)'); ax.set_ylabel('Share of stations persistent')
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda v, _: f'{v:.0%}')); ax.set_xlim(0.7, 5.8)
    panel(ax, 'a', 'Persistent share')
    ax = axs[1]; grid(ax, 'x')
    specs = [(M_W, 'log10 density (within state)', 'conditional, within state', 'o', 'cond'),
             (None, None, 'population-averaged (approx.)', 'o', 'paapprox'),
             (M_PA, 'log10 density', 'population-averaged GLM', 's', 'pa'),
             (M_IPW, 'log10 density', 'GLM, location weights', 'D', 'ipw'),
             (M_PAF, 'log10 density', 'GLM + sampling frequency', '^', 'paf'),
             ('C6 upstream GLM state+year FE: lpop5 + lq', 'lpop5', 'GLM + river discharge', 'v', 'q')]
    inds = ['ADV_BOD', 'ADV_DO', 'ADV_FC_cens', 'ADV_PH']
    for i, k in enumerate(inds):
        lab = IND[k]; y0 = (len(inds) - 1 - i) * 1.0
        for j, (m, t, nm, mk, key) in enumerate(specs):
            y = y0 + 0.36 - j * 0.14
            if key == 'paapprox':
                r = g(k, M_W, 'log10 density (within state)'); v = r.PA_OR_approx if r is not None else np.nan
                ax.scatter(v, y, marker='o', s=16, facecolor='white', edgecolor=C_[lab], lw=1, zorder=3); continue
            r = g(k, m, t)
            if r is None: continue
            ax.plot([r.lo, r.hi], [y, y], color=C_[lab], lw=0.9)
            ax.scatter(r.OR, y, marker=mk, s=16, color=C_[lab], edgecolor='white', lw=0.4, zorder=3)
    ax.axvline(1, color=INK2, lw=0.7, ls='--'); ax.set_xscale('log'); ax.set_xticks([0.5, 1, 2, 4, 8]); ax.set_xticklabels(['0.5', '1', '2', '4', '8']); ax.minorticks_off()
    ax.set_yticks([3, 2, 1, 0]); ax.set_yticklabels(['BOD', 'DO', 'FC', 'pH']); ax.set_xlabel('Odds ratio per tenfold density')
    ax.legend(handles=[Line2D([], [], marker=s[3], ls='', color=INK2 if s[4] != 'paapprox' else 'white', mec=INK2, ms=4, label=s[2]) for s in specs],
              loc='upper center', bbox_to_anchor=(0.5, -0.24), ncol=2, fontsize=5.3, columnspacing=0.8, handletextpad=0.3)
    panel(ax, 'b', 'Within-state odds ratios')
    ax = axs[2]; grid(ax, 'x')
    sl = C[(C.indicator == 'ADV_BOD') & C.model.str.startswith('state slope')].copy()
    sl['state'] = sl.model.str.replace('state slope: ', '').str.title().str.replace('&', 'and').str.replace('Pradesh', 'Pr.').str.replace('Jammu And Kashmir', 'J and K')
    sl = sl[(sl.hi / sl.lo) < 500].sort_values('OR')
    for y, (_, r) in enumerate(sl.iterrows()):
        ax.plot([max(r.lo, 0.05), min(r.hi, 200)], [y, y], color=C_['BOD'], lw=0.9)
        ax.scatter(r.OR if r.OR < 200 else 200, y, marker='o', s=14, color=C_['BOD'], edgecolor='white', lw=0.4, zorder=3)
    ax.axvline(1, color=INK2, lw=0.7, ls='--'); ax.set_xscale('log'); ax.set_xlim(0.05, 250); ax.set_xticks([0.1, 1, 10, 100]); ax.set_xticklabels(['0.1', '1', '10', '100']); ax.minorticks_off()
    ax.set_yticks(range(len(sl))); ax.set_yticklabels(sl.state, fontsize=5.4); ax.set_xlabel('BOD odds ratio per tenfold density')
    panel(ax, 'c', 'By state, BOD')
    save(fig, 'Fig09_context')

def tables():
    rows = []
    for k in ['ADV_BOD', 'ADV_DO', 'ADV_FC_cens', 'ADV_PH', 'ADV_MULTI_GE2_cens']:
        w = g(k, M_W, 'log10 density (within state)'); LG = pd.read_csv(R + 'C2_loso_glm.csv'); lo_ = LG[(LG.indicator == k) & (LG.left_out != 'none')]
        q = g(k, 'C6 upstream GLM state+year FE: lpop5 + lq', 'lpop5'); ppq = g(k, 'C6 upstream GLM state+year FE: lppq', 'lppq')
        cl = g(k, M_CL)
        rows.append([{'ADV_BOD': 'BOD', 'ADV_DO': 'DO', 'ADV_FC_cens': 'FC', 'ADV_PH': 'pH', 'ADV_MULTI_GE2_cens': '$\\ge$ 2 of 4'}[k],
                     fmt(w), f'{w.PA_OR_approx:.2f}', f'{w.AME_pp:+.1f}', fmt(g(k, M_PA, 'log10 density')), fmt(g(k, M_IPW, 'log10 density')),
                     fmt(g(k, M_PAF, 'log10 density')), f'{lo_.OR.min():.2f}-{lo_.OR.max():.2f}' if len(lo_) else '--', fmt(q), fmt(ppq), fmt(cl)])
    out = [r'''\begin{table*}[width=\textwidth,cols=11,pos=h]
\caption{Population density and adverse status within states (main spatial sample; associations, not effects). Odds ratios (95\% CI) per tenfold population density within 5 km. Conditional: within-between GLMM with year effects and state, sub-basin and station random intercepts. PA approx.: population-averaged equivalent \citep{zeger1988}. AME: average marginal effect on the probability of adverse status, percentage points. GLM: logistic regression with state and year fixed effects, station-clustered errors, unweighted, with inverse-probability-of-location weights, and with the 2023 sampling frequency. LOSO: range of the GLM OR when each state with at least 20 stations is left out. Discharge: GLM adding log modelled long-term discharge of the snapped HydroRIVERS reach. Upstream/discharge: GLM with log upstream-catchment population per unit discharge in place of density (OR per tenfold). Persistent: conditional logit by state.}\label{tbl:context}
\begin{tabular*}{\tblwidth}{@{} p{0.045\textwidth} p{0.075\textwidth} p{0.04\textwidth} p{0.04\textwidth} p{0.075\textwidth} p{0.075\textwidth} p{0.075\textwidth} p{0.055\textwidth} p{0.075\textwidth} p{0.075\textwidth} p{0.075\textwidth} @{}}
\toprule
 & Conditional & PA approx. & AME (pp) & GLM & GLM, weights & GLM + frequency & LOSO & GLM + discharge & Upstream per discharge & Persistent \\
\midrule''']
    out += [' & '.join(r) + ' \\\\' for r in rows] + ['\\bottomrule', '\\end{tabular*}', '\\end{table*}']
    open(L + 'table_context.tex', 'w').write(re.sub(r'(?<=[ {])p\{', r'>{\\raggedright\\arraybackslash}p{', '\n'.join(out))); print('\n'.join(out))  # ragged-right text columns
    # SI: LOSO, state slopes, rainfall, variance, C6 all specs
    si = C.copy(); si['indicator'] = si.indicator.map(IND)
    si.to_csv(R + 'SI_context_all_models.csv', index=False)

C_ = {'BOD': '#2a78d6', 'FC': '#eb6834', 'DO': '#1baf7a', 'pH': '#4a3aa7', '>=2 of 4': '#52514e'}
if __name__ == '__main__':
    fig9(); tables()
