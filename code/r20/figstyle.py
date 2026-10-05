"""Shared figure style for rev 20 (validated palette: BOD blue, FC orange, DO aqua, pH violet; all-pairs CVD and
normal-vision floors pass in light mode; aqua is below 3:1 contrast, so every series also carries a direct label or
a distinct marker)."""
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
mpl.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7, 'axes.titlesize': 7.5, 'axes.labelsize': 7,
                     'xtick.labelsize': 6.5, 'ytick.labelsize': 6.5, 'legend.fontsize': 6.3, 'axes.edgecolor': '#bdbcb4',
                     'axes.linewidth': 0.8, 'xtick.color': '#52514e', 'ytick.color': '#52514e', 'axes.labelcolor': '#52514e',
                     'xtick.major.size': 2.5, 'ytick.major.size': 2.5, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
                     'axes.spines.top': False, 'axes.spines.right': False, 'legend.frameon': False,
                     'savefig.dpi': 400, 'pdf.fonttype': 42, 'figure.dpi': 150})
C = {'BOD': '#2a78d6', 'FC': '#eb6834', 'DO': '#1baf7a', 'pH': '#4a3aa7', 'PH': '#4a3aa7', '>=2 of 4': '#52514e'}
MK = {'BOD': 'o', 'FC': 's', 'DO': '^', 'pH': 'D', 'PH': 'D', '>=2 of 4': 'o'}
LAB = {'BOD': 'BOD > 3 mg/L', 'FC': 'FC > 2,500 MPN/100 mL', 'DO': 'DO < 5 mg/L', 'pH': 'pH outside 6.5-8.5', 'PH': 'pH outside 6.5-8.5',
       '>=2 of 4': '>= 2 of 4 primary adverse'}
INK, INK2, MUTED, GRID, AXIS = '#0b0b0b', '#52514e', '#8a8983', '#e4e3dc', '#bdbcb4'
SEQ = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5', '#256abf', '#184f95', '#0d366b']
W2, W1 = 7.2, 3.5     # double / single column width (inches)
import os
OUT = os.path.join(os.environ.get('WQ_ROOT', os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))), 'figures') + '/'
os.makedirs(OUT, exist_ok=True)
def grid(ax, axis='x'):
    ax.grid(axis=axis, color=GRID, lw=0.6, zorder=0); ax.set_axisbelow(True)
def panel(ax, letter, title):
    ax.set_title(f'{letter}  {title}', loc='left', fontsize=7.5, fontweight='bold', color=INK, pad=6)
def save(fig, name):
    fig.savefig(OUT + name + '.pdf', bbox_inches='tight'); fig.savefig(OUT + name + '.png', bbox_inches='tight', dpi=300); plt.close(fig)
