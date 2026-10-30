#!/usr/bin/env python3
"""
make_bib_figs.py -- redraw all figures to the Briefings in Bioinformatics submission specs.

Specification highlights (BIB author instructions):
  - image resolution >= 300 dpi (bitmaps) / vector preferred
  - avoid relying on red-green alone (about 8% of male readers have red-green
    colour vision deficiency); use a triple encoding of "hue + lightness + shape"
  - font size no smaller than 6 pt at final layout (single column ~90 mm wide)
  - no pie charts (fewer than 4 or more than 12 slices); no bar charts with fewer
    than 5 observations
  - alt text must be provided

Output:
  Fig1_audit_landscape     four panels: T4 / T1 / T2 / T4r
  Fig2_patient_replication replication across 8 TAA patients
  Fig3_ssr_detection       detection rate vs. score scatter (with per-case points)
  Fig4_decision_tree       audit workflow diagram
  Fig5_t4_simulation       simulated behaviour of T4/T4r (three panels, new in v2)
  Graphical_Abstract       graphical abstract (landscape, separate file)
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Rectangle
import numpy as np, json, glob, os, re
import pandas as pd

plt.rcParams.update({
    'font.family': 'DejaVu Sans',
    'font.size': 9,
    'axes.labelsize': 9,
    'axes.titlesize': 10,
    'xtick.labelsize': 8,
    'ytick.labelsize': 8,
    'axes.linewidth': 0.8,
    'axes.unicode_minus': False,
    'savefig.dpi': 300,
    'figure.dpi': 300,
})

_HERE = os.path.dirname(os.path.abspath(__file__))

# ---- path resolution: compatible with both the working-directory and the
# submission-package layouts ----------------------------------------------
# Figures are written to the sibling 02_figures/ first (submission-package layout);
# otherwise to ../BIB论文_图表/ (directory name kept as-is).
_CAND_OUT = [
    os.path.join(os.path.dirname(_HERE), '02_figures'),
    '/workspace/BIB论文_图表',
    os.path.join(_HERE, 'figures'),
]
OUT = next((d for d in _CAND_OUT if os.path.isdir(os.path.dirname(d))), _CAND_OUT[-1])
os.makedirs(OUT, exist_ok=True)


def _findfile(name):
    """Locate a file recursively under several candidate roots."""
    roots = [
        _HERE,
        os.path.dirname(_HERE),
        os.path.dirname(os.path.dirname(_HERE)),
        '/workspace',
    ]
    for r in roots:
        if not os.path.isdir(r):
            continue
        for dp, dns, fns in os.walk(r):
            dns[:] = [d for d in dns if d != '__pycache__']
            if name in fns:
                return os.path.join(dp, name)
    raise SystemExit(f'{name} not found; please set AUDIT_ROOT and retry.')


_TABLE = os.path.join(_HERE, '论文数据表', 'paper_table.csv')   # directory name kept as-is
if not os.path.isfile(_TABLE):
    _TABLE = _findfile('paper_table.csv')
df = pd.read_csv(_TABLE)

# ---- simulation results: the reference bands in Fig 1D / Fig 5 must come from
# the same source as the text (sim_t4_t4r.json) ----
_SIM_PATH = os.path.join(_HERE, 'sim_t4_t4r.json')
if not os.path.isfile(_SIM_PATH):
    _SIM_PATH = _findfile('sim_t4_t4r.json')
with open(_SIM_PATH, encoding='utf-8') as _fh:
    SIM = json.load(_fh)
# null_arm: depth span -> mean T4/T4r
NULL_SPAN = sorted(float(k) for k in SIM['null_arm'])
T4_NULL   = {float(k): v['T4']['mean']  for k, v in SIM['null_arm'].items()}
T4R_NULL  = {float(k): v['T4r']['mean'] for k, v in SIM['null_arm'].items()}
T4_NULL_LO = {float(k): v['T4']['min']  for k, v in SIM['null_arm'].items()}
T4_NULL_HI = {float(k): v['T4']['max']  for k, v in SIM['null_arm'].items()}
# signal_arm: number of over-expressed genes -> mean T4/T4r (depth span fixed at its maximum)
SIG_K    = sorted(int(k) for k in SIM['signal_arm'])
T4_SIG   = {int(k): v['T4']['mean']  for k, v in SIM['signal_arm'].items()}
T4R_SIG  = {int(k): v['T4r']['mean'] for k, v in SIM['signal_arm'].items()}
# T4r null band (negative): min and max over all means of null_arm
T4R_NULL_LO = min(T4R_NULL.values())
T4R_NULL_HI = max(T4R_NULL.values())
# T4r signal band: range over all means of signal_arm
T4R_SIG_LO  = min(T4R_SIG.values())
T4R_SIG_HI  = max(T4R_SIG.values())
# reference anchors: the null-band value at depth span 1.0 plus the value at the
# strongest signal (the same values are cited in the figure legend of the text)
T4R_SIM_1   = T4R_NULL[1.0] if 1.0 in T4R_NULL else T4R_NULL[NULL_SPAN[0]]
T4R_SIM_3   = T4R_NULL[max(NULL_SPAN)]
T4R_SIG_10  = T4R_SIG[10] if 10 in T4R_SIG else T4R_SIG[SIG_K[0]]
T4R_SIG_40  = T4R_SIG[max(SIG_K)]

# ---- greyscale + colour-blind-friendly palette ----
C_FAIL  = '#B2182B'   # dark red: failed
C_PASS  = '#2166AC'   # dark blue: passed
C_NEUT  = '#BDBDBD'   # neutral grey
C_HI    = '#000000'   # emphasis

SHORT = {
    'AAA_GSE166676': 'AAA', 'AS_GSE159677_ours': 'AS (this study)',
    'AS_mouse_GSE214611': 'AS (mouse)', 'CAD_GSE131778': 'CAD',
    'DFU_GSE165816': 'DFU', 'PDAC_GSE212966': 'PDAC',
    'ccRCC_GSE156632': 'ccRCC-1', 'ccRCC2_GSE171306': 'ccRCC-2',
}
def shortname(n):
    if n in SHORT: return SHORT[n]
    if n.startswith('TAA_TAA'): return n.split('_')[1]
    if n == 'TAA_GSE155468': return 'TAA1'
    return n

# =====================================================================
# Figure 1 -- three-panel audit overview
# =====================================================================
def fig1():
    fig, axes = plt.subplots(1, 4, figsize=(15.0, 4.2))

    # --- Panel A: T4 rho ---
    a = axes[0]
    sub = df[df.T4_rho.notna()].copy().sort_values('T4_rho')
    labels = [shortname(x) for x in sub.dataset]
    vals = sub.T4_rho.values
    cols = [C_FAIL if v > 0.5 else C_PASS for v in vals]
    y = np.arange(len(vals))
    a.barh(y, vals, color=cols, edgecolor='black', linewidth=0.5, height=0.72)
    a.axvline(0.5, color='black', ls='--', lw=1.2)
    a.text(0.505, len(vals)-0.4, 'pass threshold\n' + r'$\rho$ = 0.5',
           fontsize=7.5, va='top', ha='left')
    a.set_yticks(y); a.set_yticklabels(labels, fontsize=7.5)
    a.set_xlim(0, 1.06); a.set_xlabel(r'Spearman $\rho$ (LRG detection rate vs. score)')
    a.set_title('A  Ranking tracks detection rate', fontsize=9.5, fontweight='bold', loc='left')
    a.grid(axis='x', alpha=0.25, lw=0.5)
    a.set_axisbelow(True)
    for s in ('top','right'): a.spines[s].set_visible(False)

    # --- Panel B: T1 SSR ---
    a = axes[1]
    sub = df[df.T1_SSR.notna()].copy().sort_values('T1_SSR')
    labels = [shortname(x) for x in sub.dataset]
    vals = sub.T1_SSR.values
    cols = [C_FAIL if v <= 1 else C_PASS for v in vals]
    y = np.arange(len(vals))
    a.barh(y, vals, color=cols, edgecolor='black', linewidth=0.5, height=0.72)
    a.axvline(1.0, color='black', ls='--', lw=1.2)
    a.text(1.03, len(vals)-0.4, 'pass threshold\nSSR = 1', fontsize=7.5, va='top', ha='left')
    a.set_yticks(y); a.set_yticklabels(labels, fontsize=7.5)
    a.set_xlabel('SSR = H(lactylation set) / H(housekeeping set)')
    a.set_title('B  Housekeeping genes discriminate as well', fontsize=9.5, fontweight='bold', loc='left')
    a.grid(axis='x', alpha=0.25, lw=0.5); a.set_axisbelow(True)
    for s in ('top','right'): a.spines[s].set_visible(False)

    # --- Panel C: T2 p ---
    a = axes[2]
    sub = df[df.T2_p.notna()].copy().sort_values('T2_p')
    labels = [shortname(x) for x in sub.dataset]
    vals = sub.T2_p.values
    cols = [C_FAIL if v >= 0.05 else C_PASS for v in vals]
    y = np.arange(len(vals))
    a.barh(y, vals, color=cols, edgecolor='black', linewidth=0.5, height=0.72)
    a.axvline(0.05, color='black', ls='--', lw=1.2)
    a.text(0.06, len(vals)-0.4, r'pass threshold' + '\n' + r'$p$ = 0.05', fontsize=7.5, va='top', ha='left')
    a.set_yticks(y); a.set_yticklabels(labels, fontsize=7.5)
    a.set_xlabel('Empirical $p$ (expression-matched null)')
    a.set_title('C  Random gene sets match the effect', fontsize=9.5, fontweight='bold', loc='left')
    a.grid(axis='x', alpha=0.25, lw=0.5); a.set_axisbelow(True)
    for s in ('top','right'): a.spines[s].set_visible(False)

    # --- Panel D: T4r (new in v2) ---
    a = axes[3]
    sub = df[df['T4r'].notna()].copy().sort_values('T4r')
    labels = [shortname(x) for x in sub.dataset]
    vals = sub['T4r'].values.astype(float)
    y = np.arange(len(vals))
    a.barh(y, vals, color=C_NEUT, edgecolor='black', linewidth=0.5, height=0.72)
    a.axvline(0, color='black', lw=1.0)
    # simulation reference band: from sim_t4_t4r.json (same source as the text)
    a.axvspan(T4R_NULL_LO, T4R_NULL_HI, color=C_PASS, alpha=0.16, zorder=0)
    a.axvline(T4R_SIG_LO, color=C_FAIL, ls=':', lw=1.2)
    a.axvline(T4R_SIG_HI, color=C_FAIL, ls=':', lw=1.2)
    a.text(T4R_SIG_HI + 0.005, len(vals) - 0.3, 'simulated\nsignal', fontsize=6.4,
           color=C_FAIL, va='top', ha='left', linespacing=1.2)
    a.set_yticks(y); a.set_yticklabels(labels, fontsize=7.5)
    a.set_xlim(-0.95, 0.30)
    a.set_xlabel(r'$\rho$ after depth residualisation (T4r)')
    a.set_title('D  Residual association is near zero', fontsize=9.5,
                fontweight='bold', loc='left')
    a.grid(axis='x', alpha=0.25, lw=0.5); a.set_axisbelow(True)
    for sp in ('top','right'): a.spines[sp].set_visible(False)

    # legend
    from matplotlib.patches import Patch
    handles = [Patch(facecolor=C_FAIL, edgecolor='black', label='Discriminating test: failed'),
               Patch(facecolor=C_PASS, edgecolor='black', label='Discriminating test: passed'),
               Patch(facecolor=C_NEUT, edgecolor='black', label='T4r (descriptive, not a criterion)')]
    fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=8.5,
               bbox_to_anchor=(0.5, -0.015))

    plt.tight_layout(rect=[0, 0.055, 1, 1])
    p = f'{OUT}/Fig1_audit_landscape.png'
    fig.savefig(p, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return p

# =====================================================================
# Figure 2 -- replication across 8 patients
# =====================================================================
def fig2():
    sub = df[df.patient.notna()].copy()
    order = ['TAA1','TAA2','TAA3','TAA4','TAA5','TAA6','TAA7','TAA8']
    sub = sub.set_index('patient').loc[order].reset_index()

    CT_COLOR = {'Neutrophil': '#2166AC', 'Macrophage': '#B2182B', 'Endothelial': '#4D4D4D'}
    CT_HATCH = {'Neutrophil': '', 'Macrophage': '///', 'Endothelial': '...'}

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.9),
                             gridspec_kw={'width_ratios': [1, 1.15]})

    # --- Panel A: the top cell type flips between patients ---
    a = axes[0]
    y = np.arange(len(order))[::-1]
    for i, (_, r) in enumerate(sub.iterrows()):
        yy = y[i]
        c = r['top']
        a.barh(yy, 1, color=CT_COLOR[c], edgecolor='black',
               linewidth=0.6, height=0.62, hatch=CT_HATCH[c])
        a.text(0.5, yy, f"{c}", ha='center', va='center', fontsize=8,
               color='white', fontweight='bold')
    a.set_yticks(y); a.set_yticklabels(order, fontsize=8.5)
    a.set_xlim(0, 1); a.set_xticks([])
    a.set_xlabel('Top-ranked cell type')
    a.set_title('A  The reported "top cell type" changes between patients',
                fontsize=9.5, fontweight='bold', loc='left')
    for s in ('top','right','bottom'): a.spines[s].set_visible(False)

    from matplotlib.patches import Patch
    h = [Patch(facecolor=CT_COLOR[k], edgecolor='black', hatch=CT_HATCH[k], label=k)
         for k in ['Neutrophil','Macrophage','Endothelial']]
    a.legend(handles=h, fontsize=7.5, frameon=False, loc='lower right',
             bbox_to_anchor=(1.02, -0.42), ncol=3)

    # --- Panel B: meanwhile the T4 rho is highly stable ---
    a = axes[1]
    rho = sub.T4_rho.values
    y = np.arange(len(order))[::-1]
    a.barh(y, rho, color=C_FAIL, edgecolor='black', linewidth=0.5, height=0.62)
    a.axvline(0.5, color='black', ls='--', lw=1.2)
    a.text(0.515, len(order)-0.5, r'pass threshold $\rho$ = 0.5', fontsize=7.5, va='top')
    a.axvline(rho.mean(), color='black', ls=':', lw=1.4)
    a.text(rho.mean()-0.015, -1.05, f'mean = {rho.mean():+.3f}', fontsize=7.5,
           ha='right', va='center')
    for i, v in enumerate(rho):
        a.text(v-0.02, y[i], f'{v:+.3f}', ha='right', va='center', fontsize=7.5,
               color='white', fontweight='bold')
    a.set_yticks(y); a.set_yticklabels(order, fontsize=8.5)
    a.set_xlim(0, 1.06)
    a.set_xlabel(r'Spearman $\rho$ (detection rate vs. score)')
    a.set_title('B  The technical signal itself is stable',
                fontsize=9.5, fontweight='bold', loc='left')
    a.grid(axis='x', alpha=0.25, lw=0.5); a.set_axisbelow(True)
    for s in ('top','right'): a.spines[s].set_visible(False)

    plt.tight_layout()
    p = f'{OUT}/Fig2_patient_replication.png'
    fig.savefig(p, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return p

# =====================================================================
# Figure 3 -- per cell type: detection rate vs. score (examples: ccRCC and TAA4)
# =====================================================================
def fig3():
    sets = [('ccRCC_GSE156632', 'A  ccRCC (GSE156632)'),
            ('TAA_TAA4_GSE155468', 'B  Thoracic aortic aneurysm, patient TAA4')]
    fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.0))
    for ax, (name, title) in zip(axes, sets):
        f = f'/root/audit/results/{name}_T4.csv'
        if not os.path.exists(f):
            ax.axis('off'); continue
        d = pd.read_csv(f)
        if len(d) < 4:
            ax.axis('off'); continue
        x = d['detect_rate_lrg'].values
        yv = d['lact_score'].values
        n = d['n'].values
        sz = 22 + 130 * (n / n.max()) ** 0.5
        ax.scatter(x, yv, s=sz, c=C_FAIL, edgecolors='black', linewidths=0.7,
                   zorder=3, alpha=0.9)
        # label every point
        for xi, yi, ci in zip(x, yv, d['celltype']):
            ax.annotate(ci, (xi, yi), fontsize=6.8,
                        xytext=(4, 3), textcoords='offset points', zorder=4)
        # linear fit as a visual guide
        if len(x) >= 3:
            z = np.polyfit(x, yv, 1)
            xs = np.linspace(x.min(), x.max(), 50)
            ax.plot(xs, np.polyval(z, xs), color='black', ls='--', lw=1.0, zorder=2)
            r = np.corrcoef(x, yv)[0, 1]
        ax.set_xlabel('Fraction of LRG detected per cell')
        ax.set_ylabel('Mean lactylation score')
        ax.set_title(title, fontsize=9.5, fontweight='bold', loc='left')
        ax.grid(alpha=0.25, lw=0.5); ax.set_axisbelow(True)
        for s in ('top','right'): ax.spines[s].set_visible(False)
        # annotate rho
        rho = d['detect_rate_lrg'].corr(d['lact_score'], method='spearman')
        ax.text(0.97, 0.06, 'Spearman ' + r'$\rho$' + f' = {rho:+.3f} (annotation only)',
                transform=ax.transAxes, ha='right', va='bottom', fontsize=7.2,
                bbox=dict(boxstyle='round,pad=0.3', fc='white', ec='black', lw=0.6))
    plt.tight_layout()
    p = f'{OUT}/Fig3_detection_vs_score.png'
    fig.savefig(p, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return p

# =====================================================================
# Figure 4 -- audit decision workflow
# =====================================================================
def fig4():
    fig, ax = plt.subplots(figsize=(7.6, 6.4))
    ax.set_xlim(0, 10); ax.set_ylim(0, 11.2); ax.axis('off')

    def box(x, y, w, h, text, fc='white', ec='black', fs=8.3, bold=False, lw=1.0):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
                     boxstyle='round,pad=0.16,rounding_size=0.18',
                     fc=fc, ec=ec, lw=lw, zorder=2))
        ax.text(x + w/2, y + h/2, text, ha='center', va='center', fontsize=fs,
                zorder=3, fontweight='bold' if bold else 'normal', linespacing=1.45)

    def arrow(x1, y1, x2, y2, color='black', style='-|>'):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2),
                     arrowstyle=style, mutation_scale=11, lw=1.1,
                     color=color, zorder=1, shrinkA=1, shrinkB=1))

    box(2.2, 10.1, 5.6, 0.85,
        'Published single-cell lactylation study\n' +
        '(LRG gene set + AddModuleScore + cell-type ranking)',
        fc='#E8E8E8', fs=8.3, bold=True)

    arrow(5.0, 10.1, 5.0, 9.55)
    box(2.6, 8.75, 4.8, 0.8,
        'STEP 1  Faithful replication\nsame data, same scoring, same ranking',
        fc='white', fs=8.3)

    arrow(5.0, 8.75, 5.0, 8.62)
    box(3.4, 7.85, 3.2, 0.62,
        'STEP 2  Apply four orthogonal controls',
        fc='#DDEEFF', fs=8.4, bold=True)

    # four branches
    ys = 6.05
    bw, bh = 2.15, 1.55
    xs = [0.35, 2.75, 5.15, 7.55]
    labels = [
        ('T1\nHousekeeping\ncontrol\n\npass if SSR > 1', 'H(real) / H(hk)'),
        ('T2\nExpression-matched\nnull\n\npass if p < 0.05', '150 permutations'),
        ('T3\nTechnical\ncovariates\n\ndiagnostic only', 'depth, nGene, HK'),
        ('T4\nCross-cell-type\ndetection rate\n\npass if ' + r'$\rho$' + ' \u2264 0.5', 'sparsity collinearity'),
    ]
    for x, (lab, sub) in zip(xs, labels):
        box(x, ys, bw, bh, lab, fc='white', fs=7.0)
        arrow(5.0, 7.85, x + bw/2, ys + bh)
        ax.text(x + bw/2, ys - 0.38, sub, ha='center', va='top',
                fontsize=6.4, style='italic', color='#444444')

    # merge
    y_merge = 4.25
    for x in xs:
        arrow(x + bw/2, ys - 0.38, x + bw/2, y_merge + 0.6)
    ax.plot([xs[0]+bw/2, xs[-1]+bw/2], [y_merge+0.6, y_merge+0.6],
            color='black', lw=1.1, zorder=1)
    arrow(5.0, y_merge + 0.6, 5.0, y_merge)

    box(2.7, 3.45, 4.6, 0.8, 'All four controls passed?', fc='#DDEEFF', fs=8.6, bold=True)

    # yes / no
    arrow(3.6, 3.45, 2.15, 2.75, color=C_PASS)
    arrow(6.4, 3.45, 7.85, 2.75, color=C_FAIL)
    ax.text(2.35, 3.22, 'YES', fontsize=7.6, color=C_PASS, fontweight='bold', ha='center')
    ax.text(7.55, 3.22, 'NO', fontsize=7.6, color=C_FAIL, fontweight='bold', ha='center')

    box(0.55, 1.55, 3.2, 1.2,
        'Signal is specific\nto the lactylation gene set.\nConclusion stands.',
        fc='#D6E6F5', ec=C_PASS, fs=8.0, bold=True)
    box(6.25, 1.55, 3.2, 1.2,
        'Ranking is driven by\ntechnical depth /\ndetection rate.',
        fc='#F5D6D6', ec=C_FAIL, fs=8.0, bold=True)

    arrow(2.15, 1.55, 2.15, 1.05)
    arrow(7.85, 1.55, 7.85, 1.05)
    box(0.75, 0.15, 8.5, 0.9,
        'For a scoring-based claim to be interpretable,\n' +
        'a housekeeping control must not perform equally well.',
        fc='#F0F0F0', fs=8.3, bold=True)

    plt.tight_layout()
    p = f'{OUT}/Fig4_audit_workflow.png'
    fig.savefig(p, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return p

# =====================================================================
# Graphical Abstract (landscape, separate file)
# =====================================================================
def gabs():
    fig, ax = plt.subplots(figsize=(9.0, 4.6))
    ax.set_xlim(0, 12); ax.set_ylim(0, 6.2); ax.axis('off')

    ax.text(6.0, 5.85, 'Housekeeping genes audit single-cell lactylation claims',
            ha='center', va='top', fontsize=12, fontweight='bold')

    # left: current practice
    ax.add_patch(FancyBboxPatch((0.25, 0.55), 3.4, 4.5,
                 boxstyle='round,pad=0.15,rounding_size=0.2',
                 fc='#F7F7F7', ec='black', lw=1.1))
    ax.text(1.95, 4.72, 'Current practice', ha='center', fontsize=9.5, fontweight='bold')
    steps = ['GEO single-cell\ndataset', '300+ lactylation-related\ngenes (LRG)',
             'AddModuleScore', 'Rank cell types', '"X cells are the most\nlactylation-active"']
    yy = 4.28
    for i, s in enumerate(steps):
        ax.add_patch(FancyBboxPatch((0.5, yy-0.42), 2.9, 0.55,
                     boxstyle='round,pad=0.08,rounding_size=0.12',
                     fc='white', ec='#666666', lw=0.8))
        ax.text(1.95, yy-0.145, s, ha='center', va='center', fontsize=6.8)
        if i < len(steps)-1:
            ax.annotate('', xy=(1.95, yy-0.44), xytext=(1.95, yy-0.48),
                        arrowprops=dict(arrowstyle='-|>', color='#666666', lw=0.9))
        yy -= 0.72

    # middle: audit
    ax.add_patch(FancyBboxPatch((4.05, 0.55), 3.5, 4.5,
                 boxstyle='round,pad=0.15,rounding_size=0.2',
                 fc='#EAF2FA', ec='black', lw=1.1))
    ax.text(5.8, 4.72, 'Four controls (T1-T4)', ha='center', fontsize=9.5, fontweight='bold')
    items = [('T1', 'Housekeeping set\nreplaces LRG'),
             ('T2', 'Expression-matched\nrandom gene sets'),
             ('T3', 'Depth / nGene\ncovariates'),
             ('T4', 'Detection rate vs.\ncell-type ranking')]
    yy = 4.15
    for k, s in items:
        ax.add_patch(FancyBboxPatch((4.3, yy-0.72), 3.0, 0.78,
                     boxstyle='round,pad=0.1,rounding_size=0.12',
                     fc='white', ec='#2166AC', lw=0.9))
        ax.text(4.55, yy-0.33, k, ha='center', va='center', fontsize=9,
                fontweight='bold', color='#2166AC')
        ax.text(5.65, yy-0.33, s, ha='left', va='center', fontsize=6.8)
        yy -= 0.88

    # right: results
    ax.add_patch(FancyBboxPatch((7.85, 0.55), 3.9, 4.5,
                 boxstyle='round,pad=0.15,rounding_size=0.2',
                 fc='#FBEFEF', ec='black', lw=1.1))
    ax.text(9.8, 4.86, '16 datasets  ·  113,419 cells', ha='center', fontsize=9.0, fontweight='bold')

    _t4r_real_med = float(df['T4r'].median())
    _n_t4r = int(df['T4r'].notna().sum())
    stats = [('9 / 16  ·  11 / 16', 'housekeeping panel and\nexpression-matched null\nfail the lactylation set'),
             (f'{T4_NULL[max(NULL_SPAN)]:+.3f}', 'T4 under a null with no\nbiology, once depth\nvaries by 3.0 log units'),
             (f'{_t4r_real_med:+.3f}', f'median residualised\nassociation (T4r),\n{_n_t4r} datasets'),
             ('4 / 3 / 1', 'top cell type across 8\npatients: macrophage /\nneutrophil / endothelial')]
    yy = 4.10
    for big, small in stats:
        ax.text(9.8, yy, big, ha='center', va='center', fontsize=11.5,
                fontweight='bold', color='#B2182B')
        ax.text(9.8, yy-0.42, small, ha='center', va='center', fontsize=6.6,
                color='#333333', linespacing=1.35)
        yy -= 0.92

    # arrows
    ax.annotate('', xy=(4.02, 2.8), xytext=(3.68, 2.8),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=1.6))
    ax.annotate('', xy=(7.82, 2.8), xytext=(7.58, 2.8),
                arrowprops=dict(arrowstyle='-|>', color='black', lw=1.6))

    ax.text(6.0, 0.22,
            'T4 cannot be failed by any gene set on realistic data; report it, do not use it as a criterion.',
            ha='center', va='bottom', fontsize=8.3, style='italic', color='#222222')

    plt.tight_layout()
    p = f'{OUT}/Graphical_Abstract.png'
    fig.savefig(p, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return p


# =====================================================================
# Figure 5 -- simulated behaviour of T4 / T4r (new in v2)
# =====================================================================
def fig5():
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 3.9))

    # ---- Panel A: T4 vs depth span ----
    a = axes[0]
    spans = NULL_SPAN
    t4m   = [T4_NULL[s] for s in spans]
    t4lo  = [T4_NULL_LO[s] for s in spans]
    t4hi  = [T4_NULL_HI[s] for s in spans]
    x = np.arange(len(spans))
    a.errorbar(x, t4m, yerr=[np.array(t4m)-np.array(t4lo),
                             np.array(t4hi)-np.array(t4m)],
               fmt='o-', color=C_FAIL, ecolor=C_FAIL, elinewidth=1.2,
               capsize=4, markersize=6, lw=1.6, label=r'T4 ($\rho$)')
    a.axhline(0.5, color='black', ls='--', lw=1.1)
    a.text(0.03, 0.53, r'original pass threshold $\rho$ = 0.5',
           fontsize=7.2, va='bottom', ha='left')
    a.plot(x, [T4R_NULL[s] for s in spans], 's--', color=C_PASS,
           markersize=6, lw=1.4, markeredgecolor='black', markeredgewidth=0.6,
           label=r'T4r ($\rho$)')
    a.set_xticks(x)
    a.set_xticklabels([f'{s:g}' for s in spans])
    a.set_xlabel('Imposed depth span across cell types\n(log units)')
    a.set_ylabel(r'Spearman $\rho$')
    a.set_ylim(min(-0.05, min(T4R_NULL.values()) - 0.10), 1.12)
    a.set_title('A  No biology; depth alone saturates T4', fontsize=9.5,
                fontweight='bold', loc='left')
    a.legend(fontsize=7.2, frameon=False, loc='center right')
    a.grid(alpha=0.25, lw=0.5); a.set_axisbelow(True)
    for sp in ('top','right'): a.spines[sp].set_visible(False)

    # ---- Panel B: T4/T4r vs signal strength ----
    a = axes[1]
    sig = SIG_K
    t4sig  = [T4_SIG[k] for k in sig]
    # T4r: null-signal anchor + means at each signal strength (same depth span)
    t4r_null   = [T4R_NULL_HI] * len(sig)
    t4r_signal = [T4R_SIG[k] for k in sig]
    a.plot(sig, t4sig, 's--', color=C_NEUT, markersize=6, lw=1.4,
           markeredgecolor='black', markeredgewidth=0.6, label=r'T4 (flat)')
    a.plot(sig, t4r_signal, 'o-', color=C_PASS, markersize=6, lw=1.8,
           markeredgecolor='black', markeredgewidth=0.6, label='T4r with signal')
    a.plot(sig, t4r_null, ':', color=C_PASS, lw=1.4, label='T4r, no signal')
    a.set_xlabel('Number of genes over-expressed\nin one cell type')
    a.set_ylabel(r'$\rho$')
    a.set_ylim(min(-1.05, min(t4r_signal) - 0.12), 1.05)
    a.set_title('B  Signal moves T4r but not T4', fontsize=9.5,
                fontweight='bold', loc='left')
    a.legend(fontsize=7.2, frameon=False, loc='center left')
    a.grid(alpha=0.25, lw=0.5); a.set_axisbelow(True)
    for sp in ('top','right'): a.spines[sp].set_visible(False)

    # ---- Panel C: real T4r vs simulated bands ----
    a = axes[2]
    sub = df[df['T4r'].notna()].copy().sort_values('T4r')
    vals = sub['T4r'].values.astype(float)
    labels = [shortname(x) for x in sub.dataset]
    y = np.arange(len(vals))
    a.axvspan(T4R_NULL_LO, T4R_NULL_HI, color=C_PASS, alpha=0.16, zorder=0,
              label='simulated null band')
    a.axvline(T4R_SIG_LO, color=C_FAIL, ls=':', lw=1.3)
    a.axvline(T4R_SIG_HI, color=C_FAIL, ls=':', lw=1.3,
              label='simulated signal (low / high)')
    a.scatter(vals, y, s=34, color=C_NEUT, edgecolor='black',
              linewidth=0.6, zorder=3)
    a.axvline(0, color='black', lw=1.0, zorder=2)
    a.set_yticks(y); a.set_yticklabels(labels, fontsize=7.2)
    a.set_xlim(-0.95, 1.35)
    a.set_xlabel(r'T4r ($\rho$)')
    a.set_title('C  Real datasets sit above the null band', fontsize=9.5,
                fontweight='bold', loc='left')
    a.legend(fontsize=6.8, frameon=False, loc='lower right')
    a.grid(axis='x', alpha=0.25, lw=0.5); a.set_axisbelow(True)
    for sp in ('top','right'): a.spines[sp].set_visible(False)

    plt.tight_layout()
    p = f'{OUT}/Fig5_t4_simulation.png'
    fig.savefig(p, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return p


def fig6():
    """Fig6 -- depth coupling and ranking disagreement under three scoring implementations.

    (A) T3 (within-group depth coupling) for 13 datasets under three estimators.
    (B) Datasets where the three estimators name different top cell types.
    """
    src = _findfile('estimator_comparison_extended.json')
    recs = json.load(open(src, encoding='utf-8'))
    df = pd.DataFrame(recs)
    df['short'] = df['dataset'].str.replace('_', ' ', regex=False)
    df.loc[df['dataset'].str.startswith('GSE155468'), 'short'] = \
        df['dataset'].str.replace('GSE155468_', '', regex=False)

    ORDER = ['score_z', 'scanpy.score_genes', 'rank_UCell']
    NICE = {'score_z': 'score_z (ours)',
            'scanpy.score_genes': 'AddModuleScore-style',
            'rank_UCell': 'UCell-style'}
    # hue + lightness double encoding, avoiding red-green contrast
    COL = {'score_z': '#1f3b73', 'scanpy.score_genes': '#5b8ac4',
           'rank_UCell': '#a8c6e8'}
    HAT = {'score_z': '', 'scanpy.score_genes': '///', 'rank_UCell': '...'}

    piv = df.pivot_table(index='dataset', columns='estimator', values='T3')
    piv = piv.reindex(columns=[c for c in ORDER if c in piv.columns])
    # sort by score_z in descending order for easier reading
    order = piv['score_z'].sort_values(ascending=False).index.tolist()
    piv = piv.loc[order]

    fig, (a, b) = plt.subplots(1, 2, figsize=(11.4, 4.4),
                               gridspec_kw={'width_ratios': [1.45, 1.0]})

    # ---- panel A: grouped bars, T3 declines across estimators ----
    y = np.arange(len(piv))
    h = 0.26
    for i, est in enumerate(piv.columns):
        a.barh(y + (i - 1) * h, piv[est].values, height=h,
               color=COL[est], edgecolor='#2b2b2b', linewidth=0.5,
               hatch=HAT[est], label=NICE[est], zorder=3)
    a.set_yticks(y)
    a.set_yticklabels([str(x).replace('GSE131778 CAD', 'CAD')
                       for x in piv.index], fontsize=5.6)
    a.invert_yaxis()
    a.set_xlabel('T3: within-cell correlation of score with detected-gene count')
    a.set_title('A  Depth coupling falls under other estimators but stays positive',
                fontsize=9.2, fontweight='bold', loc='left')
    a.axvline(0, color='#444', lw=0.8)
    a.axvline(0.5, color='#b03a2e', lw=0.9, ls='--', zorder=1)
    a.text(0.505, len(piv) - 0.3, '0.5', color='#b03a2e', fontsize=6.5, va='bottom')
    a.legend(fontsize=6.6, frameon=False, loc='lower right')
    a.grid(axis='x', alpha=0.22, lw=0.5); a.set_axisbelow(True)
    for sp in ('top', 'right'): a.spines[sp].set_visible(False)

    # ---- panel B: datasets where the three estimators disagree on the ranking ----
    flips = []
    for ds, g in df.groupby('dataset'):
        tops = list(g.set_index('estimator')['top'].reindex(
            [c for c in ORDER if c in set(g['estimator'])]).values)
        if len(set(tops)) > 1:
            flips.append((ds, tops))
    b.axis('off')
    if flips:
        rows = []
        for ds, tops in flips:
            lab = ds.replace('GSE155468_', '')
            rows.append([lab] + [str(t) for t in tops])
        tbl = b.table(cellText=rows,
                      colLabels=['Dataset', 'score_z', 'AddModule', 'UCell'],
                      cellLoc='center', loc='center')
        tbl.auto_set_font_size(False); tbl.set_fontsize(6.2)
        tbl.scale(1.0, 1.28)
        for (r, c), cell in tbl.get_celld().items():
            cell.set_linewidth(0.4)
            if r == 0:
                cell.set_facecolor('#e8eef7'); cell.set_text_props(fontweight='bold')
            elif c == 0:
                cell.set_facecolor('#f5f5f5')
    b.set_title(f'B  {len(flips)} of {df["dataset"].nunique()} datasets: '
                'the three estimators name different top cell types',
                fontsize=9.2, fontweight='bold', loc='left')

    plt.tight_layout()
    p = f'{OUT}/Fig6_estimator_coupling.png'
    fig.savefig(p, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return p


if __name__ == '__main__':
    for f in (fig1, fig2, fig3, fig4, fig5, fig6, gabs):
        print('saved:', f())
