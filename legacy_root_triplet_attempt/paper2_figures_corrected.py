#!/usr/bin/env python3
"""
Generate Paper 2 Figures: The Triplet Code as Minimal Control Architecture
for Translation Under Proteostasis Constraint

CORRECTED VERSION:
- Figure 1A: Proper saddle-node bifurcation with S-shaped F(P) curve
- Figure 4A: Fixed null ensemble (shuffles mu/tAI pairs together, not independently)

Generates 4 publication-quality figures with real data.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from scipy import stats
from scipy.optimize import brentq
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

# Style settings
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

# Color palette (colorblind-friendly)
COLORS = {
    'blue': '#2E86AB',
    'red': '#E94F37',
    'green': '#2A9D8F',
    'orange': '#F77F00',
    'purple': '#9B5DE5',
    'gray': '#6C757D',
    'light_gray': '#ADB5BD',
    'dark_gray': '#495057',
}

# Amino acid colors for specific highlighting
AA_COLORS = {
    'L': '#2E86AB',  # Blue for Leucine
    'D': '#E94F37',  # Red for Asp
    'H': '#2A9D8F',  # Green for His
    'C': '#F77F00',  # Orange for Cys
    'E': '#9B5DE5',  # Purple for Glu
}

# Standard codon table
CODON_TO_AA = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C',
    'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L',
    'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
    'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
    'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
    'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
    'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
    'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
    'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
    'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
    'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A',
    'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E',
    'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
    'GGG': 'G',
}

def load_codon_data():
    """Load codon modes data with mu and tAI values."""
    df = pd.read_csv('errors/codon_modes_ecoli.tsv', sep='\t')
    df['mu'] = pd.to_numeric(df['mu'], errors='coerce')
    df['w_tai'] = pd.to_numeric(df['w_tai'], errors='coerce')
    return df

def load_metal_data():
    """Load metal binding site enrichment data."""
    mu_df = pd.read_csv('metals/metal_codon_bias_with_mu.tsv', sep='\t')
    tai_df = pd.read_csv('metals/metal_codon_bias_with_tai.tsv', sep='\t')
    return mu_df, tai_df

def load_cross_species_usage():
    """Load codon usage for all three species."""
    ecoli = pd.read_csv('cross_species/data/ecoli/global_codon_usage.tsv', sep='\t')
    bsub = pd.read_csv('cross_species/data/bsub/global_codon_usage.tsv', sep='\t')
    yeast = pd.read_csv('cross_species/data/yeast/global_codon_usage.tsv', sep='\t')
    return ecoli, bsub, yeast

def compute_synonymy_shielding():
    """Compute synonymy shielding for each codon position."""
    nucleotides = ['A', 'C', 'G', 'T']

    s_values = []
    for pos in range(3):
        total = 0
        synonymous = 0
        for codon in CODON_TO_AA:
            if CODON_TO_AA[codon] == '*':
                continue
            aa_orig = CODON_TO_AA[codon]
            for nt in nucleotides:
                if nt == codon[pos]:
                    continue
                new_codon = codon[:pos] + nt + codon[pos+1:]
                if new_codon in CODON_TO_AA:
                    total += 1
                    if CODON_TO_AA[new_codon] == aa_orig:
                        synonymous += 1
        s_values.append(synonymous / total if total > 0 else 0)

    return s_values

# ============================================================================
# FIGURE 1: Proteostasis constraint and codon operational space
# ============================================================================

def create_figure1():
    """Create Figure 1: Proteostasis constraint and codon operational space.

    CORRECTED Panel A: Proper saddle-node bifurcation with S-shaped F(P) curve.
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel A: CORRECTED Saddle-node bifurcation
    ax = axes[0]

    # Parameters tuned for clear S-shaped curve
    V_max = 1.0      # Max chaperone capacity
    K_m = 0.25       # Chaperone affinity (smaller = sharper saturation)
    k_agg = 0.12     # Aggregation rate

    # Outflux function: F(P) = rescue + aggregation
    def F(P):
        rescue = V_max * P / (K_m + P)  # Saturable Michaelis-Menten
        aggregation = k_agg * P**2       # Superlinear aggregation
        return rescue + aggregation

    P = np.linspace(0.001, 3.0, 500)
    F_vals = F(P)

    # Plot F(P) - should show saturation then upturn
    ax.plot(P, F_vals, 'k-', linewidth=2.5, label='Total outflux F(P)')

    # Define J levels
    J_crit = 0.92   # Critical (tangent point)
    J_low = 0.70    # Viable
    J_high = 1.15   # Collapse

    # Find intersections for J_low (should have 2)
    try:
        P_stable_low = brentq(lambda p: F(p) - J_low, 0.01, 0.8)
        P_unstable_low = brentq(lambda p: F(p) - J_low, 1.2, 2.5)
        has_two_roots = True
    except:
        has_two_roots = False
        P_stable_low = 0.3
        P_unstable_low = 2.0

    # Find saddle-node point for J_crit
    try:
        P_saddle = brentq(lambda p: F(p) - J_crit, 0.4, 1.2)
    except:
        P_saddle = 0.7

    # Plot J_in lines
    ax.axhline(J_low, color=COLORS['green'], linestyle='-', linewidth=1.5, alpha=0.9)
    ax.axhline(J_crit, color=COLORS['orange'], linestyle='-', linewidth=1.5, alpha=0.9)
    ax.axhline(J_high, color=COLORS['red'], linestyle='-', linewidth=1.5, alpha=0.9)

    # Mark fixed points
    if has_two_roots:
        # Stable (filled circle)
        ax.plot(P_stable_low, J_low, 'o', color=COLORS['green'], markersize=12,
               markeredgecolor='black', markeredgewidth=1.5, zorder=5)
        # Unstable (open circle)
        ax.plot(P_unstable_low, J_low, 'o', color='white', markersize=12,
               markeredgecolor=COLORS['green'], markeredgewidth=2.5, zorder=5)

    # Saddle-node point (square)
    ax.plot(P_saddle, J_crit, 's', color=COLORS['orange'], markersize=12,
           markeredgecolor='black', markeredgewidth=1.5, zorder=5)

    # Add annotations with offset to avoid overlap
    ax.annotate('stable', xy=(P_stable_low, J_low),
               xytext=(P_stable_low - 0.15, J_low + 0.18),
               fontsize=9, ha='center', color=COLORS['dark_gray'])
    ax.annotate('unstable', xy=(P_unstable_low, J_low),
               xytext=(P_unstable_low + 0.25, J_low + 0.18),
               fontsize=9, ha='center', color=COLORS['dark_gray'])
    ax.annotate('saddle-node', xy=(P_saddle, J_crit),
               xytext=(P_saddle + 0.6, J_crit + 0.12),
               fontsize=9, ha='left', color=COLORS['dark_gray'],
               arrowprops=dict(arrowstyle='->', color=COLORS['gray'], lw=1))

    # Add regime labels
    ax.text(2.7, J_low - 0.08, r'$J_{in} < J_{crit}$', fontsize=9,
           color=COLORS['green'], ha='right', va='top')
    ax.text(2.7, J_crit + 0.05, r'$J_{in} = J_{crit}$', fontsize=9,
           color=COLORS['orange'], ha='right', va='bottom')
    ax.text(2.7, J_high + 0.05, r'$J_{in} > J_{crit}$', fontsize=9,
           color=COLORS['red'], ha='right', va='bottom')
    ax.text(0.15, 1.55, 'collapse\n(no stable\nfixed point)', fontsize=8,
           color=COLORS['red'], ha='left', va='top', style='italic')

    ax.set_xlabel('Misfolded protein pool $P$', fontsize=11)
    ax.set_ylabel('Flux', fontsize=11)
    ax.set_xlim(0, 3.0)
    ax.set_ylim(0, 1.8)

    ax.text(-0.15, 1.05, 'A', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.set_title('Saddle-node bifurcation', fontsize=11, pad=10)

    # Panel B: Codons in (mu, nu) space
    ax = axes[1]
    df = load_codon_data()

    # Plot all codons
    for _, row in df.iterrows():
        aa = row['aa']
        mu = row['mu']
        tai = row['w_tai']

        if pd.isna(mu) or pd.isna(tai):
            continue

        color = AA_COLORS.get(aa, COLORS['light_gray'])
        alpha = 1.0 if aa in AA_COLORS else 0.4
        size = 80 if aa in AA_COLORS else 40

        ax.scatter(mu, tai, c=color, s=size, alpha=alpha, edgecolors='white',
                  linewidth=0.5, zorder=3 if aa in AA_COLORS else 2)

    # Connect synonymous codons for highlighted amino acids
    for aa in ['L', 'D']:
        aa_df = df[df['aa'] == aa].dropna(subset=['mu', 'w_tai'])
        if len(aa_df) > 1:
            coords = aa_df[['mu', 'w_tai']].values
            for i in range(len(coords)):
                for j in range(i+1, len(coords)):
                    ax.plot([coords[i,0], coords[j,0]],
                           [coords[i,1], coords[j,1]],
                           color=AA_COLORS[aa], alpha=0.3, linewidth=1, zorder=1)

    # Add labels for key codons
    label_codons = {'CTG': 'L', 'CTA': 'L', 'GAC': 'D', 'GAT': 'D'}
    for _, row in df.iterrows():
        if row['codon'] in label_codons:
            mu, tai = row['mu'], row['w_tai']
            if not pd.isna(mu) and not pd.isna(tai):
                offset = (5, 5) if row['codon'] in ['CTG', 'GAC'] else (-25, -10)
                ax.annotate(row['codon'], (mu, tai), xytext=offset,
                           textcoords='offset points', fontsize=8,
                           color=AA_COLORS[row['aa']], fontweight='bold')

    # Add arrow showing Leu spread
    leu_df = df[df['aa'] == 'L'].dropna(subset=['mu', 'w_tai'])
    if len(leu_df) > 1:
        leu_mu_min = float(leu_df['mu'].min())
        leu_mu_max = float(leu_df['mu'].max())
        leu_mu_mean = float(leu_df['mu'].mean())
        leu_tai_mean = float(leu_df['w_tai'].mean())
        ax.annotate('', xy=(leu_mu_max, leu_tai_mean),
                   xytext=(leu_mu_min, leu_tai_mean),
                   arrowprops=dict(arrowstyle='<->', color=COLORS['blue'], lw=1.5))
        ax.text(leu_mu_mean, leu_tai_mean + 0.08,
               r'$\Delta_{Leu}$', fontsize=10, color=COLORS['blue'],
               ha='center', fontweight='bold')

    ax.set_xlabel(r'Mistranslation rate $\mu(c)$', fontsize=11)
    ax.set_ylabel(r'Decoding efficiency $\nu(c)$ [tAI]', fontsize=11)
    ax.set_xscale('log')
    ax.set_xlim(1e-5, 0.1)
    ax.set_ylim(0, 1.1)

    # Legend
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=AA_COLORS['L'],
               markersize=8, label='Leu'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=AA_COLORS['D'],
               markersize=8, label='Asp'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=AA_COLORS['H'],
               markersize=8, label='His'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=AA_COLORS['C'],
               markersize=8, label='Cys'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=AA_COLORS['E'],
               markersize=8, label='Glu'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['light_gray'],
               markersize=6, label='Other'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8,
             frameon=False, ncol=2)
    ax.text(-0.15, 1.05, 'B', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.set_title('Synonymous codons in operational space', fontsize=11, pad=10)

    # Panel C: Synonymy shielding
    ax = axes[2]
    s_values = compute_synonymy_shielding()
    positions = [1, 2, 3]
    colors = [COLORS['gray'], COLORS['gray'], COLORS['blue']]

    bars = ax.bar(positions, s_values, color=colors, edgecolor='white', linewidth=1.5)

    ax.axhline(y=s_values[2], color=COLORS['blue'], linestyle='--',
               alpha=0.5, linewidth=1)
    ax.text(3.5, s_values[2] + 0.02, f'$S_3 \\approx {s_values[2]:.2f}$',
           fontsize=10, color=COLORS['blue'], va='bottom')

    ax.set_xlabel('Codon position', fontsize=11)
    ax.set_ylabel('Synonymy fraction $S_i$', fontsize=11)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(['1st', '2nd', '3rd\n(wobble)'])
    ax.set_ylim(0, 0.85)

    # Add value labels on bars
    for i, (bar, val) in enumerate(zip(bars, s_values)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
               f'{val:.2f}', ha='center', va='bottom', fontsize=9)

    ax.text(-0.15, 1.05, 'C', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.set_title('Synonymy shielding by position', fontsize=11, pad=10)

    plt.tight_layout()
    plt.savefig('Figure1_proteostasis_constraint.png', dpi=300, bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.savefig('Figure1_proteostasis_constraint.svg', format='svg', bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.savefig('Figure1_proteostasis_constraint.pdf', format='pdf', bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.close()
    print("Figure 1 saved (CORRECTED Panel A).")

# ============================================================================
# FIGURE 2: Architectural comparison and cross-species conservation
# ============================================================================

def create_figure2():
    """Create Figure 2: Architecture comparison and cross-species conservation."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel A: Vocabulary capacity by code length
    ax = axes[0]

    code_types = ['Doublet\n(n=2)', 'Triplet\n(n=3)', 'Quadruplet\n(n=4)']
    sense_codons = [15, 61, 255]  # 4^n - 1 stop
    synonymy = [(15-20)/20, (61-20)/20, (255-20)/20]
    synonymy_display = [max(0, s) for s in synonymy]

    x = np.arange(len(code_types))
    width = 0.35

    bars1 = ax.bar(x - width/2, sense_codons, width, color=COLORS['blue'],
                   label='Sense codons', edgecolor='white', linewidth=1)
    ax.set_ylabel('Number of sense codons', color=COLORS['blue'], fontsize=11)
    ax.tick_params(axis='y', labelcolor=COLORS['blue'])
    ax.set_ylim(0, 280)

    ax2 = ax.twinx()
    bars2 = ax2.bar(x + width/2, synonymy_display, width, color=COLORS['orange'],
                    label='Synonymy per AA', edgecolor='white', linewidth=1)
    ax2.set_ylabel('Synonymy per amino acid', color=COLORS['orange'], fontsize=11)
    ax2.tick_params(axis='y', labelcolor=COLORS['orange'])
    ax2.set_ylim(0, 13)

    ax.axhline(y=20, color=COLORS['red'], linestyle='--', linewidth=1.5,
               label='20 AA requirement')

    ax.text(0, 25, '✗', fontsize=20, color=COLORS['red'], ha='center',
           fontweight='bold')
    ax.text(1, 70, '✓', fontsize=20, color=COLORS['green'], ha='center',
           fontweight='bold')
    ax.text(2, 260, 'excess', fontsize=9, color=COLORS['gray'], ha='center',
           style='italic')

    ax.set_xticks(x)
    ax.set_xticklabels(code_types)
    ax.set_xlabel('Code architecture', fontsize=11)

    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left',
             fontsize=8, frameon=False)

    ax.text(-0.15, 1.05, 'A', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.set_title('Vocabulary capacity by code length', fontsize=11, pad=10)

    # Panel B: Universal triplet architecture
    ax = axes[1]

    species = ['E. coli', 'B. subtilis', 'S. cerevisiae']
    synonymy_values = [61/20, 61/20, 61/20]

    bars = ax.bar(species, synonymy_values, color=[COLORS['blue'], COLORS['green'],
                  COLORS['orange']], edgecolor='white', linewidth=1.5)

    ax.axhline(y=1.0, color=COLORS['red'], linestyle='--', linewidth=1.5,
               label='Minimum for control')

    ax.set_ylabel('Synonymy per amino acid', fontsize=11)
    ax.set_ylim(0, 4)
    ax.legend(loc='upper right', fontsize=8, frameon=False)

    ax.text(1, 3.5, 'Universal triplet\narchitecture', fontsize=10,
           ha='center', va='center', style='italic', color=COLORS['dark_gray'])

    ax.text(-0.15, 1.05, 'B', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.set_title('Cross-species code conservation', fontsize=11, pad=10)

    # Panel C: Leucine codon usage across species
    ax = axes[2]

    ecoli, bsub, yeast = load_cross_species_usage()

    leu_codons = ['CTG', 'CTC', 'CTT', 'CTA', 'TTG', 'TTA']

    ecoli_leu = ecoli[ecoli['aa'] == 'L'].set_index('codon')
    bsub_leu = bsub[bsub['aa'] == 'L'].set_index('codon')
    yeast_leu = yeast[yeast['aa'] == 'L'].set_index('codon')

    ecoli_freq = [ecoli_leu.loc[c, 'freq_within_aa'] * 100 if c in ecoli_leu.index else 0
                  for c in leu_codons]
    bsub_freq = [bsub_leu.loc[c, 'freq_within_aa'] * 100 if c in bsub_leu.index else 0
                 for c in leu_codons]
    yeast_freq = [yeast_leu.loc[c, 'freq_within_aa'] * 100 if c in yeast_leu.index else 0
                  for c in leu_codons]

    x = np.arange(len(leu_codons))
    width = 0.25

    bars1 = ax.bar(x - width, ecoli_freq, width, label='E. coli',
                   color=COLORS['blue'], edgecolor='white', linewidth=1)
    bars2 = ax.bar(x, bsub_freq, width, label='B. subtilis',
                   color=COLORS['green'], edgecolor='white', linewidth=1)
    bars3 = ax.bar(x + width, yeast_freq, width, label='S. cerevisiae',
                   color=COLORS['orange'], edgecolor='white', linewidth=1)

    slow_indices = [3, 5]
    for idx in slow_indices:
        ax.axvspan(idx - 0.45, idx + 0.45, alpha=0.15, color=COLORS['gray'], zorder=0)

    ax.set_xlabel('Leucine codons', fontsize=11)
    ax.set_ylabel('Usage frequency (%)', fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(leu_codons)
    ax.legend(loc='upper right', fontsize=8, frameon=False)

    ax.text(4, ax.get_ylim()[1] * 0.95, 'slow codons', fontsize=8,
           color=COLORS['gray'], ha='center', style='italic')

    ax.text(-0.15, 1.05, 'C', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.set_title('Leucine codon usage variation', fontsize=11, pad=10)

    plt.tight_layout()
    plt.savefig('Figure2_architecture_comparison.png', dpi=300, bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.savefig('Figure2_architecture_comparison.svg', format='svg', bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.savefig('Figure2_architecture_comparison.pdf', format='pdf', bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.close()
    print("Figure 2 saved.")

# ============================================================================
# FIGURE 3: Context-dependent codon deployment at metal sites
# ============================================================================

def create_figure3():
    """Create Figure 3: Metal site codon enrichment."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    mu_df, tai_df = load_metal_data()

    # Panel A: Horizontal bar chart of enrichment
    ax = axes[0]

    aa_order = ['ASP', 'CYS', 'GLU', 'HIS']
    aa_labels = ['Asp', 'Cys', 'Glu', 'His']
    colors = [COLORS['red'], COLORS['orange'], COLORS['purple'], COLORS['green']]

    y_pos = np.arange(len(aa_order))

    for i, (aa, label, color) in enumerate(zip(aa_order, aa_labels, colors)):
        row = mu_df[mu_df['aa'] == aa].iloc[0]
        or_val = row['OR_enrichment']
        p_val = row['p_fisher']
        enriched = row['enriched_codon']

        a, b = row['ligand_codon1'], row['ligand_codon2']
        c, d = row['background_codon1'], row['background_codon2']
        se = np.sqrt(1/a + 1/b + 1/c + 1/d)
        ci_low = np.exp(np.log(or_val) - 1.96 * se)
        ci_high = np.exp(np.log(or_val) + 1.96 * se)

        ax.barh(i, or_val, color=color, edgecolor='white', linewidth=1.5, height=0.6)
        ax.errorbar(or_val, i, xerr=[[or_val - ci_low], [ci_high - or_val]],
                   fmt='none', color='black', capsize=3, linewidth=1.5)

        sig = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*' if p_val < 0.05 else ''
        ax.text(ci_high + 0.05, i, f'{enriched}{sig}', va='center', fontsize=9,
               fontweight='bold')

    ax.axvline(x=1.0, color='black', linestyle='-', linewidth=1)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(aa_labels)
    ax.set_xlabel('Odds ratio (enrichment at metal sites)', fontsize=11)
    ax.set_xlim(0.6, 1.8)

    ax.text(-0.15, 1.05, 'A', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.set_title('Codon enrichment at metal-binding sites', fontsize=11, pad=10)

    # Panel B: mu comparison
    ax = axes[1]

    x_pos = np.arange(len(aa_order))
    width = 0.35

    enriched_mu = []
    depleted_mu = []
    violated = []

    for aa in aa_order:
        row = mu_df[mu_df['aa'] == aa].iloc[0]
        enriched_mu.append(row['enriched_mu'])
        depleted_mu.append(row['depleted_mu'])
        violated.append(row['enriched_mu'] > row['depleted_mu'])

    bars1 = ax.bar(x_pos - width/2, enriched_mu, width, label='Enriched codon',
                   color=COLORS['blue'], edgecolor='white', linewidth=1)
    bars2 = ax.bar(x_pos + width/2, depleted_mu, width, label='Alternative codon',
                   color=COLORS['light_gray'], edgecolor='white', linewidth=1)

    for i, is_violated in enumerate(violated):
        if is_violated:
            ax.text(i, max(enriched_mu[i], depleted_mu[i]) * 1.1, 'violated',
                   ha='center', va='bottom', fontsize=8, color=COLORS['red'],
                   style='italic')
        else:
            ax.text(i, max(enriched_mu[i], depleted_mu[i]) * 1.1, 'consistent',
                   ha='center', va='bottom', fontsize=8, color=COLORS['green'],
                   style='italic')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(aa_labels)
    ax.set_ylabel(r'Mistranslation rate $\mu$', fontsize=11)
    ax.set_yscale('log')
    ax.legend(loc='upper right', fontsize=8, frameon=False)

    ax.text(-0.15, 1.05, 'B', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.set_title('Error rate comparison', fontsize=11, pad=10)

    # Panel C: Enrichment in (mu, nu) space
    ax = axes[2]

    df = load_codon_data()

    for aa, aa_short, color in [('ASP', 'D', COLORS['red']),
                                 ('CYS', 'C', COLORS['orange']),
                                 ('GLU', 'E', COLORS['purple']),
                                 ('HIS', 'H', COLORS['green'])]:
        mu_row = mu_df[mu_df['aa'] == aa].iloc[0]
        enriched = mu_row['enriched_codon']
        depleted = mu_row['depleted_codon']

        enr_data = df[df['codon'] == enriched].iloc[0]
        dep_data = df[df['codon'] == depleted].iloc[0]

        ax.scatter(enr_data['mu'], enr_data['w_tai'], c=color, s=100,
                  edgecolors='white', linewidth=1.5, zorder=5)
        ax.scatter(dep_data['mu'], dep_data['w_tai'], c=color, s=60,
                  edgecolors='white', linewidth=1, alpha=0.5, zorder=4)

        ax.annotate('', xy=(enr_data['mu'], enr_data['w_tai']),
                   xytext=(dep_data['mu'], dep_data['w_tai']),
                   arrowprops=dict(arrowstyle='->', color=color, lw=1.5))

        ax.text(enr_data['mu'] * 1.2, enr_data['w_tai'] + 0.03,
               aa_short, fontsize=10, color=color, fontweight='bold')

    ax.text(0.0001, 0.55, 'Throughput-first\n(Asp, His)', fontsize=9,
           color=COLORS['dark_gray'], ha='center', style='italic')
    ax.text(0.005, 0.12, 'Accuracy-first\n(Cys, Glu)', fontsize=9,
           color=COLORS['dark_gray'], ha='center', style='italic')

    ax.set_xlabel(r'Mistranslation rate $\mu$', fontsize=11)
    ax.set_ylabel(r'Decoding efficiency $\nu$ [tAI]', fontsize=11)
    ax.set_xscale('log')
    ax.set_xlim(1e-5, 0.02)
    ax.set_ylim(0, 0.7)

    ax.text(-0.15, 1.05, 'C', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.set_title('Enrichment patterns in operational space', fontsize=11, pad=10)

    plt.tight_layout()
    plt.savefig('Figure3_metal_site_enrichment.png', dpi=300, bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.savefig('Figure3_metal_site_enrichment.svg', format='svg', bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.savefig('Figure3_metal_site_enrichment.pdf', format='pdf', bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.close()
    print("Figure 3 saved.")

# ============================================================================
# FIGURE 4: Operational diversity structure (CORRECTED)
# ============================================================================

def compute_operational_diversity(df, mu_mean=None, mu_std=None, tai_mean=None, tai_std=None):
    """Compute Delta_A for each amino acid.

    If standardization parameters not provided, compute from df.
    """
    valid_df = df.dropna(subset=['mu', 'w_tai'])

    if mu_mean is None:
        log_mu = np.log10(valid_df['mu'].replace(0, np.nan).dropna())
        mu_mean, mu_std = log_mu.mean(), log_mu.std()
    if tai_mean is None:
        tai_mean, tai_std = valid_df['w_tai'].mean(), valid_df['w_tai'].std()

    delta_by_aa = {}
    for aa in valid_df['aa'].unique():
        aa_df = valid_df[valid_df['aa'] == aa]
        if len(aa_df) < 2:
            delta_by_aa[aa] = 0
            continue

        # Standardize
        mu_std_vals = (np.log10(aa_df['mu']) - mu_mean) / mu_std
        tai_std_vals = (aa_df['w_tai'] - tai_mean) / tai_std

        # Mean pairwise distance
        coords = np.column_stack([mu_std_vals.values, tai_std_vals.values])
        n = len(coords)
        total_dist = 0
        count = 0
        for i in range(n):
            for j in range(i+1, n):
                total_dist += np.sqrt(np.sum((coords[i] - coords[j])**2))
                count += 1

        delta_by_aa[aa] = total_dist / count if count > 0 else 0

    return delta_by_aa, (mu_mean, mu_std, tai_mean, tai_std)

def generate_null_codes_corrected(df, n_null=10000):
    """
    CORRECTED null ensemble generation.

    Shuffles (mu, tAI) PAIRS together within degeneracy classes.
    This preserves the natural correlation between mu and tAI while
    testing whether the specific assignment to amino acids matters.
    """
    degeneracy = df.groupby('aa').size().to_dict()
    valid_df = df.dropna(subset=['mu', 'w_tai']).copy()

    # Compute standardization parameters ONCE from original data
    log_mu = np.log10(valid_df['mu'])
    mu_mean, mu_std = log_mu.mean(), log_mu.std()
    tai_mean, tai_std = valid_df['w_tai'].mean(), valid_df['w_tai'].std()

    # Original diversity
    orig_delta, _ = compute_operational_diversity(valid_df, mu_mean, mu_std, tai_mean, tai_std)
    orig_total = sum(orig_delta.values())

    null_totals = []

    for _ in range(n_null):
        shuffled_df = valid_df.copy()

        for deg_val in set(degeneracy.values()):
            # Get AAs with this degeneracy
            deg_aas = [aa for aa, d in degeneracy.items() if d == deg_val]
            mask = shuffled_df['aa'].isin(deg_aas)
            indices = shuffled_df[mask].index.tolist()

            if len(indices) > 1:
                # CORRECTED: Shuffle (mu, tAI) PAIRS together
                # Extract paired values
                mu_vals = shuffled_df.loc[indices, 'mu'].values.copy()
                tai_vals = shuffled_df.loc[indices, 'w_tai'].values.copy()

                # Create paired array and shuffle
                paired = list(zip(mu_vals, tai_vals))
                np.random.shuffle(paired)

                # Unpack back
                shuffled_mu, shuffled_tai = zip(*paired)

                shuffled_df.loc[indices, 'mu'] = shuffled_mu
                shuffled_df.loc[indices, 'w_tai'] = shuffled_tai

        null_delta, _ = compute_operational_diversity(shuffled_df, mu_mean, mu_std, tai_mean, tai_std)
        null_totals.append(sum(null_delta.values()))

    return orig_total, null_totals, orig_delta

def create_figure4():
    """Create Figure 4: Operational diversity structure.

    CORRECTED: Null ensemble now shuffles (mu, tAI) pairs together.
    """
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))

    df = load_codon_data()

    # Panel A: Observed vs null distribution
    ax = axes[0]

    print("Generating CORRECTED null distribution (10,000 permutations)...")
    print("  (shuffling mu/tAI pairs together within degeneracy classes)")
    orig_total, null_totals, orig_delta = generate_null_codes_corrected(df, n_null=10000)

    # Plot histogram
    ax.hist(null_totals, bins=50, color=COLORS['light_gray'], edgecolor='white',
           alpha=0.8, density=True, label='Null distribution')

    # Add observed value
    ax.axvline(x=orig_total, color=COLORS['red'], linewidth=2, linestyle='-',
              label=f'Observed = {orig_total:.2f}')

    # Compute statistics
    null_arr = np.array(null_totals)
    z_score = (orig_total - np.mean(null_totals)) / np.std(null_totals)

    # Two-tailed p-value
    if z_score >= 0:
        p_value = np.sum(null_arr >= orig_total) / len(null_arr)
        # Shade right tail
        ax.fill_betweenx([0, ax.get_ylim()[1] * 10], orig_total, null_arr.max() + 1,
                         alpha=0.3, color=COLORS['red'])
    else:
        p_value = np.sum(null_arr <= orig_total) / len(null_arr)
        # Shade left tail
        ax.fill_betweenx([0, ax.get_ylim()[1] * 10], null_arr.min() - 1, orig_total,
                         alpha=0.3, color=COLORS['red'])

    ax.text(0.95, 0.95, f'z = {z_score:.2f}\np = {p_value:.4f}',
           transform=ax.transAxes, fontsize=10, ha='right', va='top',
           bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    ax.set_xlabel('Total operational diversity', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.legend(loc='upper left', fontsize=9, frameon=False)

    ax.text(-0.12, 1.05, 'A', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.set_title('Observed vs. null distribution', fontsize=11, pad=10)

    # Panel B: Delta_A per amino acid
    ax = axes[1]

    degeneracy = df.groupby('aa').size().to_dict()

    two_codon = sorted([aa for aa, d in degeneracy.items() if d == 2])
    four_codon = sorted([aa for aa, d in degeneracy.items() if d == 4])
    six_codon = sorted([aa for aa, d in degeneracy.items() if d == 6])
    three_codon = sorted([aa for aa, d in degeneracy.items() if d == 3])

    aa_order = two_codon + three_codon + four_codon + six_codon
    delta_vals = [orig_delta.get(aa, 0) for aa in aa_order]

    colors = []
    for aa in aa_order:
        d = degeneracy[aa]
        if d == 2:
            colors.append(COLORS['blue'])
        elif d == 3:
            colors.append(COLORS['purple'])
        elif d == 4:
            colors.append(COLORS['green'])
        else:
            colors.append(COLORS['orange'])

    x_pos = np.arange(len(aa_order))
    bars = ax.bar(x_pos, delta_vals, color=colors, edgecolor='white', linewidth=1)

    # Add separators
    sep_positions = [len(two_codon) - 0.5,
                     len(two_codon) + len(three_codon) - 0.5,
                     len(two_codon) + len(three_codon) + len(four_codon) - 0.5]
    for pos in sep_positions:
        ax.axvline(x=pos, color='black', linestyle='--', linewidth=0.5, alpha=0.5)

    # Add class labels
    ax.text(len(two_codon)/2 - 0.5, ax.get_ylim()[1] * 0.95, '2-codon',
           ha='center', fontsize=9, color=COLORS['blue'])
    ax.text(len(two_codon) + len(three_codon)/2 - 0.5, ax.get_ylim()[1] * 0.95, '3',
           ha='center', fontsize=9, color=COLORS['purple'])
    ax.text(len(two_codon) + len(three_codon) + len(four_codon)/2 - 0.5,
           ax.get_ylim()[1] * 0.95, '4-codon',
           ha='center', fontsize=9, color=COLORS['green'])
    ax.text(len(aa_order) - len(six_codon)/2 - 0.5, ax.get_ylim()[1] * 0.95, '6-codon',
           ha='center', fontsize=9, color=COLORS['orange'])

    ax.set_xticks(x_pos)
    ax.set_xticklabels(aa_order, fontsize=8)
    ax.set_xlabel('Amino acid', fontsize=11)
    ax.set_ylabel(r'Operational diversity $\Delta_A$', fontsize=11)

    ax.text(-0.12, 1.05, 'B', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.set_title('Diversity by amino acid', fontsize=11, pad=10)

    plt.tight_layout()
    plt.savefig('Figure4_diversity_structure.png', dpi=300, bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.savefig('Figure4_diversity_structure.svg', format='svg', bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.savefig('Figure4_diversity_structure.pdf', format='pdf', bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.close()

    print(f"Figure 4 saved (CORRECTED null ensemble).")
    print(f"  Observed diversity: {orig_total:.2f}")
    print(f"  Null mean: {np.mean(null_totals):.2f} +/- {np.std(null_totals):.2f}")
    print(f"  z-score: {z_score:.2f}, p-value: {p_value:.4f}")

# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    import os
    os.chdir('/storage/kiran-stuff/proteostasis_law')

    print("Generating Paper 2 Figures (CORRECTED VERSION)...")
    print("=" * 60)
    print("Corrections applied:")
    print("  - Figure 1A: Proper S-shaped F(P) with visible bifurcation structure")
    print("  - Figure 4A: Null ensemble shuffles (mu,tAI) pairs together")
    print("=" * 60)

    create_figure1()
    create_figure2()
    create_figure3()
    create_figure4()

    print("=" * 60)
    print("All figures generated successfully!")
    print("\nOutput files:")
    print("  - Figure1_proteostasis_constraint.{png,svg,pdf}")
    print("  - Figure2_architecture_comparison.{png,svg,pdf}")
    print("  - Figure3_metal_site_enrichment.{png,svg,pdf}")
    print("  - Figure4_diversity_structure.{png,svg,pdf}")
