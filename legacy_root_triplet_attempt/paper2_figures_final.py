#!/usr/bin/env python3
"""
Generate Paper 2 Figures: The Triplet Code as Minimal Control Architecture
for Translation Under Proteostasis Constraint

FINAL CORRECTED VERSION:
- Figure 1A: Hill-shaped F(P) curve (rises then falls due to system overload)
- Figure 4: Reframed to show wobble-driven patterns in 2-codon families

Data sources:
- errors/codon_modes_ecoli.tsv (mu from Landerer et al. 2024, tAI from stAIcalc)
- metals/metal_codon_bias_with_mu.tsv (metal site enrichment)
- cross_species/data/*/global_codon_usage.tsv (codon usage)
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

AA_COLORS = {
    'L': '#2E86AB', 'D': '#E94F37', 'H': '#2A9D8F',
    'C': '#F77F00', 'E': '#9B5DE5',
}

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
    df = pd.read_csv('errors/codon_modes_ecoli.tsv', sep='\t')
    df['mu'] = pd.to_numeric(df['mu'], errors='coerce')
    df['w_tai'] = pd.to_numeric(df['w_tai'], errors='coerce')
    return df

def load_metal_data():
    mu_df = pd.read_csv('metals/metal_codon_bias_with_mu.tsv', sep='\t')
    tai_df = pd.read_csv('metals/metal_codon_bias_with_tai.tsv', sep='\t')
    return mu_df, tai_df

def load_cross_species_usage():
    ecoli = pd.read_csv('cross_species/data/ecoli/global_codon_usage.tsv', sep='\t')
    bsub = pd.read_csv('cross_species/data/bsub/global_codon_usage.tsv', sep='\t')
    yeast = pd.read_csv('cross_species/data/yeast/global_codon_usage.tsv', sep='\t')
    return ecoli, bsub, yeast

def compute_synonymy_shielding():
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
    """
    Figure 1: Proteostasis constraint and codon operational space.

    Panel A: CORRECTED with hill-shaped F(P) curve.
    - F(P) rises at low P (chaperones handle substrate)
    - F(P) peaks at intermediate P (saturation)
    - F(P) FALLS at high P (system overload, aggregation poisoning)
    """
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel A: Hill-shaped F(P) bifurcation
    ax = axes[0]

    # Parameters for logistic-style outflux: F(P) = a * P * (1 - P/P_max)
    a = 0.8
    P_max = 12.0

    P = np.linspace(0.001, 12, 500)

    def F(P):
        """Hill-shaped outflux: rises, peaks, then falls."""
        return a * P * (1 - P / P_max)

    F_vals = F(P)

    # Peak of F(P) at P = P_max/2
    P_peak = P_max / 2
    F_peak = F(P_peak)

    # Define J levels
    J_crit = F_peak          # Exactly at peak (saddle-node bifurcation)
    J_low = F_peak * 0.65    # Well below - two intersections
    J_high = F_peak * 1.15   # Above peak - no intersection

    # Plot F(P) curve
    ax.plot(P, F_vals, 'k-', lw=2.5, label=r'Outflux $F(P)$')

    # Shade viable region (under curve, left of peak)
    ax.fill_between(P, 0, F_vals, where=(P <= P_peak), alpha=0.1, color='gray')
    # Shade collapse region (right of peak)
    ax.fill_between(P, 0, F_vals, where=(P > P_peak), alpha=0.1, color=COLORS['red'])

    # Plot J_in lines
    ax.axhline(J_low, color=COLORS['green'], linestyle='--', linewidth=1.5)
    ax.axhline(J_crit, color=COLORS['orange'], linestyle='-', linewidth=1.5)
    ax.axhline(J_high, color=COLORS['red'], linestyle='--', linewidth=1.5)

    # Find intersections for J_low: a * P * (1 - P/P_max) = J_low
    # -a/P_max * P^2 + a * P - J_low = 0
    # P = (a ± sqrt(a^2 - 4*(a/P_max)*J_low)) / (2*a/P_max)
    discriminant = a**2 - 4 * (a/P_max) * J_low
    if discriminant > 0:
        P_stable = (a - np.sqrt(discriminant)) / (2 * a / P_max)
        P_unstable = (a + np.sqrt(discriminant)) / (2 * a / P_max)

        # Stable point (filled circle)
        ax.plot(P_stable, J_low, 'o', color=COLORS['green'], markersize=12,
               markeredgecolor='black', markeredgewidth=1.5, zorder=5)
        # Unstable point (open circle)
        ax.plot(P_unstable, J_low, 'o', color='white', markersize=12,
               markeredgecolor=COLORS['green'], markeredgewidth=2.5, zorder=5)

        # Annotations
        ax.annotate('stable', xy=(P_stable, J_low),
                   xytext=(P_stable - 1.0, J_low + 0.25),
                   fontsize=9, ha='center', color=COLORS['dark_gray'])
        ax.annotate('unstable', xy=(P_unstable, J_low),
                   xytext=(P_unstable + 1.2, J_low + 0.25),
                   fontsize=9, ha='center', color=COLORS['dark_gray'])

    # Saddle-node at peak
    ax.plot(P_peak, J_crit, 's', color=COLORS['orange'], markersize=12,
           markeredgecolor='black', markeredgewidth=1.5, zorder=5)
    ax.annotate('saddle-node', xy=(P_peak, J_crit),
               xytext=(P_peak + 1.8, J_crit + 0.25),
               fontsize=9, ha='left', color=COLORS['dark_gray'],
               arrowprops=dict(arrowstyle='->', color=COLORS['gray'], lw=1))

    # Labels for J lines
    ax.text(11.5, J_low - 0.12, r'$J_{in} < J_{crit}$', fontsize=9,
           color=COLORS['green'], ha='right', va='top')
    ax.text(11.5, J_crit + 0.08, r'$J_{in} = J_{crit}$', fontsize=9,
           color=COLORS['orange'], ha='right', va='bottom')
    ax.text(11.5, J_high + 0.08, r'$J_{in} > J_{crit}$', fontsize=9,
           color=COLORS['red'], ha='right', va='bottom')

    # Collapse annotation
    ax.text(9, F_peak * 0.85, 'collapse\nregion', fontsize=9,
           color=COLORS['red'], ha='center', style='italic')

    # Equation box
    ax.text(0.03, 0.97, r'$\frac{dP}{dt} = J_{in} - F(P)$',
           transform=ax.transAxes, fontsize=10, verticalalignment='top',
           bbox=dict(boxstyle='round', facecolor='white', edgecolor='gray', alpha=0.9))

    ax.set_xlabel('Misfold burden $P$', fontsize=11)
    ax.set_ylabel('Flux', fontsize=11)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, F_peak * 1.4)

    ax.text(-0.12, 1.05, 'A', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    ax.set_title('Saddle-node bifurcation', fontsize=11, pad=10)

    # Panel B: Codons in (mu, nu) space
    ax = axes[1]
    df = load_codon_data()

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
                    ax.plot([coords[i,0], coords[j,0]], [coords[i,1], coords[j,1]],
                           color=AA_COLORS[aa], alpha=0.3, linewidth=1, zorder=1)

    # Labels for key codons
    label_codons = {'CTG': 'L', 'CTA': 'L', 'GAC': 'D', 'GAT': 'D'}
    for _, row in df.iterrows():
        if row['codon'] in label_codons:
            mu, tai = row['mu'], row['w_tai']
            if not pd.isna(mu) and not pd.isna(tai):
                offset = (5, 5) if row['codon'] in ['CTG', 'GAC'] else (-25, -10)
                ax.annotate(row['codon'], (mu, tai), xytext=offset,
                           textcoords='offset points', fontsize=8,
                           color=AA_COLORS[row['aa']], fontweight='bold')

    # Leu spread arrow
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

    legend_elements = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor=AA_COLORS['L'], markersize=8, label='Leu'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=AA_COLORS['D'], markersize=8, label='Asp'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=AA_COLORS['H'], markersize=8, label='His'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=AA_COLORS['C'], markersize=8, label='Cys'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=AA_COLORS['E'], markersize=8, label='Glu'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor=COLORS['light_gray'], markersize=6, label='Other'),
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, frameon=False, ncol=2)
    ax.text(-0.12, 1.05, 'B', transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')
    ax.set_title('Synonymous codons in operational space', fontsize=11, pad=10)

    # Panel C: Synonymy shielding
    ax = axes[2]
    s_values = compute_synonymy_shielding()
    positions = [1, 2, 3]
    colors = [COLORS['gray'], COLORS['gray'], COLORS['blue']]

    bars = ax.bar(positions, s_values, color=colors, edgecolor='white', linewidth=1.5)
    ax.axhline(y=s_values[2], color=COLORS['blue'], linestyle='--', alpha=0.5, linewidth=1)
    ax.text(3.5, s_values[2] + 0.02, f'$S_3 \\approx {s_values[2]:.2f}$',
           fontsize=10, color=COLORS['blue'], va='bottom')

    ax.set_xlabel('Codon position', fontsize=11)
    ax.set_ylabel('Synonymy fraction $S_i$', fontsize=11)
    ax.set_xticks([1, 2, 3])
    ax.set_xticklabels(['1st', '2nd', '3rd\n(wobble)'])
    ax.set_ylim(0, 0.85)

    for i, (bar, val) in enumerate(zip(bars, s_values)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.02,
               f'{val:.2f}', ha='center', va='bottom', fontsize=9)

    ax.text(-0.12, 1.05, 'C', transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')
    ax.set_title('Synonymy shielding by position', fontsize=11, pad=10)

    plt.tight_layout()
    plt.savefig('Figure1_proteostasis_constraint.png', dpi=300, bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.savefig('Figure1_proteostasis_constraint.svg', format='svg', bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.savefig('Figure1_proteostasis_constraint.pdf', format='pdf', bbox_inches='tight',
               facecolor='white', edgecolor='none')
    plt.close()
    print("Figure 1 saved (hill-shaped F(P) curve).")

# ============================================================================
# FIGURE 2: Architectural comparison
# ============================================================================

def create_figure2():
    """Figure 2: Architecture comparison and cross-species conservation."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    # Panel A: Vocabulary capacity
    ax = axes[0]
    code_types = ['Doublet\n(n=2)', 'Triplet\n(n=3)', 'Quadruplet\n(n=4)']
    sense_codons = [15, 61, 255]
    synonymy_display = [0, (61-20)/20, (255-20)/20]

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

    ax.axhline(y=20, color=COLORS['red'], linestyle='--', linewidth=1.5, label='20 AA requirement')
    ax.text(0, 25, '✗', fontsize=20, color=COLORS['red'], ha='center', fontweight='bold')
    ax.text(1, 70, '✓', fontsize=20, color=COLORS['green'], ha='center', fontweight='bold')
    ax.text(2, 260, 'excess', fontsize=9, color=COLORS['gray'], ha='center', style='italic')

    ax.set_xticks(x)
    ax.set_xticklabels(code_types)
    ax.set_xlabel('Code architecture', fontsize=11)
    lines1, labels1 = ax.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax.legend(lines1 + lines2, labels1 + labels2, loc='upper left', fontsize=8, frameon=False)
    ax.text(-0.15, 1.05, 'A', transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')
    ax.set_title('Vocabulary capacity by code length', fontsize=11, pad=10)

    # Panel B: Universal triplet architecture
    ax = axes[1]
    species = ['E. coli', 'B. subtilis', 'S. cerevisiae']
    synonymy_values = [61/20, 61/20, 61/20]
    bars = ax.bar(species, synonymy_values, color=[COLORS['blue'], COLORS['green'], COLORS['orange']],
                  edgecolor='white', linewidth=1.5)
    ax.axhline(y=1.0, color=COLORS['red'], linestyle='--', linewidth=1.5, label='Minimum for control')
    ax.set_ylabel('Synonymy per amino acid', fontsize=11)
    ax.set_ylim(0, 4)
    ax.legend(loc='upper right', fontsize=8, frameon=False)
    ax.text(1, 3.5, 'Universal triplet\narchitecture', fontsize=10, ha='center', va='center',
           style='italic', color=COLORS['dark_gray'])
    ax.text(-0.15, 1.05, 'B', transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')
    ax.set_title('Cross-species code conservation', fontsize=11, pad=10)

    # Panel C: Leucine codon usage
    ax = axes[2]
    ecoli, bsub, yeast = load_cross_species_usage()
    leu_codons = ['CTG', 'CTC', 'CTT', 'CTA', 'TTG', 'TTA']

    ecoli_leu = ecoli[ecoli['aa'] == 'L'].set_index('codon')
    bsub_leu = bsub[bsub['aa'] == 'L'].set_index('codon')
    yeast_leu = yeast[yeast['aa'] == 'L'].set_index('codon')

    ecoli_freq = [ecoli_leu.loc[c, 'freq_within_aa'] * 100 if c in ecoli_leu.index else 0 for c in leu_codons]
    bsub_freq = [bsub_leu.loc[c, 'freq_within_aa'] * 100 if c in bsub_leu.index else 0 for c in leu_codons]
    yeast_freq = [yeast_leu.loc[c, 'freq_within_aa'] * 100 if c in yeast_leu.index else 0 for c in leu_codons]

    x = np.arange(len(leu_codons))
    width = 0.25
    ax.bar(x - width, ecoli_freq, width, label='E. coli', color=COLORS['blue'], edgecolor='white')
    ax.bar(x, bsub_freq, width, label='B. subtilis', color=COLORS['green'], edgecolor='white')
    ax.bar(x + width, yeast_freq, width, label='S. cerevisiae', color=COLORS['orange'], edgecolor='white')

    for idx in [3, 5]:
        ax.axvspan(idx - 0.45, idx + 0.45, alpha=0.15, color=COLORS['gray'], zorder=0)

    ax.set_xlabel('Leucine codons', fontsize=11)
    ax.set_ylabel('Usage frequency (%)', fontsize=11)
    ax.set_xticks(x)
    ax.set_xticklabels(leu_codons)
    ax.legend(loc='upper right', fontsize=8, frameon=False)
    ax.text(4, ax.get_ylim()[1] * 0.95, 'slow codons', fontsize=8, color=COLORS['gray'], ha='center', style='italic')
    ax.text(-0.15, 1.05, 'C', transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')
    ax.set_title('Leucine codon usage variation', fontsize=11, pad=10)

    plt.tight_layout()
    plt.savefig('Figure2_architecture_comparison.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('Figure2_architecture_comparison.svg', format='svg', bbox_inches='tight', facecolor='white')
    plt.savefig('Figure2_architecture_comparison.pdf', format='pdf', bbox_inches='tight', facecolor='white')
    plt.close()
    print("Figure 2 saved.")

# ============================================================================
# FIGURE 3: Metal site enrichment
# ============================================================================

def create_figure3():
    """Figure 3: Metal site codon enrichment."""
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    mu_df, tai_df = load_metal_data()

    # Panel A: Enrichment odds ratios
    ax = axes[0]
    aa_order = ['ASP', 'CYS', 'GLU', 'HIS']
    aa_labels = ['Asp', 'Cys', 'Glu', 'His']
    colors = [COLORS['red'], COLORS['orange'], COLORS['purple'], COLORS['green']]

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
        ax.text(ci_high + 0.05, i, f'{enriched}{sig}', va='center', fontsize=9, fontweight='bold')

    ax.axvline(x=1.0, color='black', linestyle='-', linewidth=1)
    ax.set_yticks(range(len(aa_labels)))
    ax.set_yticklabels(aa_labels)
    ax.set_xlabel('Odds ratio (enrichment at metal sites)', fontsize=11)
    ax.set_xlim(0.6, 1.8)
    ax.text(-0.15, 1.05, 'A', transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')
    ax.set_title('Codon enrichment at metal-binding sites', fontsize=11, pad=10)

    # Panel B: mu comparison with strategy labels
    ax = axes[1]
    x_pos = np.arange(len(aa_order))
    width = 0.35
    enriched_mu, depleted_mu, strategies = [], [], []

    for aa in aa_order:
        mu_row = mu_df[mu_df['aa'] == aa].iloc[0]
        tai_row = tai_df[tai_df['aa'] == aa].iloc[0]
        enriched_mu.append(mu_row['enriched_mu'])
        depleted_mu.append(mu_row['depleted_mu'])

        # Determine strategy based on what enriched codon optimizes
        higher_mu = mu_row['enriched_mu'] > mu_row['depleted_mu']
        higher_tai = tai_row['enriched_tai'] > tai_row['depleted_tai']

        if higher_mu and higher_tai:
            strategies.append(('speed', COLORS['orange']))  # Throughput-first
        elif not higher_mu and not higher_tai:
            strategies.append(('accuracy', COLORS['blue']))  # Accuracy-first
        elif not higher_mu and higher_tai:
            strategies.append(('optimal', COLORS['green']))  # Better on both
        else:
            strategies.append(('trade-off', COLORS['gray']))

    ax.bar(x_pos - width/2, enriched_mu, width, label='Enriched codon', color=COLORS['blue'], edgecolor='white')
    ax.bar(x_pos + width/2, depleted_mu, width, label='Alternative codon', color=COLORS['light_gray'], edgecolor='white')

    for i, (strategy, color) in enumerate(strategies):
        ax.text(i, max(enriched_mu[i], depleted_mu[i]) * 1.15, strategy,
               ha='center', va='bottom', fontsize=8, color=color, style='italic')

    ax.set_xticks(x_pos)
    ax.set_xticklabels(aa_labels)
    ax.set_ylabel(r'Mistranslation rate $\mu$', fontsize=11)
    ax.set_yscale('log')
    ax.legend(loc='upper right', fontsize=8, frameon=False)
    ax.text(-0.15, 1.05, 'B', transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')
    ax.set_title('Error rate comparison', fontsize=11, pad=10)

    # Panel C: Enrichment in (mu, nu) space
    ax = axes[2]
    df = load_codon_data()

    for aa, aa_short, color in [('ASP', 'D', COLORS['red']), ('CYS', 'C', COLORS['orange']),
                                 ('GLU', 'E', COLORS['purple']), ('HIS', 'H', COLORS['green'])]:
        mu_row = mu_df[mu_df['aa'] == aa].iloc[0]
        enriched = mu_row['enriched_codon']
        depleted = mu_row['depleted_codon']
        enr_data = df[df['codon'] == enriched].iloc[0]
        dep_data = df[df['codon'] == depleted].iloc[0]

        ax.scatter(enr_data['mu'], enr_data['w_tai'], c=color, s=100, edgecolors='white', linewidth=1.5, zorder=5)
        ax.scatter(dep_data['mu'], dep_data['w_tai'], c=color, s=60, edgecolors='white', linewidth=1, alpha=0.5, zorder=4)
        ax.annotate('', xy=(enr_data['mu'], enr_data['w_tai']),
                   xytext=(dep_data['mu'], dep_data['w_tai']),
                   arrowprops=dict(arrowstyle='->', color=color, lw=1.5))
        ax.text(enr_data['mu'] * 1.2, enr_data['w_tai'] + 0.03, aa_short, fontsize=10, color=color, fontweight='bold')

    ax.text(0.0001, 0.55, 'Throughput-first\n(Asp, His)', fontsize=9, color=COLORS['dark_gray'], ha='center', style='italic')
    ax.text(0.005, 0.12, 'Accuracy-first\n(Cys, Glu)', fontsize=9, color=COLORS['dark_gray'], ha='center', style='italic')

    ax.set_xlabel(r'Mistranslation rate $\mu$', fontsize=11)
    ax.set_ylabel(r'Decoding efficiency $\nu$ [tAI]', fontsize=11)
    ax.set_xscale('log')
    ax.set_xlim(1e-5, 0.02)
    ax.set_ylim(0, 0.7)
    ax.text(-0.15, 1.05, 'C', transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')
    ax.set_title('Enrichment patterns in operational space', fontsize=11, pad=10)

    plt.tight_layout()
    plt.savefig('Figure3_metal_site_enrichment.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('Figure3_metal_site_enrichment.svg', format='svg', bbox_inches='tight', facecolor='white')
    plt.savefig('Figure3_metal_site_enrichment.pdf', format='pdf', bbox_inches='tight', facecolor='white')
    plt.close()
    print("Figure 3 saved.")

# ============================================================================
# FIGURE 4: Operational diversity - CORRECTED NULL CONSTRUCTION
# ============================================================================

def create_figure4():
    """
    Figure 4: Operational diversity structure.

    Panel A: Observed vs null distribution (CORRECTED null construction)
    Panel B: 2-codon family spread by wobble type
    Panel C: Delta_A by amino acid and degeneracy class

    CORRECTED NULL CONSTRUCTION:
    - Shuffles codon->AA assignments (which codons belong to which amino acid)
    - Keeps codon properties (mu, tAI) FIXED for each codon
    - Preserves degeneracy structure (9x2-codon, 1x3-codon, 5x4-codon, 3x6-codon)
    """
    import random

    # Set fixed seed for reproducibility
    random.seed(42)
    np.random.seed(42)

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))

    df = load_codon_data()
    valid_df = df.dropna(subset=['mu', 'w_tai']).copy()

    # Build codon->properties mapping (FIXED throughout)
    codon_properties = {}
    for _, row in valid_df.iterrows():
        codon_properties[row['codon']] = (row['mu'], row['w_tai'])

    # Build observed codon->AA mapping
    observed_codon_to_aa = {}
    for _, row in valid_df.iterrows():
        observed_codon_to_aa[row['codon']] = row['aa']

    # Get degeneracy structure
    aa_degeneracy = defaultdict(int)
    for aa in observed_codon_to_aa.values():
        aa_degeneracy[aa] += 1

    aa_list = sorted(aa_degeneracy.keys())
    degeneracy_list = [aa_degeneracy[aa] for aa in aa_list]
    degeneracy = dict(zip(aa_list, degeneracy_list))

    # Compute normalization parameters (fixed)
    all_mu = [codon_properties[c][0] for c in codon_properties]
    all_tai = [codon_properties[c][1] for c in codon_properties]
    log_mu = np.log10(all_mu)
    MU_MEAN, MU_STD = np.mean(log_mu), np.std(log_mu)
    TAI_MEAN, TAI_STD = np.mean(all_tai), np.std(all_tai)

    def compute_diversity_from_mapping(codon_to_aa):
        """Compute total diversity given a codon->AA mapping."""
        aa_to_codons = defaultdict(list)
        for codon, aa in codon_to_aa.items():
            if codon in codon_properties:
                aa_to_codons[aa].append(codon)

        total = 0
        delta_by_aa = {}
        for aa, codons in aa_to_codons.items():
            if len(codons) < 2:
                delta_by_aa[aa] = 0
                continue
            coords = []
            for codon in codons:
                mu, tai = codon_properties[codon]
                mu_norm = (np.log10(mu) - MU_MEAN) / MU_STD
                tai_norm = (tai - TAI_MEAN) / TAI_STD
                coords.append((mu_norm, tai_norm))

            dist_sum = 0
            count = 0
            for i in range(len(coords)):
                for j in range(i+1, len(coords)):
                    dist_sum += np.sqrt((coords[i][0] - coords[j][0])**2 +
                                       (coords[i][1] - coords[j][1])**2)
                    count += 1
            delta_by_aa[aa] = dist_sum / count if count > 0 else 0
            total += delta_by_aa[aa]

        return total, delta_by_aa

    # Compute OBSERVED diversity
    observed_total, observed_delta = compute_diversity_from_mapping(observed_codon_to_aa)

    # Generate null ensemble - CORRECT: shuffle codon->AA assignments
    all_codons = list(codon_properties.keys())
    n_null = 10000
    null_totals = []

    for _ in range(n_null):
        shuffled_codons = all_codons.copy()
        random.shuffle(shuffled_codons)

        null_codon_to_aa = {}
        idx = 0
        for aa, deg in zip(aa_list, degeneracy_list):
            for _ in range(deg):
                if idx < len(shuffled_codons):
                    null_codon_to_aa[shuffled_codons[idx]] = aa
                    idx += 1

        null_total, _ = compute_diversity_from_mapping(null_codon_to_aa)
        null_totals.append(null_total)

    null_arr = np.array(null_totals)
    null_mean = np.mean(null_arr)
    null_std = np.std(null_arr)
    z_score = (observed_total - null_mean) / null_std
    p_value = np.sum(null_arr <= observed_total) / n_null

    # ========== Panel A: Observed vs null distribution ==========
    ax = axes[0]

    ax.hist(null_totals, bins=50, color=COLORS['light_gray'], edgecolor='white',
            alpha=0.8, density=True, label='Null distribution')
    ax.axvline(x=observed_total, color=COLORS['red'], linewidth=2.5,
               label=f'Observed = {observed_total:.1f}')

    ylim = ax.get_ylim()
    if z_score < 0:
        ax.fill_betweenx([0, ylim[1]], null_arr.min() - 1, observed_total,
                         alpha=0.3, color=COLORS['red'])

    ax.text(0.95, 0.95, f'z = {z_score:.2f}\np = {p_value:.4f}',
            transform=ax.transAxes, fontsize=10, ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))

    ax.set_xlabel('Total operational diversity', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.legend(loc='upper left', fontsize=9, frameon=False)
    ax.set_title('Observed vs. null distribution', fontsize=11, pad=10)
    ax.text(-0.12, 1.05, 'A', transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')

    # ========== Panel B: 2-codon family spread by wobble type ==========
    ax = axes[1]

    two_codon_aas = sorted([aa for aa, d in degeneracy.items() if d == 2])

    purine_aas = []
    pyrimidine_aas = []
    aa_wobble_type = {}

    aa_distances = {}
    for aa in two_codon_aas:
        aa_df = valid_df[valid_df['aa'] == aa]
        if len(aa_df) != 2:
            continue

        codons = sorted(aa_df['codon'].tolist())
        pos3 = (codons[0][2], codons[1][2])

        c1 = ((np.log10(aa_df.iloc[0]['mu']) - MU_MEAN) / MU_STD,
              (aa_df.iloc[0]['w_tai'] - TAI_MEAN) / TAI_STD)
        c2 = ((np.log10(aa_df.iloc[1]['mu']) - MU_MEAN) / MU_STD,
              (aa_df.iloc[1]['w_tai'] - TAI_MEAN) / TAI_STD)

        dist = np.sqrt((c1[0] - c2[0])**2 + (c1[1] - c2[1])**2)
        aa_distances[aa] = dist

        if set(pos3) == {'A', 'G'}:
            purine_aas.append(aa)
            aa_wobble_type[aa] = 'purine'
        elif set(pos3) == {'C', 'T'}:
            pyrimidine_aas.append(aa)
            aa_wobble_type[aa] = 'pyrimidine'

    # Sort by distance within each group
    purine_sorted = sorted([(aa, aa_distances[aa]) for aa in purine_aas], key=lambda x: -x[1])
    pyrimidine_sorted = sorted([(aa, aa_distances[aa]) for aa in pyrimidine_aas], key=lambda x: -x[1])

    all_two_codon = purine_sorted + pyrimidine_sorted
    aa_labels = [x[0] for x in all_two_codon]
    spreads = [x[1] for x in all_two_codon]
    bar_colors = [COLORS['blue'] if aa_wobble_type[aa] == 'purine' else COLORS['orange']
                  for aa in aa_labels]

    x_pos = np.arange(len(aa_labels))
    ax.bar(x_pos, spreads, color=bar_colors, edgecolor='white', linewidth=1)

    # Add wobble type annotations
    for i, (aa, spread) in enumerate(all_two_codon):
        label = 'A/G' if aa_wobble_type[aa] == 'purine' else 'C/T'
        ax.text(i, spread + 0.1, label, ha='center', fontsize=7, color=COLORS['dark_gray'])

    purine_mean = np.mean([aa_distances[aa] for aa in purine_aas])
    pyrimidine_mean = np.mean([aa_distances[aa] for aa in pyrimidine_aas])

    ax.axhline(purine_mean, color=COLORS['blue'], linestyle='--', linewidth=1.5, alpha=0.7)
    ax.axhline(pyrimidine_mean, color=COLORS['orange'], linestyle='--', linewidth=1.5, alpha=0.7)

    legend_elements = [
        mpatches.Patch(facecolor=COLORS['blue'], label=f'Purine (A/G): {purine_mean:.2f}'),
        mpatches.Patch(facecolor=COLORS['orange'], label=f'Pyrimidine (C/T): {pyrimidine_mean:.2f}')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, frameon=False)

    ax.set_xticks(x_pos)
    ax.set_xticklabels(aa_labels, fontsize=9)
    ax.set_xlabel('Amino acid (2-codon families)', fontsize=11)
    ax.set_ylabel('Operational spread\n(standardized distance)', fontsize=11)
    ax.set_title('Wobble position determines operational spread', fontsize=11, pad=10)
    ax.text(-0.12, 1.05, 'B', transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')

    # ========== Panel C: Delta_A by degeneracy class ==========
    ax = axes[2]

    two_codon = sorted([aa for aa, d in degeneracy.items() if d == 2])
    three_codon = sorted([aa for aa, d in degeneracy.items() if d == 3])
    four_codon = sorted([aa for aa, d in degeneracy.items() if d == 4])
    six_codon = sorted([aa for aa, d in degeneracy.items() if d == 6])

    aa_order = two_codon + three_codon + four_codon + six_codon
    delta_vals = [observed_delta.get(aa, 0) for aa in aa_order]

    colors = []
    for aa in aa_order:
        d = degeneracy.get(aa, 0)
        if d == 2:
            colors.append(COLORS['blue'])
        elif d == 3:
            colors.append(COLORS['purple'])
        elif d == 4:
            colors.append(COLORS['green'])
        else:
            colors.append(COLORS['orange'])

    x_pos = np.arange(len(aa_order))
    ax.bar(x_pos, delta_vals, color=colors, edgecolor='white', linewidth=1)

    sep_positions = [len(two_codon) - 0.5,
                     len(two_codon) + len(three_codon) - 0.5,
                     len(two_codon) + len(three_codon) + len(four_codon) - 0.5]
    for pos in sep_positions:
        ax.axvline(x=pos, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)

    ymax = ax.get_ylim()[1]
    ax.text(len(two_codon)/2 - 0.5, ymax * 0.95, '2-codon', ha='center', fontsize=9, color=COLORS['blue'])
    ax.text(len(two_codon) + len(three_codon)/2 - 0.5, ymax * 0.95, '3', ha='center', fontsize=9, color=COLORS['purple'])
    ax.text(len(two_codon) + len(three_codon) + len(four_codon)/2 - 0.5, ymax * 0.95, '4-codon',
            ha='center', fontsize=9, color=COLORS['green'])
    ax.text(len(aa_order) - len(six_codon)/2 - 0.5, ymax * 0.95, '6-codon', ha='center', fontsize=9, color=COLORS['orange'])

    ax.set_xticks(x_pos)
    ax.set_xticklabels(aa_order, fontsize=8)
    ax.set_xlabel('Amino acid', fontsize=11)
    ax.set_ylabel(r'Mean pairwise distance $\Delta_A$', fontsize=11)
    ax.set_title('Operational diversity by amino acid', fontsize=11, pad=10)
    ax.text(-0.12, 1.05, 'C', transform=ax.transAxes, fontsize=14, fontweight='bold', va='top')

    plt.tight_layout()
    plt.savefig('Figure4_diversity_structure.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig('Figure4_diversity_structure.svg', format='svg', bbox_inches='tight', facecolor='white')
    plt.savefig('Figure4_diversity_structure.pdf', format='pdf', bbox_inches='tight', facecolor='white')
    plt.close()

    print(f"Figure 4 saved (CORRECTED null: z={z_score:.2f}, p={p_value:.4f})")

# ============================================================================
# MAIN
# ============================================================================

if __name__ == "__main__":
    import os
    os.chdir('/storage/kiran-stuff/proteostasis_law')

    print("Generating Paper 2 Figures (FINAL VERSION)...")
    print("=" * 60)
    print("Key design decisions:")
    print("  - Figure 1A: Hill-shaped F(P) - rises, peaks, falls")
    print("  - Figure 4: CORRECTED null construction (fixed seed=42)")
    print("              Shuffles codon->AA, keeps (mu,tAI) fixed")
    print("=" * 60)

    create_figure1()
    create_figure2()
    create_figure3()
    create_figure4()

    print("=" * 60)
    print("All figures generated successfully!")
