#!/usr/bin/env python3
"""
Figure 2: Metal ligands vs background - codon enrichment and μ/tAI patterns
Shows falsification of μ-only and tAI-only hypotheses.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Set style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 9

# Load data
summary_df = pd.read_csv('metals/metal_codon_bias_summary.csv')
mu_df = pd.read_csv('metals/metal_codon_bias_with_mu.tsv', sep='\t')
tai_df = pd.read_csv('metals/metal_codon_bias_with_tai.tsv', sep='\t')

# Create figure with 3 panels
fig = plt.figure(figsize=(16, 5))

# Color scheme for amino acids
aa_colors = {
    'ASP': '#E64B35',
    'CYS': '#4DBBD5',
    'GLU': '#00A087',
    'HIS': '#F39B7F'
}

amino_acids = ['ASP', 'CYS', 'GLU', 'HIS']

# Define which pattern each AA follows
mu_pattern = {'ASP': 'higher', 'CYS': 'lower', 'GLU': 'lower', 'HIS': 'higher'}
tai_pattern = {'ASP': 'higher', 'CYS': 'higher', 'GLU': 'lower', 'HIS': 'higher'}

# ============================================================================
# Panel 2A: Codon frequency comparisons (ligand vs background)
# ============================================================================
ax1 = plt.subplot(1, 3, 1)

x_positions = []
current_x = 0
group_centers = []

for aa_idx, aa in enumerate(amino_acids):
    aa_data = summary_df[summary_df['aa'] == aa].copy()

    # Calculate frequencies
    for idx, row in aa_data.iterrows():
        ligand_total = summary_df[summary_df['aa'] == aa]['ligand_count'].sum()
        background_total = summary_df[summary_df['aa'] == aa]['background_count'].sum()

        ligand_freq = row['ligand_count'] / ligand_total
        background_freq = row['background_count'] / background_total

        # Plot paired bars
        width = 0.35
        ax1.bar(current_x, ligand_freq, width, color=aa_colors[aa], alpha=0.8,
                label='Ligand' if aa_idx == 0 and idx == aa_data.index[0] else '')
        ax1.bar(current_x + width, background_freq, width, color=aa_colors[aa], alpha=0.3,
                label='Background' if aa_idx == 0 and idx == aa_data.index[0] else '')

        # Add codon label below
        ax1.text(current_x + width/2, -0.08, row['codon'],
                ha='center', va='top', fontsize=8, fontfamily='monospace')

        x_positions.append(current_x + width/2)
        current_x += 1

    # Mark group center for AA label
    group_centers.append(np.mean(x_positions[-2:]))
    current_x += 0.5  # Space between groups

# Get Fisher's exact test results from mu_df
for aa_idx, aa in enumerate(amino_acids):
    aa_mu = mu_df[mu_df['aa'] == aa].iloc[0]

    # Position for annotation (between the two codon pairs)
    group_start = aa_idx * 2.5
    group_mid = group_start + 0.5

    # Add OR and p-value annotation
    or_val = aa_mu['OR_enrichment']
    p_val = aa_mu['p_fisher']

    y_max = ax1.get_ylim()[1]

    if p_val < 0.001:
        p_text = 'p < 0.001'
    elif p_val < 0.01:
        p_text = f'p = {p_val:.3f}'
    else:
        p_text = f'p = {p_val:.4f}'

    ax1.text(group_mid, y_max * 0.95, f'OR = {or_val:.2f}\n{p_text}',
            ha='center', va='top', fontsize=7, style='italic',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.7))

# Add AA labels below x-axis
for aa_idx, aa in enumerate(amino_acids):
    group_center = aa_idx * 2.5 + 0.5
    ax1.text(group_center, -0.15, aa, ha='center', va='top',
            fontsize=10, fontweight='bold', color=aa_colors[aa])

ax1.set_ylabel('Codon Frequency', fontweight='bold')
ax1.set_xlabel('')
ax1.set_title('A. Ligand vs Background Codon Frequency', fontweight='bold', loc='left')
ax1.set_xticks([])
ax1.legend(loc='upper right', frameon=True)
ax1.set_ylim(0, ax1.get_ylim()[1] * 1.05)

# ============================================================================
# Panel 2B: μ comparison (enriched vs depleted)
# ============================================================================
ax2 = plt.subplot(1, 3, 2)

x_pos = np.arange(len(amino_acids))
width = 0.35

enriched_mu = []
depleted_mu = []

for aa in amino_acids:
    aa_mu = mu_df[mu_df['aa'] == aa].iloc[0]
    enriched_mu.append(aa_mu['enriched_mu'])
    depleted_mu.append(aa_mu['depleted_mu'])

# Create paired bars
for i, aa in enumerate(amino_acids):
    ax2.bar(i - width/2, enriched_mu[i], width,
           color=aa_colors[aa], alpha=0.8, label='Enriched' if i == 0 else '')
    ax2.bar(i + width/2, depleted_mu[i], width,
           color=aa_colors[aa], alpha=0.3, label='Depleted' if i == 0 else '')

    # Add codon labels inside or near bars
    aa_mu = mu_df[mu_df['aa'] == aa].iloc[0]

    # Enriched codon label
    ax2.text(i - width/2, enriched_mu[i]/2, aa_mu['enriched_codon'],
            ha='center', va='center', fontsize=7, fontfamily='monospace',
            rotation=90, color='white' if enriched_mu[i] > 0.002 else 'black')

    # Depleted codon label
    ax2.text(i + width/2, depleted_mu[i]/2, aa_mu['depleted_codon'],
            ha='center', va='center', fontsize=7, fontfamily='monospace',
            rotation=90, color='white' if depleted_mu[i] > 0.002 else 'black')

    # Add visual indicator for μ pattern (contradicts μ-only hypothesis)
    delta = enriched_mu[i] - depleted_mu[i]
    y_arrow = max(enriched_mu[i], depleted_mu[i]) * 1.2

    # Highlight AAs where enriched has LOWER μ (contradicts selection for low mutation)
    if delta < 0:
        # Draw box around these contradictory cases
        rect = plt.Rectangle((i - width, 0), width*2, max(enriched_mu[i], depleted_mu[i]) * 1.1,
                            fill=False, edgecolor='red', linewidth=2, linestyle='--', alpha=0.6)
        ax2.add_patch(rect)
        ax2.text(i, y_arrow * 1.05, '✗ μ-only', ha='center', fontsize=8,
                color='red', fontweight='bold')

ax2.set_ylabel('Mutation Rate (μ)', fontweight='bold')
ax2.set_xlabel('')
ax2.set_title('B. Mutation Rate: Enriched vs Depleted', fontweight='bold', loc='left')
ax2.set_xticks(x_pos)
ax2.set_xticklabels([f'{aa}' for aa in amino_acids])
ax2.legend(loc='upper left', frameon=True)
ax2.set_ylim(0, max(enriched_mu + depleted_mu) * 1.35)

# Add text box explaining the falsification
ax2.text(0.98, 0.98,
        'μ-only falsified:\nCYS & GLU enriched\ncodons have LOWER μ',
        transform=ax2.transAxes, fontsize=9, verticalalignment='top',
        horizontalalignment='right', fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFE6E6', edgecolor='red', alpha=0.9, linewidth=2))

# ============================================================================
# Panel 2C: tAI comparison (enriched vs depleted)
# ============================================================================
ax3 = plt.subplot(1, 3, 3)

enriched_tai = []
depleted_tai = []

for aa in amino_acids:
    aa_tai = tai_df[tai_df['aa'] == aa].iloc[0]
    enriched_tai.append(aa_tai['enriched_tai'])
    depleted_tai.append(aa_tai['depleted_tai'])

# Create paired bars
for i, aa in enumerate(amino_acids):
    ax3.bar(i - width/2, enriched_tai[i], width,
           color=aa_colors[aa], alpha=0.8, label='Enriched' if i == 0 else '')
    ax3.bar(i + width/2, depleted_tai[i], width,
           color=aa_colors[aa], alpha=0.3, label='Depleted' if i == 0 else '')

    # Add codon labels
    aa_tai = tai_df[tai_df['aa'] == aa].iloc[0]

    # Enriched codon label
    ax3.text(i - width/2, enriched_tai[i]/2, aa_tai['enriched_codon'],
            ha='center', va='center', fontsize=7, fontfamily='monospace',
            rotation=90, color='white' if enriched_tai[i] > 0.3 else 'black')

    # Depleted codon label
    ax3.text(i + width/2, depleted_tai[i]/2, aa_tai['depleted_codon'],
            ha='center', va='center', fontsize=7, fontfamily='monospace',
            rotation=90, color='white' if depleted_tai[i] > 0.3 else 'black')

    # Add visual indicator for tAI pattern (contradicts tAI-only hypothesis)
    delta = enriched_tai[i] - depleted_tai[i]
    y_arrow = max(enriched_tai[i], depleted_tai[i]) * 1.2

    # Highlight AAs where enriched has LOWER tAI (contradicts selection for high translation)
    if delta < 0:
        # Draw box around these contradictory cases
        rect = plt.Rectangle((i - width, 0), width*2, max(enriched_tai[i], depleted_tai[i]) * 1.1,
                            fill=False, edgecolor='red', linewidth=2, linestyle='--', alpha=0.6)
        ax3.add_patch(rect)
        ax3.text(i, y_arrow * 1.05, '✗ tAI-only', ha='center', fontsize=8,
                color='red', fontweight='bold')

ax3.set_ylabel('tRNA Adaptation Index (tAI)', fontweight='bold')
ax3.set_xlabel('')
ax3.set_title('C. tAI: Enriched vs Depleted', fontweight='bold', loc='left')
ax3.set_xticks(x_pos)
ax3.set_xticklabels([f'{aa}' for aa in amino_acids])
ax3.legend(loc='upper left', frameon=True)
ax3.set_ylim(0, max(enriched_tai + depleted_tai) * 1.35)

# Add text box explaining the falsification
ax3.text(0.98, 0.98,
        'tAI-only falsified:\nGLU enriched codon\nhas LOWER tAI',
        transform=ax3.transAxes, fontsize=9, verticalalignment='top',
        horizontalalignment='right', fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFE6E6', edgecolor='red', alpha=0.9, linewidth=2))

plt.tight_layout()
plt.savefig('figure2_metal_ligands.png', dpi=300, bbox_inches='tight')
plt.savefig('figure2_metal_ligands.pdf', bbox_inches='tight')
plt.savefig('figure2_metal_ligands.svg', bbox_inches='tight')
print("Figure 2 saved as figure2_metal_ligands.png and .pdf")
print("\nKey findings:")
print("- Panel A shows significant codon enrichment in metal ligands")
print("- Panel B demonstrates μ cannot explain all patterns (CYS & GLU enriched have LOWER μ)")
print("- Panel C demonstrates tAI cannot explain all patterns (GLU enriched has LOWER tAI)")
print("- Together, these panels falsify both μ-only and tAI-only hypotheses")
