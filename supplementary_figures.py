#!/usr/bin/env python3
"""
Supplementary Figures for Proteostasis Law Manuscript
S1: Distribution of μ and tAI across codons
S2: Metal-binding fraction per codon
S3: ROC curves for μ-only vs tAI-only models
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import csv
import math
import statsmodels.api as sm
from sklearn.metrics import roc_curve, auc

# Set style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 9

print("Loading data for supplementary figures...")

# Load codon-level data
mu_df = pd.read_csv('errors/codon_error_rates.tsv', sep='\t')
tai_df = pd.read_csv('errors/ecoli_tai_ws.tsv', sep='\t')

# Merge μ and tAI data
codon_props = mu_df.merge(tai_df, on='codon')

print(f"Loaded properties for {len(codon_props)} codons")

# ============================================================================
# Supplementary Figure S1: Distribution of μ and tAI
# ============================================================================
print("\nCreating Supplementary Figure S1...")

fig_s1 = plt.figure(figsize=(12, 5))

# S1A: μ distribution
ax1 = plt.subplot(1, 2, 1)
mu_values = codon_props['mu'].values

# Use log scale for better visualization
log_mu = np.log10(mu_values)

ax1.hist(log_mu, bins=20, color='#E64B35', alpha=0.7, edgecolor='black', linewidth=1.5)
ax1.axvline(np.median(log_mu), color='darkred', linestyle='--', linewidth=2,
           label=f'Median = {np.median(log_mu):.2f}')
ax1.axvline(np.mean(log_mu), color='orange', linestyle='--', linewidth=2,
           label=f'Mean = {np.mean(log_mu):.2f}')

ax1.set_xlabel('log₁₀(Mutation Rate μ)', fontweight='bold')
ax1.set_ylabel('Number of Codons', fontweight='bold')
ax1.set_title('S1A. Distribution of Mutation Rates', fontweight='bold', loc='left')
ax1.legend(loc='upper right', frameon=True)
ax1.grid(True, alpha=0.3)

# Add stats text
stats_text = f'Range: [{log_mu.min():.2f}, {log_mu.max():.2f}]\nStd: {log_mu.std():.2f}'
ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes,
        fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', alpha=0.8))

# S1B: tAI distribution
ax2 = plt.subplot(1, 2, 2)
tai_values = codon_props['w_tai'].values

ax2.hist(tai_values, bins=20, color='#4DBBD5', alpha=0.7, edgecolor='black', linewidth=1.5)
ax2.axvline(np.median(tai_values), color='darkblue', linestyle='--', linewidth=2,
           label=f'Median = {np.median(tai_values):.3f}')
ax2.axvline(np.mean(tai_values), color='cyan', linestyle='--', linewidth=2,
           label=f'Mean = {np.mean(tai_values):.3f}')

ax2.set_xlabel('tRNA Adaptation Index (tAI)', fontweight='bold')
ax2.set_ylabel('Number of Codons', fontweight='bold')
ax2.set_title('S1B. Distribution of Translation Efficiency', fontweight='bold', loc='left')
ax2.legend(loc='upper right', frameon=True)
ax2.grid(True, alpha=0.3)

# Add stats text
stats_text = f'Range: [{tai_values.min():.3f}, {tai_values.max():.3f}]\nStd: {tai_values.std():.3f}'
ax2.text(0.02, 0.98, stats_text, transform=ax2.transAxes,
        fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='lightcyan', alpha=0.8))

plt.tight_layout()
plt.savefig('figureS1_distributions.png', dpi=300, bbox_inches='tight')
plt.savefig('figureS1_distributions.pdf', bbox_inches='tight')
plt.savefig('figureS1_distributions.svg', bbox_inches='tight')
print("Figure S1 saved")

# ============================================================================
# Supplementary Figure S2: Metal-binding fraction per codon
# ============================================================================
print("\nCreating Supplementary Figure S2...")

# Load residue-level data
residue_df = pd.read_csv('errors/residue_kappa_table.tsv', sep='\t')

# Calculate metal fraction per codon
codon_metal_stats = residue_df.groupby('codon').agg({
    'is_metal': ['sum', 'count', 'mean']
}).reset_index()
codon_metal_stats.columns = ['codon', 'metal_count', 'total_count', 'metal_fraction']

# Merge with codon properties
codon_metal_stats = codon_metal_stats.merge(codon_props, on='codon')

# Add amino acid info
codon_to_aa = {}
with open('errors/residue_kappa_table.tsv') as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        codon_to_aa[row['codon']] = row['aa']

codon_metal_stats['aa'] = codon_metal_stats['codon'].map(codon_to_aa)

fig_s2 = plt.figure(figsize=(14, 6))

# S2A: Metal fraction distribution
ax1 = plt.subplot(1, 2, 1)

metal_fractions = codon_metal_stats['metal_fraction'].values * 100  # Convert to percentage

ax1.hist(metal_fractions, bins=25, color='#00A087', alpha=0.7, edgecolor='black', linewidth=1.5)
ax1.axvline(np.median(metal_fractions), color='darkgreen', linestyle='--', linewidth=2,
           label=f'Median = {np.median(metal_fractions):.2f}%')
ax1.axvline(np.mean(metal_fractions), color='lightgreen', linestyle='--', linewidth=2,
           label=f'Mean = {np.mean(metal_fractions):.2f}%')

ax1.set_xlabel('Metal-binding Fraction (%)', fontweight='bold')
ax1.set_ylabel('Number of Codons', fontweight='bold')
ax1.set_title('S2A. Distribution of Metal-binding Frequency', fontweight='bold', loc='left')
ax1.legend(loc='upper right', frameon=True)
ax1.grid(True, alpha=0.3)

stats_text = f'Range: [{metal_fractions.min():.2f}%, {metal_fractions.max():.2f}%]\nStd: {metal_fractions.std():.2f}%'
ax1.text(0.02, 0.98, stats_text, transform=ax1.transAxes,
        fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgreen', alpha=0.8))

# S2B: Metal fraction vs μ and tAI
ax2 = plt.subplot(1, 2, 2)

# Color by metal fraction
scatter = ax2.scatter(codon_metal_stats['mu'], codon_metal_stats['w_tai'],
                     c=codon_metal_stats['metal_fraction']*100, s=150,
                     cmap='RdYlGn', edgecolors='black', linewidth=1.5, alpha=0.8)

# Add colorbar
cbar = plt.colorbar(scatter, ax=ax2)
cbar.set_label('Metal-binding\nFraction (%)', fontweight='bold', rotation=270, labelpad=20)

ax2.set_xlabel('Mutation Rate (μ)', fontweight='bold')
ax2.set_ylabel('tRNA Adaptation Index (tAI)', fontweight='bold')
ax2.set_title('S2B. Metal-binding Fraction in μ-tAI Space', fontweight='bold', loc='left')
ax2.set_xscale('log')
ax2.grid(True, alpha=0.3)

# Highlight high metal-binding codons
high_metal = codon_metal_stats[codon_metal_stats['metal_fraction'] > 0.05]
for _, row in high_metal.iterrows():
    ax2.text(row['mu'], row['w_tai'], row['codon'], fontsize=7,
            ha='center', va='bottom', fontweight='bold')

plt.tight_layout()
plt.savefig('figureS2_metal_fractions.png', dpi=300, bbox_inches='tight')
plt.savefig('figureS2_metal_fractions.pdf', bbox_inches='tight')
plt.savefig('figureS2_metal_fractions.svg', bbox_inches='tight')
print("Figure S2 saved")

# ============================================================================
# Supplementary Figure S3: ROC curves for μ-only vs tAI-only models
# ============================================================================
print("\nCreating Supplementary Figure S3...")

TWO_CODON_AA = {
    "C": ["TGT","TGC"],
    "D": ["GAT","GAC"],
    "E": ["GAA","GAG"],
    "F": ["TTT","TTC"],
    "H": ["CAT","CAC"],
    "K": ["AAA","AAG"],
    "N": ["AAT","AAC"],
    "Q": ["CAA","CAG"],
    "Y": ["TAT","TAC"],
}

# Collect data per AA
data_by_aa = defaultdict(list)

with open('errors/residue_kappa_table.tsv') as f:
    reader = csv.DictReader(f, delimiter="\t")
    for row in reader:
        aa = row["aa"]
        if aa not in TWO_CODON_AA:
            continue
        codon = row["codon"]
        if codon not in TWO_CODON_AA[aa]:
            continue
        is_metal = int(row["is_metal"])
        try:
            mu = float(row["mu"])
            tai = float(row["w_tai"])
        except ValueError:
            continue
        if math.isnan(mu) or math.isnan(tai):
            continue
        data_by_aa[aa].append((mu, tai, is_metal))

# Create ROC plots
fig_s3 = plt.figure(figsize=(15, 10))

# Color scheme for amino acids
aa_colors = {
    'C': '#E64B35', 'D': '#4DBBD5', 'E': '#00A087', 'F': '#F39B7F',
    'H': '#3C5488', 'K': '#8491B4', 'N': '#91D1C2', 'Q': '#DC0000', 'Y': '#7E6148'
}

plot_idx = 1
roc_summary = []

for aa in sorted(TWO_CODON_AA.keys()):
    if aa not in data_by_aa or len(data_by_aa[aa]) < 50:
        continue

    ax = plt.subplot(3, 3, plot_idx)
    plot_idx += 1

    records = data_by_aa[aa]
    X_mu = np.array([r[0] for r in records])
    X_tai = np.array([r[1] for r in records])
    y = np.array([r[2] for r in records])

    # Fit models
    try:
        X_mu_ = sm.add_constant(X_mu.reshape(-1, 1))
        model_mu = sm.Logit(y, X_mu_).fit(disp=0)
        pred_mu = model_mu.predict(X_mu_)

        X_tai_ = sm.add_constant(X_tai.reshape(-1, 1))
        model_tai = sm.Logit(y, X_tai_).fit(disp=0)
        pred_tai = model_tai.predict(X_tai_)

        # Calculate ROC curves
        fpr_mu, tpr_mu, _ = roc_curve(y, pred_mu)
        auc_mu = auc(fpr_mu, tpr_mu)

        fpr_tai, tpr_tai, _ = roc_curve(y, pred_tai)
        auc_tai = auc(fpr_tai, tpr_tai)

        # Plot ROC curves
        ax.plot(fpr_mu, tpr_mu, '-', linewidth=2, color='#E64B35',
               label=f'μ-only (AUC={auc_mu:.3f})', alpha=0.8)
        ax.plot(fpr_tai, tpr_tai, '-', linewidth=2, color='#4DBBD5',
               label=f'tAI-only (AUC={auc_tai:.3f})', alpha=0.8)
        ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5, label='Random')

        ax.set_xlabel('False Positive Rate', fontsize=9)
        ax.set_ylabel('True Positive Rate', fontsize=9)
        ax.set_title(f'{aa} (n={len(records)})', fontweight='bold', fontsize=11)
        ax.legend(loc='lower right', fontsize=8, frameon=True)
        ax.grid(True, alpha=0.3)
        ax.set_xlim([0, 1])
        ax.set_ylim([0, 1])

        # Store summary
        roc_summary.append({
            'aa': aa,
            'n': len(records),
            'auc_mu': auc_mu,
            'auc_tai': auc_tai,
            'better': 'μ' if auc_mu > auc_tai else 'tAI'
        })

    except Exception as e:
        print(f"Warning: Could not fit model for {aa}: {e}")
        continue

plt.tight_layout()
plt.savefig('figureS3_roc_curves.png', dpi=300, bbox_inches='tight')
plt.savefig('figureS3_roc_curves.pdf', bbox_inches='tight')
plt.savefig('figureS3_roc_curves.svg', bbox_inches='tight')
print("Figure S3 saved")

# Print summary
print("\n=== ROC ANALYSIS SUMMARY ===")
roc_df = pd.DataFrame(roc_summary)
print(roc_df.to_string(index=False))
print(f"\nμ-dominated: {sum(roc_df['better'] == 'μ')} amino acids")
print(f"tAI-dominated: {sum(roc_df['better'] == 'tAI')} amino acids")

print("\n=== ALL SUPPLEMENTARY FIGURES COMPLETE ===")
print("Figure S1: figureS1_distributions.png/pdf")
print("Figure S2: figureS2_metal_fractions.png/pdf")
print("Figure S3: figureS3_roc_curves.png/pdf")
