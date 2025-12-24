#!/usr/bin/env python3
"""
Figure 3: Logistic regression of metal sites vs μ and tAI (residue-level)
Shows that different amino acids use different mixtures of μ and tAI.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from collections import defaultdict
import csv
import math
import statsmodels.api as sm

# Set style
sns.set_style("whitegrid")
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 9

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

# Run logistic regression analysis
print("Running logistic regression for 2-codon amino acids...")

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

# Store results
results = []

for aa, records in sorted(data_by_aa.items()):
    if len(records) < 50:
        continue

    X_mu = np.array([[r[0]] for r in records])
    X_tai = np.array([[r[1]] for r in records])
    y = np.array([r[2] for r in records])

    # Add intercept
    X_mu_ = sm.add_constant(X_mu)
    X_tai_ = sm.add_constant(X_tai)

    # μ-only model
    try:
        model_mu = sm.Logit(y, X_mu_).fit(disp=0)
        beta_mu = model_mu.params[1]
        p_mu = model_mu.pvalues[1]
        stderr_mu = model_mu.bse[1]
    except:
        beta_mu = 0
        p_mu = 1
        stderr_mu = 0

    # tAI-only model
    try:
        model_tai = sm.Logit(y, X_tai_).fit(disp=0)
        beta_tai = model_tai.params[1]
        p_tai = model_tai.pvalues[1]
        stderr_tai = model_tai.bse[1]
    except:
        beta_tai = 0
        p_tai = 1
        stderr_tai = 0

    results.append({
        'aa': aa,
        'n': len(records),
        'beta_mu': beta_mu,
        'p_mu': p_mu,
        'stderr_mu': stderr_mu,
        'beta_tai': beta_tai,
        'p_tai': p_tai,
        'stderr_tai': stderr_tai,
        'sig_mu': p_mu < 0.05,
        'sig_tai': p_tai < 0.05
    })

results_df = pd.DataFrame(results)
print(f"\nProcessed {len(results_df)} amino acids")

# Categorize amino acids based on which predictor is significant
def categorize(row):
    if row['sig_mu'] and not row['sig_tai']:
        return 'μ-dominated'
    elif row['sig_tai'] and not row['sig_mu']:
        return 'tAI-dominated'
    elif row['sig_mu'] and row['sig_tai']:
        return 'Both significant'
    else:
        return 'Weak/no signal'

results_df['category'] = results_df.apply(categorize, axis=1)

# Color scheme
category_colors = {
    'μ-dominated': '#E64B35',
    'tAI-dominated': '#4DBBD5',
    'Both significant': '#00A087',
    'Weak/no signal': '#CCCCCC'
}

# Create figure with 2 panels
fig = plt.figure(figsize=(14, 6))

# ============================================================================
# Panel 3A: Bar chart of β_μ and β_tAI with error bars
# ============================================================================
ax1 = plt.subplot(1, 2, 1)

amino_acids = results_df['aa'].tolist()
x = np.arange(len(amino_acids))
width = 0.35

# Normalize coefficients for better visualization (since magnitudes vary wildly)
# Use sign-preserving log scale
def signed_log(vals, threshold=1e-6):
    """Apply log transform while preserving sign and handling small values"""
    result = []
    for v in vals:
        if abs(v) < threshold:
            result.append(0)
        else:
            result.append(np.sign(v) * np.log10(abs(v) + 1))
    return np.array(result)

beta_mu_display = signed_log(results_df['beta_mu'].values)
beta_tai_display = results_df['beta_tai'].values  # tAI coefficients are reasonable scale

# Plot bars
for i, idx in enumerate(results_df.index):
    row = results_df.loc[idx]
    color = category_colors[row['category']]

    # β_μ bar
    bar_mu = ax1.bar(i - width/2, beta_mu_display[i], width,
                     color=color if row['sig_mu'] else '#EEEEEE',
                     alpha=0.8, edgecolor='black' if row['sig_mu'] else 'gray',
                     linewidth=2 if row['sig_mu'] else 0.5,
                     label='β_μ (μ-only)' if i == 0 else '')

    # β_tAI bar
    bar_tai = ax1.bar(i + width/2, beta_tai_display[i], width,
                      color=color if row['sig_tai'] else '#EEEEEE',
                      alpha=0.8, edgecolor='black' if row['sig_tai'] else 'gray',
                      linewidth=2 if row['sig_tai'] else 0.5,
                      label='β_tAI (tAI-only)' if i == 0 else '')

    # Add significance stars
    if row['sig_mu']:
        y_pos = beta_mu_display[i] + (0.3 if beta_mu_display[i] > 0 else -0.5)
        ax1.text(i - width/2, y_pos, '*', ha='center', va='bottom' if beta_mu_display[i] > 0 else 'top',
                fontsize=16, fontweight='bold', color='black')

    if row['sig_tai']:
        y_pos = beta_tai_display[i] + (0.3 if beta_tai_display[i] > 0 else -0.5)
        ax1.text(i + width/2, y_pos, '*', ha='center', va='bottom' if beta_tai_display[i] > 0 else 'top',
                fontsize=16, fontweight='bold', color='black')

ax1.set_ylabel('Coefficient Value', fontweight='bold')
ax1.set_xlabel('Amino Acid', fontweight='bold')
ax1.set_title('A. Logistic Regression Coefficients\n(* p < 0.05)', fontweight='bold', loc='left')
ax1.set_xticks(x)
ax1.set_xticklabels(amino_acids, fontweight='bold')
ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.8)
ax1.legend(loc='upper left', frameon=True)
ax1.grid(axis='y', alpha=0.3)

# Add note about μ coefficients
ax1.text(0.02, 0.02, 'Note: β_μ shown on log scale (sign-preserved)',
        transform=ax1.transAxes, fontsize=8, style='italic',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.7))

# ============================================================================
# Panel 3B: Scatter plot of β_tAI vs β_μ
# ============================================================================
ax2 = plt.subplot(1, 2, 2)

# Plot each amino acid
seen_categories = []
for i, idx in enumerate(results_df.index):
    row = results_df.loc[idx]
    color = category_colors[row['category']]
    marker_size = 200 if (row['sig_mu'] or row['sig_tai']) else 100

    # Use actual coefficients for scatter
    # But clip extreme values for visualization
    mu_val = np.clip(row['beta_mu'], -1e5, 1e5)
    tai_val = row['beta_tai']

    ax2.scatter(mu_val, tai_val, s=marker_size, c=color, alpha=0.7,
               edgecolors='black' if (row['sig_mu'] or row['sig_tai']) else 'gray',
               linewidth=2 if (row['sig_mu'] or row['sig_tai']) else 0.5,
               label=row['category'] if row['category'] not in seen_categories else '')

    if row['category'] not in seen_categories:
        seen_categories.append(row['category'])

    # Add amino acid label
    offset_x = mu_val * 0.05 if mu_val != 0 else 1000
    offset_y = tai_val * 0.05 if tai_val != 0 else 0.1
    ax2.text(mu_val + offset_x, tai_val + offset_y, row['aa'],
            fontsize=10, fontweight='bold', ha='center', va='center',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray', alpha=0.8))

# Add quadrant lines
ax2.axhline(y=0, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax2.axvline(x=0, color='black', linestyle='--', linewidth=1, alpha=0.5)

ax2.set_xlabel('β_μ (μ-only model)', fontweight='bold')
ax2.set_ylabel('β_tAI (tAI-only model)', fontweight='bold')
ax2.set_title('B. μ vs tAI Dominance by Amino Acid', fontweight='bold', loc='left')

# Add legend for categories
handles = []
labels = []
for cat, col in category_colors.items():
    if cat in results_df['category'].values:
        handles.append(plt.scatter([], [], s=150, c=col, edgecolors='black', linewidth=2, alpha=0.7))
        labels.append(cat)
ax2.legend(handles, labels, loc='best', frameon=True, title='Category')

ax2.grid(True, alpha=0.3)

# Add interpretation box
interpretation = (
    "Different amino acids show\n"
    "different selection pressures:\n"
    "• K: μ-dominated\n"
    "• D,E,F,H,N: tAI-dominated\n"
    "• C: Both significant\n"
    "• Y,Q: Weak signal"
)
ax2.text(0.02, 0.98, interpretation,
        transform=ax2.transAxes, fontsize=9, verticalalignment='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='lightcyan', edgecolor='blue', alpha=0.9))

plt.tight_layout()
plt.savefig('figure3_logistic_regression.png', dpi=300, bbox_inches='tight')
plt.savefig('figure3_logistic_regression.pdf', bbox_inches='tight')
plt.savefig('figure3_logistic_regression.svg', bbox_inches='tight')

# Save summary table
results_df.to_csv('figure3_logistic_coefficients.csv', index=False)

print("\nFigure 3 saved as figure3_logistic_regression.png and .pdf")
print("Coefficient table saved as figure3_logistic_coefficients.csv")
print("\nSummary by category:")
print(results_df.groupby('category').size())
