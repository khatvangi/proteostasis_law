#!/usr/bin/env python3
"""
Figure 4: μ-tAI quadrants and codon mode classifications
Shows the distribution of codons across the four proteostasis modes.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch
import seaborn as sns

# Set style
sns.set_style("white")
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 9

# Load codon modes data
df = pd.read_csv('errors/codon_modes_ecoli.tsv', sep='\t')

# Filter out stop codons and codons without mode assignments
df = df[df['mode'].notna()]
df = df[df['mode'] != 'unknown']

print(f"Loaded {len(df)} codons with mode classifications")
print(f"\nMode distribution:")
print(df['mode'].value_counts())

# Color scheme for modes
mode_colors = {
    'safe_sprinter': '#4DBBD5',      # Blue - low μ, high tAI
    'risky_sprinter': '#E64B35',     # Red - high μ, high tAI
    'safe_careful': '#00A087',       # Green - low μ, low tAI
    'risky_careful': '#F39B7F'       # Orange - high μ, low tAI
}

# Create figure with 2 panels
fig = plt.figure(figsize=(14, 6))

# ============================================================================
# Panel 4A: μ vs tAI scatter plot colored by mode
# ============================================================================
ax1 = plt.subplot(1, 2, 1)

# Use log10 scale for μ to spread out the points
df['log_mu'] = np.log10(df['mu'])

# Get median values for quadrant lines
mu_median = df['mu'].median()
tai_median = df['w_tai'].median()
log_mu_median = df['log_mu'].median()

print(f"\nMedian μ = {mu_median:.6f}, log10(μ) = {log_mu_median:.3f}")
print(f"Median tAI = {tai_median:.3f}")

# Plot each mode
for mode in ['safe_careful', 'risky_careful', 'safe_sprinter', 'risky_sprinter']:
    mode_df = df[df['mode'] == mode]
    ax1.scatter(mode_df['log_mu'], mode_df['w_tai'],
               c=mode_colors[mode], s=100, alpha=0.7,
               edgecolors='black', linewidth=0.5,
               label=mode.replace('_', ' ').title(), zorder=3)

# Add median lines to show quadrants
ax1.axvline(x=log_mu_median, color='gray', linestyle='--', linewidth=2, alpha=0.7, zorder=1)
ax1.axhline(y=tai_median, color='gray', linestyle='--', linewidth=2, alpha=0.7, zorder=1)

# Add quadrant labels
xlim = ax1.get_xlim()
ylim = ax1.get_ylim()
x_range = xlim[1] - xlim[0]
y_range = ylim[1] - ylim[0]

# Position labels in each quadrant
label_offset = 0.05
fontsize_quad = 11

# Top-left: safe_sprinter (low μ, high tAI)
ax1.text(xlim[0] + x_range * label_offset, ylim[1] - y_range * label_offset,
        'Safe Sprinter\n(low μ, high tAI)', fontsize=fontsize_quad, fontweight='bold',
        color=mode_colors['safe_sprinter'], ha='left', va='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor=mode_colors['safe_sprinter'],
                 alpha=0.8, linewidth=2), zorder=2)

# Top-right: risky_sprinter (high μ, high tAI)
ax1.text(xlim[1] - x_range * label_offset, ylim[1] - y_range * label_offset,
        'Risky Sprinter\n(high μ, high tAI)', fontsize=fontsize_quad, fontweight='bold',
        color=mode_colors['risky_sprinter'], ha='right', va='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor=mode_colors['risky_sprinter'],
                 alpha=0.8, linewidth=2), zorder=2)

# Bottom-left: safe_careful (low μ, low tAI)
ax1.text(xlim[0] + x_range * label_offset, ylim[0] + y_range * label_offset,
        'Safe Careful\n(low μ, low tAI)', fontsize=fontsize_quad, fontweight='bold',
        color=mode_colors['safe_careful'], ha='left', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor=mode_colors['safe_careful'],
                 alpha=0.8, linewidth=2), zorder=2)

# Bottom-right: risky_careful (high μ, low tAI)
ax1.text(xlim[1] - x_range * label_offset, ylim[0] + y_range * label_offset,
        'Risky Careful\n(high μ, low tAI)', fontsize=fontsize_quad, fontweight='bold',
        color=mode_colors['risky_careful'], ha='right', va='bottom',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor=mode_colors['risky_careful'],
                 alpha=0.8, linewidth=2), zorder=2)

ax1.set_xlabel('log₁₀(Mutation Rate μ)', fontweight='bold', fontsize=12)
ax1.set_ylabel('tRNA Adaptation Index (tAI)', fontweight='bold', fontsize=12)
ax1.set_title('A. E. coli Codon Distribution Across μ-tAI Space', fontweight='bold', loc='left')
ax1.grid(True, alpha=0.3, zorder=0)

# Add text explaining the median lines
ax1.text(0.98, 0.02, 'Dashed lines show median values',
        transform=ax1.transAxes, fontsize=8, ha='right', va='bottom',
        style='italic', bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.7))

# ============================================================================
# Panel 4B: Cartoon quadrant summary
# ============================================================================
ax2 = plt.subplot(1, 2, 2)
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 10)
ax2.set_aspect('equal')

# Remove axes
ax2.axis('off')

# Draw quadrants as colored boxes
quad_width = 4
quad_height = 4
quad_x_offset = 1
quad_y_offset = 1

quadrants = [
    {'name': 'Safe Sprinter', 'mode': 'safe_sprinter', 'x': quad_x_offset, 'y': quad_y_offset + quad_height,
     'desc': 'Low mutation risk\nHigh translation speed', 'symbol': '✓✓'},
    {'name': 'Risky Sprinter', 'mode': 'risky_sprinter', 'x': quad_x_offset + quad_width, 'y': quad_y_offset + quad_height,
     'desc': 'High mutation risk\nHigh translation speed', 'symbol': '⚠️→'},
    {'name': 'Safe Careful', 'mode': 'safe_careful', 'x': quad_x_offset, 'y': quad_y_offset,
     'desc': 'Low mutation risk\nLow translation speed', 'symbol': '✓⏱'},
    {'name': 'Risky Careful', 'mode': 'risky_careful', 'x': quad_x_offset + quad_width, 'y': quad_y_offset,
     'desc': 'High mutation risk\nLow translation speed', 'symbol': '⚠️⏱'}
]

for quad in quadrants:
    # Draw colored box
    rect = FancyBboxPatch((quad['x'], quad['y']), quad_width, quad_height,
                          boxstyle="round,pad=0.1",
                          facecolor=mode_colors[quad['mode']],
                          edgecolor='black', linewidth=2, alpha=0.3)
    ax2.add_patch(rect)

    # Add quadrant name
    ax2.text(quad['x'] + quad_width/2, quad['y'] + quad_height - 0.5,
            quad['name'], fontsize=13, fontweight='bold', ha='center', va='top',
            color='black')

    # Add description
    ax2.text(quad['x'] + quad_width/2, quad['y'] + quad_height/2,
            quad['desc'], fontsize=10, ha='center', va='center',
            color='black', style='italic')

    # Add count
    count = len(df[df['mode'] == quad['mode']])
    ax2.text(quad['x'] + quad_width/2, quad['y'] + 0.5,
            f'n = {count}', fontsize=11, fontweight='bold', ha='center', va='bottom',
            color='black', bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                                     edgecolor='black', linewidth=1.5))

# Add axis labels
ax2.text(5, 0.3, 'Mutation Rate (μ) →', fontsize=12, fontweight='bold', ha='center')
ax2.text(0.3, 5, 'tAI →', fontsize=12, fontweight='bold', ha='center', rotation=90)

# Add central cross
ax2.plot([quad_x_offset + quad_width, quad_x_offset + quad_width],
         [quad_y_offset, quad_y_offset + 2*quad_height],
         'k-', linewidth=3, zorder=10)
ax2.plot([quad_x_offset, quad_x_offset + 2*quad_width],
         [quad_y_offset + quad_height, quad_y_offset + quad_height],
         'k-', linewidth=3, zorder=10)

# Add title
ax2.text(5, 9.5, 'B. Proteostasis Mode Quadrants', fontsize=12, fontweight='bold', ha='center', va='top')

# Add interpretation box
interpretation = (
    "Each codon classified by its combination of:\n"
    "• Mutation pressure (μ): genetic stability\n"
    "• Translation efficiency (tAI): expression speed\n\n"
    "Different functional contexts favor different modes"
)
ax2.text(5, -0.5, interpretation, fontsize=9, ha='center', va='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='lightcyan', edgecolor='blue', alpha=0.9))

plt.tight_layout()
plt.savefig('figure4_codon_modes.png', dpi=300, bbox_inches='tight')
plt.savefig('figure4_codon_modes.pdf', bbox_inches='tight')
plt.savefig('figure4_codon_modes.svg', bbox_inches='tight')

print("\nFigure 4 saved as figure4_codon_modes.png and .pdf")
print("\nKey insight:")
print("Codons are distributed across all four quadrants, showing that E. coli")
print("uses diverse combinations of mutation resistance and translation speed.")
print("The proteostasis law predicts which mode is selected for each functional context.")
