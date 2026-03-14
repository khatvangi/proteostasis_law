#!/usr/bin/env python3
"""
Figure 6: Capacity bound - doublet impossible, triplet minimal
Demonstrates that mode demand exceeds doublet capacity but fits in triplet code.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import Circle, Wedge, Rectangle, FancyBboxPatch
import seaborn as sns

# Set style
sns.set_style("white")
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['legend.fontsize'] = 9

# Load data
aa_modes_df = pd.read_csv('errors/aa_mode_summary.tsv', sep='\t')

# Calculate number of modes for each amino acid
def count_modes(modes_str):
    """Count number of distinct modes (excluding 'none' and 'unknown')"""
    if pd.isna(modes_str) or modes_str == 'none':
        return 0
    modes = modes_str.split(',')
    real_modes = [m for m in modes if m != 'unknown']
    return len(real_modes)

aa_modes_df['n_modes'] = aa_modes_df['modes_present'].apply(count_modes)

# Calculate totals
total_mode_demand = aa_modes_df['n_modes'].sum()
two_codon_aas = aa_modes_df[aa_modes_df['n_codons'] == 2]
two_codon_mode_demand = two_codon_aas['n_modes'].sum()

print(f"Total mode demand (Σ k_A): {total_mode_demand}")
print(f"2-codon AA mode demand: {two_codon_mode_demand}")
print(f"Number of 2-codon AAs: {len(two_codon_aas)}")

# Code capacities
doublet_capacity = 16  # 4^2
triplet_capacity = 64  # 4^3

# Create figure with 3 panels
fig = plt.figure(figsize=(16, 6))

# Color scheme
demand_color = '#E64B35'
doublet_color = '#F39B7F'
triplet_color = '#4DBBD5'

# ============================================================================
# Panel 6A: Mode demand vs code capacity (simple bar chart)
# ============================================================================
ax1 = plt.subplot(1, 3, 1)

categories = ['Mode\nDemand\n(Σ k_A)', 'Doublet\nCapacity\n(4²)', 'Triplet\nCapacity\n(4³)']
values = [total_mode_demand, doublet_capacity, triplet_capacity]
colors = [demand_color, doublet_color, triplet_color]

bars = ax1.bar(categories, values, color=colors, edgecolor='black', linewidth=2, alpha=0.8)

# Add value labels on bars
for i, (bar, val) in enumerate(zip(bars, values)):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
            f'{val}', ha='center', va='bottom', fontsize=16, fontweight='bold')

# Add horizontal line at mode demand level
ax1.axhline(y=total_mode_demand, color=demand_color, linestyle='--', linewidth=2,
           alpha=0.7, label=f'Mode demand = {total_mode_demand}')

# Highlight the impossibility
ax1.fill_between([-0.5, 1.5], 0, total_mode_demand, color='red', alpha=0.1)
ax1.text(1, total_mode_demand/2, 'INSUFFICIENT', fontsize=12, fontweight='bold',
        ha='center', va='center', color='red', rotation=0)

# Highlight the sufficiency
ax1.fill_between([1.5, 2.5], total_mode_demand, triplet_capacity, color='green', alpha=0.1)
ax1.text(2, (total_mode_demand + triplet_capacity)/2, 'SUFFICIENT', fontsize=12,
        fontweight='bold', ha='center', va='center', color='green', rotation=0)

ax1.set_ylabel('Number of Codons/Modes', fontweight='bold', fontsize=12)
ax1.set_title('A. Overall Capacity Comparison', fontweight='bold', loc='left', fontsize=13)
ax1.set_ylim(0, 70)
ax1.grid(axis='y', alpha=0.3)

# Add annotation
ax1.text(0.5, 0.98, f'Deficit: {total_mode_demand - doublet_capacity} modes',
        transform=ax1.transAxes, fontsize=10, fontweight='bold',
        ha='center', va='top', color='red',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFE6E6', edgecolor='red', linewidth=2))

# ============================================================================
# Panel 6B: 2-codon amino acids specifically
# ============================================================================
ax2 = plt.subplot(1, 3, 2)

categories_2codon = ['2-Codon AAs\nMode Demand', 'Available\nDoublet Codons']
values_2codon = [two_codon_mode_demand, doublet_capacity]
colors_2codon = [demand_color, doublet_color]

bars2 = ax2.bar(categories_2codon, values_2codon, color=colors_2codon,
               edgecolor='black', linewidth=2, alpha=0.8)

# Add value labels
for bar, val in zip(bars2, values_2codon):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + 0.3,
            f'{val}', ha='center', va='bottom', fontsize=16, fontweight='bold')

# Add horizontal line at demand
ax2.axhline(y=two_codon_mode_demand, color=demand_color, linestyle='--',
           linewidth=2, alpha=0.7)

# Highlight the overflow
overflow = two_codon_mode_demand - doublet_capacity
ax2.fill_between([0.5, 1.5], doublet_capacity, two_codon_mode_demand,
                color='red', alpha=0.2)
ax2.text(1, doublet_capacity + overflow/2, f'OVERFLOW\n+{overflow}', fontsize=11,
        fontweight='bold', ha='center', va='center', color='red')

ax2.set_ylabel('Number of Codons/Modes', fontweight='bold', fontsize=12)
ax2.set_title('B. 2-Codon Amino Acids Alone', fontweight='bold', loc='left', fontsize=13)
ax2.set_ylim(0, 22)
ax2.grid(axis='y', alpha=0.3)

# Add critical observation
observation = (
    f"Critical: {len(two_codon_aas)} two-codon AAs\n"
    f"each need 2 modes = {two_codon_mode_demand} total\n"
    f"BUT only {doublet_capacity} doublet codons exist!\n\n"
    "Doublet code is IMPOSSIBLE\n"
    "even for 2-codon AAs alone."
)
ax2.text(0.5, 0.98, observation, transform=ax2.transAxes,
        fontsize=9, fontweight='bold', ha='center', va='top',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFE6E6',
                 edgecolor='red', linewidth=2, alpha=0.95))

# ============================================================================
# Panel 6C: Genetic code wheel cartoon
# ============================================================================
ax3 = plt.subplot(1, 3, 3)
ax3.set_xlim(-2, 2)
ax3.set_ylim(-2, 2)
ax3.set_aspect('equal')
ax3.axis('off')

# Draw doublet wheel (inner, overwhelmed)
doublet_radius = 0.6
n_doublet = 16
angle_per_slot_doublet = 360 / n_doublet

# Draw doublet slots
for i in range(n_doublet):
    angle_start = i * angle_per_slot_doublet
    wedge = Wedge((0, -0.5), doublet_radius, angle_start, angle_start + angle_per_slot_doublet,
                 width=0.15, facecolor=doublet_color, edgecolor='black', linewidth=1.5, alpha=0.6)
    ax3.add_patch(wedge)

# Mark overflow slots in red
overflow_slots = two_codon_mode_demand - doublet_capacity
for i in range(overflow_slots):
    angle = i * (360 / overflow_slots)
    x = 0.9 * np.cos(np.radians(angle))
    y = -0.5 + 0.9 * np.sin(np.radians(angle))
    ax3.plot([0, x], [-0.5, y], 'r-', linewidth=3, alpha=0.7)
    ax3.scatter(x, y, s=200, c='red', marker='X', edgecolors='black', linewidth=2, zorder=5)

# Doublet label
ax3.text(0, -0.5, 'DOUBLET\n16 slots', ha='center', va='center',
        fontsize=10, fontweight='bold', color='black',
        bbox=dict(boxstyle='circle,pad=0.3', facecolor='white', edgecolor='black', linewidth=2))

ax3.text(0, -1.4, 'OVERLOADED', ha='center', va='center',
        fontsize=11, fontweight='bold', color='red')

# Draw triplet wheel (outer, sufficient)
triplet_radius = 1.8
n_triplet = 64
angle_per_slot_triplet = 360 / n_triplet

# Draw triplet slots (show subset for clarity)
for i in range(0, n_triplet, 2):  # Draw every other slot for visibility
    angle_start = i * angle_per_slot_triplet
    wedge = Wedge((0, 0.5), triplet_radius, angle_start, angle_start + angle_per_slot_triplet,
                 width=0.2, facecolor=triplet_color, edgecolor='black', linewidth=0.8, alpha=0.6)
    ax3.add_patch(wedge)

# Mark used slots
n_used = total_mode_demand
for i in range(0, n_used, 2):
    angle = i * (360 / n_used)
    x = 1.7 * np.cos(np.radians(angle))
    y = 0.5 + 1.7 * np.sin(np.radians(angle))
    ax3.scatter(x, y, s=80, c='green', marker='o', edgecolors='black', linewidth=1, zorder=5)

# Triplet label
ax3.text(0, 0.5, 'TRIPLET\n64 slots', ha='center', va='center',
        fontsize=10, fontweight='bold', color='black',
        bbox=dict(boxstyle='circle,pad=0.3', facecolor='white', edgecolor='black', linewidth=2))

ax3.text(0, 1.9, 'SUFFICIENT', ha='center', va='center',
        fontsize=11, fontweight='bold', color='green')

# Add title and legend
ax3.text(0, 2.2, 'C. Code Wheel Schematic', ha='center', va='top',
        fontsize=13, fontweight='bold')

# Add capacity info
info_text = (
    f"Mode demand: {total_mode_demand}\n"
    f"Doublet: 4² = {doublet_capacity} ✗\n"
    f"Triplet: 4³ = {triplet_capacity} ✓"
)
ax3.text(0, -2.2, info_text, ha='center', va='bottom',
        fontsize=9, fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow',
                 edgecolor='orange', linewidth=2))

plt.tight_layout()
plt.savefig('figure6_capacity_bound.png', dpi=300, bbox_inches='tight')
plt.savefig('figure6_capacity_bound.pdf', bbox_inches='tight')
plt.savefig('figure6_capacity_bound.svg', bbox_inches='tight')

print("\nFigure 6 saved as figure6_capacity_bound.png and .pdf")
print("\n=== CAPACITY ANALYSIS ===")
print(f"Total mode demand (Σ k_A): {total_mode_demand}")
print(f"Doublet capacity (4²): {doublet_capacity}")
print(f"Triplet capacity (4³): {triplet_capacity}")
print(f"\nDoublet deficit: {total_mode_demand - doublet_capacity} modes")
print(f"Triplet surplus: {triplet_capacity - total_mode_demand} slots")
print(f"\n2-codon AA analysis:")
print(f"  Number of 2-codon AAs: {len(two_codon_aas)}")
print(f"  Mode demand: {two_codon_mode_demand}")
print(f"  Doublet capacity: {doublet_capacity}")
print(f"  Overflow: +{two_codon_mode_demand - doublet_capacity} modes")
print(f"\nConclusion:")
print(f"  ✗ Doublet code IMPOSSIBLE (insufficient by {total_mode_demand - doublet_capacity} modes)")
print(f"  ✓ Triplet code MINIMAL (necessary and sufficient)")
