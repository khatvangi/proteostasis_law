#!/usr/bin/env python3
"""
Figure 1: Conceptual schematic of the Proteostatic Codon Law
Panel A: Trade-off triangle between μ, tAI, and pairing safety
Panel B: κ as a "gain knob" amplifying or damping codon cost J_i(c)
Panel C: μ-tAI modes per amino acid
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle, Wedge, Rectangle, Polygon
from matplotlib.patches import ConnectionPatch
import matplotlib.lines as mlines

# Set style
plt.rcParams['font.size'] = 10
plt.rcParams['axes.labelsize'] = 11
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['font.family'] = 'sans-serif'

# Color scheme
color_mu = '#E64B35'
color_tai = '#4DBBD5'
color_safety = '#00A087'
color_kappa = '#F39B7F'
mode_colors = {
    'safe_sprinter': '#4DBBD5',
    'risky_sprinter': '#E64B35',
    'safe_careful': '#00A087',
    'risky_careful': '#F39B7F'
}

# Create figure with 3 panels
fig = plt.figure(figsize=(16, 5))

# ============================================================================
# Panel 1A: Trade-off Triangle
# ============================================================================
ax1 = plt.subplot(1, 3, 1)
ax1.set_xlim(0, 10)
ax1.set_ylim(0, 10)
ax1.set_aspect('equal')
ax1.axis('off')

# Draw triangle
triangle_points = np.array([
    [5, 8.5],   # Top vertex (μ)
    [1.5, 2],   # Bottom left vertex (tAI)
    [8.5, 2]    # Bottom right vertex (Safety)
])

triangle = Polygon(triangle_points, fill=False, edgecolor='black', linewidth=3)
ax1.add_patch(triangle)

# Add vertex labels
ax1.text(5, 9.2, 'Low Mutation\nRate (μ↓)', ha='center', va='bottom',
        fontsize=12, fontweight='bold', color=color_mu,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor=color_mu, linewidth=2))

ax1.text(1.5, 1.2, 'High Translation\nEfficiency (tAI↑)', ha='center', va='top',
        fontsize=12, fontweight='bold', color=color_tai,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor=color_tai, linewidth=2))

ax1.text(8.5, 1.2, 'Pairing\nSafety', ha='center', va='top',
        fontsize=12, fontweight='bold', color=color_safety,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white', edgecolor=color_safety, linewidth=2))

# Add edges with trade-off descriptions
# μ vs tAI edge
mid_mu_tai = (triangle_points[0] + triangle_points[1]) / 2
ax1.text(mid_mu_tai[0] - 1.2, mid_mu_tai[1], 'μ-tAI\ntrade-off', ha='center', va='center',
        fontsize=10, style='italic', rotation=60,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow', alpha=0.8))

# μ vs Safety edge
mid_mu_safety = (triangle_points[0] + triangle_points[2]) / 2
ax1.text(mid_mu_safety[0] + 1.2, mid_mu_safety[1], 'μ-safety\ntrade-off', ha='center', va='center',
        fontsize=10, style='italic', rotation=-60,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='lightcyan', alpha=0.8))

# tAI vs Safety edge
mid_tai_safety = (triangle_points[1] + triangle_points[2]) / 2
ax1.text(mid_tai_safety[0], mid_tai_safety[1] - 0.8, 'tAI-safety trade-off', ha='center', va='center',
        fontsize=10, style='italic',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.8))

# Add example codons at different positions
# High μ, low tAI (bottom right area)
ax1.scatter(7, 3.5, s=300, c=color_mu, marker='o', edgecolors='black', linewidth=2, zorder=5)
ax1.text(7, 3.5, 'TGT', ha='center', va='center', fontsize=9, fontweight='bold', color='white')
ax1.text(7.5, 3.5, '(high μ\nlow tAI)', ha='left', va='center', fontsize=8, style='italic')

# Low μ, high tAI (left area)
ax1.scatter(3, 4.5, s=300, c=color_tai, marker='o', edgecolors='black', linewidth=2, zorder=5)
ax1.text(3, 4.5, 'TGC', ha='center', va='center', fontsize=9, fontweight='bold', color='white')
ax1.text(2.5, 4.5, '(low μ\nhigh tAI)', ha='right', va='center', fontsize=8, style='italic')

# Central position (balanced)
ax1.scatter(5, 5, s=300, c='#FFCC00', marker='o', edgecolors='black', linewidth=2, zorder=5)
ax1.text(5, 5, 'GAC', ha='center', va='center', fontsize=9, fontweight='bold', color='black')
ax1.text(5, 4.2, '(balanced)', ha='center', va='top', fontsize=8, style='italic')

# Title and caption
ax1.text(5, 10.2, 'A. Proteostasis Trade-off Triangle', ha='center', va='bottom',
        fontsize=13, fontweight='bold')

caption = 'Codons occupy different positions\nin the space of competing constraints'
ax1.text(5, 0.2, caption, ha='center', va='bottom', fontsize=9, style='italic',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='lightyellow', edgecolor='orange', alpha=0.9))

# ============================================================================
# Panel 1B: κ as Gain Knob
# ============================================================================
ax2 = plt.subplot(1, 3, 2)
ax2.set_xlim(0, 10)
ax2.set_ylim(0, 10)
ax2.axis('off')

# Draw two scenarios side by side

# Left side: Low κ (κ ≈ 0) - buried, tolerant site
x_left = 2
y_center = 5

# Draw site box
buried_box = FancyBboxPatch((x_left - 1, y_center - 1.5), 2, 3,
                           boxstyle="round,pad=0.1",
                           facecolor='lightblue', edgecolor='blue', linewidth=2, alpha=0.3)
ax2.add_patch(buried_box)

ax2.text(x_left, y_center + 2, 'Low κ ≈ 0', ha='center', va='bottom',
        fontsize=11, fontweight='bold', color='blue')
ax2.text(x_left, y_center + 1.5, '(Buried/Tolerant)', ha='center', va='bottom',
        fontsize=9, style='italic', color='blue')

# Draw signal
ax2.plot([x_left - 0.5, x_left - 0.5], [y_center - 1, y_center + 1],
        'k-', linewidth=2, label='Input: J(c)')
ax2.text(x_left - 0.8, y_center - 1.5, 'J(c)', ha='center', fontsize=10, fontweight='bold')

# Amplifier symbol (low gain)
ax2.add_patch(Polygon([[x_left - 0.3, y_center - 0.8],
                       [x_left + 0.3, y_center],
                       [x_left - 0.3, y_center + 0.8]],
                      facecolor='lightgray', edgecolor='black', linewidth=2))
ax2.text(x_left, y_center, 'κ≈0', ha='center', va='center', fontsize=9, fontweight='bold')

# Output signal (damped)
ax2.plot([x_left + 0.5, x_left + 0.5], [y_center - 0.3, y_center + 0.3],
        'r-', linewidth=3, alpha=0.5)
ax2.text(x_left + 0.8, y_center - 0.8, 'κ·J(c) ≈ 0', ha='center', fontsize=10,
        fontweight='bold', color='blue')
ax2.text(x_left, y_center - 2.5, 'Codon choice\nDOESN\'T MATTER', ha='center', fontsize=9,
        fontweight='bold', color='blue')

# Right side: High κ (κ ≈ 1) - exposed, critical site
x_right = 7.5
y_center = 5

# Draw site box
exposed_box = FancyBboxPatch((x_right - 1, y_center - 1.5), 2, 3,
                            boxstyle="round,pad=0.1",
                            facecolor='#FFE6E6', edgecolor='red', linewidth=2, alpha=0.3)
ax2.add_patch(exposed_box)

ax2.text(x_right, y_center + 2, 'High κ ≈ 1', ha='center', va='bottom',
        fontsize=11, fontweight='bold', color='red')
ax2.text(x_right, y_center + 1.5, '(Exposed/Critical)', ha='center', va='bottom',
        fontsize=9, style='italic', color='red')

# Draw signal
ax2.plot([x_right - 0.5, x_right - 0.5], [y_center - 1, y_center + 1],
        'k-', linewidth=2)
ax2.text(x_right - 0.8, y_center - 1.5, 'J(c)', ha='center', fontsize=10, fontweight='bold')

# Amplifier symbol (high gain)
ax2.add_patch(Polygon([[x_right - 0.3, y_center - 0.8],
                       [x_right + 0.3, y_center],
                       [x_right - 0.3, y_center + 0.8]],
                      facecolor='#FFE6E6', edgecolor='red', linewidth=2))
ax2.text(x_right, y_center, 'κ≈1', ha='center', va='center', fontsize=9, fontweight='bold')

# Output signal (amplified)
ax2.plot([x_right + 0.5, x_right + 0.5], [y_center - 1, y_center + 1],
        'r-', linewidth=3)
ax2.text(x_right + 0.8, y_center - 1.5, '≈ J(c)', ha='center', fontsize=10,
        fontweight='bold', color='red')
ax2.text(x_right, y_center - 2.5, 'Codon choice\nMATTERS!', ha='center', fontsize=9,
        fontweight='bold', color='red')

# Title
ax2.text(5, 9.5, 'B. Site Criticality (κ) as Gain Control', ha='center', va='bottom',
        fontsize=13, fontweight='bold')

# Caption
caption = 'κ amplifies or dampens codon cost J(c)\nHigh-κ sites drive selection for optimal codons'
ax2.text(5, 0.2, caption, ha='center', va='bottom', fontsize=9, style='italic',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='#FFE6E6', edgecolor='red', alpha=0.9))

# ============================================================================
# Panel 1C: μ-tAI Modes per Amino Acid
# ============================================================================
ax3 = plt.subplot(1, 3, 3)
ax3.set_xlim(0, 10)
ax3.set_ylim(0, 10)
ax3.axis('off')

# Draw quadrant structure
center_x, center_y = 5, 5
quad_size = 3.5

# Draw axes
ax3.plot([center_x - quad_size, center_x + quad_size], [center_y, center_y],
        'k-', linewidth=2, alpha=0.5)
ax3.plot([center_x, center_x], [center_y - quad_size, center_y + quad_size],
        'k-', linewidth=2, alpha=0.5)

# Add axis labels
ax3.text(center_x + quad_size + 0.3, center_y, 'μ →', ha='left', va='center',
        fontsize=11, fontweight='bold')
ax3.text(center_x, center_y + quad_size + 0.3, 'tAI →', ha='center', va='bottom',
        fontsize=11, fontweight='bold')

# Draw four modes as colored regions with example codons

# Top-left: Safe Sprinter (low μ, high tAI)
rect_ss = Rectangle((center_x - quad_size, center_y), quad_size, quad_size,
                    facecolor=mode_colors['safe_sprinter'], alpha=0.3,
                    edgecolor=mode_colors['safe_sprinter'], linewidth=2)
ax3.add_patch(rect_ss)
ax3.text(center_x - quad_size/2, center_y + quad_size*0.7, 'Safe\nSprinter',
        ha='center', va='center', fontsize=11, fontweight='bold',
        color=mode_colors['safe_sprinter'])
ax3.scatter(center_x - quad_size/2, center_y + quad_size/3, s=200,
           c='white', marker='o', edgecolors='black', linewidth=2)
ax3.text(center_x - quad_size/2, center_y + quad_size/3, 'TGC',
        ha='center', va='center', fontsize=9, fontweight='bold')

# Top-right: Risky Sprinter (high μ, high tAI)
rect_rs = Rectangle((center_x, center_y), quad_size, quad_size,
                    facecolor=mode_colors['risky_sprinter'], alpha=0.3,
                    edgecolor=mode_colors['risky_sprinter'], linewidth=2)
ax3.add_patch(rect_rs)
ax3.text(center_x + quad_size/2, center_y + quad_size*0.7, 'Risky\nSprinter',
        ha='center', va='center', fontsize=11, fontweight='bold',
        color=mode_colors['risky_sprinter'])
ax3.scatter(center_x + quad_size/2, center_y + quad_size/3, s=200,
           c='white', marker='o', edgecolors='black', linewidth=2)
ax3.text(center_x + quad_size/2, center_y + quad_size/3, 'GAC',
        ha='center', va='center', fontsize=9, fontweight='bold')

# Bottom-left: Safe Careful (low μ, low tAI)
rect_sc = Rectangle((center_x - quad_size, center_y - quad_size), quad_size, quad_size,
                    facecolor=mode_colors['safe_careful'], alpha=0.3,
                    edgecolor=mode_colors['safe_careful'], linewidth=2)
ax3.add_patch(rect_sc)
ax3.text(center_x - quad_size/2, center_y - quad_size*0.3, 'Safe\nCareful',
        ha='center', va='center', fontsize=11, fontweight='bold',
        color=mode_colors['safe_careful'])
ax3.scatter(center_x - quad_size/2, center_y - quad_size*0.65, s=200,
           c='white', marker='o', edgecolors='black', linewidth=2)
ax3.text(center_x - quad_size/2, center_y - quad_size*0.65, 'GAT',
        ha='center', va='center', fontsize=9, fontweight='bold')

# Bottom-right: Risky Careful (high μ, low tAI)
rect_rc = Rectangle((center_x, center_y - quad_size), quad_size, quad_size,
                    facecolor=mode_colors['risky_careful'], alpha=0.3,
                    edgecolor=mode_colors['risky_careful'], linewidth=2)
ax3.add_patch(rect_rc)
ax3.text(center_x + quad_size/2, center_y - quad_size*0.3, 'Risky\nCareful',
        ha='center', va='center', fontsize=11, fontweight='bold',
        color=mode_colors['risky_careful'])
ax3.scatter(center_x + quad_size/2, center_y - quad_size*0.65, s=200,
           c='white', marker='o', edgecolors='black', linewidth=2)
ax3.text(center_x + quad_size/2, center_y - quad_size*0.65, 'TGT',
        ha='center', va='center', fontsize=9, fontweight='bold')

# Title
ax3.text(5, 9.5, 'C. Four Proteostasis Modes', ha='center', va='bottom',
        fontsize=13, fontweight='bold')

# Caption
caption = ('Each amino acid\'s codons span multiple modes\n'
          'Different functional contexts select different modes')
ax3.text(5, 0.2, caption, ha='center', va='bottom', fontsize=9, style='italic',
        bbox=dict(boxstyle='round,pad=0.5', facecolor='lightcyan', edgecolor='blue', alpha=0.9))

# Add annotation showing same AA (Cys) has multiple modes
ax3.annotate('', xy=(center_x - quad_size/2, center_y + quad_size/3),
            xytext=(center_x + quad_size/2, center_y - quad_size*0.65),
            arrowprops=dict(arrowstyle='<->', lw=2, color='purple', linestyle='--'))
ax3.text(center_x + 1.5, center_y - 0.5, 'Same AA (Cys)\ndifferent modes',
        ha='center', va='center', fontsize=8, color='purple', fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='purple', linewidth=1.5))

plt.tight_layout()
plt.savefig('figure1_conceptual.png', dpi=300, bbox_inches='tight')
plt.savefig('figure1_conceptual.pdf', bbox_inches='tight')
plt.savefig('figure1_conceptual.svg', bbox_inches='tight')

print("Figure 1 (conceptual schematic) saved as figure1_conceptual.png and .pdf")
print("\nPanel A: Trade-off triangle showing competing constraints (μ, tAI, safety)")
print("Panel B: κ as gain knob amplifying or damping codon cost")
print("Panel C: Four proteostasis modes in μ-tAI space with example codons")
