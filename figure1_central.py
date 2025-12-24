#!/usr/bin/env python3
"""
Figure 1: Central conceptual schematic of the Proteostatic Codon Law
Three horizontal panels (A-C) with consistent visual language
Panel A: Proteostatic trade-off triangle
Panel B: κ as gain knob
Panel C: μ-tAI modes per amino acid
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, Circle, Polygon, Rectangle, Wedge
from matplotlib.patches import FancyArrowPatch
import matplotlib.lines as mlines

# Set style - minimal, clean
plt.rcParams['font.size'] = 9
plt.rcParams['axes.labelsize'] = 10
plt.rcParams['axes.titlesize'] = 11
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

# Color scheme - single accent color (dark blue) + light grey
accent_color = '#2E5090'  # Dark blue
light_grey = '#CCCCCC'
very_light_grey = '#F0F0F0'

# Mode colors for Panel C
mode_color_cys = '#00A087'  # Green
mode_color_asp = '#4DBBD5'  # Blue

# Create figure - 16:9 aspect ratio
fig = plt.figure(figsize=(16, 9))

# ============================================================================
# Panel A: Proteostatic Trade-off Triangle
# ============================================================================
ax1 = plt.subplot(1, 3, 1)
ax1.set_xlim(-0.2, 1.2)
ax1.set_ylim(-0.2, 1.2)
ax1.set_aspect('equal')
ax1.axis('off')

# Panel label
ax1.text(0.05, 1.15, 'A', fontsize=16, fontweight='bold', ha='left', va='top')

# Draw equilateral triangle
triangle_height = 0.866  # sqrt(3)/2
triangle_points = np.array([
    [0.5, 0.9],      # Top vertex
    [0.1, 0.1],      # Bottom left
    [0.9, 0.1]       # Bottom right
])

triangle = Polygon(triangle_points, fill=False, edgecolor=light_grey, linewidth=1.5)
ax1.add_patch(triangle)

# Vertex circles and labels with different shades
vertex_colors = ['#4A6FA5', '#6B8FC2', '#8CA8D8']

# Top vertex - Accuracy (μ)
ax1.scatter(triangle_points[0, 0], triangle_points[0, 1], s=200,
           c=vertex_colors[0], edgecolors='black', linewidth=1, zorder=5)
ax1.text(triangle_points[0, 0], triangle_points[0, 1] + 0.08,
        'Accuracy (μ)', ha='center', va='bottom', fontsize=10, fontweight='bold')
ax1.text(triangle_points[0, 0], triangle_points[0, 1] + 0.04,
        'Misincorporation', ha='center', va='bottom', fontsize=8, style='italic', color='grey')

# Bottom-left vertex - Decoding performance (tAI)
ax1.scatter(triangle_points[1, 0], triangle_points[1, 1], s=200,
           c=vertex_colors[1], edgecolors='black', linewidth=1, zorder=5)
ax1.text(triangle_points[1, 0], triangle_points[1, 1] - 0.08,
        'Decoding performance (tAI)', ha='center', va='top', fontsize=10, fontweight='bold')
ax1.text(triangle_points[1, 0], triangle_points[1, 1] - 0.04,
        'Speed / throughput', ha='center', va='top', fontsize=8, style='italic', color='grey')

# Bottom-right vertex - Pairing safety
ax1.scatter(triangle_points[2, 0], triangle_points[2, 1], s=200,
           c=vertex_colors[2], edgecolors='black', linewidth=1, zorder=5)
ax1.text(triangle_points[2, 0], triangle_points[2, 1] - 0.08,
        'Pairing safety', ha='center', va='top', fontsize=10, fontweight='bold')
ax1.text(triangle_points[2, 0], triangle_points[2, 1] - 0.04,
        'Watson–Crick vs wobble', ha='center', va='top', fontsize=8, style='italic', color='grey')

# Codon point (slightly off-center)
codon_point = np.array([0.5, 0.45])
ax1.scatter(codon_point[0], codon_point[1], s=150, c=accent_color,
           marker='o', edgecolors='black', linewidth=1.5, zorder=6)
ax1.text(codon_point[0] + 0.08, codon_point[1], 'codon c',
        fontsize=9, style='italic', ha='left', va='center')

# Barycentric arrows (faint)
for vertex in triangle_points:
    ax1.annotate('', xy=codon_point, xytext=vertex,
                arrowprops=dict(arrowstyle='->', lw=0.8, color=light_grey, alpha=0.5))

# Equation at center
equation = r'$\Delta J_i^{(0)}(c) = \alpha \mu(c) + \beta f(a(c)) + \gamma g(b(c))$'
ax1.text(codon_point[0], codon_point[1] - 0.12, equation,
        fontsize=9, ha='center', va='top',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor=light_grey, linewidth=1))

# Intuition label at bottom
intuition = ('Each codon c has a fixed proteostatic profile:\n'
            'accuracy μ, decoding performance (tAI), and pairing safety.')
ax1.text(0.5, -0.15, intuition, ha='center', va='top', fontsize=8.5,
        style='italic', color='#555555')

# ============================================================================
# Panel B: κ as Gain Knob
# ============================================================================
ax2 = plt.subplot(1, 3, 2)
ax2.set_xlim(0, 1)
ax2.set_ylim(0, 1)
ax2.axis('off')

# Panel label
ax2.text(0.05, 0.98, 'B', fontsize=16, fontweight='bold', ha='left', va='top')

# --- Top mini-panel: High κ ---
y_top = 0.75

# Label
ax2.text(0.5, y_top + 0.18, 'High-κ site (metal ligand / active site)',
        ha='center', va='bottom', fontsize=9, fontweight='bold')

# Slider/dial
slider_x = 0.15
slider_height = 0.25
# Draw slider track
ax2.plot([slider_x, slider_x], [y_top - 0.05, y_top + slider_height],
        'k-', linewidth=2, alpha=0.3)
# Draw slider knob at high position (0.9)
knob_y = y_top + slider_height * 0.9
ax2.scatter(slider_x, knob_y, s=300, c=accent_color, marker='s',
           edgecolors='black', linewidth=2, zorder=5)
ax2.text(slider_x - 0.08, knob_y, r'κ$_i$ ≈ 1', ha='right', va='center',
        fontsize=9, fontweight='bold')

# Labels on slider
ax2.text(slider_x + 0.02, y_top + slider_height, '1', ha='left', va='center', fontsize=8)
ax2.text(slider_x + 0.02, y_top - 0.05, '0', ha='left', va='center', fontsize=8)

# J bar (tall)
bar_x = 0.35
bar_width = 0.08
bar_height_high = 0.23
bar_rect_high = Rectangle((bar_x, y_top), bar_width, bar_height_high,
                         facecolor=accent_color, edgecolor='black', linewidth=1.5)
ax2.add_patch(bar_rect_high)
ax2.text(bar_x + bar_width/2, y_top + bar_height_high/2, r'$J_i(c)$',
        ha='center', va='center', fontsize=9, fontweight='bold', color='white')

# Formula
ax2.text(bar_x + bar_width + 0.08, y_top + bar_height_high/2,
        r'$J_i(c) = \kappa_i \Delta J_i^{(0)}(c)$',
        ha='left', va='center', fontsize=9)

# --- Bottom mini-panel: Low κ ---
y_bottom = 0.25

# Label
ax2.text(0.5, y_bottom + 0.18, 'Low-κ site (solvent-exposed / flexible loop)',
        ha='center', va='bottom', fontsize=9, fontweight='bold')

# Slider/dial
# Draw slider track
ax2.plot([slider_x, slider_x], [y_bottom - 0.05, y_bottom + slider_height],
        'k-', linewidth=2, alpha=0.3)
# Draw slider knob at low position (0.1)
knob_y_low = y_bottom + slider_height * 0.1
ax2.scatter(slider_x, knob_y_low, s=300, c='#CCCCCC', marker='s',
           edgecolors='black', linewidth=2, zorder=5)
ax2.text(slider_x - 0.08, knob_y_low, r'κ$_i$ ≈ 0', ha='right', va='center',
        fontsize=9, fontweight='bold')

# Labels on slider
ax2.text(slider_x + 0.02, y_bottom + slider_height, '1', ha='left', va='center', fontsize=8)
ax2.text(slider_x + 0.02, y_bottom - 0.05, '0', ha='left', va='center', fontsize=8)

# J bar (short, pale)
bar_height_low = 0.03
bar_rect_low = Rectangle((bar_x, y_bottom), bar_width, bar_height_low,
                         facecolor='#E0E0E0', edgecolor='black', linewidth=1.5)
ax2.add_patch(bar_rect_low)
ax2.text(bar_x + bar_width/2, y_bottom + bar_height_low/2, r'$J_i(c)$',
        ha='center', va='center', fontsize=7)

# Formula
ax2.text(bar_x + bar_width + 0.08, y_bottom + bar_height_low/2,
        r'$J_i(c) = \kappa_i \Delta J_i^{(0)}(c)$',
        ha='left', va='center', fontsize=9)

# Note about same codon
ax2.text(0.5, 0.52, '(same μ, tAI, pairing; different κ$_i$)',
        ha='center', va='center', fontsize=8, style='italic', color='#666666')

# Right side annotation
annotation_x = 0.6
annotation_text = (
    r'κ$_i$ weights how much local errors' + '\n' +
    'matter for proteostasis.\n\n' +
    'The same codon c produces orders-of-\n' +
    'magnitude higher burden at κ≈1 sites\n' +
    'than at κ≈0 sites.'
)
ax2.text(annotation_x, 0.5, annotation_text, ha='left', va='center',
        fontsize=8.5, style='italic', color='#555555',
        bbox=dict(boxstyle='round,pad=0.5', facecolor=very_light_grey,
                 edgecolor=light_grey, linewidth=1))

# ============================================================================
# Panel C: μ-tAI Modes per Amino Acid
# ============================================================================
ax3 = plt.subplot(1, 3, 3)
ax3.set_xlim(-0.1, 1.1)
ax3.set_ylim(-0.1, 1.1)
ax3.set_aspect('equal')

# Panel label
ax3.text(-0.05, 1.05, 'C', fontsize=16, fontweight='bold', ha='left', va='top',
        transform=ax3.transData)

# Axes
ax3.axhline(0.5, color=light_grey, linewidth=1.5, zorder=1)
ax3.axvline(0.5, color=light_grey, linewidth=1.5, zorder=1)

ax3.set_xlabel('μ(c) →', fontsize=10, fontweight='bold')
ax3.set_ylabel('tAI(c) →', fontsize=10, fontweight='bold')
ax3.tick_params(left=False, bottom=False, labelleft=False, labelbottom=False)
ax3.spines['top'].set_visible(False)
ax3.spines['right'].set_visible(False)
ax3.spines['bottom'].set_visible(False)
ax3.spines['left'].set_visible(False)

# Quadrant labels
label_offset_x = 0.25
label_offset_y = 0.75

# Top-left: safe sprinter
ax3.text(label_offset_x, label_offset_y, 'safe sprinter\n(low μ, high tAI)',
        ha='center', va='center', fontsize=9, style='italic', color='#555555')

# Top-right: risky sprinter
ax3.text(1 - label_offset_x, label_offset_y, 'risky sprinter\n(high μ, high tAI)',
        ha='center', va='center', fontsize=9, style='italic', color='#555555')

# Bottom-left: safe careful
ax3.text(label_offset_x, 1 - label_offset_y, 'safe careful\n(low μ, low tAI)',
        ha='center', va='center', fontsize=9, style='italic', color='#555555')

# Bottom-right: risky careful
ax3.text(1 - label_offset_x, 1 - label_offset_y, 'risky careful\n(high μ, low tAI)',
        ha='center', va='center', fontsize=9, style='italic', color='#555555')

# Cys example (green)
# TGC: safe sprinter (low μ, high tAI)
tgc_x, tgc_y = 0.2, 0.75
# TGT: risky careful (high μ, low tAI)
tgt_x, tgt_y = 0.85, 0.2

ax3.scatter(tgc_x, tgc_y, s=250, c=mode_color_cys, marker='o',
           edgecolors='black', linewidth=2, zorder=5)
ax3.text(tgc_x, tgc_y, 'TGC', ha='center', va='center',
        fontsize=9, fontweight='bold', color='white')

ax3.scatter(tgt_x, tgt_y, s=250, c=mode_color_cys, marker='o',
           edgecolors='black', linewidth=2, zorder=5)
ax3.text(tgt_x, tgt_y, 'TGT', ha='center', va='center',
        fontsize=9, fontweight='bold', color='white')

# Enclosing bracket for Cys
from matplotlib.patches import FancyBboxPatch
cys_box = FancyBboxPatch((tgc_x - 0.1, tgt_y - 0.05),
                         (tgt_x - tgc_x + 0.2), (tgc_y - tgt_y + 0.1),
                         boxstyle="round,pad=0.05",
                         fill=False, edgecolor=mode_color_cys,
                         linewidth=2, linestyle='--', alpha=0.6)
ax3.add_patch(cys_box)
ax3.text((tgc_x + tgt_x)/2, tgt_y - 0.12, 'Cys modes',
        ha='center', va='top', fontsize=9, fontweight='bold', color=mode_color_cys)

# Asp example (blue)
# GAC: risky sprinter (high μ, high tAI)
gac_x, gac_y = 0.75, 0.8
# GAT: safe careful (low μ, low tAI)
gat_x, gat_y = 0.15, 0.25

ax3.scatter(gac_x, gac_y, s=250, c=mode_color_asp, marker='o',
           edgecolors='black', linewidth=2, zorder=5)
ax3.text(gac_x, gac_y, 'GAC', ha='center', va='center',
        fontsize=9, fontweight='bold', color='white')

ax3.scatter(gat_x, gat_y, s=250, c=mode_color_asp, marker='o',
           edgecolors='black', linewidth=2, zorder=5)
ax3.text(gat_x, gat_y, 'GAT', ha='center', va='center',
        fontsize=9, fontweight='bold', color='white')

# Enclosing bracket for Asp
asp_box = FancyBboxPatch((gat_x - 0.05, gat_y - 0.05),
                         (gac_x - gat_x + 0.1), (gac_y - gat_y + 0.1),
                         boxstyle="round,pad=0.05",
                         fill=False, edgecolor=mode_color_asp,
                         linewidth=2, linestyle='--', alpha=0.6)
ax3.add_patch(asp_box)
ax3.text((gac_x + gat_x)/2 + 0.12, gac_y + 0.08, 'Asp modes',
        ha='center', va='bottom', fontsize=9, fontweight='bold', color=mode_color_asp)

# Mode count summary in corner
summary_text = 'AA   modes (k$_A$)\nCys      2\nAsp      2\n…         …'
ax3.text(0.98, 0.02, summary_text, ha='right', va='bottom',
        fontsize=8, family='monospace',
        transform=ax3.transAxes,
        bbox=dict(boxstyle='round,pad=0.5', facecolor='white',
                 edgecolor=light_grey, linewidth=1))

# Caption
caption = ('Each amino acid uses a small set of μ–tAI modes;\n'
          'the code must supply enough codons to realize them.')
ax3.text(0.5, -0.05, caption, ha='center', va='top', fontsize=8.5,
        style='italic', color='#555555', transform=ax3.transAxes)

# Overall title
fig.suptitle('The Proteostatic Codon Law', fontsize=14, fontweight='bold', y=0.98)

plt.tight_layout(rect=[0, 0.02, 1, 0.96])
plt.savefig('figure1_central.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig('figure1_central.pdf', bbox_inches='tight', facecolor='white')
plt.savefig('figure1_central.svg', bbox_inches='tight', facecolor='white')

print("Figure 1 (central/redesigned) saved as:")
print("  - figure1_central.png (300 DPI raster)")
print("  - figure1_central.pdf (vector, fully editable)")
print("  - figure1_central.svg (vector, web/Inkscape friendly)")
print("\nPanel A: Proteostatic trade-off triangle with barycentric formula")
print("Panel B: κ as gain knob (dual slider/bar comparison)")
print("Panel C: μ-tAI quadrant with Cys and Asp mode examples")
