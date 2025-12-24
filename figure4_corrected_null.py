#!/usr/bin/env python3
"""
Figure 4: Operational Diversity Structure

CORRECTED NULL CONSTRUCTION:
- Shuffles codon→AA assignments (which codons belong to which amino acid)
- Keeps codon properties (μ, tAI) FIXED for each codon
- Preserves degeneracy structure (9×2-codon, 1×3-codon, 5×4-codon, 3×6-codon)

This is the correct null hypothesis: "Does the SPECIFIC assignment of codons
to amino acids matter, or would any random assignment give similar diversity?"
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from collections import defaultdict
import random

# Style settings
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.linewidth'] = 1.2
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False

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

def load_codon_data(path='errors/codon_modes_ecoli.tsv'):
    """Load codon data with mu and tAI values."""
    df = pd.read_csv(path, sep='\t')
    df['mu'] = pd.to_numeric(df['mu'], errors='coerce')
    df['w_tai'] = pd.to_numeric(df['w_tai'], errors='coerce')
    return df


def compute_operational_diversity(codon_to_aa, codon_properties, norm_params):
    """
    Compute operational diversity Δ_A for each amino acid.
    
    Parameters:
    -----------
    codon_to_aa : dict
        Mapping from codon to amino acid
    codon_properties : dict
        Mapping from codon to (mu, tAI) tuple
    norm_params : tuple
        (mu_mean, mu_std, tai_mean, tai_std) for standardization
    
    Returns:
    --------
    dict : aa -> Δ_A (mean pairwise distance in standardized space)
    """
    mu_mean, mu_std, tai_mean, tai_std = norm_params
    
    # Group codons by amino acid
    aa_to_codons = defaultdict(list)
    for codon, aa in codon_to_aa.items():
        if codon in codon_properties:
            aa_to_codons[aa].append(codon)
    
    delta_by_aa = {}
    
    for aa, codons in aa_to_codons.items():
        if len(codons) < 2:
            delta_by_aa[aa] = 0.0
            continue
        
        # Get standardized coordinates for each codon
        coords = []
        for codon in codons:
            mu, tai = codon_properties[codon]
            if mu is None or tai is None or np.isnan(mu) or np.isnan(tai):
                continue
            # Standardize
            mu_norm = (np.log10(mu) - mu_mean) / mu_std
            tai_norm = (tai - tai_mean) / tai_std
            coords.append((mu_norm, tai_norm))
        
        if len(coords) < 2:
            delta_by_aa[aa] = 0.0
            continue
        
        # Compute mean pairwise Euclidean distance
        total_dist = 0.0
        count = 0
        for i in range(len(coords)):
            for j in range(i + 1, len(coords)):
                dist = np.sqrt((coords[i][0] - coords[j][0])**2 + 
                              (coords[i][1] - coords[j][1])**2)
                total_dist += dist
                count += 1
        
        delta_by_aa[aa] = total_dist / count if count > 0 else 0.0
    
    return delta_by_aa


def generate_null_ensemble_correct(df, n_null=10000, seed=42):
    """
    CORRECT null ensemble generation.
    
    Shuffles codon→AA assignments while:
    - Keeping codon properties (μ, tAI) FIXED for each codon
    - Preserving degeneracy structure
    
    This asks: "Given the distribution of codon properties across all 61 codons,
    does the specific assignment to amino acids produce unusual diversity?"
    """
    random.seed(seed)
    np.random.seed(seed)
    
    # Get valid codons with both mu and tAI
    valid_df = df.dropna(subset=['mu', 'w_tai']).copy()
    
    # Build codon→properties mapping (FIXED throughout)
    codon_properties = {}
    for _, row in valid_df.iterrows():
        codon_properties[row['codon']] = (row['mu'], row['w_tai'])
    
    # Build observed codon→AA mapping
    observed_codon_to_aa = {}
    for _, row in valid_df.iterrows():
        observed_codon_to_aa[row['codon']] = row['aa']
    
    # Get degeneracy structure from observed code
    aa_degeneracy = defaultdict(int)
    for aa in observed_codon_to_aa.values():
        aa_degeneracy[aa] += 1
    
    # Sort AAs by name for consistent ordering
    aa_list = sorted(aa_degeneracy.keys())
    degeneracy_list = [aa_degeneracy[aa] for aa in aa_list]
    
    # Compute normalization parameters from ALL codons (fixed)
    all_mu = [codon_properties[c][0] for c in codon_properties]
    all_tai = [codon_properties[c][1] for c in codon_properties]
    
    log_mu = np.log10(all_mu)
    mu_mean, mu_std = np.mean(log_mu), np.std(log_mu)
    tai_mean, tai_std = np.mean(all_tai), np.std(all_tai)
    norm_params = (mu_mean, mu_std, tai_mean, tai_std)
    
    # Compute OBSERVED diversity
    observed_delta = compute_operational_diversity(
        observed_codon_to_aa, codon_properties, norm_params
    )
    observed_total = sum(observed_delta.values())
    
    # Generate null distribution
    all_codons = list(codon_properties.keys())
    null_totals = []
    
    for iteration in range(n_null):
        # Shuffle codons
        shuffled_codons = all_codons.copy()
        random.shuffle(shuffled_codons)
        
        # Assign to AAs preserving degeneracy
        null_codon_to_aa = {}
        idx = 0
        for aa, deg in zip(aa_list, degeneracy_list):
            for _ in range(deg):
                if idx < len(shuffled_codons):
                    null_codon_to_aa[shuffled_codons[idx]] = aa
                    idx += 1
        
        # Compute diversity for this null code
        null_delta = compute_operational_diversity(
            null_codon_to_aa, codon_properties, norm_params
        )
        null_totals.append(sum(null_delta.values()))
    
    return observed_total, null_totals, observed_delta, norm_params


def get_wobble_type(codon1, codon2):
    """Determine if two synonymous codons differ at purine (A/G) or pyrimidine (C/T) position."""
    # Find differing position
    for i in range(3):
        if codon1[i] != codon2[i]:
            nts = {codon1[i], codon2[i]}
            if nts == {'A', 'G'}:
                return 'purine'
            elif nts == {'C', 'T'} or nts == {'U', 'T'}:
                return 'pyrimidine'
    return 'other'


def analyze_wobble_pattern(df, observed_delta):
    """Analyze the purine vs pyrimidine pattern for 2-codon families."""
    # Get 2-codon families
    degeneracy = df.groupby('aa').size().to_dict()
    two_codon_aas = [aa for aa, d in degeneracy.items() if d == 2]
    
    purine_spreads = []
    pyrimidine_spreads = []
    aa_wobble_type = {}
    
    for aa in two_codon_aas:
        codons = df[df['aa'] == aa]['codon'].tolist()
        if len(codons) != 2:
            continue
        
        # Determine wobble type
        wobble = get_wobble_type(codons[0], codons[1])
        aa_wobble_type[aa] = wobble
        
        spread = observed_delta.get(aa, 0)
        
        if wobble == 'purine':
            purine_spreads.append((aa, spread))
        elif wobble == 'pyrimidine':
            pyrimidine_spreads.append((aa, spread))
    
    return purine_spreads, pyrimidine_spreads, aa_wobble_type


def create_figure4_corrected(data_path='errors/codon_modes_ecoli.tsv', 
                              output_prefix='Figure4_diversity_CORRECTED',
                              n_null=10000):
    """
    Create Figure 4 with CORRECT null construction.
    
    Panel A: Observed vs null distribution
    Panel B: Wobble pattern (purine vs pyrimidine)
    Panel C: Diversity by amino acid
    """
    print("=" * 70)
    print("GENERATING FIGURE 4 WITH CORRECT NULL CONSTRUCTION")
    print("=" * 70)
    print(f"\nNull construction: Shuffle codon→AA assignments")
    print(f"                   Keep codon properties (μ, tAI) FIXED")
    print(f"                   Preserve degeneracy structure")
    print(f"\nNumber of null samples: {n_null}")
    
    # Load data
    df = load_codon_data(data_path)
    print(f"\nLoaded {len(df)} codons from {data_path}")
    
    valid_df = df.dropna(subset=['mu', 'w_tai'])
    print(f"Valid codons (with both μ and tAI): {len(valid_df)}")
    
    # Generate null ensemble
    print(f"\nGenerating null ensemble...")
    observed_total, null_totals, observed_delta, norm_params = \
        generate_null_ensemble_correct(df, n_null=n_null)
    
    null_arr = np.array(null_totals)
    null_mean = np.mean(null_arr)
    null_std = np.std(null_arr)
    z_score = (observed_total - null_mean) / null_std
    
    # Compute p-value (one-tailed, testing if observed < null)
    p_value_lower = np.sum(null_arr <= observed_total) / len(null_arr)
    p_value_upper = np.sum(null_arr >= observed_total) / len(null_arr)
    p_value = min(p_value_lower, p_value_upper)  # One-tailed
    
    print(f"\n{'='*50}")
    print("RESULTS:")
    print(f"{'='*50}")
    print(f"  Observed total diversity: {observed_total:.2f}")
    print(f"  Null mean ± std: {null_mean:.2f} ± {null_std:.2f}")
    print(f"  z-score: {z_score:.2f}")
    print(f"  p-value (one-tailed, lower): {p_value_lower:.4f}")
    print(f"  p-value (one-tailed, upper): {p_value_upper:.4f}")
    print(f"{'='*50}")
    
    # Analyze wobble pattern
    purine_spreads, pyrimidine_spreads, aa_wobble_type = \
        analyze_wobble_pattern(df, observed_delta)
    
    purine_vals = [x[1] for x in purine_spreads]
    pyrimidine_vals = [x[1] for x in pyrimidine_spreads]
    
    purine_mean = np.mean(purine_vals) if purine_vals else 0
    pyrimidine_mean = np.mean(pyrimidine_vals) if pyrimidine_vals else 0
    
    # Statistical test for wobble difference
    if len(purine_vals) >= 2 and len(pyrimidine_vals) >= 2:
        t_stat, wobble_p = stats.ttest_ind(purine_vals, pyrimidine_vals)
        # Also do permutation test
        combined = purine_vals + pyrimidine_vals
        observed_diff = purine_mean - pyrimidine_mean
        perm_diffs = []
        for _ in range(10000):
            random.shuffle(combined)
            perm_purine = combined[:len(purine_vals)]
            perm_pyrimidine = combined[len(purine_vals):]
            perm_diffs.append(np.mean(perm_purine) - np.mean(perm_pyrimidine))
        wobble_p_perm = np.sum(np.abs(perm_diffs) >= np.abs(observed_diff)) / len(perm_diffs)
    else:
        wobble_p = np.nan
        wobble_p_perm = np.nan
    
    print(f"\nWOBBLE PATTERN (2-codon families):")
    print(f"  Purine pairs (A/G):     {purine_mean:.2f} (n={len(purine_vals)})")
    print(f"  Pyrimidine pairs (C/T): {pyrimidine_mean:.2f} (n={len(pyrimidine_vals)})")
    print(f"  Ratio: {purine_mean/pyrimidine_mean:.2f}x" if pyrimidine_mean > 0 else "  Ratio: N/A")
    print(f"  t-test p-value: {wobble_p:.4f}")
    print(f"  Permutation p-value: {wobble_p_perm:.4f}")
    
    # Create figure
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    
    # Panel A: Observed vs null distribution
    ax = axes[0]
    
    ax.hist(null_totals, bins=50, color=COLORS['light_gray'], edgecolor='white',
            alpha=0.8, density=True, label='Null distribution')
    
    ax.axvline(x=observed_total, color=COLORS['red'], linewidth=2.5, 
               label=f'Observed = {observed_total:.1f}')
    
    # Shade appropriate tail
    ylim = ax.get_ylim()
    if z_score < 0:
        ax.fill_betweenx([0, ylim[1]], null_arr.min() - 1, observed_total,
                         alpha=0.3, color=COLORS['red'])
    else:
        ax.fill_betweenx([0, ylim[1]], observed_total, null_arr.max() + 1,
                         alpha=0.3, color=COLORS['red'])
    
    ax.text(0.95, 0.95, f'z = {z_score:.2f}\np = {p_value:.4f}',
            transform=ax.transAxes, fontsize=10, ha='right', va='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8, edgecolor='gray'))
    
    ax.set_xlabel('Total operational diversity', fontsize=11)
    ax.set_ylabel('Density', fontsize=11)
    ax.legend(loc='upper left', fontsize=9, frameon=False)
    ax.set_title('Observed vs. null distribution', fontsize=11, pad=10)
    ax.text(-0.12, 1.05, 'A', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    
    # Panel B: Wobble pattern
    ax = axes[1]
    
    # Sort by spread within each group
    purine_sorted = sorted(purine_spreads, key=lambda x: -x[1])
    pyrimidine_sorted = sorted(pyrimidine_spreads, key=lambda x: -x[1])
    
    all_two_codon = purine_sorted + pyrimidine_sorted
    aa_labels = [x[0] for x in all_two_codon]
    spreads = [x[1] for x in all_two_codon]
    bar_colors = [COLORS['blue'] if aa_wobble_type[aa] == 'purine' else COLORS['orange'] 
                  for aa in aa_labels]
    
    x_pos = np.arange(len(aa_labels))
    bars = ax.bar(x_pos, spreads, color=bar_colors, edgecolor='white', linewidth=1)
    
    # Add wobble type annotations
    for i, (aa, spread) in enumerate(all_two_codon):
        wobble = aa_wobble_type[aa]
        label = 'A/G' if wobble == 'purine' else 'C/T'
        ax.text(i, spread + 0.1, label, ha='center', fontsize=7, color=COLORS['dark_gray'])
    
    # Add mean lines
    ax.axhline(purine_mean, color=COLORS['blue'], linestyle='--', linewidth=1.5, alpha=0.7)
    ax.axhline(pyrimidine_mean, color=COLORS['orange'], linestyle='--', linewidth=1.5, alpha=0.7)
    
    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor=COLORS['blue'], label=f'Purine pair (A/G): mean={purine_mean:.2f}'),
        Patch(facecolor=COLORS['orange'], label=f'Pyrimidine pair (C/T): mean={pyrimidine_mean:.2f}')
    ]
    ax.legend(handles=legend_elements, loc='upper right', fontsize=8, frameon=False)
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels(aa_labels, fontsize=9)
    ax.set_xlabel('Amino acid (2-codon families)', fontsize=11)
    ax.set_ylabel('Operational spread\n(standardized distance)', fontsize=11)
    ax.set_title('Wobble position determines operational spread', fontsize=11, pad=10)
    ax.text(-0.12, 1.05, 'B', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    
    # Panel C: All amino acids by degeneracy
    ax = axes[2]
    
    degeneracy = df.groupby('aa').size().to_dict()
    
    two_codon = sorted([aa for aa, d in degeneracy.items() if d == 2])
    three_codon = sorted([aa for aa, d in degeneracy.items() if d == 3])
    four_codon = sorted([aa for aa, d in degeneracy.items() if d == 4])
    six_codon = sorted([aa for aa, d in degeneracy.items() if d == 6])
    
    aa_order = two_codon + three_codon + four_codon + six_codon
    delta_vals = [observed_delta.get(aa, 0) for aa in aa_order]
    
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
    ax.bar(x_pos, delta_vals, color=colors, edgecolor='white', linewidth=1)
    
    # Add separators
    sep_positions = [len(two_codon) - 0.5,
                     len(two_codon) + len(three_codon) - 0.5,
                     len(two_codon) + len(three_codon) + len(four_codon) - 0.5]
    for pos in sep_positions:
        ax.axvline(x=pos, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)
    
    # Add class labels at top
    ymax = ax.get_ylim()[1]
    ax.text(len(two_codon)/2 - 0.5, ymax * 0.95, '2-codon',
            ha='center', fontsize=9, color=COLORS['blue'])
    ax.text(len(two_codon) + len(three_codon)/2 - 0.5, ymax * 0.95, '3',
            ha='center', fontsize=9, color=COLORS['purple'])
    ax.text(len(two_codon) + len(three_codon) + len(four_codon)/2 - 0.5,
            ymax * 0.95, '4-codon', ha='center', fontsize=9, color=COLORS['green'])
    ax.text(len(aa_order) - len(six_codon)/2 - 0.5, ymax * 0.95, '6-codon',
            ha='center', fontsize=9, color=COLORS['orange'])
    
    ax.set_xticks(x_pos)
    ax.set_xticklabels(aa_order, fontsize=8)
    ax.set_xlabel('Amino acid', fontsize=11)
    ax.set_ylabel(r'Mean pairwise distance $\Delta_A$', fontsize=11)
    ax.set_title('Operational diversity by amino acid', fontsize=11, pad=10)
    ax.text(-0.12, 1.05, 'C', transform=ax.transAxes, fontsize=14,
            fontweight='bold', va='top')
    
    plt.tight_layout()
    
    # Save
    for ext in ['png', 'pdf', 'svg']:
        plt.savefig(f'{output_prefix}.{ext}', dpi=300, bbox_inches='tight',
                    facecolor='white', edgecolor='none')
    plt.close()
    
    print(f"\nFigure saved as {output_prefix}.{{png,pdf,svg}}")
    
    # Return summary for verification
    return {
        'observed_total': observed_total,
        'null_mean': null_mean,
        'null_std': null_std,
        'z_score': z_score,
        'p_value': p_value,
        'purine_mean': purine_mean,
        'pyrimidine_mean': pyrimidine_mean,
        'wobble_ratio': purine_mean / pyrimidine_mean if pyrimidine_mean > 0 else np.nan,
        'wobble_p': wobble_p_perm,
        'observed_delta': observed_delta,
    }


def verify_null_construction(df, n_samples=5):
    """
    Verification: Print examples of null codes to confirm correct construction.
    """
    print("\n" + "=" * 70)
    print("VERIFICATION: Showing examples of null code construction")
    print("=" * 70)
    
    valid_df = df.dropna(subset=['mu', 'w_tai']).copy()
    
    # Observed assignment
    print("\nOBSERVED CODE (first 10 codons):")
    print("-" * 40)
    for _, row in valid_df.head(10).iterrows():
        print(f"  {row['codon']} → {row['aa']}  (μ={row['mu']:.4f}, tAI={row['w_tai']:.3f})")
    
    # Build structures
    codon_properties = {}
    for _, row in valid_df.iterrows():
        codon_properties[row['codon']] = (row['mu'], row['w_tai'])
    
    observed_codon_to_aa = {}
    for _, row in valid_df.iterrows():
        observed_codon_to_aa[row['codon']] = row['aa']
    
    aa_degeneracy = defaultdict(int)
    for aa in observed_codon_to_aa.values():
        aa_degeneracy[aa] += 1
    
    aa_list = sorted(aa_degeneracy.keys())
    degeneracy_list = [aa_degeneracy[aa] for aa in aa_list]
    
    all_codons = list(codon_properties.keys())
    
    print(f"\nDegeneracy structure: {dict(zip(aa_list, degeneracy_list))}")
    print(f"Total codons: {len(all_codons)}")
    
    # Generate a few null samples
    for i in range(n_samples):
        print(f"\nNULL SAMPLE {i+1} (first 10 codons):")
        print("-" * 40)
        
        shuffled_codons = all_codons.copy()
        random.shuffle(shuffled_codons)
        
        null_codon_to_aa = {}
        idx = 0
        for aa, deg in zip(aa_list, degeneracy_list):
            for _ in range(deg):
                if idx < len(shuffled_codons):
                    null_codon_to_aa[shuffled_codons[idx]] = aa
                    idx += 1
        
        # Show first 10
        for codon in all_codons[:10]:
            mu, tai = codon_properties[codon]
            obs_aa = observed_codon_to_aa[codon]
            null_aa = null_codon_to_aa[codon]
            changed = "← CHANGED" if obs_aa != null_aa else ""
            print(f"  {codon}: {obs_aa} → {null_aa}  (μ={mu:.4f}, tAI={tai:.3f}) {changed}")
        
        # Verify degeneracy preserved
        null_degeneracy = defaultdict(int)
        for aa in null_codon_to_aa.values():
            null_degeneracy[aa] += 1
        
        if dict(null_degeneracy) == dict(zip(aa_list, degeneracy_list)):
            print("  ✓ Degeneracy structure PRESERVED")
        else:
            print("  ✗ ERROR: Degeneracy structure changed!")


if __name__ == "__main__":
    import os
    
    # Try to find data directory
    possible_paths = [
        'errors/codon_modes_ecoli.tsv',
        '/storage/kiran-stuff/proteostasis_law/errors/codon_modes_ecoli.tsv',
        '../errors/codon_modes_ecoli.tsv',
    ]
    
    data_path = None
    for path in possible_paths:
        if os.path.exists(path):
            data_path = path
            break
    
    if data_path is None:
        print("ERROR: Could not find codon_modes_ecoli.tsv")
        print("Please provide the path as argument or ensure file exists")
        exit(1)
    
    print(f"Using data from: {data_path}")
    
    # Load data
    df = load_codon_data(data_path)
    
    # Run verification first
    verify_null_construction(df, n_samples=3)
    
    # Generate figure
    results = create_figure4_corrected(
        data_path=data_path,
        output_prefix='Figure4_diversity_CORRECTED',
        n_null=10000
    )
    
    print("\n" + "=" * 70)
    print("FINAL SUMMARY")
    print("=" * 70)
    print(f"Observed diversity:     {results['observed_total']:.2f}")
    print(f"Null mean ± std:        {results['null_mean']:.2f} ± {results['null_std']:.2f}")
    print(f"z-score:                {results['z_score']:.2f}")
    print(f"p-value:                {results['p_value']:.4f}")
    print(f"")
    print(f"Purine mean spread:     {results['purine_mean']:.2f}")
    print(f"Pyrimidine mean spread: {results['pyrimidine_mean']:.2f}")
    print(f"Ratio:                  {results['wobble_ratio']:.2f}x")
    print(f"Wobble p-value:         {results['wobble_p']:.4f}")
