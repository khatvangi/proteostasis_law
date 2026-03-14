#!/bin/bash
# Add SVG export to all figure scripts

echo "Adding SVG export to all figures..."

# Function to add SVG save line after PDF save
add_svg_export() {
    local file=$1
    # Check if file already has .svg export
    if grep -q "\.svg" "$file"; then
        echo "  $file - already has SVG export ✓"
    else
        echo "  $file - adding SVG export..."
        # Add SVG export after each .savefig PDF line
        sed -i "s/plt\.savefig('\([^']*\)\.pdf'/plt.savefig('\1.pdf'\nplt.savefig('\1.svg'/g" "$file"
    fi
}

# Process all figure files
add_svg_export "figure1_conceptual.py"
add_svg_export "figure2_metal_ligands.py"
add_svg_export "figure3_logistic_regression.py"
add_svg_export "figure4_codon_modes.py"
add_svg_export "figure5_mode_richness.py"
add_svg_export "figure6_capacity_bound.py"
add_svg_export "supplementary_figures.py"

echo "Done! Now regenerating all figures with SVG export..."
